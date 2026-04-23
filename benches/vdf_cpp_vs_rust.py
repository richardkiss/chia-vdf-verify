#!/usr/bin/env python3
"""
Real-world C++ chiavdf vs Rust chia_vdf_verify benchmark.

Extracts actual VDF proofs from the mainnet blockchain DB (end-of-slot proofs,
which always use the identity element as input and are self-contained), then
measures:
  1. Single-threaded latency: chiavdf vs chia_vdf_verify
  2. Multi-threaded throughput: chiavdf (GIL-held) vs chia_vdf_verify (GIL-released)

Usage:
  python benches/vdf_cpp_vs_rust.py [--db PATH] [--proofs N] [--threads T]
"""
from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import zstd
import chia_rs
import chiavdf
from chia_vdf_verify import (
    create_discriminant_bytes,
    verify_n_wesolowski_bytes,
)

DISC_BITS = 1024
# Identity element for class group (input to EOS VDFs)
IDENTITY = bytes([0x08]) + bytes(99)


@dataclass
class VDFCase:
    challenge: bytes          # 32-byte challenge seed
    proof_blob: bytes         # output_form (100 bytes) + witness bytes
    iters: int
    witness_type: int
    height: int

    # Pre-computed discriminant representations (set after extraction)
    disc_str: str = ""        # decimal string for chiavdf
    disc_bytes: bytes = b""   # sign+magnitude bytes for chia_vdf_verify


# ── Extraction ─────────────────────────────────────────────────────────────────

def extract_proofs(db_path: str, n: int) -> list[VDFCase]:
    """Pull end-of-slot CC VDF proofs from the block store."""
    conn = sqlite3.connect(db_path)
    cases: list[VDFCase] = []
    # Walk blocks from the tip backwards to find EOS blocks quickly
    rows = conn.execute(
        "SELECT block FROM full_blocks ORDER BY height DESC LIMIT 5000"
    ).fetchall()
    for (raw,) in rows:
        block = chia_rs.FullBlock.from_bytes(zstd.decompress(raw))
        for ss in (block.finished_sub_slots or []):
            vdf = ss.challenge_chain.challenge_chain_end_of_slot_vdf
            proof = ss.proofs.challenge_chain_slot_proof
            blob = vdf.output.data + bytes(proof.witness)
            cases.append(VDFCase(
                challenge=bytes(vdf.challenge),
                proof_blob=blob,
                iters=vdf.number_of_iterations,
                witness_type=proof.witness_type,
                height=block.height,
            ))
        if len(cases) >= n:
            break
    conn.close()
    if not cases:
        sys.exit("No EOS VDF proofs found in DB")
    return cases[:n]


def precompute_discriminants(cases: list[VDFCase]) -> None:
    """Pre-generate discriminants so timing doesn't include disc creation."""
    seen: dict[bytes, tuple[str, bytes]] = {}
    for c in cases:
        if c.challenge not in seen:
            disc_hex = chiavdf.create_discriminant(c.challenge, DISC_BITS)
            disc_int = int(disc_hex, 16)
            disc_str = str(disc_int)
            disc_bytes = create_discriminant_bytes(c.challenge, DISC_BITS)
            seen[c.challenge] = (disc_str, disc_bytes)
        c.disc_str, c.disc_bytes = seen[c.challenge]


# ── Verification wrappers ──────────────────────────────────────────────────────

def verify_cpp(c: VDFCase) -> bool:
    return chiavdf.verify_n_wesolowski(
        c.disc_str, IDENTITY, c.proof_blob,
        c.iters, DISC_BITS, c.witness_type,
    )


def verify_rust(c: VDFCase) -> bool:
    return verify_n_wesolowski_bytes(
        c.disc_bytes, IDENTITY, c.proof_blob,
        c.iters, DISC_BITS, c.witness_type,
    )


# ── Benchmarks ─────────────────────────────────────────────────────────────────

def bench_single_threaded(cases: list[VDFCase], rounds: int = 3) -> dict:
    """Run each verifier over all cases, return median throughput."""
    results = {}
    for label, fn in [("chiavdf (C++)", verify_cpp), ("chia_vdf_verify (Rust)", verify_rust)]:
        times = []
        for _ in range(rounds):
            t0 = time.perf_counter()
            ok = sum(1 for c in cases if fn(c))
            elapsed = time.perf_counter() - t0
            times.append(elapsed)
        times.sort()
        median = times[len(times) // 2]
        results[label] = {
            "total_s": median,
            "per_proof_ms": median / len(cases) * 1000,
            "throughput": len(cases) / median,
            "verified": ok,
        }
    return results


def bench_multithreaded(cases: list[VDFCase], threads: int, per_thread: int = 10) -> dict:
    """Submit all cases to a thread pool and measure wall-clock throughput."""
    # Ensure at least per_thread proofs per thread
    total = threads * per_thread
    workload = (cases * (total // len(cases) + 1))[:total]
    results = {}
    for label, fn in [("chiavdf (C++)", verify_cpp), ("chia_vdf_verify (Rust)", verify_rust)]:
        with ThreadPoolExecutor(max_workers=threads) as pool:
            t0 = time.perf_counter()
            futs = [pool.submit(fn, c) for c in workload]
            ok = sum(1 for f in as_completed(futs) if f.result())
            elapsed = time.perf_counter() - t0
        results[label] = {
            "total_s": elapsed,
            "throughput": len(workload) / elapsed,
            "verified": ok,
            "n": len(workload),
        }
    return results


# ── Reporting ──────────────────────────────────────────────────────────────────

def fmt(n: float, unit: str = "") -> str:
    return f"{n:.2f}{unit}"


def report(cases: list[VDFCase], threads: int) -> None:
    print()
    print("=" * 60)
    print("VDF Benchmark: C++ chiavdf vs Rust chia_vdf_verify")
    print("=" * 60)
    print(f"Proofs: {len(cases)}  |  Threads: {threads}")
    iters_set = sorted({c.iters for c in cases})
    wt_set = sorted({c.witness_type for c in cases})
    print(f"Iterations: {iters_set}  |  Witness types: {wt_set}")
    print(f"Height range: {min(c.height for c in cases):,} – {max(c.height for c in cases):,}")
    print()

    print("── Single-threaded ─────────────────────────────────────")
    st = bench_single_threaded(cases)
    for label, r in st.items():
        print(f"  {label}")
        print(f"    {fmt(r['per_proof_ms'], ' ms/proof')}  |  "
              f"{fmt(r['throughput'], ' proofs/s')}  |  "
              f"verified {r['verified']}/{len(cases)}")
    labels = list(st.keys())
    if len(labels) == 2:
        ratio = st[labels[0]]["per_proof_ms"] / st[labels[1]]["per_proof_ms"]
        faster = labels[0] if ratio < 1 else labels[1]
        print(f"  → C++ is {max(ratio, 1/ratio):.2f}x {'faster' if ratio < 1 else 'slower'} single-threaded")
    print()

    print(f"── Multi-threaded ({threads} threads) ─────────────────────────")
    mt = bench_multithreaded(cases, threads)
    for label, r in mt.items():
        print(f"  {label}")
        print(f"    {fmt(r['throughput'], ' proofs/s')} wall-clock  "
              f"({r['n']} proofs in {fmt(r['total_s'], 's')})")
    labels = list(mt.keys())
    if len(labels) == 2:
        ratio = mt[labels[1]]["throughput"] / mt[labels[0]]["throughput"]
        print(f"  → Rust is {ratio:.2f}x {'faster' if ratio > 1 else 'slower'} multi-threaded "
              f"(GIL release vs GIL hold)")
    print()

    # Scaling efficiency
    print("── GIL effect (multi/single throughput ratio) ──────────")
    for label in labels:
        st_tp = st[label]["throughput"]
        mt_tp = mt[label]["throughput"]
        scale = mt_tp / st_tp
        print(f"  {label}: {scale:.1f}x scaling with {threads} threads "
              f"(ideal: {threads}x)")
    print()


# ── Entry point ────────────────────────────────────────────────────────────────

def main() -> None:
    default_db = os.path.expanduser(
        "~/.chia/rust-vdf-test/db/blockchain_v2_mainnet.sqlite"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=default_db, help="Path to blockchain SQLite DB")
    parser.add_argument("--proofs", type=int, default=20, help="Number of EOS proofs to extract")
    parser.add_argument("--threads", type=int, default=8, help="Thread count for multi-threaded bench")
    args = parser.parse_args()

    if not Path(args.db).exists():
        sys.exit(f"DB not found: {args.db}")

    print(f"Extracting {args.proofs} EOS VDF proofs from {args.db} ...")
    cases = extract_proofs(args.db, args.proofs)
    print(f"  Got {len(cases)} proofs. Pre-computing discriminants ...")
    precompute_discriminants(cases)
    print("  Done.")

    report(cases, args.threads)


if __name__ == "__main__":
    main()
