#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "chia-rs>=0.38",
#   "zstd>=1.5",
# ]
# ///
"""Differential fuzzing: compare chia_vdf_verify (Rust) vs chiavdf (C++) on real and mutated proofs.

Extracts VDF proofs from the Chia mainnet SQLite DB, verifies each with both
verifiers to confirm agreement on valid proofs, then mutates each proof and
checks that both verifiers agree on the rejection verdict.

Usage:
    python scripts/differential_fuzz.py
    python scripts/differential_fuzz.py --db /path/to/db.sqlite
    python scripts/differential_fuzz.py --sample 200 --mutations 20
    python scripts/differential_fuzz.py --seed 42

Environment: must have both chiavdf and chia_vdf_verify importable.
Recommended: run under /home/kiss/projects/chia-blockchain/.venv
"""
from __future__ import annotations

import argparse
import random
import sqlite3
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import zstd
from chia_rs import FullBlock
from chia.consensus.default_constants import DEFAULT_CONSTANTS
from chia.types.blockchain_format.classgroup import ClassgroupElement

try:
    import chiavdf as _chiavdf_module
    from chiavdf import create_discriminant as _chiavdf_disc
    from chiavdf import verify_n_wesolowski as _chiavdf_verify
    CHIAVDF_OK = True
except ImportError:
    CHIAVDF_OK = False
    print("WARNING: chiavdf not importable — only testing chia_vdf_verify", file=sys.stderr)

try:
    import chia_vdf_verify as _rust_module
    from chia_vdf_verify import create_discriminant as _rust_disc
    from chia_vdf_verify import verify_n_wesolowski as _rust_verify
    RUST_OK = True
except ImportError:
    RUST_OK = False
    print("ERROR: chia_vdf_verify not importable", file=sys.stderr)
    sys.exit(1)

IDENTITY = ClassgroupElement.get_default_element()
CONSTANTS = DEFAULT_CONSTANTS


@dataclass
class ProofCase:
    height: int
    label: str
    disc_str: str       # decimal string, negative
    input_el: bytes     # 100-byte input form
    output: bytes       # output_form_bytes + witness_bytes
    n_iters: int
    disc_size: int
    witness_type: int


@dataclass
class Stats:
    valid_proofs: int = 0
    valid_agree: int = 0
    valid_disagree: int = 0
    mutations_tested: int = 0
    mutations_both_reject: int = 0
    mutations_rust_only_accept: int = 0
    mutations_chiavdf_only_accept: int = 0
    mutations_both_accept: int = 0
    rust_crashes: int = 0
    chiavdf_crashes: int = 0
    disagreements: list[str] = field(default_factory=list)


def call_rust(disc: str, input_el: bytes, output: bytes, n_iters: int, disc_size: int, witness_type: int) -> Optional[bool]:
    try:
        return bool(_rust_verify(disc, input_el, output, n_iters, disc_size, witness_type))
    except Exception:
        return None  # unexpected exception → treated as rejection


def call_chiavdf(disc: str, input_el: bytes, output: bytes, n_iters: int, disc_size: int, witness_type: int) -> Optional[bool]:
    if not CHIAVDF_OK:
        return None
    try:
        return bool(_chiavdf_verify(disc, input_el, output, n_iters, disc_size, witness_type))
    except Exception:
        # chiavdf raises exceptions (bqfc_export overflow etc.) on malformed inputs
        # treat as rejection — not a crash in the traditional sense
        return None


def extract_proofs(db_path: str, sample_size: int, rng: random.Random) -> list[ProofCase]:
    """Extract up to sample_size valid VDF proofs from the blockchain DB."""
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    conn.execute("pragma query_only = ON")

    total_rows = conn.execute(
        "SELECT COUNT(*) FROM full_blocks WHERE in_main_chain = 1"
    ).fetchone()[0]

    # Sample heights across the entire chain
    max_height = conn.execute(
        "SELECT MAX(height) FROM full_blocks WHERE in_main_chain = 1"
    ).fetchone()[0]
    min_height = conn.execute(
        "SELECT MIN(height) FROM full_blocks WHERE in_main_chain = 1"
    ).fetchone()[0]

    print(f"DB has {total_rows} main-chain blocks, heights {min_height}..{max_height}")

    # Build a list of candidate heights spread across the chain
    step = max(1, (max_height - min_height) // (sample_size * 3))
    candidates = list(range(min_height, max_height + 1, step))
    rng.shuffle(candidates)
    candidates = candidates[: sample_size * 5]  # oversample to find enough proofs

    cases: list[ProofCase] = []

    disc_size = CONSTANTS.DISCRIMINANT_SIZE_BITS

    for height in sorted(candidates):
        if len(cases) >= sample_size:
            break
        row = conn.execute(
            "SELECT block FROM full_blocks WHERE in_main_chain = 1 AND height = ?",
            (height,),
        ).fetchone()
        if row is None:
            continue

        try:
            block = FullBlock.from_bytes(zstd.decompress(row[0]))
        except Exception:
            continue

        rc = block.reward_chain_block

        def _disc_str(challenge_bytes: bytes) -> str:
            return str(int(_rust_disc(challenge_bytes, disc_size), 16))

        candidate_proofs: list[ProofCase] = []

        if block.challenge_chain_sp_proof is not None and rc.challenge_chain_sp_vdf is not None:
            p = block.challenge_chain_sp_proof
            if p.normalized_to_identity:
                info = rc.challenge_chain_sp_vdf
                candidate_proofs.append(ProofCase(
                    height=height,
                    label="cc_sp",
                    disc_str=_disc_str(bytes(info.challenge)),
                    input_el=bytes(IDENTITY.data),
                    output=bytes(info.output.data) + bytes(p.witness),
                    n_iters=int(info.number_of_iterations),
                    disc_size=disc_size,
                    witness_type=int(p.witness_type),
                ))

        if block.challenge_chain_ip_proof is not None:
            p = block.challenge_chain_ip_proof
            if p.normalized_to_identity:
                info = rc.challenge_chain_ip_vdf
                candidate_proofs.append(ProofCase(
                    height=height,
                    label="cc_ip",
                    disc_str=_disc_str(bytes(info.challenge)),
                    input_el=bytes(IDENTITY.data),
                    output=bytes(info.output.data) + bytes(p.witness),
                    n_iters=int(info.number_of_iterations),
                    disc_size=disc_size,
                    witness_type=int(p.witness_type),
                ))

        if block.reward_chain_sp_proof is not None and rc.reward_chain_sp_vdf is not None:
            p = block.reward_chain_sp_proof
            info = rc.reward_chain_sp_vdf
            candidate_proofs.append(ProofCase(
                height=height,
                label="rc_sp",
                disc_str=_disc_str(bytes(info.challenge)),
                input_el=bytes(IDENTITY.data),
                output=bytes(info.output.data) + bytes(p.witness),
                n_iters=int(info.number_of_iterations),
                disc_size=disc_size,
                witness_type=int(p.witness_type),
            ))

        if block.reward_chain_ip_proof is not None:
            p = block.reward_chain_ip_proof
            info = rc.reward_chain_ip_vdf
            candidate_proofs.append(ProofCase(
                height=height,
                label="rc_ip",
                disc_str=_disc_str(bytes(info.challenge)),
                input_el=bytes(IDENTITY.data),
                output=bytes(info.output.data) + bytes(p.witness),
                n_iters=int(info.number_of_iterations),
                disc_size=disc_size,
                witness_type=int(p.witness_type),
            ))

        if block.infused_challenge_chain_ip_proof is not None and rc.infused_challenge_chain_ip_vdf is not None:
            p = block.infused_challenge_chain_ip_proof
            if p.normalized_to_identity:
                info = rc.infused_challenge_chain_ip_vdf
                candidate_proofs.append(ProofCase(
                    height=height,
                    label="icc_ip",
                    disc_str=_disc_str(bytes(info.challenge)),
                    input_el=bytes(IDENTITY.data),
                    output=bytes(info.output.data) + bytes(p.witness),
                    n_iters=int(info.number_of_iterations),
                    disc_size=disc_size,
                    witness_type=int(p.witness_type),
                ))

        cases.extend(candidate_proofs)

    conn.close()
    return cases[:sample_size]


def mutate_bytes(data: bytes, rng: random.Random, n_flips: int = 1) -> bytes:
    """Flip n_flips random bytes at random positions."""
    arr = bytearray(data)
    for _ in range(n_flips):
        pos = rng.randrange(len(arr))
        arr[pos] ^= rng.randint(1, 255)  # guaranteed flip (no-op avoided)
    return bytes(arr)


def run_valid_check(case: ProofCase, stats: Stats) -> None:
    """Verify a known-valid proof; check both verifiers agree on True."""
    rust_result = call_rust(
        case.disc_str, case.input_el, case.output,
        case.n_iters, case.disc_size, case.witness_type,
    )
    chiavdf_result = call_chiavdf(
        case.disc_str, case.input_el, case.output,
        case.n_iters, case.disc_size, case.witness_type,
    ) if CHIAVDF_OK else None

    stats.valid_proofs += 1

    rust_ok = rust_result is True
    chiavdf_ok = chiavdf_result is True or not CHIAVDF_OK

    if rust_result is None:
        stats.rust_crashes += 1
    if chiavdf_result is None and CHIAVDF_OK:
        stats.chiavdf_crashes += 1

    if CHIAVDF_OK and rust_result != chiavdf_result:
        msg = (
            f"VALID PROOF DISAGREE height={case.height} label={case.label} "
            f"rust={rust_result} chiavdf={chiavdf_result}"
        )
        stats.disagreements.append(msg)
        stats.valid_disagree += 1
        print(f"  !! {msg}", flush=True)
    else:
        stats.valid_agree += 1
        if not rust_ok:
            # Both rejected a supposedly-valid proof — report but not a "disagreement"
            msg = (
                f"BOTH REJECTED valid proof height={case.height} label={case.label} "
                f"rust={rust_result} chiavdf={chiavdf_result}"
            )
            stats.disagreements.append(msg)
            print(f"  !! {msg}", flush=True)


def run_mutation(case: ProofCase, mutated_output: bytes, stats: Stats, tag: str) -> None:
    """Run both verifiers on a mutated proof and record agreement."""
    rust_result = call_rust(
        case.disc_str, case.input_el, mutated_output,
        case.n_iters, case.disc_size, case.witness_type,
    )
    chiavdf_result = call_chiavdf(
        case.disc_str, case.input_el, mutated_output,
        case.n_iters, case.disc_size, case.witness_type,
    ) if CHIAVDF_OK else None

    if rust_result is None:
        stats.rust_crashes += 1
    if chiavdf_result is None and CHIAVDF_OK:
        stats.chiavdf_crashes += 1

    stats.mutations_tested += 1

    if CHIAVDF_OK:
        if rust_result is True and chiavdf_result is not True:
            stats.mutations_rust_only_accept += 1
            msg = (
                f"MUTATION DISAGREE (rust=True, chiavdf={chiavdf_result}) "
                f"height={case.height} label={case.label} mutation={tag}"
            )
            stats.disagreements.append(msg)
            print(f"  !! {msg}", flush=True)
        elif rust_result is not True and chiavdf_result is True:
            stats.mutations_chiavdf_only_accept += 1
            msg = (
                f"MUTATION DISAGREE (rust={rust_result}, chiavdf=True) "
                f"height={case.height} label={case.label} mutation={tag}"
            )
            stats.disagreements.append(msg)
            print(f"  !! {msg}", flush=True)
        elif rust_result is True and chiavdf_result is True:
            stats.mutations_both_accept += 1
        else:
            stats.mutations_both_reject += 1
    else:
        # No chiavdf — just count rust rejections
        if rust_result is False or rust_result is None:
            stats.mutations_both_reject += 1
        else:
            stats.mutations_both_accept += 1


def main() -> None:
    default_db = str(Path.home() / ".chia/mainnet/db/blockchain_v2_mainnet.sqlite")

    parser = argparse.ArgumentParser(description="Differential VDF fuzzer: Rust vs chiavdf")
    parser.add_argument("--db", default=default_db, help="Path to blockchain_v2 SQLite DB")
    parser.add_argument("--sample", type=int, default=100, help="Number of valid proofs to sample")
    parser.add_argument("--mutations", type=int, default=20, help="Mutations per proof")
    parser.add_argument("--seed", type=int, default=0, help="Random seed")
    parser.add_argument("--verbose", action="store_true", help="Print each proof check")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    stats = Stats()
    t0 = time.time()

    print("=" * 60)
    print("Differential VDF Fuzzer")
    print(f"  DB:         {args.db}")
    print(f"  Sample:     {args.sample} proofs")
    print(f"  Mutations:  {args.mutations} per proof")
    print(f"  Seed:       {args.seed}")
    print(f"  chiavdf:    {'available' if CHIAVDF_OK else 'NOT AVAILABLE (single-verifier mode)'}")
    print(f"  rust vdf:   available")
    print("=" * 60)
    sys.stdout.flush()

    print(f"\nExtracting proof cases from DB...")
    cases = extract_proofs(args.db, args.sample, rng)
    print(f"Extracted {len(cases)} proof cases\n")

    if not cases:
        print("ERROR: no proof cases found", file=sys.stderr)
        sys.exit(1)

    # Phase 1: valid proof checks
    print(f"Phase 1: Checking {len(cases)} valid proofs...")
    for i, case in enumerate(cases):
        run_valid_check(case, stats)
        if args.verbose or (i + 1) % 20 == 0:
            print(f"  [{i+1}/{len(cases)}] height={case.height} label={case.label} "
                  f"agree={stats.valid_agree} disagree={stats.valid_disagree}")
            sys.stdout.flush()

    elapsed1 = time.time() - t0
    print(f"\nPhase 1 done in {elapsed1:.1f}s")
    print(f"  Valid proofs: {stats.valid_proofs}")
    print(f"  Agree:        {stats.valid_agree}")
    print(f"  Disagree:     {stats.valid_disagree}")
    print(f"  Rust exceptions: {stats.rust_crashes}")
    if CHIAVDF_OK:
        print(f"  Chiavdf exceptions: {stats.chiavdf_crashes}")
    sys.stdout.flush()

    # Phase 2: mutation testing
    print(f"\nPhase 2: Mutation testing ({len(cases)} proofs × {args.mutations} mutations)...")
    t2 = time.time()
    total_mutations = len(cases) * args.mutations

    for i, case in enumerate(cases):
        for m in range(args.mutations):
            # Vary mutation intensity: mostly 1 flip, sometimes 2-4
            n_flips = rng.choices([1, 2, 4, 8], weights=[60, 25, 10, 5])[0]
            mutated = mutate_bytes(case.output, rng, n_flips)
            tag = f"m{m}_f{n_flips}"
            run_mutation(case, mutated, stats, tag)

        if (i + 1) % 10 == 0:
            done = (i + 1) * args.mutations
            elapsed = time.time() - t2
            rate = done / elapsed if elapsed > 0 else 0
            print(
                f"  [{done}/{total_mutations}] "
                f"both_reject={stats.mutations_both_reject} "
                f"rust_only_accept={stats.mutations_rust_only_accept} "
                f"chiavdf_only_accept={stats.mutations_chiavdf_only_accept} "
                f"both_accept={stats.mutations_both_accept} "
                f"rate={rate:.0f}/s"
            )
            sys.stdout.flush()

    elapsed2 = time.time() - t2
    total_elapsed = time.time() - t0

    print(f"\nPhase 2 done in {elapsed2:.1f}s")
    print()
    print("=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print(f"Valid proofs sampled:      {stats.valid_proofs}")
    print(f"  Verifier agreement:      {stats.valid_agree}")
    print(f"  Verifier disagreement:   {stats.valid_disagree}")
    print()
    print(f"Mutations tested:          {stats.mutations_tested}")
    print(f"  Both rejected:           {stats.mutations_both_reject} "
          f"({100*stats.mutations_both_reject/max(1,stats.mutations_tested):.1f}%)")
    if CHIAVDF_OK:
        print(f"  Rust=True, chiavdf=False: {stats.mutations_rust_only_accept}  ← BUG if > 0")
        print(f"  Rust=False, chiavdf=True: {stats.mutations_chiavdf_only_accept}  ← BUG if > 0")
    print(f"  Both accepted (surprise): {stats.mutations_both_accept}  ← unexpected if > 0")
    print()
    print(f"Rust exceptions:           {stats.rust_crashes}")
    if CHIAVDF_OK:
        print(f"Chiavdf exceptions:        {stats.chiavdf_crashes}  (bqfc_export overflow on malformed inputs, counted as rejection)")
    print(f"Total time:                {total_elapsed:.1f}s")
    print("=" * 60)

    total_disagree = stats.valid_disagree + stats.mutations_rust_only_accept + stats.mutations_chiavdf_only_accept
    if total_disagree == 0 and stats.rust_crashes == 0:
        print("\nOK: No disagreements found. Rust verifier matches chiavdf.")
    else:
        print(f"\nFAIL: {total_disagree} disagreements found!")
        for d in stats.disagreements[:20]:
            print(f"  {d}")
        if len(stats.disagreements) > 20:
            print(f"  ... and {len(stats.disagreements) - 20} more")

    sys.exit(0 if total_disagree == 0 else 1)


if __name__ == "__main__":
    main()
