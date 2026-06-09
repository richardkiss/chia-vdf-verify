#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "chia-blockchain",
#   "zstd",
#   "tqdm",
# ]
# ///
"""
Benchmark VDF proof verification against a synced Chia mainnet blockchain DB.

Measures the lower bound on sync time (pure proof-checking cost) and lets you
compare the C++ chiavdf backend against the pure-Rust chia-vdf-verify backend.

Usage:
    # C++ backend (default)
    uv run scripts/bench-sync-time.py --threads 8

    # Rust backend (from local repo)
    uv run --with . scripts/bench-sync-time.py --backend rust --threads 8

    # Rust backend (from any git URL — note the git+https:// syntax)
    uv run --with "chia-vdf-verify @ git+https://github.com/your-fork/chia-vdf-verify" \
        scripts/bench-sync-time.py --backend rust --threads 8
"""

import argparse
import collections
import sqlite3
import sys
import time
from concurrent.futures import ThreadPoolExecutor, Future
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

try:
    from chia_rs import FullBlock, VDFInfo, VDFProof, ConsensusConstants
    from chia_rs.sized_bytes import bytes32, bytes100
    from chia.consensus.default_constants import DEFAULT_CONSTANTS
    from chia.types.blockchain_format.classgroup import ClassgroupElement
except ImportError as e:
    print(f"Error: chia-blockchain not found: {e}", file=sys.stderr)
    sys.exit(1)

try:
    import zstd
except ImportError:
    print("Error: zstd not found. Install with: pip install zstd", file=sys.stderr)
    sys.exit(1)

try:
    from tqdm import tqdm
    HAVE_TQDM = True
except ImportError:
    HAVE_TQDM = False

# --- backends ---

try:
    from chiavdf import create_discriminant as _cpp_disc, verify_n_wesolowski as _cpp_verify
    HAVE_CPP = True
except ImportError:
    HAVE_CPP = False

try:
    from chia_vdf_verify import create_discriminant as _rust_disc, verify_n_wesolowski as _rust_verify
    HAVE_RUST = True
except ImportError:
    HAVE_RUST = False


@dataclass
class VDFTask:
    proof: VDFProof
    info: VDFInfo
    input_el: bytes100
    height: int
    name: str


def _disc_cache_key(challenge: bytes, size_bits: int) -> tuple:
    return (bytes(challenge), size_bits)


_disc_cache: dict = {}


def get_discriminant(challenge: bytes, size_bits: int, backend: str) -> int:
    key = _disc_cache_key(challenge, size_bits)
    if key not in _disc_cache:
        fn = _cpp_disc if backend == "cpp" else _rust_disc
        _disc_cache[key] = int(fn(challenge, size_bits), 16)
    return _disc_cache[key]


def verify_task(task: VDFTask, backend: str, primes_only: bool = False) -> tuple[bool, Optional[str]]:
    try:
        disc = get_discriminant(bytes(task.info.challenge), DEFAULT_CONSTANTS.DISCRIMINANT_SIZE_BITS, backend)
        if primes_only:
            return True, None
        verify_fn = _cpp_verify if backend == "cpp" else _rust_verify
        ok = verify_fn(
            str(disc),
            bytes(task.input_el),
            bytes(task.proof.witness),
            task.info.number_of_iterations,
            DEFAULT_CONSTANTS.DISCRIMINANT_SIZE_BITS,
            task.proof.witness_type,
        )
        return bool(ok), None
    except Exception as e:
        return False, str(e)


def extract_tasks(block: FullBlock) -> list[VDFTask]:
    tasks = []
    h = int(block.height)
    el = ClassgroupElement.get_default_element().data

    for i, ss in enumerate(block.finished_sub_slots):
        p = ss.proofs
        cc = ss.challenge_chain
        rc = ss.reward_chain
        if p.challenge_chain_slot_proof is not None:
            tasks.append(VDFTask(p.challenge_chain_slot_proof, cc.challenge_chain_end_of_slot_vdf, el, h, f"cc_slot_{i}"))
        if p.reward_chain_slot_proof is not None:
            tasks.append(VDFTask(p.reward_chain_slot_proof, rc.end_of_slot_vdf, el, h, f"rc_slot_{i}"))
        if ss.infused_challenge_chain is not None and p.infused_challenge_chain_slot_proof is not None:
            tasks.append(VDFTask(p.infused_challenge_chain_slot_proof,
                                 ss.infused_challenge_chain.infused_challenge_chain_end_of_slot_vdf, el, h, f"icc_slot_{i}"))

    rc = block.reward_chain_block
    if rc.challenge_chain_sp_vdf is not None and block.challenge_chain_sp_proof is not None:
        tasks.append(VDFTask(block.challenge_chain_sp_proof, rc.challenge_chain_sp_vdf, el, h, "cc_sp"))
    if block.challenge_chain_ip_proof is not None:
        tasks.append(VDFTask(block.challenge_chain_ip_proof, rc.challenge_chain_ip_vdf, el, h, "cc_ip"))
    if rc.reward_chain_sp_vdf is not None and block.reward_chain_sp_proof is not None:
        tasks.append(VDFTask(block.reward_chain_sp_proof, rc.reward_chain_sp_vdf, el, h, "rc_sp"))
    if block.reward_chain_ip_proof is not None:
        tasks.append(VDFTask(block.reward_chain_ip_proof, rc.reward_chain_ip_vdf, el, h, "rc_ip"))
    if rc.infused_challenge_chain_ip_vdf is not None and block.infused_challenge_chain_ip_proof is not None:
        tasks.append(VDFTask(block.infused_challenge_chain_ip_proof, rc.infused_challenge_chain_ip_vdf, el, h, "icc_ip"))

    return tasks


def count_blocks(db: Path, min_height: int, max_height: Optional[int]) -> int:
    conn = sqlite3.connect(str(db))
    q = "SELECT COUNT(*) FROM full_blocks WHERE in_main_chain = 1 AND height >= ?"
    params: list = [min_height]
    if max_height is not None:
        q += " AND height <= ?"
        params.append(max_height)
    (n,) = conn.execute(q, params).fetchone()
    conn.close()
    return n


def iter_blocks(db: Path, min_height: int, max_height: Optional[int]):
    """Stream (height, FullBlock) in chain order — O(1) memory."""
    conn = sqlite3.connect(str(db))
    conn.row_factory = sqlite3.Row
    q = "SELECT height, block FROM full_blocks WHERE in_main_chain = 1 AND height >= ? ORDER BY height"
    params: list = [min_height]
    if max_height is not None:
        q = q.replace("ORDER BY", "AND height <= ? ORDER BY")
        params.insert(1, max_height)

    for row in conn.execute(q, params):
        try:
            yield int(row["height"]), FullBlock.from_bytes(zstd.decompress(row["block"]))
        except Exception as e:
            print(f"Warning: height {row['height']}: {e}", file=sys.stderr)
    conn.close()


def iter_tasks(db: Path, min_height: int, max_height: Optional[int]):
    """Stream (height, VDFTask) pairs in chain order."""
    for height, block in iter_blocks(db, min_height, max_height):
        for task in extract_tasks(block):
            yield task


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", type=Path,
                        default=Path.home() / ".chia/mainnet/db/blockchain_v2_mainnet.sqlite")
    parser.add_argument("--backend", choices=["cpp", "rust"], default="cpp",
                        help="cpp=chiavdf (default), rust=chia-vdf-verify")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--start-height", type=int, default=0, metavar="N",
                        help="First block height to process (default: 0)")
    parser.add_argument("--max-height", type=int, metavar="N",
                        help="Last block height to process (default: chain tip)")
    parser.add_argument("--primes-only", action="store_true",
                        help="Only compute discriminants (prime finding), skip proof verification")
    args = parser.parse_args()

    if not args.db.exists():
        sys.exit(f"Database not found: {args.db}")
    if args.backend == "cpp" and not HAVE_CPP:
        sys.exit("C++ backend not available: pip install chiavdf")
    if args.backend == "rust" and not HAVE_RUST:
        sys.exit(
            "Rust backend not available. Run with one of:\n"
            "  uv run --with . scripts/bench-sync-time.py --backend rust\n"
            '  uv run --with "chia-vdf-verify @ git+https://github.com/OWNER/chia-vdf-verify"'
            " scripts/bench-sync-time.py --backend rust"
        )

    nblocks = count_blocks(args.db, args.start_height, args.max_height)
    height_range = f"{args.start_height:,} – {args.max_height:,}" if args.max_height else f"{args.start_height:,} – tip"
    backend_label = "chiavdf C++" if args.backend == "cpp" else "chia-vdf-verify Rust"
    mode_label = "primes only (no verification)" if args.primes_only else "full verification"
    print(f"Blocks in range [{height_range}]: {nblocks:,}")
    print(f"Backend: {backend_label}   Threads: {args.threads}   Mode: {mode_label}\n")

    ok = fail = 0
    blocks_seen = 0
    last_height = args.start_height
    errors: list[tuple[int, str]] = []
    t0 = time.perf_counter()

    # Sliding window: keep at most `window` futures in flight so we never
    # buffer the whole chain's worth of tasks in memory.
    window = args.threads * 4
    in_flight: collections.deque[tuple[Future, VDFTask]] = collections.deque()

    progress = tqdm(total=nblocks, desc="blocks", unit="block") if HAVE_TQDM else None
    have_progress = progress is not None
    prev_height = -1

    def drain_one() -> None:
        nonlocal ok, fail
        fut, task = in_flight.popleft()
        success, err = fut.result()
        if success:
            ok += 1
        else:
            fail += 1
            if err:
                errors.append((task.height, err))
                if len(errors) <= 5:
                    print(f"\n  error at height {task.height}: {err}", file=sys.stderr)

    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        for task in iter_tasks(args.db, args.start_height, args.max_height):
            if have_progress and task.height != prev_height:
                progress.update(1)
                progress.set_postfix(height=f"{task.height:,}", proofs=ok + fail)
                prev_height = task.height
            last_height = task.height

            if len(in_flight) >= window:
                drain_one()
            in_flight.append((ex.submit(verify_task, task, args.backend, args.primes_only), task))

        while in_flight:
            drain_one()

    if have_progress:
        progress.close()

    elapsed = time.perf_counter() - t0
    total = ok + fail
    pps = total / elapsed if elapsed else 0
    single_thread_est = elapsed * args.threads

    print(f"\n{'='*60}")
    print(f"Heights:         {args.start_height:>10,} – {last_height:,}")
    print(f"Proofs checked:  {total:>10,}  ({ok:,} ok, {fail:,} failed)")
    print(f"Wall time:       {elapsed:>10.2f}s")
    print(f"Proofs/sec:      {pps:>10.1f}")
    print(f"Single-thread ≈  {single_thread_est:>10.1f}s  (lower bound on sequential sync)")
    print(f"Backend:         {args.backend:>10}  ({backend_label})")
    print(f"Mode:            {'primes only':>10}" if args.primes_only else f"Mode:            {'full verify':>10}")
    print(f"Threads:         {args.threads:>10}")
    print(f"{'='*60}")

    if len(errors) > 5:
        print(f"({len(errors) - 5} more errors suppressed)", file=sys.stderr)


if __name__ == "__main__":
    main()
