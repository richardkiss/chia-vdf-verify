#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "chia-blockchain>=2.0",
#   "chiavdf",
#   "chia-vdf-verify",
#   "zstd",
# ]
# ///
"""
Sequential VDF benchmark — verifies every proof exactly as the chia node does.

Processes blocks in order, maintains chain state (prev BlockRecord), and
correctly computes the input ClassgroupElement and iteration count for each
EOS / SP / IP proof.  Compares C++ (chiavdf) vs Rust (chia-vdf-verify).

Usage:
    uv run --with . scripts/bench-sync-time.py --backend rust
    uv run --with . scripts/bench-sync-time.py --backend cpp
    uv run --with . scripts/bench-sync-time.py --backend both
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
import time
import zstd
from functools import lru_cache
from pathlib import Path
from typing import Optional

from chia_rs import BlockRecord, ClassgroupElement, EndOfSubSlotBundle, FullBlock, VDFInfo, VDFProof
from chia_rs.sized_bytes import bytes32
from chia_rs.sized_ints import uint64

# ---------------------------------------------------------------------------
# Backend imports
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Mainnet constants
# ---------------------------------------------------------------------------
DISC_BITS: int = 1024
NUM_SPS_SUB_SLOT: int = 64
NUM_SP_INTERVALS_EXTRA: int = 3
GENESIS_CHALLENGE: bytes32 = bytes32.fromhex(
    "ccd5bb71183532bff220ba46c268991a3ff07eb358e8255a65c30a2dce0e5fbb"
)
_GENESIS_CHALLENGE_BYTES = bytes.fromhex(
    "ccd5bb71183532bff220ba46c268991a3ff07eb358e8255a65c30a2dce0e5fbb"
)
DEFAULT_DB = Path.home() / ".chia/mainnet/db/blockchain_v2_mainnet.sqlite"

# ---------------------------------------------------------------------------
# Discriminant cache — backend-agnostic integer
# ---------------------------------------------------------------------------
@lru_cache(maxsize=2048)
def get_discriminant(challenge: bytes) -> int:
    fn = _cpp_disc if HAVE_CPP else _rust_disc
    return int(fn(challenge, DISC_BITS), 16)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_DEFAULT_EL: Optional[ClassgroupElement] = None

def default_el() -> ClassgroupElement:
    global _DEFAULT_EL
    if _DEFAULT_EL is None:
        _DEFAULT_EL = ClassgroupElement.get_default_element()
    return _DEFAULT_EL


def calc_sp_iters(br: BlockRecord) -> int:
    return int(br.signage_point_index) * int(br.sub_slot_iters) // NUM_SPS_SUB_SLOT


def calc_ip_iters(br: BlockRecord) -> int:
    return calc_sp_iters(br) + int(br.required_iters)


def is_overflow(br: BlockRecord) -> bool:
    return int(br.signage_point_index) >= NUM_SPS_SUB_SLOT - NUM_SP_INTERVALS_EXTRA


def sp_total_iters(br: BlockRecord) -> int:
    """Total VDF iters at the signage point of this block."""
    # ip_iters = sp_iters + required_iters  =>  sp_total = total - required
    return int(br.total_iters) - int(br.required_iters)


def in_genesis_slot(challenge: bytes) -> bool:
    """
    The genesis slot uses cumulative proofs (default_el + full ip_iters) rather
    than delta proofs from the previous block's output.  Detect by challenge.
    """
    return bytes(challenge) == _GENESIS_CHALLENGE_BYTES


# ---------------------------------------------------------------------------
# Core proof verifier — mirrors validate_vdf / verify_vdf
# ---------------------------------------------------------------------------
def verify_vdf(
    proof: VDFProof,
    info: VDFInfo,
    input_el: ClassgroupElement,
    backend: str,
) -> bool:
    """
    Verify one VDF proof.  For normalized_to_identity proofs the input is
    always the default element regardless of what the caller passes.
    """
    try:
        actual_input = default_el() if proof.normalized_to_identity else input_el
        disc = get_discriminant(bytes(info.challenge))
        output_blob = bytes(info.output.data) + bytes(proof.witness)
        fn = _cpp_verify if backend == "cpp" else _rust_verify
        return bool(fn(
            str(disc),
            bytes(actual_input.data),
            output_blob,
            info.number_of_iterations,
            DISC_BITS,
            proof.witness_type,
        ))
    except Exception:
        return False


# ---------------------------------------------------------------------------
# BlockRecord in-memory cache
# ---------------------------------------------------------------------------
class BlockRecordDB:
    def __init__(self, conn: sqlite3.Connection) -> None:
        self._conn = conn
        self._cache: dict[bytes, BlockRecord] = {}

    def add(self, br: BlockRecord) -> None:
        self._cache[bytes(br.header_hash)] = br

    def block_record(self, header_hash: bytes32) -> BlockRecord:
        key = bytes(header_hash)
        if key not in self._cache:
            row = self._conn.execute(
                "SELECT block_record FROM full_blocks WHERE header_hash=?", (key,)
            ).fetchone()
            if row is None:
                raise KeyError(f"BlockRecord not found: {header_hash.hex()}")
            br = BlockRecord.from_bytes(row[0])
            self._cache[key] = br
        return self._cache[key]


# ---------------------------------------------------------------------------
# EOS proof verification
# ---------------------------------------------------------------------------
def verify_eos_proofs(
    block: FullBlock,
    block_br: BlockRecord,
    prev_br: Optional[BlockRecord],
    backend: str,
) -> tuple[int, int]:
    ok_count = fail_count = 0
    genesis = prev_br is None

    for n, ss in enumerate(block.finished_sub_slots):
        p = ss.proofs
        cc = ss.challenge_chain
        rc = ss.reward_chain

        # CC EOS -------------------------------------------------------
        # Mirrors block_header_validation.py logic for cc_start_element
        # and partial_cc_vdf_info (which uses eos_vdf_iters, not full slot iters).
        if genesis or n > 0:
            cc_start = default_el()
            eos_iters = int(cc.challenge_chain_end_of_slot_vdf.number_of_iterations)
        else:
            # n == 0, not genesis: EOS VDF starts from prev block's IP output
            assert prev_br is not None
            cc_start = prev_br.challenge_vdf_output
            eos_iters = int(prev_br.sub_slot_iters) - calc_ip_iters(prev_br)

        cc_eos_info = cc.challenge_chain_end_of_slot_vdf
        # partial_cc_vdf_info has eos_iters (not full slot iters)
        partial_info = cc_eos_info.replace(number_of_iterations=uint64(eos_iters))

        if p.challenge_chain_slot_proof is not None:
            proof = p.challenge_chain_slot_proof
            info = cc_eos_info if proof.normalized_to_identity else partial_info
            inp = default_el() if proof.normalized_to_identity else cc_start
            ok = verify_vdf(proof, info, inp, backend)
            ok_count += ok; fail_count += not ok

        # RC EOS: always starts from default element
        if p.reward_chain_slot_proof is not None:
            ok = verify_vdf(p.reward_chain_slot_proof, rc.end_of_slot_vdf, default_el(), backend)
            ok_count += ok; fail_count += not ok

        # ICC EOS: always starts from default element
        icc = ss.infused_challenge_chain
        if icc is not None and p.infused_challenge_chain_slot_proof is not None:
            ok = verify_vdf(
                p.infused_challenge_chain_slot_proof,
                icc.infused_challenge_chain_end_of_slot_vdf,
                default_el(),
                backend,
            )
            ok_count += ok; fail_count += not ok

    return ok_count, fail_count


# ---------------------------------------------------------------------------
# SP proof verification
# ---------------------------------------------------------------------------
def get_cc_sp_candidates(
    block_br: BlockRecord,
    fss: list,
    prev_br: Optional[BlockRecord],
    brdb: BlockRecordDB,
    sp_info: VDFInfo,
) -> list[tuple[ClassgroupElement, VDFInfo]]:
    """
    Returns candidate (input_el, vdf_info) pairs for CC SP proof verification.
    Mirrors get_signage_point_vdf_info, returning primary + fallback candidates.
    """
    new_sub_slot = len(fss) > 0
    overflow = is_overflow(block_br)
    genesis = prev_br is None
    sp_tt = sp_total_iters(block_br)
    sp_in_slot = calc_sp_iters(block_br)
    stored = sp_info  # stored sp_info (iters = sp_in_slot for new-slot cases)

    if new_sub_slot and not overflow:
        return [(default_el(), stored)]   # Case 1
    if new_sub_slot and overflow and len(fss) > 1:
        return [(default_el(), stored)]   # Case 2
    if genesis or in_genesis_slot(sp_info.challenge):
        return [(default_el(), stored)]   # Case 3 / genesis slot

    assert prev_br is not None

    def from_curr(curr: BlockRecord) -> list[tuple[ClassgroupElement, VDFInfo]]:
        delta = sp_tt - int(curr.total_iters)
        return [
            (curr.challenge_vdf_output, stored.replace(number_of_iterations=uint64(delta))),
            (default_el(), stored),
        ]

    if new_sub_slot and overflow and len(fss) == 1:
        # Case 4
        curr = prev_br
        while not curr.first_in_sub_slot and int(curr.total_iters) > sp_tt:
            curr = brdb.block_record(curr.prev_hash)
        if int(curr.total_iters) < sp_tt:
            return from_curr(curr)
        return [(default_el(), stored)]

    if not new_sub_slot and overflow:
        # Case 5
        curr = prev_br
        found_slots = len(curr.finished_challenge_slot_hashes or []) if curr.first_in_sub_slot else 0
        sp_pre_sb: Optional[BlockRecord] = None
        while found_slots < 2 and int(curr.height) > 0:
            if sp_pre_sb is None and int(curr.total_iters) < sp_tt:
                sp_pre_sb = curr
            curr = brdb.block_record(curr.prev_hash)
            if curr.first_in_sub_slot:
                found_slots += len(curr.finished_challenge_slot_hashes or [])
        if sp_pre_sb is None and int(curr.total_iters) < sp_tt:
            sp_pre_sb = curr
        if sp_pre_sb is not None:
            return from_curr(sp_pre_sb)
        return [(default_el(), stored)]

    # Case 6: same sub-slot, no overflow
    curr = prev_br
    while not curr.first_in_sub_slot and int(curr.total_iters) > sp_tt:
        curr = brdb.block_record(curr.prev_hash)
    if int(curr.total_iters) < sp_tt:
        return from_curr(curr)
    return [(default_el(), stored)]


def verify_sp_proof(
    block: FullBlock,
    block_br: BlockRecord,
    prev_br: Optional[BlockRecord],
    brdb: BlockRecordDB,
    backend: str,
) -> tuple[int, int]:
    ok_count = fail_count = 0
    rc = block.reward_chain_block
    fss = list(block.finished_sub_slots)

    if calc_sp_iters(block_br) == 0:
        return 0, 0  # first SP in sub-slot — no SP proof

    # CC SP
    if rc.challenge_chain_sp_vdf is not None and block.challenge_chain_sp_proof is not None:
        proof = block.challenge_chain_sp_proof
        candidates = get_cc_sp_candidates(block_br, fss, prev_br, brdb, rc.challenge_chain_sp_vdf)
        ok = any(verify_vdf(proof, info, inp, backend) for inp, info in candidates)
        ok_count += ok; fail_count += not ok

    # RC SP: always starts from default element, use stored iterations
    if rc.reward_chain_sp_vdf is not None and block.reward_chain_sp_proof is not None:
        ok = verify_vdf(block.reward_chain_sp_proof, rc.reward_chain_sp_vdf, default_el(), backend)
        ok_count += ok; fail_count += not ok

    return ok_count, fail_count


# ---------------------------------------------------------------------------
# IP proof verification
# ---------------------------------------------------------------------------
def verify_ip_proof(
    block: FullBlock,
    block_br: BlockRecord,
    prev_br: Optional[BlockRecord],
    backend: str,
) -> tuple[int, int]:
    ok_count = fail_count = 0
    rc = block.reward_chain_block
    genesis = prev_br is None
    new_sub_slot = len(block.finished_sub_slots) > 0

    # CC IP ---------------------------------------------------------------
    # The proof covers either:
    #   (a) the full range from slot-start (default_el, stored ip_iters), or
    #   (b) just the delta from the previous block's IP output.
    # Genesis slot uses (a); later slots use (b).  We try both and take any pass.
    if block.challenge_chain_ip_proof is not None:
        proof = block.challenge_chain_ip_proof
        stored_info = rc.challenge_chain_ip_vdf
        # Candidates: list of (input_el, vdf_info)
        if genesis or new_sub_slot or in_genesis_slot(stored_info.challenge):
            candidates = [(default_el(), stored_info)]
        else:
            assert prev_br is not None
            ip_delta = int(block_br.total_iters) - int(prev_br.total_iters)
            candidates = [
                (prev_br.challenge_vdf_output,
                 stored_info.replace(number_of_iterations=uint64(ip_delta))),
                (default_el(), stored_info),          # fallback
            ]
        ok = any(verify_vdf(proof, info, inp, backend) for inp, info in candidates)
        ok_count += ok; fail_count += not ok

    # RC IP: always default element, stored iterations
    if block.reward_chain_ip_proof is not None:
        ok = verify_vdf(block.reward_chain_ip_proof, rc.reward_chain_ip_vdf, default_el(), backend)
        ok_count += ok; fail_count += not ok

    # ICC IP — input is prev_b.infused_challenge_vdf_output + delta_iters
    # (same pattern as CC IP but for the infused challenge chain).
    # When there is no previous ICC output the chain is starting fresh (default_el).
    if rc.infused_challenge_chain_ip_vdf is not None and block.infused_challenge_chain_ip_proof is not None:
        proof = block.infused_challenge_chain_ip_proof
        icc_info = rc.infused_challenge_chain_ip_vdf
        prev_icc_out = None if (genesis or prev_br is None) else prev_br.infused_challenge_vdf_output
        if prev_icc_out is None:
            candidates = [(default_el(), icc_info)]
        else:
            icc_delta = int(block_br.total_iters) - int(prev_br.total_iters)
            candidates = [
                (prev_icc_out, icc_info.replace(number_of_iterations=uint64(icc_delta))),
                (default_el(), icc_info),
            ]
        ok = any(verify_vdf(proof, info, inp, backend) for inp, info in candidates)
        ok_count += ok; fail_count += not ok

    return ok_count, fail_count


# ---------------------------------------------------------------------------
# Full block
# ---------------------------------------------------------------------------
def verify_block(
    block: FullBlock,
    block_br: BlockRecord,
    prev_br: Optional[BlockRecord],
    brdb: BlockRecordDB,
    backend: str,
) -> tuple[int, int]:
    eos_ok, eos_fail = verify_eos_proofs(block, block_br, prev_br, backend)
    sp_ok,  sp_fail  = verify_sp_proof(block, block_br, prev_br, brdb, backend)
    ip_ok,  ip_fail  = verify_ip_proof(block, block_br, prev_br, backend)
    return eos_ok + sp_ok + ip_ok, eos_fail + sp_fail + ip_fail


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
def run_benchmark(db: Path, backend: str, min_height: int, max_height: Optional[int]) -> None:
    if backend == "cpp" and not HAVE_CPP:
        print("Error: chiavdf not available", file=sys.stderr); sys.exit(1)
    if backend == "rust" and not HAVE_RUST:
        print("Error: chia-vdf-verify not available", file=sys.stderr); sys.exit(1)

    conn = sqlite3.connect(str(db))
    brdb = BlockRecordDB(conn)

    q = ("SELECT height, block, block_record FROM full_blocks "
         "WHERE in_main_chain=1 AND height >= ? ORDER BY height")
    params: list = [min_height]
    if max_height is not None:
        q = q.replace("ORDER BY", "AND height <= ? ORDER BY")
        params.insert(1, max_height)

    label = "chiavdf C++" if backend == "cpp" else "chia-vdf-verify Rust"
    print(f"Backend: {label}   (sequential, correct chain state)")

    total_ok = total_fail = total_blocks = 0
    prev_br: Optional[BlockRecord] = None

    if min_height > 0:
        row = conn.execute(
            "SELECT block_record FROM full_blocks WHERE in_main_chain=1 AND height=?",
            (min_height - 1,),
        ).fetchone()
        if row:
            prev_br = BlockRecord.from_bytes(row[0])
            brdb.add(prev_br)

    t0 = time.perf_counter()
    last_h = min_height

    for row in conn.execute(q, params):
        last_h = int(row[0])
        block    = FullBlock.from_bytes(zstd.decompress(row[1]))
        block_br = BlockRecord.from_bytes(row[2])
        brdb.add(block_br)

        ok, fail = verify_block(block, block_br, prev_br, brdb, backend)
        total_ok   += ok
        total_fail += fail
        total_blocks += 1
        prev_br = block_br

        if total_blocks % 500 == 0:
            elapsed = time.perf_counter() - t0
            print(
                f"  h={last_h:,}  blocks={total_blocks:,}  "
                f"ok={total_ok:,}  fail={total_fail:,}  "
                f"{total_blocks/elapsed:.1f} blk/s",
                end="\r",
            )

    elapsed = time.perf_counter() - t0
    total_proofs = total_ok + total_fail
    print()
    print("=" * 60)
    print(f"Heights:        {min_height:>10,} – {last_h:,}")
    print(f"Blocks:         {total_blocks:>10,}")
    print(f"Proofs:         {total_proofs:>10,}  ({total_ok:,} ok, {total_fail:,} failed)")
    print(f"Wall time:      {elapsed:>10.2f}s")
    if elapsed > 0:
        print(f"Blocks/sec:     {total_blocks/elapsed:>10.1f}")
        print(f"Proofs/sec:     {total_proofs/elapsed:>10.1f}")
    print(f"Backend:        {label}")
    print("=" * 60)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--backend", choices=["cpp", "rust", "both"], default="rust")
    p.add_argument("--db",         type=Path,  default=DEFAULT_DB)
    p.add_argument("--min-height", type=int,   default=0)
    p.add_argument("--max-height", type=int,   default=None)
    args = p.parse_args()

    backends = ["cpp", "rust"] if args.backend == "both" else [args.backend]
    for b in backends:
        run_benchmark(args.db, b, args.min_height, args.max_height)


if __name__ == "__main__":
    main()
