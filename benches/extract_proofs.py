#!/usr/bin/env python3
"""Extract end-of-slot VDF proofs from a Chia blockchain DB.

Usage:
    python benches/extract_proofs.py [--db PATH] [--count N] [--output PATH]

Requires: chia_rs, zstd (pip install chia-rs zstd)
"""
from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
from pathlib import Path

import zstd
import chia_rs


def extract(db_path: str, n: int) -> list[dict]:
    conn = sqlite3.connect(db_path)
    cases: list[dict] = []
    rows = conn.execute(
        "SELECT block FROM full_blocks ORDER BY height DESC LIMIT ?",
        (n * 100,),
    ).fetchall()
    for (raw,) in rows:
        block = chia_rs.FullBlock.from_bytes(zstd.decompress(raw))
        for ss in (block.finished_sub_slots or []):
            vdf = ss.challenge_chain.challenge_chain_end_of_slot_vdf
            proof = ss.proofs.challenge_chain_slot_proof
            blob = vdf.output.data + bytes(proof.witness)
            cases.append({
                "challenge": bytes(vdf.challenge).hex(),
                "proof_blob": blob.hex(),
                "iters": vdf.number_of_iterations,
                "witness_type": proof.witness_type,
                "height": block.height,
            })
        if len(cases) >= n:
            break
    conn.close()
    return cases[:n]


def main() -> None:
    default_db = os.path.expanduser("~/.chia/mainnet/db/blockchain_v2_mainnet.sqlite")
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--db", default=default_db, help="blockchain SQLite DB path")
    parser.add_argument("--count", type=int, default=20, help="number of proofs to extract")
    parser.add_argument("--output", default="benches/proofs.json", help="output JSON path")
    args = parser.parse_args()

    if not Path(args.db).exists():
        sys.exit(f"DB not found: {args.db}")

    cases = extract(args.db, args.count)
    if not cases:
        sys.exit("No EOS VDF proofs found")

    with open(args.output, "w") as f:
        json.dump(cases, f, indent=2)

    heights = sorted(c["height"] for c in cases)
    wtypes = sorted(set(c["witness_type"] for c in cases))
    print(f"Extracted {len(cases)} proofs (heights {heights[0]:,}–{heights[-1]:,}, "
          f"witness types {wtypes})")
    print(f"Wrote {args.output} ({os.path.getsize(args.output):,} bytes)")


if __name__ == "__main__":
    main()
