#!/usr/bin/env python3
"""
Read packing_ver3_batch_summary.tsv (from batch_pack_ver3_all_pdbs.py) and write
packing_scores_ver3.tsv in the 6QZH project root (parent of `scripts/`, default).

Columns: basename, pdb_path, finalscore, n_events, status, error
Rows sorted by finalscore descending (missing scores last).
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    here = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--summary",
        type=Path,
        default=here / "packing_ver3_batch_summary.tsv",
        help="Input TSV from batch pack",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=here / "packing_scores_ver3.tsv",
        help="Output score file",
    )
    args = ap.parse_args()

    rows = []
    with open(args.summary, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            p = (row.get("pdb_path") or "").strip()
            fs = (row.get("finalscore") or "").strip()
            ne = (row.get("n_events") or "").strip()
            err = (row.get("error") or "").strip()
            base = Path(p).name if p else ""
            rows.append(
                {
                    "basename": base,
                    "pdb_path": p,
                    "finalscore": fs,
                    "n_events": ne,
                    "status": "error" if err else "ok",
                    "error": err,
                }
            )

    def sort_key(row: dict) -> tuple:
        fs = row["finalscore"]
        if not fs:
            return (1, 1e9, row["basename"])
        try:
            return (0, -float(fs), row["basename"])
        except ValueError:
            return (1, 1e9, row["basename"])

    rows.sort(key=sort_key)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["basename", "pdb_path", "finalscore", "n_events", "status", "error"],
            delimiter="\t",
        )
        w.writeheader()
        for row in rows:
            w.writerow(row)

    print(f"Wrote {args.out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
