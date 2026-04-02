#!/usr/bin/env python
"""
Copy Rosetta relaxed PDBs whose KIH knob_score is above a threshold.

Inputs:
  --kih_csv        outputs_4/kih_packing_results_opt.csv
  --rosetta_pdb_dir outputs_4/rosetta_fastrelax
  --out_dir        outputs_4/knob_gt20
  --threshold      default 20 (knob_score > threshold)
"""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kih_csv", required=True)
    ap.add_argument("--rosetta_pdb_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--threshold", type=float, default=20.0)
    args = ap.parse_args()

    kih_csv = Path(args.kih_csv)
    pdb_dir = Path(args.rosetta_pdb_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    missing = 0
    rows = 0

    with open(kih_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        if "pdb_file" not in reader.fieldnames or "knob_score" not in reader.fieldnames:
            raise SystemExit("kih_csv must contain columns: pdb_file, knob_score")

        for row in reader:
            rows += 1
            pdb_file = row["pdb_file"]
            knob_score = float(row["knob_score"])
            if knob_score > args.threshold:
                src = pdb_dir / pdb_file
                if not src.exists():
                    missing += 1
                    continue
                shutil.copy2(src, out_dir / pdb_file)
                copied += 1

    print(f"Total rows in kih csv: {rows}")
    print(f"knob_score > {args.threshold}: {copied} copied, {missing} missing")


if __name__ == "__main__":
    main()

