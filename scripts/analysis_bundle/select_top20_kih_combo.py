#!/usr/bin/env python
"""
Select Top20 designs by a composite score from:
  - sc_value (higher better)
  - knob_score (higher better, KIH)
  - i_ptm (higher better)
  - rmsd (lower better)
  - holes_chA_chB (lower better)

Then copy corresponding Rosetta relaxed PDBs into a target folder.
"""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


def minmax_norm(series: pd.Series, higher_is_better: bool) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    mn = np.nanmin(s.values)
    mx = np.nanmax(s.values)
    if not np.isfinite(mn) or not np.isfinite(mx) or mx == mn:
        return pd.Series(np.full(len(s), 0.5), index=s.index)
    norm = (s - mn) / (mx - mn)
    if higher_is_better:
        return norm
    return 1.0 - norm


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kih_csv", required=True, help="KIH packing results CSV (knob_score)")
    ap.add_argument("--merged_csv", required=True, help="Merged MPNN+Rosetta scores CSV (sc_value, rmsd, i_ptm, holes)")
    ap.add_argument("--rosetta_pdb_dir", required=True, help="Directory with Rosetta relaxed PDBs")
    ap.add_argument("--out_dir", required=True, help="Output directory for selected Top20")
    ap.add_argument(
        "--weights",
        default="1,1,1,1,1",
        help="Weights for: sc_value, knob_score, i_ptm, rmsd, holes_chA_chB (comma separated)",
    )
    ap.add_argument("--top_n", type=int, default=20)
    args = ap.parse_args()

    w_sc, w_knob, w_ptm, w_rmsd, w_holes = [float(x) for x in args.weights.split(",")]

    kih = pd.read_csv(args.kih_csv)
    merged = pd.read_csv(args.merged_csv)

    # KIH CSV uses pdb_file like "..._0001.pdb"
    kih["description"] = kih["pdb_file"].astype(str).map(lambda x: os.path.splitext(os.path.basename(x))[0])

    # merged_scores has "description" column already
    need_cols_merged = {"description", "sc_value", "i_ptm", "rmsd", "holes_chA_chB"}
    missing = need_cols_merged - set(merged.columns)
    if missing:
        raise SystemExit(f"merged_csv missing columns: {sorted(missing)}")

    need_cols_kih = {"description", "knob_score"}
    missing_k = need_cols_kih - set(kih.columns)
    if missing_k:
        raise SystemExit(f"kih_csv missing columns: {sorted(missing_k)}")

    df = merged.merge(kih[["description", "knob_score"]], on="description", how="inner")
    if df.empty:
        raise SystemExit("Merge resulted in 0 rows. Check description/pdb filename alignment.")

    # Composite features
    df["sc_norm"] = minmax_norm(df["sc_value"], higher_is_better=True)
    df["knob_norm"] = minmax_norm(df["knob_score"], higher_is_better=True)
    df["ptm_norm"] = minmax_norm(df["i_ptm"], higher_is_better=True)
    df["rmsd_norm"] = minmax_norm(df["rmsd"], higher_is_better=False)
    df["holes_norm"] = minmax_norm(df["holes_chA_chB"], higher_is_better=False)

    denom = (w_sc + w_knob + w_ptm + w_rmsd + w_holes) or 1.0
    df["combo_score"] = (
        w_sc * df["sc_norm"]
        + w_knob * df["knob_norm"]
        + w_ptm * df["ptm_norm"]
        + w_rmsd * df["rmsd_norm"]
        + w_holes * df["holes_norm"]
    ) / denom

    top = df.sort_values("combo_score", ascending=False).head(args.top_n).copy()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_csv = out_dir / f"top{args.top_n}_kih_combo.csv"
    keep_cols = [
        "description",
        "combo_score",
        "sc_value",
        "knob_score",
        "i_ptm",
        "rmsd",
        "holes_chA_chB",
        "dG_separated",
        "packstat",
    ]
    keep_cols = [c for c in keep_cols if c in top.columns]
    top[keep_cols].to_csv(out_csv, index=False)

    # Copy PDBs
    missing_pdb = []
    copied = 0
    for desc in top["description"].astype(str).tolist():
        src = Path(args.rosetta_pdb_dir) / f"{desc}.pdb"
        if not src.exists():
            missing_pdb.append(str(src))
            continue
        dst = out_dir / f"{desc}.pdb"
        shutil.copy2(src, dst)
        copied += 1

    print(f"Selected: {len(top)} rows")
    print(f"Copied PDBs: {copied} into {out_dir}")
    if missing_pdb:
        print(f"WARNING: missing {len(missing_pdb)} PDBs. Example: {missing_pdb[0]}")
    print(f"Wrote: {out_csv}")


if __name__ == "__main__":
    main()

