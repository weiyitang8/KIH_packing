#!/usr/bin/env python
"""
Plot a correlation heatmap between Rosetta/AF2 metrics and KIH knob metrics.

Inputs:
  --merged_csv : outputs_4/merged_scores.csv
  --kih_csv     : outputs_4/kih_packing_results_opt.csv
Output:
  --out_png     : correlation heatmap with knobs included
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def extract_design_n_from_pdb_filename(pdb_file: str) -> tuple[int | None, int | None]:
    """
    Expected pdb_file like: design13_n38_design13_n38_0001.pdb
    Returns (design, n) from the first designX_nY occurrence.
    """
    base = Path(pdb_file).name
    m = re.search(r"design(\d+)_n(\d+)", base)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged_csv", required=True, help="outputs_4/merged_scores.csv")
    ap.add_argument("--kih_csv", required=True, help="outputs_4/kih_packing_results_opt.csv")
    ap.add_argument("--out_png", required=True, help="Output png path")
    args = ap.parse_args()

    merged = pd.read_csv(args.merged_csv)
    kih = pd.read_csv(args.kih_csv)

    # Drop extra unnamed index column if present
    merged = merged.loc[:, ~merged.columns.str.match(r"^Unnamed:")]

    needed_merged = {"design", "n", "sc_value", "i_ptm", "rmsd", "holes_chA_chB"}
    missing_merged = needed_merged - set(merged.columns)
    if missing_merged:
        raise SystemExit(f"merged_csv missing columns: {sorted(missing_merged)}")

    needed_kih = {"pdb_file", "knob_score", "total_knobs", "knob_density"}
    missing_kih = needed_kih - set(kih.columns)
    if missing_kih:
        raise SystemExit(f"kih_csv missing columns: {sorted(missing_kih)}")

    # Add design/n to kih
    design_n = kih["pdb_file"].map(extract_design_n_from_pdb_filename)
    kih[["design", "n"]] = pd.DataFrame(design_n.tolist(), index=kih.index)
    kih = kih.dropna(subset=["design", "n"]).copy()
    kih["design"] = kih["design"].astype(int)
    kih["n"] = kih["n"].astype(int)

    # Merge on design/n
    df = merged.merge(kih[["design", "n", "knob_score", "total_knobs", "knob_density"]],
                      on=["design", "n"], how="inner")
    if df.empty:
        raise SystemExit("After merging merged_csv and kih_csv, got 0 rows.")

    # Small correlation subset requested by user
    # (Note: user text may contain typos; use exact column names from CSVs)
    key_cols = [
        "mpnn",
        "plddt",
        "i_ptm",
        "rmsd",
        "dG_separated",
        "sc_value",
        "packstat",
        "dSASA_int",
        "holes_chA_chB",
        "delta_unsatHbonds",
        "knob_score",
    ]

    available_cols = [c for c in key_cols if c in df.columns]
    if len(available_cols) < 3:
        raise SystemExit("Not enough columns available to compute correlation.")

    corr = df[available_cols].corr(numeric_only=True)

    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(
        corr,
        annot=True,
        fmt=".2f",
        cmap="RdBu_r",
        center=0,
        square=True,
        ax=ax,
    )
    ax.set_title("Correlation (with KIH knobs)")
    plt.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

    print(f"Wrote: {out_png}")
    print(f"Rows used: {len(df)}, columns used: {len(available_cols)}")


if __name__ == "__main__":
    main()

