#!/usr/bin/env python
"""
Plot histograms for key KIH knob metrics.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kih_csv", required=True, help="outputs_4/kih_packing_results_opt.csv")
    ap.add_argument("--out_png", required=True, help="Output histogram png")
    args = ap.parse_args()

    kih = pd.read_csv(args.kih_csv)

    # Key knob-related metrics (all exist in current KIH csv)
    cols = [
        ("knob_score", "Knob Score (higher better)"),
        ("total_knobs", "Total Knobs (higher better)"),
        ("knob_density", "Knob Density (higher better)"),
        ("avg_knob_contacts", "Avg Contacts per Knob (higher better)"),
        ("hydrophobic_contacts", "Hydrophobic Contacts (higher better)"),
        ("hydrophobic_fraction", "Hydrophobic Fraction (higher better)"),
    ]

    existing = [(c, label) for c, label in cols if c in kih.columns]
    if len(existing) < 3:
        raise SystemExit(f"Not enough columns found in kih_csv. Found: {list(kih.columns)}")

    n = len(existing)
    nrows = 2 if n <= 4 else 3
    ncols = int(np.ceil(n / nrows))

    sns.set_style("whitegrid")
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 3.8 * nrows))
    axes = np.array(axes).reshape(-1)

    for i, (col, label) in enumerate(existing):
        ax = axes[i]
        vals = pd.to_numeric(kih[col], errors="coerce").dropna()
        if len(vals) == 0:
            continue
        ax.hist(vals, bins=30, edgecolor="black", alpha=0.75)
        ax.axvline(vals.median(), color="r", linestyle="--", linewidth=2, label=f"Median: {vals.median():.2f}")
        ax.set_xlabel(label)
        ax.set_ylabel("Count")
        ax.legend()

    # Turn off unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.suptitle("KIH Knob Metrics Histograms", fontsize=16)
    plt.tight_layout(rect=[0, 0.02, 1, 0.98])

    out = Path(args.out_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)

    print(f"Wrote: {out}")


if __name__ == "__main__":
    main()

