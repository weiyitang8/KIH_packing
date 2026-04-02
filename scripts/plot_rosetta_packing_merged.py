#!/usr/bin/env python3
"""
Plot Rosetta (InterfaceAnalyzer) metrics vs pack_analysis_muOR_ver3 scores from
`rosetta_packing_merged_scores.tsv`.

Outputs (default under 6QZH/figures_rosetta_packing/):
  - correlation_heatmap.png   Pearson r matrix for numeric columns
  - scatter_packing_vs_rosetta.png  KIH_score vs each Rosetta metric (2x3)
  - hist_packing_finalscore.png     distribution of KIH_score (TSV: packing_finalscore)
  - hist_rosetta_metrics.png        histograms of dG_separated, sc_value, Holes, … (2x3)

Requires: pandas, matplotlib, seaborn, numpy

Usage:
  cd .../6QZH
  python3 scripts/plot_rosetta_packing_merged.py
  python3 scripts/plot_rosetta_packing_merged.py --tsv rosetta_packing_merged_scores.tsv --out-dir figures_custom
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# TSV column vs label on figures
PACKING_SCORE_COL = "packing_finalscore"
PACKING_SCORE_LABEL = "KIH_score"

NUMERIC_COLS = [
    "dG_separated",
    "sc_value",
    "Holes",
    "pack_stat",
    "dSASA_int",
    "delta_unsatHbonds",
    PACKING_SCORE_COL,
    "packing_n_events",
]

# Pearson heatmap: omit packing_n_events (event count, not comparable as continuous r)
HEATMAP_COLS = [c for c in NUMERIC_COLS if c != "packing_n_events"]


def load_table(tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv, sep="\t")
    for c in NUMERIC_COLS:
        if c not in df.columns:
            continue
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def plot_heatmap(df: pd.DataFrame, out_png: Path) -> None:
    sub = df[[c for c in HEATMAP_COLS if c in df.columns]].copy()
    # Pairwise Pearson (pandas skips NaN per pair)
    corr = sub.corr(numeric_only=True, min_periods=30)
    if PACKING_SCORE_COL in corr.columns:
        corr = corr.rename(
            columns={PACKING_SCORE_COL: PACKING_SCORE_LABEL},
            index={PACKING_SCORE_COL: PACKING_SCORE_LABEL},
        )
    if corr.empty or len(corr) < 2:
        print("heatmap: not enough numeric overlap, skip")
        return
    # Smaller figure; matrix row/column names (e.g. Holes) larger than colorbar ticks
    fs_names = 18
    fs_cbar_ticks = 18
    fs_annot = 16
    fig, ax = plt.subplots(figsize=(9, 7.2))
    sns.heatmap(
        corr,
        annot=True,
        fmt=".2f",
        cmap="RdBu_r",
        center=0,
        square=True,
        ax=ax,
        vmin=-1,
        vmax=1,
        annot_kws={"size": fs_annot},
    )
    ax.tick_params(axis="both", labelsize=fs_names)
    for a in fig.axes:
        if a is not ax:
            a.tick_params(labelsize=fs_cbar_ticks)
    plt.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_png}")


def plot_scatter_grid(df: pd.DataFrame, out_png: Path) -> None:
    rosetta_x = [
        "dG_separated",
        "sc_value",
        "Holes",
        "pack_stat",
        "dSASA_int",
        "delta_unsatHbonds",
    ]
    y = PACKING_SCORE_COL
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    axes = axes.ravel()
    for i, xc in enumerate(rosetta_x):
        ax = axes[i]
        pair = df[[xc, y]].dropna()
        if pair.empty:
            ax.set_visible(False)
            continue
        ax.scatter(pair[xc], pair[y], alpha=0.25, s=12, edgecolors="none")
        ax.set_xlabel(xc)
        ax.set_ylabel(PACKING_SCORE_LABEL)
        if len(pair) >= 3:
            z = np.polyfit(pair[xc].values, pair[y].values, 1)
            xs = np.linspace(pair[xc].min(), pair[xc].max(), 50)
            ax.plot(xs, np.poly1d(z)(xs), "r-", lw=1, alpha=0.8)
    plt.suptitle(f"{PACKING_SCORE_LABEL} vs Rosetta / InterfaceAnalyzer metrics")
    plt.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_png}")


def plot_hist_packing(df: pd.DataFrame, out_png: Path) -> None:
    s = df[PACKING_SCORE_COL].dropna()
    if s.empty:
        print(f"hist: no {PACKING_SCORE_COL}, skip")
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(s, bins=40, color="steelblue", edgecolor="white", alpha=0.9)
    ax.set_xlabel(f"{PACKING_SCORE_LABEL} (ver3)")
    ax.set_ylabel("count")
    ax.set_title(f"Distribution of {PACKING_SCORE_LABEL}")
    plt.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_png}")


def plot_hist_rosetta_metrics(df: pd.DataFrame, out_png: Path) -> None:
    """2x3 histograms for InterfaceAnalyzer / Rosetta columns (excluding packing columns)."""
    rosetta_cols = [
        "dG_separated",
        "sc_value",
        "Holes",
        "pack_stat",
        "dSASA_int",
        "delta_unsatHbonds",
    ]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    axes = axes.ravel()
    for i, col in enumerate(rosetta_cols):
        ax = axes[i]
        s = df[col].dropna() if col in df.columns else pd.Series(dtype=float)
        if s.empty:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(col)
            continue
        n_bins = min(40, max(10, int(np.sqrt(len(s)))))
        ax.hist(s, bins=n_bins, color="coral", edgecolor="white", alpha=0.9)
        ax.set_xlabel(col)
        ax.set_ylabel("count")
        ax.set_title(f"{col} (n={len(s)})")
    plt.suptitle("Rosetta / InterfaceAnalyzer score distributions")
    plt.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_png}")


def main() -> None:
    here = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--tsv",
        type=Path,
        default=here / "rosetta_packing_merged_scores.tsv",
        help="Merged TSV from merge_rosetta_and_packing_scores.py",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=here / "figures_rosetta_packing",
        help="Output directory for PNGs",
    )
    ap.add_argument(
        "--only-rosetta-ok",
        action="store_true",
        help="Drop rows where rosetta_parse_note is non-empty",
    )
    args = ap.parse_args()

    df = load_table(args.tsv)
    if args.only_rosetta_ok and "rosetta_parse_note" in df.columns:
        df = df[df["rosetta_parse_note"].fillna("") == ""].copy()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        try:
            plt.style.use("seaborn-whitegrid")
        except OSError:
            pass

    plot_heatmap(df, args.out_dir / "correlation_heatmap.png")
    plot_scatter_grid(df, args.out_dir / "scatter_packing_vs_rosetta.png")
    plot_hist_packing(df, args.out_dir / "hist_packing_finalscore.png")
    plot_hist_rosetta_metrics(df, args.out_dir / "hist_rosetta_metrics.png")

    print(f"Done. Figures in: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
