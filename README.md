# 6QZH μOR ensemble — Rosetta interface metrics & three-helix packing scores

Thesis-side pipeline for **6QZH** designs: batch **pack_analysis_muOR_ver3** (three-helix knob–hole–style packing + SASA filters), merge with **Rosetta InterfaceAnalyzer** fields parsed from PDB footers, and produce publication-style figures.

**This directory is the Git repository root** (not the whole `thesis/` tree).

## What the packing score is

- **`packing_finalscore` in TSV** = `finalscore` from **`pack_analysis_muOR_ver3`** (`muOR/bin/pack_analysis_muOR_ver3.py`).
- It is **not** the two-chain script in `scripts/analysis_bundle/analyze_kih_packing.py`; ver3 explicitly models **three helices** as residue segments **TM5, TM6, TMB** and aggregates scored packing events.
- Figures label this quantity **`KIH_score`** for readability; the merged table column remains `packing_finalscore`.

## Repository layout

| Path | Role |
|------|------|
| **`scripts/`** | All Python entry points |
| `scripts/batch_pack_ver3_all_pdbs.py` | Batch-run ver3 on `pdbs/**/*.pdb` → `packing_ver3_batch_summary.tsv` |
| `scripts/generate_packing_scores_from_batch_summary.py` | Sortable `packing_scores_ver3.tsv` |
| `scripts/merge_rosetta_and_packing_scores.py` | PDB footer Rosetta + packing → `rosetta_packing_merged_scores.tsv` |
| `scripts/plot_rosetta_packing_merged.py` | Correlation heatmap, scatter grid, histograms |
| `scripts/plot_sequence_logo.py` | Sequence logo from ensemble (needs `pdbs/`) |
| `scripts/dedupe_pdbs_by_sequence.py` | Optional deduplication helper |
| `scripts/analysis_bundle/` | Extra KIH / knob analysis scripts |
| `figures_rosetta_packing/` | **Exported PNGs** (checked in for thesis/GitHub preview) |
| `README_SCORES.md` | TSV column cheat sheet & command snippets (Chinese) |

**`pdbs/`** (~1k+ structures) is **gitignored**. Place your ensemble locally or distribute via Zenodo/lab storage; the committed TSVs reproduce plots without PDBs.

## Dependencies

- **Batch packing**: Python 3 with BioPython, SciPy, etc., plus a **muOR checkout** above this repo (or `PYTHONPATH` to `muOR/bin`) so `bin/pack_analysis_muOR_ver3.py` is found — see `scripts/batch_pack_ver3_all_pdbs.py` header.
- **Plots**: `pip install -r requirements_plots.txt`

## Quick start (figures only)

```bash
cd /path/to/this/6QZH   # repository root
pip install -r requirements_plots.txt
python3 scripts/plot_rosetta_packing_merged.py
# → figures_rosetta_packing/correlation_heatmap.png, scatter_*, hist_*, …
```

Optional: `python3 scripts/plot_rosetta_packing_merged.py --only-rosetta-ok`

Full score rebuild: see **`README_SCORES.md`**.

## Figure outputs (`figures_rosetta_packing/`)

| File | Content |
|------|---------|
| `correlation_heatmap.png` | Pearson **r** between Rosetta metrics and **KIH_score** (no title; thesis-tuned fonts) |
| `scatter_packing_vs_rosetta.png` | **KIH_score** vs each Rosetta metric (2×3) |
| `hist_packing_finalscore.png` | **KIH_score** distribution |
| `hist_rosetta_metrics.png` | Rosetta metric histograms |
| `sequence_logo.png` | Logo from chain **B** (regenerate with `scripts/plot_sequence_logo.py` when `pdbs/` is present) |

## Data files in this folder (small, tracked)

- `rosetta_packing_merged_scores.tsv` — main merged table for plotting  
- `packing_scores_ver3.tsv`, `packing_ver3_batch_summary.tsv` — packing-only tables  
- `pdbs_sequence_dedup_report.tsv` — deduplication report (if used)

## License

If not specified elsewhere, follow your lab / institution policy for structures and code.
