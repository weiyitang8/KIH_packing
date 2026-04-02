#!/usr/bin/env python3
"""
Build a **sequence logo** from an ensemble of PDB files (same chain, PDB numbering).

- For each structure, reads amino acids on the chosen chain (default ``B``) in model 0,
  keyed by residue number (resSeq + insertion code).
- Tallies counts per position and plots with **logomaker** (Schneider–Stephens style).

Requires: biopython, pandas, matplotlib, logomaker

Examples::

  cd .../thesis/figuremaking/6QZH
  pip install logomaker
  python3 plot_sequence_logo.py
  python3 plot_sequence_logo.py --dedupe-by-sequence
  python3 plot_sequence_logo.py --region 100 250 --out figures_rosetta_packing/sequence_logo_TM.png
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio.Data.PDBData import protein_letters_3to1, protein_letters_3to1_extended
from Bio.PDB import PDBParser, PPBuilder

# Optional: skip redundant structures (same AA sequence as another file)
def _first_model(structure):
    return next(structure.get_models())


def _chain_sequences(structure) -> list[tuple[str, str]]:
    model = _first_model(structure)
    ppb = PPBuilder()
    out: list[tuple[str, str]] = []
    for chain in sorted(model.get_chains(), key=lambda c: str(c.id)):
        parts: list[str] = []
        for pp in ppb.build_peptides(chain):
            parts.append(str(pp.get_sequence()))
        seq = "".join(parts)
        out.append((str(chain.id), seq))
    return out


def _sequence_fingerprint_key(structure) -> str:
    parts = _chain_sequences(structure)
    return "|".join(f"{cid}:{seq}" for cid, seq in parts if seq)


def _fp_hash(key: str) -> str:
    return hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]


def residue_to_one(res) -> str | None:
    """Map PDB residue to one-letter code (20 standard + common HET like MSE→M)."""
    name = res.get_resname().strip().upper()
    if name in protein_letters_3to1:
        return protein_letters_3to1[name]
    if name in protein_letters_3to1_extended:
        return protein_letters_3to1_extended[name]
    return None


def get_chain(structure, chain_id: str):
    model = _first_model(structure)
    cid = chain_id.strip()
    for ch in model.get_chains():
        if str(ch.id).strip() == cid:
            return ch
    raise KeyError(f"chain {chain_id!r} not found")


def chain_position_counts(chain) -> dict[tuple[int, str], str]:
    """Map (resseq, icode) -> one-letter aa for standard protein residues."""
    out: dict[tuple[int, str], str] = {}
    for res in chain:
        if res.id[0] != " ":
            continue
        aa = residue_to_one(res)
        if aa is None:
            continue
        rseq = res.id[1]
        icode = res.id[2]
        if not isinstance(icode, str):
            icode = str(icode)
        out[(int(rseq), icode)] = aa
    return out


def pos_label(key: tuple[int, str]) -> str:
    rseq, icode = key
    s = str(rseq)
    if icode and icode != " ":
        s += icode.strip()
    return s


def build_count_matrix(
    pdb_paths: list[Path],
    chain_id: str,
    region: tuple[int | None, int | None],
) -> tuple[pd.DataFrame, int, int]:
    """
    Rows = positions (sorted), columns = 20 AAs (A,C,...,Y,W order used by logomaker).
    Returns (count_df, n_used_pdbs, n_skipped).
    """
    letters = list("ACDEFGHIKLMNPQRSTVWY")
    pos_totals: dict[tuple[int, str], Counter[str]] = defaultdict(Counter)
    used = 0
    skipped = 0
    parser = PDBParser(QUIET=True)

    lo, hi = region
    for p in pdb_paths:
        try:
            st = parser.get_structure(p.stem, str(p))
            ch = get_chain(st, chain_id)
            pmap = chain_position_counts(ch)
        except Exception:
            skipped += 1
            continue
        if not pmap:
            skipped += 1
            continue
        for pos, aa in pmap.items():
            rseq = pos[0]
            if lo is not None and rseq < lo:
                continue
            if hi is not None and rseq > hi:
                continue
            if aa in letters:
                pos_totals[pos][aa] += 1
        used += 1

    if not pos_totals:
        raise SystemExit("No residues tallied; check --chain, --region, or PDB contents.")

    sorted_pos = sorted(pos_totals.keys(), key=lambda x: (x[0], x[1]))
    rows = []
    index_labels = []
    for pos in sorted_pos:
        c = pos_totals[pos]
        total = sum(c.values())
        if total == 0:
            continue
        row = {a: c.get(a, 0) for a in letters}
        rows.append(row)
        index_labels.append(pos_label(pos))

    df = pd.DataFrame(rows, index=index_labels, columns=letters)
    return df, used, skipped


def plot_logo(
    count_df: pd.DataFrame,
    out_png: Path,
    title: str,
    xtick_every: int,
) -> None:
    import logomaker

    # logomaker requires integer row index 0..N-1
    res_labels = list(count_df.index)
    prob_df = count_df.div(count_df.sum(axis=1).replace(0, 1), axis=0).copy()
    prob_df.index = range(len(prob_df))

    n = len(prob_df)
    width_in = max(8.0, min(48.0, 0.12 * n))
    fig, ax = plt.subplots(figsize=(width_in, 3.5))
    logomaker.Logo(
        prob_df,
        ax=ax,
        color_scheme="chemistry",
        font_name="DejaVu Sans",  # Arial often missing on Linux
        show_spines=False,
    )
    ax.set_ylabel("probability")
    ax.set_xlabel("residue (PDB numbering)")
    ax.set_title(title)
    if n > 80 and xtick_every > 1:
        ticks = list(range(0, n, xtick_every))
        ax.set_xticks(ticks)
        ax.set_xticklabels([res_labels[i] for i in ticks], rotation=65, ha="right", fontsize=6)
    else:
        ax.set_xticks(range(n))
        ax.set_xticklabels(res_labels, rotation=65, ha="right", fontsize=7)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    here = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description="Sequence logo from PDB ensemble")
    ap.add_argument("--pdbs-root", type=Path, default=here / "pdbs", help="Root for *.pdb (recursive)")
    ap.add_argument("--chain", type=str, default="B", help="Chain ID (default B)")
    ap.add_argument(
        "--dedupe-by-sequence",
        action="store_true",
        help="Keep one PDB per distinct full-structure sequence (fewer redundant stacks)",
    )
    ap.add_argument(
        "--region",
        type=int,
        nargs=2,
        metavar=("START", "END"),
        help="Inclusive residue number range (PDB resSeq); omit for full length",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=here / "figures_rosetta_packing" / "sequence_logo.png",
        help="Output PNG path",
    )
    ap.add_argument(
        "--xtick-every",
        type=int,
        default=10,
        help="Label every Nth column when many positions (default 10)",
    )
    args = ap.parse_args()

    try:
        import logomaker  # noqa: F401
    except ImportError:
        print("Install logomaker:  pip install logomaker", file=sys.stderr)
        sys.exit(1)

    root = args.pdbs_root.resolve()
    pdb_paths = sorted(root.rglob("*.pdb"))
    if not pdb_paths:
        sys.exit(f"No .pdb under {root}")

    if args.dedupe_by_sequence:
        parser = PDBParser(QUIET=True)
        seen: set[str] = set()
        unique_paths: list[Path] = []
        for p in pdb_paths:
            try:
                st = parser.get_structure(p.stem, str(p))
                key = _sequence_fingerprint_key(st)
                if not key.strip():
                    continue
                fp = _fp_hash(key)
                if fp in seen:
                    continue
                seen.add(fp)
                unique_paths.append(p)
            except Exception:
                continue
        print(
            f"Dedupe: {len(pdb_paths)} PDBs -> {len(unique_paths)} unique sequences",
            flush=True,
        )
        pdb_paths = unique_paths

    region = (None, None)
    if args.region:
        region = (args.region[0], args.region[1])

    count_df, used, skipped = build_count_matrix(pdb_paths, args.chain.strip(), region)
    title = (
        f"Sequence logo (chain {args.chain}, n={used} structures"
        + (f", res {args.region[0]}–{args.region[1]}" if args.region else "")
        + ")"
    )
    plot_logo(count_df, args.out, title, args.xtick_every if len(count_df) > 80 else 1)
    print(f"Used {used} PDBs, skipped {skipped} (parse/empty). Positions: {len(count_df)}")
    print(f"Wrote {args.out.resolve()}")


if __name__ == "__main__":
    main()
