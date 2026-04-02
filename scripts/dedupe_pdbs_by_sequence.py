#!/usr/bin/env python3
"""
Deduplicate PDB files under a directory by **amino-acid sequence** (all chains, model 0).

- Builds a canonical key: sorted by chain id, format ``A:SEQ|B:SEQ|...`` (one-letter codes).
- Default: write a report only (dry-run). Use ``--move`` to move duplicate files aside.

Requires: biopython

Example:
  cd .../6QZH
  python3 scripts/dedupe_pdbs_by_sequence.py
  python3 scripts/dedupe_pdbs_by_sequence.py --move
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
from collections import defaultdict
from pathlib import Path

from Bio.PDB import PDBParser, PPBuilder


def first_model(structure):
    """First MODEL in file; do not use structure[0] — that requires model id == 0."""
    return next(structure.get_models())


def chain_sequences(structure) -> list[tuple[str, str]]:
    """Return [(chain_id, one_letter_seq), ...] sorted by chain id."""
    model = first_model(structure)
    ppb = PPBuilder()
    out: list[tuple[str, str]] = []
    chains = sorted(model.get_chains(), key=lambda c: str(c.id))
    for chain in chains:
        parts: list[str] = []
        for pp in ppb.build_peptides(chain):
            parts.append(str(pp.get_sequence()))
        seq = "".join(parts)
        out.append((str(chain.id), seq))
    return out


def sequence_key(structure) -> str:
    parts = chain_sequences(structure)
    return "|".join(f"{cid}:{seq}" for cid, seq in parts if seq)


def sequence_fingerprint(key: str) -> str:
    return hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]


def main() -> None:
    here = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pdbs-root",
        type=Path,
        default=here / "pdbs",
        help="Root directory to scan for *.pdb (recursive)",
    )
    ap.add_argument(
        "--out-report",
        type=Path,
        default=here / "pdbs_sequence_dedup_report.tsv",
        help="TSV report path",
    )
    ap.add_argument(
        "--dup-dir",
        type=Path,
        default=here / "pdbs_duplicates_by_sequence",
        help="Where to move duplicate PDBs when --move is set",
    )
    ap.add_argument(
        "--move",
        action="store_true",
        help="Move duplicate files (not the kept representative) to --dup-dir",
    )
    ap.add_argument(
        "--keep",
        choices=("first", "last"),
        default="first",
        help="Which file to keep per sequence group (sorted by rel path)",
    )
    args = ap.parse_args()

    root = args.pdbs_root.resolve()
    parser = PDBParser(QUIET=True)

    pdb_paths = sorted(root.rglob("*.pdb"))
    if not pdb_paths:
        print(f"No .pdb under {root}", flush=True)
        return

    # seq_hash -> list of (relpath, full_key)
    groups: dict[str, list[tuple[str, str]]] = defaultdict(list)
    errors: list[tuple[str, str]] = []

    for p in pdb_paths:
        rel = str(p.relative_to(root))
        try:
            st = parser.get_structure(p.stem, str(p))
            key = sequence_key(st)
            if not key or key.replace("|", "").replace(":", "") == "":
                errors.append((rel, "empty_sequence"))
                continue
            fp = sequence_fingerprint(key)
            groups[fp].append((rel, key))
        except Exception as e:
            errors.append((rel, repr(e)))

    # pick representative per group
    keep_rows: list[tuple[str, str, str, str]] = []
    dup_rows: list[tuple[str, str, str, str]] = []
    # columns: role, seq_fp, pdb_rel, sequence_key_preview

    for fp, items in sorted(groups.items(), key=lambda x: x[0]):
        items_sorted = sorted(items, key=lambda t: t[0])
        if args.keep == "last":
            items_sorted = list(reversed(items_sorted))
        rep_rel, rep_key = items_sorted[0]
        keep_rows.append(("keep", fp, rep_rel, rep_key[:200] + ("..." if len(rep_key) > 200 else "")))
        for rel, key in items_sorted[1:]:
            dup_rows.append(("duplicate", fp, rel, key[:200] + ("..." if len(key) > 200 else "")))

    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_report, "w", encoding="utf-8") as f:
        f.write("role\tseq_fingerprint\tpdb_path\tsequence_key_preview\n")
        for row in keep_rows + dup_rows:
            f.write("\t".join(row) + "\n")
        for rel, err in errors:
            f.write(f"error\t\t{rel}\t{err}\n")

    n_groups = len(groups)
    n_dup = len(dup_rows)
    print(f"Scanned {len(pdb_paths)} PDBs, {n_groups} unique sequences, {n_dup} duplicates, {len(errors)} errors.")
    print(f"Report: {args.out_report}")

    if args.move:
        if not dup_rows:
            print("Nothing to move.")
            return
        args.dup_dir.mkdir(parents=True, exist_ok=True)
        moved = 0
        for role, fp, rel, _ in dup_rows:
            if role != "duplicate":
                continue
            src = root / rel
            dst = args.dup_dir / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dst))
            moved += 1
        print(f"Moved {moved} duplicate files under {args.dup_dir}")
        print("Representatives remain under:", root)


if __name__ == "__main__":
    main()
