#!/usr/bin/env python3
"""
Merge InterfaceAnalyzer / Rosetta score lines at the end of each PDB with
packing_scores_ver3.tsv (pack_analysis_muOR_ver3).

Default: --root = 6QZH/pdbs, --packing / --out = files in 6QZH/ (project root, parent of `scripts/`).

Rosetta keys parsed (last occurrence in file wins):
  dG_separated, sc_value, holes-ABX_interface, packstat, dSASA_int, delta_unsatHbonds
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

ROSETTA_KEYS = {
    "dG_separated": "dG_separated",
    "sc_value": "sc_value",
    "holes-ABX_interface": "Holes",
    "packstat": "pack_stat",
    "dSASA_int": "dSASA_int",
    "delta_unsatHbonds": "delta_unsatHbonds",
}

_LINE_RE = re.compile(
    r"^(dG_separated|sc_value|holes-ABX_interface|packstat|dSASA_int|delta_unsatHbonds)\s+(.+?)\s*$"
)


def parse_rosetta_tail(pdb_path: Path) -> tuple[dict[str, str], str]:
    out: dict[str, str] = {}
    note = ""
    try:
        text = pdb_path.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        return {}, f"read_error:{e}"

    for line in text.splitlines():
        m = _LINE_RE.match(line.strip())
        if not m:
            continue
        key, val = m.group(1), m.group(2).strip()
        col = ROSETTA_KEYS[key]
        out[col] = val

    if not out:
        note = "no_rosetta_score_lines"
    return out, note


def main() -> None:
    here = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root",
        type=Path,
        default=here / "pdbs",
        help="Directory tree containing PDBs (paths in packing TSV are relative to this)",
    )
    ap.add_argument(
        "--packing",
        type=Path,
        default=here / "packing_scores_ver3.tsv",
        help="Input packing score table",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=here / "rosetta_packing_merged_scores.tsv",
        help="Output merged TSV",
    )
    args = ap.parse_args()
    root = args.root.resolve()

    packing_rows: list[dict] = []
    with open(args.packing, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        packing_rows = list(r)

    fieldnames = [
        "basename",
        "pdb_path",
        "dG_separated",
        "sc_value",
        "Holes",
        "pack_stat",
        "dSASA_int",
        "delta_unsatHbonds",
        "packing_finalscore",
        "packing_n_events",
        "packing_status",
        "packing_error",
        "rosetta_parse_note",
    ]

    merged = []
    for row in packing_rows:
        rel = (row.get("pdb_path") or "").strip()
        base = (row.get("basename") or "").strip()
        pfs = (row.get("finalscore") or "").strip()
        pne = (row.get("n_events") or "").strip()
        st = (row.get("status") or "").strip()
        err = (row.get("error") or "").strip()

        pdb_abs = root / rel if rel else root / base
        rosetta, note = parse_rosetta_tail(pdb_abs)

        merged.append(
            {
                "basename": base,
                "pdb_path": rel,
                "dG_separated": rosetta.get("dG_separated", ""),
                "sc_value": rosetta.get("sc_value", ""),
                "Holes": rosetta.get("Holes", ""),
                "pack_stat": rosetta.get("pack_stat", ""),
                "dSASA_int": rosetta.get("dSASA_int", ""),
                "delta_unsatHbonds": rosetta.get("delta_unsatHbonds", ""),
                "packing_finalscore": pfs,
                "packing_n_events": pne,
                "packing_status": st,
                "packing_error": err,
                "rosetta_parse_note": note,
            }
        )

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in merged:
            w.writerow(row)

    n_miss = sum(1 for m in merged if m["rosetta_parse_note"])
    print(f"Wrote {args.out} ({len(merged)} rows). Rosetta lines missing: {n_miss}")


if __name__ == "__main__":
    main()
