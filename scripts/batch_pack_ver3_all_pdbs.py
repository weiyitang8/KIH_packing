#!/usr/bin/env python3
"""
Batch-run pack_analysis_muOR_ver3 on every *.pdb under this script's target directory
(default: same folder as this file's parent .../6QZH/pdbs).

Writes (default: **6QZH project root**, parent of `scripts/`):
  packing_ver3_batch_summary.tsv   (tab-separated, one row per PDB)
  packing_ver3_batch_errors.log    (tracebacks for failed files)

Usage:
  cd /path/to/6QZH
  python3 scripts/batch_pack_ver3_all_pdbs.py

Optional env:
  MUSOR_BATCH_JOBS=8   number of parallel workers (default: min(8, CPU count))
  MUSOR_PDB_ROOT=/path/to/pdbs   override scan root (default: 6QZH/pdbs)
  MUSOR_SCORE_OUT=/path/to/dir   override output dir for TSV/log (default: 6QZH/)
"""
from __future__ import annotations

import concurrent.futures
import os
import sys
import traceback
from pathlib import Path

# pack_analysis_muOR_ver3 and sasa_sourcecode live in muOR/bin (walk up until found)


def _find_muor_root() -> Path:
    start = Path(__file__).resolve().parent
    for ancestor in [start, *start.parents]:
        if (ancestor / "bin" / "pack_analysis_muOR_ver3.py").is_file():
            return ancestor
    raise RuntimeError(
        "Cannot find muOR root (need bin/pack_analysis_muOR_ver3.py in a parent directory). "
        "Use a full muOR checkout or extend PYTHONPATH to muOR/bin."
    )


_MUOR_ROOT = _find_muor_root()
_BIN = _MUOR_ROOT / "bin"
if str(_BIN) not in sys.path:
    sys.path.insert(0, str(_BIN))  # worker processes need this for `pack_analysis_muOR_ver3`


def _process_one(pdb_path: str) -> tuple:
    """Returns (relpath, finalscore, n_events, err_msg)."""
    from pack_analysis_muOR_ver3 import excute, parser

    rel = pdb_path
    try:
        st = parser.get_structure(Path(pdb_path).stem, pdb_path)[0]
        chainA = st["A"]
        chainB = st["B"]

        def rr(chain, start, end):
            out = []
            for r in range(start, end + 1):
                try:
                    out.append(chain[int(r)])
                except KeyError:
                    pass
            return out

        TM5 = rr(chainA, 225, 263)
        TM6 = rr(chainA, 268, 309)
        TMB = rr(chainB, 346, 384)
        if not TM5 or not TM6 or not TMB:
            return (
                rel,
                "",
                "",
                "empty_TM5_TM6_or_TMB",
            )

        fs, res = excute(protein=st, TM5=TM5, TM6=TM6, TMB=TMB, x=0.0, y=7.0)
        return (rel, str(fs), str(len(res)), "")
    except Exception:
        return (rel, "", "", traceback.format_exc())


def main() -> None:
    _root_6qzh = Path(__file__).resolve().parents[1]
    root = Path(os.environ.get("MUSOR_PDB_ROOT", _root_6qzh / "pdbs"))
    out_dir = Path(os.environ.get("MUSOR_SCORE_OUT", _root_6qzh))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "packing_ver3_batch_summary.tsv"
    out_err = out_dir / "packing_ver3_batch_errors.log"

    pdbs = sorted(root.rglob("*.pdb"))
    n = len(pdbs)
    if n == 0:
        print(f"No .pdb under {root}", file=sys.stderr)
        sys.exit(1)

    jobs = int(os.environ.get("MUSOR_BATCH_JOBS", min(8, os.cpu_count() or 4)))
    jobs = max(1, jobs)

    print(f"Root: {root}")
    print(f"PDB count: {n}")
    print(f"Workers: {jobs}")
    print(f"Output: {out_tsv}")
    print(f"Errors: {out_err}")

    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("pdb_path\tfinalscore\tn_events\terror\n")

    err_f = open(out_err, "w", encoding="utf-8")

    done = 0
    root_res = root.resolve()
    with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(_process_one, str(p)): p for p in pdbs}
        for fut in concurrent.futures.as_completed(futs):
            rel, fs, ne, err = fut.result()
            try:
                rel_out = str(Path(rel).resolve().relative_to(root_res))
            except ValueError:
                rel_out = rel
            done += 1
            err_cell = err.replace("\t", " ").replace("\n", " | ") if err else ""
            line = f"{rel_out}\t{fs}\t{ne}\t{err_cell}\n"
            with open(out_tsv, "a", encoding="utf-8") as f:
                f.write(line)
            if err:
                err_f.write(f"=== {rel} ===\n{err}\n")
                err_f.flush()
            if done % 50 == 0 or done == n:
                print(f"progress {done}/{n}", flush=True)

    err_f.close()
    print("Done.")


if __name__ == "__main__":
    main()
