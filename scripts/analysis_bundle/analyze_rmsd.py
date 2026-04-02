#!/usr/bin/env python

import os
import glob
import argparse
import csv

from Bio.PDB import PDBParser
import MDAnalysis as mda
from MDAnalysis.analysis import rms


def find_relax_pdb(relax_dir, base_name):
    """
    Find the Rosetta Relax PDB corresponding to base_name in relax_dir.
    Prefers base_name_relax_*.pdb, otherwise any PDB containing base_name.
    """
    pattern1 = os.path.join(relax_dir, f"{base_name}_relax_*.pdb")
    files = glob.glob(pattern1)
    if files:
        return sorted(files)[0]

    pattern2 = os.path.join(relax_dir, f"*{base_name}*.pdb")
    files = glob.glob(pattern2)
    if files:
        return sorted(files)[0]

    return None


def sequences_compatible(pdb1, pdb2):
    """
    Quick sanity check: same number of CA atoms in two structures.
    """
    parser = PDBParser(QUIET=True)
    s1 = parser.get_structure("af3", pdb1)
    s2 = parser.get_structure("relax", pdb2)

    def count_ca(struct):
        n = 0
        for model in struct:
            for chain in model:
                for res in chain:
                    if "CA" in res:
                        n += 1
        return n

    n1 = count_ca(s1)
    n2 = count_ca(s2)

    return n1 == n2, n1, n2


def compute_backbone_rmsd(af3_pdb, relax_pdb):
    """
    Compute backbone (N, CA, C) RMSD using MDAnalysis.
    """
    u_ref = mda.Universe(af3_pdb)
    u_mob = mda.Universe(relax_pdb)

    sel = "name N CA C"
    ref_atoms = u_ref.select_atoms(sel)
    mob_atoms = u_mob.select_atoms(sel)

    if len(ref_atoms) != len(mob_atoms):
        raise ValueError(
            f"Backbone atom count mismatch: {len(ref_atoms)} vs {len(mob_atoms)}"
        )

    R = rms.RMSD(mob_atoms, ref_atoms, center=True, superposition=True)
    R.run()
    # R.rmsd: (n_frames, 4), third column is RMSD (Å)
    rmsd_value = R.rmsd[-1, 2]
    return rmsd_value


def main():
    parser = argparse.ArgumentParser(
        description="Compute backbone RMSD between AF3 and Rosetta FastRelax structures."
    )
    parser.add_argument("--af3_dir", required=True, help="Directory with AF3 PDB files")
    parser.add_argument(
        "--relax_dir", required=True, help="Directory with Rosetta FastRelax PDB files"
    )
    parser.add_argument(
        "--out_csv", required=True, help="Output CSV file for RMSD results"
    )
    args = parser.parse_args()

    af3_dir = args.af3_dir
    relax_dir = args.relax_dir
    out_csv = args.out_csv

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)

    af3_pdbs = sorted(glob.glob(os.path.join(af3_dir, "*.pdb")))
    if not af3_pdbs:
        print(f"No AF3 PDB files found in {af3_dir}")
        return

    rows = []

    for af3_pdb in af3_pdbs:
        base = os.path.basename(af3_pdb)
        base_name = os.path.splitext(base)[0]

        relax_pdb = find_relax_pdb(relax_dir, base_name)
        if relax_pdb is None:
            print(f"[WARN] No relax PDB found for {base_name}, skip")
            continue

        compatible, n1, n2 = sequences_compatible(af3_pdb, relax_pdb)
        if not compatible:
            print(
                f"[WARN] CA count mismatch for {base_name}: AF3={n1}, Relax={n2}, skip"
            )
            continue

        try:
            rmsd_value = compute_backbone_rmsd(af3_pdb, relax_pdb)
            print(f"{base_name}: RMSD = {rmsd_value:.3f} Å")
            rows.append(
                [
                    base_name,
                    os.path.basename(af3_pdb),
                    os.path.basename(relax_pdb),
                    n1,
                    rmsd_value,
                ]
            )
        except Exception as e:
            print(f"[ERROR] Failed to compute RMSD for {base_name}: {e}")

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["design_name", "af3_pdb", "relax_pdb", "n_CA", "backbone_rmsd_A"]
        )
        for r in rows:
            writer.writerow(r)

    print(f"RMSD results written to {out_csv}")


if __name__ == "__main__":
    main()

