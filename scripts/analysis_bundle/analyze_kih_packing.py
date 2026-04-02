#!/usr/bin/env python
"""
Simplified Knob-into-Holes (KIH) Packing Analysis for Helix-Helix Interfaces.

Based on the concept that good helix packing involves side chains (knobs) 
fitting into holes formed by backbone atoms of the opposing helix.

This script analyzes:
1. Interface contacts between two helical chains
2. Side chain packing quality (knob-hole interactions)
3. Helix-helix geometry (crossing angle, distance)
"""

import argparse
import os
import glob
import csv
import numpy as np
from collections import defaultdict

try:
    from Bio.PDB import PDBParser, NeighborSearch, Selection
    from Bio.PDB.vectors import calc_dihedral, Vector
except ImportError:
    print("ERROR: BioPython required. Install with: pip install biopython")
    exit(1)


# Amino acid classifications
SMALL_AA = {'G', 'A', 'S', 'C', 'P'}
MEDIUM_AA = {'V', 'T', 'D', 'N', 'I', 'L'}
LARGE_AA = {'M', 'F', 'Y', 'W', 'K', 'R', 'E', 'Q', 'H'}
HYDROPHOBIC = {'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P', 'G'}
POLAR = {'S', 'T', 'N', 'Q', 'Y', 'C'}
CHARGED = {'D', 'E', 'K', 'R', 'H'}

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'H', 'HA', '1HA', '2HA', 'HA2', 'HA3'}


class KIHAnalyzer:
    def __init__(self, pdb_path, chain_a='A', chain_b='B'):
        self.pdb_path = pdb_path
        self.chain_a_id = chain_a
        self.chain_b_id = chain_b
        self.parser = PDBParser(QUIET=True)
        self.structure = None
        self.chain_a = None
        self.chain_b = None
        
    def load_structure(self):
        """Load PDB structure."""
        try:
            self.structure = self.parser.get_structure('protein', self.pdb_path)
            model = self.structure[0]
            
            # Get chains
            chains = {c.id: c for c in model.get_chains()}
            
            if self.chain_a_id not in chains:
                # Try to find first two chains
                chain_ids = sorted(chains.keys())
                if len(chain_ids) >= 2:
                    self.chain_a_id = chain_ids[0]
                    self.chain_b_id = chain_ids[1]
            
            self.chain_a = chains.get(self.chain_a_id)
            self.chain_b = chains.get(self.chain_b_id)
            
            if self.chain_a is None or self.chain_b is None:
                return False
            return True
        except Exception as e:
            print(f"Error loading {self.pdb_path}: {e}")
            return False
    
    def get_residues(self, chain):
        """Get list of residues (excluding heteroatoms)."""
        return [r for r in chain.get_residues() 
                if r.id[0] == ' ' and r.has_id('CA')]
    
    def get_ca_coords(self, residues):
        """Get CA coordinates for residues."""
        coords = []
        for r in residues:
            if r.has_id('CA'):
                coords.append(r['CA'].get_coord())
        return np.array(coords)
    
    def calc_helix_axis(self, ca_coords):
        """Fit a line through CA atoms to get helix axis."""
        if len(ca_coords) < 3:
            return None, None
        
        # Center the points
        centroid = np.mean(ca_coords, axis=0)
        centered = ca_coords - centroid
        
        # SVD to find principal axis
        _, _, vh = np.linalg.svd(centered)
        axis = vh[0]  # First principal component
        
        return axis, centroid
    
    def calc_crossing_angle(self):
        """Calculate helix-helix crossing angle."""
        res_a = self.get_residues(self.chain_a)
        res_b = self.get_residues(self.chain_b)
        
        ca_a = self.get_ca_coords(res_a)
        ca_b = self.get_ca_coords(res_b)
        
        if len(ca_a) < 5 or len(ca_b) < 5:
            return None
        
        axis_a, _ = self.calc_helix_axis(ca_a)
        axis_b, _ = self.calc_helix_axis(ca_b)
        
        if axis_a is None or axis_b is None:
            return None
        
        # Calculate angle between axes
        dot = np.abs(np.dot(axis_a, axis_b))
        angle = np.degrees(np.arccos(np.clip(dot, -1, 1)))
        
        # Return smaller angle (0-90)
        return min(angle, 180 - angle)
    
    def calc_interface_distance(self):
        """Calculate average CA-CA distance at interface."""
        res_a = self.get_residues(self.chain_a)
        res_b = self.get_residues(self.chain_b)
        
        ca_a = self.get_ca_coords(res_a)
        ca_b = self.get_ca_coords(res_b)
        
        if len(ca_a) == 0 or len(ca_b) == 0:
            return None, None
        
        # Calculate all pairwise distances
        dists = []
        for ca1 in ca_a:
            for ca2 in ca_b:
                d = np.linalg.norm(ca1 - ca2)
                if d < 12:  # Only interface residues
                    dists.append(d)
        
        if not dists:
            return None, None
        
        return np.mean(dists), np.min(dists)
    
    def identify_interface_residues(self, cutoff=8.0):
        """
        Find residues at the interface using CA-CA distance cutoff.

        Accelerated with Bio.PDB.NeighborSearch to avoid O(N^2) residue-pair loops.
        Returns:
          (interface_a_residues, interface_b_residues) as lists of residue objects.
        """
        res_a = self.get_residues(self.chain_a)
        res_b = self.get_residues(self.chain_b)

        ca_b_atoms = []
        ca_b_atom_to_residue = {}
        for rb in res_b:
            if rb.has_id('CA'):
                atom = rb['CA']
                ca_b_atoms.append(atom)
                ca_b_atom_to_residue[id(atom)] = rb

        if not ca_b_atoms:
            return [], []

        ns = NeighborSearch(ca_b_atoms)

        interface_a_map = {}
        interface_b_map = {}

        for ra in res_a:
            if not ra.has_id('CA'):
                continue
            ca_a = ra['CA']
            neigh_atoms = ns.search(ca_a.get_coord(), cutoff)
            if not neigh_atoms:
                continue

            # Mark chain A residue
            interface_a_map[id(ra)] = ra

            # Mark chain B residues
            for nb in neigh_atoms:
                rb = ca_b_atom_to_residue.get(id(nb))
                if rb is not None:
                    interface_b_map[id(rb)] = rb

        return list(interface_a_map.values()), list(interface_b_map.values())
    
    def get_sidechain_atoms(self, residue):
        """Get sidechain (non-backbone) atoms."""
        return [a for a in residue.get_atoms() 
                if a.get_name() not in BACKBONE_ATOMS]
    
    def get_backbone_atoms(self, residue):
        """Get backbone atoms."""
        return [a for a in residue.get_atoms() 
                if a.get_name() in BACKBONE_ATOMS]
    
    def calc_sidechain_to_backbone_contacts(self, sc_atoms, bb_atoms, cutoff=4.0):
        """Count contacts between sidechain and backbone atoms."""
        contacts = 0
        min_dist = float('inf')
        
        for sc in sc_atoms:
            sc_coord = sc.get_coord()
            for bb in bb_atoms:
                bb_coord = bb.get_coord()
                d = np.linalg.norm(sc_coord - bb_coord)
                if d < cutoff:
                    contacts += 1
                if d < min_dist:
                    min_dist = d
        
        return contacts, min_dist if min_dist != float('inf') else None
    
    def analyze_knob_hole_packing(self, interface_a, interface_b):
        """
        Analyze KIH packing between chains.
        
        A "knob" is a side chain that packs into a "hole" formed by 
        backbone atoms of the opposing helix.
        
        Returns:
            dict with packing metrics
        """
        if not interface_a or not interface_b:
            return {
                'n_knobs_a_to_b': 0,
                'n_knobs_b_to_a': 0,
                'total_knobs': 0,
                'knob_score': 0,
                'avg_knob_contacts': 0,
                'interface_res_a': 0,
                'interface_res_b': 0,
            }
        
        # Contact cutoffs (Å)
        CONTACT_CUTOFF = 4.0
        MIN_DIST_CUTOFF = 3.5

        # Build NeighborSearch indices for backbone atoms (much faster than residue-pair loops)
        bb_b_atoms = []
        for res_b in interface_b:
            bb_b_atoms.extend(self.get_backbone_atoms(res_b))
        ns_bb_b = NeighborSearch(bb_b_atoms) if bb_b_atoms else None
        bb_a_atoms = []
        for res_a in interface_a:
            bb_a_atoms.extend(self.get_backbone_atoms(res_a))
        ns_bb_a = NeighborSearch(bb_a_atoms) if bb_a_atoms else None

        # Analyze knobs from A going into holes in B
        knobs_a_to_b = []
        if ns_bb_b is not None:
            for res_a in interface_a:
                sc_a = self.get_sidechain_atoms(res_a)
                if not sc_a:  # GLY
                    continue

                contacted_bb_ids = set()
                min_dist_to_bb = float('inf')

                for sc_atom in sc_a:
                    neigh_bb = ns_bb_b.search(sc_atom.get_coord(), CONTACT_CUTOFF)
                    for nb in neigh_bb:
                        contacted_bb_ids.add(id(nb))
                        d = np.linalg.norm(sc_atom.get_coord() - nb.get_coord())
                        if d < MIN_DIST_CUTOFF and d < min_dist_to_bb:
                            min_dist_to_bb = d

                total_contacts = len(contacted_bb_ids)
                if total_contacts >= 2 or min_dist_to_bb < MIN_DIST_CUTOFF:
                    knobs_a_to_b.append({
                        'residue': f"{res_a.get_resname()}{res_a.id[1]}",
                        'contacts': total_contacts,
                        'min_dist': min_dist_to_bb,
                    })

        # Analyze knobs from B going into holes in A
        knobs_b_to_a = []
        if ns_bb_a is not None:
            for res_b in interface_b:
                sc_b = self.get_sidechain_atoms(res_b)
                if not sc_b:  # GLY
                    continue

                contacted_bb_ids = set()
                min_dist_to_bb = float('inf')

                for sc_atom in sc_b:
                    neigh_bb = ns_bb_a.search(sc_atom.get_coord(), CONTACT_CUTOFF)
                    for nb in neigh_bb:
                        contacted_bb_ids.add(id(nb))
                        d = np.linalg.norm(sc_atom.get_coord() - nb.get_coord())
                        if d < MIN_DIST_CUTOFF and d < min_dist_to_bb:
                            min_dist_to_bb = d

                total_contacts = len(contacted_bb_ids)
                if total_contacts >= 2 or min_dist_to_bb < MIN_DIST_CUTOFF:
                    knobs_b_to_a.append({
                        'residue': f"{res_b.get_resname()}{res_b.id[1]}",
                        'contacts': total_contacts,
                        'min_dist': min_dist_to_bb,
                    })
        
        # Calculate scores
        n_knobs_a = len(knobs_a_to_b)
        n_knobs_b = len(knobs_b_to_a)
        total_knobs = n_knobs_a + n_knobs_b
        
        all_contacts = ([k['contacts'] for k in knobs_a_to_b] + 
                       [k['contacts'] for k in knobs_b_to_a])
        avg_contacts = np.mean(all_contacts) if all_contacts else 0
        
        # Knob score: weighted by contacts and distance
        knob_score = 0
        for k in knobs_a_to_b + knobs_b_to_a:
            if k['min_dist'] < 3.1:
                knob_score += 1.0
            elif k['contacts'] >= 3:
                knob_score += 0.75
            elif k['contacts'] >= 2:
                knob_score += 0.5
        
        return {
            'n_knobs_a_to_b': n_knobs_a,
            'n_knobs_b_to_a': n_knobs_b,
            'total_knobs': total_knobs,
            'knob_score': round(knob_score, 2),
            'avg_knob_contacts': round(avg_contacts, 2),
            'interface_res_a': len(interface_a),
            'interface_res_b': len(interface_b),
            'knobs_a_detail': knobs_a_to_b,
            'knobs_b_detail': knobs_b_to_a,
        }
    
    def calc_hydrophobic_packing(self, interface_a, interface_b):
        """Calculate hydrophobic contacts at interface."""
        if not interface_a or not interface_b:
            return 0, 0.0

        hp_a_set = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO'}

        # Use CA atoms and a neighbor search to avoid N^2 distance loops
        ca_b_atoms = []
        ca_b_atom_to_is_hp = {}
        for res_b in interface_b:
            if not res_b.has_id('CA'):
                continue
            atom = res_b['CA']
            ca_b_atoms.append(atom)
            ca_b_atom_to_is_hp[id(atom)] = res_b.get_resname() in hp_a_set

        if not ca_b_atoms:
            return 0, 0.0

        ns = NeighborSearch(ca_b_atoms)

        hp_contacts = 0
        total_contacts = 0

        CONTACT_CA_CUTOFF = 8.0
        for res_a in interface_a:
            if not res_a.has_id('CA'):
                continue
            aa_a = res_a.get_resname()
            is_hp_a = aa_a in hp_a_set
            ca_a = res_a['CA']

            neigh_atoms = ns.search(ca_a.get_coord(), CONTACT_CA_CUTOFF)
            for nb in neigh_atoms:
                total_contacts += 1
                is_hp_b = ca_b_atom_to_is_hp.get(id(nb), False)
                if is_hp_a and is_hp_b:
                    hp_contacts += 1

        hp_fraction = hp_contacts / total_contacts if total_contacts > 0 else 0.0
        return hp_contacts, hp_fraction
    
    def analyze(self):
        """Run full analysis and return results dict."""
        if not self.load_structure():
            return None
        
        results = {
            'pdb_file': os.path.basename(self.pdb_path),
            'chain_a': self.chain_a_id,
            'chain_b': self.chain_b_id,
        }
        
        # Chain lengths
        res_a = self.get_residues(self.chain_a)
        res_b = self.get_residues(self.chain_b)
        results['len_chain_a'] = len(res_a)
        results['len_chain_b'] = len(res_b)
        
        # Crossing angle
        crossing_angle = self.calc_crossing_angle()
        results['crossing_angle'] = round(crossing_angle, 2) if crossing_angle else None
        
        # Interface distance
        avg_dist, min_dist = self.calc_interface_distance()
        results['avg_interface_dist'] = round(avg_dist, 2) if avg_dist else None
        results['min_interface_dist'] = round(min_dist, 2) if min_dist else None
        
        # Interface set once (used by both KIH packing and hydrophobic packing)
        interface_a, interface_b = self.identify_interface_residues(cutoff=10.0)

        # KIH packing analysis
        kih = self.analyze_knob_hole_packing(interface_a, interface_b)
        results.update({
            'n_knobs_a_to_b': kih['n_knobs_a_to_b'],
            'n_knobs_b_to_a': kih['n_knobs_b_to_a'],
            'total_knobs': kih['total_knobs'],
            'knob_score': kih['knob_score'],
            'avg_knob_contacts': kih['avg_knob_contacts'],
            'interface_res_a': kih['interface_res_a'],
            'interface_res_b': kih['interface_res_b'],
        })
        
        # Hydrophobic packing
        hp_contacts, hp_fraction = self.calc_hydrophobic_packing(interface_a, interface_b)
        results['hydrophobic_contacts'] = hp_contacts
        results['hydrophobic_fraction'] = round(hp_fraction, 3)
        
        # Normalized scores
        max_possible_knobs = min(results['interface_res_a'], results['interface_res_b'])
        results['knob_density'] = round(
            results['total_knobs'] / max_possible_knobs, 3
        ) if max_possible_knobs > 0 else 0
        
        return results


def analyze_directory(
    pdb_dir,
    chain_a='A',
    chain_b='B',
    out_csv='kih_packing_results.csv',
    max_pdbs=0,
):
    """Analyze all PDBs in a directory."""
    pdb_files = glob.glob(os.path.join(pdb_dir, '*.pdb'))
    
    if not pdb_files:
        print(f"No PDB files found in {pdb_dir}")
        return
    
    print(f"Found {len(pdb_files)} PDB files to analyze...")
    
    results = []
    for i, pdb_path in enumerate(sorted(pdb_files), 1):
        if max_pdbs and max_pdbs > 0 and i > max_pdbs:
            break
        if i % 50 == 0:
            print(f"  Processing {i}/{len(pdb_files)}...")
        
        analyzer = KIHAnalyzer(pdb_path, chain_a, chain_b)
        result = analyzer.analyze()
        
        if result:
            results.append(result)
    
    if not results:
        print("No valid results to write")
        return
    
    # Write CSV
    fieldnames = [
        'pdb_file', 'chain_a', 'chain_b', 'len_chain_a', 'len_chain_b',
        'crossing_angle', 'avg_interface_dist', 'min_interface_dist',
        'interface_res_a', 'interface_res_b',
        'n_knobs_a_to_b', 'n_knobs_b_to_a', 'total_knobs',
        'knob_score', 'avg_knob_contacts', 'knob_density',
        'hydrophobic_contacts', 'hydrophobic_fraction'
    ]
    
    with open(out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r.get(k, '') for k in fieldnames})
    
    print(f"\nResults written to {out_csv}")
    print(f"Analyzed {len(results)} structures")
    
    # Print summary
    if results:
        knob_scores = [r['knob_score'] for r in results if r['knob_score'] is not None]
        total_knobs = [r['total_knobs'] for r in results]
        
        print(f"\n--- Summary ---")
        print(f"Knob Score:  min={min(knob_scores):.1f}, max={max(knob_scores):.1f}, "
              f"mean={np.mean(knob_scores):.2f}")
        print(f"Total Knobs: min={min(total_knobs)}, max={max(total_knobs)}, "
              f"mean={np.mean(total_knobs):.1f}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze Knob-into-Holes (KIH) packing in helix-helix interfaces'
    )
    parser.add_argument('--pdb_dir', default=None,
                       help='Directory containing PDB files (required when --single is not set)')
    parser.add_argument('--chain_a', default='A',
                       help='Chain ID for first helix (default: A)')
    parser.add_argument('--chain_b', default='B',
                       help='Chain ID for second helix (default: B)')
    parser.add_argument('--out_csv', default='kih_packing_results.csv',
                       help='Output CSV file')
    parser.add_argument('--single', default=None,
                       help='Analyze a single PDB file (verbose output)')
    parser.add_argument('--max_pdbs', type=int, default=0,
                       help='Limit number of PDBs to analyze (0 = all)')
    
    args = parser.parse_args()
    
    if args.single:
        # Analyze single PDB with verbose output
        analyzer = KIHAnalyzer(args.single, args.chain_a, args.chain_b)
        result = analyzer.analyze()
        
        if result:
            print(f"\n=== KIH Packing Analysis: {result['pdb_file']} ===\n")
            
            print(f"Chains: {result['chain_a']} ({result['len_chain_a']} res) vs "
                  f"{result['chain_b']} ({result['len_chain_b']} res)")
            print(f"Crossing Angle: {result['crossing_angle']}°")
            print(f"Interface Distance: avg={result['avg_interface_dist']} Å, "
                  f"min={result['min_interface_dist']} Å")
            
            print(f"\n--- Interface ---")
            print(f"Interface residues: A={result['interface_res_a']}, B={result['interface_res_b']}")
            
            print(f"\n--- KIH Packing ---")
            print(f"Knobs A→B: {result['n_knobs_a_to_b']}")
            print(f"Knobs B→A: {result['n_knobs_b_to_a']}")
            print(f"Total Knobs: {result['total_knobs']}")
            print(f"Knob Score: {result['knob_score']}")
            print(f"Knob Density: {result['knob_density']}")
            print(f"Avg Contacts per Knob: {result['avg_knob_contacts']}")
            
            print(f"\n--- Hydrophobic ---")
            print(f"HP Contacts: {result['hydrophobic_contacts']}")
            print(f"HP Fraction: {result['hydrophobic_fraction']}")
    else:
        # Analyze directory
        if not args.pdb_dir:
            raise SystemExit("ERROR: --pdb_dir is required when --single is not set")
        analyze_directory(
            args.pdb_dir,
            args.chain_a,
            args.chain_b,
            args.out_csv,
            max_pdbs=args.max_pdbs,
        )


if __name__ == '__main__':
    main()
