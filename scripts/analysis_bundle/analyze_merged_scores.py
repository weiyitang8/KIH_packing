#!/usr/bin/env python
"""
Merge MPNN results with Rosetta scores and generate analysis plots.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re

def parse_rosetta_scores(score_file):
    """Parse Rosetta score file into DataFrame."""
    rows = []
    header = None
    
    with open(score_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("SEQUENCE:"):
                continue
            if line.startswith("SCORE:") and header is None:
                parts = line.split()
                header = parts[1:]  # skip "SCORE:"
                continue
            if line.startswith("SCORE:"):
                parts = line.split()
                values = parts[1:]  # skip "SCORE:"
                if len(values) == len(header):
                    rows.append(values)
    
    df = pd.DataFrame(rows, columns=header)
    
    # Convert numeric columns
    for col in df.columns:
        if col != 'description':
            try:
                df[col] = pd.to_numeric(df[col])
            except:
                pass
    
    # Extract design and n from description (e.g., "design0_n0_design0_n0_0001")
    def extract_design_n(desc):
        match = re.match(r'design(\d+)_n(\d+)', desc)
        if match:
            return int(match.group(1)), int(match.group(2))
        return None, None
    
    df['design'], df['n'] = zip(*df['description'].apply(extract_design_n))
    
    return df

def load_mpnn_results(csv_file):
    """Load MPNN results CSV."""
    df = pd.read_csv(csv_file)
    return df

def merge_data(mpnn_df, ros_df):
    """Merge MPNN and Rosetta data on design/n."""
    merged = pd.merge(mpnn_df, ros_df, on=['design', 'n'], how='inner')
    return merged

def calculate_composite_score(df):
    """
    Calculate a composite score for ranking designs.
    Higher is better.
    """
    # Normalize each metric to 0-1 range
    def normalize(series, higher_is_better=True):
        min_val, max_val = series.min(), series.max()
        if max_val == min_val:
            return pd.Series([0.5] * len(series))
        normalized = (series - min_val) / (max_val - min_val)
        return normalized if higher_is_better else 1 - normalized
    
    scores = pd.DataFrame()
    
    # MPNN metrics (from AF2/ColabFold screening)
    if 'plddt' in df.columns:
        scores['plddt_norm'] = normalize(df['plddt'], higher_is_better=True)
    if 'i_ptm' in df.columns:
        scores['i_ptm_norm'] = normalize(df['i_ptm'], higher_is_better=True)
    if 'i_pae' in df.columns:
        scores['i_pae_norm'] = normalize(df['i_pae'], higher_is_better=False)
    if 'rmsd' in df.columns:
        scores['rmsd_norm'] = normalize(df['rmsd'], higher_is_better=False)
    
    # Rosetta metrics
    if 'total_score' in df.columns:
        scores['total_score_norm'] = normalize(df['total_score'], higher_is_better=False)
    if 'dG_separated' in df.columns:
        scores['dG_norm'] = normalize(df['dG_separated'], higher_is_better=False)
    if 'sc_value' in df.columns:
        scores['sc_norm'] = normalize(df['sc_value'], higher_is_better=True)
    if 'packstat' in df.columns:
        scores['packstat_norm'] = normalize(df['packstat'], higher_is_better=True)
    if 'hbonds_int' in df.columns:
        scores['hbonds_norm'] = normalize(df['hbonds_int'], higher_is_better=True)
    if 'dSASA_int' in df.columns:
        scores['dSASA_norm'] = normalize(df['dSASA_int'], higher_is_better=True)
    if 'holes_chA_chB' in df.columns:
        scores['holes_norm'] = normalize(df['holes_chA_chB'], higher_is_better=False)
    if 'delta_unsatHbonds' in df.columns:
        scores['unsat_norm'] = normalize(df['delta_unsatHbonds'], higher_is_better=False)
    
    # Weighted composite score
    weights = {
        'plddt_norm': 1.0,
        'i_ptm_norm': 1.5,      # Interface pTM is important
        'i_pae_norm': 1.0,
        'rmsd_norm': 1.0,       # RMSD from template
        'total_score_norm': 0.5,
        'dG_norm': 2.0,         # Binding energy is key
        'sc_norm': 1.5,         # Shape complementarity important
        'packstat_norm': 1.0,
        'hbonds_norm': 1.0,
        'dSASA_norm': 0.5,
        'holes_norm': 1.0,      # Holes/voids (lower is better)
        'unsat_norm': 1.0,      # Unsatisfied H-bonds (lower is better)
    }
    
    composite = pd.Series([0.0] * len(df))
    total_weight = 0
    for col, weight in weights.items():
        if col in scores.columns:
            composite += scores[col] * weight
            total_weight += weight
    
    df['composite_score'] = composite / total_weight if total_weight > 0 else 0
    
    return df

def generate_plots(df, out_dir):
    """Generate analysis plots."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # 1. Correlation heatmap of key metrics
    key_cols = ['mpnn', 'plddt', 'i_ptm', 'i_pae', 'rmsd', 
                'total_score', 'dG_separated', 'sc_value', 'packstat', 
                'hbonds_int', 'dSASA_int', 'holes_chA_chB', 'delta_unsatHbonds',
                'composite_score']
    available_cols = [c for c in key_cols if c in df.columns]
    
    if len(available_cols) > 2:
        fig, ax = plt.subplots(figsize=(12, 10))
        corr = df[available_cols].corr()
        sns.heatmap(corr, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                    square=True, ax=ax)
        plt.title('Correlation between MPNN and Rosetta Metrics')
        plt.tight_layout()
        plt.savefig(out_dir / 'correlation_heatmap.png', dpi=150)
        plt.close()
    
    # 2. dG_separated vs i_ptm (key binding metrics)
    if 'dG_separated' in df.columns and 'i_ptm' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(df['dG_separated'], df['i_ptm'], 
                            c=df['composite_score'], cmap='viridis',
                            alpha=0.6, s=30)
        plt.colorbar(scatter, label='Composite Score')
        ax.set_xlabel('dG_separated (Rosetta, kcal/mol)')
        ax.set_ylabel('i_pTM (AF2)')
        ax.set_title('Binding Energy vs Interface pTM')
        
        # Mark top 10
        top10 = df.nlargest(10, 'composite_score')
        ax.scatter(top10['dG_separated'], top10['i_ptm'], 
                  c='red', s=100, marker='*', label='Top 10', edgecolors='black')
        ax.legend()
        plt.tight_layout()
        plt.savefig(out_dir / 'dG_vs_iptm.png', dpi=150)
        plt.close()
    
    # 3. sc_value vs packstat (structure quality)
    if 'sc_value' in df.columns and 'packstat' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(df['sc_value'], df['packstat'], 
                            c=df['dG_separated'] if 'dG_separated' in df.columns else df['composite_score'],
                            cmap='viridis_r', alpha=0.6, s=30)
        plt.colorbar(scatter, label='dG_separated' if 'dG_separated' in df.columns else 'Composite')
        ax.set_xlabel('Shape Complementarity (sc_value)')
        ax.set_ylabel('Packing Quality (packstat)')
        ax.axvline(x=0.65, color='r', linestyle='--', alpha=0.5, label='sc > 0.65')
        ax.axhline(y=0.55, color='r', linestyle='--', alpha=0.5, label='pack > 0.55')
        ax.set_title('Shape Complementarity vs Packing')
        plt.tight_layout()
        plt.savefig(out_dir / 'sc_vs_packstat.png', dpi=150)
        plt.close()
    
    # 4. Distribution of key metrics
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    metrics_to_plot = [
        ('dG_separated', 'Binding Energy (kcal/mol)'),
        ('i_ptm', 'Interface pTM'),
        ('sc_value', 'Shape Complementarity'),
        ('packstat', 'Packing Quality'),
        ('plddt', 'pLDDT'),
        ('rmsd', 'RMSD (Å)'),
        ('holes_chA_chB', 'Holes Score'),
        ('delta_unsatHbonds', 'Unsatisfied H-bonds'),
        ('composite_score', 'Composite Score')
    ]
    
    for ax, (col, label) in zip(axes.flat, metrics_to_plot):
        if col in df.columns:
            ax.hist(df[col], bins=30, edgecolor='black', alpha=0.7)
            ax.axvline(df[col].median(), color='r', linestyle='--', label=f'Median: {df[col].median():.2f}')
            ax.set_xlabel(label)
            ax.set_ylabel('Count')
            ax.legend()
    
    plt.suptitle('Distribution of Key Metrics')
    plt.tight_layout()
    plt.savefig(out_dir / 'distributions.png', dpi=150)
    plt.close()
    
    # 5. RMSD vs dG plot
    if 'rmsd' in df.columns and 'dG_separated' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(df['rmsd'], df['dG_separated'], 
                            c=df['composite_score'], cmap='viridis',
                            alpha=0.6, s=30)
        plt.colorbar(scatter, label='Composite Score')
        ax.set_xlabel('RMSD from template (Å)')
        ax.set_ylabel('dG_separated (kcal/mol)')
        ax.set_title('RMSD vs Binding Energy')
        
        # Mark top 10
        top10 = df.nlargest(10, 'composite_score')
        ax.scatter(top10['rmsd'], top10['dG_separated'], 
                  c='red', s=100, marker='*', label='Top 10', edgecolors='black')
        ax.legend()
        plt.tight_layout()
        plt.savefig(out_dir / 'rmsd_vs_dG.png', dpi=150)
        plt.close()
    
    # 6. Holes vs Unsatisfied H-bonds plot
    if 'holes_chA_chB' in df.columns and 'delta_unsatHbonds' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(df['holes_chA_chB'], df['delta_unsatHbonds'], 
                            c=df['dG_separated'] if 'dG_separated' in df.columns else df['composite_score'],
                            cmap='viridis_r', alpha=0.6, s=30)
        plt.colorbar(scatter, label='dG_separated' if 'dG_separated' in df.columns else 'Composite')
        ax.set_xlabel('Holes Score (lower is better)')
        ax.set_ylabel('Unsatisfied H-bonds (lower is better)')
        ax.axvline(x=3.0, color='r', linestyle='--', alpha=0.5, label='holes < 3')
        ax.axhline(y=4.0, color='r', linestyle='--', alpha=0.5, label='unsat < 4')
        ax.set_title('Interface Quality: Holes vs Unsatisfied H-bonds')
        ax.legend()
        plt.tight_layout()
        plt.savefig(out_dir / 'holes_vs_unsatHbonds.png', dpi=150)
        plt.close()
    
    # 7. Top candidates table
    if 'composite_score' in df.columns:
        top20 = df.nlargest(20, 'composite_score')
        cols_to_show = ['design', 'n', 'plddt', 'i_ptm', 'rmsd',
                        'dG_separated', 'sc_value', 'packstat', 
                        'holes_chA_chB', 'delta_unsatHbonds', 'composite_score']
        cols_available = [c for c in cols_to_show if c in top20.columns]
        
        fig, ax = plt.subplots(figsize=(18, 8))
        ax.axis('off')
        table = ax.table(cellText=top20[cols_available].round(3).values,
                        colLabels=cols_available,
                        loc='center',
                        cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.2, 1.5)
        plt.title('Top 20 Candidates by Composite Score', fontsize=14, pad=20)
        plt.tight_layout()
        plt.savefig(out_dir / 'top20_table.png', dpi=150, bbox_inches='tight')
        plt.close()
    
    print(f"Plots saved to {out_dir}")

def print_summary(df):
    """Print summary statistics and recommendations."""
    print("\n" + "="*70)
    print("ANALYSIS SUMMARY")
    print("="*70)
    
    print(f"\nTotal designs analyzed: {len(df)}")
    
    # Key metrics summary
    key_metrics = {
        'dG_separated': ('Binding Energy', 'kcal/mol', 'lower is better, < -30 good'),
        'i_ptm': ('Interface pTM', '', 'higher is better, > 0.5 good'),
        'sc_value': ('Shape Complementarity', '', 'higher is better, > 0.65 good'),
        'packstat': ('Packing Quality', '', 'higher is better, > 0.55 good'),
        'plddt': ('pLDDT', '', 'higher is better, > 0.7 good'),
        'rmsd': ('RMSD', 'Å', 'lower is better, < 5 good'),
        'holes_chA_chB': ('Holes Score', '', 'lower is better, < 3 good'),
        'delta_unsatHbonds': ('Unsatisfied H-bonds', '', 'lower is better, < 4 good'),
    }
    
    print("\n--- Key Metrics Summary ---")
    for col, (name, unit, interpret) in key_metrics.items():
        if col in df.columns:
            print(f"\n{name}:")
            print(f"  Range: {df[col].min():.2f} to {df[col].max():.2f} {unit}")
            print(f"  Mean:  {df[col].mean():.2f}, Median: {df[col].median():.2f}")
            print(f"  Interpretation: {interpret}")
    
    # Good candidates
    print("\n--- Filtering Good Candidates ---")
    
    good_mask = pd.Series([True] * len(df))
    criteria = []
    
    if 'dG_separated' in df.columns:
        good_mask &= df['dG_separated'] < -30
        criteria.append('dG < -30')
    if 'sc_value' in df.columns:
        good_mask &= df['sc_value'] > 0.60
        criteria.append('sc > 0.60')
    if 'packstat' in df.columns:
        good_mask &= df['packstat'] > 0.50
        criteria.append('packstat > 0.50')
    if 'i_ptm' in df.columns:
        good_mask &= df['i_ptm'] > 0.15
        criteria.append('i_ptm > 0.15')
    if 'rmsd' in df.columns:
        good_mask &= df['rmsd'] < 15
        criteria.append('rmsd < 15')
    if 'holes_chA_chB' in df.columns:
        good_mask &= df['holes_chA_chB'] < 4
        criteria.append('holes < 4')
    if 'delta_unsatHbonds' in df.columns:
        good_mask &= df['delta_unsatHbonds'] < 6
        criteria.append('unsatHbonds < 6')
    
    n_good = good_mask.sum()
    print(f"Criteria: {', '.join(criteria)}")
    print(f"Designs passing all criteria: {n_good} / {len(df)} ({100*n_good/len(df):.1f}%)")
    
    # Top 10 by composite score
    if 'composite_score' in df.columns:
        print("\n--- Top 10 by Composite Score ---")
        top10 = df.nlargest(10, 'composite_score')
        cols = ['design', 'n', 'composite_score', 'dG_separated', 'i_ptm', 'sc_value', 
                'packstat', 'rmsd', 'holes_chA_chB', 'delta_unsatHbonds']
        cols = [c for c in cols if c in top10.columns]
        print(top10[cols].to_string(index=False))
    
    print("\n" + "="*70)
    print("HOW TO INTERPRET / WHAT TO LOOK FOR:")
    print("="*70)
    print("""
1. MPNN Metrics (from AF2 screening):
   - pLDDT: Confidence in structure prediction (>0.7 is good)
   - i_pTM: Interface predicted TM-score (>0.5 suggests good interface)
   - i_PAE: Interface predicted aligned error (lower is better, <15 good)
   - RMSD: Deviation from template (context-dependent)

2. Rosetta Metrics (from FastRelax):
   - dG_separated: Binding energy. More negative = stronger binding
     * < -30 kcal/mol: Promising binder
     * < -50 kcal/mol: Strong binder
   - sc_value: Shape complementarity (0-1). Higher = better fit
     * > 0.65: Good complementarity
   - packstat: Packing quality (0-1). Higher = tighter packing
     * > 0.55: Well-packed interface
   - hbonds_int: H-bonds at interface. More = more specific binding

3. Composite Score:
   - Weighted combination of normalized metrics
   - Higher is better
   - Use for initial ranking, then examine top candidates manually

4. Recommended Workflow:
   a) Sort by composite_score to get top candidates
   b) Check that top candidates have reasonable individual metrics
   c) Visually inspect top 10-20 structures in PyMOL
   d) Consider experimental validation for top 3-5
""")

def main():
    parser = argparse.ArgumentParser(description='Merge and analyze MPNN + Rosetta scores')
    parser.add_argument('--mpnn', required=True, help='MPNN results CSV file')
    parser.add_argument('--rosetta', required=True, help='Rosetta score file')
    parser.add_argument('--out_dir', default='./analysis_plots', help='Output directory for plots')
    parser.add_argument('--out_csv', default='./merged_scores.csv', help='Output merged CSV')
    args = parser.parse_args()
    
    print("Loading MPNN results...")
    mpnn_df = load_mpnn_results(args.mpnn)
    print(f"  Loaded {len(mpnn_df)} MPNN entries")
    
    print("Loading Rosetta scores...")
    ros_df = parse_rosetta_scores(args.rosetta)
    print(f"  Loaded {len(ros_df)} Rosetta entries")
    
    print("Merging data...")
    merged = merge_data(mpnn_df, ros_df)
    print(f"  Merged: {len(merged)} entries")
    
    if len(merged) == 0:
        print("ERROR: No matching entries found. Check design/n naming.")
        return
    
    print("Calculating composite scores...")
    merged = calculate_composite_score(merged)
    
    print("Generating plots...")
    generate_plots(merged, args.out_dir)
    
    print(f"Saving merged data to {args.out_csv}...")
    merged.to_csv(args.out_csv, index=False)
    
    print_summary(merged)
    
    print(f"\nDone! Check {args.out_dir} for plots and {args.out_csv} for merged data.")

if __name__ == "__main__":
    main()
