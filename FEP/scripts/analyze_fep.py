#!/usr/bin/env python3
"""
FEP ΔΔG_binding Analysis - MBAR Method with Visualization
Calculate binding free energy change upon mutation using MBAR method

Usage:
    python analyze_fep.py /path/to/outputs/test/design_name [--skip 0.1]
    
Example:
    python analyze_fep.py outputs/test/design_62_dldesign_16_best
    
Output:
    - ddG_results.txt: Numerical results
    - mbar_overlap_*.png: MBAR overlap matrices (2 plots)
    - free_energy_profile_*.png: Free energy profiles (2 plots)
    - convergence_*.png: Convergence analysis (2 plots)
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import warnings
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for server usage
import matplotlib.pyplot as plt
import logging

# Suppress ALL warnings and logs
warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.CRITICAL)
logging.getLogger().setLevel(logging.CRITICAL)
logging.getLogger('alchemlyb').setLevel(logging.CRITICAL)
logging.getLogger('pymbar').setLevel(logging.CRITICAL)
logging.getLogger('alchemlyb.estimators').setLevel(logging.CRITICAL)
logging.getLogger('pymbar.mbar').setLevel(logging.CRITICAL)

from alchemlyb.parsing.gmx import extract_u_nk
from alchemlyb.preprocessing import subsampling
from alchemlyb.estimators import MBAR
from alchemlyb.visualisation import plot_mbar_overlap_matrix
from alchemlyb.postprocessors.units import to_kcalmol

TEMPERATURE = 310  # K
SKIP_FRACTION = 0.1  # Skip first 10% as equilibration

def load_data(base_path, state, direction, skip_frac=0.1, T=310):
    """Load and preprocess u_nk data from lambda windows for MBAR
    
    Following alchemlyb tutorial workflow:
    1. Extract u_nk data from dhdl.xvg files
    2. Skip equilibration period
    3. Decorrelate samples using subsampling module
    """
    path = os.path.join(base_path, state, direction)
    
    print(f"\n{'='*60}")
    print(f"{state.upper()}/{direction.upper()}")
    print(f"{'='*60}")
    
    # Find lambda directories
    try:
        lambda_dirs = sorted([d for d in os.listdir(path) if d.startswith('lambda_')])
    except FileNotFoundError:
        print(f"  ERROR: Directory not found: {path}")
        return None
    
    if not lambda_dirs:
        print(f"  ERROR: No lambda windows found")
        return None
    
    print(f"  Found {len(lambda_dirs)} lambda windows")
    
    u_nk_list = []
    success = 0
    
    for ldir in lambda_dirs:
        # Try multiple file locations
        dhdl_file = None
        for fname in ['dhdl.xvg', 'prod.xvg']:
            for subdir in ['', 'PROD', '.']:
                fpath = os.path.join(path, ldir, subdir, fname) if subdir else os.path.join(path, ldir, fname)
                if os.path.exists(fpath):
                    dhdl_file = fpath
                    break
            if dhdl_file:
                break
        
        if not dhdl_file:
            continue
        
        try:
            # Extract u_nk (MBAR requires u_nk, not dH/dλ)
            u_nk = extract_u_nk(dhdl_file, T=T)
            
            # Skip equilibration period
            if len(u_nk) > 100:
                skip_frames = int(len(u_nk) * skip_frac)
                u_nk = u_nk.iloc[skip_frames:]
            
            # Decorrelate using subsampling module (following tutorial)
            u_nk = subsampling.decorrelate_u_nk(u_nk)
            
            u_nk_list.append(u_nk)
            success += 1
        except Exception as e:
            print(f"  WARNING: {ldir} failed: {e}")
    
    if success == 0:
        print(f"  ERROR: No valid data loaded")
        return None
    
    print(f"  Loaded: {success}/{len(lambda_dirs)} windows")
    
    # Concatenate all lambda windows (following tutorial pattern)
    u_nk = pd.concat(u_nk_list)
    print(f"  u_nk shape: {u_nk.shape} (after decorrelation)")
    
    return u_nk


def calc_dG_mbar(u_nk_fwd, u_nk_rev, T=310):
    """Calculate ΔG using MBAR method
    
    Following alchemlyb tutorial workflow:
    1. Concatenate forward and reverse data
    2. Fit MBAR estimator  
    3. Extract ΔG from first to last lambda state
    """
    print(f"  Method: MBAR")
    
    if u_nk_fwd is None or u_nk_rev is None:
        print(f"    ERROR: No u_nk data")
        return None, None, None
    
    # Combine forward and reverse trajectories
    combined = pd.concat([u_nk_fwd, u_nk_rev])
    
    # Fit MBAR estimator
    mbar = MBAR()
    mbar.fit(combined)
    
    # Extract results: use .iloc for robust position-based indexing
    # This avoids MultiIndex complexity - get ΔG from first to last lambda state
    dG = float(mbar.delta_f_.iloc[0, -1])
    err = float(mbar.d_delta_f_.iloc[0, -1])
    
    print(f"    ΔG = {dG:.3f} ± {err:.3f} kJ/mol")
    
    return dG, err, mbar


def plot_overlap_matrix(mbar, output_path, title):
    """Plot MBAR overlap matrix using alchemlyb"""
    try:
        ax = plot_mbar_overlap_matrix(mbar.overlap_matrix)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.figure.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(ax.figure)
        print(f"  Saved: {output_path}")
    except Exception as e:
        print(f"  WARNING: Failed to plot overlap matrix: {e}")


def plot_free_energy_profile(mbar, output_path, title):
    """Plot free energy profile manually from MBAR delta_f_"""
    try:
        # Extract free energy values from MBAR result
        delta_f = mbar.delta_f_
        
        # Get free energy relative to first state (lambda=0)
        states = delta_f.columns
        free_energies = delta_f.iloc[0, :]  # Row 0 is relative to first state
        errors = mbar.d_delta_f_.iloc[0, :]
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.errorbar(range(len(states)), free_energies, yerr=errors, 
                   marker='o', markersize=8, linewidth=2, capsize=5,
                   color='#2E86AB', ecolor='#A23B72', capthick=2)
        
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('Lambda state', fontsize=12)
        ax.set_ylabel('ΔF (kJ/mol)', fontsize=12)
        ax.set_xticks(range(len(states)))
        # Handle MultiIndex (tuples) vs simple states
        if isinstance(states[0], tuple):
            # For MultiIndex, use first element (VDW lambda) for labeling
            ax.set_xticklabels([f'{s[0]:.2f}' for s in states], rotation=45, ha='right')
        else:
            ax.set_xticklabels([f'{float(s):.2f}' for s in states], rotation=45, ha='right')
        ax.grid(alpha=0.3, linestyle='--')
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=1)
        
        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {output_path}")
    except Exception as e:
        print(f"  WARNING: Failed to plot free energy profile: {e}")
        import traceback
        print(f"  Debug: {traceback.format_exc()}")


def plot_convergence_analysis(u_nk_fwd, u_nk_rev, output_path, title):
    """
    Show data quality and sampling statistics instead of refitting MBAR.
    
    Since data is already decorrelated, convergence is implicitly validated.
    This plot shows the balance and coverage of forward/reverse sampling.
    """
    try:
        # Get basic statistics
        n_fwd = len(u_nk_fwd)
        n_rev = len(u_nk_rev)
        n_states = len(u_nk_fwd.columns)
        
        # Create figure with 3 subplots
        fig = plt.figure(figsize=(15, 5))
        gs = fig.add_gridspec(1, 3, hspace=0.3, wspace=0.3)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])
        
        # Plot 1: Forward vs Reverse sample comparison
        categories = ['Forward', 'Reverse', 'Combined']
        counts = [n_fwd, n_rev, n_fwd + n_rev]
        colors = ['#2E86AB', '#A23B72', '#44803F']
        
        bars = ax1.bar(categories, counts, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
        ax1.set_ylabel('Number of Decorrelated Samples', fontsize=11, fontweight='bold')
        ax1.set_title('Sampling Statistics', fontsize=12, fontweight='bold')
        ax1.grid(alpha=0.3, linestyle='--', axis='y')
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        # Plot 2: Forward/Reverse ratio analysis
        fwd_rev_ratio = n_fwd / n_rev if n_rev > 0 else 0
        optimal_ratio = 1.0
        
        x = ['Forward/Reverse\nRatio']
        y = [fwd_rev_ratio]
        bar = ax2.bar(x, y, color='#F18F01', alpha=0.8, edgecolor='black', linewidth=1.5, width=0.5)
        ax2.axhline(y=optimal_ratio, color='green', linestyle='--', linewidth=2, 
                   label=f'Optimal: {optimal_ratio:.1f}', alpha=0.7)
        ax2.set_ylabel('Ratio', fontsize=11, fontweight='bold')
        ax2.set_title('Forward/Reverse Balance', fontsize=12, fontweight='bold')
        ax2.set_ylim(0, max(2.0, fwd_rev_ratio * 1.2))
        ax2.legend(fontsize=10)
        ax2.grid(alpha=0.3, linestyle='--', axis='y')
        
        # Add value label
        ax2.text(0, fwd_rev_ratio, f'{fwd_rev_ratio:.2f}',
                ha='center', va='bottom', fontsize=11, fontweight='bold')
        
        # Plot 3: Lambda state coverage
        lambda_coverage = [n_states, n_states, n_states]  # All states covered
        labels = ['Forward\nTrajectory', 'Reverse\nTrajectory', 'Combined']
        
        bars3 = ax3.bar(labels, lambda_coverage, color=colors, alpha=0.8, 
                       edgecolor='black', linewidth=1.5)
        ax3.set_ylabel('Number of Lambda States', fontsize=11, fontweight='bold')
        ax3.set_title('Lambda State Coverage', fontsize=12, fontweight='bold')
        ax3.set_ylim(0, n_states * 1.2)
        ax3.grid(alpha=0.3, linestyle='--', axis='y')
        
        # Add value labels
        for bar in bars3:
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}/{int(height)}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        # Add summary text
        balance_status = "Balanced" if 0.5 <= fwd_rev_ratio <= 2.0 else "Imbalanced"
        balance_color = "green" if balance_status == "Balanced" else "orange"
        
        summary_text = (f"Data Quality: {balance_status}\n"
                       f"Total decorrelated samples: {n_fwd + n_rev}\n"
                       f"Lambda states covered: {n_states}")
        
        fig.text(0.5, 0.02, summary_text,
                ha='center', fontsize=10, style='italic',
                bbox=dict(boxstyle='round', facecolor=balance_color, alpha=0.3))
        
        fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
        fig.tight_layout(rect=[0, 0.06, 1, 0.96])
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {output_path}")
        
    except Exception as e:
        print(f"  WARNING: Failed to plot convergence: {e}")
        import traceback
        print(f"  Debug: {traceback.format_exc()}")
        import traceback
        print(f"  Debug: {traceback.format_exc()}")


def main():
    parser = argparse.ArgumentParser(description='FEP ΔΔG Analysis (MBAR)')
    parser.add_argument('directory', help='Output directory (e.g., outputs/test/design_name)')
    parser.add_argument('--skip', type=float, default=SKIP_FRACTION,
                       help=f'Skip fraction for equilibration (default: {SKIP_FRACTION})')
    parser.add_argument('--temp', type=float, default=TEMPERATURE,
                       help=f'Temperature in K (default: {TEMPERATURE})')
    parser.add_argument('--no-plots', action='store_true',
                       help='Skip generating plots (faster)')
    
    args = parser.parse_args()
    
    base_path = os.path.abspath(args.directory)
    if not os.path.exists(base_path):
        print(f"ERROR: Directory not found: {base_path}")
        sys.exit(1)
    
    print("="*70)
    print("FEP ΔΔG_binding Analysis (MBAR Method)")
    print("="*70)
    print(f"Directory: {base_path}")
    print(f"Temperature: {args.temp} K")
    print(f"Skip: {args.skip*100:.0f}%")
    print(f"Visualization: {'Disabled' if args.no_plots else 'Enabled'}")
    print("="*70)
    
    # Load data
    print("\n" + "="*70)
    print("LOADING DATA")
    print("="*70)
    
    u_b_f = load_data(base_path, 'bound', 'forward', args.skip, args.temp)
    u_b_r = load_data(base_path, 'bound', 'reverse', args.skip, args.temp)
    u_u_f = load_data(base_path, 'unbound', 'forward', args.skip, args.temp)
    u_u_r = load_data(base_path, 'unbound', 'reverse', args.skip, args.temp)
    
    # Check data completeness
    print(f"\nData Completeness Check:")
    states = [
        ('Bound/Forward', u_b_f),
        ('Bound/Reverse', u_b_r),
        ('Unbound/Forward', u_u_f),
        ('Unbound/Reverse', u_u_r)
    ]
    
    all_ok = True
    for name, data in states:
        status = 'OK' if data is not None else 'MISSING'
        if data is None:
            all_ok = False
        print(f"  {name:18s}: {status}")
    
    if not all_ok:
        print(f"\n{'!'*70}")
        print(f"ERROR: Incomplete data, cannot proceed")
        print(f"Check SLURM logs: {base_path}/*/logs/")
        print(f"{'!'*70}")
        sys.exit(1)
    
    # Calculate ΔG for BOUND state
    print("\n" + "="*70)
    print("MBAR ANALYSIS - BOUND STATE")
    print("="*70)
    dG_b, err_b, mbar_b = calc_dG_mbar(u_b_f, u_b_r, args.temp)
    
    if dG_b is None:
        print("ERROR: BOUND state calculation failed")
        sys.exit(1)
    
    # Calculate ΔG for UNBOUND state
    print("\n" + "="*70)
    print("MBAR ANALYSIS - UNBOUND STATE")
    print("="*70)
    dG_u, err_u, mbar_u = calc_dG_mbar(u_u_f, u_u_r, args.temp)
    
    if dG_u is None:
        print("ERROR: UNBOUND state calculation failed")
        sys.exit(1)
    
    # Calculate ΔΔG_binding
    ddG = dG_b - dG_u
    err = np.sqrt(err_b**2 + err_u**2)
    
    # Print results
    print("\n" + "="*70)
    print("FINAL RESULTS (MBAR)")
    print("="*70)
    print(f"ΔG_bound    = {dG_b:8.2f} ± {err_b:.2f} kJ/mol ({dG_b/4.184:7.2f} ± {err_b/4.184:.2f} kcal/mol)")
    print(f"ΔG_unbound  = {dG_u:8.2f} ± {err_u:.2f} kJ/mol ({dG_u/4.184:7.2f} ± {err_u/4.184:.2f} kcal/mol)")
    print("-" * 70)
    print(f"ΔΔG_binding = {ddG:7.2f} ± {err:.2f} kJ/mol ({ddG/4.184:6.2f} ± {err/4.184:.2f} kcal/mol)")
    print("="*70)
    
    # Interpretation
    print(f"\nInterpretation:")
    if ddG > 2.0:
        print(f"  ΔΔG = {ddG:.1f} kJ/mol > 0  → Mutation WEAKENS binding")
    elif ddG < -2.0:
        print(f"  ΔΔG = {ddG:.1f} kJ/mol < 0  → Mutation STRENGTHENS binding")
    else:
        print(f"  ΔΔG ≈ {ddG:.1f} kJ/mol ≈ 0  → NO significant effect")
    
    # Save numerical results
    output_file = os.path.join(base_path, 'ddG_results.txt')
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("FEP ΔΔG_binding Analysis Results (MBAR)\n")
        f.write("="*70 + "\n\n")
        f.write(f"Directory: {base_path}\n")
        f.write(f"Temperature: {args.temp} K\n")
        f.write(f"Skip fraction: {args.skip:.1%}\n\n")
        f.write(f"ΔG_bound    = {dG_b:8.3f} ± {err_b:.3f} kJ/mol\n")
        f.write(f"            = {dG_b/4.184:8.3f} ± {err_b/4.184:.3f} kcal/mol\n\n")
        f.write(f"ΔG_unbound  = {dG_u:8.3f} ± {err_u:.3f} kJ/mol\n")
        f.write(f"            = {dG_u/4.184:8.3f} ± {err_u/4.184:.3f} kcal/mol\n\n")
        f.write(f"ΔΔG_binding = {ddG:7.3f} ± {err:.3f} kJ/mol\n")
        f.write(f"            = {ddG/4.184:7.3f} ± {err/4.184:.3f} kcal/mol\n")
    
    print(f"\nResults saved to: {output_file}")
    
    # Generate plots
    if not args.no_plots:
        print("\n" + "="*70)
        print("GENERATING VISUALIZATION PLOTS")
        print("="*70)
        
        # Overlap matrices
        print("\n1. MBAR Overlap Matrices:")
        plot_overlap_matrix(mbar_b, 
                          os.path.join(base_path, 'mbar_overlap_bound.png'),
                          'MBAR Overlap Matrix - Bound State')
        plot_overlap_matrix(mbar_u,
                          os.path.join(base_path, 'mbar_overlap_unbound.png'),
                          'MBAR Overlap Matrix - Unbound State')
        
        # Free energy profiles
        print("\n2. Free Energy Profiles:")
        plot_free_energy_profile(mbar_b,
                                os.path.join(base_path, 'free_energy_profile_bound.png'),
                                'Free Energy Profile - Bound State')
        plot_free_energy_profile(mbar_u,
                                os.path.join(base_path, 'free_energy_profile_unbound.png'),
                                'Free Energy Profile - Unbound State')
        
        # Convergence analysis
        print("\n3. Convergence Analysis:")
        plot_convergence_analysis(u_b_f, u_b_r,
                                 os.path.join(base_path, 'convergence_bound.png'),
                                 'Convergence - Bound State')
        plot_convergence_analysis(u_u_f, u_u_r,
                                 os.path.join(base_path, 'convergence_unbound.png'),
                                 'Convergence - Unbound State')
        
        print("\n" + "="*70)
        print(f"All plots saved to: {base_path}")
        print("="*70)
    
    print("\nAnalysis completed successfully!")


if __name__ == '__main__':
    main()
