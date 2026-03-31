#!/usr/bin/env python3
"""
Filter RF2 PDB analysis results and copy selected PDB files to AF3/RFantibody_best.

This script:
1. Copies rf2_pdb_analysis.csv from a specified RFantibody output directory
2. Filters rows meeting criteria: Framework_CDR_RMSD < 2 AND sc_<=5A >= 4
3. Copies corresponding PDB files from rf2/ subdirectory to AF3/RFantibody_best/

Usage:
    python3 filter_and_copy_rf2_best.py --rf-output-dir /path/to/RFantibody/output

Example:
    python3 filter_and_copy_rf2_best.py \
        --rf-output-dir /public/home/xuziyi/RFantibody/outputs/200T_700_32_L1_10-13_H3_11-14_T82-159_replica2
"""

import os
import sys
import argparse
import shutil
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description='Filter RF2 analysis results and copy best PDBs to AF3/RFantibody_best'
    )
    parser.add_argument(
        '--rf-output-dir',
        required=True,
        help='Path to RFantibody output directory (e.g., /public/home/xuziyi/RFantibody/outputs/200T_700_32_L1_10-13_H3_11-14_T82-159_replica2)'
    )
    parser.add_argument(
        '--framework-cdr-rmsd-cutoff',
        type=float,
        default=2.0,
        help='Maximum Framework_CDR_RMSD (default: 2.0)'
    )
    parser.add_argument(
        '--sc-5a-cutoff',
        type=int,
        default=4,
        help='Minimum sc_<=5A (default: 4)'
    )
    parser.add_argument(
        '--output-dir',
        default='RFantibody_best',
        help='Output directory relative to AF3/ (default: RFantibody_best)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print actions without copying files'
    )
    parser.add_argument(
        '--debug-csv',
        action='store_true',
        help='Print first few lines of CSV for debugging'
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Validate RF output directory
    rf_output_dir = os.path.abspath(args.rf_output_dir)
    if not os.path.isdir(rf_output_dir):
        print(f"❌ Error: RF output directory not found: {rf_output_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Locate rf2_pdb_analysis.csv
    csv_path = os.path.join(rf_output_dir, 'rf2_pdb_analysis.csv')
    if not os.path.isfile(csv_path):
        print(f"❌ Error: rf2_pdb_analysis.csv not found at: {csv_path}", file=sys.stderr)
        sys.exit(1)
    
    # Locate rf2/ subdirectory
    rf2_dir = os.path.join(rf_output_dir, 'rf2')
    if not os.path.isdir(rf2_dir):
        print(f"❌ Error: rf2/ subdirectory not found at: {rf2_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Setup output directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, args.output_dir)
    
    if not args.dry_run:
        os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 80)
    print("RF2 Best PDB Filter and Copy")
    print("=" * 80)
    print(f"RF output directory: {rf_output_dir}")
    print(f"CSV file: {csv_path}")
    print(f"RF2 PDB directory: {rf2_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Framework_CDR_RMSD cutoff: < {args.framework_cdr_rmsd_cutoff}")
    print(f"sc_<=5A cutoff: >= {args.sc_5a_cutoff}")
    print()
    
    # Debug mode: print raw CSV lines
    if args.debug_csv:
        print("=" * 80)
        print("DEBUG: First 20 lines of CSV file:")
        print("=" * 80)
        try:
            with open(csv_path, 'r', encoding='utf-8') as f:
                for i, line in enumerate(f):
                    if i < 20:
                        print(f"Line {i+1}: {repr(line[:200])}")  # Show first 200 chars
                    else:
                        break
        except Exception as e:
            print(f"Error reading file: {e}")
        print("=" * 80)
        print()
    
    # Step 1: Copy CSV file
    csv_dest = os.path.join(output_dir, 'rf2_pdb_analysis.csv')
    if not args.dry_run:
        shutil.copy2(csv_path, csv_dest)
        print(f"✅ Copied CSV to: {csv_dest}")
    else:
        print(f"[DRY RUN] Would copy CSV to: {csv_dest}")
    
    # Step 2: Load and filter CSV
    print()
    print("Loading and filtering CSV...")
    try:
        # Try reading with different parameters to handle edge cases
        df = pd.read_csv(csv_path, encoding='utf-8', on_bad_lines='skip', engine='python')
        print(f"Total rows in CSV: {len(df)}")
    except Exception as e:
        print(f"❌ Error reading CSV with utf-8: {e}", file=sys.stderr)
        # Try alternative encoding
        try:
            df = pd.read_csv(csv_path, encoding='latin1', on_bad_lines='skip', engine='python')
            print(f"Total rows in CSV (with latin1 encoding): {len(df)}")
        except Exception as e2:
            print(f"❌ Error reading CSV with latin1: {e2}", file=sys.stderr)
            # Last resort: try with different error handling
            try:
                df = pd.read_csv(csv_path, encoding='utf-8', error_bad_lines=False, warn_bad_lines=True)
                print(f"Total rows in CSV (skipped bad lines): {len(df)}")
            except:
                print(f"❌ Failed to read CSV with all methods", file=sys.stderr)
                sys.exit(1)
    
    # Check required columns
    required_cols = ['Filename', 'Framework_CDR_RMSD', 'sc_<=5A']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"❌ Error: Missing required columns: {missing_cols}", file=sys.stderr)
        print(f"Available columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)
    
    # Apply filters
    filtered_df = df[
        (df['Framework_CDR_RMSD'] < args.framework_cdr_rmsd_cutoff) &
        (df['sc_<=5A'] >= args.sc_5a_cutoff)
    ]
    
    print(f"Rows after filtering: {len(filtered_df)}")
    print()
    
    if len(filtered_df) == 0:
        print("⚠️  Warning: No rows meet the filtering criteria")
        print("Exiting without copying any PDB files")
        sys.exit(0)
    
    # Show statistics
    print("Filtering Statistics:")
    print(f"  Framework_CDR_RMSD range: {filtered_df['Framework_CDR_RMSD'].min():.2f} - {filtered_df['Framework_CDR_RMSD'].max():.2f}")
    print(f"  sc_<=5A range: {filtered_df['sc_<=5A'].min()} - {filtered_df['sc_<=5A'].max()}")
    print()
    
    # Step 3: Copy PDB files
    print("Copying PDB files...")
    print("-" * 80)
    
    copied_count = 0
    missing_count = 0
    failed_count = 0
    
    for idx, row in filtered_df.iterrows():
        filename = row['Filename']
        
        # Construct source and destination paths
        # Filename format: design_488_dldesign_10 -> design_488_dldesign_10_best.pdb
        pdb_basename = f"{filename}_best.pdb"
        src_pdb = os.path.join(rf2_dir, pdb_basename)
        dest_pdb = os.path.join(output_dir, pdb_basename)
        
        # Check if source exists
        if not os.path.isfile(src_pdb):
            print(f"  ⚠️  [{idx+1}/{len(filtered_df)}] Missing: {pdb_basename}")
            missing_count += 1
            continue
        
        # Copy file
        if not args.dry_run:
            try:
                shutil.copy2(src_pdb, dest_pdb)
                copied_count += 1
                if copied_count <= 10 or copied_count % 50 == 0:
                    print(f"  ✅ [{copied_count}/{len(filtered_df)}] Copied: {pdb_basename}")
            except Exception as e:
                print(f"  ❌ [{idx+1}/{len(filtered_df)}] Failed to copy {pdb_basename}: {e}")
                failed_count += 1
        else:
            print(f"  [DRY RUN] [{idx+1}/{len(filtered_df)}] Would copy: {pdb_basename}")
            copied_count += 1
    
    # Summary
    print()
    print("=" * 80)
    print("Summary:")
    print("=" * 80)
    print(f"Total filtered rows: {len(filtered_df)}")
    print(f"Successfully copied: {copied_count}")
    print(f"Missing PDB files: {missing_count}")
    print(f"Failed to copy: {failed_count}")
    print()
    
    if not args.dry_run:
        print(f"✅ PDB files copied to: {output_dir}")
        print()
        print("Next steps:")
        print("  1. Review the copied PDB files")
        print("  2. Run batch_prepare_af3.py to prepare AF3 jobs")
        print("  3. Submit jobs to HPC using alphafold3.sh")
    else:
        print("[DRY RUN] No files were actually copied")
    
    print()


if __name__ == '__main__':
    main()
