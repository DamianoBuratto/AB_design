#!/bin/bash
# Batch compute RMSD for all jobs

# Configuration
MUT_ONLY=false  # Set to true to skip wild-type (WT) jobs

echo "=========================================="
echo "Batch RMSD Computation"
if [ "$MUT_ONLY" = true ]; then
    echo "Mode: Mutant only (skipping WT jobs)"
else
    echo "Mode: All jobs (mut + wt)"
fi
echo "Start: $(date)"
echo "=========================================="
echo ""

success_count=0
fail_count=0
total_count=0

# Search ranking_scores.csv job directories
for csv_file in jobs/*/outputs/Antibody-pHLA_Complex/Antibody-pHLA_Complex_ranking_scores.csv; do
    if [ ! -f "$csv_file" ]; then
        continue
    fi
    
    # Extract job name
    job_dir=$(dirname $(dirname $(dirname "$csv_file")))
    job_name=$(basename "$job_dir")
    
    # Skip WT jobs if MUT_ONLY is enabled
    if [ "$MUT_ONLY" = true ] && [[ "$job_name" == *_wt ]]; then
        continue
    fi
    
    total_count=$((total_count + 1))
    
    echo "----------------------------------------"
    echo "[$total_count] Processing Job: $job_name"
    echo "----------------------------------------"
    
    # Find corresponding RF2 PDB file
    rf2_pdb=$(find "$job_dir" -maxdepth 1 -name "*_best.pdb" | head -1)
    
    if [ -z "$rf2_pdb" ]; then
        echo "  ERROR: RF2 PDB not found, skipping"
        fail_count=$((fail_count + 1))
        echo ""
        continue
    fi
    
    output_dir="$job_dir/outputs/Antibody-pHLA_Complex"
    
    echo "  RF2 PDB: $(basename "$rf2_pdb")"
    echo "  Output dir: $output_dir"
    
    # Check if antibody_rmsd column exists
    if head -1 "$csv_file" | grep -q "antibody_rmsd"; then
        echo "  WARNING:  RMSD data exists, recalculating..."
    fi
    
    # Run RMSD computation
    python3 compute_rmsd.py "$output_dir" "$rf2_pdb" "$csv_file" 2>&1 | grep -E "Found|Processing|Successfully|RMSD range|Error" | sed 's/^/    /'
    
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "  OK: Success"
        success_count=$((success_count + 1))
    else
        echo "  ERROR: Failed"
        fail_count=$((fail_count + 1))
    fi
    
    echo ""
done

echo "=========================================="
echo "Batch Computation Complete"
echo "Total: $total_count jobs"
echo "Success: $success_count"
echo "Failed: $fail_count"
echo "End: $(date)"
echo "=========================================="
echo ""

# Merge all CSVs and calculate delta metrics
echo "=========================================="
echo "Analyzing All Results with Delta Metrics..."
echo "=========================================="

python3 << 'PYTHON_EOF'
import pandas as pd
import glob
import os
import re

def extract_job_info(job_dir_name):
    """
    Extract job name and type from directory name.
    Examples: pHLA_4_1 -> (pHLA_4_1, mut)
              pHLA_4_1_wt -> (pHLA_4_1, wt)
    """
    if job_dir_name.endswith('_wt'):
        job_name = job_dir_name[:-3]
        job_type = 'wt'
    else:
        job_name = job_dir_name
        job_type = 'mut'
    return job_name, job_type

def calculate_deltas_and_stats(df_group):
    """
    Calculate metrics for antibody and CDR for a job+type group.
    Antibody metrics displayed in the row with minimum antibody_rmsd.
    CDR metrics displayed in the row with minimum cdr_rmsd (independent).
    
    Antibody metrics: delta_score, delta_rmsd, std_score, std_rmsd
    CDR metrics: cdr_delta_score, cdr_delta_rmsd, cdr_std_score, cdr_std_rmsd
    """
    # Find minimum antibody RMSD row
    min_ab_rmsd_idx = df_group['antibody_rmsd'].idxmin()
    max_ab_rmsd = df_group['antibody_rmsd'].max()
    max_score = df_group['ranking_score'].max()
    min_ab_rmsd_value = df_group.loc[min_ab_rmsd_idx, 'antibody_rmsd']
    score_at_min_ab_rmsd = df_group.loc[min_ab_rmsd_idx, 'ranking_score']
    
    # Antibody Delta metrics
    delta_rmsd = min_ab_rmsd_value - max_ab_rmsd
    delta_score = score_at_min_ab_rmsd - max_score
    
    # Antibody Average
    avg_rmsd = df_group['antibody_rmsd'].mean()
    avg_score = df_group['ranking_score'].mean()
    
    # Antibody Standard deviation (population std using 1/n)
    std_rmsd = df_group['antibody_rmsd'].std(ddof=0)  # 1/n
    std_score = df_group['ranking_score'].std(ddof=0)  # 1/n
    
    # CDR metrics (if cdr_rmsd column exists)
    if 'cdr_rmsd' in df_group.columns:
        # Find minimum CDR RMSD row (independent from antibody)
        min_cdr_rmsd_idx = df_group['cdr_rmsd'].idxmin()
        max_cdr_rmsd = df_group['cdr_rmsd'].max()
        min_cdr_rmsd_value = df_group.loc[min_cdr_rmsd_idx, 'cdr_rmsd']
        score_at_min_cdr_rmsd = df_group.loc[min_cdr_rmsd_idx, 'ranking_score']
        
        cdr_delta_rmsd = min_cdr_rmsd_value - max_cdr_rmsd
        cdr_delta_score = score_at_min_cdr_rmsd - max_score
        
        # CDR Average
        avg_cdr_rmsd = df_group['cdr_rmsd'].mean()
        
        # CDR Standard deviation (population std using 1/n)
        cdr_std_rmsd = df_group['cdr_rmsd'].std(ddof=0)  # 1/n
        cdr_std_score = std_score  # Same as antibody
    else:
        min_cdr_rmsd_idx = None
        avg_cdr_rmsd = pd.NA
        cdr_delta_rmsd = pd.NA
        cdr_delta_score = pd.NA
        cdr_std_rmsd = pd.NA
        cdr_std_score = pd.NA
    
    # Initialize all metric columns to NA
    df_group['avg_rmsd'] = pd.NA
    df_group['delta_score'] = pd.NA
    df_group['delta_rmsd'] = pd.NA
    df_group['std_score'] = pd.NA
    df_group['std_rmsd'] = pd.NA
    df_group['avg_cdr_rmsd'] = pd.NA
    df_group['cdr_delta_score'] = pd.NA
    df_group['cdr_delta_rmsd'] = pd.NA
    df_group['cdr_std_score'] = pd.NA
    df_group['cdr_std_rmsd'] = pd.NA
    
    # Assign antibody metrics to minimum antibody RMSD row
    df_group.loc[min_ab_rmsd_idx, 'avg_rmsd'] = avg_rmsd
    df_group.loc[min_ab_rmsd_idx, 'delta_score'] = delta_score
    df_group.loc[min_ab_rmsd_idx, 'delta_rmsd'] = delta_rmsd
    df_group.loc[min_ab_rmsd_idx, 'std_score'] = std_score
    df_group.loc[min_ab_rmsd_idx, 'std_rmsd'] = std_rmsd
    
    # Assign CDR metrics to minimum CDR RMSD row (independent)
    if min_cdr_rmsd_idx is not None:
        df_group.loc[min_cdr_rmsd_idx, 'avg_cdr_rmsd'] = avg_cdr_rmsd
        df_group.loc[min_cdr_rmsd_idx, 'cdr_delta_score'] = cdr_delta_score
        df_group.loc[min_cdr_rmsd_idx, 'cdr_delta_rmsd'] = cdr_delta_rmsd
        df_group.loc[min_cdr_rmsd_idx, 'cdr_std_score'] = cdr_std_score
        df_group.loc[min_cdr_rmsd_idx, 'cdr_std_rmsd'] = cdr_std_rmsd
    
    return df_group

# Collect all CSV files
all_data = []
csv_files = sorted(glob.glob("jobs/*/outputs/Antibody-pHLA_Complex/Antibody-pHLA_Complex_ranking_scores.csv"))

print(f"Found {len(csv_files)} CSV files\n")

for csv_file in csv_files:
    job_dir = csv_file.split('/')[1]
    job_name, job_type = extract_job_info(job_dir)
    
    try:
        df = pd.read_csv(csv_file)
        
        if 'antibody_rmsd' not in df.columns:
            print(f"  WARNING:  {job_dir}: Missing antibody_rmsd column, skipping")
            continue
        
        # Check if seed and sample columns already exist (from compute_rmsd.py)
        if 'seed' not in df.columns or 'sample' not in df.columns:
            # Try to extract from seed_sample or model_name columns
            if 'seed_sample' in df.columns:
                df[['seed', 'sample']] = df['seed_sample'].str.extract(r'seed-(\d+)_sample-(\d+)')
                df['seed'] = df['seed'].astype(int)
                df['sample'] = df['sample'].astype(int)
            elif 'model_name' in df.columns:
                df[['seed', 'sample']] = df['model_name'].str.extract(r'seed-(\d+)_sample-(\d+)')
                df['seed'] = df['seed'].astype(int)
                df['sample'] = df['sample'].astype(int)
            else:
                print(f"  WARNING:  {job_dir}: No seed/sample columns found, using indices")
                df['seed'] = range(len(df))
                df['sample'] = 0
        # If columns exist, ensure they are integers
        else:
            df['seed'] = df['seed'].astype(int)
            df['sample'] = df['sample'].astype(int)
        
        # Add job info
        df.insert(0, 'job', job_name)
        df.insert(1, 'type', job_type)
        
        # Check for abnormal cdr_rmsd values (should not exceed antibody_rmsd significantly)
        if 'cdr_rmsd' in df.columns:
            abnormal_cdr = df[df['cdr_rmsd'] > df['antibody_rmsd'] * 1.2]
            if len(abnormal_cdr) > 0:
                print(f"  WARNING:  {job_dir}: {len(abnormal_cdr)} rows with abnormally high cdr_rmsd (>1.2x antibody_rmsd)")
                print(f"      This may indicate CDR atom count mismatch between RF2 and AF3")
        
        all_data.append(df)
        print(f"  OK: {job_dir} -> {job_name} ({job_type}): {len(df)} rows")
        
    except Exception as e:
        print(f"  ERROR: {job_dir}: Failed - {e}")
        continue

if not all_data:
    print("\nERROR: No data to process!")
    exit(1)

# Combine all data
print(f"\nCombining all data...")
combined_df = pd.concat(all_data, ignore_index=True)
print(f"Total rows before delta calculation: {len(combined_df)}")

# Calculate deltas and stats for each job+type group
print(f"Calculating antibody and CDR metrics (delta and std)...\n")
result_df = combined_df.groupby(['job', 'type'], group_keys=False).apply(calculate_deltas_and_stats)

# Sort by job, type (wt before mut), then RMSD
result_df = result_df.sort_values(
    ['job', 'type', 'antibody_rmsd'],
    ascending=[True, False, True]
).reset_index(drop=True)

# Reorder columns
column_order = ['job', 'type', 'seed', 'sample', 'ranking_score', 'antibody_rmsd', 'cdr_rmsd',
                'avg_rmsd', 'delta_score', 'delta_rmsd', 'std_score', 'std_rmsd',
                'avg_cdr_rmsd', 'cdr_delta_score', 'cdr_delta_rmsd', 'cdr_std_score', 'cdr_std_rmsd']
other_cols = [col for col in result_df.columns if col not in column_order]
final_columns = [col for col in column_order + other_cols if col in result_df.columns]
result_df = result_df[final_columns]

# Save final analysis
output_file = "af3_results_analysis.csv"
result_df.to_csv(output_file, index=False)

print("=" * 80)
print("Summary Statistics")
print("=" * 80)
print(f"Total rows: {len(result_df)}")
print(f"Unique jobs: {result_df['job'].nunique()}")
print(f"Jobs with mut: {result_df[result_df['type']=='mut']['job'].nunique()}")
print(f"Jobs with wt: {result_df[result_df['type']=='wt']['job'].nunique()}")

delta_rows = result_df[result_df['delta_rmsd'].notna()]
print(f"\nRows with delta/std values: {len(delta_rows)}")
if len(delta_rows) > 0:
    print(f"\nAntibody metrics:")
    print(f"  Delta RMSD range: {delta_rows['delta_rmsd'].min():.2f} to {delta_rows['delta_rmsd'].max():.2f}")
    print(f"  Delta score range: {delta_rows['delta_score'].min():.4f} to {delta_rows['delta_score'].max():.4f}")
    print(f"  Std RMSD range: {delta_rows['std_rmsd'].min():.2f} to {delta_rows['std_rmsd'].max():.2f}")
    print(f"  Std score range: {delta_rows['std_score'].min():.4f} to {delta_rows['std_score'].max():.4f}")
    
    # CDR metrics (if available)
    if 'cdr_delta_rmsd' in delta_rows.columns:
        cdr_rows = delta_rows[delta_rows['cdr_delta_rmsd'].notna()]
        if len(cdr_rows) > 0:
            print(f"\nCDR metrics:")
            print(f"  CDR Delta RMSD range: {cdr_rows['cdr_delta_rmsd'].min():.2f} to {cdr_rows['cdr_delta_rmsd'].max():.2f}")
            print(f"  CDR Delta score range: {cdr_rows['cdr_delta_score'].min():.4f} to {cdr_rows['cdr_delta_score'].max():.4f}")
            print(f"  CDR Std RMSD range: {cdr_rows['cdr_std_rmsd'].min():.2f} to {cdr_rows['cdr_std_rmsd'].max():.2f}")
            print(f"  CDR Std score range: {cdr_rows['cdr_std_score'].min():.4f} to {cdr_rows['cdr_std_score'].max():.4f}")

print(f"\nOK: Analysis complete!")
print(f"📊 Results saved to: {output_file}")

# Show preview
print(f"\n{'-'*80}")
print("Preview (first 20 rows):")
print(f"{'-'*80}")
print(result_df.head(20).to_string(index=False, max_colwidth=15))

# Show best results per job
print(f"\n{'-'*80}")
print("🏆 Best Results per Job (Top 10, by antibody RMSD):")
print(f"{'-'*80}")
best_per_job_type = result_df.loc[result_df.groupby(['job', 'type'])['antibody_rmsd'].idxmin()]
best_per_job_type = best_per_job_type.sort_values('antibody_rmsd').head(10)
display_cols = ['job', 'type', 'seed', 'sample', 'ranking_score', 'antibody_rmsd', 'cdr_rmsd',
                'avg_rmsd', 'delta_score', 'delta_rmsd', 'std_score', 'std_rmsd']
# Add CDR columns if they exist
if 'cdr_delta_score' in best_per_job_type.columns:
    display_cols += ['avg_cdr_rmsd', 'cdr_delta_score', 'cdr_delta_rmsd', 'cdr_std_score', 'cdr_std_rmsd']
display_cols = [col for col in display_cols if col in best_per_job_type.columns]
print(best_per_job_type[display_cols].to_string(index=False))

PYTHON_EOF

echo ""
echo "=========================================="
echo "All Done!"
echo "=========================================="
