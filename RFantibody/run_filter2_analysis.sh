#!/bin/bash
# ==============================================================================
# Standalone Filter2 Sidechain Contact Analysis
# ==============================================================================
# Purpose: Analyze sidechain contacts between ProteinMPNN outputs and HIS hotspot
# Usage:   bash run_filter2_analysis.sh <proteinmpnn_dir> <output_dir>
#
# Example:
#   bash run_filter2_analysis.sh \
#       outputs/replica3/proteinmpnn \
#       outputs/replica3
#
# Output:
#   - filter2_sidechain_contacts.csv (contact statistics)
#   - filter2_analysis_summary.txt (statistical summary)
# ==============================================================================

set -e

# ==============================================================================
# Input validation
# ==============================================================================
if [ "$#" -ne 2 ]; then
    echo "ERROR: Incorrect number of arguments"
    echo ""
    echo "Usage: bash run_filter2_analysis.sh <proteinmpnn_dir> <output_dir>"
    echo ""
    echo "Arguments:"
    echo "  proteinmpnn_dir  : Directory containing ProteinMPNN PDB files"
    echo "  output_dir       : Directory for output CSV and summary"
    echo ""
    echo "Example:"
    echo "  bash run_filter2_analysis.sh outputs/replica3/proteinmpnn outputs/replica3"
    exit 1
fi

PROTEINMPNN_DIR="$1"
OUTPUT_DIR="$2"

# Validate inputs
if [ ! -d "$PROTEINMPNN_DIR" ]; then
    echo "ERROR: ProteinMPNN directory does not exist: $PROTEINMPNN_DIR"
    exit 1
fi

pdb_count=$(ls "$PROTEINMPNN_DIR"/*.pdb 2>/dev/null | wc -l)
if [ "$pdb_count" -eq 0 ]; then
    echo "ERROR: No PDB files found in $PROTEINMPNN_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# ==============================================================================
# Configuration
# ==============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${SCRIPT_DIR}"

# Load config if available (for GLIBC and Python setup)
CONFIG_FILE="${SCRIPT_DIR}/run_config.txt"
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    echo "WARNING: Config file not found, using defaults"
    PROJECT_DIR="$(pwd)"
fi

# ==============================================================================
# Environment setup
# ==============================================================================
# Try to setup GLIBC wrapper if available
if [ -n "${RFANTIBODY_PYTHON_BIN:-}" ]; then
    PYTHON_BIN="${RFANTIBODY_PYTHON_BIN}"
    USE_GLIBC_WRAPPER=true
else
    PYTHON_BIN="python"
    USE_GLIBC_WRAPPER=false
    echo "INFO: Using system Python (no GLIBC wrapper configured)"
fi

# Define run_with_glibc wrapper
if [ "$USE_GLIBC_WRAPPER" = true ] && [ -n "${GLIBC_LOADER:-}" ]; then
    run_with_glibc() {
        if [ "$#" -eq 0 ]; then
            echo "run_with_glibc: missing command" >&2
            return 1
        fi
        local cmd="$1"
        shift
        if [[ "${cmd}" != /* ]]; then
            cmd=$(command -v "${cmd}")
        fi
        "${GLIBC_LOADER}" --library-path "${RFANTIBODY_GLIBC_LIBRARY_PATH}" "${cmd}" "$@"
    }
else
    run_with_glibc() {
        "$@"
    }
fi

# ==============================================================================
# Run Filter2 analysis
# ==============================================================================
echo "=========================================="
echo "Filter2 Contact Analysis"
echo "=========================================="
echo "ProteinMPNN directory: $PROTEINMPNN_DIR"
echo "Output directory:      $OUTPUT_DIR"
echo "PDB files to analyze:  $pdb_count"
echo ""

echo "Step 1/3: Sidechain contact analysis..."
run_with_glibc "${PYTHON_BIN}" "${PROJECT_DIR}/count_contacts_with_hotspot.py" \
    --pdb-dir "$PROTEINMPNN_DIR" \
    --workdir "${OUTPUT_DIR}/.tmp_sc" \
    --mode sidechain \
    --cutoffs 5 6 > /dev/null

if [ $? -ne 0 ]; then
    echo "ERROR: Sidechain analysis failed"
    exit 1
fi

echo "Step 2/3: Backbone contact analysis..."
run_with_glibc "${PYTHON_BIN}" "${PROJECT_DIR}/count_contacts_with_hotspot.py" \
    --pdb-dir "$PROTEINMPNN_DIR" \
    --workdir "${OUTPUT_DIR}/.tmp_bb" \
    --mode backbone \
    --cutoffs 5 6 > /dev/null

if [ $? -ne 0 ]; then
    echo "ERROR: Backbone analysis failed"
    exit 1
fi

echo "Step 3/3: Merging results..."
python3 << PYEOF
import csv
import sys

sc_file = "${OUTPUT_DIR}/.tmp_sc/contacts_with_hotspot_summary.csv"
bb_file = "${OUTPUT_DIR}/.tmp_bb/contacts_with_hotspot_summary.csv"
output_file = "${OUTPUT_DIR}/filter2_sidechain_contacts.csv"

# Read both CSVs
with open(sc_file) as f:
    sc_data = {row['PDB File']: row for row in csv.DictReader(f)}

with open(bb_file) as f:
    bb_data = {row['PDB File']: row for row in csv.DictReader(f)}

# Merge data
merged = []
for pdb_name in sorted(sc_data.keys()):
    merged_row = {
        'PDB File': pdb_name,
        'sc_<=5A': sc_data[pdb_name]['sc_<=5A'],
        'bb_<=5A': bb_data[pdb_name]['bb_<=5A'],
        'sc_<=6A': sc_data[pdb_name]['sc_<=6A'],
        'bb_<=6A': bb_data[pdb_name]['bb_<=6A'],
    }
    merged.append(merged_row)

# Write merged CSV
with open(output_file, 'w', newline='') as f:
    fieldnames = ['PDB File', 'sc_<=5A', 'bb_<=5A', 'sc_<=6A', 'bb_<=6A']
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(merged)

print(f"Merged {len(merged)} designs")
PYEOF

# Cleanup temp directories
rm -rf "${OUTPUT_DIR}/.tmp_sc" "${OUTPUT_DIR}/.tmp_bb"

echo "✓ Filter2 analysis complete"
echo ""

# ==============================================================================
# Generate statistical summary
# ==============================================================================
CSV_FILE="$OUTPUT_DIR/filter2_sidechain_contacts.csv"
SUMMARY_FILE="$OUTPUT_DIR/filter2_analysis_summary.txt"

if [ ! -f "$CSV_FILE" ]; then
    echo "ERROR: Output CSV not found: $CSV_FILE"
    exit 1
fi

echo "Generating statistical summary..."

run_with_glibc "${PYTHON_BIN}" -c "
import csv
import sys
from collections import Counter

csv_file = '${CSV_FILE}'

# Read data
data_5a = []
data_6a = []

with open(csv_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        val_5a = row.get('(<= 5A)', '0')
        val_6a = row.get('(<= 6A)', '0')
        
        # Handle non-numeric values
        try:
            data_5a.append(int(val_5a))
        except ValueError:
            continue
        
        try:
            data_6a.append(int(val_6a))
        except ValueError:
            continue

if not data_5a or not data_6a:
    print('ERROR: No valid contact data found', file=sys.stderr)
    sys.exit(1)

# Calculate statistics
def get_stats(data):
    data_sorted = sorted(data)
    n = len(data_sorted)
    
    return {
        'count': n,
        'min': min(data),
        'max': max(data),
        'mean': sum(data) / n,
        'median': data_sorted[n//2] if n % 2 == 1 else (data_sorted[n//2-1] + data_sorted[n//2]) / 2,
        'q1': data_sorted[n//4],
        'q3': data_sorted[3*n//4],
    }

stats_5a = get_stats(data_5a)
stats_6a = get_stats(data_6a)

# Count distributions
count_5a_dist = Counter(data_5a)
count_6a_dist = Counter(data_6a)

# Write summary
summary_lines = []
summary_lines.append('='*70)
summary_lines.append('Filter2 Sidechain Contact Analysis Summary')
summary_lines.append('='*70)
summary_lines.append('')
summary_lines.append(f'Total sequences analyzed: {stats_5a[\"count\"]}')
summary_lines.append('')
summary_lines.append('-'*70)
summary_lines.append('Contact Statistics (<= 5A)')
summary_lines.append('-'*70)
summary_lines.append(f'  Minimum:        {stats_5a[\"min\"]:6.0f}')
summary_lines.append(f'  Maximum:        {stats_5a[\"max\"]:6.0f}')
summary_lines.append(f'  Mean:           {stats_5a[\"mean\"]:6.2f}')
summary_lines.append(f'  Median:         {stats_5a[\"median\"]:6.0f}')
summary_lines.append(f'  Q1 (25%):       {stats_5a[\"q1\"]:6.0f}')
summary_lines.append(f'  Q3 (75%):       {stats_5a[\"q3\"]:6.0f}')
summary_lines.append('')
summary_lines.append('-'*70)
summary_lines.append('Contact Statistics (<= 6A)')
summary_lines.append('-'*70)
summary_lines.append(f'  Minimum:        {stats_6a[\"min\"]:6.0f}')
summary_lines.append(f'  Maximum:        {stats_6a[\"max\"]:6.0f}')
summary_lines.append(f'  Mean:           {stats_6a[\"mean\"]:6.2f}')
summary_lines.append(f'  Median:         {stats_6a[\"median\"]:6.0f}')
summary_lines.append(f'  Q1 (25%):       {stats_6a[\"q1\"]:6.0f}')
summary_lines.append(f'  Q3 (75%):       {stats_6a[\"q3\"]:6.0f}')
summary_lines.append('')
summary_lines.append('-'*70)
summary_lines.append('Distribution (<= 6A contacts)')
summary_lines.append('-'*70)

# Show top 10 most common contact counts
for count, freq in sorted(count_6a_dist.items(), key=lambda x: -x[1])[:10]:
    pct = 100.0 * freq / stats_6a['count']
    summary_lines.append(f'  {count:3d} contacts: {freq:5d} sequences ({pct:5.1f}%)')

summary_lines.append('')
summary_lines.append('-'*70)
summary_lines.append('Filter Criteria Suggestions')
summary_lines.append('-'*70)
summary_lines.append('')
summary_lines.append('Based on the distribution, consider these filtering thresholds:')
summary_lines.append('')
summary_lines.append(f'1. Conservative (top 25%):  >= {stats_6a[\"q3\"]:0.0f} contacts at 6A')
summary_lines.append(f'2. Moderate (top 50%):      >= {stats_6a[\"median\"]:0.0f} contacts at 6A')
summary_lines.append(f'3. Permissive (mean):       >= {stats_6a[\"mean\"]:0.0f} contacts at 6A')
summary_lines.append('')
summary_lines.append('You can also combine with RF2 quality metrics:')
summary_lines.append('  - Example: (<= 6A) >= threshold AND pLDDT > 80')
summary_lines.append('')
summary_lines.append('='*70)

# Print to screen and file
output_text = '\\n'.join(summary_lines)
print(output_text)

with open('${SUMMARY_FILE}', 'w') as f:
    f.write(output_text + '\\n')

print('')
print(f'Summary saved to: ${SUMMARY_FILE}')
"

STATS_EXIT=$?

if [ $STATS_EXIT -ne 0 ]; then
    echo ""
    echo "WARNING: Statistical analysis failed"
else
    echo ""
    echo "=========================================="
    echo "Output Files"
    echo "=========================================="
    echo "Contact data:       $CSV_FILE"
    echo "Statistical summary: $SUMMARY_FILE"
    echo ""
fi

echo "=========================================="
echo "Analysis Complete"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Review $SUMMARY_FILE for statistics"
echo "2. Review $CSV_FILE for individual sequence contacts"
echo "3. Decide on filtering threshold based on distribution"
echo "4. Optionally combine with RF2 quality metrics for final selection"
echo ""
