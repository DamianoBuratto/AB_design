#!/bin/bash
# ============================================================================
# Continue RFantibody Pipeline from Filtering Stage
# ============================================================================
# Purpose: Resume pipeline after RFdiffusion completed but job was killed
# Runs: Filtering → ProteinMPNN → RF2 → Analysis
#
# Usage:
#   bash continue_from_filtering.sh <output_dir>
#   
# Example:
#   bash continue_from_filtering.sh ${PROJECT_DIR}/outputs/200T_700_32_L1_10-13_H3_11-14_T82-159_replica3
# ============================================================================

set -e

# ============================================================================
# Configuration
# ============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/run_config.txt"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

echo "Loading configuration from: $CONFIG_FILE"
source "$CONFIG_FILE"

# Get output directory from command line
if [ -z "$1" ]; then
    echo "ERROR: Please provide output directory"
    echo "Usage: bash continue_from_filtering.sh <output_dir>"
    echo ""
    echo "Example:"
    echo "  bash continue_from_filtering.sh \${PROJECT_DIR}/outputs/200T_700_32_L1_10-13_H3_11-14_T82-159_replica3"
    exit 1
fi

CONFIG_OUTPUT_DIR="$1"

if [ ! -d "$CONFIG_OUTPUT_DIR" ]; then
    echo "ERROR: Output directory does not exist: $CONFIG_OUTPUT_DIR"
    exit 1
fi

if [ ! -d "$CONFIG_OUTPUT_DIR/rfdiffusion" ]; then
    echo "ERROR: RFdiffusion directory not found: $CONFIG_OUTPUT_DIR/rfdiffusion"
    echo "Make sure RFdiffusion has completed and generated designs"
    exit 1
fi

# Extract configuration name from path
CONFIG_NAME=$(basename "$CONFIG_OUTPUT_DIR")

# Create batch ID for this submission
BATCH_ID=$(date +%Y%m%d_%H%M%S)
BATCH_LOG_DIR="${PROJECT_DIR}/logs/continue_filtering_${BATCH_ID}"
mkdir -p "$BATCH_LOG_DIR"

echo "=========================================="
echo "Continue RFantibody Pipeline"
echo "From: Filtering Stage"
echo "=========================================="
echo "Configuration: ${CONFIG_NAME}"
echo "Output Directory: ${CONFIG_OUTPUT_DIR}"
echo "Batch ID: ${BATCH_ID}"
echo "Log Directory: ${BATCH_LOG_DIR}"
echo ""

# Verify RFdiffusion outputs exist
design_count=$(ls ${CONFIG_OUTPUT_DIR}/rfdiffusion/design_*.pdb 2>/dev/null | wc -l)
echo "Found ${design_count} RFdiffusion designs"

if [ ${design_count} -eq 0 ]; then
    echo "ERROR: No RFdiffusion designs found in ${CONFIG_OUTPUT_DIR}/rfdiffusion/"
    exit 1
fi

# ============================================================================
# Create SLURM Job Script
# ============================================================================
JOB_SCRIPT="${BATCH_LOG_DIR}/${CONFIG_NAME}_continue.slurm"

cat > "$JOB_SCRIPT" <<'EOF'
#!/bin/bash
#SBATCH --job-name=RFAB_Continue_NAME_PLACEHOLDER
#SBATCH --output=BATCH_LOG_DIR_PLACEHOLDER/NAME_PLACEHOLDER_continue_%j.out
#SBATCH --error=BATCH_LOG_DIR_PLACEHOLDER/NAME_PLACEHOLDER_continue_%j.err
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node20,node21,node24,node26,node27,node28
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00

# Setup environment
source /public/home/xuziyi/miniconda/etc/profile.d/conda.sh
conda activate RFantibody_HPC

source /public/software/profile.d/compiler_intel-2021.3.0.sh
source /public/software/profile.d/cuda12.2-env.sh
source /public/software/profile.d/fftw3_3.3.10_float-sse2-avx2.sh

# GLIBC setup
GLIBC_ROOT=/public/software/mathlib/glibc/2.34
GLIBC_LOADER=${GLIBC_ROOT}/lib/ld-linux-x86-64.so.2
GLIBC_LIB_PATH=${GLIBC_ROOT}/lib:${GLIBC_ROOT}/lib64

CONDA_LIB_DIR="${CONDA_PREFIX}/lib"
GLIBC_RUNTIME_PATH="${GLIBC_LIB_PATH}:${CONDA_LIB_DIR}"
if [ -n "${LD_LIBRARY_PATH:-}" ]; then
    GLIBC_RUNTIME_PATH="${GLIBC_RUNTIME_PATH}:${LD_LIBRARY_PATH}"
fi

export GLIBC_ROOT
export GLIBC_LOADER
export RFANTIBODY_GLIBC_LIBRARY_PATH="${GLIBC_RUNTIME_PATH}"
export RFANTIBODY_PYTHON_BIN="${CONDA_PREFIX}/bin/python"

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

export -f run_with_glibc

export PYTHONPATH=PROJECT_DIR_PLACEHOLDER/src:PROJECT_DIR_PLACEHOLDER/src/rfantibody/rfdiffusion:PROJECT_DIR_PLACEHOLDER/src/rfantibody:PROJECT_DIR_PLACEHOLDER/include/SE3Transformer:$PYTHONPATH
export PATH=PROJECT_DIR_PLACEHOLDER/include/USalign:$PATH
export WEIGHTS_DIR=WEIGHTS_DIR_PLACEHOLDER
export ICECREAM_COLORS=never
export HYDRA_FULL_ERROR=1
export PYTHONUNBUFFERED=1
export NO_COLOR=1
export TERM=dumb

cd PROJECT_DIR_PLACEHOLDER || exit 1

echo "=========================================="
echo "Continue Pipeline: Filtering → ProteinMPNN → RF2 → Analysis"
echo "=========================================="
echo "Configuration: NAME_PLACEHOLDER"
echo "Started: $(date)"
echo "Output: CONFIG_OUTPUT_DIR_PLACEHOLDER"
echo "=========================================="
echo ""

# ============================================================================
# Step 1: Filtering
# ============================================================================
echo "==> Step 1: Filtering by RHC-HIS contacts (<6A)"
echo ""

run_with_glibc "${RFANTIBODY_PYTHON_BIN}" "PROJECT_DIR_PLACEHOLDER/count_contacts_with_hotspot.py" \
    --pdb-dir "CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion" \
    --workdir "CONFIG_OUTPUT_DIR_PLACEHOLDER" \
    --mode backbone \
    --primary-cutoff 6 \
    --cutoffs 5 6 7 8

CONTACT_EXIT=$?
if [ $CONTACT_EXIT -ne 0 ]; then
    echo "✗ Contact counting failed with exit code $CONTACT_EXIT"
    exit 1
fi

# Extract passing designs
CSV_FILE="CONFIG_OUTPUT_DIR_PLACEHOLDER/contacts_with_hotspot_summary.csv"
MAPFILE=$(mktemp)

run_with_glibc "${RFANTIBODY_PYTHON_BIN}" -c "
import csv
csv_path = '${CSV_FILE}'
primary_col = 'bb_<=6A'

with open(csv_path) as handle:
    reader = csv.DictReader(handle)
    if primary_col not in reader.fieldnames:
        raise ValueError(f'Column {primary_col} not found')
    for row in reader:
        count = int(row[primary_col])
        if count > 0:
            print(row['PDB File'])
" > "${MAPFILE}"

MATCHING_PDBS=()
while IFS= read -r line; do
    MATCHING_PDBS+=("$line")
done < "${MAPFILE}"
rm -f "${MAPFILE}"

echo "Found ${#MATCHING_PDBS[@]} designs with hotspot contacts"

if [ ${#MATCHING_PDBS[@]} -eq 0 ]; then
    echo "✗ No designs passed filtering"
    exit 1
fi

# Move original to rfdiffusion_raw
ORIG_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion"
RAW_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion_raw"
NEW_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion"

if [ ! -d "${RAW_DIR}" ]; then
    mv "${ORIG_DIR}" "${RAW_DIR}"
    mkdir -p "${NEW_DIR}"
    
    for pdb in "${MATCHING_PDBS[@]}"; do
        cp "${RAW_DIR}/${pdb}" "${NEW_DIR}/"
    done
fi

filtered_count=$(ls ${NEW_DIR}/design_*.pdb 2>/dev/null | wc -l)
echo "✓ Filtering complete: ${filtered_count} designs retained"
echo ""

# ============================================================================
# Step 2: ProteinMPNN
# ============================================================================
echo "==> Step 2: ProteinMPNN (32 sequences per design)"
echo ""

bash PROJECT_DIR_PLACEHOLDER/run_antibody_design.sh --proteinmpnn-only --output-dir CONFIG_OUTPUT_DIR_PLACEHOLDER </dev/null

MPNN_EXIT=$?
if [ $MPNN_EXIT -ne 0 ]; then
    echo "✗ ProteinMPNN failed with exit code $MPNN_EXIT"
    exit 1
fi

mpnn_count=$(ls CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn/*.pdb 2>/dev/null | wc -l)
echo "✓ ProteinMPNN completed: ${mpnn_count} PDB files"
echo ""

# ============================================================================
# Step 2.5: Filter2 - Sidechain and Backbone Contact Analysis
# ============================================================================
echo "==> Step 2.5: Filter2 - Contact analysis (sidechain + backbone)"
echo ""

echo "  Running sidechain analysis..."
run_with_glibc "${RFANTIBODY_PYTHON_BIN}" "PROJECT_DIR_PLACEHOLDER/count_contacts_with_hotspot.py" \
    --pdb-dir "CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn" \
    --workdir "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_sc" \
    --mode sidechain \
    --cutoffs 5 6 > /dev/null

FILTER2_SC_EXIT=$?
if [ $FILTER2_SC_EXIT -ne 0 ]; then
    echo "✗ Filter2 sidechain analysis failed"
    exit 1
fi

echo "  Running backbone analysis..."
run_with_glibc "${RFANTIBODY_PYTHON_BIN}" "PROJECT_DIR_PLACEHOLDER/count_contacts_with_hotspot.py" \
    --pdb-dir "CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn" \
    --workdir "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_bb" \
    --mode backbone \
    --cutoffs 5 6 > /dev/null

FILTER2_BB_EXIT=$?
if [ $FILTER2_BB_EXIT -ne 0 ]; then
    echo "✗ Filter2 backbone analysis failed"
    exit 1
fi

echo "  Merging results..."
python3 << PYEOF
import csv

sc_file = "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_sc/contacts_with_hotspot_summary.csv"
bb_file = "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_bb/contacts_with_hotspot_summary.csv"
output_file = "CONFIG_OUTPUT_DIR_PLACEHOLDER/filter2_sidechain_contacts.csv"

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
rm -rf "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_sc" "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_bb"

filter2_count=$(wc -l < "CONFIG_OUTPUT_DIR_PLACEHOLDER/filter2_sidechain_contacts.csv")
filter2_count=$((filter2_count - 1))  # Subtract header
echo "✓ Filter2 complete: ${filter2_count} designs analyzed"
echo "  Results: CONFIG_OUTPUT_DIR_PLACEHOLDER/filter2_sidechain_contacts.csv"
echo ""

# ============================================================================
# Step 3: RF2
# ============================================================================
echo "==> Step 3: RF2 structure prediction"
echo ""

echo "Processing ${mpnn_count} sequences from ProteinMPNN"
echo ""

bash PROJECT_DIR_PLACEHOLDER/run_antibody_design.sh --rf2-only --output-dir CONFIG_OUTPUT_DIR_PLACEHOLDER </dev/null

RF2_EXIT=$?
if [ $RF2_EXIT -ne 0 ]; then
    echo "✗ RF2 failed with exit code $RF2_EXIT"
    exit 1
fi

rf2_count=$(ls CONFIG_OUTPUT_DIR_PLACEHOLDER/rf2/*.pdb 2>/dev/null | wc -l)
echo "✓ RF2 completed: ${rf2_count} structures"
echo ""

# ============================================================================
# Step 4: Analysis
# ============================================================================
echo "==> Step 4: Analysis"
echo ""

run_with_glibc "${RFANTIBODY_PYTHON_BIN}" "PROJECT_DIR_PLACEHOLDER/analyze_rf2_pdb.py" \
    --rf2-dir "CONFIG_OUTPUT_DIR_PLACEHOLDER/rf2" \
    --output "CONFIG_OUTPUT_DIR_PLACEHOLDER/rf2_pdb_analysis.csv"

if [ $? -eq 0 ]; then
    echo "✓ Analysis complete"
else
    echo "⚠ Analysis completed with warnings"
fi

echo ""
echo "=========================================="
echo "Pipeline Completed: $(date)"
echo "=========================================="
echo "Summary:"
echo "  Initial designs: DESIGN_COUNT_PLACEHOLDER"
echo "  After filtering: ${filtered_count}"
echo "  ProteinMPNN outputs: ${mpnn_count}"
echo "  RF2 structures: ${rf2_count}"
echo "  Results: CONFIG_OUTPUT_DIR_PLACEHOLDER/"
echo "=========================================="
EOF

# Replace placeholders
sed -i "s@NAME_PLACEHOLDER@${CONFIG_NAME}@g" "$JOB_SCRIPT"
sed -i "s@CONFIG_OUTPUT_DIR_PLACEHOLDER@${CONFIG_OUTPUT_DIR}@g" "$JOB_SCRIPT"
sed -i "s@BATCH_LOG_DIR_PLACEHOLDER@${BATCH_LOG_DIR}@g" "$JOB_SCRIPT"
sed -i "s@PROJECT_DIR_PLACEHOLDER@${PROJECT_DIR}@g" "$JOB_SCRIPT"
sed -i "s@WEIGHTS_DIR_PLACEHOLDER@${WEIGHTS_DIR}@g" "$JOB_SCRIPT"
sed -i "s@DESIGN_COUNT_PLACEHOLDER@${design_count}@g" "$JOB_SCRIPT"

# ============================================================================
# Submit Job
# ============================================================================
echo "Submitting SLURM job..."
JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")

echo ""
echo "=========================================="
echo "Job Submitted Successfully"
echo "=========================================="
echo "Job ID: ${JOB_ID}"
echo "Configuration: ${CONFIG_NAME}"
echo "Output Directory: ${CONFIG_OUTPUT_DIR}"
echo "Log Files:"
echo "  Output: ${BATCH_LOG_DIR}/${CONFIG_NAME}_continue_${JOB_ID}.out"
echo "  Error:  ${BATCH_LOG_DIR}/${CONFIG_NAME}_continue_${JOB_ID}.err"
echo ""
echo "To monitor job:"
echo "  squeue -j ${JOB_ID}"
echo "  tail -f ${BATCH_LOG_DIR}/${CONFIG_NAME}_continue_${JOB_ID}.out"
echo ""
echo "To cancel job:"
echo "  scancel ${JOB_ID}"
echo ""
echo "=========================================="
