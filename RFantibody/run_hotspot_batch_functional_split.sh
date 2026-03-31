#!/bin/bash

# ============================================================================
# Stage 1
#   - RFdiffusion
#   - Filter1
#   - ProteinMPNN
#
# Stage 2
#   - RF2
#   - Analysis
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

# Validate required variables
if [ -z "$PROJECT_DIR" ] || [ ! -d "$PROJECT_DIR" ]; then
    echo "ERROR: PROJECT_DIR not set or does not exist: $PROJECT_DIR"
    exit 1
fi

# Create batch ID for this submission
BATCH_ID=$(date +%Y%m%d_%H%M%S)
BATCH_LOG_DIR="${PROJECT_DIR}/logs/functional_split_${BATCH_ID}"
mkdir -p "$BATCH_LOG_DIR"

echo "=========================================="
echo "Two-Stage RFantibody Pipeline"
echo "Functional Split Mode"
echo "=========================================="
echo "Stage 1: RFdiffusion → Filter1 → ProteinMPNN"
echo "Stage 2: RF2 → Analysis"
echo ""
echo "Batch ID: ${BATCH_ID}"
echo "Batch log directory: ${BATCH_LOG_DIR}"
echo ""

# ============================================================================
# Define Configurations
# ============================================================================
# Format: "NAME|CDR_LOOPS|HOTSPOTS|NUM_DESIGNS|T_VALUE"
declare -a REPLICA_CONFIGS=(
    "200T_600_32_L1_10-13_H3_11-14_T82-159_replica4|[L1:10-13,L2:7,L3:9-11,H1:7,H2:6,H3:11-14]|[T82,T85,T155,T156,T159]|600|200"
    "200T_600_32_L1_10-13_H3_11-14_T82-159_replica5|[L1:10-13,L2:7,L3:9-11,H1:7,H2:6,H3:11-14]|[T82,T85,T155,T156,T159]|600|200"
)

# ============================================================================
# Job tracking
# ============================================================================
declare -a STAGE1_JOB_IDS=()
declare -a STAGE2_JOB_IDS=()
declare -a JOB_NAMES=()

# ============================================================================
# Submit jobs for each configuration
# ============================================================================
for config in "${REPLICA_CONFIGS[@]}"; do
    IFS='|' read -r NAME CDR_LOOPS HOTSPOTS NUM_DESIGNS T_VALUE <<< "$config"
    
    echo "----------------------------------------"
    echo "Configuration: $NAME"
    echo "  Hotspots:  $HOTSPOTS"
    echo "  CDR Loops: $CDR_LOOPS"
    echo "  T:         $T_VALUE"
    echo "  Designs:   $NUM_DESIGNS"
    
    CONFIG_OUTPUT_DIR="${PROJECT_DIR}/outputs/${NAME}"
    echo "  Output: $CONFIG_OUTPUT_DIR"
    
    mkdir -p "${CONFIG_OUTPUT_DIR}"
    
    # Create temporary config for this task
    TEMP_CONFIG="${BATCH_LOG_DIR}/${NAME}_run_config.txt"
    cat > "$TEMP_CONFIG" <<CFGEOF
INPUT_PDB=${INPUT_PDB}
FRAMEWORK_PDB=${FRAMEWORK_PDB}
OUTPUT_DIR=${CONFIG_OUTPUT_DIR}
WEIGHTS_DIR=${WEIGHTS_DIR}
PROJECT_DIR=${PROJECT_DIR}
CFGEOF
    
    # ========================================================================
    # Stage 1 Job Script: RFdiffusion → Filtering → ProteinMPNN
    # ========================================================================
    STAGE1_SCRIPT="${BATCH_LOG_DIR}/${NAME}_stage1.slurm"
    
    cat > "$STAGE1_SCRIPT" <<'EOF'
#!/bin/bash
#SBATCH --job-name=RFAB_S1_NAME_PLACEHOLDER
#SBATCH --output=BATCH_LOG_DIR_PLACEHOLDER/NAME_PLACEHOLDER_stage1_%j.out
#SBATCH --error=BATCH_LOG_DIR_PLACEHOLDER/NAME_PLACEHOLDER_stage1_%j.err
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

# Generate random seed
JOB_SEED=$RANDOM

# Environment setup
export PYTHONPATH=PROJECT_DIR_PLACEHOLDER/src:PROJECT_DIR_PLACEHOLDER/src/rfantibody/rfdiffusion:PROJECT_DIR_PLACEHOLDER/src/rfantibody:PROJECT_DIR_PLACEHOLDER/include/SE3Transformer:$PYTHONPATH
export PATH=PROJECT_DIR_PLACEHOLDER/include/USalign:$PATH
export WEIGHTS_DIR=WEIGHTS_DIR_PLACEHOLDER
export ICECREAM_COLORS=never
export HYDRA_FULL_ERROR=1
export PYTHONHASHSEED=${JOB_SEED}
export PYTHONUNBUFFERED=1
export NO_COLOR=1
export TERM=dumb
export CONFIG_FILE=TEMP_CONFIG_PLACEHOLDER

cd PROJECT_DIR_PLACEHOLDER || exit 1

echo "=========================================="
echo "STAGE 1: RFdiffusion → Filtering → ProteinMPNN"
echo "=========================================="
echo "Configuration: NAME_PLACEHOLDER"
echo "Started: $(date)"
echo "Random Seed: ${JOB_SEED}"
echo ""
echo "Parameters:"
echo "  Hotspots:  HOTSPOTS_PLACEHOLDER"
echo "  CDR Loops: CDR_LOOPS_PLACEHOLDER"
echo "  T:         T_VALUE_PLACEHOLDER"
echo "  Designs:   NUM_DESIGNS_PLACEHOLDER"
echo "  Output:    CONFIG_OUTPUT_DIR_PLACEHOLDER"
echo "=========================================="
echo ""

# ============================================================================
# Step 1: RFdiffusion
# ============================================================================
echo "==> Step 1: RFdiffusion (NUM_DESIGNS_PLACEHOLDER designs)"
echo ""

# Create symlink for config (handle race condition for parallel jobs)
mkdir -p scripts/config

# Check if symlink already exists and is correct
EXPECTED_TARGET="PROJECT_DIR_PLACEHOLDER/src/rfantibody/rfdiffusion/config/inference"
if [ -L "scripts/config/inference" ]; then
    CURRENT_TARGET=$(readlink "scripts/config/inference")
    if [ "${CURRENT_TARGET}" = "${EXPECTED_TARGET}" ]; then
        echo "Config symlink already exists and is correct"
    else
        echo "Updating config symlink to correct target"
        rm -f scripts/config/inference
        ln -sfn "${EXPECTED_TARGET}" scripts/config/inference
    fi
else
    echo "Creating config symlink"
    rm -rf scripts/config/inference
    ln -sfn "${EXPECTED_TARGET}" scripts/config/inference
fi

# Verify symlink (use -e to check if target exists)
if [ ! -e "scripts/config/inference" ]; then
    echo "ERROR: Config symlink verification failed" >&2
    echo "  Expected: scripts/config/inference -> ${EXPECTED_TARGET}" >&2
    echo "  Check if target directory exists: ${EXPECTED_TARGET}" >&2
    exit 1
fi

echo "Config symlink verified: scripts/config/inference"

# Create output directory structure
mkdir -p "CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion"

# Run RFdiffusion with random seed (non-deterministic mode)
run_with_glibc "${RFANTIBODY_PYTHON_BIN}" scripts/rfdiffusion_inference.py \
    --config-name antibody \
    antibody.target_pdb="INPUT_PDB_PLACEHOLDER" \
    antibody.framework_pdb="FRAMEWORK_PDB_PLACEHOLDER" \
    inference.ckpt_override_path="WEIGHTS_DIR_PLACEHOLDER/RFdiffusion_Ab.pt" \
    'ppi.hotspot_res=HOTSPOTS_PLACEHOLDER' \
    'antibody.design_loops=CDR_LOOPS_PLACEHOLDER' \
    inference.num_designs=NUM_DESIGNS_PLACEHOLDER \
    diffuser.T=T_VALUE_PLACEHOLDER \
    ++inference.seed=${JOB_SEED} \
    inference.output_prefix="CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion/design"

RFDIFF_EXIT=$?

if [ $RFDIFF_EXIT -ne 0 ]; then
    echo "✗ RFdiffusion failed with exit code $RFDIFF_EXIT"
    exit $RFDIFF_EXIT
fi

design_count=$(ls CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion/design_*.pdb 2>/dev/null | wc -l)
echo "✓ RFdiffusion completed: ${design_count} designs"
echo ""

# ============================================================================
# Step 2: Filtering
# ============================================================================
echo "==> Step 2: Filtering by RHC-HIS contacts (<6A)"
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
# Step 3: ProteinMPNN
# ============================================================================
echo "==> Step 3: ProteinMPNN (32 sequences per design)"
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
# Step 3.5: Filter2 - DISABLED (COMMENTED OUT)
# ============================================================================
# Filter2 has been disabled per user request (2026-02-06)
# Final CSV outputs: contacts_with_hotspot_summary.csv + rf2_pdb_analysis.csv only
#
# echo "==> Step 3.5: Filter2 - Contact analysis (sidechain + backbone)"
# echo ""
#
# echo "  Running sidechain analysis..."
# run_with_glibc "${RFANTIBODY_PYTHON_BIN}" "PROJECT_DIR_PLACEHOLDER/count_contacts_with_hotspot.py" \
#     --pdb-dir "CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn" \
#     --workdir "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_sc" \
#     --mode sidechain \
#     --cutoffs 5 6 > /dev/null
#
# FILTER2_SC_EXIT=$?
# if [ $FILTER2_SC_EXIT -ne 0 ]; then
#     echo "✗ Filter2 sidechain analysis failed"
#     exit 1
# fi
#
# echo "  Running backbone analysis..."
# run_with_glibc "${RFANTIBODY_PYTHON_BIN}" "PROJECT_DIR_PLACEHOLDER/count_contacts_with_hotspot.py" \
#     --pdb-dir "CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn" \
#     --workdir "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_bb" \
#     --mode backbone \
#     --cutoffs 5 6 > /dev/null
#
# FILTER2_BB_EXIT=$?
# if [ $FILTER2_BB_EXIT -ne 0 ]; then
#     echo "✗ Filter2 backbone analysis failed"
#     exit 1
# fi
#
# echo "  Merging results..."
# run_with_glibc "${RFANTIBODY_PYTHON_BIN}" << PYEOF
# import csv
#
# sc_file = "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_sc/contacts_with_hotspot_summary.csv"
# bb_file = "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_bb/contacts_with_hotspot_summary.csv"
# output_file = "CONFIG_OUTPUT_DIR_PLACEHOLDER/filter2_sidechain_contacts.csv"
#
# # Read both CSVs
# with open(sc_file) as f:
#     sc_data = {row['PDB File']: row for row in csv.DictReader(f)}
#
# with open(bb_file) as f:
#     bb_data = {row['PDB File']: row for row in csv.DictReader(f)}
#
# # Merge data
# merged = []
# for pdb_name in sorted(sc_data.keys()):
#     merged_row = {
#         'PDB File': pdb_name,
#         'sc_<=5A': sc_data[pdb_name]['sc_<=5A'],
#         'bb_<=5A': bb_data[pdb_name]['bb_<=5A'],
#         'sc_<=6A': sc_data[pdb_name]['sc_<=6A'],
#         'bb_<=6A': bb_data[pdb_name]['bb_<=6A'],
#     }
#     merged.append(merged_row)
#
# # Write merged CSV
# with open(output_file, 'w', newline='') as f:
#     fieldnames = ['PDB File', 'sc_<=5A', 'bb_<=5A', 'sc_<=6A', 'bb_<=6A']
#     writer = csv.DictWriter(f, fieldnames=fieldnames)
#     writer.writeheader()
#     writer.writerows(merged)
#
# print(f"Merged {len(merged)} designs")
# PYEOF
#
# # Cleanup temp directories
# rm -rf "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_sc" "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_bb"
#
# filter2_count=$(wc -l < "CONFIG_OUTPUT_DIR_PLACEHOLDER/filter2_sidechain_contacts.csv")
# filter2_count=$((filter2_count - 1))  # Subtract header
# echo "✓ Filter2 complete: ${filter2_count} designs analyzed"
# echo "  Results: CONFIG_OUTPUT_DIR_PLACEHOLDER/filter2_sidechain_contacts.csv"
# echo "  Columns: sc_<=5A, bb_<=5A, sc_<=6A, bb_<=6A"
# echo ""

echo "=========================================="
echo "STAGE 1 Completed: $(date)"
echo "=========================================="
echo "Summary:"
echo "  Initial designs: ${design_count}"
echo "  After Filter1 (backbone): ${filtered_count}"
echo "  ProteinMPNN outputs: ${mpnn_count}"
echo "  Filter2: DISABLED (skipped)"
echo ""
echo "Stage 2 will automatically start to process ${mpnn_count} sequences"
echo "========================================="
EOF

    # Replace placeholders
    sed -i "s@NAME_PLACEHOLDER@${NAME}@g" "$STAGE1_SCRIPT"
    sed -i "s@HOTSPOTS_PLACEHOLDER@${HOTSPOTS}@g" "$STAGE1_SCRIPT"
    sed -i "s@CDR_LOOPS_PLACEHOLDER@${CDR_LOOPS}@g" "$STAGE1_SCRIPT"
    sed -i "s@NUM_DESIGNS_PLACEHOLDER@${NUM_DESIGNS}@g" "$STAGE1_SCRIPT"
    sed -i "s@T_VALUE_PLACEHOLDER@${T_VALUE}@g" "$STAGE1_SCRIPT"
    sed -i "s@CONFIG_OUTPUT_DIR_PLACEHOLDER@${CONFIG_OUTPUT_DIR}@g" "$STAGE1_SCRIPT"
    sed -i "s@INPUT_PDB_PLACEHOLDER@${INPUT_PDB}@g" "$STAGE1_SCRIPT"
    sed -i "s@FRAMEWORK_PDB_PLACEHOLDER@${FRAMEWORK_PDB}@g" "$STAGE1_SCRIPT"
    sed -i "s@BATCH_LOG_DIR_PLACEHOLDER@${BATCH_LOG_DIR}@g" "$STAGE1_SCRIPT"
    sed -i "s@PROJECT_DIR_PLACEHOLDER@${PROJECT_DIR}@g" "$STAGE1_SCRIPT"
    sed -i "s@WEIGHTS_DIR_PLACEHOLDER@${WEIGHTS_DIR}@g" "$STAGE1_SCRIPT"
    sed -i "s@TEMP_CONFIG_PLACEHOLDER@${TEMP_CONFIG}@g" "$STAGE1_SCRIPT"
    
    # ========================================================================
    # Stage 2 Job Script: RF2 → Analysis
    # ========================================================================
    STAGE2_SCRIPT="${BATCH_LOG_DIR}/${NAME}_stage2.slurm"
    
    cat > "$STAGE2_SCRIPT" <<'EOF'
#!/bin/bash
#SBATCH --job-name=RFAB_S2_NAME_PLACEHOLDER
#SBATCH --output=BATCH_LOG_DIR_PLACEHOLDER/NAME_PLACEHOLDER_stage2_%j.out
#SBATCH --error=BATCH_LOG_DIR_PLACEHOLDER/NAME_PLACEHOLDER_stage2_%j.err
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
export CONFIG_FILE=TEMP_CONFIG_PLACEHOLDER

cd PROJECT_DIR_PLACEHOLDER || exit 1

echo "=========================================="
echo "STAGE 2: RF2 → Analysis"
echo "=========================================="
echo "Configuration: NAME_PLACEHOLDER"
echo "Started: $(date)"
echo ""
echo "Output: CONFIG_OUTPUT_DIR_PLACEHOLDER"
echo "=========================================="
echo ""

# ============================================================================
# Step 3.5: Filter2 - DISABLED (COMMENTED OUT)
# ============================================================================
# Filter2 fallback has been disabled per user request (2026-02-06)
# Final CSV outputs: contacts_with_hotspot_summary.csv + rf2_pdb_analysis.csv only
#
# if [ ! -f "CONFIG_OUTPUT_DIR_PLACEHOLDER/filter2_sidechain_contacts.csv" ]; then
#     echo "==> Step 3.5: Filter2 - Contact analysis (sidechain + backbone)"
#     echo ""
#
#     echo "  Running sidechain analysis..."
#     run_with_glibc "${RFANTIBODY_PYTHON_BIN}" "PROJECT_DIR_PLACEHOLDER/count_contacts_with_hotspot.py" \
#         --pdb-dir "CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn" \
#         --workdir "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_sc" \
#         --mode sidechain \
#         --cutoffs 5 6 > /dev/null
#
#     FILTER2_SC_EXIT=$?
#     if [ $FILTER2_SC_EXIT -ne 0 ]; then
#         echo "✗ Filter2 sidechain analysis failed"
#         exit 1
#     fi
#
#     echo "  Running backbone analysis..."
#     run_with_glibc "${RFANTIBODY_PYTHON_BIN}" "PROJECT_DIR_PLACEHOLDER/count_contacts_with_hotspot.py" \
#         --pdb-dir "CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn" \
#         --workdir "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_bb" \
#         --mode backbone \
#         --cutoffs 5 6 > /dev/null
#
#     FILTER2_BB_EXIT=$?
#     if [ $FILTER2_BB_EXIT -ne 0 ]; then
#         echo "✗ Filter2 backbone analysis failed"
#         exit 1
#     fi
#
#     echo "  Merging results..."
#     run_with_glibc "${RFANTIBODY_PYTHON_BIN}" << PYEOF
# import csv
#
# sc_file = "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_sc/contacts_with_hotspot_summary.csv"
# bb_file = "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_bb/contacts_with_hotspot_summary.csv"
# output_file = "CONFIG_OUTPUT_DIR_PLACEHOLDER/filter2_sidechain_contacts.csv"
#
# # Read both CSVs
# with open(sc_file) as f:
#     sc_data = {row['PDB File']: row for row in csv.DictReader(f)}
#
# with open(bb_file) as f:
#     bb_data = {row['PDB File']: row for row in csv.DictReader(f)}
#
# # Merge data
# merged = []
# for pdb_name in sorted(sc_data.keys()):
#     merged_row = {
#         'PDB File': pdb_name,
#         'sc_<=5A': sc_data[pdb_name]['sc_<=5A'],
#         'bb_<=5A': bb_data[pdb_name]['bb_<=5A'],
#         'sc_<=6A': sc_data[pdb_name]['sc_<=6A'],
#         'bb_<=6A': bb_data[pdb_name]['bb_<=6A'],
#     }
#     merged.append(merged_row)
#
# # Write merged CSV
# with open(output_file, 'w', newline='') as f:
#     fieldnames = ['PDB File', 'sc_<=5A', 'bb_<=5A', 'sc_<=6A', 'bb_<=6A']
#     writer = csv.DictWriter(f, fieldnames=fieldnames)
#     writer.writeheader()
#     writer.writerows(merged)
#
# print(f"Merged {len(merged)} designs")
# PYEOF
#
#     # Cleanup temp directories
#     rm -rf "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_sc" "CONFIG_OUTPUT_DIR_PLACEHOLDER/.tmp_bb"
#
#     filter2_count=$(wc -l < "CONFIG_OUTPUT_DIR_PLACEHOLDER/filter2_sidechain_contacts.csv")
#     filter2_count=$((filter2_count - 1))  # Subtract header
#     echo "✓ Filter2 complete: ${filter2_count} designs analyzed"
#     echo "  Results: CONFIG_OUTPUT_DIR_PLACEHOLDER/filter2_sidechain_contacts.csv"
#     echo "  Columns: sc_<=5A, bb_<=5A, sc_<=6A, bb_<=6A"
#     echo ""
# else
#     echo "==> Step 3.5: Filter2 already completed (skipping)"
#     echo ""
# fi

# ============================================================================
# Step 4: RF2
# ============================================================================

# ============================================================================
# Step 4: RF2
# ============================================================================
echo "==> Step 4: RF2 structure prediction"
echo ""

mpnn_count=$(ls CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn/*.pdb 2>/dev/null | wc -l)
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
# Step 5: Analysis
# ============================================================================
echo "==> Step 5: Analysis"
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
echo "STAGE 2 Completed: $(date)"
echo "=========================================="
echo "Summary:"
echo "  ProteinMPNN inputs: ${mpnn_count}"
echo "  RF2 structures: ${rf2_count}"
echo "  Filter2: DISABLED (skipped)"
echo "  Results: CONFIG_OUTPUT_DIR_PLACEHOLDER/"
echo "=========================================="
EOF

    # Replace placeholders
    sed -i "s@NAME_PLACEHOLDER@${NAME}@g" "$STAGE2_SCRIPT"
    sed -i "s@CONFIG_OUTPUT_DIR_PLACEHOLDER@${CONFIG_OUTPUT_DIR}@g" "$STAGE2_SCRIPT"
    sed -i "s@BATCH_LOG_DIR_PLACEHOLDER@${BATCH_LOG_DIR}@g" "$STAGE2_SCRIPT"
    sed -i "s@PROJECT_DIR_PLACEHOLDER@${PROJECT_DIR}@g" "$STAGE2_SCRIPT"
    sed -i "s@WEIGHTS_DIR_PLACEHOLDER@${WEIGHTS_DIR}@g" "$STAGE2_SCRIPT"
    sed -i "s@TEMP_CONFIG_PLACEHOLDER@${TEMP_CONFIG}@g" "$STAGE2_SCRIPT"
    
    # ========================================================================
    # Submit jobs with dependency
    # ========================================================================
    echo "  Submitting Stage 1 job..."
    STAGE1_JOB_ID=$(sbatch --parsable "$STAGE1_SCRIPT")
    STAGE1_JOB_IDS+=("$STAGE1_JOB_ID")
    echo "  ✓ Stage 1 job: $STAGE1_JOB_ID"
    
    echo "  Submitting Stage 2 job (depends on Stage 1)..."
    STAGE2_JOB_ID=$(sbatch --parsable --dependency=afterany:$STAGE1_JOB_ID "$STAGE2_SCRIPT")
    STAGE2_JOB_IDS+=("$STAGE2_JOB_ID")
    echo "  ✓ Stage 2 job: $STAGE2_JOB_ID (waits for $STAGE1_JOB_ID)"
    
    JOB_NAMES+=("$NAME")
    echo ""
done

# ============================================================================
# Summary
# ============================================================================
echo "=========================================="
echo "Two-Stage Batch Submission Complete"
echo "=========================================="
echo ""
echo "Submitted job pairs:"
for i in "${!JOB_NAMES[@]}"; do
    echo "  ${JOB_NAMES[$i]}:"
    echo "    Stage 1: ${STAGE1_JOB_IDS[$i]} (RFdiff→Filtering→MPNN)"
    echo "    Stage 2: ${STAGE2_JOB_IDS[$i]} (RF2→Analysis, depends on Stage 1)"
done
echo ""
echo "Log directory: ${BATCH_LOG_DIR}"
echo ""
echo "To monitor jobs:"
echo "  squeue -u \$USER"
echo ""
echo "To cancel all jobs:"
echo "  scancel ${STAGE1_JOB_IDS[@]} ${STAGE2_JOB_IDS[@]}"
echo ""
echo "=========================================="
