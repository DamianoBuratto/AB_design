#!/bin/bash

# ============================================================================
# Resume Replica Jobs - Smart Recovery
# ============================================================================
# Purpose: Resume interrupted replica jobs from the appropriate stage
#
# Replica 1 & 2:
#   - RFdiffusion complete (700 → filtered to 404/422)
#   - Filtering complete
#   - ProteinMPNN complete (12,928 and 13,504 sequences)
#   - RF2 incomplete (9% and 20%)
#   → Resume from: RF2 step
#
# Replica 3:
#   - RFdiffusion incomplete (436/700 = 62%)
#   - Filtering not run
#   - ProteinMPNN not run
#   - RF2 not run
#   → Resume from: Filtering step (using 436 existing designs)
#
# Usage:
#   bash resume_replicas_from_rf2.sh
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

source "$CONFIG_FILE"

if [ -z "$PROJECT_DIR" ] || [ ! -d "$PROJECT_DIR" ]; then
    echo "ERROR: PROJECT_DIR not set or does not exist: $PROJECT_DIR"
    exit 1
fi

# Create batch ID
BATCH_ID=$(date +%Y%m%d_%H%M%S)
RESUME_LOG_DIR="${PROJECT_DIR}/logs/resume_rf2_${BATCH_ID}"
mkdir -p "$RESUME_LOG_DIR"

echo "=========================================="
echo "Resume Replica Jobs - Smart Recovery"
echo "=========================================="
echo "Batch ID: ${BATCH_ID}"
echo "Resume log directory: $RESUME_LOG_DIR"
echo ""

# ============================================================================
# Define replica configurations
# ============================================================================
# Format: NAME|START_STEP|DESCRIPTION
# START_STEP: 2=Filtering, 4=RF2, 5=Analysis
# Note: Replica2 completed, Replica3 excluded (handled by run_hotspot_batch_functional_split.sh)
declare -a REPLICA_CONFIGS=(
    "200T_700_32_L1_10-13_H3_11-14_T82-159_replica1|4|RF2 incomplete - resume only"
)

# ============================================================================
# Job tracking
# ============================================================================
declare -a JOB_IDS=()
declare -a JOB_NAMES=()

# ============================================================================
# Submit resume jobs
# ============================================================================
for config in "${REPLICA_CONFIGS[@]}"; do
    IFS='|' read -r NAME START_STEP DESCRIPTION <<< "$config"
    
    echo "----------------------------------------"
    echo "Configuration: $NAME"
    echo "  Status: $DESCRIPTION"
    echo "  Resume from: Step $START_STEP"
    
    CONFIG_OUTPUT_DIR="${PROJECT_DIR}/outputs/${NAME}"
    echo "  Output: $CONFIG_OUTPUT_DIR"
    
    if [ ! -d "$CONFIG_OUTPUT_DIR" ]; then
        echo "  WARNING: Output directory not found, skipping"
        echo ""
        continue
    fi
    
    # Create temp config
    TEMP_CONFIG="${RESUME_LOG_DIR}/${NAME}_run_config.txt"
    cat > "$TEMP_CONFIG" <<CFGEOF
INPUT_PDB=${INPUT_PDB}
FRAMEWORK_PDB=${FRAMEWORK_PDB}
OUTPUT_DIR=${CONFIG_OUTPUT_DIR}
WEIGHTS_DIR=${WEIGHTS_DIR}
PROJECT_DIR=${PROJECT_DIR}
CFGEOF
    
    # Create job script
    JOB_SCRIPT="${RESUME_LOG_DIR}/${NAME}_resume.slurm"
    
    cat > "$JOB_SCRIPT" <<'EOF'
#!/bin/bash
#SBATCH --job-name=Resume_NAME_PLACEHOLDER
#SBATCH --output=BATCH_LOG_DIR_PLACEHOLDER/NAME_PLACEHOLDER_%j.out
#SBATCH --error=BATCH_LOG_DIR_PLACEHOLDER/NAME_PLACEHOLDER_%j.err
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node20,node21,node24,node26,node27,node28
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00

START_STEP=START_STEP_PLACEHOLDER

# Setup environment
source /public/home/xuziyi/miniconda/etc/profile.d/conda.sh
conda activate RFantibody_HPC

source /public/software/profile.d/compiler_intel-2021.3.0.sh
source /public/software/profile.d/cuda12.2-env.sh
source /public/software/profile.d/fftw3_3.3.10_float-sse2-avx2.sh

# Setup GLIBC environment
GLIBC_ROOT=/public/software/mathlib/glibc/2.34
GLIBC_LOADER=${GLIBC_ROOT}/lib/ld-linux-x86-64.so.2
GLIBC_LIB_PATH=${GLIBC_ROOT}/lib:${GLIBC_ROOT}/lib64

CONDA_LIB_DIR="${CONDA_PREFIX}/lib"
if [ ! -d "${CONDA_LIB_DIR}" ]; then
    echo "Expected conda lib directory not found: ${CONDA_LIB_DIR}" >&2
    exit 1
fi

GLIBC_RUNTIME_PATH="${GLIBC_LIB_PATH}:${CONDA_LIB_DIR}"
if [ -n "${LD_LIBRARY_PATH:-}" ]; then
    GLIBC_RUNTIME_PATH="${GLIBC_RUNTIME_PATH}:${LD_LIBRARY_PATH}"
fi

export GLIBC_ROOT
export GLIBC_LOADER
export RFANTIBODY_GLIBC_LIBRARY_PATH="${GLIBC_RUNTIME_PATH}"
export RFANTIBODY_PYTHON_BIN="${CONDA_PREFIX}/bin/python"

if [ ! -x "${RFANTIBODY_PYTHON_BIN}" ]; then
    echo "Python interpreter not found at ${RFANTIBODY_PYTHON_BIN}" >&2
    exit 1
fi

if [ ! -x "${GLIBC_LOADER}" ]; then
    echo "Custom glibc loader not found at ${GLIBC_LOADER}" >&2
    exit 1
fi

run_with_glibc() {
    if [ "$#" -eq 0 ]; then
        echo "run_with_glibc: missing command" >&2
        return 1
    fi
    if [ -z "${GLIBC_LOADER:-}" ]; then
        echo "run_with_glibc: GLIBC_LOADER is not set" >&2
        return 1
    fi
    if [ ! -x "${GLIBC_LOADER}" ]; then
        echo "run_with_glibc: loader not executable: ${GLIBC_LOADER}" >&2
        return 1
    fi
    local cmd="$1"
    shift
    if [[ "${cmd}" != /* ]]; then
        cmd=$(command -v "${cmd}")
    fi
    if [ -z "${cmd}" ]; then
        echo "run_with_glibc: command not found" >&2
        return 1
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

cd PROJECT_DIR_PLACEHOLDER || {
    echo "ERROR: Failed to change to project directory" >&2
    exit 1
}

echo "=========================================="
echo "Resume Configuration: NAME_PLACEHOLDER"
echo "Started: $(date)"
echo "=========================================="
echo "Output Directory: CONFIG_OUTPUT_DIR_PLACEHOLDER"
echo "Working directory: $(pwd)"
echo "Resume from step: ${START_STEP}"
echo "=========================================="
echo ""

# ============================================================================
# Step 2: Filtering (only for replica3)
# ============================================================================
if [ ${START_STEP} -le 2 ]; then
    echo "==> Step 2: Filtering by RHC-HIS contacts (<6A)"
    echo ""
    
    RFDIFF_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion"
    
    if [ ! -d "$RFDIFF_DIR" ]; then
        echo "ERROR: RFdiffusion directory not found: $RFDIFF_DIR"
        exit 1
    fi
    
    design_count=$(ls ${RFDIFF_DIR}/*.pdb 2>/dev/null | wc -l)
    echo "Found ${design_count} RFdiffusion designs"
    
    run_with_glibc "${RFANTIBODY_PYTHON_BIN}" "PROJECT_DIR_PLACEHOLDER/count_contacts_with_hotspot.py" \
        --pdb-dir "${RFDIFF_DIR}" \
        --workdir "CONFIG_OUTPUT_DIR_PLACEHOLDER" \
        --mode backbone \
        --primary-cutoff 6 \
        --cutoffs 5 6 7 8
    
    if [ $? -ne 0 ]; then
        echo "✗ Contact counting failed"
        exit 1
    fi
    
    # Extract passing designs
    CSV_FILE="CONFIG_OUTPUT_DIR_PLACEHOLDER/contacts_with_hotspot_summary.csv"
    MAPFILE=$(mktemp)
    
    run_with_glibc "${RFANTIBODY_PYTHON_BIN}" -c "
import csv
csv_path = '${CSV_FILE}'
primary_col = '(<= 6A)'

with open(csv_path) as handle:
    reader = csv.DictReader(handle)
    if primary_col not in reader.fieldnames:
        raise ValueError(f'Column {primary_col} not found')
    for row in reader:
        count = int(row[primary_col])
        if count > 0:
            print(row['PDB'])
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
    
    # Move to rfdiffusion_raw if not already done
    ORIG_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion"
    RAW_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion_raw"
    NEW_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion"
    
    if [ ! -d "${RAW_DIR}" ]; then
        echo "Moving original designs to rfdiffusion_raw..."
        mv "${ORIG_DIR}" "${RAW_DIR}"
        mkdir -p "${NEW_DIR}"
        
        for pdb in "${MATCHING_PDBS[@]}"; do
            cp "${RAW_DIR}/${pdb}" "${NEW_DIR}/"
        done
    else
        echo "Note: rfdiffusion_raw already exists, skipping move"
    fi
    
    filtered_count=$(ls ${NEW_DIR}/design_*.pdb 2>/dev/null | wc -l)
    echo "✓ Filtering complete: ${filtered_count} designs retained"
    echo ""
fi

# ============================================================================
# Step 3: ProteinMPNN
# ============================================================================
if [ ${START_STEP} -le 3 ]; then
    echo "==> Step 3: ProteinMPNN (32 sequences per design)"
    echo ""
    
    PROTEINMPNN_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn"
    
    if [ -d "$PROTEINMPNN_DIR" ]; then
        existing_count=$(ls ${PROTEINMPNN_DIR}/*.pdb 2>/dev/null | wc -l)
        echo "Found ${existing_count} existing ProteinMPNN outputs"
        
        if [ $existing_count -gt 0 ]; then
            echo "ProteinMPNN already run, skipping"
        else
            bash PROJECT_DIR_PLACEHOLDER/run_antibody_design.sh --proteinmpnn-only --output-dir CONFIG_OUTPUT_DIR_PLACEHOLDER </dev/null
        fi
    else
        bash PROJECT_DIR_PLACEHOLDER/run_antibody_design.sh --proteinmpnn-only --output-dir CONFIG_OUTPUT_DIR_PLACEHOLDER </dev/null
    fi
    
    if [ $? -ne 0 ]; then
        echo "✗ ProteinMPNN failed"
        exit 1
    fi
    
    mpnn_count=$(ls CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn/*.pdb 2>/dev/null | wc -l)
    echo "✓ ProteinMPNN completed: ${mpnn_count} PDB files"
    echo ""
fi

# ============================================================================
# Step 4: RF2 (check completion status)
# ============================================================================
if [ ${START_STEP} -le 4 ]; then
    echo "==> Step 4: Checking RF2 completion status"
    echo ""
    
    RF2_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/rf2"
    PROTEINMPNN_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn"
    
    if [ ! -d "$PROTEINMPNN_DIR" ]; then
        echo "ERROR: ProteinMPNN directory not found: $PROTEINMPNN_DIR"
        exit 1
    fi
    
    mpnn_count=$(ls ${PROTEINMPNN_DIR}/*.pdb 2>/dev/null | wc -l)
    echo "ProteinMPNN PDB files: $mpnn_count"
    
    if [ -d "$RF2_DIR" ]; then
        rf2_count=$(ls ${RF2_DIR}/*.pdb 2>/dev/null | wc -l)
        echo "Existing RF2 structures: $rf2_count"
    else
        rf2_count=0
        echo "RF2 directory not found, will create"
        mkdir -p "$RF2_DIR"
    fi
    
    if [ $rf2_count -lt $mpnn_count ]; then
        echo ""
        echo "RF2 incomplete: $rf2_count / $mpnn_count structures"
        echo "Continuing RF2 predictions..."
        echo ""
        
        bash PROJECT_DIR_PLACEHOLDER/run_antibody_design.sh --rf2-only --output-dir CONFIG_OUTPUT_DIR_PLACEHOLDER </dev/null
        
        if [ $? -ne 0 ]; then
            echo "⚠ RF2 completed with errors or timeout"
            echo "You can re-run this script to continue from where it stopped"
        else
            rf2_count=$(ls ${RF2_DIR}/*.pdb 2>/dev/null | wc -l)
            echo "✓ RF2 completed: $rf2_count structures"
        fi
    else
        echo "✓ RF2 already complete: $rf2_count / $mpnn_count structures"
    fi
    echo ""
fi

# ============================================================================
# Step 5: Analysis
# ============================================================================
if [ ${START_STEP} -le 5 ]; then
    echo "==> Step 5: Running analysis to generate rf2_pdb_analysis.csv"
    echo ""
    
    ANALYSIS_CSV="CONFIG_OUTPUT_DIR_PLACEHOLDER/rf2_pdb_analysis.csv"
    RF2_DIR="CONFIG_OUTPUT_DIR_PLACEHOLDER/rf2"
    
    if [ ! -d "$RF2_DIR" ]; then
        echo "ERROR: RF2 directory not found: $RF2_DIR"
        exit 1
    fi
    
    rf2_count=$(ls ${RF2_DIR}/*.pdb 2>/dev/null | wc -l)
    
    if [ $rf2_count -eq 0 ]; then
        echo "WARNING: No RF2 structures found, skipping analysis"
    else
        echo "Analyzing ${rf2_count} RF2 structures..."
        
        run_with_glibc "${RFANTIBODY_PYTHON_BIN}" "PROJECT_DIR_PLACEHOLDER/analyze_rf2_pdb.py" \
            --rf2-dir "${RF2_DIR}" \
            --output "${ANALYSIS_CSV}"
        
        if [ $? -eq 0 ]; then
            csv_rows=$(wc -l < "${ANALYSIS_CSV}")
            echo "✓ Analysis completed: ${ANALYSIS_CSV}"
            echo "  CSV rows: ${csv_rows} (including header)"
        else
            echo "⚠ Analysis completed with warnings"
        fi
    fi
fi

echo ""
echo "=========================================="
echo "Resume Job Completed: $(date)"
echo "=========================================="
echo "Results location: CONFIG_OUTPUT_DIR_PLACEHOLDER/"
if [ -d "CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion_raw" ]; then
    raw_count=$(ls CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion_raw/*.pdb 2>/dev/null | wc -l)
    echo "  - rfdiffusion_raw/ : ${raw_count} original designs"
fi
if [ -d "CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion" ]; then
    filt_count=$(ls CONFIG_OUTPUT_DIR_PLACEHOLDER/rfdiffusion/*.pdb 2>/dev/null | wc -l)
    echo "  - rfdiffusion/     : ${filt_count} filtered designs"
fi
if [ -d "CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn" ]; then
    mpnn_count=$(ls CONFIG_OUTPUT_DIR_PLACEHOLDER/proteinmpnn/*.pdb 2>/dev/null | wc -l)
    echo "  - proteinmpnn/     : ${mpnn_count} sequences"
fi
if [ -d "CONFIG_OUTPUT_DIR_PLACEHOLDER/rf2" ]; then
    rf2_count=$(ls CONFIG_OUTPUT_DIR_PLACEHOLDER/rf2/*.pdb 2>/dev/null | wc -l)
    echo "  - rf2/             : ${rf2_count} final structures"
fi
if [ -f "CONFIG_OUTPUT_DIR_PLACEHOLDER/rf2_pdb_analysis.csv" ]; then
    echo "  - rf2_pdb_analysis.csv : Quality metrics"
fi
echo "=========================================="
EOF

    # Replace placeholders
    sed -i "s@NAME_PLACEHOLDER@${NAME}@g" "$JOB_SCRIPT"
    sed -i "s@START_STEP_PLACEHOLDER@${START_STEP}@g" "$JOB_SCRIPT"
    sed -i "s@CONFIG_OUTPUT_DIR_PLACEHOLDER@${CONFIG_OUTPUT_DIR}@g" "$JOB_SCRIPT"
    sed -i "s@BATCH_LOG_DIR_PLACEHOLDER@${RESUME_LOG_DIR}@g" "$JOB_SCRIPT"
    sed -i "s@PROJECT_DIR_PLACEHOLDER@${PROJECT_DIR}@g" "$JOB_SCRIPT"
    sed -i "s@WEIGHTS_DIR_PLACEHOLDER@${WEIGHTS_DIR}@g" "$JOB_SCRIPT"
    sed -i "s@TEMP_CONFIG_PLACEHOLDER@${TEMP_CONFIG}@g" "$JOB_SCRIPT"
    
    # Submit job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")
    JOB_IDS+=("$JOB_ID")
    JOB_NAMES+=("$NAME")
    
    echo "  ✓ Submitted SLURM job: $JOB_ID"
    echo ""
done

# ============================================================================
# Summary
# ============================================================================
echo "=========================================="
echo "Resume Batch Submission Complete"
echo "=========================================="
echo ""
echo "Submitted jobs:"
for i in "${!JOB_NAMES[@]}"; do
    echo "  [${JOB_IDS[$i]}] ${JOB_NAMES[$i]}"
done
echo ""
echo "Log directory: $RESUME_LOG_DIR"
echo ""
echo "To monitor jobs:"
echo "  squeue -u \$USER"
echo ""
echo "To check specific job:"
for i in "${!JOB_NAMES[@]}"; do
    echo "  squeue -j ${JOB_IDS[$i]}  # ${JOB_NAMES[$i]}"
done
echo ""
echo "To cancel jobs:"
echo "  scancel ${JOB_IDS[@]}"
echo ""
echo "=========================================="
