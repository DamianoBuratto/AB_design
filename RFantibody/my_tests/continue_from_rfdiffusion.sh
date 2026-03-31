#!/bin/bash

# ============================================================================
# Continue Pipeline from Existing RFdiffusion Results
# ============================================================================
# Usage:
#   bash continue_from_rfdiffusion.sh <output_dir> [start_step]
#
# Arguments:
#   output_dir   - Directory containing rfdiffusion/ folder with existing results
#   start_step   - Optional: 2=Filtering, 3=ProteinMPNN, 4=RF2, 5=Analysis (default: 2)
#
# Example:
#   # Continue from filtering step
#   bash continue_from_rfdiffusion.sh /public/home/xuziyi/outputs/region1
#
#   # Continue from ProteinMPNN step (skip filtering)
#   bash continue_from_rfdiffusion.sh /public/home/xuziyi/outputs/region1 3
# ============================================================================

set -e

if [ $# -lt 1 ]; then
    echo "ERROR: Output directory required"
    echo "Usage: bash continue_from_rfdiffusion.sh <output_dir> [start_step]"
    echo ""
    echo "Example:"
    echo "  bash continue_from_rfdiffusion.sh /public/home/xuziyi/outputs/region1"
    echo "  bash continue_from_rfdiffusion.sh /public/home/xuziyi/outputs/region1 3"
    exit 1
fi

OUTPUT_DIR="$1"
START_STEP="${2:-2}"  # Default to step 2 (filtering)

# Validate output directory
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "ERROR: Output directory not found: $OUTPUT_DIR"
    exit 1
fi

if [ ! -d "$OUTPUT_DIR/rfdiffusion" ]; then
    echo "ERROR: rfdiffusion subdirectory not found in: $OUTPUT_DIR"
    echo "Expected: $OUTPUT_DIR/rfdiffusion/"
    exit 1
fi

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/run_config.txt"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

source "$CONFIG_FILE"

# Determine configuration name from output directory
CONFIG_NAME=$(basename "$OUTPUT_DIR")

echo "=========================================="
echo "Continue Pipeline from RFdiffusion Results"
echo "=========================================="
echo "Configuration: $CONFIG_NAME"
echo "Output directory: $OUTPUT_DIR"
echo "Start from step: $START_STEP"
echo "Project directory: $PROJECT_DIR"
echo ""

# Count existing designs
design_count=$(ls $OUTPUT_DIR/rfdiffusion/design_*.pdb 2>/dev/null | wc -l)
echo "Found $design_count existing RFdiffusion designs"
echo ""

# Create job script
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
JOB_SCRIPT="/tmp/continue_${CONFIG_NAME}_${TIMESTAMP}.slurm"

cat > "$JOB_SCRIPT" <<'EOFSCRIPT'
#!/bin/bash
#SBATCH --job-name=Continue_CONF_NAME
#SBATCH --output=OUT_DIR/continue_%j.out
#SBATCH --error=OUT_DIR/continue_%j.err
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node20,node21,node24,node26,node27,node28
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00

START_STEP=START_STEP_VAL
OUTPUT_DIR="OUT_DIR"
PROJECT_DIR="PROJ_DIR"
WEIGHTS_DIR="WEIGHTS_DIR_VAL"

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

# Setup environment
export PYTHONPATH=${PROJECT_DIR}/src:${PROJECT_DIR}/src/rfantibody/rfdiffusion:${PROJECT_DIR}/src/rfantibody:${PROJECT_DIR}/include/SE3Transformer:$PYTHONPATH
export PATH=${PROJECT_DIR}/include/USalign:$PATH
export WEIGHTS_DIR
export ICECREAM_COLORS=never
export PYTHONUNBUFFERED=1
export NO_COLOR=1
export TERM=dumb

cd ${PROJECT_DIR} || exit 1

echo "=========================================="
echo "Continuing Pipeline: CONF_NAME"
echo "Started: $(date)"
echo "=========================================="
echo "Start from step: ${START_STEP}"
echo "Output: ${OUTPUT_DIR}"
echo "=========================================="
echo ""

# Step 2: Filtering
if [ ${START_STEP} -le 2 ]; then
    echo "==> Step 2: Filtering designs by RHC-HIS contacts (<6A)"
    echo ""
    
    ${RFANTIBODY_PYTHON_BIN} "${PROJECT_DIR}/count_contacts_with_hotspot.py" \
        --pdb-dir "${OUTPUT_DIR}/rfdiffusion" \
        --workdir "${OUTPUT_DIR}" \
        --mode backbone \
        --primary-cutoff 6 \
        --cutoffs 6 7 8
    
    if [ $? -ne 0 ]; then
        echo "✗ Contact counting failed"
        exit 1
    fi
    
    # Extract filtered PDBs
    CSV_FILE="${OUTPUT_DIR}/contacts_with_hotspot_summary.csv"
    MAPFILE=$(mktemp)
    ${RFANTIBODY_PYTHON_BIN} -c "
import csv
with open('${CSV_FILE}') as f:
    reader = csv.DictReader(f)
    for row in reader:
        val = row.get('(<= 6A)', '').strip()
        try:
            if float(val) > 0:
                print(row.get('PDB File', ''))
        except ValueError:
            pass
" > "${MAPFILE}"
    
    MATCHING_PDBS=()
    while IFS= read -r line; do
        MATCHING_PDBS+=("$line")
    done < "${MAPFILE}"
    rm -f "${MAPFILE}"
    
    echo "Found ${#MATCHING_PDBS[@]} designs with contacts"
    
    if [ ${#MATCHING_PDBS[@]} -eq 0 ]; then
        echo "✗ No designs passed filtering"
        exit 0
    fi
    
    # Move to rfdiffusion_raw if not already done
    if [ ! -d "${OUTPUT_DIR}/rfdiffusion_raw" ]; then
        mv "${OUTPUT_DIR}/rfdiffusion" "${OUTPUT_DIR}/rfdiffusion_raw"
        mkdir -p "${OUTPUT_DIR}/rfdiffusion"
        
        for pdb in "${MATCHING_PDBS[@]}"; do
            cp "${OUTPUT_DIR}/rfdiffusion_raw/${pdb}" "${OUTPUT_DIR}/rfdiffusion/"
        done
    fi
    
    filtered_count=$(ls ${OUTPUT_DIR}/rfdiffusion/design_*.pdb 2>/dev/null | wc -l)
    echo "✓ Filtering complete: ${filtered_count} designs retained"
    echo ""
fi

# Step 3: ProteinMPNN
if [ ${START_STEP} -le 3 ]; then
    echo "==> Step 3: Running ProteinMPNN"
    echo ""
    
    bash ${PROJECT_DIR}/run_antibody_design.sh --proteinmpnn-only --output-dir ${OUTPUT_DIR} </dev/null
    
    if [ $? -ne 0 ]; then
        echo "✗ ProteinMPNN failed"
        exit 1
    fi
    
    mpnn_count=$(ls ${OUTPUT_DIR}/proteinmpnn/design_*.fa 2>/dev/null | wc -l)
    echo "✓ ProteinMPNN completed: ${mpnn_count} files"
    echo ""
fi

# Step 4: RF2
if [ ${START_STEP} -le 4 ]; then
    echo "==> Step 4: Running RF2"
    echo ""
    
    bash ${PROJECT_DIR}/run_antibody_design.sh --rf2-only --output-dir ${OUTPUT_DIR} </dev/null
    
    if [ $? -ne 0 ]; then
        echo "✗ RF2 failed"
        exit 1
    fi
    
    rf2_count=$(ls ${OUTPUT_DIR}/rf2/design_*.pdb 2>/dev/null | wc -l)
    echo "✓ RF2 completed: ${rf2_count} structures"
    echo ""
fi

# Step 5: Analysis
if [ ${START_STEP} -le 5 ]; then
    echo "==> Step 5: Analyzing RF2 structures"
    echo ""
    
    ${RFANTIBODY_PYTHON_BIN} "${PROJECT_DIR}/analyze_rf2_pdb.py" \
        --rf2-dir "${OUTPUT_DIR}/rf2" \
        --output "${OUTPUT_DIR}/rf2_pdb_analysis.csv"
    
    echo "✓ Analysis complete"
fi

echo ""
echo "=========================================="
echo "Pipeline Completed: $(date)"
echo "=========================================="
EOFSCRIPT

# Replace placeholders
sed -i "s|CONF_NAME|${CONFIG_NAME}|g" "$JOB_SCRIPT"
sed -i "s|OUT_DIR|${OUTPUT_DIR}|g" "$JOB_SCRIPT"
sed -i "s|PROJ_DIR|${PROJECT_DIR}|g" "$JOB_SCRIPT"
sed -i "s|START_STEP_VAL|${START_STEP}|g" "$JOB_SCRIPT"
sed -i "s|WEIGHTS_DIR_VAL|${WEIGHTS_DIR}|g" "$JOB_SCRIPT"

echo "Generated SLURM job script: $JOB_SCRIPT"
echo ""
echo "To submit the job:"
echo "  sbatch $JOB_SCRIPT"
echo ""
echo "Or run directly (non-SLURM):"
echo "  bash $JOB_SCRIPT"
echo ""

# Ask if user wants to submit
read -p "Submit job now? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")
    echo "✓ Job submitted: $JOB_ID"
    echo ""
    echo "Monitor with:"
    echo "  squeue -j $JOB_ID"
    echo "  tail -f ${OUTPUT_DIR}/continue_${JOB_ID}.out"
fi
