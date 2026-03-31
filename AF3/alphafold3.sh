#!/bin/bash
#SBATCH -J AF3
#SBATCH --partition=quick
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node20,node21,node24,node26
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=./logs/af3_%j_%Y%m%d_%H%M%S.out
#SBATCH --error=./logs/af3_%j.err

set -euo pipefail

date

echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: ${SLURM_JOB_NODELIST:-N/A}"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-N/A}"
echo "Job directory: $(pwd)"
echo "CPUS_used: ${SLURM_CPUS_PER_TASK:-1}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Minimal interface: follow AF3v3-template.sh style.
# Usage: sbatch alphafold3.sh /path/to/input.json [/path/to/output_dir]
usage() {
    cat <<'USAGE'
Usage: sbatch alphafold3.sh /path/to/input.json [/path/to/output_dir]

This script strictly follows the AF3v3 template style: it activates AF3v3,
cds to the system tools directory and runs run_alphafold.py with --json_path
and --output_dir only. If OUTPUT_DIR is omitted, outputs will be written to
<json_parent>/outputs.
USAGE
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    usage
    exit 0
fi

if [ "$#" -lt 1 ]; then
    echo "ERROR: missing JSON_PATH" >&2
    usage
    exit 2
fi

JSON_PATH="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
OUTPUT_PATH="${2:-$(dirname "$JSON_PATH")/outputs}"

echo "Using JSON_PATH=$JSON_PATH"
echo "Using OUTPUT_PATH=$OUTPUT_PATH"

mkdir -p "$OUTPUT_PATH" ./logs

# Strictly require AF3v3 conda environment. If activation fails, exit with error.
echo "Activating conda environment AF3v3..."
if source /public/software/apps/conda/2024.06-1/bin/activate AF3v3 2>/dev/null; then
    echo "Activated AF3v3"
else
    echo "ERROR: could not activate AF3v3 conda environment. Please ensure AF3v3 is installed and accessible." >&2
    exit 2
fi

PYTHON_EXEC=$(which python || true)
if [ -z "${PYTHON_EXEC}" ]; then
    echo "ERROR: python not found in AF3v3 environment (activation succeeded but no python)." >&2
    exit 3
fi

export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

# PYTHON_EXEC set by AF3v3 activation
echo "Using python: ${PYTHON_EXEC}"

# Run Alphafold3 from the system alphafold3 tool directory (template style).
TOOLS_DIR="/public/software/apps/tools/alphafold3"
if [ ! -d "${TOOLS_DIR}" ]; then
    echo "ERROR: tools dir not found: ${TOOLS_DIR}" >&2
    exit 1
fi

echo "Running Alphafold3 using tool dir: ${TOOLS_DIR}"
cd "${TOOLS_DIR}" || exit 1

echo "Preparing to run: run_alphafold.py --json_path=${JSON_PATH} --output_dir=${OUTPUT_PATH}"

# Some clusters provide a newer glibc in a shared location; compiled C extensions
# in the AF3 python package may require GLIBC_2.29+. If such a loader exists,
# run python via that loader and include the loader lib dir and the conda env lib
# in --library-path (see AF3v3-example.sh). This avoids system-wide glibc upgrades.
GLIBC_LOADER="/public/software/mathlib/glibc/2.34/lib/ld-linux-x86-64.so.2"
CONDA_PREFIX="$(dirname "$(dirname "${PYTHON_EXEC}")")"

if [ -x "${GLIBC_LOADER}" ]; then
    echo "Using GLIBC loader: ${GLIBC_LOADER}"
    LIB_PATHS="/public/software/mathlib/glibc/2.34/lib:${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
    echo "Library path: ${LIB_PATHS}"
    time CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0} "${GLIBC_LOADER}" --library-path "${LIB_PATHS}" "${PYTHON_EXEC}" run_alphafold.py \
        --json_path="${JSON_PATH}" \
        --output_dir="${OUTPUT_PATH}"
else
    echo "GLIBC loader not found at ${GLIBC_LOADER}; invoking python directly (may fail if system glibc is too old)"
    time "${PYTHON_EXEC}" run_alphafold.py \
        --json_path="${JSON_PATH}" \
        --output_dir="${OUTPUT_PATH}"
fi

echo "End time: $(date)"

# ============================================================================
# Post-processing: Rename best model and compute RMSD
# ============================================================================
echo ""
echo "Starting post-processing..."

# Extract job name from JSON path (e.g., pHLA_34_2 from /path/to/jobs/pHLA_34_2/example.json)
JOB_DIR="$(dirname "$JSON_PATH")"
JOB_NAME="$(basename "$JOB_DIR")"

echo "Job directory: $JOB_DIR"
echo "Job name: $JOB_NAME"

# Find the best model CIF file in output directory
BEST_MODEL="${OUTPUT_PATH}/Antibody-pHLA_Complex/Antibody-pHLA_Complex_model.cif"
if [ -f "$BEST_MODEL" ]; then
    # Rename to AF_<job_name>.cif in job directory
    RENAMED_MODEL="${JOB_DIR}/AF_${JOB_NAME}.cif"
    cp "$BEST_MODEL" "$RENAMED_MODEL"
    echo "Copied best model to: $RENAMED_MODEL"
else
    echo "Warning: Best model not found at $BEST_MODEL"
fi

# Find corresponding RF2 PDB file in job directory
RF2_PDB=$(find "$JOB_DIR" -maxdepth 1 -name "*_best.pdb" | head -1)
if [ -z "$RF2_PDB" ]; then
    echo "Warning: No RF2 *_best.pdb file found in $JOB_DIR"
    echo "Skipping RMSD calculation"
else
    echo "Found RF2 reference: $RF2_PDB"
    
    # Find ranking CSV
    RANKING_CSV="${OUTPUT_PATH}/Antibody-pHLA_Complex/Antibody-pHLA_Complex_ranking_scores.csv"
    if [ ! -f "$RANKING_CSV" ]; then
        echo "Warning: Ranking CSV not found at $RANKING_CSV"
        echo "Skipping RMSD calculation"
    else
        # Compute RMSD using our script
        # Try multiple possible locations for compute_rmsd.py
        COMPUTE_RMSD_SCRIPT=""
        for possible_path in \
            "${JOB_DIR}/../compute_rmsd.py" \
            "/public/home/xuziyi/AF3/compute_rmsd.py" \
            "$(pwd)/compute_rmsd.py" \
            "$(dirname "$JSON_PATH")/../compute_rmsd.py"; do
            if [ -f "$possible_path" ]; then
                COMPUTE_RMSD_SCRIPT="$possible_path"
                break
            fi
        done
        
        if [ -n "$COMPUTE_RMSD_SCRIPT" ] && [ -f "$COMPUTE_RMSD_SCRIPT" ] && [ -f "$BEST_MODEL" ]; then
            echo "Computing antibody RMSD for all seed-samples..."
            echo "  AF3 output dir: ${OUTPUT_PATH}/Antibody-pHLA_Complex"
            echo "  RF2 PDB: $RF2_PDB"
            echo "  Ranking CSV: $RANKING_CSV"
            
            # Use the same Python from AF3 environment
            # New signature: compute_rmsd.py <af3_output_dir> <rf2_pdb> <ranking_csv>
            "${PYTHON_EXEC}" "$COMPUTE_RMSD_SCRIPT" \
                "${OUTPUT_PATH}/Antibody-pHLA_Complex" \
                "$RF2_PDB" \
                "$RANKING_CSV"
            
            if [ $? -eq 0 ]; then
                echo "RMSD calculation completed successfully"
            else
                echo "Warning: RMSD calculation failed"
            fi
        else
            echo "Warning: compute_rmsd.py not found or best model missing"
            echo "Skipping RMSD calculation"
        fi
    fi
fi

echo "Post-processing completed"
echo "Final end time: $(date)"

