#!/bin/bash

# ============================================================================
# Hotspot Batch Testing Script for HPC
# ============================================================================
# Usage:
#   bash run_hotspot_batch.sh                    # Run complete pipeline
#   bash run_hotspot_batch.sh --skip-rfdiffusion # Skip RFdiffusion, use existing results
#   bash run_hotspot_batch.sh --from-step 2      # Start from step 2 (filtering)
#
# This script submits multiple SLURM jobs to test different hotspot combinations
# Each combination gets its own output directory and job
#
# Edit the HOTSPOT_CONFIGS array below to define your test cases
# ============================================================================

set -e

# ============================================================================
# Parse command line arguments
# ============================================================================
SKIP_RFDIFFUSION=false
START_STEP=1
SELECTED_CONFIGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-rfdiffusion)
            SKIP_RFDIFFUSION=true
            START_STEP=2
            shift
            ;;
        --from-step)
            START_STEP="$2"
            if [ "$START_STEP" -ge 2 ]; then
                SKIP_RFDIFFUSION=true
            fi
            shift 2
            ;;
        --config)
            SELECTED_CONFIGS+=("$2")
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: bash run_hotspot_batch.sh [--skip-rfdiffusion] [--from-step N] [--config NAME]"
            echo "  --config NAME: Run only specified configuration (can be used multiple times)"
            echo "  --from-step N: Start from step N (1=RFdiffusion, 2=Filtering, 3=ProteinMPNN, 4=RF2, 5=Analysis)"
            echo "  --skip-rfdiffusion: Skip RFdiffusion step (same as --from-step 2)"
            echo ""
            echo "Examples:"
            echo "  bash run_hotspot_batch.sh                    # Run all configs from step 1"
            echo "  bash run_hotspot_batch.sh --from-step 2      # Run all configs from step 2"
            echo "  bash run_hotspot_batch.sh --config region2   # Run only region2 from step 1"
            echo "  bash run_hotspot_batch.sh --config region2 --config combined --from-step 3  # Run region2 and combined from step 3"
            exit 1
            ;;
    esac
done

# ============================================================================
# Configuration: Define your test cases here
# ============================================================================
# 🔧 BATCH TEST MODE: Choose what to vary
# - "HOTSPOTS": Vary hotspots, fix CDR loops
# - "CDR_LOOPS": Vary CDR loops, fix hotspots
# - "MIXED": Explicitly specify both hotspots and CDR loops for each config
#
BATCH_MODE="MIXED"  # Options: "HOTSPOTS", "CDR_LOOPS", or "MIXED"

# 🎯 Fixed parameters (used when BATCH_MODE is HOTSPOTS or CDR_LOOPS)
FIXED_HOTSPOTS="[T82,T85,T155,T156,T159]"
FIXED_CDR_LOOPS="[L1:10-13,L2:7,L3:9-11,H1:7,H2:6,H3:11-14]"
FIXED_T=200  # RFdiffusion T parameter (diffusion steps)
FIXED_NUM_DESIGNS=700  # Number of designs per configuration

# Note: Each job will generate its own random seed for non-deterministic runs

# ============================================================================
# Test Configurations
# ============================================================================
# Format depends on BATCH_MODE:
#
# When BATCH_MODE="HOTSPOTS" or "CDR_LOOPS":
#   "NAME|VARIABLE_VALUE|NUM_DESIGNS|T_VALUE"
#   VARIABLE_VALUE: The varying parameter (hotspots or CDR loops)
#
# When BATCH_MODE="MIXED":
#   "NAME|CDR_LOOPS|HOTSPOTS|NUM_DESIGNS|T_VALUE"
#   Both CDR_LOOPS and HOTSPOTS must be explicitly specified
#
# NUM_DESIGNS and T_VALUE are optional and default to FIXED_NUM_DESIGNS and FIXED_T
# NOTE: Use | as delimiter to avoid conflict with : in CDR loop notation
#
# Examples:
#   MIXED mode: "optimal_small|[L1:10-12,L2:7,L3:9-11,H1:7,H2:6,H3:11-13]|[T82,T85,T155,T156,T159]|100|200"

declare -a ALL_CONFIGS=(
    # ============================================================================
    # Parallel Runs with Identical Parameters (Random Seeds)
    # All runs use: Hotspots=[T82,T85,T155,T156,T159]
    #               CDR_LOOPS=[L1:10-13,L2:7,L3:9-11,H1:7,H2:6,H3:11-14]
    #               T=200, 700 designs, ProteinMPNN=32 seqs
    #               Deterministic=False (each run uses different random seed)
    # Goal: Explore stochastic variation with identical input parameters
    # Naming: 200T_700_32_L1_10-13_H3_11-14_T82-159_replica[1-3]
    #   200T = T parameter (200 diffusion steps)
    #   700 = Number of designs
    #   32 = ProteinMPNN sequences per design
    #   L1_10-13 = L1 loop range
    #   H3_11-14 = H3 loop range
    #   T82-159 = Hotspots (T82,T85,T155,T156,T159)
    #   replica[1-3] = Independent replicate number
    # ============================================================================
    
    # Replica 1: Random seed will be generated at runtime
    "200T_700_32_L1_10-13_H3_11-14_T82-159_replica1|[L1:10-13,L2:7,L3:9-11,H1:7,H2:6,H3:11-14]|[T82,T85,T155,T156,T159]|700|200"
    
    # Replica 2: Random seed will be generated at runtime
    "200T_700_32_L1_10-13_H3_11-14_T82-159_replica2|[L1:10-13,L2:7,L3:9-11,H1:7,H2:6,H3:11-14]|[T82,T85,T155,T156,T159]|700|200"
    
    # Replica 3: Random seed will be generated at runtime
    "200T_700_32_L1_10-13_H3_11-14_T82-159_replica3|[L1:10-13,L2:7,L3:9-11,H1:7,H2:6,H3:11-14]|[T82,T85,T155,T156,T159]|700|200"
)

# ============================================================================
# Select configurations to run
# ============================================================================
if [ ${#SELECTED_CONFIGS[@]} -gt 0 ]; then
    # Filter configurations based on selection
    declare -a HOTSPOT_CONFIGS=()
    for config in "${ALL_CONFIGS[@]}"; do
        config_name=$(echo "$config" | cut -d'|' -f1)
        if [[ " ${SELECTED_CONFIGS[@]} " =~ " ${config_name} " ]]; then
            HOTSPOT_CONFIGS+=("$config")
        fi
    done

    if [ ${#HOTSPOT_CONFIGS[@]} -eq 0 ]; then
        echo "ERROR: No matching configurations found for: ${SELECTED_CONFIGS[*]}"
        echo "Available configurations: $(echo "${ALL_CONFIGS[@]}" | sed 's/|[^|]*//g' | tr ' ' '\n' | sort | tr '\n' ' ')"
        exit 1
    fi
else
    # Run all configurations
    HOTSPOT_CONFIGS=("${ALL_CONFIGS[@]}")
fi

# ============================================================================
# Load configuration from run_config.txt
# ============================================================================
# Determine script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/run_config.txt"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    echo "Please ensure run_config.txt exists in the same directory as this script"
    exit 1
fi

echo "Loading configuration from: $CONFIG_FILE"
source "$CONFIG_FILE"

# Validate PROJECT_DIR was loaded
if [ -z "$PROJECT_DIR" ]; then
    echo "ERROR: PROJECT_DIR not set in $CONFIG_FILE"
    exit 1
fi

if [ ! -d "$PROJECT_DIR" ]; then
    echo "ERROR: PROJECT_DIR does not exist: $PROJECT_DIR"
    exit 1
fi

echo "Project directory: $PROJECT_DIR"

# Backup original run_config.txt to restore after each configuration
ORIGINAL_CONFIG_BACKUP="${PROJECT_DIR}/run_config.txt.original_backup"
cp "$CONFIG_FILE" "$ORIGINAL_CONFIG_BACKUP"

# CDR loops configuration (from your test_my_input.slurm)
CDR_LOOPS="[L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13]"

# Validate input files exist
if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: Input PDB not found: $INPUT_PDB"
    exit 1
fi

if [ ! -f "$FRAMEWORK_PDB" ]; then
    echo "ERROR: Framework PDB not found: $FRAMEWORK_PDB"
    exit 1
fi

# Create batch ID and batch logs directory
BATCH_ID=$(date +%Y%m%d_%H%M%S)
# Batch log dir includes the batch id so outputs can be grouped per submission
BATCH_LOG_DIR="${PROJECT_DIR}/logs/hotspot_batch_${BATCH_ID}"
mkdir -p "$BATCH_LOG_DIR"

echo "=========================================="
echo "RFantibody Parallel Runs - Random Seeds"
echo "=========================================="
echo "📋 BATCH MODE: ${BATCH_MODE}"
echo "🎲 Seed Strategy: Random (each job generates unique seed)"
echo ""
echo "🔧 Fixed Parameters:"
if [ "$BATCH_MODE" = "CDR_LOOPS" ]; then
    echo "  Hotspots:    ${FIXED_HOTSPOTS}"
    echo "  T:           ${FIXED_T}"
    echo "  Num Designs: ${FIXED_NUM_DESIGNS}"
    echo ""
    echo "🔄 Varying: CDR Loops"
elif [ "$BATCH_MODE" = "HOTSPOTS" ]; then
    echo "  CDR Loops:   ${FIXED_CDR_LOOPS}"
    echo "  T:           ${FIXED_T}"
    echo "  Num Designs: ${FIXED_NUM_DESIGNS}"
    echo ""
    echo "🔄 Varying: Hotspots"
else
    echo "  Mode: MIXED (both CDR loops and hotspots specified per config)"
    echo "  Default T:           ${FIXED_T}"
    echo "  Default Num Designs: ${FIXED_NUM_DESIGNS}"
    echo ""
    echo "🔄 Varying: Both CDR Loops and Hotspots"
fi
echo ""
echo "Configurations to test: ${#HOTSPOT_CONFIGS[@]}"
echo "Start from step: $START_STEP"
if [ "$SKIP_RFDIFFUSION" = true ]; then
    echo "Mode: Using existing RFdiffusion results"
else
    echo "Mode: Running complete pipeline"
fi
echo "Batch log directory: $BATCH_LOG_DIR"
echo ""

# ============================================================================
# Job tracking
# ============================================================================
declare -a JOB_IDS=()
declare -a JOB_NAMES=()

# ============================================================================
# Submit SLURM jobs for each configuration
# ============================================================================
for config in "${HOTSPOT_CONFIGS[@]}"; do
    # Parse configuration based on BATCH_MODE
    if [ "$BATCH_MODE" = "MIXED" ]; then
        # Format: NAME|CDR_LOOPS|HOTSPOTS|NUM_DESIGNS|T_VALUE
        IFS='|' read -r NAME CDR_LOOPS HOTSPOTS CONFIG_NUM_DESIGNS CONFIG_T_VALUE <<< "$config"
    else
        # Format: NAME|VARIABLE_VALUE|NUM_DESIGNS|T_VALUE
        IFS='|' read -r NAME VARIABLE_VALUE CONFIG_NUM_DESIGNS CONFIG_T_VALUE <<< "$config"
        
        # Determine actual hotspots and CDR loops based on batch mode
        if [ "$BATCH_MODE" = "CDR_LOOPS" ]; then
            HOTSPOTS="$FIXED_HOTSPOTS"
            CDR_LOOPS="$VARIABLE_VALUE"
        else
            HOTSPOTS="$VARIABLE_VALUE"
            CDR_LOOPS="$FIXED_CDR_LOOPS"
        fi
    fi
    
    # Set defaults if not specified
    NUM_DESIGNS="${CONFIG_NUM_DESIGNS:-$FIXED_NUM_DESIGNS}"
    T_VALUE="${CONFIG_T_VALUE:-$FIXED_T}"
    
    echo "----------------------------------------"
    echo "Configuration: $NAME"
    echo "  Hotspots:  $HOTSPOTS"
    echo "  CDR Loops: $CDR_LOOPS"
    echo "  T:         $T_VALUE"
    echo "  Designs:   $NUM_DESIGNS"
    
    # Create unique output directory for this configuration
    BASE_OUTPUT_DIR="${PROJECT_DIR}/outputs"
    CONFIG_OUTPUT_DIR="${BASE_OUTPUT_DIR}/${NAME}"
    echo "  Output: $CONFIG_OUTPUT_DIR"
    
    # Create output directories
    mkdir -p "${CONFIG_OUTPUT_DIR}"
    
    # Create a temporary config file for this hotspot configuration
    TEMP_CONFIG="${BATCH_LOG_DIR}/${NAME}_run_config.txt"
    cat > "$TEMP_CONFIG" <<CFGEOF
INPUT_PDB=${INPUT_PDB}
FRAMEWORK_PDB=${FRAMEWORK_PDB}
OUTPUT_DIR=${CONFIG_OUTPUT_DIR}
WEIGHTS_DIR=${WEIGHTS_DIR}
PROJECT_DIR=${PROJECT_DIR}
CFGEOF
    
    # Create SLURM job script for this configuration
    JOB_SCRIPT="${BATCH_LOG_DIR}/${NAME}_job.slurm"
    
    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=RFantibody_${NAME}
#SBATCH --output=${BATCH_LOG_DIR}/${NAME}_%j.out
#SBATCH --error=${BATCH_LOG_DIR}/${NAME}_%j.err
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node20,node21,node24,node26,node27,node28
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00

# Pipeline start step (1=RFdiffusion, 2=Filtering, 3=ProteinMPNN, 4=RF2, 5=Analysis)
START_STEP=${START_STEP}

source /public/home/xuziyi/miniconda/etc/profile.d/conda.sh
conda activate RFantibody_HPC

source /public/software/profile.d/compiler_intel-2021.3.0.sh
source /public/software/profile.d/cuda12.2-env.sh
source /public/software/profile.d/fftw3_3.3.10_float-sse2-avx2.sh

# Setup GLIBC environment
GLIBC_ROOT=/public/software/mathlib/glibc/2.34
GLIBC_LOADER=\${GLIBC_ROOT}/lib/ld-linux-x86-64.so.2
GLIBC_LIB_PATH=\${GLIBC_ROOT}/lib:\${GLIBC_ROOT}/lib64

CONDA_LIB_DIR="\${CONDA_PREFIX}/lib"
if [ ! -d "\${CONDA_LIB_DIR}" ]; then
    echo "Expected conda lib directory not found: \${CONDA_LIB_DIR}" >&2
    exit 1
fi

GLIBC_RUNTIME_PATH="\${GLIBC_LIB_PATH}:\${CONDA_LIB_DIR}"
if [ -n "\${LD_LIBRARY_PATH:-}" ]; then
    GLIBC_RUNTIME_PATH="\${GLIBC_RUNTIME_PATH}:\${LD_LIBRARY_PATH}"
fi

export GLIBC_ROOT
export GLIBC_LOADER
export RFANTIBODY_GLIBC_LIBRARY_PATH="\${GLIBC_RUNTIME_PATH}"
export RFANTIBODY_PYTHON_BIN="\${CONDA_PREFIX}/bin/python"

if [ ! -x "\${RFANTIBODY_PYTHON_BIN}" ]; then
    echo "Python interpreter not found at \${RFANTIBODY_PYTHON_BIN}" >&2
    exit 1
fi

if [ ! -x "\${GLIBC_LOADER}" ]; then
    echo "Custom glibc loader not found at \${GLIBC_LOADER}" >&2
    exit 1
fi

run_with_glibc() {
    if [ "\$#" -eq 0 ]; then
        echo "run_with_glibc: missing command" >&2
        return 1
    fi
    if [ -z "\${GLIBC_LOADER:-}" ]; then
        echo "run_with_glibc: GLIBC_LOADER is not set" >&2
        return 1
    fi
    if [ ! -x "\${GLIBC_LOADER}" ]; then
        echo "run_with_glibc: loader not executable: \${GLIBC_LOADER}" >&2
        return 1
    fi
    local cmd="\$1"
    shift
    if [[ "\${cmd}" != /* ]]; then
        cmd=\$(command -v "\${cmd}")
    fi
    if [ -z "\${cmd}" ]; then
        echo "run_with_glibc: command not found" >&2
        return 1
    fi
    "\${GLIBC_LOADER}" --library-path "\${RFANTIBODY_GLIBC_LIBRARY_PATH}" "\${cmd}" "\$@"
}

export -f run_with_glibc

# Generate random seed for this job
JOB_SEED=\$RANDOM
echo "🎲 Generated random seed for this job: \${JOB_SEED}"
echo ""

# Setup environment
export PYTHONPATH=${PROJECT_DIR}/src:${PROJECT_DIR}/src/rfantibody/rfdiffusion:${PROJECT_DIR}/src/rfantibody:${PROJECT_DIR}/include/SE3Transformer:\$PYTHONPATH
export PATH=${PROJECT_DIR}/include/USalign:\$PATH
export WEIGHTS_DIR=${WEIGHTS_DIR}
export ICECREAM_COLORS=never
export HYDRA_FULL_ERROR=1

# Set random seed environment variables (for reproducibility tracking)
export PYTHONHASHSEED=\${JOB_SEED}

# Disable ANSI color codes
export PYTHONUNBUFFERED=1
export NO_COLOR=1
export TERM=dumb

# Change to project directory
cd ${PROJECT_DIR} || {
    echo "ERROR: Failed to change to project directory: ${PROJECT_DIR}" >&2
    exit 1
}

echo "=========================================="
echo "Configuration: ${NAME}"
echo "Started: \$(date)"
echo "=========================================="
echo "📋 Batch Mode: ${BATCH_MODE}"
echo "🎲 Random Seed: \${JOB_SEED}"
echo ""
echo "🔧 RFdiffusion Parameters:"
echo "  Hotspots:   ${HOTSPOTS}"
echo "  CDR Loops:  ${CDR_LOOPS}"
echo "  T:          ${T_VALUE} (diffusion steps)"
echo "  Designs:    ${NUM_DESIGNS}"
echo "  Deterministic: False (random seed=\${JOB_SEED})"
echo ""
echo "🧬 ProteinMPNN:"
echo "  Sequences per design: 32"
echo ""
echo "📁 Output:"
echo "  Directory: ${CONFIG_OUTPUT_DIR}"
echo ""
echo "🔍 Filtering:"
echo "  Contact distance: <6A to RHC-HIS (ARG-HIS-CYS motif)"
echo ""
echo "⚙️  Pipeline:"
echo "  Start from step: \${START_STEP}"
echo "  Working directory: \$(pwd)"
echo "=========================================="
echo ""
echo "💡 Note: Random seed mode explores stochastic variation"
echo "         Each job uses unique seed (\${JOB_SEED}) for independent sampling"
echo "         Filtering step targets RHC motif HIS (ARG-HIS-CYS motif)"
echo "         RFdiffusion hotspot constraint: ${HOTSPOTS}"
echo ""

# Use task-specific config file to avoid race conditions between parallel jobs
# Set environment variable for run_antibody_design.sh to use
export CONFIG_FILE="${TEMP_CONFIG}"

echo "Using task-specific config: \${CONFIG_FILE}"
echo ""

# ============================================================================
# Step 1: RFdiffusion
# ============================================================================
if [ \${START_STEP} -le 1 ]; then
    echo "==> Step 1: Running RFdiffusion (${NUM_DESIGNS} designs)"
    echo ""

    # Create symlink for config (handle race condition for parallel jobs)
    mkdir -p scripts/config
    
    # Check if symlink already exists and is correct
    EXPECTED_TARGET="${PROJECT_DIR}/src/rfantibody/rfdiffusion/config/inference"
    if [ -L "scripts/config/inference" ]; then
        CURRENT_TARGET=\$(readlink "scripts/config/inference")
        if [ "\${CURRENT_TARGET}" = "\${EXPECTED_TARGET}" ]; then
            echo "Config symlink already exists and is correct"
        else
            echo "Updating config symlink to correct target"
            rm -f scripts/config/inference
            ln -sfn "\${EXPECTED_TARGET}" scripts/config/inference
        fi
    else
        echo "Creating config symlink"
        rm -rf scripts/config/inference
        ln -sfn "\${EXPECTED_TARGET}" scripts/config/inference
    fi
    
    # Verify symlink (use -e to check if target exists)
    if [ ! -e "scripts/config/inference" ]; then
        echo "ERROR: Config symlink verification failed" >&2
        echo "  Expected: scripts/config/inference -> \${EXPECTED_TARGET}" >&2
        echo "  Check if target directory exists: \${EXPECTED_TARGET}" >&2
        exit 1
    fi
    
    echo "Config symlink verified: scripts/config/inference"

    # Create output directory structure
    mkdir -p "${CONFIG_OUTPUT_DIR}/rfdiffusion"

    # Run RFdiffusion with random seed (non-deterministic mode)
    # Seed is set for reproducibility tracking but deterministic mode is disabled
    run_with_glibc "\${RFANTIBODY_PYTHON_BIN}" scripts/rfdiffusion_inference.py \\
        --config-name antibody \\
        antibody.target_pdb="${INPUT_PDB}" \\
        antibody.framework_pdb="${FRAMEWORK_PDB}" \\
        inference.ckpt_override_path="${WEIGHTS_DIR}/RFdiffusion_Ab.pt" \\
        'ppi.hotspot_res=${HOTSPOTS}' \\
        'antibody.design_loops=${CDR_LOOPS}' \\
        inference.num_designs=${NUM_DESIGNS} \\
        diffuser.T=${T_VALUE} \\
        ++inference.seed=\${JOB_SEED} \\
        inference.output_prefix="${CONFIG_OUTPUT_DIR}/rfdiffusion/design"

    RFDIFF_EXIT=\$?

    if [ \$RFDIFF_EXIT -ne 0 ]; then
        echo "✗ RFdiffusion failed with exit code \$RFDIFF_EXIT"
        exit \$RFDIFF_EXIT
    fi

    design_count=\$(ls ${CONFIG_OUTPUT_DIR}/rfdiffusion/design_*.pdb 2>/dev/null | wc -l)
    echo "✓ RFdiffusion completed: \$design_count designs"
    echo ""
else
    echo "==> Step 1: Skipping RFdiffusion (using existing results)"
    if [ ! -d "${CONFIG_OUTPUT_DIR}/rfdiffusion" ]; then
        echo "ERROR: RFdiffusion output directory not found: ${CONFIG_OUTPUT_DIR}/rfdiffusion"
        echo "Cannot skip RFdiffusion step without existing results!"
        exit 1
    fi
    design_count=\$(ls ${CONFIG_OUTPUT_DIR}/rfdiffusion/design_*.pdb 2>/dev/null | wc -l)
    echo "Found \$design_count existing designs in ${CONFIG_OUTPUT_DIR}/rfdiffusion"
    echo ""
fi

# ============================================================================
# Step 2: Filtering
# ============================================================================
if [ \${START_STEP} -le 2 ]; then
    echo "==> Step 2: Filtering designs by RHC-HIS contacts (<6A)"
    echo "            (RHC = ARG-HIS-CYS motif, reference HIS for all configs)"
    echo ""

    # Run contact counting (use simple python call without run_with_glibc)
    CONTACT_SCRIPT="${PROJECT_DIR}/count_contacts_with_hotspot.py"
    RFDIFF_OUTPUT="${CONFIG_OUTPUT_DIR}/rfdiffusion"
    FILTER_WORKDIR="${CONFIG_OUTPUT_DIR}"

    echo "Running: \${RFANTIBODY_PYTHON_BIN} \${CONTACT_SCRIPT}"
    echo "  --pdb-dir \${RFDIFF_OUTPUT}"
    echo "  --workdir \${FILTER_WORKDIR}"
    echo ""

    \${RFANTIBODY_PYTHON_BIN} "\${CONTACT_SCRIPT}" \\
        --pdb-dir "\${RFDIFF_OUTPUT}" \\
        --workdir "\${FILTER_WORKDIR}" \\
        --mode backbone \\
        --primary-cutoff 6 \\
        --cutoffs 5 6 7 8

    CONTACT_EXIT=\$?
    if [ \$CONTACT_EXIT -ne 0 ]; then
        echo "✗ Contact counting failed with exit code \$CONTACT_EXIT"
        exit 1
    fi

    # Filter and move designs
    CSV_FILE="${CONFIG_OUTPUT_DIR}/contacts_with_hotspot_summary.csv"
    if [ ! -f "\${CSV_FILE}" ]; then
        echo "ERROR: Contact summary CSV not found: \${CSV_FILE}"
        exit 1
    fi

    # Extract PDB filenames with contacts
    MAPFILE=\$(mktemp)
    CSV_PATH="\${CSV_FILE}"
    \${RFANTIBODY_PYTHON_BIN} -c "
import csv
csv_path = '\${CSV_PATH}'
primary_col = '(<= 6A)'

try:
    with open(csv_path) as handle:
        reader = csv.DictReader(handle)
        if primary_col not in reader.fieldnames:
            import sys
            print('ERROR: Primary contact column not found', file=sys.stderr)
            sys.exit(1)
        
        for row in reader:
            value = row.get(primary_col, '').strip()
            try:
                if float(value) > 0:
                    print(row.get('PDB File', row.get('Filename', '')))
            except ValueError:
                pass
except FileNotFoundError as e:
    import sys
    print(f'ERROR: CSV file not found: {csv_path}', file=sys.stderr)
    sys.exit(1)
except Exception as e:
    import sys
    print(f'ERROR: {e}', file=sys.stderr)
    sys.exit(1)
" > "\${MAPFILE}"

    CSV_FILTER_EXIT=\$?

    if [ \$CSV_FILTER_EXIT -ne 0 ]; then
        echo "✗ CSV filtering failed"
        rm -f "\${MAPFILE}"
        exit 1
    fi

    MATCHING_PDBS=()
    while IFS= read -r line; do
        MATCHING_PDBS+=("\$line")
    done < "\${MAPFILE}"
    rm -f "\${MAPFILE}"

    echo "Found \${#MATCHING_PDBS[@]} designs with hotspot contacts"

    if [ \${#MATCHING_PDBS[@]} -eq 0 ]; then
        echo "✗ No designs passed filtering. Pipeline stopped."
        exit 0
    fi

    # Move original rfdiffusion to rfdiffusion_raw (only if not already moved)
    ORIG_DIR="${CONFIG_OUTPUT_DIR}/rfdiffusion"
    RAW_DIR="${CONFIG_OUTPUT_DIR}/rfdiffusion_raw"
    NEW_DIR="${CONFIG_OUTPUT_DIR}/rfdiffusion"

    if [ ! -d "\${RAW_DIR}" ]; then
        mv "\${ORIG_DIR}" "\${RAW_DIR}"
        mkdir -p "\${NEW_DIR}"

        # Copy filtered designs
        for pdb in "\${MATCHING_PDBS[@]}"; do
            src="\${RAW_DIR}/\${pdb}"
            if [ -f "\${src}" ]; then
                cp "\${src}" "\${NEW_DIR}/"
            fi
        done
    else
        echo "Note: rfdiffusion_raw already exists, skipping move operation"
        echo "      Using existing filtered designs in rfdiffusion/"
    fi

    filtered_count=\$(ls \${NEW_DIR}/design_*.pdb 2>/dev/null | wc -l)
    echo "✓ Filtering complete: \${filtered_count} designs retained"
    echo ""
else
    echo "==> Step 2: Skipping filtering (using existing filtered results)"
    if [ ! -d "${CONFIG_OUTPUT_DIR}/rfdiffusion" ]; then
        echo "ERROR: Filtered rfdiffusion directory not found!"
        exit 1
    fi
    filtered_count=\$(ls ${CONFIG_OUTPUT_DIR}/rfdiffusion/design_*.pdb 2>/dev/null | wc -l)
    echo "Found \$filtered_count filtered designs"
    echo ""
fi

# ============================================================================
# Step 3: ProteinMPNN
# ============================================================================
if [ \${START_STEP} -le 3 ]; then
    echo "==> Step 3: Running ProteinMPNN (32 sequences per design)"
    echo ""

    # Run ProteinMPNN
    bash ${PROJECT_DIR}/run_antibody_design.sh --proteinmpnn-only --output-dir ${CONFIG_OUTPUT_DIR} </dev/null

    if [ \$? -ne 0 ]; then
        echo "✗ ProteinMPNN failed"
        exit 1
    fi

    mpnn_count=\$(ls ${CONFIG_OUTPUT_DIR}/proteinmpnn/design_*.fa 2>/dev/null | wc -l)
    echo "✓ ProteinMPNN completed: \$mpnn_count FASTA files"
    echo ""
else
    echo "==> Step 3: Skipping ProteinMPNN (using existing results)"
    mpnn_count=\$(ls ${CONFIG_OUTPUT_DIR}/proteinmpnn/design_*.fa 2>/dev/null | wc -l)
    echo "Found \$mpnn_count FASTA files"
    echo ""
fi

# ============================================================================
# Step 4: RF2
# ============================================================================
if [ \${START_STEP} -le 4 ]; then
    echo "==> Step 4: Running RF2 structure prediction"
    echo ""

    # Run RF2
    bash ${PROJECT_DIR}/run_antibody_design.sh --rf2-only --output-dir ${CONFIG_OUTPUT_DIR} </dev/null

    if [ \$? -ne 0 ]; then
        echo "✗ RF2 failed"
        exit 1
    fi

    rf2_count=\$(ls ${CONFIG_OUTPUT_DIR}/rf2/design_*.pdb 2>/dev/null | wc -l)
    echo "✓ RF2 completed: \$rf2_count structures"
    echo ""
else
    echo "==> Step 4: Skipping RF2 (using existing results)"
    rf2_count=\$(ls ${CONFIG_OUTPUT_DIR}/rf2/design_*.pdb 2>/dev/null | wc -l)
    echo "Found \$rf2_count structures"
    echo ""
fi

# ============================================================================
# Step 5: Analysis
# ============================================================================
if [ \${START_STEP} -le 5 ]; then
    echo "==> Step 5: Analyzing RF2 structures"
    echo ""

    # Run analysis
    ANALYZE_SCRIPT="${PROJECT_DIR}/analyze_rf2_pdb.py"
    RF2_OUTPUT="${CONFIG_OUTPUT_DIR}/rf2"
    ANALYSIS_CSV="${CONFIG_OUTPUT_DIR}/rf2_pdb_analysis.csv"

    echo "Running: \${RFANTIBODY_PYTHON_BIN} \${ANALYZE_SCRIPT}"
    echo "  --rf2-dir \${RF2_OUTPUT}"
    echo "  --output \${ANALYSIS_CSV}"
    echo ""

    \${RFANTIBODY_PYTHON_BIN} "\${ANALYZE_SCRIPT}" \\
        --rf2-dir "\${RF2_OUTPUT}" \\
        --output "\${ANALYSIS_CSV}"

    if [ \$? -eq 0 ]; then
        echo "✓ Analysis complete: \${ANALYSIS_CSV}"
    else
        echo "⚠ Analysis script completed with warnings"
    fi
else
    echo "==> Step 5: Skipping analysis"
fi

echo ""
echo "=========================================="
echo "Pipeline Completed: \$(date)"
echo "=========================================="
echo "Summary:"
if [ \${START_STEP} -le 1 ]; then
    echo "  Initial designs: \$design_count"
fi
if [ \${START_STEP} -le 2 ]; then
    echo "  After filtering: \$filtered_count"
fi
if [ \${START_STEP} -le 3 ]; then
    echo "  ProteinMPNN FASTA: \$mpnn_count"
fi
if [ \${START_STEP} -le 4 ]; then
    echo "  RF2 structures: \$rf2_count"
fi
echo ""
echo "Results location: ${CONFIG_OUTPUT_DIR}/"
if [ -d "${CONFIG_OUTPUT_DIR}/rfdiffusion_raw" ]; then
    raw_count=\$(ls ${CONFIG_OUTPUT_DIR}/rfdiffusion_raw/design_*.pdb 2>/dev/null | wc -l)
    echo "  - rfdiffusion_raw/ : All \${raw_count} initial designs"
fi
if [ -d "${CONFIG_OUTPUT_DIR}/rfdiffusion" ]; then
    echo "  - rfdiffusion/     : Filtered designs"
fi
if [ -d "${CONFIG_OUTPUT_DIR}/proteinmpnn" ]; then
    echo "  - proteinmpnn/     : FASTA sequences"
fi
if [ -d "${CONFIG_OUTPUT_DIR}/rf2" ]; then
    echo "  - rf2/             : Final structures"
fi
if [ -f "${CONFIG_OUTPUT_DIR}/rf2_pdb_analysis.csv" ]; then
    echo "  - rf2_pdb_analysis.csv : Quality metrics"
fi
echo "=========================================="
EOF
    
    # Replace placeholders with actual values (escape special characters for sed)
    # Use @ as delimiter to avoid conflicts with paths containing /
    sed -i "s@NAME_PLACEHOLDER@${NAME}@g" "$JOB_SCRIPT"
    sed -i "s@HOTSPOTS_PLACEHOLDER@${HOTSPOTS}@g" "$JOB_SCRIPT"
    sed -i "s@CDR_LOOPS_PLACEHOLDER@${CDR_LOOPS}@g" "$JOB_SCRIPT"
    sed -i "s@NUM_DESIGNS_PLACEHOLDER@${NUM_DESIGNS}@g" "$JOB_SCRIPT"
    sed -i "s@CONFIG_OUTPUT_DIR_PLACEHOLDER@${CONFIG_OUTPUT_DIR}@g" "$JOB_SCRIPT"
    sed -i "s@INPUT_PDB_PLACEHOLDER@${INPUT_PDB}@g" "$JOB_SCRIPT"
    sed -i "s@FRAMEWORK_PDB_PLACEHOLDER@${FRAMEWORK_PDB}@g" "$JOB_SCRIPT"
    sed -i "s@BATCH_LOG_DIR_PLACEHOLDER@${BATCH_LOG_DIR}@g" "$JOB_SCRIPT"
    sed -i "s@PROJECT_DIR_PLACEHOLDER@${PROJECT_DIR}@g" "$JOB_SCRIPT"
    sed -i "s@WEIGHTS_DIR_PLACEHOLDER@${WEIGHTS_DIR}@g" "$JOB_SCRIPT"
    sed -i "s@TEMP_CONFIG_PLACEHOLDER@${TEMP_CONFIG}@g" "$JOB_SCRIPT"
    
    # Submit the job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")
    JOB_IDS+=("$JOB_ID")
    JOB_NAMES+=("$NAME")
    
    echo "  ✓ Submitted SLURM job: $JOB_ID"
    echo ""
    
    # Restore original run_config.txt immediately after job submission
    cp "$ORIGINAL_CONFIG_BACKUP" "$CONFIG_FILE"
done

# ============================================================================
# Summary
# ============================================================================
echo "=========================================="
echo "Batch Submission Complete"
echo "=========================================="
echo ""
echo "Submitted jobs:"
for i in "${!JOB_NAMES[@]}"; do
    echo "  [${JOB_IDS[$i]}] ${JOB_NAMES[$i]}"
done
echo ""
echo "Log directory: $BATCH_LOG_DIR"
echo ""
echo "To monitor jobs:"
echo "  squeue -u \$USER"
echo ""
echo "To check specific job:"
for i in "${!JOB_NAMES[@]}"; do
    echo "  squeue -j ${JOB_IDS[$i]}  # ${JOB_NAMES[$i]}"
done
echo ""
echo "To view job output:"
for i in "${!JOB_NAMES[@]}"; do
    echo "  tail -f ${BATCH_LOG_DIR}/${JOB_NAMES[$i]}_*.out"
done
echo ""
echo "To cancel jobs:"
echo "  scancel ${JOB_IDS[@]}"
echo ""
echo "Results will be in:"
for config in "${HOTSPOT_CONFIGS[@]}"; do
    IFS='|' read -r NAME VARIABLE_VALUE CONFIG_NUM_DESIGNS CONFIG_T_VALUE <<< "$config"
    echo "  ${BASE_OUTPUT_DIR}/${NAME}/"
done
echo ""
echo "=========================================="

# Clean up backup file
rm -f "$ORIGINAL_CONFIG_BACKUP"
