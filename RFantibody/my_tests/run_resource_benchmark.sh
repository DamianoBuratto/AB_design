#!/bin/bash

# ============================================================================
# Resource Benchmark Script for RFantibody
# ============================================================================
# This script tests different GPU/CPU configurations with identical parameters
# to benchmark resource requirements and execution time
#
# Fixed Parameters:
#   - Hotspots: [T82,T85,T155,T156,T159]
#   - Designs: 50
#   - ProteinMPNN: 16 sequences per design
#   - T: 200
#   - CDR: [L1:10-13,L2:7,L3:9-11,H1:7,H2:6,H3:11-14]
#   - Deterministic: True
#   - Seed: 12818 (fixed for reproducibility)
# Usage: bash run_resource_benchmark.sh
# ============================================================================

set -e

# ============================================================================
# Fixed Parameters
# ============================================================================
FIXED_HOTSPOTS="[T82,T85,T155,T156,T159]"
FIXED_CDR_LOOPS="[L1:10-13,L2:7,L3:9-11,H1:7,H2:6,H3:11-14]"
FIXED_T=200
FIXED_NUM_DESIGNS=25
FIXED_MPNN_SEQS=2
FIXED_SEED=12818
DETERMINISTIC=True

# ============================================================================
# Resource Configurations to Test
# ============================================================================
# Format: "NAME|GPU_COUNT|CPU_COUNT"
# Memory is not specified - SLURM will use default allocation
# CPU configuration: test different GPU:CPU ratios (1:1 to 1:16)
# Per HPC guide, recommended: 16 CPUs per GPU (128/8=16)
declare -a RESOURCE_CONFIGS=(
    # Single GPU configurations (test 1:1, 1:2, 1:4, 1:8, 1:16)
    "1GPU_1CPU|1|1"      # 1:1 - minimal CPU, test CPU bottleneck
    "1GPU_2CPU|1|2"      # 1:2
    "1GPU_4CPU|1|4"      # 1:4
    "1GPU_8CPU|1|8"      # 1:8
    "1GPU_16CPU|1|16"    # 1:16 - recommended ratio
    
    # Dual GPU configurations (test 1:1, 1:2, 1:4, 1:8, 1:16)
    "2GPU_2CPU|2|2"      # 1:1
    "2GPU_4CPU|2|4"      # 1:2
    "2GPU_8CPU|2|8"      # 1:4
    "2GPU_16CPU|2|16"    # 1:8
    "2GPU_32CPU|2|32"    # 1:16 - recommended ratio
    
    # Quad GPU configurations (test 1:1, 1:2, 1:4, 1:8, 1:16)
    "4GPU_4CPU|4|4"      # 1:1
    "4GPU_8CPU|4|8"      # 1:2
    "4GPU_16CPU|4|16"    # 1:4
    "4GPU_32CPU|4|32"    # 1:8
    "4GPU_64CPU|4|64"    # 1:16 - recommended ratio
    
    # Note: 8GPU configurations removed - nodes 21-23,25,27-28 do not support 8 GPUs
)

# ============================================================================
# Load configuration from run_config.txt
# ============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/run_config.txt"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

echo "Loading configuration from: $CONFIG_FILE"
source "$CONFIG_FILE"

# Validate PROJECT_DIR
if [ -z "$PROJECT_DIR" ] || [ ! -d "$PROJECT_DIR" ]; then
    echo "ERROR: PROJECT_DIR not set or does not exist: $PROJECT_DIR"
    exit 1
fi

echo "Project directory: $PROJECT_DIR"

# Validate input files
if [ ! -f "$INPUT_PDB" ] || [ ! -f "$FRAMEWORK_PDB" ]; then
    echo "ERROR: Input files not found"
    exit 1
fi

# Create batch ID and logs directory
BATCH_ID=$(date +%Y%m%d_%H%M%S)
BENCHMARK_LOG_DIR="${PROJECT_DIR}/logs/resource_benchmark_${BATCH_ID}"
mkdir -p "$BENCHMARK_LOG_DIR"

# Create results summary file
SUMMARY_FILE="${BENCHMARK_LOG_DIR}/benchmark_summary.txt"
cat > "$SUMMARY_FILE" <<SUMMARY_HEADER
========================================
RFantibody Resource Benchmark
========================================
Batch ID: ${BATCH_ID}
Started: $(date)

Fixed Parameters:
  Hotspots:     ${FIXED_HOTSPOTS}
  CDR Loops:    ${FIXED_CDR_LOOPS}
  T:            ${FIXED_T}
  Designs:      ${FIXED_NUM_DESIGNS}
  MPNN Seqs:    ${FIXED_MPNN_SEQS}
  Seed:         ${FIXED_SEED}
  Deterministic: ${DETERMINISTIC}

Resource Configurations: ${#RESOURCE_CONFIGS[@]}

========================================
Results (will be updated as jobs complete):
========================================

SUMMARY_HEADER

echo "=========================================="
echo "RFantibody Resource Benchmark"
echo "=========================================="
echo "Batch ID: ${BATCH_ID}"
echo "Benchmark log directory: $BENCHMARK_LOG_DIR"
echo ""
echo "Fixed Parameters:"
echo "  Hotspots:     ${FIXED_HOTSPOTS}"
echo "  CDR Loops:    ${FIXED_CDR_LOOPS}"
echo "  T:            ${FIXED_T}"
echo "  Designs:      ${FIXED_NUM_DESIGNS}"
echo "  MPNN Seqs:    ${FIXED_MPNN_SEQS}"
echo "  Seed:         ${FIXED_SEED} (Deterministic=${DETERMINISTIC})"
echo ""
echo "Testing ${#RESOURCE_CONFIGS[@]} resource configurations"
echo ""

# ============================================================================
# Job tracking
# ============================================================================
declare -a JOB_IDS=()
declare -a JOB_NAMES=()

# ============================================================================
# Submit benchmark jobs
# ============================================================================
for config in "${RESOURCE_CONFIGS[@]}"; do
    # Parse configuration
    IFS='|' read -r NAME GPU_COUNT CPU_COUNT <<< "$config"
    
    echo "----------------------------------------"
    echo "Configuration: $NAME"
    echo "  GPU:    $GPU_COUNT"
    echo "  CPU:    $CPU_COUNT"
    echo "  Memory: (SLURM default)"
    
    # Create output directory with batch ID to avoid conflicts
    BASE_OUTPUT_DIR="${PROJECT_DIR}/outputs/benchmark_${BATCH_ID}"
    CONFIG_OUTPUT_DIR="${BASE_OUTPUT_DIR}/${NAME}"
    echo "  Output: $CONFIG_OUTPUT_DIR"
    
    mkdir -p "${CONFIG_OUTPUT_DIR}"
    
    # Create temporary config file
    TEMP_CONFIG="${BENCHMARK_LOG_DIR}/${NAME}_run_config.txt"
    cat > "$TEMP_CONFIG" <<CFGEOF
INPUT_PDB=${INPUT_PDB}
FRAMEWORK_PDB=${FRAMEWORK_PDB}
OUTPUT_DIR=${CONFIG_OUTPUT_DIR}
WEIGHTS_DIR=${WEIGHTS_DIR}
PROJECT_DIR=${PROJECT_DIR}
CFGEOF
    
    # Create SLURM job script
    JOB_SCRIPT="${BENCHMARK_LOG_DIR}/${NAME}_benchmark.slurm"
    
    cat > "$JOB_SCRIPT" <<'EOF'
#!/bin/bash
#SBATCH --job-name=Bench_BENCH_NAME
#SBATCH --output=BENCH_LOG_DIR/BENCH_NAME_%j.out
#SBATCH --error=BENCH_LOG_DIR/BENCH_NAME_%j.err
#SBATCH --partition=quick
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node20,node21,node24,node26,node27,node28
#SBATCH --gres=gpu:GPU_COUNT_VALUE
#SBATCH --cpus-per-task=CPU_COUNT_VALUE
#SBATCH --time=12:00:00

# Record start time
START_TIME=$(date +%s)
START_TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

# Record node information
RUNNING_NODE=$(hostname)
SLURM_NODE_INFO="${SLURM_NODELIST:-$RUNNING_NODE}"

echo "=========================================="
echo "Resource Benchmark: BENCH_NAME"
echo "=========================================="
echo "Job ID:   ${SLURM_JOB_ID}"
echo "Node:     ${SLURM_NODE_INFO}"
echo "Hostname: ${RUNNING_NODE}"
echo "Started:  ${START_TIMESTAMP}"
echo ""
echo "Resource Allocation:"
echo "  GPUs:   GPU_COUNT_VALUE"
echo "  CPUs:   CPU_COUNT_VALUE"
echo "  Memory: (SLURM default)"
echo ""
echo "Fixed Parameters:"
echo "  Hotspots:     HOTSPOTS_VALUE"
echo "  CDR Loops:    CDR_LOOPS_VALUE"
echo "  T:            T_VALUE"
echo "  Designs:      NUM_DESIGNS_VALUE"
echo "  MPNN Seqs:    MPNN_SEQS_VALUE"
echo "  Seed:         SEED_VALUE (Deterministic=DETERMINISTIC_VALUE)"
echo ""
echo "Output Directory: CONFIG_OUTPUT_DIR_VALUE"
echo "=========================================="
echo ""

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

# Setup environment variables
export PYTHONPATH=PROJECT_DIR_VALUE/src:PROJECT_DIR_VALUE/src/rfantibody/rfdiffusion:PROJECT_DIR_VALUE/src/rfantibody:PROJECT_DIR_VALUE/include/SE3Transformer:$PYTHONPATH
export PATH=PROJECT_DIR_VALUE/include/USalign:$PATH
export WEIGHTS_DIR=WEIGHTS_DIR_VALUE
export ICECREAM_COLORS=never
export HYDRA_FULL_ERROR=1
export PYTHONHASHSEED=SEED_VALUE
export PYTHONUNBUFFERED=1
export NO_COLOR=1
export TERM=dumb
export CONFIG_FILE="TEMP_CONFIG_VALUE"

cd PROJECT_DIR_VALUE || exit 1

# ============================================================================
# Step 1: RFdiffusion
# ============================================================================
STEP_START=$(date +%s)
echo "==> Step 1: RFdiffusion (NUM_DESIGNS_VALUE designs, deterministic)"
echo ""

mkdir -p scripts/config
ln -sfn "PROJECT_DIR_VALUE/src/rfantibody/rfdiffusion/config/inference" scripts/config/inference
mkdir -p "CONFIG_OUTPUT_DIR_VALUE/rfdiffusion"

run_with_glibc "${RFANTIBODY_PYTHON_BIN}" scripts/rfdiffusion_inference.py \
    --config-name antibody \
    antibody.target_pdb="INPUT_PDB_VALUE" \
    antibody.framework_pdb="FRAMEWORK_PDB_VALUE" \
    inference.ckpt_override_path="WEIGHTS_DIR_VALUE/RFdiffusion_Ab.pt" \
    'ppi.hotspot_res=HOTSPOTS_VALUE' \
    'antibody.design_loops=CDR_LOOPS_VALUE' \
    inference.num_designs=NUM_DESIGNS_VALUE \
    diffuser.T=T_VALUE \
    ++inference.seed=SEED_VALUE \
    inference.output_prefix="CONFIG_OUTPUT_DIR_VALUE/rfdiffusion/design" \
    2>&1 | grep -v "ic|" || true

RFDIFF_EXIT=$?
STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

if [ $RFDIFF_EXIT -ne 0 ]; then
    echo "✗ RFdiffusion failed with exit code $RFDIFF_EXIT"
    exit $RFDIFF_EXIT
fi

design_count=$(ls CONFIG_OUTPUT_DIR_VALUE/rfdiffusion/design_*.pdb 2>/dev/null | wc -l)
echo "✓ RFdiffusion completed: $design_count designs (${STEP_DURATION}s)"
echo ""

# ============================================================================
# Step 2: Filtering (DISABLED for benchmark - using all designs)
# ============================================================================
STEP_START=$(date +%s)
echo "==> Step 2: Filtering DISABLED - using all designs for benchmark"
echo ""

# Skip filtering - use all RFdiffusion designs directly
filtered_count=$(ls CONFIG_OUTPUT_DIR_VALUE/rfdiffusion/design_*.pdb 2>/dev/null | wc -l)

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

echo "✓ Using all ${filtered_count} designs (filtering disabled) (${STEP_DURATION}s)"
echo ""

# ============================================================================
# Step 3: ProteinMPNN
# ============================================================================
STEP_START=$(date +%s)
echo "==> Step 3: ProteinMPNN (MPNN_SEQS_VALUE sequences per design)"
echo ""

bash PROJECT_DIR_VALUE/run_antibody_design.sh --proteinmpnn-only --output-dir CONFIG_OUTPUT_DIR_VALUE </dev/null

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

mpnn_count=$(ls CONFIG_OUTPUT_DIR_VALUE/proteinmpnn/design_*.fa 2>/dev/null | wc -l)
echo "✓ ProteinMPNN completed: $mpnn_count FASTA files (${STEP_DURATION}s)"
echo ""

# ============================================================================
# Step 4: RF2
# ============================================================================
# Record time before RF2 for phase timing analysis
PRE_RF2_END=$(date +%s)
PRE_RF2_DURATION=$((PRE_RF2_END - START_TIME))

STEP_START=$(date +%s)
echo "==> Step 4: RF2 structure prediction"
echo ""

bash PROJECT_DIR_VALUE/run_antibody_design.sh --rf2-only --output-dir CONFIG_OUTPUT_DIR_VALUE </dev/null

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

rf2_count=$(ls CONFIG_OUTPUT_DIR_VALUE/rf2/design_*.pdb 2>/dev/null | wc -l)
echo "✓ RF2 completed: $rf2_count structures (${STEP_DURATION}s)"
echo ""

# ============================================================================
# Step 5: Analysis
# ============================================================================
STEP_START=$(date +%s)
echo "==> Step 5: Analysis"
echo ""

run_with_glibc "${RFANTIBODY_PYTHON_BIN}" analyze_rf2_pdb.py \
    --rf2-dir "CONFIG_OUTPUT_DIR_VALUE/rf2" \
    --output "CONFIG_OUTPUT_DIR_VALUE/rf2_pdb_analysis.csv"

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

echo "✓ Analysis completed (${STEP_DURATION}s)"
echo ""

# ============================================================================
# Final Summary
# ============================================================================
END_TIME=$(date +%s)
END_TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
TOTAL_DURATION=$((END_TIME - START_TIME))
RF2_DURATION=$((END_TIME - PRE_RF2_END))

HOURS=$((TOTAL_DURATION / 3600))
MINUTES=$(((TOTAL_DURATION % 3600) / 60))
SECONDS=$((TOTAL_DURATION % 60))

PRE_RF2_HOURS=$((PRE_RF2_DURATION / 3600))
PRE_RF2_MINUTES=$(((PRE_RF2_DURATION % 3600) / 60))
PRE_RF2_SECONDS=$((PRE_RF2_DURATION % 60))

RF2_HOURS=$((RF2_DURATION / 3600))
RF2_MINUTES=$(((RF2_DURATION % 3600) / 60))
RF2_SECONDS=$((RF2_DURATION % 60))

echo "=========================================="
echo "Benchmark Completed: BENCH_NAME"
echo "=========================================="
echo "Started:  ${START_TIMESTAMP}"
echo "Finished: ${END_TIMESTAMP}"
echo ""
echo "Total Time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${TOTAL_DURATION} seconds)"
echo "  Pre-RF2 (RFdiff+Filter+MPNN): ${PRE_RF2_HOURS}h ${PRE_RF2_MINUTES}m ${PRE_RF2_SECONDS}s (${PRE_RF2_DURATION}s)"
echo "  RF2 Prediction:               ${RF2_HOURS}h ${RF2_MINUTES}m ${RF2_SECONDS}s (${RF2_DURATION}s)"
echo ""
echo "Pipeline Summary:"
echo "  Initial designs:   $design_count"
echo "  After filtering:   $filtered_count"
echo "  ProteinMPNN FASTA: $mpnn_count"
echo "  RF2 structures:    $rf2_count"
echo ""
echo "Resource Configuration:"
echo "  GPUs:   GPU_COUNT_VALUE"
echo "  CPUs:   CPU_COUNT_VALUE"
echo "  Memory: (SLURM default)"
echo "=========================================="

# Write results to summary file
{
    echo "----------------------------------------"
    echo "Configuration: BENCH_NAME"
    echo "  Job ID:    ${SLURM_JOB_ID}"
    echo "  Node:      ${SLURM_NODE_INFO}"
    echo "  Resources: GPU_COUNT_VALUE GPU, CPU_COUNT_VALUE CPU"
    echo "  Total Time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${TOTAL_DURATION}s)"
    echo "    ├─ Pre-RF2:  ${PRE_RF2_HOURS}h ${PRE_RF2_MINUTES}m ${PRE_RF2_SECONDS}s (${PRE_RF2_DURATION}s)"
    echo "    └─ RF2:      ${RF2_HOURS}h ${RF2_MINUTES}m ${RF2_SECONDS}s (${RF2_DURATION}s)"
    echo "  Designs: $design_count -> $filtered_count -> $rf2_count"
    echo "  Started:  ${START_TIMESTAMP}"
    echo "  Finished: ${END_TIMESTAMP}"
    echo ""
} >> SUMMARY_FILE_VALUE
EOF

    # Replace placeholders
    sed -i "s|BENCH_NAME|${NAME}|g" "$JOB_SCRIPT"
    sed -i "s|BENCH_LOG_DIR|${BENCHMARK_LOG_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|GPU_COUNT_VALUE|${GPU_COUNT}|g" "$JOB_SCRIPT"
    sed -i "s|CPU_COUNT_VALUE|${CPU_COUNT}|g" "$JOB_SCRIPT"
    sed -i "s|PROJECT_DIR_VALUE|${PROJECT_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|CONFIG_OUTPUT_DIR_VALUE|${CONFIG_OUTPUT_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|TEMP_CONFIG_VALUE|${TEMP_CONFIG}|g" "$JOB_SCRIPT"
    sed -i "s|WEIGHTS_DIR_VALUE|${WEIGHTS_DIR}|g" "$JOB_SCRIPT"
    sed -i "s|INPUT_PDB_VALUE|${INPUT_PDB}|g" "$JOB_SCRIPT"
    sed -i "s|FRAMEWORK_PDB_VALUE|${FRAMEWORK_PDB}|g" "$JOB_SCRIPT"
    sed -i "s|HOTSPOTS_VALUE|${FIXED_HOTSPOTS}|g" "$JOB_SCRIPT"
    sed -i "s|CDR_LOOPS_VALUE|${FIXED_CDR_LOOPS}|g" "$JOB_SCRIPT"
    sed -i "s|NUM_DESIGNS_VALUE|${FIXED_NUM_DESIGNS}|g" "$JOB_SCRIPT"
    sed -i "s|T_VALUE|${FIXED_T}|g" "$JOB_SCRIPT"
    sed -i "s|MPNN_SEQS_VALUE|${FIXED_MPNN_SEQS}|g" "$JOB_SCRIPT"
    sed -i "s|SEED_VALUE|${FIXED_SEED}|g" "$JOB_SCRIPT"
    sed -i "s|DETERMINISTIC_VALUE|${DETERMINISTIC}|g" "$JOB_SCRIPT"
    sed -i "s|SUMMARY_FILE_VALUE|${SUMMARY_FILE}|g" "$JOB_SCRIPT"
    
    # Submit job
    JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")
    JOB_IDS+=("$JOB_ID")
    JOB_NAMES+=("$NAME")
    
    echo "  ✓ Submitted: Job ID ${JOB_ID}"
    echo ""
done

# ============================================================================
# Summary
# ============================================================================
echo "=========================================="
echo "Benchmark Submission Complete"
echo "=========================================="
echo ""
echo "Submitted ${#JOB_IDS[@]} benchmark jobs:"
for i in "${!JOB_NAMES[@]}"; do
    echo "  [${JOB_IDS[$i]}] ${JOB_NAMES[$i]}"
done
echo ""
echo "Log directory: $BENCHMARK_LOG_DIR"
echo "Summary file:  $SUMMARY_FILE"
echo ""
echo "To monitor jobs:"
echo "  squeue -u \$USER"
echo ""
echo "To view summary:"
echo "  cat $SUMMARY_FILE"
echo ""
echo "To cancel all benchmark jobs:"
echo "  scancel ${JOB_IDS[@]}"
echo ""
echo "Results will be in: ${PROJECT_DIR}/outputs/benchmark/"
echo "=========================================="
