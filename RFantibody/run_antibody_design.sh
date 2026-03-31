#!/bin/bash

# Usage:
#   poetry run bash run_antibody_design.sh                    # Run all phases (default)
#
# Options:
#   --clean                : Automatically clean previous results without prompting
#   --no-clean             : Skip cleanup without prompting
#   --rfdiffusion-only     : Only run RFdiffusion (Phase 1)
#   --proteinmpnn-only     : Only run ProteinMPNN (Phase 2, requires RFdiffusion results)
#   --rf2-only             : Only run RF2 (Phase 3, requires ProteinMPNN results)
#   --phases X,Y,Z         : Run specific phases (e.g., --phases 1,2 or --phases 2,3)  
#
# Before use:
# Have run prepare_inputs.sh to complete environment check and input preparation

set -e  # Exit script immediately on any error

# Remember how many args were passed originally so we can detect interactive runs
ORIGINAL_ARGC="$#"

# ============================================================================
# Logging Setup
# ============================================================================
# Determine the project directory first
if [ "$(pwd)" = "/home" ] && [ -d "/home/weights" ]; then
    PROJECT_DIR="/home"
else
    PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Create logs directory if it doesn't exist
LOGS_DIR="${PROJECT_DIR}/logs"
mkdir -p "$LOGS_DIR"

# Generate timestamp for log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOGS_DIR}/rfantibody_design_${TIMESTAMP}.log"

# Function to log messages with timestamp
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Redirect all output to log file while still showing on screen
exec > >(tee -a "$LOG_FILE") 2>&1

log_message "=== RFantibody Design Pipeline Started ==="
log_message "Log file: $LOG_FILE"

# Parse command line arguments
CLEAN_REQUESTED=false
RUN_PHASES="1,2,3"  # Default: run all phases

while [ "$#" -gt 0 ]; do
    arg="$1"
    case $arg in
        --clean)
            log_message "Auto-cleaning previous results (--clean flag provided)..."
            CLEAN_REQUESTED=true
            shift
            ;;
        --no-clean)
            log_message "Skipping cleanup (--no-clean flag provided)..."
            CLEAN_REQUESTED=false
            shift
            ;;
        --rfdiffusion-only)
            log_message "RFdiffusion-only mode enabled..."
            RUN_PHASES="1"
            CLEAN_REQUESTED=false
            shift
            ;;
        --proteinmpnn-only)
            log_message "ProteinMPNN-only mode enabled..."
            RUN_PHASES="2"
            CLEAN_REQUESTED=false
            shift
            ;;
        --rf2-only)
            log_message "RF2-only mode enabled..."
            RUN_PHASES="3"
            CLEAN_REQUESTED=false
            shift
            ;;
        --phases)
            shift
            if [ "$#" -eq 0 ]; then
                log_message "ERROR: --phases requires an argument (e.g. --phases 1,2)"
                exit 1
            fi
            RUN_PHASES="$1"
            log_message "Custom phases specified: $RUN_PHASES"
            CLEAN_REQUESTED=false
            shift
            ;;
        --output-dir)
            shift
            if [ "$#" -eq 0 ]; then
                log_message "ERROR: --output-dir requires a directory argument"
                exit 1
            fi
            OUTPUT_DIR_OVERRIDE="$1"
            log_message "Overriding OUTPUT_DIR from cmdline: $OUTPUT_DIR_OVERRIDE"
            shift
            ;;
        --phases=*)
            RUN_PHASES="${arg#*=}"
            log_message "Custom phases specified: $RUN_PHASES"
            CLEAN_REQUESTED=false
            shift
            ;;
        *)
            log_message "Unknown argument: $arg"
            log_message "Use --help for usage information"
            exit 1
            ;;
    esac
done

if [ "$ORIGINAL_ARGC" -eq 0 ]; then
    if [ -t 0 ]; then
        read -p "Do you want to clean previous results? (y/N): " -r
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            log_message "Cleaning previous results..."
            CLEAN_REQUESTED=true
        fi
    else
        log_message "No CLI args provided and stdin is non-interactive; skipping cleanup prompt."
    fi
fi

# ============================================================================
# Global Variable Definitions
# ============================================================================
# PROJECT_DIR is already set above in logging setup

# Design parameters
# IMPORTANT: These hotspot residues correspond to the T-chain numbering in the converted PDB
HOTSPOTS="[T7,T8,T9]" 
CDR_LOOPS="[L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13]"
NUM_DESIGNS=100
NUM_SEQ_PER_TARGET=32

# Set Python path to include all necessary directories (based on example_test.slurm)
export PYTHONPATH="${PROJECT_DIR}/src:${PROJECT_DIR}/src/rfantibody/rfdiffusion:${PROJECT_DIR}/src/rfantibody:${PROJECT_DIR}/include/SE3Transformer:${PYTHONPATH}"

# Disable ANSI color codes in logs
export PYTHONUNBUFFERED=1
export NO_COLOR=1
export TERM=dumb
export ICECREAM_COLORS=never

# ============================================================================
# Detect environment and set Python command
# ============================================================================
# Define run_with_glibc function if in HPC environment
if [ -n "${RFANTIBODY_PYTHON_BIN:-}" ] && [ -n "${GLIBC_LOADER:-}" ]; then
    # HPC environment - define the wrapper function
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
    
    log_message "HPC environment detected, using GLIBC wrapper"
    PYTHON_CMD="run_with_glibc ${RFANTIBODY_PYTHON_BIN}"
    PYTHON_SIMPLE="${RFANTIBODY_PYTHON_BIN}"
else
    log_message "Docker/local environment detected, using standard python3"
    PYTHON_CMD="python3"
    PYTHON_SIMPLE="python3"
fi

# Read parameters from config file
# Allow CONFIG_FILE to be set via environment variable (for parallel batch jobs)
# Otherwise use default location
if [ -z "${CONFIG_FILE}" ]; then
    CONFIG_FILE="${PROJECT_DIR}/run_config.txt"
fi
log_message "Looking for config file at: $CONFIG_FILE"
log_message "Current PROJECT_DIR: $PROJECT_DIR"
log_message "Current working directory: $(pwd)"

if [ ! -f "$CONFIG_FILE" ]; then
    log_message "ERROR: Config file does not exist: $CONFIG_FILE"
    log_message "Please run prepare_inputs.sh first"
    log_message "Available files in PROJECT_DIR:"
    ls -la "$PROJECT_DIR" | head -10
    exit 1
fi

# Load configuration
source "$CONFIG_FILE"
log_message "Configuration loaded from: $CONFIG_FILE"

# If output dir provided on command line, override the config value
if [ -n "${OUTPUT_DIR_OVERRIDE:-}" ]; then
    OUTPUT_DIR="${OUTPUT_DIR_OVERRIDE}"
    log_message "OUTPUT_DIR overridden to: ${OUTPUT_DIR}"
fi

# Clean previous results if requested and now OUTPUT_DIR is available
if [ "$CLEAN_REQUESTED" = true ] && [ -n "$OUTPUT_DIR" ]; then
    log_message "Now cleaning previous results from: $OUTPUT_DIR"
    rm -rf "${OUTPUT_DIR}/rfdiffusion" "${OUTPUT_DIR}/proteinmpnn" "${OUTPUT_DIR}/rf2"
    log_message "Previous results cleaned."
fi

# Ensure output directories exist
log_message "Creating output directories..."
mkdir -p "${OUTPUT_DIR}/rfdiffusion" "${OUTPUT_DIR}/proteinmpnn" "${OUTPUT_DIR}/rf2"
log_message "Output structure:"
log_message "  RFdiffusion results: ${OUTPUT_DIR}/rfdiffusion/"
log_message "  ProteinMPNN results: ${OUTPUT_DIR}/proteinmpnn/"
log_message "  RF2 results: ${OUTPUT_DIR}/rf2/"

# ============================================================================
# Function: check_prerequisites - Check if required inputs exist for each phase
# ============================================================================
check_prerequisites() {
    local phase="$1"
    
    case $phase in
        1)
            # Phase 1 (RFdiffusion) requires input PDB files
            if [ ! -f "$INPUT_PDB" ]; then
                log_message "ERROR: Input PDB file not found: $INPUT_PDB"
                return 1
            fi
            if [ ! -f "$FRAMEWORK_PDB" ]; then
                log_message "ERROR: Framework PDB file not found: $FRAMEWORK_PDB"
                return 1
            fi
            ;;
        2)
            # Phase 2 (ProteinMPNN) requires RFdiffusion results
            if [ ! -d "${OUTPUT_DIR}/rfdiffusion" ] || [ -z "$(ls -A ${OUTPUT_DIR}/rfdiffusion/*.pdb 2>/dev/null)" ]; then
                log_message "ERROR: No RFdiffusion results found in ${OUTPUT_DIR}/rfdiffusion/"
                log_message "Please run Phase 1 (RFdiffusion) first."
                return 1
            fi
            ;;
        3)
            # Phase 3 (RF2) requires ProteinMPNN results
            if [ ! -d "${OUTPUT_DIR}/proteinmpnn" ] || [ -z "$(ls -A ${OUTPUT_DIR}/proteinmpnn/*.pdb 2>/dev/null)" ]; then
                log_message "ERROR: No ProteinMPNN results found in ${OUTPUT_DIR}/proteinmpnn/"
                log_message "Please run Phase 2 (ProteinMPNN) first."
                return 1
            fi
            ;;
    esac
    return 0
}
# ============================================================================
# Function: run_rfdiffusion - Run RFdiffusion for antibody backbone design
# ============================================================================
run_rfdiffusion() {
    log_message "=========================================="
    log_message "Phase 1: RFdiffusion backbone design"
    log_message "=========================================="
    log_message "Design parameters:"
    log_message "  - Target PDB: $INPUT_PDB"
    log_message "  - Framework PDB: $FRAMEWORK_PDB"
    log_message "  - CDR loops: $CDR_LOOPS"
    log_message "  - Number of designs: $NUM_DESIGNS"
    
    # Use the successful pattern from example_test.slurm:
    # 1. Stay in project root directory (don't cd to rfdiffusion)
    # 2. Set PYTHONPATH with rfdiffusion paths
    # 3. Call rfdiffusion_inference.py directly (not wrapper)
    
    # Set PYTHONPATH exactly as in example_test.slurm
    local pythonpath="${PROJECT_DIR}/src:${PROJECT_DIR}/src/rfantibody/rfdiffusion:${PROJECT_DIR}/src/rfantibody:${PROJECT_DIR}/include/SE3Transformer"
    if [ -n "${PYTHONPATH:-}" ]; then
        pythonpath="${pythonpath}:${PYTHONPATH}"
    fi
    
    # Build command array - use rfdiffusion_inference.py directly like example
    local -a cmd=(
        python3 "${PROJECT_DIR}/scripts/rfdiffusion_inference.py"
        --config-name antibody
        "antibody.target_pdb=${INPUT_PDB}"
        "antibody.framework_pdb=${FRAMEWORK_PDB}"
        "inference.ckpt_override_path=${WEIGHTS_DIR}/RFdiffusion_Ab.pt"
        "antibody.design_loops=${CDR_LOOPS}"
        "inference.num_designs=${NUM_DESIGNS}"
        "inference.output_prefix=${OUTPUT_DIR}/rfdiffusion/design"
        "ppi.hotspot_res=${HOTSPOTS}"
        "diffuser.T=200"
        "inference.cautious=False"
    )

    log_message "Executing RFdiffusion (using example_test.slurm pattern)..."
    log_message "Working directory: $(pwd)"

    # Create symlink for config directory (required by Hydra)
    mkdir -p "${PROJECT_DIR}/scripts/config"
    ln -sfn "${PROJECT_DIR}/src/rfantibody/rfdiffusion/config/inference" "${PROJECT_DIR}/scripts/config/inference"
    
    # Verify symlink was created
    if [ ! -e "${PROJECT_DIR}/scripts/config/inference" ]; then
        log_message "ERROR: Failed to create config symlink at ${PROJECT_DIR}/scripts/config/inference"
        log_message "Source: ${PROJECT_DIR}/src/rfantibody/rfdiffusion/config/inference"
        exit 1
    fi
    
    log_message "Config symlink verified: scripts/config/inference -> src/rfantibody/rfdiffusion/config/inference"

    # Execute with proper PYTHONPATH
    if ! PYTHONPATH="${pythonpath}" "${cmd[@]}"; then
        log_message "Error: RFdiffusion execution failed"
        exit 1
    fi
    
    # Verify results
    local design_count=$(ls ${OUTPUT_DIR}/rfdiffusion/design_*.pdb 2>/dev/null | wc -l)
    if [ "$design_count" -lt 1 ]; then
        log_message "Error: RFdiffusion generated no designs"
        exit 1
    fi
    
    log_message "RFdiffusion completed, generated $design_count designs"
    log_message "Phase 1 results location: ${OUTPUT_DIR}/rfdiffusion/"
}

# ============================================================================
# Function: run_proteinmpnn - Run ProteinMPNN for sequence optimization
# ============================================================================
run_proteinmpnn() {
    log_message "=========================================="
    log_message "Phase 2: ProteinMPNN sequence optimization"
    log_message "=========================================="
    
    # Export WEIGHTS_DIR for ProteinMPNN script to find checkpoint
    export WEIGHTS_DIR
    
    PYTHONPATH="${PROJECT_DIR}/src:$PYTHONPATH" ${PYTHON_SIMPLE} ${PROJECT_DIR}/scripts/proteinmpnn_interface_design.py \
        -pdbdir "${OUTPUT_DIR}/rfdiffusion" \
        -outpdbdir "${OUTPUT_DIR}/proteinmpnn" \
        -seqs_per_struct $NUM_SEQ_PER_TARGET \
        -temperature 0.1
    
    local seq_count=$(ls ${OUTPUT_DIR}/proteinmpnn/*.pdb 2>/dev/null | wc -l)
    log_message "ProteinMPNN completed, generated $seq_count sequence designs"
    log_message "Phase 2 results location: ${OUTPUT_DIR}/proteinmpnn/"
}

# ============================================================================
# Function: run_rf2 - Run RF2 for structure prediction
# ============================================================================
run_rf2() {
    log_message "=========================================="
    log_message "Phase 3: RF2 structure prediction"
    log_message "=========================================="
    
    PYTHONPATH="${PROJECT_DIR}/src:$PYTHONPATH" ${PYTHON_CMD} ${PROJECT_DIR}/scripts/rf2_predict.py \
        input.pdb_dir="${OUTPUT_DIR}/proteinmpnn" \
        output.pdb_dir="${OUTPUT_DIR}/rf2" \
        model.model_weights="${WEIGHTS_DIR}/RF2_ab.pt"
    
    local pred_count=$(ls ${OUTPUT_DIR}/rf2/*.pdb 2>/dev/null | wc -l)
    log_message "RF2 completed, predicted $pred_count structures"
    log_message "Phase 3 results location: ${OUTPUT_DIR}/rf2/"
    
    # Analyze RF2 results automatically
    log_message "Running RF2 results analysis..."
    if [ -f "${PROJECT_DIR}/analyze_rf2_pdb.py" ]; then
        PYTHONPATH="${PROJECT_DIR}/src:$PYTHONPATH" ${PYTHON_SIMPLE} "${PROJECT_DIR}/analyze_rf2_pdb.py" \
            --rf2-dir "${OUTPUT_DIR}/rf2" \
            --output "${OUTPUT_DIR}/rf2/rf2_analysis_results.csv"
        
        if [ -f "${OUTPUT_DIR}/rf2/rf2_analysis_results.csv" ]; then
            log_message "RF2 analysis completed. Results saved to: ${OUTPUT_DIR}/rf2/rf2_analysis_results.csv"
        else
            log_message "Warning: RF2 analysis failed to generate CSV output"
        fi
    else
        log_message "Warning: RF2 analysis script not found, skipping analysis"
        log_message "You can manually analyze results later using: python3 analyze_rf2_pdb.py"
    fi
}

# ============================================================================
# Main Function
# ============================================================================
main() {
    log_message "=========================================="
    log_message "RFantibody p53_R175H Antibody Design Pipeline"
    log_message "=========================================="
    log_message "Phases to run: $RUN_PHASES"
    
    # Parse phases to run
    IFS=',' read -ra PHASES <<< "$RUN_PHASES"
    
    # Validate and run each requested phase
    for phase in "${PHASES[@]}"; do
        # Trim whitespace
        phase=$(echo "$phase" | tr -d ' ')
        
        case $phase in
            1)
                log_message "Checking prerequisites for Phase 1 (RFdiffusion)..."
                if check_prerequisites 1; then
                    run_rfdiffusion
                else
                    log_message "Phase 1 prerequisites not met. Exiting."
                    exit 1
                fi
                ;;
            2)
                log_message "Checking prerequisites for Phase 2 (ProteinMPNN)..."
                if check_prerequisites 2; then
                    run_proteinmpnn
                else
                    log_message "Phase 2 prerequisites not met. Exiting."
                    exit 1
                fi
                ;;
            3)
                log_message "Checking prerequisites for Phase 3 (RF2)..."
                if check_prerequisites 3; then
                    run_rf2
                else
                    log_message "Phase 3 prerequisites not met. Exiting."
                    exit 1
                fi
                ;;
            *)
                log_message "ERROR: Invalid phase number: $phase"
                log_message "Valid phases are: 1 (RFdiffusion), 2 (ProteinMPNN), 3 (RF2)"
                exit 1
                ;;
        esac
    done
        
    log_message "=========================================="
    log_message "Pipeline completed successfully!"
    log_message "Executed phases: $RUN_PHASES"
    log_message "Results location: ${OUTPUT_DIR}/"
    log_message "=========================================="
}


# Run main function
main "$@"