#!/bin/bash
# Batch submit AF3 jobs to HPC with automatic queue monitoring
# Handles cases where there are more than 66 jobs (HPC quick partition limit)

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Maximum number of jobs to submit in one batch (HPC quick partition limit)
MAX_JOBS_PER_BATCH=66

# Directory containing job folders (each with a JSON file)
JOBS_DIR="${1:-jobs}"

# Path to alphafold3.sh script
AF3_SCRIPT="${2:-./alphafold3.sh}"

# Wait time between batches (seconds) - time for jobs to start running
BATCH_WAIT_TIME=60

# Log directory
LOG_DIR="./logs"

# ============================================================================
# Validate inputs
# ============================================================================

if [ ! -d "$JOBS_DIR" ]; then
    echo "❌ Error: Jobs directory not found: $JOBS_DIR" >&2
    exit 1
fi

if [ ! -f "$AF3_SCRIPT" ]; then
    echo "❌ Error: AlphaFold3 script not found: $AF3_SCRIPT" >&2
    exit 1
fi

mkdir -p "$LOG_DIR"

# ============================================================================
# Find all JSON files in job directories
# ============================================================================

echo "=============================================="
echo "AF3 Batch Job Submission"
echo "=============================================="
echo "Jobs directory: $JOBS_DIR"
echo "AF3 script: $AF3_SCRIPT"
echo "Max jobs per batch: $MAX_JOBS_PER_BATCH"
echo "Wait between batches: ${BATCH_WAIT_TIME}s"
echo ""

# Find job directories and select one JSON file per directory
# Prefer pHLA_*.json over example.json to avoid duplicates
echo "Scanning for job directories..."
declare -a JSON_FILES
duplicate_count=0

for job_dir in "$JOBS_DIR"/*/; do
    if [ ! -d "$job_dir" ]; then
        continue
    fi
    
    job_name=$(basename "$job_dir")
    
    # Find all JSON files in this directory
    mapfile -t json_list < <(find "$job_dir" -maxdepth 1 -name "*.json")
    json_count=${#json_list[@]}
    
    if [ "$json_count" -eq 0 ]; then
        echo "  ⚠️  Warning: No JSON file found in $job_name, skipping"
        continue
    elif [ "$json_count" -gt 1 ]; then
        duplicate_count=$((duplicate_count + 1))
        echo "  ⚠️  Warning: Multiple JSON files in $job_name:"
        for jf in "${json_list[@]}"; do
            echo "      - $(basename "$jf")"
        done
        
        # Prefer pHLA_*.json over example.json
        json_file=$(printf '%s\n' "${json_list[@]}" | grep "pHLA_.*\.json" | head -1)
        if [ -z "$json_file" ]; then
            json_file="${json_list[0]}"
        fi
        echo "      Using: $(basename "$json_file")"
    else
        json_file="${json_list[0]}"
    fi
    
    # Add to array
    JSON_FILES+=("$json_file")
done

# Sort the array
IFS=$'\n' JSON_FILES=($(sort <<<"${JSON_FILES[*]}"))
unset IFS

if [ "$duplicate_count" -gt 0 ]; then
    echo ""
    echo "⚠️  Found $duplicate_count directories with multiple JSON files"
    echo "   Consider cleaning up duplicate files before submission"
    echo ""
fi

TOTAL_JOBS=${#JSON_FILES[@]}

if [ "$TOTAL_JOBS" -eq 0 ]; then
    echo "❌ No JSON files found in $JOBS_DIR subdirectories"
    exit 1
fi

echo "Found $TOTAL_JOBS job(s) to submit"
echo ""

# ============================================================================
# Calculate batching strategy
# ============================================================================

NUM_BATCHES=$(( (TOTAL_JOBS + MAX_JOBS_PER_BATCH - 1) / MAX_JOBS_PER_BATCH ))

if [ "$TOTAL_JOBS" -le "$MAX_JOBS_PER_BATCH" ]; then
    echo "✅ Total jobs ($TOTAL_JOBS) within single batch limit"
    echo "   Submitting all jobs in one batch..."
else
    echo "⚠️  Total jobs ($TOTAL_JOBS) exceeds single batch limit ($MAX_JOBS_PER_BATCH)"
    echo "   Will submit in $NUM_BATCHES batches"
fi

echo ""
echo "=============================================="
echo "Batch Submission Plan:"
echo "=============================================="

for (( batch=0; batch<NUM_BATCHES; batch++ )); do
    START_IDX=$((batch * MAX_JOBS_PER_BATCH))
    END_IDX=$(( START_IDX + MAX_JOBS_PER_BATCH ))
    if [ "$END_IDX" -gt "$TOTAL_JOBS" ]; then
        END_IDX=$TOTAL_JOBS
    fi
    BATCH_SIZE=$((END_IDX - START_IDX))
    echo "Batch $((batch + 1))/$NUM_BATCHES: Jobs $((START_IDX + 1))-$END_IDX ($BATCH_SIZE jobs)"
done

echo ""

echo "Proceed with submission? (y/n) "
read -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted by user"
    exit 0
fi

# ============================================================================
# Submit jobs in batches
# ============================================================================

submitted_count=0
failed_count=0

for (( batch=0; batch<NUM_BATCHES; batch++ )); do
    START_IDX=$((batch * MAX_JOBS_PER_BATCH))
    END_IDX=$(( START_IDX + MAX_JOBS_PER_BATCH ))
    if [ "$END_IDX" -gt "$TOTAL_JOBS" ]; then
        END_IDX=$TOTAL_JOBS
    fi
    
    echo ""
    echo "=============================================="
    echo "Submitting Batch $((batch + 1))/$NUM_BATCHES"
    echo "=============================================="
    echo ""
    
    for (( i=START_IDX; i<END_IDX; i++ )); do
        JSON_FILE="${JSON_FILES[$i]}"
        JOB_DIR="$(dirname "$JSON_FILE")"
        JOB_NAME="$(basename "$JOB_DIR")"
        
        echo "[$((i + 1))/$TOTAL_JOBS] Submitting job: $JOB_NAME"
        echo "  JSON: $JSON_FILE"
        
        # Submit with sbatch
        OUTPUT_LOG="${LOG_DIR}/${JOB_NAME}.out"
        ERROR_LOG="${LOG_DIR}/${JOB_NAME}.err"
        
        SUBMIT_CMD="sbatch --output=\"$OUTPUT_LOG\" --error=\"$ERROR_LOG\" \"$AF3_SCRIPT\" \"$JSON_FILE\""
        
        if eval "$SUBMIT_CMD"; then
            submitted_count=$((submitted_count + 1))
            echo "  ✅ Submitted"
        else
            failed_count=$((failed_count + 1))
            echo "  ❌ Failed to submit"
        fi
        echo ""
    done
    
    # If this is not the last batch, wait between submissions
    if [ $((batch + 1)) -lt "$NUM_BATCHES" ]; then
        echo ""
        echo "=============================================="
        echo "Batch $((batch + 1))/$NUM_BATCHES Submitted"
        echo "=============================================="
        echo ""
        echo "Waiting ${BATCH_WAIT_TIME}s before submitting next batch..."
        echo "Press Ctrl+C to stop (jobs already submitted will continue)"
        sleep "$BATCH_WAIT_TIME"
        echo "Continuing to next batch..."
        echo ""
    fi
done

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "Submission Complete"
echo "=============================================="
echo "Total jobs: $TOTAL_JOBS"
echo "Successfully submitted: $submitted_count"
echo "Failed: $failed_count"
echo ""
echo "Monitor job status with:"
echo "  squeue -u \$USER"
echo ""
echo "Check logs in:"
echo "  $LOG_DIR/"
echo ""
