#!/bin/bash
# Smart queue manager for AF3 jobs
# Monitors the queue and automatically submits new jobs to keep queue full
# Run with: screen -S af3_submit or tmux

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

JOBS_DIR="${1:-jobs}"
AF3_SCRIPT="${2:-./alphafold3.sh}"
LOG_DIR="./logs"

# Queue management settings
MAX_QUEUE_SIZE=50          # Maximum total jobs in queue (RUNNING + PENDING)
MIN_QUEUE_SIZE=20          # Minimum jobs before submitting more
SUBMIT_BATCH_SIZE=30       # How many jobs to submit each time
CHECK_INTERVAL=180         # Check queue every 3 minutes (180 seconds)
MUT_ONLY=true             # Set to true to skip wild-type (WT) jobs

# Tracking file
SUBMITTED_JOBS_FILE="./submitted_jobs.txt"
PROGRESS_LOG="./submission_progress.log"

mkdir -p "$LOG_DIR"

# ============================================================================
# Helper Functions
# ============================================================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$PROGRESS_LOG"
}

get_queue_count() {
    # Count AF3 jobs (RUNNING + PENDING) for current user
    squeue -u "$USER" -n "AF3" -h 2>/dev/null | wc -l || echo "0"
}

get_running_count() {
    squeue -u "$USER" -n "AF3" -t RUNNING -h 2>/dev/null | wc -l || echo "0"
}

get_pending_count() {
    squeue -u "$USER" -n "AF3" -t PENDING -h 2>/dev/null | wc -l || echo "0"
}

# ============================================================================
# Find all JSON files
# ============================================================================

log_message "=========================================="
log_message "AF3 Smart Queue Manager Started"
log_message "=========================================="
log_message "Jobs directory: $JOBS_DIR"
log_message "AF3 script: $AF3_SCRIPT"
log_message "Max queue size: $MAX_QUEUE_SIZE"
log_message "Min queue size: $MIN_QUEUE_SIZE"
log_message "Submit batch size: $SUBMIT_BATCH_SIZE"
log_message "Check interval: ${CHECK_INTERVAL}s"
log_message ""

# Validate inputs
if [ ! -d "$JOBS_DIR" ]; then
    log_message "❌ Error: Jobs directory not found: $JOBS_DIR"
    exit 1
fi

if [ ! -f "$AF3_SCRIPT" ]; then
    log_message "❌ Error: AlphaFold3 script not found: $AF3_SCRIPT"
    exit 1
fi

# Find all JSON files
log_message "Scanning for job directories..."
declare -a ALL_JSON_FILES

for job_dir in "$JOBS_DIR"/*/; do
    if [ ! -d "$job_dir" ]; then
        continue
    fi
    
    # Skip WT jobs if MUT_ONLY is enabled
    job_name=$(basename "$job_dir")
    if [ "$MUT_ONLY" = true ] && [[ "$job_name" == *_wt ]]; then
        continue
    fi
    
    # Find JSON files
    mapfile -t json_list < <(find "$job_dir" -maxdepth 1 -name "*.json")
    json_count=${#json_list[@]}
    
    if [ "$json_count" -eq 0 ]; then
        continue
    elif [ "$json_count" -gt 1 ]; then
        # Prefer pHLA_*.json
        json_file=$(printf '%s\n' "${json_list[@]}" | grep "pHLA_.*\.json" | head -1)
        if [ -z "$json_file" ]; then
            json_file="${json_list[0]}"
        fi
    else
        json_file="${json_list[0]}"
    fi
    
    ALL_JSON_FILES+=("$json_file")
done

# Sort files
IFS=$'\n' ALL_JSON_FILES=($(sort <<<"${ALL_JSON_FILES[*]}"))
unset IFS

TOTAL_JOBS=${#ALL_JSON_FILES[@]}
log_message "Found $TOTAL_JOBS total jobs to process"
log_message ""

if [ "$TOTAL_JOBS" -eq 0 ]; then
    log_message "❌ No jobs found"
    exit 1
fi

# ============================================================================
# Load already submitted jobs
# ============================================================================

declare -A SUBMITTED_JOBS

if [ -f "$SUBMITTED_JOBS_FILE" ]; then
    log_message "Loading previously submitted jobs from $SUBMITTED_JOBS_FILE"
    while IFS= read -r json_file; do
        SUBMITTED_JOBS["$json_file"]=1
    done < "$SUBMITTED_JOBS_FILE"
    log_message "Previously submitted: ${#SUBMITTED_JOBS[@]} jobs"
    log_message ""
fi

# ============================================================================
# Main monitoring loop
# ============================================================================

log_message "Starting queue monitoring loop..."
log_message "Press Ctrl+C to stop (already submitted jobs will continue)"
log_message ""

iteration=0

while true; do
    iteration=$((iteration + 1))
    
    # Get current queue status
    QUEUE_COUNT=$(get_queue_count)
    RUNNING_COUNT=$(get_running_count)
    PENDING_COUNT=$(get_pending_count)
    SUBMITTED_COUNT=${#SUBMITTED_JOBS[@]}
    REMAINING=$((TOTAL_JOBS - SUBMITTED_COUNT))
    
    log_message "=========================================="
    log_message "Iteration #$iteration - Queue Status"
    log_message "=========================================="
    log_message "Total jobs: $TOTAL_JOBS"
    log_message "Already submitted: $SUBMITTED_COUNT"
    log_message "Remaining: $REMAINING"
    log_message "Current queue: $QUEUE_COUNT (Running: $RUNNING_COUNT, Pending: $PENDING_COUNT)"
    log_message ""
    
    # Check if all jobs submitted
    if [ "$REMAINING" -eq 0 ]; then
        log_message "✅ All jobs have been submitted!"
        log_message "Monitor remaining jobs with: squeue -u $USER -n AF3"
        log_message ""
        log_message "Final statistics:"
        log_message "  Total submitted: $SUBMITTED_COUNT"
        log_message "  Currently in queue: $QUEUE_COUNT"
        log_message ""
        log_message "Queue manager finished at $(date)"
        exit 0
    fi
    
    # Check if we should submit more jobs
    if [ "$QUEUE_COUNT" -lt "$MIN_QUEUE_SIZE" ]; then
        # Calculate how many to submit
        AVAILABLE_SLOTS=$((MAX_QUEUE_SIZE - QUEUE_COUNT))
        JOBS_TO_SUBMIT=$SUBMIT_BATCH_SIZE
        
        if [ "$JOBS_TO_SUBMIT" -gt "$AVAILABLE_SLOTS" ]; then
            JOBS_TO_SUBMIT=$AVAILABLE_SLOTS
        fi
        
        if [ "$JOBS_TO_SUBMIT" -gt "$REMAINING" ]; then
            JOBS_TO_SUBMIT=$REMAINING
        fi
        
        log_message "Queue size ($QUEUE_COUNT) below minimum ($MIN_QUEUE_SIZE)"
        log_message "Submitting $JOBS_TO_SUBMIT new jobs..."
        log_message ""
        
        submitted_now=0
        
        for json_file in "${ALL_JSON_FILES[@]}"; do
            # Skip if already submitted
            if [ -n "${SUBMITTED_JOBS[$json_file]:-}" ]; then
                continue
            fi
            
            # Submit job
            JOB_DIR="$(dirname "$json_file")"
            JOB_NAME="$(basename "$JOB_DIR")"
            OUTPUT_LOG="${LOG_DIR}/${JOB_NAME}.out"
            ERROR_LOG="${LOG_DIR}/${JOB_NAME}.err"
            
            if sbatch --output="$OUTPUT_LOG" --error="$ERROR_LOG" "$AF3_SCRIPT" "$json_file" > /dev/null 2>&1; then
                SUBMITTED_JOBS["$json_file"]=1
                echo "$json_file" >> "$SUBMITTED_JOBS_FILE"
                submitted_now=$((submitted_now + 1))
                log_message "  ✅ [$submitted_now/$JOBS_TO_SUBMIT] Submitted: $JOB_NAME"
                
                if [ "$submitted_now" -ge "$JOBS_TO_SUBMIT" ]; then
                    break
                fi
            else
                log_message "  ❌ Failed to submit: $JOB_NAME"
            fi
        done
        
        log_message ""
        log_message "Submitted $submitted_now jobs in this batch"
        log_message ""
    else
        log_message "Queue size ($QUEUE_COUNT) is adequate (min: $MIN_QUEUE_SIZE, max: $MAX_QUEUE_SIZE)"
        log_message "No submission needed"
        log_message ""
    fi
    
    # Wait before next check
    log_message "Next check in ${CHECK_INTERVAL}s..."
    log_message ""
    sleep "$CHECK_INTERVAL"
done
