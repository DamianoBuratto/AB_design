#!/bin/bash
#####################################################
# submit_fep.sh - Submit FEP HREX jobs to SLURM
#
# Run ON HPC after setup_fep.sh has completed.
#
# Usage:
#   bash scripts/submit_fep.sh                         # forward only
#   bash scripts/submit_fep.sh reverse                 # reverse only
#   bash scripts/submit_fep.sh forward reverse         # both
#####################################################

#===============================================================
# CONFIGURATION - Must match setup_fep.sh settings
#===============================================================

# --- HPC paths ---
FEP_ROOT="/public/home/xuziyi/FEP"
GMXLIB_PATH="/public/home/xuziyi/FEP/force_fields/mutff"

# --- Design ---
DESIGN_NAME="design_373_dldesign_13_best"
NR_REPLICAS=3

# --- SLURM resources ---
SLURM_PARTITION="quick"
SLURM_EXCLUDE="node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node24,node26,node23"
SLURM_TIME="12:00:00"
SLURM_NODES=1
SLURM_NTASKS=1
SLURM_CPUS=32
SLURM_GPUS=2

# --- GROMACS / environment (written into each SLURM script) ---
GROMACS_PROFILE_1="/public/software/profile.d/apps_gromacs_2024.2_mpi.sh"
GROMACS_PROFILE_2="/public/software/profile.d/apps_gromacs_2023.2_mpi.sh"
MPI_PROFILE="/public/software/profile.d/mpi_openmpi-intel-2.1.2.sh"
OMP_THREADS=1    # 32 MPI ranks on 32 CPUs (1 rank/core); 1 OMP thread per rank

#===============================================================
# DERIVED PATHS
#===============================================================
OUTPUTS_DIR="${FEP_ROOT}/outputs/${DESIGN_NAME}"
MDP_DIR="${FEP_ROOT}/mdps"
MDP_VERSION="_fep"   # suffix: em_fep.mdp, nvt_fep.mdp, etc.

#===============================================================
# ARGUMENT HANDLING
# Default: forward only.
# --with-reverse   shortcut to submit both forward and reverse
# Positional args: e.g.  forward reverse
#===============================================================
DIRECTIONS="forward"
for _arg in "$@"; do
    case "$_arg" in
        --with-reverse)    DIRECTIONS="forward reverse" ;;
        forward|reverse)   [[ "$DIRECTIONS" != *"$_arg"* ]] && DIRECTIONS="$DIRECTIONS $_arg" ;;
        *) echo "WARNING: Unknown argument '$_arg', ignoring." ;;
    esac
done
DIRECTIONS="${DIRECTIONS## }"  # strip leading space

[[ -d "$OUTPUTS_DIR" ]] || { echo "ERROR: Outputs dir not found: $OUTPUTS_DIR"; exit 1; }
[[ -d "$MDP_DIR"     ]] || { echo "ERROR: MDP dir not found: $MDP_DIR"; exit 1; }

echo "============================================"
echo " FEP Submit"
echo " Design:    ${DESIGN_NAME}"
echo " Outputs:   ${OUTPUTS_DIR}"
echo " Direction: ${DIRECTIONS}"
echo " Replicas:  ${NR_REPLICAS}"
echo "============================================"

start_dir="$PWD"

for direction in ${DIRECTIONS}; do
    [[ "$direction" == "forward" || "$direction" == "reverse" ]] \
        || { echo "WARNING: Unknown direction '$direction', skipping."; continue; }

    for state in bound unbound; do
        for replica in $(seq 1 ${NR_REPLICAS}); do

            state_dir="${OUTPUTS_DIR}/${direction}/${state}_R${replica}"

            if [[ ! -d "$state_dir" ]]; then
                echo "SKIP: ${direction}/${state}_R${replica} (directory not found)"
                continue
            fi

            cd "$state_dir" || { echo "ERROR: Cannot cd to ${state_dir}"; continue; }

            # Skip if production already complete
            if [[ -e "lambda31/PROD/prod.gro" ]]; then
                echo "DONE: ${direction}/${state}_R${replica} (prod.gro exists)"
                cd "$start_dir"; continue
            fi

            # Find starting GRO (*ions.gro written by setup_fep.sh)
            gro_file=$(ls *ions.gro 2>/dev/null | head -1)
            if [[ -z "$gro_file" ]]; then
                echo "SKIP: ${direction}/${state}_R${replica} (no *ions.gro found)"
                cd "$start_dir"; continue
            fi

            # Job name (max 15 chars): R<n><d><S>  e.g. R1fB
            dir_code="${direction:0:1}"
            [[ "$state" == "bound" ]] && state_code="B" || state_code="U"
            job_name="R${replica}${dir_code}${state_code}"
            output_log="${direction}_${state}_R${replica}.log"

            #-----------------------------------------------------------
            # Write SLURM script inline
            # Note: \$ escapes prevent expansion NOW; they expand at job runtime.
            #-----------------------------------------------------------
            cat > run_fep.slurm << SLURM_EOF
#!/bin/bash
#SBATCH -J ${job_name}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --exclude=${SLURM_EXCLUDE}
#SBATCH --nodes=${SLURM_NODES}
#SBATCH --ntasks=${SLURM_NTASKS}
#SBATCH --cpus-per-task=${SLURM_CPUS}
#SBATCH --gres=gpu:${SLURM_GPUS}
#SBATCH --time=${SLURM_TIME}
#SBATCH --output="${output_log}"
#SBATCH --open-mode=append

set -o pipefail    # pipefail only; set -u added AFTER sourcing profile scripts (profile.d scripts
                   # often reference unset vars; -u during source causes silent job exit)

#####################################################
# Environment
#####################################################
source ${GROMACS_PROFILE_1} 2>/dev/null \\
    || source ${GROMACS_PROFILE_2} 2>/dev/null \\
    || { echo "ERROR: Cannot load GROMACS. Adjust GROMACS_PROFILE in submit_fep.sh"; exit 1; }
source ${MPI_PROFILE} 2>/dev/null || true
export GMXLIB="${GMXLIB_PATH}"
export OMP_NUM_THREADS=${OMP_THREADS}

set -u   # Unbound-variable check enabled HERE (after all source commands)

#####################################################
# Job variables
#####################################################
NUMBER_OF_WINDOWS=32
STARTING_FOLDER="\$(pwd)"
STRUCTURE_FILE="\${STARTING_FOLDER}/${gro_file}"
TOPOLOGY_FILE="\${STARTING_FOLDER}/topol_hybrid.top"
INDEX_FILE="\${STARTING_FOLDER}/index.ndx"
MDPs_FOLDER="${MDP_DIR}"
mdp_version="${MDP_VERSION}"

echo "========================================================"
echo "Job:        ${job_name}  [${direction}/${state}_R${replica}]"
echo "Start:      \$(date)"
echo "Node:       \$(hostname)"
echo "GPUs:       \${CUDA_VISIBLE_DEVICES:-NOT_SET}"
echo "Dir:        \$(pwd)"
echo "GMXLIB:     \$GMXLIB"
echo "GMX binary: \$(which gmx_mpi 2>/dev/null || echo NOT FOUND)"
echo "GMX version:\$(gmx_mpi --version 2>&1 | grep 'GROMACS version' || true)"
echo "MPI tasks:  \${NUMBER_OF_WINDOWS} (oversubscribed on ${SLURM_CPUS} CPUs)"
echo "OMP/rank:   \${OMP_NUM_THREADS}"
echo "Structure:  \${STRUCTURE_FILE}"
echo "Topology:   \${TOPOLOGY_FILE}"
echo "Index:      \${INDEX_FILE}"
echo "MDPs dir:   \${MDPs_FOLDER}"
echo "Disk (start):\$(du -sh . 2>/dev/null | cut -f1)"
echo "========================================================"

window_start=0
window_end=\$((NUMBER_OF_WINDOWS - 1))
directories="\$(seq "\$window_start" "\$window_end" | awk '{printf "lambda%s ", \$0}')"
num_gpus=\$(echo "\${CUDA_VISIBLE_DEVICES:-0}" | tr ',' '\n' | wc -l)
num_gpus_minus_one=\$((num_gpus - 1))
gpu_ids="\$(seq 0 "\$num_gpus_minus_one" | awk '{printf "%s", \$0}')"
echo "GPU ids for mdrun: \${gpu_ids}"

# Helper: verify all windows have a given file, print missing ones
check_outputs() {
    local _label="\$1" _file="\$2" _missing=0
    for _i in \$(seq "\$window_start" "\$window_end"); do
        [[ -s "lambda\${_i}/\${_file}" ]] || { echo "  MISSING: lambda\${_i}/\${_file}"; ((_missing++)); }
    done
    if [[ "\$_missing" -eq 0 ]]; then
        echo "\$(date) CHECK OK: all \${NUMBER_OF_WINDOWS} windows have \${_label}"
    else
        echo "\$(date) CHECK FAIL: \${_missing} window(s) missing \${_label}" >&2
        exit 1
    fi
}

# Helper: run grompp for one window; log to per-window file; dump tail on error
run_grompp() {
    local _tag="\$1" _win="\$2"; shift 2
    local _log="lambda\${_win}/grompp_\${_tag}.log"
    if gmx_mpi grompp "\$@" &> "\${_log}"; then
        :   # success: log is kept for reference but not echoed
    else
        echo "ERROR: grompp \${_tag} lambda\${_win} failed. Last 30 lines of \${_log}:" >&2
        tail -30 "\${_log}" >&2
        exit 1
    fi
}

#####################################################
# STEP 1: Lambda directories + MDP files
#####################################################
echo ""
echo "--- STEP 1: Lambda dirs and MDP files  \$(date) ---"
for window_idx in \$(seq "\$window_start" "\$window_end"); do
    mkdir -p "lambda\${window_idx}/EM" "lambda\${window_idx}/NVT" \\
             "lambda\${window_idx}/NPT" "lambda\${window_idx}/PROD"
    sed "s/KEY_win_idx/\${window_idx}/g" "\${MDPs_FOLDER}/em\${mdp_version}.mdp" \\
        > "lambda\${window_idx}/em.mdp"    || { echo "ERROR: em mdp \${window_idx}";    exit 1; }
    sed "s/steep/cg/g" "lambda\${window_idx}/em.mdp" \\
        > "lambda\${window_idx}/em_cg.mdp" || { echo "ERROR: em_cg mdp \${window_idx}"; exit 1; }
    sed "s/KEY_win_idx/\${window_idx}/g" "\${MDPs_FOLDER}/nvt\${mdp_version}.mdp" \\
        > "lambda\${window_idx}/nvt.mdp"   || { echo "ERROR: nvt mdp \${window_idx}";   exit 1; }
    sed "s/KEY_win_idx/\${window_idx}/g" "\${MDPs_FOLDER}/npt\${mdp_version}.mdp" \\
        > "lambda\${window_idx}/npt.mdp"   || { echo "ERROR: npt mdp \${window_idx}";   exit 1; }
    sed "s/KEY_win_idx/\${window_idx}/g" "\${MDPs_FOLDER}/prod\${mdp_version}.mdp" \\
        > "lambda\${window_idx}/prod.mdp"  || { echo "ERROR: prod mdp \${window_idx}";  exit 1; }
done
echo "\$(date) STEP 1 done: MDP files written for \${NUMBER_OF_WINDOWS} windows"

#####################################################
# STEP 2: Energy Minimization (3-pass: steep -> CG -> steep)
#####################################################
echo ""
echo "--- STEP 2: Energy Minimization  \$(date) ---"
if [[ -e "lambda\${window_end}/EM/em.gro" ]]; then
    echo "\$(date) EM already complete, skipping."
else
    # Pass 1: steep
    echo "\$(date) EM pass 1/3 (steep): grompp..."
    for window_idx in \$(seq "\$window_start" "\$window_end"); do
        run_grompp em1 "\${window_idx}" \\
            -f "lambda\${window_idx}/em.mdp" \\
            -c "\${STRUCTURE_FILE}" -r "\${STRUCTURE_FILE}" \\
            -p "\${TOPOLOGY_FILE}" -n "\${INDEX_FILE}" \\
            -o "lambda\${window_idx}/em1.tpr" -maxwarn 3
    done
    echo "\$(date) EM pass 1/3 (steep): mdrun..."
    mpirun -oversubscribe -np \${NUMBER_OF_WINDOWS} gmx_mpi mdrun \\
        -v -deffnm em1 -multidir \$(echo "\${directories}") -gpu_id "\${gpu_ids}" \\
        &>> em1.out || { echo "ERROR: mdrun EM pass 1. See em1.out"; tail -40 em1.out >&2; exit 1; }
    check_outputs "em1.gro" "em1.gro"

    # Pass 2: CG
    echo "\$(date) EM pass 2/3 (CG): grompp..."
    for window_idx in \$(seq "\$window_start" "\$window_end"); do
        run_grompp em_cg "\${window_idx}" \\
            -f "lambda\${window_idx}/em_cg.mdp" \\
            -c "lambda\${window_idx}/em1.gro" -r "lambda\${window_idx}/em1.gro" \\
            -p "\${TOPOLOGY_FILE}" -n "\${INDEX_FILE}" \\
            -o "lambda\${window_idx}/em_cg.tpr" -maxwarn 3
    done
    echo "\$(date) EM pass 2/3 (CG): mdrun..."
    mpirun -oversubscribe -np \${NUMBER_OF_WINDOWS} gmx_mpi mdrun \\
        -v -deffnm em_cg -multidir \$(echo "\${directories}") -gpu_id "\${gpu_ids}" \\
        &>> em_cg.out || { echo "ERROR: mdrun EM pass 2. See em_cg.out"; tail -40 em_cg.out >&2; exit 1; }
    check_outputs "em_cg.gro" "em_cg.gro"

    # Pass 3: steep from CG
    echo "\$(date) EM pass 3/3 (steep from CG): grompp..."
    for window_idx in \$(seq "\$window_start" "\$window_end"); do
        run_grompp em3 "\${window_idx}" \\
            -f "lambda\${window_idx}/em.mdp" \\
            -c "lambda\${window_idx}/em_cg.gro" -r "lambda\${window_idx}/em_cg.gro" \\
            -p "\${TOPOLOGY_FILE}" -n "\${INDEX_FILE}" \\
            -o "lambda\${window_idx}/em.tpr" -maxwarn 3
    done
    echo "\$(date) EM pass 3/3 (steep): mdrun..."
    mpirun -oversubscribe -np \${NUMBER_OF_WINDOWS} gmx_mpi mdrun \\
        -v -deffnm em -multidir \$(echo "\${directories}") -gpu_id "\${gpu_ids}" \\
        &>> em.out || { echo "ERROR: mdrun EM pass 3. See em.out"; tail -40 em.out >&2; exit 1; }
    check_outputs "em.gro" "em.gro"

    for window_idx in \$(seq "\$window_start" "\$window_end"); do
        mv ./lambda\${window_idx}/em* ./lambda\${window_idx}/EM/
    done
    echo "\$(date) STEP 2 done: EM complete. Disk: \$(du -sh . | cut -f1)"
fi

#####################################################
# STEP 3: NVT Equilibration (300 ps)
#####################################################
echo ""
echo "--- STEP 3: NVT Equilibration  \$(date) ---"
if [[ -e "lambda\${window_end}/NVT/nvt.gro" ]]; then
    echo "\$(date) NVT already complete, skipping."
else
    echo "\$(date) NVT: grompp..."
    for window_idx in \$(seq "\$window_start" "\$window_end"); do
        run_grompp nvt "\${window_idx}" \\
            -f "lambda\${window_idx}/nvt.mdp" \\
            -c "lambda\${window_idx}/EM/em.gro" -r "lambda\${window_idx}/EM/em.gro" \\
            -p "\${TOPOLOGY_FILE}" -n "\${INDEX_FILE}" \\
            -o "lambda\${window_idx}/nvt.tpr" -maxwarn 3
    done
    echo "\$(date) NVT: mdrun..."
    mpirun -oversubscribe -np \${NUMBER_OF_WINDOWS} gmx_mpi mdrun \\
        -v -deffnm nvt -nb gpu -bonded gpu \\
        -multidir \$(echo "\${directories}") -gpu_id "\${gpu_ids}" \\
        &>> nvt.out || { echo "ERROR: mdrun NVT. See nvt.out"; tail -40 nvt.out >&2; exit 1; }
    check_outputs "nvt.gro" "nvt.gro"
    for window_idx in \$(seq "\$window_start" "\$window_end"); do
        mv ./lambda\${window_idx}/nvt* ./lambda\${window_idx}/NVT/
    done
    # Quick temperature sanity: print last Temp line from lambda0 log
    echo "\$(date) NVT temperature check (lambda0):"
    grep -E "^\s+[0-9]" lambda0/NVT/nvt.log 2>/dev/null | awk 'END{printf "  Last log line: %s\n", \$0}' || true
    echo "\$(date) STEP 3 done: NVT complete. Disk: \$(du -sh . | cut -f1)"
fi

#####################################################
# STEP 4: NPT Equilibration (500 ps, Berendsen)
#####################################################
echo ""
echo "--- STEP 4: NPT Equilibration  \$(date) ---"
if [[ -e "lambda\${window_end}/NPT/npt.gro" ]]; then
    echo "\$(date) NPT already complete, skipping."
else
    echo "\$(date) NPT: grompp..."
    for window_idx in \$(seq "\$window_start" "\$window_end"); do
        run_grompp npt "\${window_idx}" \\
            -f "lambda\${window_idx}/npt.mdp" \\
            -c "lambda\${window_idx}/NVT/nvt.gro" -r "lambda\${window_idx}/NVT/nvt.gro" \\
            -p "\${TOPOLOGY_FILE}" -n "\${INDEX_FILE}" \\
            -o "lambda\${window_idx}/npt.tpr" -maxwarn 4
    done
    echo "\$(date) NPT: mdrun..."
    mpirun -oversubscribe -np \${NUMBER_OF_WINDOWS} gmx_mpi mdrun \\
        -v -deffnm npt -nb gpu -bonded gpu \\
        -multidir \$(echo "\${directories}") -gpu_id "\${gpu_ids}" \\
        &>> npt.out || { echo "ERROR: mdrun NPT. See npt.out"; tail -40 npt.out >&2; exit 1; }
    check_outputs "npt.gro" "npt.gro"
    for window_idx in \$(seq "\$window_start" "\$window_end"); do
        mv ./lambda\${window_idx}/npt* ./lambda\${window_idx}/NPT/
    done
    echo "\$(date) STEP 4 done: NPT complete. Disk: \$(du -sh . | cut -f1)"
fi

#####################################################
# STEP 5: FEP Production HREX (1.5 ns, C-rescale, replex 1000)
#####################################################
echo ""
echo "--- STEP 5: PROD HREX  \$(date) ---"
if [[ -e "lambda\${window_end}/PROD/prod.gro" ]]; then
    echo "\$(date) PROD already complete."
else
    echo "\$(date) PROD: grompp..."
    for window_idx in \$(seq "\$window_start" "\$window_end"); do
        run_grompp prod "\${window_idx}" \\
            -f "lambda\${window_idx}/prod.mdp" \\
            -c "lambda\${window_idx}/NPT/npt.gro" -r "lambda\${window_idx}/NVT/nvt.gro" \\
            -p "\${TOPOLOGY_FILE}" -n "\${INDEX_FILE}" \\
            -o "lambda\${window_idx}/prod.tpr" -maxwarn 4
    done
    echo "\$(date) PROD: mdrun HREX (replex=1000)..."
    mpirun -oversubscribe -np \${NUMBER_OF_WINDOWS} gmx_mpi mdrun \\
        -v -deffnm prod -cpi prod.cpt -nb gpu -bonded gpu \\
        -replex 1000 -multidir \$(echo "\${directories}") \\
        -gpu_id "\${gpu_ids}" -pin on \\
        &>> prod.out || { echo "ERROR: mdrun PROD. See prod.out"; tail -40 prod.out >&2; exit 1; }
    check_outputs "prod.gro" "prod.gro"
    for window_idx in \$(seq "\$window_start" "\$window_end"); do
        mkdir -p ./lambda\${window_idx}/PROD
        mv ./lambda\${window_idx}/prod* ./lambda\${window_idx}/dhdl* \\
           ./lambda\${window_idx}/PROD/ 2>/dev/null || true
    done
    echo "\$(date) STEP 5 done: PROD complete."
fi

#####################################################
# Final summary
#####################################################
echo ""
echo "========================================================"
echo "FEP HREX finished: \$(date)"
echo "Total disk usage:  \$(du -sh . | cut -f1)"
echo ""
echo "Output file check (prod.gro per lambda):"
ok=0; missing=0
for _i in \$(seq "\$window_start" "\$window_end"); do
    _f="lambda\${_i}/PROD/prod.gro"
    if [[ -s "\${_f}" ]]; then
        printf "  lambda%2d: OK  (%s)\n" "\${_i}" "\$(ls -lh "\${_f}" | awk '{print \$5}')"
        ((ok++))
    else
        printf "  lambda%2d: MISSING\n" "\${_i}"
        ((missing++))
    fi
done
echo "  -> \${ok}/\${NUMBER_OF_WINDOWS} windows have prod.gro  (\${missing} missing)"
echo ""
echo "dhdl file check (dhdl.xvg per lambda):"
ok_dh=0
for _i in \$(seq "\$window_start" "\$window_end"); do
    _f="lambda\${_i}/PROD/dhdl.xvg"
    [[ -s "\${_f}" ]] && ((ok_dh++)) || printf "  lambda%2d: dhdl.xvg MISSING\n" "\${_i}"
done
echo "  -> \${ok_dh}/\${NUMBER_OF_WINDOWS} windows have dhdl.xvg"
echo "========================================================"
SLURM_EOF
            #-----------------------------------------------------------
            # Submit
            #-----------------------------------------------------------
            sbatch_out=$(sbatch ./run_fep.slurm) \
                || { echo "ERROR: sbatch failed for ${direction}/${state}_R${replica}"; cd "$start_dir"; continue; }
            job_id=$(echo "$sbatch_out" | awk '{print $NF}')
            echo "SUBMIT: ${direction}/${state}_R${replica}  [job=${job_name}, gro=${gro_file}] -> SLURM job ${job_id}"

            cd "$start_dir"
        done
    done
done

echo ""
echo "All submissions done. Monitor: squeue -u \$USER"
