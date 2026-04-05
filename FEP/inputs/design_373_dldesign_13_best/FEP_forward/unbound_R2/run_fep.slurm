#!/bin/bash
#SBATCH -J R2fU
#SBATCH --partition=quick
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node24,node26,node23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:2
#SBATCH --time=12:00:00
#SBATCH --output="forward_unbound_R2.log"
#SBATCH --open-mode=append

set -o pipefail    # pipefail only; set -u added AFTER sourcing profile scripts (profile.d scripts
                   # often reference unset vars; -u during source causes silent job exit)

#####################################################
# Environment
#####################################################
source /public/software/profile.d/apps_gromacs_2024.2_mpi.sh 2>/dev/null \
    || source /public/software/profile.d/apps_gromacs_2023.2_mpi.sh 2>/dev/null \
    || { echo "ERROR: Cannot load GROMACS. Adjust GROMACS_PROFILE in submit_fep.sh"; exit 1; }
source /public/software/profile.d/mpi_openmpi-intel-2.1.2.sh 2>/dev/null || true
export GMXLIB="/public/home/xuziyi/FEP/force_fields/mutff"
export OMP_NUM_THREADS=1

set -u   # Unbound-variable check enabled HERE (after all source commands)

#####################################################
# Job variables
#####################################################
NUMBER_OF_WINDOWS=32
STARTING_FOLDER="$(pwd)"
STRUCTURE_FILE="${STARTING_FOLDER}/pHLA_forw_mutated_ff_ions.gro"
TOPOLOGY_FILE="${STARTING_FOLDER}/topol_hybrid.top"
INDEX_FILE="${STARTING_FOLDER}/index.ndx"
MDPs_FOLDER="/public/home/xuziyi/FEP/mdps"
mdp_version="_fep"

echo "========================================================"
echo "Job:        R2fU  [forward/unbound_R2]"
echo "Start:      $(date)"
echo "Node:       $(hostname)"
echo "GPUs:       ${CUDA_VISIBLE_DEVICES:-NOT_SET}"
echo "Dir:        $(pwd)"
echo "GMXLIB:     $GMXLIB"
echo "GMX binary: $(which gmx_mpi 2>/dev/null || echo NOT FOUND)"
echo "GMX version:$(gmx_mpi --version 2>&1 | grep 'GROMACS version' || true)"
echo "MPI tasks:  ${NUMBER_OF_WINDOWS} (oversubscribed on 32 CPUs)"
echo "OMP/rank:   ${OMP_NUM_THREADS}"
echo "Structure:  ${STRUCTURE_FILE}"
echo "Topology:   ${TOPOLOGY_FILE}"
echo "Index:      ${INDEX_FILE}"
echo "MDPs dir:   ${MDPs_FOLDER}"
echo "Disk (start):$(du -sh . 2>/dev/null | cut -f1)"
echo "========================================================"

window_start=0
window_end=$((NUMBER_OF_WINDOWS - 1))
directories="$(seq "$window_start" "$window_end" | awk '{printf "lambda%s ", $0}')"
num_gpus=$(echo "${CUDA_VISIBLE_DEVICES:-0}" | tr ',' '\n' | wc -l)
num_gpus_minus_one=$((num_gpus - 1))
gpu_ids="$(seq 0 "$num_gpus_minus_one" | awk '{printf "%s", $0}')"
echo "GPU ids for mdrun: ${gpu_ids}"

# Helper: verify all windows have a given file, print missing ones
check_outputs() {
    local _label="$1" _file="$2" _missing=0
    for _i in $(seq "$window_start" "$window_end"); do
        [[ -s "lambda${_i}/${_file}" ]] || { echo "  MISSING: lambda${_i}/${_file}"; ((_missing++)); }
    done
    if [[ "$_missing" -eq 0 ]]; then
        echo "$(date) CHECK OK: all ${NUMBER_OF_WINDOWS} windows have ${_label}"
    else
        echo "$(date) CHECK FAIL: ${_missing} window(s) missing ${_label}" >&2
        exit 1
    fi
}

# Helper: run grompp for one window; log to per-window file; dump tail on error
run_grompp() {
    local _tag="$1" _win="$2"; shift 2
    local _log="lambda${_win}/grompp_${_tag}.log"
    if gmx_mpi grompp "$@" &> "${_log}"; then
        :   # success: log is kept for reference but not echoed
    else
        echo "ERROR: grompp ${_tag} lambda${_win} failed. Last 30 lines of ${_log}:" >&2
        tail -30 "${_log}" >&2
        exit 1
    fi
}

#####################################################
# STEP 1: Lambda directories + MDP files
#####################################################
echo ""
echo "--- STEP 1: Lambda dirs and MDP files  $(date) ---"
for window_idx in $(seq "$window_start" "$window_end"); do
    mkdir -p "lambda${window_idx}/EM" "lambda${window_idx}/NVT" \
             "lambda${window_idx}/NPT" "lambda${window_idx}/PROD"
    sed "s/KEY_win_idx/${window_idx}/g" "${MDPs_FOLDER}/em${mdp_version}.mdp" \
        > "lambda${window_idx}/em.mdp"    || { echo "ERROR: em mdp ${window_idx}";    exit 1; }
    sed "s/steep/cg/g" "lambda${window_idx}/em.mdp" \
        > "lambda${window_idx}/em_cg.mdp" || { echo "ERROR: em_cg mdp ${window_idx}"; exit 1; }
    sed "s/KEY_win_idx/${window_idx}/g" "${MDPs_FOLDER}/nvt${mdp_version}.mdp" \
        > "lambda${window_idx}/nvt.mdp"   || { echo "ERROR: nvt mdp ${window_idx}";   exit 1; }
    sed "s/KEY_win_idx/${window_idx}/g" "${MDPs_FOLDER}/npt${mdp_version}.mdp" \
        > "lambda${window_idx}/npt.mdp"   || { echo "ERROR: npt mdp ${window_idx}";   exit 1; }
    sed "s/KEY_win_idx/${window_idx}/g" "${MDPs_FOLDER}/prod${mdp_version}.mdp" \
        > "lambda${window_idx}/prod.mdp"  || { echo "ERROR: prod mdp ${window_idx}";  exit 1; }
done
echo "$(date) STEP 1 done: MDP files written for ${NUMBER_OF_WINDOWS} windows"

#####################################################
# STEP 2: Energy Minimization (3-pass: steep -> CG -> steep)
#####################################################
echo ""
echo "--- STEP 2: Energy Minimization  $(date) ---"
if [[ -e "lambda${window_end}/EM/em.gro" ]]; then
    echo "$(date) EM already complete, skipping."
else
    # Pass 1: steep
    echo "$(date) EM pass 1/3 (steep): grompp..."
    for window_idx in $(seq "$window_start" "$window_end"); do
        run_grompp em1 "${window_idx}" \
            -f "lambda${window_idx}/em.mdp" \
            -c "${STRUCTURE_FILE}" -r "${STRUCTURE_FILE}" \
            -p "${TOPOLOGY_FILE}" -n "${INDEX_FILE}" \
            -o "lambda${window_idx}/em1.tpr" -maxwarn 3
    done
    echo "$(date) EM pass 1/3 (steep): mdrun..."
    mpirun -oversubscribe -np ${NUMBER_OF_WINDOWS} gmx_mpi mdrun \
        -v -deffnm em1 -multidir $(echo "${directories}") -gpu_id "${gpu_ids}" \
        &>> em1.out || { echo "ERROR: mdrun EM pass 1. See em1.out"; tail -40 em1.out >&2; exit 1; }
    check_outputs "em1.gro" "em1.gro"

    # Pass 2: CG
    echo "$(date) EM pass 2/3 (CG): grompp..."
    for window_idx in $(seq "$window_start" "$window_end"); do
        run_grompp em_cg "${window_idx}" \
            -f "lambda${window_idx}/em_cg.mdp" \
            -c "lambda${window_idx}/em1.gro" -r "lambda${window_idx}/em1.gro" \
            -p "${TOPOLOGY_FILE}" -n "${INDEX_FILE}" \
            -o "lambda${window_idx}/em_cg.tpr" -maxwarn 3
    done
    echo "$(date) EM pass 2/3 (CG): mdrun..."
    mpirun -oversubscribe -np ${NUMBER_OF_WINDOWS} gmx_mpi mdrun \
        -v -deffnm em_cg -multidir $(echo "${directories}") -gpu_id "${gpu_ids}" \
        &>> em_cg.out || { echo "ERROR: mdrun EM pass 2. See em_cg.out"; tail -40 em_cg.out >&2; exit 1; }
    check_outputs "em_cg.gro" "em_cg.gro"

    # Pass 3: steep from CG
    echo "$(date) EM pass 3/3 (steep from CG): grompp..."
    for window_idx in $(seq "$window_start" "$window_end"); do
        run_grompp em3 "${window_idx}" \
            -f "lambda${window_idx}/em.mdp" \
            -c "lambda${window_idx}/em_cg.gro" -r "lambda${window_idx}/em_cg.gro" \
            -p "${TOPOLOGY_FILE}" -n "${INDEX_FILE}" \
            -o "lambda${window_idx}/em.tpr" -maxwarn 3
    done
    echo "$(date) EM pass 3/3 (steep): mdrun..."
    mpirun -oversubscribe -np ${NUMBER_OF_WINDOWS} gmx_mpi mdrun \
        -v -deffnm em -multidir $(echo "${directories}") -gpu_id "${gpu_ids}" \
        &>> em.out || { echo "ERROR: mdrun EM pass 3. See em.out"; tail -40 em.out >&2; exit 1; }
    check_outputs "em.gro" "em.gro"

    for window_idx in $(seq "$window_start" "$window_end"); do
        mv ./lambda${window_idx}/em* ./lambda${window_idx}/EM/
    done
    echo "$(date) STEP 2 done: EM complete. Disk: $(du -sh . | cut -f1)"
fi

#####################################################
# STEP 3: NVT Equilibration (300 ps)
#####################################################
echo ""
echo "--- STEP 3: NVT Equilibration  $(date) ---"
if [[ -e "lambda${window_end}/NVT/nvt.gro" ]]; then
    echo "$(date) NVT already complete, skipping."
else
    echo "$(date) NVT: grompp..."
    for window_idx in $(seq "$window_start" "$window_end"); do
        run_grompp nvt "${window_idx}" \
            -f "lambda${window_idx}/nvt.mdp" \
            -c "lambda${window_idx}/EM/em.gro" -r "lambda${window_idx}/EM/em.gro" \
            -p "${TOPOLOGY_FILE}" -n "${INDEX_FILE}" \
            -o "lambda${window_idx}/nvt.tpr" -maxwarn 3
    done
    echo "$(date) NVT: mdrun..."
    mpirun -oversubscribe -np ${NUMBER_OF_WINDOWS} gmx_mpi mdrun \
        -v -deffnm nvt -nb gpu -bonded gpu \
        -multidir $(echo "${directories}") -gpu_id "${gpu_ids}" \
        &>> nvt.out || { echo "ERROR: mdrun NVT. See nvt.out"; tail -40 nvt.out >&2; exit 1; }
    check_outputs "nvt.gro" "nvt.gro"
    for window_idx in $(seq "$window_start" "$window_end"); do
        mv ./lambda${window_idx}/nvt* ./lambda${window_idx}/NVT/
    done
    # Quick temperature sanity: print last Temp line from lambda0 log
    echo "$(date) NVT temperature check (lambda0):"
    grep -E "^\s+[0-9]" lambda0/NVT/nvt.log 2>/dev/null | awk 'END{printf "  Last log line: %s\n", $0}' || true
    echo "$(date) STEP 3 done: NVT complete. Disk: $(du -sh . | cut -f1)"
fi

#####################################################
# STEP 4: NPT Equilibration (500 ps, Berendsen)
#####################################################
echo ""
echo "--- STEP 4: NPT Equilibration  $(date) ---"
if [[ -e "lambda${window_end}/NPT/npt.gro" ]]; then
    echo "$(date) NPT already complete, skipping."
else
    echo "$(date) NPT: grompp..."
    for window_idx in $(seq "$window_start" "$window_end"); do
        run_grompp npt "${window_idx}" \
            -f "lambda${window_idx}/npt.mdp" \
            -c "lambda${window_idx}/NVT/nvt.gro" -r "lambda${window_idx}/NVT/nvt.gro" \
            -p "${TOPOLOGY_FILE}" -n "${INDEX_FILE}" \
            -o "lambda${window_idx}/npt.tpr" -maxwarn 4
    done
    echo "$(date) NPT: mdrun..."
    mpirun -oversubscribe -np ${NUMBER_OF_WINDOWS} gmx_mpi mdrun \
        -v -deffnm npt -nb gpu -bonded gpu \
        -multidir $(echo "${directories}") -gpu_id "${gpu_ids}" \
        &>> npt.out || { echo "ERROR: mdrun NPT. See npt.out"; tail -40 npt.out >&2; exit 1; }
    check_outputs "npt.gro" "npt.gro"
    for window_idx in $(seq "$window_start" "$window_end"); do
        mv ./lambda${window_idx}/npt* ./lambda${window_idx}/NPT/
    done
    echo "$(date) STEP 4 done: NPT complete. Disk: $(du -sh . | cut -f1)"
fi

#####################################################
# STEP 5: FEP Production HREX (1.5 ns, C-rescale, replex 1000)
#####################################################
echo ""
echo "--- STEP 5: PROD HREX  $(date) ---"
if [[ -e "lambda${window_end}/PROD/prod.gro" ]]; then
    echo "$(date) PROD already complete."
else
    echo "$(date) PROD: grompp..."
    for window_idx in $(seq "$window_start" "$window_end"); do
        run_grompp prod "${window_idx}" \
            -f "lambda${window_idx}/prod.mdp" \
            -c "lambda${window_idx}/NPT/npt.gro" -r "lambda${window_idx}/NVT/nvt.gro" \
            -p "${TOPOLOGY_FILE}" -n "${INDEX_FILE}" \
            -o "lambda${window_idx}/prod.tpr" -maxwarn 4
    done
    echo "$(date) PROD: mdrun HREX (replex=1000)..."
    mpirun -oversubscribe -np ${NUMBER_OF_WINDOWS} gmx_mpi mdrun \
        -v -deffnm prod -cpi prod.cpt -nb gpu -bonded gpu \
        -replex 1000 -multidir $(echo "${directories}") \
        -gpu_id "${gpu_ids}" -pin on \
        &>> prod.out || { echo "ERROR: mdrun PROD. See prod.out"; tail -40 prod.out >&2; exit 1; }
    check_outputs "prod.gro" "prod.gro"
    for window_idx in $(seq "$window_start" "$window_end"); do
        mkdir -p ./lambda${window_idx}/PROD
        mv ./lambda${window_idx}/prod* ./lambda${window_idx}/dhdl* \
           ./lambda${window_idx}/PROD/ 2>/dev/null || true
    done
    echo "$(date) STEP 5 done: PROD complete."
fi

#####################################################
# Final summary
#####################################################
echo ""
echo "========================================================"
echo "FEP HREX finished: $(date)"
echo "Total disk usage:  $(du -sh . | cut -f1)"
echo ""
echo "Output file check (prod.gro per lambda):"
ok=0; missing=0
for _i in $(seq "$window_start" "$window_end"); do
    _f="lambda${_i}/PROD/prod.gro"
    if [[ -s "${_f}" ]]; then
        printf "  lambda%2d: OK  (%s)\n" "${_i}" "$(ls -lh "${_f}" | awk '{print $5}')"
        ((ok++))
    else
        printf "  lambda%2d: MISSING\n" "${_i}"
        ((missing++))
    fi
done
echo "  -> ${ok}/${NUMBER_OF_WINDOWS} windows have prod.gro  (${missing} missing)"
echo ""
echo "dhdl file check (dhdl.xvg per lambda):"
ok_dh=0
for _i in $(seq "$window_start" "$window_end"); do
    _f="lambda${_i}/PROD/dhdl.xvg"
    [[ -s "${_f}" ]] && ((ok_dh++)) || printf "  lambda%2d: dhdl.xvg MISSING\n" "${_i}"
done
echo "  -> ${ok_dh}/${NUMBER_OF_WINDOWS} windows have dhdl.xvg"
echo "========================================================"
