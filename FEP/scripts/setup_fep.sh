#!/bin/bash
#####################################################
# setup_fep.sh - FEP Binding Free Energy Setup
# System:    Antibody-pHLA complex
# Mutation:  HIS239 -> ARG (peptide chain P, resid 8)
# Method:    Alchemical FEP (pmx + GROMACS)
# States:    BOUND (Ab+pHLA) / UNBOUND (pHLA only)
#
# Run ON HPC:
#   bash scripts/setup_fep.sh                          # forward only
#   bash scripts/setup_fep.sh --direction both         # forward + reverse skeleton
#   bash scripts/setup_fep.sh --direction reverse \
#       --rev-gro /path/bound_R1_last.gro \
#       --rev-gro-unbound /path/unbound_R1_last.gro    # reverse with starting structures
#####################################################

#===============================================================
# CONFIGURATION - Edit before running on HPC
#===============================================================

# --- HPC paths ---
FEP_ROOT="/public/home/xuziyi/FEP"               # Root of FEP directory on HPC
GMXLIB_PATH="/public/home/xuziyi/FEP/force_fields/mutff"   # mutff force field directory
CONDA_ENV="fep_env"                                # Conda environment name

# --- Design / mutation ---
DESIGN_NAME="design_373_dldesign_13_best"         # Subdirectory under inputs/
INPUT_STRUCTURE="replica_1/md_r1.gro"             # Final MD GRO (relative to inputs/DESIGN_NAME)
INPUT_TPR="replica_1/md_r1.tpr"                   # TPR matching the GRO (for gmx trjconv)
INPUT_PDB="setup/system.pdb"                       # pdb2gmx OUTPUT from original MD setup (with H, 6504 atoms = matches GRO)
MUTATION="HIS239ARG"                              # Label only (for output naming)
PEPTIDE_CHAIN="P"                                 # Chain ID of the mutated peptide
PEPTIDE_RESID="8"                                 # Residue index within the peptide chain
MUTATION_TO="R"                                   # Target amino acid (single-letter)

# --- Force field / simulation ---
FORCE_FIELD="charmm36m-mut"
WATER_MODEL="tip3p"
SIMULATION_TEMP=310                               # K
NR_REPLICAS=3

# Lambda schedule (32 windows, dense at endpoints)
LAMBDAS=(0 0.0001 0.001 0.01 0.02 0.03 0.04 0.07 0.10 0.15 0.20 0.25 0.30 0.35 \
         0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.93 0.96 0.98 \
         0.99 0.999 0.9999 1.0)

# --- Position restraints (HLA alpha-1/alpha-2 floor, chain A) ---
POSRE_CHAINS="A"
POSRE_RESIDUES="7 39 41 75 77 79 170 172 174 209 211 213"

# --- Ions (must match original MD setup) ---
ION_POS="K"
ION_NEG="CL"
ION_CONC=0.15    # M

# --- Chains to remove for UNBOUND state (antibody H+L) ---
CHAINS_TO_REMOVE="H L"

#===============================================================
# ARGUMENT PARSING
#===============================================================
DIRECTION="forward"
REV_GRO=""
REV_GRO_UNBOUND=""

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

while [[ $# -gt 0 ]]; do
    case "$1" in
        --with-reverse)     DIRECTION="both";       shift ;;   # shortcut: forward + reverse
        --direction)        DIRECTION="$2";        shift 2 ;;
        --rev-gro)          REV_GRO="$2";           shift 2 ;;
        --rev-gro-unbound)  REV_GRO_UNBOUND="$2";  shift 2 ;;
        -h|--help)
            cat << 'HELP'
Usage: bash setup_fep.sh [OPTIONS]

  --with-reverse
      Run both forward AND reverse directions (shortcut for --direction both).
      Default: forward only.
      For reverse, provide --rev-gro and --rev-gro-unbound after forward PROD is done;
      omitting them creates empty reverse dirs to be filled later.

  --direction forward|reverse|both
      Fine-grained control (use --with-reverse for the common case).

  --rev-gro <file>
      ARG-end GRO for BOUND reverse (OPTIONAL).
      If omitted, the script auto-extracts the last frame of
      forward/bound_R1/lambda{N-1}/PROD/prod.xtc. Only needed if
      you want to supply a custom starting structure.

  --rev-gro-unbound <file>
      ARG-end GRO for UNBOUND reverse (OPTIONAL; atom count DIFFERS from
      BOUND -- must be a separate file if provided manually).
      If omitted, auto-extracted from forward/unbound_R1/lambda{N-1}/PROD/prod.xtc.
HELP
            exit 0 ;;
        *) echo "ERROR: Unknown argument: $1. Use --help." >&2; exit 1 ;;
    esac
done

[[ "$DIRECTION" =~ ^(forward|reverse|both)$ ]] \
    || { echo "ERROR: --direction must be: forward, reverse, or both" >&2; exit 1; }

if [[ "$DIRECTION" == "reverse" || "$DIRECTION" == "both" ]]; then
    [[ -z "$REV_GRO" ]] && echo -e "${YELLOW}NOTE: --rev-gro not provided; will auto-extract from forward PROD lambda31.${NC}"
    [[ -z "$REV_GRO_UNBOUND" ]] && echo -e "${YELLOW}NOTE: --rev-gro-unbound not provided; will auto-extract from forward PROD lambda31.${NC}"
fi

#===============================================================
# DERIVED PATHS
#===============================================================
INPUT_DIR="${FEP_ROOT}/inputs/${DESIGN_NAME}"
WORK_DIR="${FEP_ROOT}/outputs/${DESIGN_NAME}"
MDP_DIR="${FEP_ROOT}/mdps"

cd "$FEP_ROOT" || { echo "ERROR: Cannot access FEP_ROOT: $FEP_ROOT"; exit 1; }

#===============================================================
# UTILITY FUNCTIONS
#===============================================================
fatal() { echo -e "${RED}FATAL (code $1): $2${NC}" >&2; exit "$1"; }

addTERpdb() {
    local _pdbFileName="${1%.pdb}"
    local _O_terminal_lines _flag _increment=0 _temp="addTERpdb_temp"
    cp "${_pdbFileName}.pdb" "./${_temp}.pdb"
    _O_terminal_lines="$(grep -n ' OC2 ' "${_temp}.pdb" | cut -d: -f1)"
    for lineOC2 in $_O_terminal_lines; do
        ((lineOC2 += _increment))
        local lineTER_check=$(( lineOC2 + 1 ))
        _flag=$(sed -n "${lineTER_check}p" "${_temp}.pdb" | cut -b1-3)
        if [[ "$_flag" != "TER" ]]; then
            ((_increment++))
            sed -i "${lineOC2}a\\TER" "${_temp}.pdb"
        fi
    done
    mv "${_temp}.pdb" "${_pdbFileName}.pdb"
    rm -f ./addTERpdb_temp*.pdb
}

make_posres() {
    local _pdb_name="${1%.pdb}" _posres_residues="$2" _chain="$3" _itp_glob="$4" _which_atoms="${5:-4}"
    local _posres_F="200 200 200"
    local _itp_name; _itp_name=$(ls ${_itp_glob}); _itp_name="${_itp_name%.itp}"

    echo -e "${BLUE}  Creating POSRES: chain=${_chain} residues=${_posres_residues}${NC}"
    [[ -r "${_pdb_name}.pdb" ]] || fatal 1 "make_posres: ${_pdb_name}.pdb not readable"
    grep "^ATOM.\{17\}${_chain} " "${_pdb_name}.pdb" > "temp_posres_ch${_chain}.pdb"
    [[ -s "temp_posres_ch${_chain}.pdb" ]] \
        || fatal 1 "make_posres: No ATOM records for chain ${_chain}"

    echo -e "keep ${_which_atoms}\nr ${_posres_residues}\n0 & 1\nname 2 posres-custom\n\nq\n" | \
        gmx make_ndx -f "temp_posres_ch${_chain}.pdb" -o "index_posres_ch${_chain}.ndx" \
        &> "index_posres_ch${_chain}.out" \
        || fatal 1 "make_ndx failed for chain ${_chain}"

    echo "2" | gmx genrestr \
        -f "temp_posres_ch${_chain}.pdb" \
        -n "index_posres_ch${_chain}.ndx" \
        -fc ${_posres_F} \
        -o "posre_custom_ch${_chain}.itp" \
        &> "genrestr_ch${_chain}.out" \
        || fatal 1 "genrestr failed for chain ${_chain}"

    printf "\n; POSRES_custom -> pdb=%s resid=%s chain=%s\n" "${_pdb_name}" "${_posres_residues}" "${_chain}" >> "${_itp_name}.itp"
    printf "#ifdef POSRES_custom\n#include \"posre_custom_ch%s.itp\"\n#endif\n\n" "${_chain}" >> "${_itp_name}.itp"
    rm -f "temp_posres_ch${_chain}.pdb"
    echo -e "${GREEN}  ✓ POSRES created for chain ${_chain}${NC}"
}

write_minim_mdp() {
    local _lambda_vector _zeros _n
    _lambda_vector=$(printf " %s" "${LAMBDAS[@]}"); _lambda_vector="${_lambda_vector:1}"
    _n="${#LAMBDAS[@]}"
    _zeros=$(python3 -c "print(' '.join(['0.0']*${_n}))")
    cat > minim.mdp << EOF
; Minimal MDP for genion TPR generation (pmx hybrid topology)
integrator = steep
nsteps = 0
emtol = 1000.0
emstep = 0.01
nstcomm = 1
pbc = xyz
cutoff-scheme = Verlet
ns-type = grid
nstlist = 10
rlist = 1.2
coulombtype = PME
rcoulomb = 1.2
pme-order = 4
fourierspacing = 0.125
vdw-type = cut-off
rvdw = 1.2
vdw-modifier = Force-switch
rvdw-switch = 1.0
DispCorr = no
constraints = h-bonds
Tcoupl = no
Pcoupl = no
gen_vel = no
free-energy = yes
sc-coul = yes
sc-alpha = 0.5
sc-power = 1
sc-sigma = 0.3
nstdhdl = 0
calc-lambda-neighbors = -1
init-lambda-state = 0
fep_lambdas = ${_lambda_vector}
restraint-lambdas = ${_zeros}
mass-lambdas = ${_zeros}
temperature-lambdas = ${_zeros}
EOF
}

his_input_string() {
    local _pdb="$1" _num_his _string=""
    _num_his=$(grep -Ec "(HIS     CA)|(CA  HIS)" "${_pdb}" 2>/dev/null || echo 0)
    for ((i=0; i<_num_his; i++)); do _string+="1\n"; done
    printf "%s" "${_string}"
}

normalize_his() {
    local _f="$1"
    for his_type in "HISE" "HISD" "HISH" "HIE " "HID " "HIH " "HSE " "HSD " "HSP "; do
        sed -i "s/${his_type}/HIS /g" "${_f}"
    done
}


#===============================================================
# Print configuration
#===============================================================
echo -e "${YELLOW}=== FEP Setup: ${MUTATION} ===${NC}"
echo -e "${BLUE}FEP_ROOT:    ${FEP_ROOT}${NC}"
echo -e "${BLUE}Design:      ${DESIGN_NAME}${NC}"
echo -e "${BLUE}Input GRO:   ${INPUT_DIR}/${INPUT_STRUCTURE}${NC}"
echo -e "${BLUE}Temperature: ${SIMULATION_TEMP} K${NC}"
echo -e "${GREEN}Direction:   ${DIRECTION}  |  ${#LAMBDAS[@]} lambda windows  |  ${NR_REPLICAS} replicas/state${NC}"

#===============================================================
# Step 1: Load environment
#===============================================================
echo -e "\n${BLUE}[Step 1] Loading environment...${NC}"

# Load GROMACS (serial version for setup; MPI used by SLURM jobs)
GROMACS_LOADED=false
if command -v module &>/dev/null; then
    for mod in gromacs/2024.2 gromacs/2024 gromacs; do
        if module load "$mod" 2>/dev/null && command -v gmx &>/dev/null; then
            echo -e "${GREEN}✓ GROMACS loaded via module: $mod${NC}"; GROMACS_LOADED=true; break
        fi
    done
fi
if [[ "$GROMACS_LOADED" = false ]]; then
    for gmx_path in /public/software/apps/gromacs/2024.2/bin/gmx \
                    /public/software/apps/gromacs/2024/bin/gmx \
                    /usr/local/gromacs/bin/gmx; do
        if [[ -x "$gmx_path" ]]; then
            export PATH="$(dirname "${gmx_path}"):${PATH}"
            echo -e "${GREEN}✓ GROMACS from: $gmx_path${NC}"; GROMACS_LOADED=true; break
        fi
    done
fi
command -v gmx &>/dev/null && GROMACS_LOADED=true
[[ "$GROMACS_LOADED" = false ]] && fatal 1 "GROMACS not found. Set PATH or use 'module load gromacs'."

# Activate conda env
CONDA_BASE=$(conda info --base 2>/dev/null) || fatal 1 "conda not found"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}" || fatal 1 "Failed to activate conda env: ${CONDA_ENV}"
echo -e "${GREEN}✓ conda env: ${CONDA_ENV}${NC}"

# Set GMXLIB for mutff force field
export GMXLIB="$GMXLIB_PATH"
[[ -d "$GMXLIB" ]] || fatal 1 "GMXLIB not found: $GMXLIB"
echo -e "${GREEN}✓ GMXLIB: ${GMXLIB}${NC}"

# Verify tools (pip show only reads metadata - no slow NFS package import)
command -v gmx &>/dev/null || fatal 1 "gmx not found"
# pmx is packaged as pmx_biobb on bioconda; both names checked for robustness
pip show pmx_biobb > /dev/null 2>&1 \
    || pip show pmx > /dev/null 2>&1 \
    || fatal 1 "pmx/pmx_biobb not installed in ${CONDA_ENV}"
echo -e "${GREEN}✓ gmx, pmx all available${NC}"

#===============================================================
# Step 2: Check input files
#===============================================================
echo -e "\n${BLUE}[Step 2] Checking input files...${NC}"
[[ -f "${INPUT_DIR}/${INPUT_STRUCTURE}" ]] \
    || fatal 2 "GRO not found: ${INPUT_DIR}/${INPUT_STRUCTURE}"
if [[ "$DIRECTION" != "reverse" ]]; then
    [[ -f "${INPUT_DIR}/${INPUT_TPR}" ]] \
        || fatal 2 "TPR not found: ${INPUT_DIR}/${INPUT_TPR}"
    [[ -f "${INPUT_DIR}/${INPUT_PDB}" ]] \
        || fatal 2 "PDB not found: ${INPUT_DIR}/${INPUT_PDB}"
fi
echo -e "${GREEN}✓ Input files OK${NC}"

#===============================================================
# Step 3: Create directory structure
#===============================================================
echo -e "\n${BLUE}[Step 3] Creating output directories...${NC}"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR" || fatal 3 "Cannot enter $WORK_DIR"
[[ "$DIRECTION" == "forward" || "$DIRECTION" == "both" ]] && mkdir -p "forward/bound" "forward/unbound"
[[ "$DIRECTION" == "reverse" || "$DIRECTION" == "both" ]] && mkdir -p "reverse/bound" "reverse/unbound"
echo -e "${GREEN}✓ Output dir: ${WORK_DIR}${NC}"

#===============================================================
# FORWARD DIRECTION (Steps 4-7)
# Skipped entirely when --direction reverse
#===============================================================
if [[ "$DIRECTION" != "reverse" ]]; then

#---------------------------------------------------------------
# Step 4: GRO -> PDB (gmx trjconv, chain IDs from TPR) + pdb2gmx pass 1
#---------------------------------------------------------------
echo -e "\n${BLUE}[Step 4] GRO -> PDB + pdb2gmx pass 1...${NC}"
cd "${WORK_DIR}/forward" || fatal 4 "Cannot enter forward/"

# Use System group (0) - "Protein" index in this TPR only covers antibody chains;
# pHLA chains (A, P) belong to separate molecule groups and would be missed.
echo "0" | gmx trjconv \
    -f "${INPUT_DIR}/${INPUT_STRUCTURE}" \
    -s "${INPUT_DIR}/${INPUT_TPR}" \
    -o "bound_forw_structure_raw.pdb" \
    &> trjconv_gro2pdb.out \
    || fatal 4 "gmx trjconv (GRO->PDB) failed. See: ${WORK_DIR}/forward/trjconv_gro2pdb.out"

# Inject chain IDs from system.pdb (INPUT_PDB) by POSITIONAL mapping (atom index order).
# system.pdb is the pdb2gmx OUTPUT (with H, 6504 atoms) — same atom count/order as md_r1.gro.
# Resid-based mapping CANNOT work because each chain restarts resid numbering from 1.
awk '
    NR==FNR {
        if (/^ATOM/) canon[++n] = substr($0,22,1)
        next
    }
    /^ATOM/ {
        rn = substr($0,18,3); gsub(/ /,"",rn)
        if (rn=="SOL"||rn=="HOH"||rn=="WAT"||rn~/^TIP/ \
            ||rn~/^K\+?$/||rn~/^CL-?$/||rn~/^NA\+?$/ \
            ||rn=="MG"||rn=="ZN"||rn=="CA") next
        $0 = substr($0,1,21) canon[++m] substr($0,23)
        print
    }
    END {
        if (n != m) {
            printf "ERROR: atom count mismatch: system.pdb=%d  gro_protein=%d\n"\
                   "  Ensure INPUT_PDB is the pdb2gmx output (with H) that produced INPUT_STRUCTURE.\n", \
                   n, m > "/dev/stderr"
            exit 1
        }
        print "END"
    }
' "${INPUT_DIR}/${INPUT_PDB}" "bound_forw_structure_raw.pdb" \
    > "bound_forw_structure.pdb" \
    || fatal 4 "Chain ID injection failed (see atom count error above)"
rm -f "bound_forw_structure_raw.pdb"

normalize_his "bound_forw_structure.pdb"

HIS_INPUT=$(his_input_string "bound_forw_structure.pdb")
echo -e "${HIS_INPUT}" | gmx pdb2gmx \
    -f  "bound_forw_structure.pdb" \
    -o  "bound_forw_structure_ff.pdb" \
    -ignh \
    -ff "${FORCE_FIELD}" \
    -water "${WATER_MODEL}" \
    -his \
    &> pdb2gmx_pass1.out \
    || fatal 4 "pdb2gmx pass 1 failed. See: ${WORK_DIR}/forward/pdb2gmx_pass1.out"

normalize_his "bound_forw_structure_ff.pdb"
echo -e "${GREEN}✓ Step 4 done: bound_forw_structure_ff.pdb${NC}"

#---------------------------------------------------------------
# Step 5: BOUND state setup
#---------------------------------------------------------------
echo -e "\n${BLUE}[Step 5] BOUND state setup (Ab+pHLA)...${NC}"
cd "${WORK_DIR}/forward/bound" || fatal 5 "Cannot enter forward/bound"

BASE="structure_forw_mutated"

# 5a: pmx mutate - introduce hybrid HIS(X)->ARG residue
echo -e "${BLUE}  [5a] pmx mutate: chain${PEPTIDE_CHAIN} res${PEPTIDE_RESID} -> ${MUTATION_TO}${NC}"
echo "${PEPTIDE_CHAIN} ${PEPTIDE_RESID} ${MUTATION_TO}" > mutations_Forw.txt
pmx mutate \
    -f ../"bound_forw_structure_ff.pdb" \
    -o "${BASE}.pdb" \
    -ff "${FORCE_FIELD}" \
    --script ./mutations_Forw.txt \
    --keep_resid \
    > pmx_mutate.out 2>&1 \
    || { tail -20 pmx_mutate.out >&2
         fatal 5 "pmx mutate failed. See: ${WORK_DIR}/forward/bound/pmx_mutate.out"; }
normalize_his "${BASE}.pdb"
echo -e "${GREEN}  ✓ pmx mutate done${NC}"

# 5b: Add TER records at chain C-termini (required by pdb2gmx for multi-chain PDB)
addTERpdb "${BASE}.pdb"
echo -e "${GREEN}  ✓ TER records added${NC}"

# 5c: pdb2gmx pass 2 - build topology for hybrid structure
echo -e "${BLUE}  [5c] pdb2gmx pass 2: building hybrid topology...${NC}"
HIS_INPUT=$(his_input_string "${BASE}.pdb")
echo -e "${HIS_INPUT}" | gmx pdb2gmx \
    -f "${BASE}.pdb" \
    -o "${BASE}_ff.pdb" \
    -ff "${FORCE_FIELD}" \
    -water "${WATER_MODEL}" \
    -his \
    &> pdb2gmx_pass2.out \
    || fatal 5 "pdb2gmx pass 2 failed. See: ${WORK_DIR}/forward/bound/pdb2gmx_pass2.out"
normalize_his "${BASE}_ff.pdb"
echo -e "${GREEN}  ✓ pdb2gmx pass 2 done${NC}"

# 5d: editconf - define simulation box (triclinic, 1.2 nm margin, principal-axis aligned)
echo -e "${BLUE}  [5d] editconf: setting box...${NC}"
echo "1" | gmx editconf \
    -f "${BASE}_ff.pdb" \
    -o "${BASE}_ff_box.pdb" \
    -d 1.2 -bt triclinic -c -princ \
    &> editconf.out \
    || fatal 5 "editconf failed. See: ${WORK_DIR}/forward/bound/editconf.out"
echo -e "${GREEN}  ✓ editconf done${NC}"

# 5e: pmx gentop - generate dual (hybrid) topology
echo -e "${BLUE}  [5e] pmx gentop: generating hybrid topology...${NC}"
pmx gentop -p topol.top -o topol_hybrid.top -ff "${FORCE_FIELD}" \
    > pmx_gentop.out 2>&1 \
    || { tail -20 pmx_gentop.out >&2
         fatal 5 "pmx gentop failed. See: ${WORK_DIR}/forward/bound/pmx_gentop.out"; }
echo -e "${GREEN}  ✓ topol_hybrid.top created${NC}"

# 5f: Position restraints for HLA alpha chain
echo -e "${BLUE}  [5f] POSRES for chain(s): ${POSRE_CHAINS}...${NC}"
for chain in ${POSRE_CHAINS}; do
    make_posres "${BASE}_ff_box.pdb" "${POSRE_RESIDUES}" "${chain}" \
        "pmx_topol_Protein_chain_${chain}.itp"
done

# 5g: Copy topology files to unbound directory (shared before chain removal)
echo -e "${BLUE}  [5g] Copying files to forward/unbound/...${NC}"
cp ./pmx_topol_Protein_chain_*.itp \
   ./posre_custom_ch*.itp \
   ./topol_hybrid.top \
   "./${BASE}_ff_box.pdb" \
   "${WORK_DIR}/forward/unbound/" \
    || fatal 5 "Failed to copy files to unbound/"
rm -f ./*\#

# 5h: Solvate + genion + index for BOUND
echo -e "${BLUE}  [5h] Solvating BOUND...${NC}"
gmx solvate \
    -cp "${BASE}_ff_box.pdb" -cs spc216.gro \
    -o "${BASE}_ff_water.gro" -p topol_hybrid.top \
    &> solvate.out \
    || fatal 5 "solvate failed. See: ${WORK_DIR}/forward/bound/solvate.out"

write_minim_mdp
gmx grompp \
    -f minim.mdp -c "${BASE}_ff_water.gro" -p topol_hybrid.top \
    -o "${BASE}_ff_ions.tpr" -maxwarn 2 \
    &> grompp_genion.out \
    || fatal 5 "grompp (genion) failed"

echo "SOL" | gmx genion \
    -s "${BASE}_ff_ions.tpr" -o "${BASE}_ff_ions.gro" -p topol_hybrid.top \
    -nname "${ION_NEG}" -pname "${ION_POS}" -neutral -conc "${ION_CONC}" \
    &> genion.out \
    || fatal 5 "genion failed. See: ${WORK_DIR}/forward/bound/genion.out"

echo -e "keep 0\nr SOL ${ION_POS} ${ION_NEG}\n0 & ! 1\nname 2 SOLU\nname 1 SOLV\n\nq\n" | \
    gmx make_ndx -f "${BASE}_ff_ions.gro" -o index.ndx &> make_ndx.out \
    || fatal 5 "make_ndx (bound) failed"
echo -e "${GREEN}✓ BOUND setup complete${NC}"

#---------------------------------------------------------------
# Step 6: UNBOUND state setup (remove Ab chains, re-solvate)
#---------------------------------------------------------------
echo -e "\n${BLUE}[Step 6] UNBOUND state setup (pHLA only)...${NC}"
cd "${WORK_DIR}/forward/unbound" || fatal 6 "Cannot enter forward/unbound"

UNB_BASE="pHLA_forw_mutated"
mv "${BASE}_ff_box.pdb" "${UNB_BASE}_ff_box.pdb" || fatal 6 "Cannot rename structure in unbound/"

# Remove antibody chains from PDB + topology + itp files
for chain in ${CHAINS_TO_REMOVE}; do
    sed -i "/^ATOM.\{17\}${chain} /d; /^HETATM.\{15\}${chain} /d" "${UNB_BASE}_ff_box.pdb"
    sed -i "/Protein_chain_${chain}/d" topol_hybrid.top
    rm -f "pmx_topol_Protein_chain_${chain}.itp" "posre_custom_ch${chain}.itp"
    echo -e "${GREEN}  ✓ Removed chain ${chain}${NC}"
done
rm -f ./*\#

gmx solvate \
    -cp "${UNB_BASE}_ff_box.pdb" -cs spc216.gro \
    -o "${UNB_BASE}_ff_water.gro" -p topol_hybrid.top \
    &> solvate.out \
    || fatal 6 "solvate (unbound) failed. See: ${WORK_DIR}/forward/unbound/solvate.out"

write_minim_mdp
gmx grompp \
    -f minim.mdp -c "${UNB_BASE}_ff_water.gro" -p topol_hybrid.top \
    -o "${UNB_BASE}_ff_ions.tpr" -maxwarn 2 \
    &> grompp_genion.out \
    || fatal 6 "grompp (genion, unbound) failed"

echo "SOL" | gmx genion \
    -s "${UNB_BASE}_ff_ions.tpr" -o "${UNB_BASE}_ff_ions.gro" -p topol_hybrid.top \
    -nname "${ION_NEG}" -pname "${ION_POS}" -neutral -conc "${ION_CONC}" \
    &> genion.out \
    || fatal 6 "genion (unbound) failed. See: ${WORK_DIR}/forward/unbound/genion.out"

echo -e "keep 0\nr SOL ${ION_POS} ${ION_NEG}\n0 & ! 1\nname 2 SOLU\nname 1 SOLV\n\nq\n" | \
    gmx make_ndx -f "${UNB_BASE}_ff_ions.gro" -o index.ndx &> make_ndx.out \
    || fatal 6 "make_ndx (unbound) failed"
echo -e "${GREEN}✓ UNBOUND setup complete${NC}"

#---------------------------------------------------------------
# Step 7: Generate NR_REPLICAS independent replicas
#---------------------------------------------------------------
echo -e "\n${BLUE}[Step 7] Generating ${NR_REPLICAS} replicas...${NC}"
cd "${WORK_DIR}" || fatal 7 "Cannot enter $WORK_DIR"

mv "forward/bound"   "forward/bound_R1"   || fatal 7 "Cannot rename forward/bound -> bound_R1"
mv "forward/unbound" "forward/unbound_R1" || fatal 7 "Cannot rename forward/unbound -> unbound_R1"

for replica in $(seq 2 "${NR_REPLICAS}"); do
    cp -r "forward/bound_R1"   "forward/bound_R${replica}"   || fatal 7 "Copy bound_R${replica} failed"
    cp -r "forward/unbound_R1" "forward/unbound_R${replica}" || fatal 7 "Copy unbound_R${replica} failed"
    echo -e "${GREEN}  ✓ Replica ${replica} created${NC}"
done
echo -e "${GREEN}✓ Replicas R1..R${NR_REPLICAS} ready${NC}"

fi  # end forward section

#===============================================================
# Step 8: Reverse direction setup
# Reuses forward hybrid topology; only starting GRO differs (ARG-end)
#===============================================================
if [[ "$DIRECTION" == "reverse" || "$DIRECTION" == "both" ]]; then
    echo -e "\n${BLUE}[Step 8] Setting up REVERSE directories...${NC}"
    cd "$WORK_DIR" || fatal 8 "Cannot enter $WORK_DIR"

    [[ -d "forward/bound_R1" ]] \
        || fatal 8 "forward/bound_R1 not found. Run --direction forward first."

    # Auto-extract reverse starting structures from the last lambda window of forward PROD
    # if --rev-gro / --rev-gro-unbound were not supplied.
    N_LAMBDA="${#LAMBDAS[@]}"
    LAST_WIN=$(( N_LAMBDA - 1 ))  # 31 for 32-window schedule

    if [[ -z "$REV_GRO" ]]; then
        FWD_XTC="${WORK_DIR}/forward/bound_R1/lambda${LAST_WIN}/PROD/prod.xtc"
        FWD_TPR="${WORK_DIR}/forward/bound_R1/lambda${LAST_WIN}/PROD/prod.tpr"
        if [[ -f "$FWD_XTC" && -f "$FWD_TPR" ]]; then
            echo -e "${BLUE}  [8-extract] BOUND: extracting last frame from lambda${LAST_WIN}...${NC}"
            REV_GRO="${WORK_DIR}/rev_start_bound.gro"
            echo "0" | gmx trjconv \
                -f "$FWD_XTC" -s "$FWD_TPR" \
                -o "$REV_GRO" -dump -1 \
                &> "${WORK_DIR}/rev_extract_bound.out" \
                || fatal 8 "Auto-extract BOUND GRO failed. See: ${WORK_DIR}/rev_extract_bound.out"
            echo -e "${GREEN}  ✓ rev_start_bound.gro extracted${NC}"
        else
            echo -e "${YELLOW}  WARNING: forward PROD data not found; cannot auto-extract BOUND GRO.${NC}"
            echo -e "${YELLOW}    Expected: ${FWD_XTC}${NC}"
            echo -e "${YELLOW}    Provide via --rev-gro or run forward PROD first.${NC}"
        fi
    fi

    if [[ -z "$REV_GRO_UNBOUND" ]]; then
        FWD_XTC="${WORK_DIR}/forward/unbound_R1/lambda${LAST_WIN}/PROD/prod.xtc"
        FWD_TPR="${WORK_DIR}/forward/unbound_R1/lambda${LAST_WIN}/PROD/prod.tpr"
        if [[ -f "$FWD_XTC" && -f "$FWD_TPR" ]]; then
            echo -e "${BLUE}  [8-extract] UNBOUND: extracting last frame from lambda${LAST_WIN}...${NC}"
            REV_GRO_UNBOUND="${WORK_DIR}/rev_start_unbound.gro"
            echo "0" | gmx trjconv \
                -f "$FWD_XTC" -s "$FWD_TPR" \
                -o "$REV_GRO_UNBOUND" -dump -1 \
                &> "${WORK_DIR}/rev_extract_unbound.out" \
                || fatal 8 "Auto-extract UNBOUND GRO failed. See: ${WORK_DIR}/rev_extract_unbound.out"
            echo -e "${GREEN}  ✓ rev_start_unbound.gro extracted${NC}"
        else
            echo -e "${YELLOW}  WARNING: forward PROD data not found; cannot auto-extract UNBOUND GRO.${NC}"
            echo -e "${YELLOW}    Expected: ${FWD_XTC}${NC}"
            echo -e "${YELLOW}    Provide via --rev-gro-unbound or run forward PROD first.${NC}"
        fi
    fi

    for replica in $(seq 1 "${NR_REPLICAS}"); do
        for state in bound unbound; do
            src="forward/${state}_R${replica}"
            dst="reverse/${state}_R${replica}"
            mkdir -p "$dst"

            cp "${src}/topol_hybrid.top" "$dst/"
            cp "${src}"/*.itp "$dst/" 2>/dev/null || true
            cp "${src}/index.ndx" "$dst/"

            if [[ "$state" == "bound" ]]; then
                gro_src="$REV_GRO"; gro_flag="--rev-gro"
            else
                gro_src="$REV_GRO_UNBOUND"; gro_flag="--rev-gro-unbound"
            fi

            if [[ -n "$gro_src" && -f "$gro_src" ]]; then
                fwd_gro=$(ls "${src}"/*ions.gro 2>/dev/null | head -1)
                gro_basename=$(basename "${fwd_gro:-structure_forw_mutated_ff_ions.gro}")
                cp "$gro_src" "${dst}/${gro_basename}"
                echo -e "${GREEN}  ✓ reverse/${state}_R${replica}: ions.gro placed${NC}"
            else
                echo -e "${YELLOW}  NOTE: reverse/${state}_R${replica} has no ions.gro yet (provide via ${gro_flag}).${NC}"
            fi
        done
    done
    echo -e "${GREEN}✓ Reverse directories ready${NC}"
fi

#===============================================================
# Summary
#===============================================================
echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}  FEP SETUP COMPLETE${NC}"
echo -e "${GREEN}  Design:    ${DESIGN_NAME}${NC}"
echo -e "${GREEN}  Mutation:  ${MUTATION}${NC}"
echo -e "${GREEN}  Direction: ${DIRECTION}${NC}"
echo -e "${GREEN}  Output:    ${WORK_DIR}${NC}"
[[ "$DIRECTION" != "reverse" ]] && \
    echo -e "${GREEN}  Forward:   bound_R1..R${NR_REPLICAS} + unbound_R1..R${NR_REPLICAS}${NC}"
[[ "$DIRECTION" == "reverse" || "$DIRECTION" == "both" ]] && \
    echo -e "${GREEN}  Reverse:   ${WORK_DIR}/reverse/${NC}"
echo -e "${GREEN}  Next:      bash scripts/submit_fep.sh${NC}"
echo -e "${GREEN}========================================${NC}"
