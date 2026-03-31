#!/bin/bash
#SBATCH --job-name=pHLA_setup_charmm
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node20,node21,node24,node26
#SBATCH --dependency=singleton
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=3-00:00:00
#SBATCH --output=/public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/setup/slurm_setup_%j.out
#SBATCH --error=/public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/setup/slurm_setup_%j.err

# ============================================================================
# Conda Environment Activation
# ============================================================================
CONDA_BASE=/public/home/xuziyi/miniconda
CONDA_ENV=pHLA_MD_env

source $CONDA_BASE/etc/profile.d/conda.sh
conda activate $CONDA_ENV

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate conda environment: $CONDA_ENV"
    echo "Please run setup_conda_env.sh first to create the environment"
    exit 1
fi

echo "Conda environment activated: $(conda info --envs | grep '*' | awk '{print $1}')"
echo "Python version: $(python3 --version)"
echo ""

# ============================================================================
# Load HPC Modules
# ============================================================================
module load gromacs/2024.2
module load cuda

# Set custom CHARMM36m force field path
export GMXLIB=/public/home/xuziyi/FEP/force_fields/mutff

# GPU optimization flags
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_MAXCONSTRWARN=-1

# GROMACS path
GMX="gmx"

# MDP files directory
MDP_DIR="/public/home/xuziyi/MD/mdp_files_charmm"

# Python scripts path
SCRIPT_DIR="/public/home/xuziyi/MD/scripts_charmm"

# Input PDB file
PDB_FILE="/public/home/xuziyi/MD/data/design_373_dldesign_13_best.pdb"

echo "======================================"
echo "Start processing: design_373_dldesign_13_best (CHARMM36m)"
echo "======================================"

# ========================================
# Step 1: System Setup
# ========================================
echo "Step 1: System preparation (topology, solvation, ions, minimization)..."
cd /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/setup

python3 $SCRIPT_DIR/setup_system_charmm.py \
    -i "$PDB_FILE" \
    -w /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/setup \
    -ions "$MDP_DIR/ions.mdp" \
    -em "$MDP_DIR/minim.mdp" \
    -gmx "$GMX" \
    -ff charmm36m-mut \
    -water tip3p \
    -box triclinic \
    -d 1.2 \
    -conc 0.10

if [ $? -ne 0 ]; then
    echo "ERROR: System setup failed!"
    exit 1
fi

# ========================================
# Step 2: Equilibration (NVT + NPT)
# ========================================
echo "Step 2: Running NVT and NPT equilibration..."

python3 $SCRIPT_DIR/run_equilibration_charmm.py \
    -w /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/setup \
    -nvt "$MDP_DIR/nvt.mdp" \
    -npt "$MDP_DIR/npt.mdp" \
    -gmx "$GMX"

if [ $? -ne 0 ]; then
    echo "ERROR: Equilibration failed!"
    exit 1
fi

echo "======================================"
echo "Setup and equilibration completed!"
echo "======================================"
echo "Now submit production jobs:"
echo "  sbatch /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/submit_charmm_design_373_dldesign_13_best_r1.slurm"
echo "  sbatch /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/submit_charmm_design_373_dldesign_13_best_r2.slurm"
echo "  sbatch /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/submit_charmm_design_373_dldesign_13_best_r3.slurm"
