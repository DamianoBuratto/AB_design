#!/bin/bash
#SBATCH --job-name=charmm_design_373_dldesign_13_best_r2
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node20,node21,node24,node26
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=/public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/replica_2/slurm_r2_%j.out
#SBATCH --error=/public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/replica_2/slurm_r2_%j.err

# ============================================================================
# Conda Environment Activation
# ============================================================================
CONDA_BASE=/public/home/xuziyi/miniconda
CONDA_ENV=pHLA_MD_env

source $CONDA_BASE/etc/profile.d/conda.sh
conda activate $CONDA_ENV

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate conda environment: $CONDA_ENV"
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

echo "======================================"
echo "pHLA MD Production - Replica 2"
echo "======================================"
echo "System: design_373_dldesign_13_best"
echo "Replica: 2"
echo "Working directory: /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/replica_2"
echo "Setup directory: /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/setup"
echo ""

# Check if equilibration is complete
if [ ! -f "/public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/setup/npt.gro" ]; then
    echo "ERROR: Equilibration not complete! Missing npt.gro"
    echo "Please run setup script first: sbatch /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/setup/submit_design_373_dldesign_13_best_setup.slurm"
    exit 1
fi

# ========================================
# Production MD
# ========================================
cd /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/replica_2

python3 $SCRIPT_DIR/run_production_charmm.py \
    -w /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/replica_2 \
    -md "$MDP_DIR/production.mdp" \
    -r 2 \
    -setup /public/home/xuziyi/MD/results/production_run_20260315/design_373_dldesign_13_best/setup \
    -time 150 \
    -gmx "$GMX" \
    --no-analysis

if [ $? -ne 0 ]; then
    echo "ERROR: Production MD failed for replica 2!"
    exit 1
fi

echo "======================================"
echo "Replica 2 completed!"
echo "======================================"
