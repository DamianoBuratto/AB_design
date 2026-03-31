#!/bin/bash
#SBATCH --job-name=RFantibody_test
#SBATCH --output=logs/rfantibody_test_%j.out
#SBATCH --error=logs/rfantibody_test_%j.err
#SBATCH --partition=debug
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node10
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00

source /public/home/xuziyi/miniconda/etc/profile.d/conda.sh
conda activate RFantibody_HPC

# Load required modules
source /public/software/profile.d/compiler_intel-2021.3.0.sh
source /public/software/profile.d/cuda12.2-env.sh
source /public/software/profile.d/fftw3_3.3.10_float-sse2-avx2.sh

# Setup glibc environment
GLIBC_ROOT=/public/software/mathlib/glibc/2.34
GLIBC_LOADER=${GLIBC_ROOT}/lib/ld-linux-x86-64.so.2
GLIBC_LIB_PATH=${GLIBC_ROOT}/lib:${GLIBC_ROOT}/lib64

CONDA_LIB_DIR="${CONDA_PREFIX}/lib"
if [ ! -d "${CONDA_LIB_DIR}" ]; then
	echo "Expected conda lib directory not found: ${CONDA_LIB_DIR}" >&2
	exit 1
fi

GLIBC_RUNTIME_PATH="${GLIBC_LIB_PATH}:${CONDA_LIB_DIR}"
if [ -n "${LD_LIBRARY_PATH:-}" ]; then
	GLIBC_RUNTIME_PATH="${GLIBC_RUNTIME_PATH}:${LD_LIBRARY_PATH}"
fi

export GLIBC_ROOT
export GLIBC_LOADER
export RFANTIBODY_GLIBC_LIBRARY_PATH="${GLIBC_RUNTIME_PATH}"
export RFANTIBODY_PYTHON_BIN="${CONDA_PREFIX}/bin/python"
if [ ! -x "${RFANTIBODY_PYTHON_BIN}" ]; then
	echo "Python interpreter not found at ${RFANTIBODY_PYTHON_BIN}" >&2
	exit 1
fi

if [ ! -x "${GLIBC_LOADER}" ]; then
	echo "Custom glibc loader not found at ${GLIBC_LOADER}" >&2
	exit 1
fi

# Helper function to run commands with custom glibc
run_with_glibc() {
	if [ "$#" -eq 0 ]; then
		echo "run_with_glibc: missing command" >&2
		return 1
	fi
	if [ -z "${GLIBC_LOADER:-}" ]; then
		echo "run_with_glibc: GLIBC_LOADER is not set" >&2
		return 1
	fi
	if [ ! -x "${GLIBC_LOADER}" ]; then
		echo "run_with_glibc: loader not executable: ${GLIBC_LOADER}" >&2
		return 1
	fi
	local cmd="$1"
	shift
	if [[ "${cmd}" != /* ]]; then
		cmd=$(command -v "${cmd}")
	fi
	if [ -z "${cmd}" ]; then
		echo "run_with_glibc: command not found" >&2
		return 1
	fi
	"${GLIBC_LOADER}" --library-path "${RFANTIBODY_GLIBC_LIBRARY_PATH}" "${cmd}" "$@"
}

export -f run_with_glibc

# Setup environment variables
export PYTHONPATH=$HOME/RFantibody/src:$HOME/RFantibody/src/rfantibody/rfdiffusion:$HOME/RFantibody/src/rfantibody:$HOME/RFantibody/include/SE3Transformer:$PYTHONPATH
export PATH=$HOME/RFantibody/include/USalign:$PATH
export ICECREAM_COLORS=never
export HYDRA_FULL_ERROR=1

cd $HOME/RFantibody

# Create output directories
mkdir -p scripts/examples/example_outputs
mkdir -p logs
mkdir -p scripts/config
ln -sfn "$HOME/RFantibody/src/rfantibody/rfdiffusion/config/inference" scripts/config/inference

echo "Test job started: $(date)"

# Run RFdiffusion inference using custom glibc-wrapped python
run_with_glibc "${RFANTIBODY_PYTHON_BIN}" scripts/rfdiffusion_inference.py \
    --config-name antibody \
    antibody.target_pdb=$HOME/RFantibody/scripts/examples/example_inputs/rsv_site3.pdb \
    antibody.framework_pdb=$HOME/RFantibody/scripts/examples/example_inputs/hu-4D5-8_Fv.pdb \
    inference.ckpt_override_path=$HOME/RFantibody/weights/RFdiffusion_Ab.pt \
    'ppi.hotspot_res=[T305,T456]' \
    'antibody.design_loops=[L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13]' \
    inference.num_designs=2 \
    inference.final_step=48 \
    diffuser.T=50 \
    inference.deterministic=True \
    inference.output_prefix=$HOME/RFantibody/scripts/examples/example_outputs/ab_des

echo "Test job completed: $(date)"
