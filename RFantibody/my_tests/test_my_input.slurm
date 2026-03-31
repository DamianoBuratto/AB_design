#!/bin/bash
#SBATCH --job-name=RFantibody_p53
#SBATCH --output=logs/rfantibody_p53_%j.out
#SBATCH --error=logs/rfantibody_p53_%j.err
#SBATCH --partition=debug
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node10
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00

# ============================================================================
# 简化测试脚本 - 直接调用RFdiffusion（仿照example_test.slurm的成功模式）
# ============================================================================

source /public/home/xuziyi/miniconda/etc/profile.d/conda.sh
conda activate RFantibody_HPC

# 加载必需的模块
source /public/software/profile.d/compiler_intel-2021.3.0.sh
source /public/software/profile.d/cuda12.2-env.sh
source /public/software/profile.d/fftw3_3.3.10_float-sse2-avx2.sh

# ============================================================================
# 设置GLIBC环境（与example_test.slurm完全相同）
# ============================================================================
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

# ============================================================================
# 设置环境变量
# ============================================================================
export PYTHONPATH=$HOME/RFantibody/src:$HOME/RFantibody/src/rfantibody/rfdiffusion:$HOME/RFantibody/src/rfantibody:$HOME/RFantibody/include/SE3Transformer:$PYTHONPATH
export PATH=$HOME/RFantibody/include/USalign:$PATH
export ICECREAM_COLORS=never
export HYDRA_FULL_ERROR=1

cd $HOME/RFantibody

# ============================================================================
# 创建输出目录
# ============================================================================
mkdir -p outputs/p53_test
mkdir -p logs
mkdir -p scripts/config
ln -sfn "$HOME/RFantibody/src/rfantibody/rfdiffusion/config/inference" scripts/config/inference

echo "========================================"
echo "测试任务开始: $(date)"
echo "========================================"

# ============================================================================
# 准备输入文件（使用HLT格式文件）
# ============================================================================
echo "步骤1: 检查HLT格式输入文件..."

# 检查必需的HLT格式文件是否存在
if [ ! -f "$HOME/RFantibody/data/6w51_antibody_HLT.pdb" ]; then
    echo "错误: 找不到 6w51_antibody_HLT.pdb"
    echo "请先运行: python3 scripts/util/chothia2HLT.py data/6w51_antibody.pdb --heavy H --light L -o data/6w51_antibody_HLT.pdb"
    exit 1
fi

if [ ! -f "$HOME/RFantibody/data/p53_R175H_cut_HLT.pdb" ]; then
    echo "错误: 找不到 p53_R175H_cut_HLT.pdb"
    exit 1
fi

echo "✓ HLT格式文件已准备就绪"

# ============================================================================
# 运行RFdiffusion（使用HLT格式文件）
# ============================================================================
echo ""
echo "步骤2: 运行RFdiffusion..."
echo "----------------------------------------"
echo "输入文件:"
echo "  - 目标蛋白: $HOME/RFantibody/data/p53_R175H_cut_HLT.pdb (HLT格式,连续编号)"
echo "  - 抗体框架: $HOME/RFantibody/data/6w51_antibody_HLT.pdb (HLT格式,含CDR标注)"
echo "  - 输出目录: $HOME/RFantibody/outputs/p53_test"
echo "----------------------------------------"

# 使用HLT格式文件运行RFdiffusion
# HLT格式要求：
#  - 抗体框架包含H和L链，带REMARK标注CDR loop位置
#  - 目标蛋白使用T链标识
#  - 残基编号连续（H链从1开始，L链接着H链，T链独立编号）
run_with_glibc "${RFANTIBODY_PYTHON_BIN}" scripts/rfdiffusion_inference.py \
    --config-name antibody \
    antibody.target_pdb=$HOME/RFantibody/data/p53_R175H_cut_HLT.pdb \
    antibody.framework_pdb=$HOME/RFantibody/data/6w51_antibody_HLT.pdb \
    inference.ckpt_override_path=$HOME/RFantibody/weights/RFdiffusion_Ab.pt \
    'ppi.hotspot_res=[T7,T8,T9]' \
    'antibody.design_loops=[L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13]' \
    inference.num_designs=10 \
    inference.final_step=48 \
    diffuser.T=50 \
    inference.deterministic=False \
    inference.output_prefix=$HOME/RFantibody/outputs/p53_test/design

echo ""
echo "========================================"
echo "测试任务完成: $(date)"
echo "========================================"

# 检查结果
if [ -d "$HOME/RFantibody/outputs/p53_test" ]; then
    design_count=$(ls $HOME/RFantibody/outputs/p53_test/design_*.pdb 2>/dev/null | wc -l)
    echo "生成的设计数量: $design_count"
    if [ $design_count -gt 0 ]; then
        echo "✓ RFdiffusion成功完成！"
        echo "结果文件位于: $HOME/RFantibody/outputs/p53_test/"
    else
        echo "✗ 未生成任何设计文件，请检查日志"
    fi
fi

# 自动显示最新日志内容
jobid=${SLURM_JOB_ID:-$(ls logs/rfantibody_p53_*.out 2>/dev/null | sed 's/.*_//;s/.out//' | sort -n | tail -1)}
if [ -n "$jobid" ]; then
    echo "========== SLURM STDOUT =========="
    cat logs/rfantibody_p53_${jobid}.out
    echo "========== SLURM STDERR =========="
    cat logs/rfantibody_p53_${jobid}.err
fi