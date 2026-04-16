#!/bin/bash
# ============================================================================
# setup_mmpbsa_env.sh — 创建 gmxMMPBSA conda 环境
#
# HPC 集群 (tc6000) 无直接外网，所有包通过 ZJU 镜像安装：
#   - conda-forge (由 ~/.condarc 中 ZJU 镜像覆盖)
#   - pip (浙大 PyPI 镜像)
#
# 运行前确认 ~/.condarc 已配置：
#   channels:
#     - defaults
#   default_channels:
#     - https://mirrors.zju.edu.cn/anaconda/pkgs/main
#     - https://mirrors.zju.edu.cn/anaconda/pkgs/r
#   custom_channels:
#     conda-forge: https://mirrors.zju.edu.cn/anaconda/cloud
#
# 用法:
#   bash setup_mmpbsa_env.sh
# ============================================================================

set -euo pipefail

ENV_NAME="gmxMMPBSA"
PYTHON_VER="3.10"
PIP_MIRROR="https://mirrors.zju.edu.cn/pypi/web/simple"
CONDA_BASE="${HOME}/miniconda"

echo "========================================"
echo "  Creating conda environment: ${ENV_NAME}"
echo "========================================"

# ── 清除索引缓存（避免旧频道残留）──────────────────────────────────────────
conda clean -i -y

# ── 创建环境 ─────────────────────────────────────────────────────────────────
# ambertools + mpi4py 均来自 conda-forge 预编译二进制
conda create -y -n "${ENV_NAME}" \
    -c conda-forge --override-channels \
    python="${PYTHON_VER}" \
    ambertools \
    compilers \
    "mpi4py>=3.1"

echo ""
echo "Activating ${ENV_NAME}..."
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

echo "Python: $(python3 --version)"
echo "AmberTools: $(python3 -c 'import pytraj; print(pytraj.__version__)' 2>/dev/null || echo 'N/A')"

# ── pip 安装 gmx_MMPBSA ────────────────────────────────────────────────────
# --no-deps 避免 pip 重复安装 ambertools/mpi4py（已有 conda 版）
echo ""
echo "Installing gmx_MMPBSA via pip (no-deps)..."
pip install -i "${PIP_MIRROR}" --no-deps gmx_MMPBSA

# 安装 gmx_MMPBSA 的 Python 依赖（排除 ambertools/mpi4py）
echo "Installing gmx_MMPBSA Python dependencies..."
pip install -i "${PIP_MIRROR}" \
    pandas scipy seaborn tqdm

# ── 验证 ─────────────────────────────────────────────────────────────────────
echo ""
echo "========================================"
echo "  Verification"
echo "========================================"
echo "gmx_MMPBSA version:"
gmx_MMPBSA --version 2>&1 | head -3 || true
echo ""
echo "mpi4py:"
python3 -c "import mpi4py; print(f'  mpi4py {mpi4py.__version__}')"
echo ""

echo "========================================"
echo "  Environment ${ENV_NAME} is ready!"
echo "  Activate with: conda activate ${ENV_NAME}"
echo "========================================"
#!/bin/bash
# ==============================================================================
# Setup gmx_MMPBSA Conda Environment
# ==============================================================================
# Creates a conda environment named "gmxMMPBSA" with:
#   - AmberTools (via conda-forge — no ambermd channel needed)
#   - gmx_MMPBSA (via pip / PyPI mirror)
#   - Python analysis dependencies
#
# HPC Notes:
#   - Requires ~/.condarc configured with ZJU (Zhejiang University) mirrors.
#   - Does NOT use the ambermd channel (no external internet required).
#   - Run once on the login node or an interactive node.
#   - Usage: bash setup_mmpbsa_env.sh
#
# Required ~/.condarc (ZJU mirror):
#   channels:
#     - defaults
#   show_channel_urls: true
#   default_channels:
#     - https://mirrors.zju.edu.cn/anaconda/pkgs/main
#     - https://mirrors.zju.edu.cn/anaconda/pkgs/r
#     - https://mirrors.zju.edu.cn/anaconda/pkgs/msys2
#   custom_channels:
#     conda-forge: https://mirrors.zju.edu.cn/anaconda/cloud
#     bioconda: https://mirrors.zju.edu.cn/anaconda/cloud
#
# If conda-forge packages are still not found, run:
#   conda clean -i   (manually clear index cache)
# ==============================================================================

set -euo pipefail

CONDA_BASE="/public/home/xuziyi/miniconda"
ENV_NAME="gmxMMPBSA"
LOG_FILE="setup_mmpbsa_env.log"

# pip mirror (ZJU PyPI mirror for gmx_MMPBSA and analysis packages)
PIP_INDEX="https://mirrors.zju.edu.cn/pypi/web/simple"
PIP_TRUSTED="mirrors.zju.edu.cn"
PIP_FLAGS="-i ${PIP_INDEX} --trusted-host ${PIP_TRUSTED}"

echo "========================================"
echo " gmx_MMPBSA Environment Setup"
echo " $(date)"
echo "========================================"
echo ""

# ── Activate conda ──────────────────────────────────────────────────────────
echo "[1/6] Initializing conda..."
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# ── Print current channel configuration ─────────────────────────────────────
echo "    Conda channels (from ~/.condarc):"
conda config --show channels 2>/dev/null | head -10 || true
echo ""

# ── Clear stale index cache (ensures ZJU mirror index is used) ───────────────
echo "[2/6] Clearing conda index cache (conda clean -i)..."
conda clean -i -y 2>&1 | tee -a "${LOG_FILE}"
echo "    Index cache cleared."

# ── Check if environment already exists ─────────────────────────────────────
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "    WARNING: Environment '${ENV_NAME}' already exists."
    echo "    To reinstall:  conda env remove -n ${ENV_NAME}  then re-run this script."
    echo "    To update:     conda activate ${ENV_NAME} && pip install -U gmx_MMPBSA ${PIP_FLAGS}"
    exit 0
fi

# ── Create environment ────────────────────────────────────────────────────────
# AmberTools is on conda-forge → served by ZJU mirror (no ambermd channel needed).
# --override-channels ensures only the channels in ~/.condarc are used.
echo "[3/6] Creating conda environment '${ENV_NAME}' (AmberTools from conda-forge)..."
echo "    Mirror: https://mirrors.zju.edu.cn/anaconda/cloud/conda-forge"
echo "    This may take 10-20 minutes..."
conda create -n "${ENV_NAME}" \
    -c conda-forge \
    --override-channels \
    python=3.10 \
    ambertools \
    compilers \
    "mpi4py>=3.1" \
    -y \
    2>&1 | tee "${LOG_FILE}"

echo "    Environment created."

# ── Activate and install gmx_MMPBSA via pip (ZJU PyPI mirror) ───────────────
# mpi4py is already installed by conda above (pre-built with MPICH).
# --no-deps skips re-installing dependencies that are already satisfied,
# then we explicitly install remaining pip-only deps.
echo "[4/6] Installing gmx_MMPBSA via pip (ZJU mirror)..."
conda activate "${ENV_NAME}"

# Install gmx_MMPBSA without pulling in mpi4py from pip (conda version preferred)
pip install gmx_MMPBSA --no-deps ${PIP_FLAGS} 2>&1 | tee -a "${LOG_FILE}"

# Install remaining gmx_MMPBSA runtime dependencies (excluding mpi4py and parmed)
pip install "numpy==1.26.4" "pandas==1.5.3" "matplotlib==3.7.3" \
    "scipy==1.14.1" seaborn tqdm ${PIP_FLAGS} 2>&1 | tee -a "${LOG_FILE}"

# ── Install additional Python analysis packages ──────────────────────────────
echo "[5/6] Installing extra analysis tools..."
pip install scikit-learn ${PIP_FLAGS} 2>&1 | tee -a "${LOG_FILE}"

# ── Verify installation ──────────────────────────────────────────────────────
echo "[6/6] Verifying installation..."
echo ""
echo "  Python    : $(python3 --version)"
echo "  gmx_MMPBSA: $(gmx_MMPBSA --version 2>&1 | head -1)"
echo "  tleap     : $(which tleap 2>/dev/null || echo 'NOT FOUND')"
echo ""

# ── Load GROMACS to verify compatibility ────────────────────────────────────
echo "  Testing GROMACS module..."
module load gromacs/2024.2 2>/dev/null \
    && echo "  GROMACS loaded successfully." \
    || echo "  WARNING: 'module load gromacs/2024.2' not available on this node."

echo ""
echo "========================================"
echo " Installation complete!"
echo " Log saved to: ${LOG_FILE}"
echo ""
echo " To activate:"
echo "   source ${CONDA_BASE}/etc/profile.d/conda.sh"
echo "   conda activate ${ENV_NAME}"
echo ""
echo " Quick test:"
echo "   gmx_MMPBSA --help"
echo "========================================"
