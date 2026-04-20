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
