#!/bin/bash

# This script is responsible for:
# 1. Checking and preparing input files
# 2. Creating output directory structure
# 3. Standardizing PDB file format
# 4. Checking environment and dependencies

set -e  # Exit script immediately on any error

# ============================================================================
# Global Variable Definitions
# ============================================================================
# Check if running from /home directory (container environment)
if [ "$(pwd)" = "/home" ] && [ -d "/home/weights" ]; then
    PROJECT_DIR="/home"
else
    PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Set Python path to include src directory for proper module imports
# Add both the general src directory and the specific rfdiffusion directory
export PYTHONPATH="${PROJECT_DIR}/src/rfantibody/rfdiffusion:${PROJECT_DIR}/src:${PYTHONPATH}"

INPUT_PDB="${PROJECT_DIR}/data/p53_R175H_cut.pdb"
OUTPUT_DIR="${PROJECT_DIR}/outputs"
FRAMEWORK_PDB="${PROJECT_DIR}/data/6w51_antibody.pdb"
WEIGHTS_DIR="${PROJECT_DIR}/weights"

# ============================================================================
# Function Definitions
# ============================================================================
check_environment() {
    echo "Checking runtime environment..."
    
    # Check Python environment
    if ! command -v python3 >/dev/null 2>&1; then
        echo "Error: Python3 is not available"
        exit 1
    fi
    
    # Check Python packages
    echo "Checking Python dependencies..."
    python3 -c "
import sys
missing_packages = []
critical_deps = {
    'torch': 'PyTorch deep learning framework',
    'numpy': 'NumPy numerical computing',
    'hydra': 'Hydra configuration management',
    'omegaconf': 'OmegaConf configuration processing',
    'pandas': 'Pandas data processing',
    'scipy': 'SciPy scientific computing'
}

for package, description in critical_deps.items():
    try:
        __import__(package)
        print(f'OK {package}')
    except ImportError:
        missing_packages.append(package)
        print(f'MISSING {package} - {description}')

# Check GPU availability
try:
    import torch
    if torch.cuda.is_available():
        gpu_name = torch.cuda.get_device_properties(0).name
        print(f'GPU available: {gpu_name}')
    else:
        print('WARNING: No GPU detected, inference will be slow')
except:
    pass

if missing_packages:
    print(f'Missing critical packages: {missing_packages}')
    print('Install with: pip install ' + ' '.join(missing_packages))
    sys.exit(1)
else:
    print('All critical dependencies satisfied')
"
    
    # Check weight files
    if [ ! -f "${WEIGHTS_DIR}/RFdiffusion_Ab.pt" ]; then
        echo "ERROR: Model weight file missing"
        echo "Run: bash include/download_weights.sh"
        exit 1
    fi
    
    # Check input files
    if [ ! -f "$INPUT_PDB" ]; then
        echo "ERROR: Input PDB file does not exist: $INPUT_PDB"
        exit 1
    fi
    
    if [ ! -f "$FRAMEWORK_PDB" ]; then
        echo "ERROR: Framework PDB file does not exist: $FRAMEWORK_PDB"
        exit 1
    fi
    
    echo "Environment check passed"
}

setup_directories() {
    echo "Creating output directory structure..."
    mkdir -p "$OUTPUT_DIR"/{rfdiffusion,proteinmpnn,rf2}
    echo "Output directories created"
}

prepare_pdb_files() {
    echo "Preparing input files..."
    
    # Create configuration file
    CONFIG_FILE="${PROJECT_DIR}/run_config.txt"
    cat > "$CONFIG_FILE" << EOF
# RFantibody Runtime Configuration
INPUT_PDB=${INPUT_PDB}
FRAMEWORK_PDB=${FRAMEWORK_PDB}
OUTPUT_DIR=${OUTPUT_DIR}
WEIGHTS_DIR=${WEIGHTS_DIR}
PROJECT_DIR=${PROJECT_DIR}
EOF
    
    echo "Configuration saved to $CONFIG_FILE"
    echo "PROJECT_DIR is set to: $PROJECT_DIR"
}

# ============================================================================
# Main Function
# ============================================================================
main() {
    echo "RFantibody Input Preparation"
    echo "============================"
    
    check_environment
    setup_directories
    prepare_pdb_files
    
    echo "============================"
    echo "Preparation completed"
    echo "Run: bash run_antibody_design.sh"
}

# Run main function
main "$@"
