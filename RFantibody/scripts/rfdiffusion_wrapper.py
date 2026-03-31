#!/usr/bin/env python3
"""
RFdiffusion wrapper - Simplified version based on successful example_test.slurm pattern.
Sets up proper environment and directly calls rfdiffusion_inference.py.
"""

import sys
import os
from pathlib import Path

# Get the project root directory
project_root = Path(__file__).parent.parent
rfdiffusion_path = project_root / "src" / "rfantibody" / "rfdiffusion"

# Set up PYTHONPATH similar to example_test.slurm
# Order matters: rfdiffusion dir first for relative imports
sys.path.insert(0, str(rfdiffusion_path))
sys.path.insert(0, str(project_root / "src" / "rfantibody"))
sys.path.insert(0, str(project_root / "src"))

# Change to rfdiffusion directory for config path resolution
# This is critical - example_test.slurm doesn't change dir, but config_path in @hydra.main is relative
os.chdir(str(rfdiffusion_path))

# Now import and run the original inference script
if __name__ == "__main__":
    # The script uses hydra which handles command line args automatically
    # Just need to make it importable and executable
    inference_script = project_root / "scripts" / "rfdiffusion_inference.py"
    
    # Execute the script in the current namespace to preserve command line args
    with open(inference_script) as f:
        exec(f.read())
