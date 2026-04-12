#!/usr/bin/env python3
"""
pHLA MD Automation - Batch Processing Script (CHARMM36m)
Batch process multiple PDB files with 3 independent replicas each

1. Scan all PDB files in data/ directory
2. Create independent working directory for each PDB (with charmm prefix)
3. Run complete MD workflow using CHARMM36m force field
4. Generate SLURM batch submission scripts
"""

import os
import sys
import argparse
import shutil
from pathlib import Path
from datetime import datetime


class BatchProcessor:
    """Batch processing manager class"""
    
    def __init__(self, data_dir, output_dir, mdp_dir, 
                 gmx_path="gmx", num_replicas=3):
        """
        Args:
            data_dir (str): PDB files directory
            output_dir (str): Output root directory
            mdp_dir (str): MDP files directory
            gmx_path (str): GROMACS executable path
            num_replicas (int): Number of replicas
        """
        self.data_dir = Path(data_dir).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.mdp_dir = Path(mdp_dir).resolve()
        self.gmx_path = gmx_path
        self.num_replicas = num_replicas
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Log file
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.batch_log = self.output_dir / f"batch_process_{timestamp}.log"
        
    def find_pdb_files(self):
        """Scan all PDB files in data/ directory"""
        pdb_files = sorted(self.data_dir.glob("*.pdb"))
        
        print(f"Found {len(pdb_files)} PDB file(s):")
        for pdb in pdb_files:
            print(f"  - {pdb.name}")
        
        return pdb_files
    
    def create_work_structure(self, pdb_file):
        """Create directory structure for MD workflow"""
        pdb_name = pdb_file.stem
        pdb_dir = self.output_dir / pdb_name
        
        # Create setup and replica directories
        setup_dir = pdb_dir / "setup"
        setup_dir.mkdir(parents=True, exist_ok=True)
        
        replica_dirs = []
        for i in range(1, self.num_replicas + 1):
            replica_dir = pdb_dir / f"replica_{i}"
            replica_dir.mkdir(parents=True, exist_ok=True)
            replica_dirs.append(replica_dir)
        
        return pdb_dir, setup_dir, replica_dirs
    
    def generate_setup_script(self, pdb_name, setup_dir):
        """Generate SLURM script for system setup and equilibration"""
        slurm_script = f"""#!/bin/bash
#SBATCH --job-name=pHLA_setup_charmm
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node15,node20,node21,node24,node26
#SBATCH --dependency=singleton
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=3-00:00:00
#SBATCH --output={setup_dir}/slurm_setup_%j.out
#SBATCH --error={setup_dir}/slurm_setup_%j.err

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

echo "Conda environment activated: $(conda info --envs | grep '*' | awk '{{print $1}}')"
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
GMX="{self.gmx_path}"

# MDP files directory
MDP_DIR="{self.mdp_dir}"

# Python scripts path
SCRIPT_DIR="{Path(__file__).parent.resolve()}"

# Input PDB file
PDB_FILE="{self.data_dir / (pdb_name + '.pdb')}"

echo "======================================"
echo "Start processing: {pdb_name} (CHARMM36m)"
echo "======================================"

# ========================================
# Step 1: System Setup
# ========================================
echo "Step 1: System preparation (topology, solvation, ions, minimization)..."
cd {setup_dir}

python3 $SCRIPT_DIR/setup_system_charmm.py \\
    -i "$PDB_FILE" \\
    -w {setup_dir} \\
    -ions "$MDP_DIR/ions.mdp" \\
    -em "$MDP_DIR/minim.mdp" \\
    -gmx "$GMX" \\
    -ff charmm36m-mut \\
    -water tip3p \\
    -box triclinic \\
    -d 1.2 \\
    -conc 0.10

if [ $? -ne 0 ]; then
    echo "ERROR: System setup failed!"
    exit 1
fi

# ========================================
# Step 2: Equilibration (NVT + NPT)
# ========================================
echo "Step 2: Running NVT and NPT equilibration..."

python3 $SCRIPT_DIR/run_equilibration_charmm.py \\
    -w {setup_dir} \\
    -nvt "$MDP_DIR/nvt.mdp" \\
    -npt "$MDP_DIR/npt.mdp" \\
    -gmx "$GMX"

if [ $? -ne 0 ]; then
    echo "ERROR: Equilibration failed!"
    exit 1
fi

echo "======================================"
echo "Setup and equilibration completed!"
echo "======================================"
echo "Now submit production jobs:"
echo "  sbatch {setup_dir.parent}/submit_charmm_{pdb_name}_r1.slurm"
echo "  sbatch {setup_dir.parent}/submit_charmm_{pdb_name}_r2.slurm"
echo "  sbatch {setup_dir.parent}/submit_charmm_{pdb_name}_r3.slurm"
"""
        
        # Write setup SLURM script
        slurm_file = setup_dir / f"submit_charmm_{pdb_name}_setup.slurm"
        with open(slurm_file, 'w') as f:
            f.write(slurm_script)
        slurm_file.chmod(0o755)
        
        print(f"  Generated setup script: {slurm_file}")
        return slurm_file
    
    def generate_replica_script(self, pdb_name, setup_dir, replica_dir, replica_num):
        """Generate SLURM script for a single production replica"""
        slurm_script = f"""#!/bin/bash
#SBATCH --job-name=charmm_{pdb_name}_r{replica_num}
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node15,node20,node21,node24,node26
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output={replica_dir}/slurm_r{replica_num}_%j.out
#SBATCH --error={replica_dir}/slurm_r{replica_num}_%j.err

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

echo "Conda environment activated: $(conda info --envs | grep '*' | awk '{{print $1}}')"
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
GMX="{self.gmx_path}"

# MDP files directory
MDP_DIR="{self.mdp_dir}"

# Python scripts path
SCRIPT_DIR="{Path(__file__).parent.resolve()}"

echo "======================================"
echo "pHLA MD Production - Replica {replica_num}"
echo "======================================"
echo "System: {pdb_name}"
echo "Replica: {replica_num}"
echo "Working directory: {replica_dir}"
echo "Setup directory: {setup_dir}"
echo ""

# Check if equilibration is complete
if [ ! -f "{setup_dir}/npt.gro" ]; then
    echo "ERROR: Equilibration not complete! Missing npt.gro"
    echo "Please run setup script first: sbatch {setup_dir}/submit_{pdb_name}_setup.slurm"
    exit 1
fi

# ========================================
# Production MD
# ========================================
cd {replica_dir}

python3 $SCRIPT_DIR/run_production_charmm.py \\
    -w {replica_dir} \\
    -md "$MDP_DIR/production.mdp" \\
    -r {replica_num} \\
    -setup {setup_dir} \\
    -time 150 \\
    -gmx "$GMX" \\
    --no-analysis

if [ $? -ne 0 ]; then
    echo "ERROR: Production MD failed for replica {replica_num}!"
    exit 1
fi

echo "======================================"
echo "Replica {replica_num} completed!"
echo "======================================"
"""
        
        # Write replica SLURM script
        slurm_file = setup_dir.parent / f"submit_charmm_{pdb_name}_r{replica_num}.slurm"
        with open(slurm_file, 'w') as f:
            f.write(slurm_script)
        slurm_file.chmod(0o755)
        
        print(f"  Generated replica {replica_num} script: {slurm_file}")
        return slurm_file
    
    def generate_slurm_scripts(self, pdb_name, setup_dir, replica_dirs):
        """Generate all SLURM scripts for one system"""
        scripts = []
        
        # Generate setup script
        setup_script = self.generate_setup_script(pdb_name, setup_dir)
        scripts.append(setup_script)
        
        # Generate replica scripts
        for i, replica_dir in enumerate(replica_dirs, 1):
            replica_script = self.generate_replica_script(pdb_name, setup_dir, replica_dir, i)
            scripts.append(replica_script)
        
        return scripts
    
    def generate_system_submit_script(self, pdb_name, setup_dir, scripts):
        """Generate convenience script to submit all jobs for one system"""
        submit_script = setup_dir.parent / f"submit_charmm_{pdb_name}_all.sh"
        
        # Calculate number of replicas from scripts list (scripts[0] is setup, rest are replicas)
        num_replicas = len(scripts) - 1
        
        with open(submit_script, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(f"# Submit all jobs for {pdb_name}\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("\n")
            f.write("# Step 1: Submit setup job\n")
            f.write(f"SETUP_JOB=$(sbatch --parsable {scripts[0]})\n")
            f.write('echo "Setup job submitted: $SETUP_JOB"\n')
            f.write("\n")
            f.write("# Step 2: Submit 3 replica jobs in parallel (all depend on setup completion)\n")
            for i in range(1, num_replicas + 1):
                f.write(f"REPLICA{i}_JOB=$(sbatch --parsable --dependency=afterok:$SETUP_JOB {scripts[i]})\n")
                f.write(f'echo "Replica {i} job submitted (parallel): $REPLICA{i}_JOB"\n')
            f.write("\n")
            f.write(f'echo "All jobs submitted for {pdb_name}"\n')
            f.write('echo "3 replicas will run in parallel after setup completes"\n')
            f.write('echo "Monitor with: squeue -u $USER"\n')
        
        submit_script.chmod(0o755)
        print(f"  Generated submit-all script: {submit_script}")
        return submit_script
    
    def generate_master_submit_script(self, system_scripts):
        """Generate master submit script for all systems"""
        master_script = self.output_dir / "submit_charmm_all.sh"
        
        with open(master_script, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# Master submission script for all systems\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("\n")
            
            for system_script in system_scripts:
                f.write(f"bash {system_script}\n")
        
        # Make executable
        master_script.chmod(0o755)
        
        print("\n" + "=" * 70)
        print("Generated master submission script:")
        print(f"  {master_script}")
        print("\nTo submit all systems, run:")
        print(f"  bash {master_script}")
        print("\nOr submit individual systems:")
        for script in system_scripts:
            print(f"  bash {script}")
        print("=" * 70)
        
        return master_script
    
    def run_batch_processing(self):
        """Run batch processing for all PDB files"""
        print("=" * 60)
        print("pHLA MD Batch Processing")
        print("=" * 60)
        print(f"Data directory: {self.data_dir}")
        print(f"Output directory: {self.output_dir}")
        print(f"MDP directory: {self.mdp_dir}")
        print(f"Number of replicas: {self.num_replicas}")
        print("=" * 60)
        
        # Find all PDB files
        pdb_files = self.find_pdb_files()
        
        if not pdb_files:
            print("No PDB files found in data directory!")
            return
        
        # Generate SLURM scripts for each PDB
        all_system_scripts = []
        
        for pdb_file in pdb_files:
            print(f"\nProcessing: {pdb_file.name}")
            
            pdb_dir, setup_dir, replica_dirs = self.create_work_structure(pdb_file)
            print(f"  Working directory: {pdb_dir}")
            
            # Generate all scripts for this system (1 setup + 3 replicas)
            scripts = self.generate_slurm_scripts(
                pdb_file.stem, 
                setup_dir, 
                replica_dirs
            )
            
            # Generate system-specific submit script
            system_submit = self.generate_system_submit_script(
                pdb_file.stem,
                setup_dir,
                scripts
            )
            all_system_scripts.append(system_submit)
        
        # Generate master submission script
        master_script = self.generate_master_submit_script(all_system_scripts)
        
        # Write batch processing log
        with open(self.batch_log, 'w') as f:
            f.write(f"Batch Processing Log\n")
            f.write(f"Data directory: {self.data_dir}\n")
            f.write(f"Output directory: {self.output_dir}\n")
            f.write(f"MDP directory: {self.mdp_dir}\n")
            f.write(f"Number of replicas: {self.num_replicas}\n\n")
            
            f.write("Processed PDB files:\n")
            for pdb in pdb_files:
                f.write(f"  - {pdb.name}\n")
        
        print("=" * 60)
        print("Batch processing completed!")
        print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="pHLA MD Batch Processing (CHARMM36m) - Generate SLURM scripts for multiple PDB files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Example:
    python3 batch_submit_charmm.py -data ./data -output ./simulations -mdp ./mdp_files_charmm -n 3
        """)
    
    parser.add_argument('-data', '--data-dir', required=True,
                       help='Directory containing PDB files')
    parser.add_argument('-output', '--output-dir', required=True,
                       help='Output directory for all simulations')
    parser.add_argument('-mdp', '--mdp-dir', required=True,
                       help='Directory containing MDP files')
    parser.add_argument('-n', '--num-replicas', type=int, default=3,
                       help='Number of independent replicas per PDB (default: 3)')
    parser.add_argument('-gmx', '--gromacs', default='gmx',
                       help='GROMACS executable path (default: gmx)')
    
    args = parser.parse_args()
    
    # Validate directories
    if not Path(args.data_dir).exists():
        print(f"Error: Data directory not found: {args.data_dir}")
        sys.exit(1)
    
    if not Path(args.mdp_dir).exists():
        print(f"Error: MDP directory not found: {args.mdp_dir}")
        sys.exit(1)
    
    processor = BatchProcessor(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        mdp_dir=args.mdp_dir,
        gmx_path=args.gromacs,
        num_replicas=args.num_replicas
    )
    
    processor.run_batch_processing()


if __name__ == "__main__":
    main()
