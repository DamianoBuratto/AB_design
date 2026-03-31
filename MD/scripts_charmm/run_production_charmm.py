#!/usr/bin/env python3
"""
pHLA MD Production Run (CHARMM36m Force Field)
Run production MD simulation
Based on ANUBI framework
"""

import argparse
import os
import socket
import subprocess
import logging
from pathlib import Path
import shutil
import sys

class ProductionMD:
    def __init__(self, work_dir, replica_num, gmx_path="gmx"):
        self.work_dir = Path(work_dir)
        self.replica_num = replica_num
        self.gmx_path = gmx_path
        
        # Create working directory
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.work_dir / f'production_r{replica_num}.log'),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def run_gmx(self, cmd, description):
        """Run GROMACS command with logging"""
        self.logger.info(f"Running: {description}")
        self.logger.debug(f"Command: {cmd}")
        
        log_file = self.work_dir / f"{description.replace(' ', '_')}.log"
        
        # Ensure environment variables are passed to subprocess
        env = os.environ.copy()
        
        try:
            with open(log_file, 'w') as log:
                subprocess.run(
                    cmd,
                    shell=True,
                    check=True,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    cwd=self.work_dir,
                    env=env
                )
            
            self.logger.info(f"✓ {description} completed")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"✗ {description} failed with code {e.returncode}")
            self.logger.error(f"Check log: {log_file}")
            raise
    
    def prepare_system(self, setup_dir):
        """Copy necessary files from setup directory"""
        self.logger.info("Preparing system files...")
        
        setup_path = Path(setup_dir)
        required_files = {
            'npt.gro': 'npt.gro',
            'topol.top': 'topol.top',
            'state_npt.cpt': 'state_npt.cpt'
        }
        
        # Copy all topology include files (.itp)
        for itp_file in setup_path.glob('*.itp'):
            required_files[itp_file.name] = itp_file.name
        
        for src_name, dst_name in required_files.items():
            src = setup_path / src_name
            dst = self.work_dir / dst_name
            
            if src.exists():
                shutil.copy(src, dst)
                self.logger.info(f"  Copied: {src_name}")
            else:
                if src_name in ['npt.gro', 'topol.top']:
                    raise FileNotFoundError(f"Required file missing: {src}")
                else:
                    self.logger.warning(f"  Optional file missing: {src_name}")
    
    def run_production(self, md_mdp, sim_time=150):
        """Run production MD simulation"""
        self.logger.info("="*60)
        self.logger.info(f"PRODUCTION MD - Replica {self.replica_num}")
        self.logger.info(f"Simulation time: {sim_time} ns")
        self.logger.info("="*60)
        
        # Log node information
        hostname = socket.gethostname()
        cuda_visible = os.environ.get('CUDA_VISIBLE_DEVICES', 'Not set')
        
        self.logger.info(f"Running on node: {hostname}")
        self.logger.info(f"CUDA_VISIBLE_DEVICES: {cuda_visible}")
        
        # Check required files
        if not (self.work_dir / "npt.gro").exists():
            raise FileNotFoundError("Missing npt.gro (run preparation first)")
        if not (self.work_dir / "topol.top").exists():
            raise FileNotFoundError("Missing topol.top")
        
        # Prepare TPR
        cmd_grompp = (
            f"{self.gmx_path} grompp "
            f"-f {md_mdp} "
            f"-c npt.gro "
            f"-r npt.gro "
            f"-t state_npt.cpt "
            f"-p topol.top "
            f"-o md_r{self.replica_num}.tpr "
            f"-maxwarn 1"
        )
        
        self.run_gmx(cmd_grompp, f"grompp_md_r{self.replica_num}")
        
        # Run production MD with GPU acceleration
        cmd_mdrun = (
            f"{self.gmx_path} mdrun -v -deffnm md_r{self.replica_num} "
            f"-ntmpi 1 -ntomp 16 -gpu_id 0 "
            f"-nb gpu -bonded gpu -pme gpu "
            f"-pin on -pinstride 1"
        )
        
        self.run_gmx(cmd_mdrun, f"mdrun_r{self.replica_num}")
        
        # Verify outputs
        expected_outputs = [
            f"md_r{self.replica_num}.gro",
            f"md_r{self.replica_num}.xtc",
            f"md_r{self.replica_num}.edr"
        ]
        
        for output in expected_outputs:
            if not (self.work_dir / output).exists():
                self.logger.warning(f"Expected output missing: {output}")
    
    def analyze_trajectory(self):
        """Basic trajectory analysis"""
        self.logger.info("="*60)
        self.logger.info("TRAJECTORY ANALYSIS")
        self.logger.info("="*60)
        
        tpr = f"md_r{self.replica_num}.tpr"
        xtc = f"md_r{self.replica_num}.xtc"
        
        # RMSD analysis
        cmd_rmsd = (
            f'echo "4 4" | {self.gmx_path} rms '
            f'-s {tpr} -f {xtc} '
            f'-o rmsd_md_r{self.replica_num}.xvg -tu ns'
        )
        
        try:
            self.run_gmx(cmd_rmsd, f"rmsd_r{self.replica_num}")
        except:
            self.logger.warning("RMSD analysis failed (optional)")
        
        # Energy analysis
        cmd_energy = (
            f'echo "10 0" | {self.gmx_path} energy '
            f'-f md_r{self.replica_num}.edr '
            f'-o energy_md_r{self.replica_num}.xvg'
        )
        
        try:
            self.run_gmx(cmd_energy, f"energy_r{self.replica_num}")
        except:
            self.logger.warning("Energy analysis failed (optional)")
    
    def run_all(self, setup_dir, md_mdp, sim_time=150, run_analysis=True):
        """Run complete production MD workflow"""
        self.logger.info("="*60)
        self.logger.info("pHLA MD PRODUCTION RUN")
        self.logger.info(f"Replica: {self.replica_num}")
        self.logger.info(f"Working directory: {self.work_dir}")
        self.logger.info("="*60)
        
        try:
            self.prepare_system(setup_dir)
            self.run_production(md_mdp, sim_time)
            
            if run_analysis:
                self.analyze_trajectory()
            
            self.logger.info("="*60)
            self.logger.info("✓ PRODUCTION MD COMPLETED SUCCESSFULLY")
            self.logger.info("="*60)
            self.logger.info(f"Output files in: {self.work_dir}")
            self.logger.info(f"  - md_r{self.replica_num}.xtc (trajectory)")
            self.logger.info(f"  - md_r{self.replica_num}.edr (energies)")
            self.logger.info(f"  - md_r{self.replica_num}.gro (final structure)")
            
        except Exception as e:
            self.logger.error("="*60)
            self.logger.error("✗ PRODUCTION MD FAILED")
            self.logger.error(str(e))
            self.logger.error("="*60)
            raise

def main():
    parser = argparse.ArgumentParser(
        description='pHLA MD Production Run - Run production MD simulation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example:
    python run_production.py -w ./replica_1 -md production.mdp -r 1 -setup ./setup
        '''
    )
    
    parser.add_argument('-w', '--workdir', required=True,
                       help='Working directory for this replica')
    parser.add_argument('-md', '--md-mdp', required=True,
                       help='MDP file for production MD')
    parser.add_argument('-r', '--replica', type=int, required=True,
                       help='Replica number (for output naming)')
    parser.add_argument('-setup', '--setup-dir', 
                       help='Setup directory containing npt.gro and topol.top (default: ../setup)')
    parser.add_argument('-time', '--time', type=int, default=150,
                       help='Simulation time in nanoseconds (default: 150)')
    parser.add_argument('-gmx', '--gromacs', default='gmx',
                       help='GROMACS executable path (default: gmx)')
    parser.add_argument('--no-analysis', action='store_true',
                       help='Skip analysis after MD run')
    
    args = parser.parse_args()
    
    # Determine setup directory
    if args.setup_dir:
        setup_dir = args.setup_dir
    else:
        setup_dir = Path(args.workdir).parent / "setup"
    
    # Run production MD
    prod = ProductionMD(
        work_dir=args.workdir,
        replica_num=args.replica,
        gmx_path=args.gromacs
    )
    
    prod.run_all(
        setup_dir=setup_dir,
        md_mdp=args.md_mdp,
        sim_time=args.time,
        run_analysis=not args.no_analysis
    )

if __name__ == "__main__":
    main()
