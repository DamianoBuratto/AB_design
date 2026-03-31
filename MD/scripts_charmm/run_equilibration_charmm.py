#!/usr/bin/env python3
"""
pHLA MD Equilibration (CHARMM36m Force Field)
Run NVT and NPT equilibration phases
Based on ANUBI framework
"""

import argparse
import os
import shutil
import socket
import subprocess
import logging
from pathlib import Path
import sys

class Equilibration:
    def __init__(self, work_dir, gmx_path="gmx"):
        self.work_dir = Path(work_dir)
        self.gmx_path = gmx_path
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.work_dir / 'equilibration.log'),
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
    
    def nvt_equilibration(self, nvt_mdp, input_gro="system_minim.gro"):
        """NVT (constant Number, Volume, Temperature) equilibration"""
        self.logger.info("="*60)
        self.logger.info("NVT EQUILIBRATION")
        self.logger.info("="*60)
        
        # Check required files
        if not (self.work_dir / input_gro).exists():
            raise FileNotFoundError(f"Missing input: {input_gro}")
        if not (self.work_dir / "topol.top").exists():
            raise FileNotFoundError("Missing topology: topol.top")
        
        # Prepare TPR
        cmd_grompp = (
            f"{self.gmx_path} grompp "
            f"-f {nvt_mdp} "
            f"-c {input_gro} "
            f"-r {input_gro} "
            f"-p topol.top "
            f"-o nvt.tpr "
            f"-maxwarn 1"
        )
        
        self.run_gmx(cmd_grompp, "grompp_nvt")
        
        # Run NVT with GPU acceleration
        cmd_mdrun = (
            f"{self.gmx_path} mdrun -v -deffnm nvt "
            f"-ntmpi 1 -ntomp 16 -gpu_id 0 "
            f"-nb gpu -bonded gpu -pme gpu "
            f"-pin on -pinstride 1"
        )
        
        self.run_gmx(cmd_mdrun, "mdrun_nvt")
        
        if not (self.work_dir / "nvt.gro").exists():
            raise FileNotFoundError("Missing output: nvt.gro")
    
    def npt_equilibration(self, npt_mdp):
        """NPT (constant Number, Pressure, Temperature) equilibration"""
        self.logger.info("="*60)
        self.logger.info("NPT EQUILIBRATION")
        self.logger.info("="*60)
        
        # Check required files
        if not (self.work_dir / "nvt.gro").exists():
            raise FileNotFoundError("Missing input: nvt.gro (run NVT first)")
        if not (self.work_dir / "nvt.cpt").exists():
            self.logger.warning("Missing nvt.cpt - will start NPT without NVT checkpoint")
        
        # Prepare TPR
        cmd_grompp = (
            f"{self.gmx_path} grompp "
            f"-f {npt_mdp} "
            f"-c nvt.gro "
            f"-r nvt.gro "
            f"-t nvt.cpt "
            f"-p topol.top "
            f"-o npt.tpr "
            f"-maxwarn 1"
        )
        
        self.run_gmx(cmd_grompp, "grompp_npt")
        
        # Run NPT with GPU acceleration
        cmd_mdrun = (
            f"{self.gmx_path} mdrun -v -deffnm npt "
            f"-ntmpi 1 -ntomp 16 -gpu_id 0 "
            f"-nb gpu -bonded gpu -pme gpu "
            f"-pin on -pinstride 1"
        )
        
        self.run_gmx(cmd_mdrun, "mdrun_npt")
        
        if not (self.work_dir / "npt.gro").exists():
            raise FileNotFoundError("Missing output: npt.gro")
        
        # Rename checkpoint for production MD
        if (self.work_dir / "npt.cpt").exists():
            shutil.copy(self.work_dir / "npt.cpt", self.work_dir / "state_npt.cpt")
    
    def run_all(self, nvt_mdp, npt_mdp):
        """Run complete equilibration"""
        hostname = socket.gethostname()
        cuda_visible = os.environ.get('CUDA_VISIBLE_DEVICES', 'Not set')

        self.logger.info("="*60)
        self.logger.info("pHLA MD EQUILIBRATION")
        self.logger.info(f"Working directory: {self.work_dir}")
        self.logger.info(f"Node: {hostname}")
        self.logger.info(f"CUDA_VISIBLE_DEVICES: {cuda_visible}")
        self.logger.info("="*60)
        
        try:
            self.nvt_equilibration(nvt_mdp)
            self.npt_equilibration(npt_mdp)
            
            self.logger.info("="*60)
            self.logger.info("✓ EQUILIBRATION COMPLETED SUCCESSFULLY")
            self.logger.info("="*60)
            self.logger.info(f"Output files in: {self.work_dir}")
            self.logger.info("  - nvt.gro (NVT equilibrated)")
            self.logger.info("  - npt.gro (NPT equilibrated)")
            self.logger.info("  - state_npt.cpt (checkpoint for production)")
            
        except Exception as e:
            self.logger.error("="*60)
            self.logger.error("✗ EQUILIBRATION FAILED")
            self.logger.error(str(e))
            self.logger.error("="*60)
            raise

def main():
    parser = argparse.ArgumentParser(
        description='pHLA MD Equilibration - Run NVT and NPT equilibration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example:
    python run_equilibration.py -w ./work -nvt nvt.mdp -npt npt.mdp
        '''
    )
    
    parser.add_argument('-w', '--workdir', required=True,
                       help='Working directory (must contain minim.gro and topol.top)')
    parser.add_argument('-nvt', '--nvt-mdp', required=True,
                       help='MDP file for NVT equilibration')
    parser.add_argument('-npt', '--npt-mdp', required=True,
                       help='MDP file for NPT equilibration')
    parser.add_argument('-gmx', '--gromacs', default='gmx',
                       help='GROMACS executable path (default: gmx)')
    
    args = parser.parse_args()
    
    # Run equilibration
    eq = Equilibration(
        work_dir=args.workdir,
        gmx_path=args.gromacs
    )
    
    eq.run_all(
        nvt_mdp=args.nvt_mdp,
        npt_mdp=args.npt_mdp
    )

if __name__ == "__main__":
    main()
