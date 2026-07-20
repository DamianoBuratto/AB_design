#!/usr/bin/env python3
"""
pHLA MD System Setup (CHARMM36m Force Field)
- Generates topology with CHARMM36m
- Generates LOCAL position restraints for HLA chain (Fixes grompp errors)
- Solvates, adds ions, and minimizes energy
"""

import argparse
import subprocess
import logging
from pathlib import Path
import shutil
import os

class SystemSetup:
    def __init__(self, pdb_file, work_dir, gmx_path="gmx", 
                 forcefield="charmm36m-mut", water="tip3p"):
        self.pdb_file = Path(pdb_file).resolve()
        self.work_dir = Path(work_dir).resolve()
        self.gmx_path = gmx_path
        self.forcefield = forcefield
        self.water = water
        
        # Create working directory
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.work_dir / 'setup_system.log', mode='w'),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)

    def run_gmx(self, cmd, description, input_str=None):
        """Run GROMACS command with logging and error handling"""
        self.logger.info(f"Running: {description}")
        self.logger.debug(f"Command: {cmd}")
        
        log_file = self.work_dir / f"{description.replace(' ', '_')}.log"
        
        # Ensure environment variables are passed to subprocess
        env = os.environ.copy()
        
        try:
            with open(log_file, 'w') as log:
                if input_str:
                    process = subprocess.Popen(
                        cmd,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        cwd=self.work_dir,
                        env=env,
                        text=True
                    )
                    process.communicate(input=input_str)
                    if process.returncode != 0:
                        raise subprocess.CalledProcessError(process.returncode, cmd)
                else:
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
            self.logger.error(f"Check log for details: {log_file}")
            raise

    def get_last_group_index(self, ndx_file):
        """Helper to find the index of the last group in an ndx file"""
        last_idx = 0
        with open(self.work_dir / ndx_file, 'r') as f:
            for line in f:
                if line.strip().startswith('[') and line.strip().endswith(']'):
                    last_idx += 1
        return last_idx - 1 if last_idx > 0 else 0

    def step0_preprocess_pdb(self, peptide_sequence="HMTEVVRHC"):
        """
        Preprocess PDB before topology generation:
        1. Split chain T into P (peptide) and A (HLA)
        2. Renumber residues starting from 1 for each chain
        3. Add TER records between chains
        4. Convert all HIS to HSE (CHARMM36m ε-protonated histidine)
        
        Args:
            peptide_sequence: Expected peptide sequence (default: HMTEVVRHC)
        """
        self.logger.info("="*60)
        self.logger.info("STEP 0: PDB Preprocessing")
        self.logger.info("="*60)
        
        input_pdb = self.work_dir / self.pdb_file.name
        output_pdb = self.work_dir / "processed.pdb"
        
        # Copy original PDB to work directory
        if not input_pdb.exists():
            shutil.copy(self.pdb_file, input_pdb)
        
        self.logger.info(f"Processing PDB: {input_pdb.name}")
        self.logger.info(f"Peptide sequence: {peptide_sequence} ({len(peptide_sequence)} residues)")
        
        # Parse PDB file
        with open(input_pdb, 'r') as f:
            lines = f.readlines()
        
        # Collect ATOM/HETATM records by chain
        chain_t_atoms = []
        other_atoms = []
        header_lines = []
        
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21:22].strip()
                if chain_id == 'T':
                    chain_t_atoms.append(line)
                else:
                    other_atoms.append(line)
            elif line.startswith(('CRYST', 'MODEL', 'REMARK', 'TITLE', 'HEADER')):
                header_lines.append(line)
        
        if not chain_t_atoms:
            raise ValueError("No chain T found in PDB file!")
        
        self.logger.info(f"Found {len(chain_t_atoms)} atoms in chain T")
        
        # Identify peptide length from sequence
        peptide_length = len(peptide_sequence)
        
        # Group atoms by residue to split peptide/HLA
        residues_t = {}
        for line in chain_t_atoms:
            resnum = int(line[22:26].strip())
            if resnum not in residues_t:
                residues_t[resnum] = []
            residues_t[resnum].append(line)
        
        sorted_resnums = sorted(residues_t.keys())
        self.logger.info(f"Chain T contains {len(sorted_resnums)} residues: {sorted_resnums[0]}-{sorted_resnums[-1]}")
        
        # Split into peptide (first N residues) and HLA (remaining)
        peptide_resnums = sorted_resnums[:peptide_length]
        hla_resnums = sorted_resnums[peptide_length:]
        
        self.logger.info(f"Peptide: {len(peptide_resnums)} residues → Chain P")
        self.logger.info(f"HLA: {len(hla_resnums)} residues → Chain A")
        
        # Write processed PDB
        with open(output_pdb, 'w') as f:
            # Write header
            for line in header_lines:
                f.write(line)
            
            # Write Peptide chain (P)
            atom_serial = 1
            for new_resnum, old_resnum in enumerate(peptide_resnums, start=1):
                for atom_line in residues_t[old_resnum]:
                    # Convert HIS to HSE (CHARMM naming)
                    resname = atom_line[17:20].strip()
                    if resname == 'HIS':
                        resname = 'HSE'
                    
                    # Rewrite line with new chain, residue number, atom serial
                    new_line = (
                        f"{atom_line[:6]}"           # ATOM/HETATM
                        f"{atom_serial:>5} "          # Atom serial
                        f"{atom_line[12:17]}"         # Atom name + altLoc
                        f"{resname:>3} "              # Residue name (HSE conversion)
                        f"P"                          # Chain ID = P (peptide)
                        f"{new_resnum:>4}    "        # Residue number (renumbered)
                        f"{atom_line[30:54]}"         # Coordinates
                        f"{atom_line[54:76]}"         # Occupancy, temp factor
                        f"{atom_line[76:]}"           # Element, charge
                    )
                    f.write(new_line)
                    atom_serial += 1
            
            f.write("TER\n")
            
            # Write HLA chain (A)
            for new_resnum, old_resnum in enumerate(hla_resnums, start=1):
                for atom_line in residues_t[old_resnum]:
                    # Convert HIS to HSE (CHARMM naming)
                    resname = atom_line[17:20].strip()
                    if resname == 'HIS':
                        resname = 'HSE'
                    
                    # Rewrite line with new chain, residue number, atom serial
                    new_line = (
                        f"{atom_line[:6]}"           # ATOM/HETATM
                        f"{atom_serial:>5} "          # Atom serial
                        f"{atom_line[12:17]}"         # Atom name + altLoc
                        f"{resname:>3} "              # Residue name (HSE conversion)
                        f"A"                          # Chain ID = A (HLA)
                        f"{new_resnum:>4}    "        # Residue number (renumbered from 1)
                        f"{atom_line[30:54]}"         # Coordinates
                        f"{atom_line[54:76]}"         # Occupancy, temp factor
                        f"{atom_line[76:]}"           # Element, charge
                    )
                    f.write(new_line)
                    atom_serial += 1
            
            f.write("TER\n")
            
            # Write other chains if present
            for line in other_atoms:
                f.write(line)
            
            f.write("END\n")
        
        self.logger.info(f"✓ Preprocessed PDB written: {output_pdb.name}")
        self.logger.info(f"  - Chain P (peptide): {len(peptide_resnums)} residues")
        self.logger.info(f"  - Chain A (HLA): {len(hla_resnums)} residues")
        self.logger.info(f"  - All HIS converted to HSE")
        
        # Update self.pdb_file to use processed version
        self.pdb_file = output_pdb
        
        return output_pdb

    def step1_generate_topology(self):
        """Generate topology using pdb2gmx with CHARMM36m force field."""
        self.logger.info("="*60)
        self.logger.info("STEP 1: Generate Topology")
        self.logger.info("="*60)
        
        local_pdb = self.work_dir / self.pdb_file.name
        if not local_pdb.exists():
            shutil.copy(self.pdb_file, local_pdb)
        
        # -ignh ignores H atoms in input and rebuilds them (good for naming consistency)
        cmd = (
            f"{self.gmx_path} pdb2gmx "
            f"-f {local_pdb.name} "
            f"-o system.pdb "
            f"-p topol.top "
            f"-ignh "
            f"-ff {self.forcefield} "
            f"-water {self.water}"
        )
        
        self.run_gmx(cmd, "pdb2gmx")
        
        if not (self.work_dir / "system.pdb").exists():
            raise FileNotFoundError("pdb2gmx failed: system.pdb not created")

    def step1b_generate_hla_restraints(self, hla_chain="A", 
                                        hla_residues="4-11,94-102",
                                        force_constant=200):
        """
        Generate position restraints on the ISOLATED HLA chain to match local topology indices.
        
        After preprocessing, HLA is chain A with residues renumbered from 1.
        Restraints target the two beta-sheet floor regions of the binding groove:
          A:4-11   (beta-sheet region 1)
          A:94-102 (beta-sheet region 2)
        
        Args:
            hla_chain: HLA chain identifier (default: 'A' after preprocessing)
            hla_residues: Comma-separated residue ranges in chain A numbering (default: '4-11,94-102')
            force_constant: Restraint force constant in kJ/(mol·nm²)
        """
        self.logger.info("="*60)
        self.logger.info("STEP 1B: Generate HLA Position Restraints (Local Indices)")
        self.logger.info("="*60)
        
        # 1. Create index file to identify Chain T in the system
        # "chain T" might differ in naming, but pdb2gmx usually preserves chain IDs if possible
        # or assigns them sequentially. We assume standard pdb2gmx output.
        
        self.logger.info(f"Isolating Chain {hla_chain}...")
        
        # Create an index for the whole system, keep only system group + target chain
        cmd_ndx_sys = f"{self.gmx_path} make_ndx -f system.pdb -o index_sys.ndx"
        # Keep group 0, add chain, quit
        self.run_gmx(cmd_ndx_sys, "make_ndx_system", input_str=f"keep 0\nchain {hla_chain}\nq\n")
        
        # 2. Extract Chain T to a separate PDB
        # The group name for chain T is typically "chT" or similar.
        # We find the last group created by make_ndx.

        # Use "keep 0" to ensure chain group is always at position 1
        chain_group_idx = 1
        
        isolated_pdb = f"chain_{hla_chain}_isolated.pdb"
        cmd_extract = (
            f"{self.gmx_path} trjconv "
            f"-f system.pdb "
            f"-n index_sys.ndx "
            f"-o {isolated_pdb} "
            f"-s system.pdb"
        )
        
        self.run_gmx(cmd_extract, f"extract_chain_{hla_chain}", input_str=f"{chain_group_idx}\n")
        
        # 3. Generate restraints on the ISOLATED PDB
        # This ensures atom indices start from 1, matching the .itp file.
        self.logger.info(f"Generating restraints on isolated {isolated_pdb}...")
        
        residue_ranges = [r.strip() for r in hla_residues.split(',')]
        
        # Strategy for make_ndx on isolated chain:
        # 1. Create group for target residues
        # 2. Create group for CA atoms
        # 3. Intersect them
        
        # Construct input for make_ndx
        # We first check default groups to know what to delete or ignore, 
        # but it's safer to just create new groups at the end.
        
        ndx_cmds = []
        # Union of all residue ranges: "r 4-11 | r 94-102"
        res_cmd = " | ".join([f"r {r}" for r in residue_ranges])
        ndx_cmds.append(res_cmd) 
        
        # Robust Clean & Build Sequence:
        # 1. keep 0 -> Group 0 = System (all atoms)
        # 2. r 4-11 | r 94-102 -> Group 1 = target residues
        # 3. 1 & a CA -> Group 2 = target CA atoms (used for genrestr)

        clean_cmds = ["keep 0"]  # Keep only group 0 (System), remove all others
        build_cmds = [
            res_cmd,      # Creates Group 1 (target residues)
            "1 & a CA",   # Creates Group 2 (target CA atoms)
            "q"
        ]
        
        full_input = "\n".join(clean_cmds + build_cmds) + "\n"
        
        self.run_gmx(
            f"{self.gmx_path} make_ndx -f {isolated_pdb} -o index_hla_local.ndx", 
            "make_ndx_local", 
            input_str=full_input
        )
        
        # The target group should be the last one (index 1 if cleaning worked, or just last)
        target_group_idx = self.get_last_group_index("index_hla_local.ndx")
        
        # 4. Generate the .itp file
        cmd_genrestr = (
            f"{self.gmx_path} genrestr "
            f"-f {isolated_pdb} "
            f"-n index_hla_local.ndx "
            f"-o posre_HLA_CA.itp "
            f"-fc {force_constant} {force_constant} {force_constant}"
        )
        
        self.run_gmx(cmd_genrestr, "genrestr_hla", input_str=f"{target_group_idx}\n")
        
        if not (self.work_dir / "posre_HLA_CA.itp").exists():
            raise FileNotFoundError("posre_HLA_CA.itp was not generated.")
            
        # 5. Modify the specific CHAIN topology file
        self.modify_topology_for_hla_restraints(hla_chain)

    def modify_topology_for_hla_restraints(self, hla_chain):
        """Append the restraint include to the specific chain .itp file"""
        # pdb2gmx typically names it topol_Protein_chain_T.itp
        chain_itp = self.work_dir / f"topol_Protein_chain_{hla_chain}.itp"
        
        if not chain_itp.exists():
            # Try finding it if name differs
            candidates = list(self.work_dir.glob(f"*chain_{hla_chain}.itp"))
            if candidates:
                chain_itp = candidates[0]
            else:
                self.logger.error(f"Could not find topology file for chain {hla_chain}")
                raise FileNotFoundError(f"Missing topol_Protein_chain_{hla_chain}.itp")
        
        self.logger.info(f"Adding restraints to {chain_itp.name}...")
        
        with open(chain_itp, 'a') as f:
            f.write("\n")
            f.write("; Custom HLA position restraints (CA atoms only)\n")
            f.write("#ifdef POSRES_HLA\n")
            f.write('#include "posre_HLA_CA.itp"\n')
            f.write("#endif\n")

    def step2_build_box(self, box_type="cubic", distance=1.2):
        """Build simulation box"""
        self.logger.info("="*60)
        self.logger.info("STEP 2: Build Box")
        self.logger.info("="*60)
        
        cmd = (
            f"{self.gmx_path} editconf "
            f"-f system.pdb "
            f"-o system_box.gro "
            f"-c -bt {box_type} -d {distance}"
        )
        self.run_gmx(cmd, "editconf")

    def step3_solvate(self):
        """Add water molecules"""
        self.logger.info("="*60)
        self.logger.info("STEP 3: Solvate System")
        self.logger.info("="*60)
        
        cmd = (
            f"{self.gmx_path} solvate "
            f"-cp system_box.gro "
            f"-cs spc216.gro "
            f"-o system_water.gro "
            f"-p topol.top"
        )
        self.run_gmx(cmd, "solvate")

    def step4_add_ions(self, ions_mdp, concentration=0.10):
        """Add ions for neutralization and salt concentration"""
        self.logger.info("="*60)
        self.logger.info("STEP 4: Add Ions")
        self.logger.info("="*60)
        
        # Prepare TPR for genion
        cmd_grompp = (
            f"{self.gmx_path} grompp "
            f"-f {ions_mdp} "
            f"-c system_water.gro "
            f"-p topol.top "
            f"-o system_ions.tpr "
            f"-maxwarn 2"
        )
        self.run_gmx(cmd_grompp, "grompp_ions")
        
        # Add ions (replace SOL molecules)
        cmd_genion = (
            f"{self.gmx_path} genion "
            f"-s system_ions.tpr "
            f"-o system_ions.gro "
            f"-p topol.top "
            f"-pname K "
            f"-nname CL "
            f"-neutral "
            f"-conc {concentration}"
        )
        # Select Group 13 (SOL) usually, but safer to use name
        self.run_gmx(cmd_genion, "genion", input_str="SOL\n")

    def step5_energy_minimization(self, minim_mdp, nsteps=50000):
        """Energy minimization"""
        self.logger.info("="*60)
        self.logger.info("STEP 5: Energy Minimization")
        self.logger.info("="*60)
        
        cmd_grompp = (
            f"{self.gmx_path} grompp "
            f"-f {minim_mdp} "
            f"-c system_ions.gro "
            f"-p topol.top "
            f"-o em.tpr "
            f"-maxwarn 1"
        )
        self.run_gmx(cmd_grompp, "grompp_em")
        
        # Run minimization
        cmd_mdrun = f"{self.gmx_path} mdrun -v -deffnm em"
        self.run_gmx(cmd_mdrun, "mdrun_em")
        
        # Rename for consistency
        if (self.work_dir / "em.gro").exists():
            shutil.copy(self.work_dir / "em.gro", self.work_dir / "system_minim.gro")
            self.logger.info("✓ Energy minimization completed")
        else:
            raise FileNotFoundError("Missing output: em.gro")
        
        # Check convergence
        self.check_minimization()

    def check_minimization(self):
        """Check if energy minimization converged"""
        em_log = self.work_dir / "em.log"
        if not em_log.exists(): 
            return
            
        with open(em_log, 'r') as f:
            for line in f:
                if "Norm of force" in line:
                    self.logger.info(f"Minimization status: {line.strip()}")

    def run_all(self, ions_mdp, minim_mdp, box_type="triclinic", 
                distance=1.2, concentration=0.10, em_steps=50000):
        try:
            self.step0_preprocess_pdb()  # Split chains, renumber, HIS→HSE
            self.step1_generate_topology()
            # Generate HLA restraints BEFORE building the box, but using isolated logic
            self.step1b_generate_hla_restraints()
            self.step2_build_box(box_type, distance)
            self.step3_solvate()
            self.step4_add_ions(ions_mdp, concentration)
            self.step5_energy_minimization(minim_mdp, em_steps)
            
            self.logger.info("="*60)
            self.logger.info("✓ SYSTEM SETUP COMPLETED SUCCESSFULLY")
            self.logger.info("="*60)
            
        except Exception as e:
            self.logger.error("="*60)
            self.logger.error("✗ SYSTEM SETUP FAILED")
            self.logger.error(str(e))
            self.logger.error("="*60)
            raise

def main():
    parser = argparse.ArgumentParser(
        description='pHLA MD System Setup (CHARMM36m)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-i', '--input', required=True, help='Input PDB file')
    parser.add_argument('-w', '--workdir', required=True, help='Working directory')
    parser.add_argument('-ions', '--ions-mdp', required=True, help='MDP file for ion addition')
    parser.add_argument('-em', '--minim-mdp', required=True, help='MDP file for energy minimization')
    parser.add_argument('-gmx', '--gromacs', default='gmx', help='GROMACS executable path')
    parser.add_argument('-ff', '--forcefield', default='charmm36m-mut', help='Force field (default: charmm36m-mut from custom path)')
    parser.add_argument('-water', '--water-model', default='tip3p', help='Water model (TIP3P for CHARMM)')
    parser.add_argument('-box', '--box-type', default='triclinic', help='Box type (default: triclinic)')
    parser.add_argument('-d', '--distance', type=float, default=1.2, help='Distance in nm (default: 1.2)')
    parser.add_argument('-conc', '--concentration', type=float, default=0.10, help='Salt concentration')
    parser.add_argument('--em-steps', type=int, default=50000, help='EM steps')
    
    args = parser.parse_args()
    
    setup = SystemSetup(
        pdb_file=args.input,
        work_dir=args.workdir,
        gmx_path=args.gromacs,
        forcefield=args.forcefield,
        water=args.water_model
    )
    
    setup.run_all(
        ions_mdp=args.ions_mdp,
        minim_mdp=args.minim_mdp,
        box_type=args.box_type,
        distance=args.distance,
        concentration=args.concentration,
        em_steps=args.em_steps
    )

if __name__ == "__main__":
    main()