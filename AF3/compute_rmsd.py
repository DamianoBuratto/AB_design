#!/usr/bin/env python3
"""
Compute antibody and CDR RMSD between AF3 CIF results and RF2 PDB result.

This script:
1. Reads all seed-sample CIF files from AF3 output
2. Reads RF2 PDB file (chains: H, L, T)
3. Parses CDR regions from PDB REMARK PDBinfo-LABEL annotations
4. For each seed-sample: aligns using fixed chains AF3(A+C) vs RF2(T)
5. Computes RMSD of mobile chains: AF3(H+L) vs RF2(H+L) (full antibody)
6. Computes RMSD of CDR regions only: H1/H2/H3/L1/L2/L3
7. Updates the ranking_scores.csv with antibody_rmsd and cdr_rmsd columns

Usage:
    python3 compute_rmsd.py <output_base_dir> <rf2_pdb> <ranking_csv>
    
Example:
    python3 compute_rmsd.py \\
        jobs/pHLA_34_2/outputs/Antibody-pHLA_Complex \\
        jobs/pHLA_34_2/design_34_dldesign_2_best.pdb \\
        jobs/pHLA_34_2/outputs/Antibody-pHLA_Complex/Antibody-pHLA_Complex_ranking_scores.csv
"""

import sys
import os
import csv
import glob
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio.PDB.Polypeptide import is_aa
import numpy as np
import re


def parse_cdr_regions_from_pdb(pdb_file):
    """
    Parse CDR regions from PDB REMARK PDBinfo-LABEL annotations.
    
    Returns:
        dict: {'H': {1: [res_nums], 2: [res_nums], 3: [res_nums]},
               'L': {1: [res_nums], 2: [res_nums], 3: [res_nums]}}
    """
    cdr_regions = {'H': {1: [], 2: [], 3: []}, 'L': {1: [], 2: [], 3: []}}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('REMARK PDBinfo-LABEL:'):
                # Format: REMARK PDBinfo-LABEL:   25 H1
                match = re.search(r'REMARK PDBinfo-LABEL:\s+(\d+)\s+([HL])(\d)', line)
                if match:
                    res_num = int(match.group(1))
                    chain = match.group(2)  # H or L
                    cdr_num = int(match.group(3))  # 1, 2, or 3
                    cdr_regions[chain][cdr_num].append(res_num)
    
    return cdr_regions


def get_ca_atoms_from_chains(structure, model_id, chain_ids):
    """Extract CA atoms from specified chains for alignment."""
    atoms = []
    try:
        model = structure[model_id]
        for chain_id in chain_ids:
            if chain_id not in model:
                print(f"Warning: Chain {chain_id} not found in structure", file=sys.stderr)
                continue
            chain = model[chain_id]
            for residue in chain:
                if is_aa(residue, standard=True):
                    if 'CA' in residue:
                        atoms.append(residue['CA'])
    except Exception as e:
        print(f"Error extracting CA atoms from chains {chain_ids}: {e}", file=sys.stderr)
    return atoms


def get_pHLA_atoms_for_alignment(rf2_structure, af3_structure, model_id=0):
    """
    Extract pHLA CA atoms for alignment, properly mapping RF2 T chain to AF3 C+A chains.
    
    RF2 T chain structure: [peptide + HLA] concatenated
    AF3 structure: C chain (peptide), A chain (HLA)
    
    Returns:
        rf2_phla_atoms: List of CA atoms from RF2 T chain (peptide + HLA)
        af3_phla_atoms: List of CA atoms from AF3 C+A chains (peptide + HLA)
    """
    # Get AF3 C and A chains first to determine split point
    af3_model = af3_structure[model_id]
    af3_c_residues = [res for res in af3_model['C'] if is_aa(res, standard=True)] if 'C' in af3_model else []
    af3_a_residues = [res for res in af3_model['A'] if is_aa(res, standard=True)] if 'A' in af3_model else []
    
    peptide_len = len(af3_c_residues)
    
    # Get RF2 T chain and split
    rf2_model = rf2_structure[model_id]
    rf2_t_residues = [res for res in rf2_model['T'] if is_aa(res, standard=True)] if 'T' in rf2_model else []
    
    # Split RF2 T chain: first peptide_len residues = peptide, rest = HLA
    rf2_peptide_residues = rf2_t_residues[:peptide_len]
    rf2_hla_residues = rf2_t_residues[peptide_len:]
    
    # Extract CA atoms
    rf2_peptide_ca = [res['CA'] for res in rf2_peptide_residues if 'CA' in res]
    rf2_hla_ca = [res['CA'] for res in rf2_hla_residues if 'CA' in res]
    af3_peptide_ca = [res['CA'] for res in af3_c_residues if 'CA' in res]
    af3_hla_ca = [res['CA'] for res in af3_a_residues if 'CA' in res]
    
    # Combine peptide + HLA
    rf2_phla_atoms = rf2_peptide_ca + rf2_hla_ca
    af3_phla_atoms = af3_peptide_ca + af3_hla_ca
    
    return rf2_phla_atoms, af3_phla_atoms


def get_all_atoms_from_chains(structure, model_id, chain_ids):
    """Extract all heavy atoms from specified chains."""
    atoms = []
    try:
        model = structure[model_id]
        for chain_id in chain_ids:
            if chain_id not in model:
                print(f"Warning: Chain {chain_id} not found in structure", file=sys.stderr)
                continue
            chain = model[chain_id]
            for residue in chain:
                if is_aa(residue, standard=True):
                    for atom in residue:
                        # Skip hydrogen atoms
                        if not atom.element.startswith('H'):
                            atoms.append(atom)
    except Exception as e:
        print(f"Error extracting atoms from chains {chain_ids}: {e}", file=sys.stderr)
    return atoms


def get_ca_atoms_from_cdr_regions(structure, model_id, cdr_regions, is_af3=False):
    """
    Extract CA atoms from CDR regions (H1/H2/H3/L1/L2/L3).
    
    Args:
        structure: Biopython structure
        model_id: Model index
        cdr_regions: Dict from parse_cdr_regions_from_pdb() with RF2 residue numbering
        is_af3: If True, convert L chain residue numbers (RF2 L chain starts higher than AF3)
    
    Returns:
        List of CA atoms from all CDR regions
    """
    atoms = []
    try:
        model = structure[model_id]
        
        # Determine L chain offset for AF3
        # RF2: L chain typically starts at 125 or 126
        # AF3: L chain always starts at 1
        # Need to calculate offset dynamically
        l_chain_offset = 0
        if is_af3 and 'L' in model:
            # Get first residue number of L chain in this structure
            l_chain_residues = [res for res in model['L'] if is_aa(res, standard=True)]
            if l_chain_residues:
                af3_l_first = min(res.id[1] for res in l_chain_residues)
                # Get RF2 L chain first residue from CDR regions
                all_l_residues = []
                for cdr_num in [1, 2, 3]:
                    all_l_residues.extend(cdr_regions['L'][cdr_num])
                if all_l_residues:
                    rf2_l_min = min(all_l_residues)
                    # offset = RF2_first - AF3_first (e.g., 125 - 1 = 124)
                    l_chain_offset = rf2_l_min - 1 - (rf2_l_min % 100)  # Approximate
                    # Better: assume AF3 starts at 1, RF2 starts at H_chain_length + 1
                    # For most cases: offset = number of H chain residues
                    if 'H' in model:
                        h_chain_residues = [res for res in model['H'] if is_aa(res, standard=True)]
                        if h_chain_residues:
                            l_chain_offset = len(h_chain_residues)
        
        for chain_id in ['H', 'L']:
            if chain_id not in model:
                continue
            
            chain = model[chain_id]
            
            # Collect all CDR residue numbers for this chain
            cdr_res_nums = set()
            for cdr_num in [1, 2, 3]:
                res_nums = cdr_regions[chain_id][cdr_num]
                if chain_id == 'L' and is_af3 and l_chain_offset > 0:
                    # Convert L chain residue numbers for AF3
                    res_nums = [r - l_chain_offset for r in res_nums]
                cdr_res_nums.update(res_nums)
            
            # Extract CA atoms from CDR residues
            for residue in chain:
                res_id = residue.id[1]  # Get residue number
                if res_id in cdr_res_nums and is_aa(residue, standard=True):
                    if 'CA' in residue:
                        atoms.append(residue['CA'])
    
    except Exception as e:
        print(f"Error extracting CDR CA atoms: {e}", file=sys.stderr)
    
    return atoms


def compute_rmsd_after_superposition(fixed_ref, fixed_mob, mobile_ref, mobile_mob, verbose=False):
    """
    Superimpose structures using fixed atoms, then compute RMSD of mobile atoms.
    
    Args:
        fixed_ref: List of atoms from reference structure (fixed part)
        fixed_mob: List of atoms from mobile structure (fixed part)
        mobile_ref: List of atoms from reference structure (mobile part)
        mobile_mob: List of atoms from mobile structure (mobile part)
        verbose: Print warnings for atom count mismatches
    
    Returns:
        rmsd: RMSD of mobile part after alignment
    """
    if len(fixed_ref) == 0 or len(fixed_mob) == 0:
        raise ValueError("Fixed atom lists cannot be empty")
    
    if len(mobile_ref) == 0 or len(mobile_mob) == 0:
        raise ValueError("Mobile atom lists cannot be empty")
    
    if len(fixed_ref) != len(fixed_mob):
        min_len = min(len(fixed_ref), len(fixed_mob))
        if verbose:
            print(f"  Warning: Fixed atom count mismatch (ref={len(fixed_ref)}, mob={len(fixed_mob)}), using first {min_len}")
        fixed_ref = fixed_ref[:min_len]
        fixed_mob = fixed_mob[:min_len]
    
    if len(mobile_ref) != len(mobile_mob):
        min_len = min(len(mobile_ref), len(mobile_mob))
        if verbose:
            print(f"  Warning: Mobile atom count mismatch (ref={len(mobile_ref)}, mob={len(mobile_mob)}), using first {min_len}")
        mobile_ref = mobile_ref[:min_len]
        mobile_mob = mobile_mob[:min_len]
    
    # Perform superposition using fixed atoms
    sup = Superimposer()
    sup.set_atoms(fixed_ref, fixed_mob)
    
    # Extract mobile atom coordinates WITHOUT modifying the structure
    # BioPython's Superimposer uses right-multiplication: coord @ R + t
    R = sup.rotran[0].astype('f')  # Rotation matrix (3x3)
    t = sup.rotran[1].astype('f')  # Translation vector (3,)
    
    mobile_ref_coords = np.array([atom.coord for atom in mobile_ref])
    mobile_mob_coords_original = np.array([atom.coord for atom in mobile_mob])
    
    # Apply transformation: coords_transformed = coords @ R + t (right multiplication!)
    mobile_mob_coords_transformed = np.dot(mobile_mob_coords_original, R) + t
    
    # Compute RMSD
    diff = mobile_ref_coords - mobile_mob_coords_transformed
    rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
    
    return rmsd


def parse_seed_sample_from_path(cif_path):
    """Extract seed and sample from path like seed-7570_sample-0/model.cif"""
    import re
    match = re.search(r'seed-(\d+)_sample-(\d+)', cif_path)
    if match:
        return int(match.group(1)), int(match.group(2))
    return None, None


def main():
    if len(sys.argv) != 4:
        print("Usage: python3 compute_rmsd.py <af3_output_dir> <rf2_pdb> <ranking_csv>", file=sys.stderr)
        print("\nExample:", file=sys.stderr)
        print("  python3 compute_rmsd.py \\", file=sys.stderr)
        print("    jobs/pHLA_34_2/outputs/Antibody-pHLA_Complex \\", file=sys.stderr)
        print("    jobs/pHLA_34_2/design_34_dldesign_2_best.pdb \\", file=sys.stderr)
        print("    jobs/pHLA_34_2/outputs/Antibody-pHLA_Complex/ranking_scores.csv", file=sys.stderr)
        sys.exit(1)
    
    af3_output_dir = sys.argv[1]
    rf2_pdb = sys.argv[2]
    ranking_csv = sys.argv[3]
    
    if not os.path.isdir(af3_output_dir):
        print(f"Error: AF3 output directory not found: {af3_output_dir}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.isfile(rf2_pdb):
        print(f"Error: RF2 PDB file not found: {rf2_pdb}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.isfile(ranking_csv):
        print(f"Error: Ranking CSV not found: {ranking_csv}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Reading RF2 PDB: {rf2_pdb}")
    pdb_parser = PDBParser(QUIET=True)
    rf2_structure = pdb_parser.get_structure('rf2', rf2_pdb)
    
    # Parse CDR regions from RF2 PDB
    print("Parsing CDR regions from RF2 PDB...")
    cdr_regions = parse_cdr_regions_from_pdb(rf2_pdb)
    total_cdr_residues = sum(len(cdr_regions[chain][cdr]) for chain in ['H', 'L'] for cdr in [1, 2, 3])
    print(f"  Found CDR residues: H1={len(cdr_regions['H'][1])}, H2={len(cdr_regions['H'][2])}, H3={len(cdr_regions['H'][3])}, "
          f"L1={len(cdr_regions['L'][1])}, L2={len(cdr_regions['L'][2])}, L3={len(cdr_regions['L'][3])} (total: {total_cdr_residues})")
    
    # Extract RF2 reference atoms (once)
    print("Extracting RF2 reference atoms...")
    rf2_mobile = get_ca_atoms_from_chains(rf2_structure, 0, ['H', 'L'])
    rf2_cdr = get_ca_atoms_from_cdr_regions(rf2_structure, 0, cdr_regions, is_af3=False)
    print(f"  RF2 mobile CA atoms (H+L chains): {len(rf2_mobile)}")
    print(f"  RF2 CDR CA atoms: {len(rf2_cdr)}")
    
    if len(rf2_mobile) == 0:
        print("Error: No CA atoms found in RF2 H/L chains (antibody)", file=sys.stderr)
        print("  Available chains:", list(rf2_structure[0].child_dict.keys()), file=sys.stderr)
        sys.exit(1)
    
    if len(rf2_cdr) == 0:
        print("Warning: No CDR CA atoms found in RF2. CDR RMSD will not be calculated.")
    
    # Find all seed-sample CIF files (try multiple patterns)
    cif_pattern1 = os.path.join(af3_output_dir, "seed-*_sample-*/*_model.cif")
    cif_pattern2 = os.path.join(af3_output_dir, "seed-*_sample-*/model.cif")
    cif_files = glob.glob(cif_pattern1)
    if not cif_files:
        cif_files = glob.glob(cif_pattern2)
    
    if not cif_files:
        print(f"Error: No seed-sample CIF files found", file=sys.stderr)
        print(f"  Tried pattern 1: {cif_pattern1}", file=sys.stderr)
        print(f"  Tried pattern 2: {cif_pattern2}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(cif_files)} CIF files to process")
    
    # Compute RMSD for each seed-sample
    rmsd_dict = {}  # key: (seed, sample), value: antibody_rmsd
    cdr_rmsd_dict = {}  # key: (seed, sample), value: cdr_rmsd
    cif_parser = MMCIFParser(QUIET=True)
    
    for cif_file in cif_files:
        seed, sample = parse_seed_sample_from_path(cif_file)
        if seed is None or sample is None:
            print(f"Warning: Could not parse seed/sample from {cif_file}, skipping")
            continue
        
        print(f"Processing seed={seed}, sample={sample}...")
        
        try:
            af3_structure = cif_parser.get_structure('af3', cif_file)
            
            # Check available chains
            available_chains = list(af3_structure[0].child_dict.keys())
            
            # Extract pHLA atoms using proper mapping (RF2 T chain -> AF3 C+A chains)
            rf2_phla, af3_phla = get_pHLA_atoms_for_alignment(rf2_structure, af3_structure, model_id=0)
            
            # Extract antibody atoms
            af3_mobile = get_ca_atoms_from_chains(af3_structure, 0, ['H', 'L'])
            af3_cdr = get_ca_atoms_from_cdr_regions(af3_structure, 0, cdr_regions, is_af3=True)
            
            if len(rf2_phla) == 0 or len(af3_phla) == 0:
                print(f"  Warning: No pHLA CA atoms (RF2={len(rf2_phla)}, AF3={len(af3_phla)}), available: {available_chains}, skipping")
                continue
            
            if len(af3_mobile) == 0:
                print(f"  Warning: No CA atoms in AF3 H/L chains (antibody), available: {available_chains}, skipping")
                continue
            
            print(f"  pHLA CA atoms: RF2={len(rf2_phla)}, AF3={len(af3_phla)}")
            print(f"  Antibody CA atoms: RF2={len(rf2_mobile)}, AF3={len(af3_mobile)}")
            print(f"  CDR CA atoms: RF2={len(rf2_cdr)}, AF3={len(af3_cdr)}")
            
            # Compute antibody RMSD (align on pHLA, measure antibody)
            rmsd = compute_rmsd_after_superposition(
                fixed_ref=rf2_phla,
                fixed_mob=af3_phla,
                mobile_ref=rf2_mobile,
                mobile_mob=af3_mobile,
                verbose=True
            )
            
            rmsd_dict[(seed, sample)] = rmsd
            print(f"  Antibody RMSD = {rmsd:.4f} Å")
            
            # Compute CDR RMSD (if CDR atoms available)
            if len(rf2_cdr) > 0 and len(af3_cdr) > 0:
                # Check for significant atom count mismatch
                atom_ratio = len(af3_cdr) / len(rf2_cdr) if len(rf2_cdr) > 0 else 0
                if atom_ratio < 0.8 or atom_ratio > 1.2:
                    print(f"  Warning: Large CDR atom count mismatch (RF2={len(rf2_cdr)}, AF3={len(af3_cdr)}, ratio={atom_ratio:.2f})")
                    print(f"  CDR RMSD may be unreliable - skipping")
                    cdr_rmsd_dict[(seed, sample)] = None  # Mark as unreliable
                else:
                    cdr_rmsd = compute_rmsd_after_superposition(
                        fixed_ref=rf2_phla,
                        fixed_mob=af3_phla,
                        mobile_ref=rf2_cdr,
                        mobile_mob=af3_cdr,
                        verbose=True  # Show warnings for atom count adjustments
                    )
                    # Sanity check: CDR RMSD should not be much higher than antibody RMSD
                    if cdr_rmsd > rmsd * 1.5:
                        print(f"  Warning: CDR RMSD ({cdr_rmsd:.2f}) >> Antibody RMSD ({rmsd:.2f}) - may indicate alignment issue")
                    cdr_rmsd_dict[(seed, sample)] = cdr_rmsd
                    print(f"  CDR RMSD = {cdr_rmsd:.4f} Å")
            else:
                print(f"  CDR RMSD = N/A (CDR atoms: RF2={len(rf2_cdr)}, AF3={len(af3_cdr)})")
            
        except Exception as e:
            print(f"  Error processing {cif_file}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    if not rmsd_dict:
        print("Error: No RMSD values computed successfully", file=sys.stderr)
        sys.exit(1)
    
    # Update CSV
    print(f"\nUpdating {ranking_csv}...")
    rows = []
    with open(ranking_csv, 'r') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        
        if not fieldnames:
            print("Error: CSV file has no header", file=sys.stderr)
            sys.exit(1)
        
        # Check for required columns
        if 'seed' not in fieldnames or 'sample' not in fieldnames:
            print(f"Error: CSV missing required columns. Found: {fieldnames}", file=sys.stderr)
            sys.exit(1)
        
        # Add new columns if not present
        new_columns = []
        if 'antibody_rmsd' not in fieldnames:
            new_columns.append('antibody_rmsd')
        if 'cdr_rmsd' not in fieldnames:
            new_columns.append('cdr_rmsd')
        
        if new_columns:
            fieldnames = list(fieldnames) + new_columns
        
        for row in reader:
            try:
                seed = int(row['seed'])
                sample = int(row['sample'])
                
                # Update antibody RMSD
                if (seed, sample) in rmsd_dict:
                    row['antibody_rmsd'] = f"{rmsd_dict[(seed, sample)]:.4f}"
                else:
                    row['antibody_rmsd'] = "N/A"
                
                # Update CDR RMSD
                if (seed, sample) in cdr_rmsd_dict:
                    cdr_val = cdr_rmsd_dict[(seed, sample)]
                    if cdr_val is None:
                        row['cdr_rmsd'] = "N/A"  # Unreliable due to atom mismatch
                    else:
                        row['cdr_rmsd'] = f"{cdr_val:.4f}"
                else:
                    row['cdr_rmsd'] = "N/A"
                    
            except (KeyError, ValueError) as e:
                print(f"Warning: Error processing row {row}: {e}")
                row['antibody_rmsd'] = "ERROR"
                row['cdr_rmsd'] = "ERROR"
            
            rows.append(row)
    
    # Write updated CSV
    with open(ranking_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"Successfully updated {len(rows)} rows with RMSD values")
    print(f"Antibody RMSD range: {min(rmsd_dict.values()):.4f} - {max(rmsd_dict.values()):.4f} Å")
    if cdr_rmsd_dict:
        valid_cdr = [v for v in cdr_rmsd_dict.values() if v is not None]
        if valid_cdr:
            print(f"CDR RMSD range: {min(valid_cdr):.4f} - {max(valid_cdr):.4f} Å")
            print(f"Valid CDR RMSD: {len(valid_cdr)}/{len(cdr_rmsd_dict)}")


if __name__ == '__main__':
    main()
