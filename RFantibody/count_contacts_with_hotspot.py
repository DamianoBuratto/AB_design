import os
import glob
import math
import csv
import argparse

DEFAULT_CUTOFFS = [5, 6, 7, 8]
BACKBONE_EXCLUDE = {'N', 'C', 'O', 'OXT'}
BACKBONE_INCLUDE = {'N', 'CA', 'C', 'O', 'OXT'}

def parse_pdb_structure(pdb_file):
    """Parse PDB structure by chain, storing residue order and atom coordinates."""
    structure = {}
    with open(pdb_file) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue  # Skip non-primary conformations
            resname = line[17:20].strip()
            chain = line[21].strip()
            resseq = int(line[22:26])
            icode = line[26].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip()

            chain_records = structure.setdefault(chain, [])
            if chain_records and chain_records[-1]['resseq'] == resseq and chain_records[-1]['icode'] == icode:
                residue = chain_records[-1]
            else:
                residue = {'resname': resname, 'resseq': resseq, 'icode': icode, 'atoms': []}
                chain_records.append(residue)
            residue['atoms'].append({'name': atom_name, 'coord': (x, y, z), 'element': element})
    return structure

def calc_dist(coord1, coord2):
    return math.sqrt(sum((a-b)**2 for a, b in zip(coord1, coord2)))

def collect_hotspot_atoms(residue, mode):
    """Collect HIS atom coordinates based on mode."""
    coords = []
    names = []
    added = set()

    def add_atom(atom_name, coord):
        if atom_name in added:
            return
        names.append(atom_name)
        coords.append(coord)
        added.add(atom_name)

    for atom in residue['atoms']:
        name = atom['name']
        element = atom.get('element', '')
        clean_name = name.strip().lstrip('0123456789')

        if mode == 'ca_only':
            if name == 'CA':
                add_atom(name, atom['coord'])
            continue

        if mode == 'backbone':
            if name in BACKBONE_INCLUDE:
                add_atom(name, atom['coord'])
            continue

        # Sidechain mode: include ALL heavy atoms (backbone + sidechain)
        # Skip only hydrogen atoms
        if element.upper() == 'H' or clean_name.startswith(('H', 'D')):
            continue
        add_atom(name, atom['coord'])

    if not coords:
        for atom in residue['atoms']:
            if atom['name'] in BACKBONE_INCLUDE:
                add_atom(atom['name'], atom['coord'])

    return names, coords


def gather_antibody_ca(structure):
    residues = []
    for chain_id, chain_records in structure.items():
        if chain_id == 'T':
            continue
        for residue in chain_records:
            for atom in residue['atoms']:
                if atom['name'] == 'CA':
                    residues.append({'chain': chain_id, 'resseq': residue['resseq'], 'resname': residue['resname'], 'coord': atom['coord']})
                    break
    return residues


def count_contacts_multi_cutoffs(pdb_file, cutoffs, mode='sidechain', structure=None, verbose=True):
    structure = structure or parse_pdb_structure(pdb_file)
    cutoffs = sorted({int(round(c)) for c in cutoffs}) if cutoffs else DEFAULT_CUTOFFS
    if not cutoffs:
        cutoffs = DEFAULT_CUTOFFS

    t_chain = structure.get('T', [])
    if not t_chain:
        if verbose:
            print(f"{os.path.basename(pdb_file)}: T chain not found")
        return None, {c: [] for c in cutoffs}

    antibody_ca = gather_antibody_ca(structure)
    if not antibody_ca:
        if verbose:
            print(f"{os.path.basename(pdb_file)}: Antibody CA atoms not found")
        return None, {c: [] for c in cutoffs}

    # Find the target HIS in the RHC motif (R-H-C pattern)
    hotspot_res = None
    prev_name = None
    next_name = None
    t_index = None
    
    for i in range(1, len(t_chain) - 1):
        if (t_chain[i]['resname'] == 'HIS' and 
            t_chain[i-1]['resname'] == 'ARG' and 
            t_chain[i+1]['resname'] == 'CYS'):
            hotspot_res = t_chain[i]
            t_index = i
            prev_name = t_chain[i-1]['resname']
            next_name = t_chain[i+1]['resname']
            break
    
    if hotspot_res is None:
        if verbose:
            print(f"{os.path.basename(pdb_file)}: HIS in RHC motif not found (requires ARG-HIS-CYS)")
        return None, {c: [] for c in cutoffs}
    
    context_msg = f"prev={prev_name}, next={next_name}" if verbose else ""
    atom_names, hotspot_coords = collect_hotspot_atoms(hotspot_res, mode)
    if not hotspot_coords:
        if verbose:
            print(f"{os.path.basename(pdb_file)}: HIS missing usable atom coordinates")
        return None, {c: [] for c in cutoffs}
    if verbose and mode == 'sidechain':
        has_sidechain = any(name not in BACKBONE_INCLUDE for name in atom_names)
        if not has_sidechain:
            print("HIS missing sidechain atoms, fallback to backbone heavy atoms.")

    contacts_by_cutoff = {c: [] for c in cutoffs}
    min_distance = float('inf')

    for residue in antibody_ca:
        min_dist = min(calc_dist(coord, residue['coord']) for coord in hotspot_coords)
        if min_dist < min_distance:
            min_distance = min_dist
        for c in cutoffs:
            if min_dist <= c:
                contacts_by_cutoff[c].append((residue['chain'], residue['resseq'], residue['resname']))

    if verbose:
        ca_coord = next((atom['coord'] for atom in hotspot_res['atoms'] if atom['name'] == 'CA'), None)
        if ca_coord:
            print(
                f"{os.path.basename(pdb_file)} Found HIS in RHC: resseq={hotspot_res['resseq']} CA={ca_coord}"
            )
        if context_msg:
            print(f"HIS neighbors: {context_msg} (RHC motif)")
        print(f"HIS atoms used: {', '.join(atom_names)}")
        if math.isfinite(min_distance):
            print(f"Min HIS-antibody CA distance: {min_distance:.2f} Å")

    hotspot_info = {
        'chain': 'T',
        'resname': hotspot_res['resname'],
        'resseq': hotspot_res['resseq'],
        'atoms': atom_names,
    }
    return hotspot_info, contacts_by_cutoff

def main():
    parser = argparse.ArgumentParser(description='Count contacts between HIS in RHC motif (T-chain) and antibody CA atoms')
    parser.add_argument('--workdir', '-w', default=os.getcwd(),
                        help='Working directory for output and default PDB search')
    parser.add_argument('--pdb-dir', '-p', default=None,
                        help='Directory containing PDB files (overrides default)')
    parser.add_argument('--mode', choices=['sidechain', 'backbone', 'ca_only'], default='sidechain',
                        help='sidechain: HIS CA+non-H sidechain (fallback to backbone); backbone: N/CA/C/O only; ca_only: CA only')
    parser.add_argument('--cutoffs', '-c', type=float, nargs='+', default=None,
                        help='Contact distance cutoffs in Angstroms (default: 5 6 7 8)')
    parser.add_argument('--primary-cutoff', type=float, default=None,
                        help='Cutoff for residue detail column (default: minimum cutoff)')
    args = parser.parse_args()

    cutoffs = sorted({int(round(c)) for c in (args.cutoffs or DEFAULT_CUTOFFS)})
    if not cutoffs:
        cutoffs = DEFAULT_CUTOFFS.copy()
    primary_cutoff = int(round(args.primary_cutoff)) if args.primary_cutoff is not None else cutoffs[0]
    if primary_cutoff not in cutoffs:
        cutoffs.append(primary_cutoff)
        cutoffs = sorted(set(cutoffs))

    if args.pdb_dir:
        pdb_dir = args.pdb_dir
    else:
        pdb_dir = os.path.join(args.workdir, 'test2_outputs', 'rfdiffusion')

    if not os.path.isdir(pdb_dir):
        print(f"PDB directory does not exist: {pdb_dir}")
        return

    pdb_files = glob.glob(os.path.join(pdb_dir, '*.pdb'))
    results = []
    for pdb_file in pdb_files:
        hotspot, contacts_by_cutoff = count_contacts_multi_cutoffs(pdb_file, cutoffs, mode=args.mode)
        row = {'PDB File': os.path.basename(pdb_file)}
        
        # Use rf2_analysis.csv naming convention
        for c in cutoffs:
            if args.mode == 'sidechain':
                label = f'sc_<={int(c)}A'
            elif args.mode == 'backbone':
                label = f'bb_<={int(c)}A'
            else:
                label = f'(<= {c}A)'
            row[label] = len(set(contacts_by_cutoff[c])) if hotspot is not None else 0
        results.append(row)

    csv_path = os.path.join(args.workdir, 'contacts_with_hotspot_summary.csv')
    try:
        os.makedirs(args.workdir, exist_ok=True)
    except Exception as e:
        print(f"Failed to create workdir {args.workdir}: {e}")
        return

    # Column names match rf2_analysis.csv format
    if args.mode == 'sidechain':
        fieldnames = ['PDB File'] + [f'sc_<={int(c)}A' for c in cutoffs]
    elif args.mode == 'backbone':
        fieldnames = ['PDB File'] + [f'bb_<={int(c)}A' for c in cutoffs]
    else:
        fieldnames = ['PDB File'] + [f'(<= {c}A)' for c in cutoffs]
    
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"Results written to: {csv_path}\n")
    header = f"{'PDB File':40} | " + ' | '.join([f'{fn}'.ljust(10) for fn in fieldnames[1:]])
    print(header)
    print('-'*len(header))
    for row in results:
        counts = ' | '.join([str(row[fn]).ljust(10) for fn in fieldnames[1:]])
        print(f"{row['PDB File']:40} | {counts}")

if __name__ == '__main__':
    main()