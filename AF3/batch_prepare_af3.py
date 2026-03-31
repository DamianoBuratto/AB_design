#!/usr/bin/env python3
"""
Batch prepare AF3 jobs from RF2 best pdbs.

Usage (from AF3/):
  python3 batch_prepare_af3.py --input-dir ../test2_outputs/rf2 --output-base ../test2_outputs/af3_jobs --run-alphafold --sbatch

"""

import os
import sys
import argparse
import re
import shutil
import json
from collections import OrderedDict

# Three-letter to one-letter map
THREE_TO_ONE = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I',
    'LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S',
    'THR':'T','VAL':'V','TRP':'W','TYR':'Y','SEC':'U','PYL':'O'
}

def parse_args():
    p = argparse.ArgumentParser(description='Prepare AF3 jobs from RF2 best pdbs (run from AF3/)')
    p.add_argument('--input-dir', default='RFantibody_best', help='Directory with *_best.pdb files (default: AF3/RFantibody_best)')
    p.add_argument('--templates-dir', default='templates', help='Directory containing pHLA_A.fasta and pHLA_C.fasta (default: templates/)')
    p.add_argument('--output-base', default='jobs', help='Base dir to create pHLA_* job folders (default: AF3/jobs)')
    p.add_argument('--chains', default='H,L', help='Comma-separated chain IDs to keep (default H,L)')
    p.add_argument('--build-script', default=None, help='Path to build_af3_json.py to use (optional)')
    p.add_argument('--alphafold-script', default='./alphafold3.sh', help='Path to alphafold3.sh relative to AF3/ (default: ./alphafold3.sh)')
    p.add_argument('--run-alphafold', action='store_true', help='Run alphafold3.sh after job creation')
    p.add_argument('--sbatch', action='store_true', help='When running alphafold, submit with sbatch instead of bash')
    p.add_argument('--mut-only', action='store_true', help='Only create mutant jobs, skip wild-type (WT) jobs')
    p.add_argument('--dry-run', action='store_true', help='Only print actions without making changes')
    return p.parse_args()

def extract_design_seq(filename):
    # match design_166_dldesign_1_best.pdb
    m = re.search(r'design_(\d+)_dldesign_(\d+)', filename)
    if m:
        return m.group(1), m.group(2)
    # fallback: use numeric groups or basename
    base = os.path.splitext(os.path.basename(filename))[0]
    return base, '0'

def extract_chains(pdb_path, keep_chains):
    """Return dict chain_id -> list of ATOM/HETATM lines for that chain"""
    chains = {c: [] for c in keep_chains}
    other_lines = []
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith(('ATOM','HETATM')):
                    chain = line[21].strip() or ' '
                    if chain in keep_chains:
                        chains[chain].append(line.rstrip('\n'))
                else:
                    other_lines.append(line.rstrip('\n'))
    except Exception:
        raise
    return chains, other_lines

def write_antibody_pdb(jobdir, outname, chains_dict, other_lines):
    path = os.path.join(jobdir, outname)
    with open(path, 'w') as out:
        # write non-ATOM header lines that look useful
        for l in other_lines:
            if l.startswith(('HEADER','REMARK','TITLE','MODEL')):
                out.write(l + '\n')
        # write atoms for chains
        for chain in sorted(chains_dict.keys()):
            for l in chains_dict[chain]:
                out.write(l + '\n')
            out.write('TER\n')
        out.write('END\n')
    return path

def pdb_chain_to_sequence(chain_lines):
    # collect residues in order using (resseq, iCode) to avoid duplicates
    residues = OrderedDict()
    for l in chain_lines:
        # columns: resname 17-20, chain 21, resseq 22-26, iCode 26, atom name cols
        resname = l[17:20].strip()
        resseq = l[22:26].strip()
        icode = l[26].strip()
        key = (resseq, icode)
        if key not in residues:
            residues[key] = resname
    # convert to 1-letter (unknown->X)
    seq = ''.join(THREE_TO_ONE.get(r.upper(),'X') for r in residues.values())
    return seq

def write_fasta(path, header, seq):
    with open(path, 'w') as f:
        f.write(f'>{header}\n')
        # wrap 80
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')

def find_build_script(repo_af3_dir):
    # search AF3/build_af3_json.py or AF3/*/build_af3_json.py
    import glob
    candidate = os.path.join(repo_af3_dir, 'build_af3_json.py')
    if os.path.isfile(candidate):
        return candidate
    for p in glob.glob(os.path.join(repo_af3_dir, '*', 'build_af3_json.py')):
        return p
    return None

def main():
    args = parse_args()
    # When running from AF3/, most paths should be specified relative to cwd
    in_dir = args.input_dir
    tpl_dir = args.templates_dir
    out_base = args.output_base
    keep_chains = [c.strip() for c in args.chains.split(',') if c.strip()]

    if not os.path.isdir(in_dir):
        print('Input dir not found:', in_dir, file=sys.stderr); sys.exit(1)
    if not os.path.isdir(tpl_dir):
        print('Templates dir not found:', tpl_dir, file=sys.stderr); sys.exit(1)

    a_template = os.path.join(tpl_dir, 'pHLA_A.fasta')
    c_template = os.path.join(tpl_dir, 'pHLA_C.fasta')
    if not os.path.isfile(a_template) or not os.path.isfile(c_template):
        print('Template files pHLA_A.fasta and/or pHLA_C.fasta missing in', tpl_dir, file=sys.stderr); sys.exit(1)

    os.makedirs(out_base, exist_ok=True)

    # repo root is parent of AF3/ when cwd=AF3
    repo_root = os.path.abspath(os.path.join(os.getcwd(), '..'))
    repo_af3 = os.path.join(repo_root, 'AF3')
    build_script = args.build_script or find_build_script(repo_af3)
    if build_script:
        print('Using build script:', build_script)
    else:
        print('No build_af3_json.py found under AF3/; will create minimal JSON files. You may want to supply --build-script', file=sys.stderr)

    pdbs = sorted([os.path.join(in_dir, f) for f in os.listdir(in_dir) if f.endswith('_best.pdb')])
    if not pdbs:
        print('No *_best.pdb files found in', in_dir); sys.exit(1)

    for pdbfile in pdbs:
        fname = os.path.basename(pdbfile)
        design, seq = extract_design_seq(fname)
        jobname = f'pHLA_{design}_{seq}'
        jobdir = os.path.join(out_base, jobname)
        print('\nPreparing job', jobname)
        if args.dry_run:
            print('DRY RUN: would create', jobdir)
        else:
            os.makedirs(jobdir, exist_ok=True)

        # copy original pdb for record
        if not args.dry_run:
            shutil.copy2(pdbfile, os.path.join(jobdir, fname))

        # extract chains
        chains_dict, other_lines = extract_chains(pdbfile, keep_chains)
        antibody_pdb_name = f'antibody_{design}_{seq}.pdb'
        antibody_pdb_path = os.path.join(jobdir, antibody_pdb_name)
        if not args.dry_run:
            write_antibody_pdb(jobdir, antibody_pdb_name, chains_dict, other_lines)
            print('Wrote', antibody_pdb_path)
        else:
            print('DRY RUN: would write', antibody_pdb_path)

        # write per-chain fasta
        for chain in keep_chains:
            seqres = pdb_chain_to_sequence(chains_dict.get(chain, []))
            if not seqres:
                print(f'Warning: no residues found for chain {chain} in {fname}')
            header = f'antibody_{design}_{seq}_{chain}'
            fasta_name = f'antibody_{chain}.fasta'  # antibody_H.fasta, antibody_L.fasta
            fasta_path = os.path.join(jobdir, fasta_name)
            if not args.dry_run:
                write_fasta(fasta_path, header, seqres)
                print('Wrote', fasta_path)
            else:
                print('DRY RUN: would write', fasta_path)

        # copy templates
        if not args.dry_run:
            shutil.copy2(a_template, os.path.join(jobdir, 'pHLA_A.fasta'))
            shutil.copy2(c_template, os.path.join(jobdir, 'pHLA_C.fasta'))
            print('Copied pHLA_A.fasta and pHLA_C.fasta to', jobdir)
        else:
            print('DRY RUN: would copy templates to', jobdir)

        # attempt to create JSON via build script
        json_out = os.path.join(jobdir, f'{jobname}.json')
        if build_script and not args.dry_run:
            try:
                # copy build script into jobdir and run it there
                bs_dest = os.path.join(jobdir, os.path.basename(build_script))
                shutil.copy2(build_script, bs_dest)
                print('Copied build script to', bs_dest)
                # run it in jobdir
                rc = os.system(f'cd "{jobdir}" && python3 "{os.path.basename(build_script)}"')
                if rc != 0:
                    print('Warning: build script returned non-zero code; creating minimal JSON instead')
                    raise RuntimeError('build script failed')
                # try to find produced json file in jobdir
                candidates = [f for f in os.listdir(jobdir) if f.endswith('.json')]
                if candidates:
                    # prefer one named same as job
                    for c in candidates:
                        if jobname in c:
                            json_out = os.path.join(jobdir, c); break
                    else:
                        json_out = os.path.join(jobdir, candidates[0])
                print('JSON created (by build script):', json_out)
            except Exception as e:
                print('Build script failed or not compatible:', e)
                # fall through to minimal JSON

        if not os.path.isfile(json_out) and not args.dry_run:
            # create minimal JSON with paths (user should check compatibility)
            minimal = {
                'antibody_h_fasta': os.path.abspath(os.path.join(jobdir, 'antibody_H.fasta')),
                'antibody_l_fasta': os.path.abspath(os.path.join(jobdir, 'antibody_L.fasta')),
                'pHLA_A_fasta': os.path.abspath(os.path.join(jobdir, 'pHLA_A.fasta')),
                'pHLA_C_fasta': os.path.abspath(os.path.join(jobdir, 'pHLA_C.fasta')),
                'output_dir': os.path.abspath(os.path.join(jobdir, 'outputs'))
            }
            with open(json_out, 'w') as jf:
                json.dump(minimal, jf, indent=2)
            print('Wrote minimal JSON to', json_out)
        elif args.dry_run:
            print('DRY RUN: would write minimal JSON to', json_out)

        # optionally run alphafold on this job
        if args.run_alphafold and not args.dry_run:
            af_script = os.path.abspath(args.alphafold_script)
            if not os.path.isfile(af_script):
                print('alphafold script not found at', af_script, 'skip running')
            else:
                cmd = f'"{af_script}" "{json_out}"'
                # ensure logs dir exists at AF3/logs
                logs_dir = os.path.join(os.getcwd(), 'logs')
                if args.sbatch:
                    full = f'sbatch --output="{logs_dir}/{jobname}.out" --error="{logs_dir}/{jobname}.err" {af_script} "{json_out}"'
                    print('Submitting with sbatch:', full)
                    # create logs dir before submitting
                    os.makedirs(logs_dir, exist_ok=True)
                    os.system(full)
                else:
                    print('Running alphafold script (blocking):', cmd)
                    os.system(f'bash {af_script} "{json_out}"')
        elif args.run_alphafold and args.dry_run:
            print('DRY RUN: would run alphafold script on', json_out)

        # --- create a parallel _wt jobdir with modified JSON where protein C sequence H->R ---
        # Skip WT creation if --mut-only flag is set
        if args.mut_only:
            if args.dry_run:
                print('DRY RUN: --mut-only flag set, skipping WT job creation')
            else:
                print('--mut-only flag set, skipping WT job creation')
            continue
        
        # create wt jobdir name
        wt_jobname = f'{jobname}_wt'
        wt_jobdir = os.path.join(out_base, wt_jobname)
        if args.dry_run:
            print('DRY RUN: would create WT jobdir', wt_jobdir)
        else:
            # copy jobdir into wt_jobdir (if not existing)
            if not os.path.isdir(jobdir):
                print('Warning: original jobdir not found to create WT:', jobdir)
            else:
                if os.path.isdir(wt_jobdir):
                    print('WT jobdir already exists, will overwrite JSON only:', wt_jobdir)
                else:
                    shutil.copytree(jobdir, wt_jobdir)

        # locate json in wt_jobdir and modify it
        wt_json = os.path.join(wt_jobdir, f'{wt_jobname}.json')
        # if not present, try the non-wt json name in wt_jobdir
        if not os.path.isfile(wt_json):
            # look for any json in wt_jobdir
            candidates = [f for f in os.listdir(wt_jobdir) if f.endswith('.json')] if os.path.isdir(wt_jobdir) else []
            if candidates:
                wt_json = os.path.join(wt_jobdir, candidates[0])

        if args.dry_run:
            print('DRY RUN: would modify WT JSON at', wt_json)
        else:
            if os.path.isfile(wt_json):
                try:
                    with open(wt_json, 'r') as jf:
                        data = json.load(jf)
                    # find protein entry with id == 'C'
                    if isinstance(data, dict):
                        # support either single protein or list under 'proteins' or 'protein'
                        modified = False
                        # try common key names
                        if 'protein' in data and isinstance(data['protein'], dict):
                            p = data['protein']
                            if p.get('id') == 'C' and isinstance(p.get('sequence'), str):
                                p['sequence'] = p['sequence'].replace('H', 'R')
                                modified = True
                        # also try a 'proteins' list
                        if 'proteins' in data and isinstance(data['proteins'], list):
                            for p in data['proteins']:
                                if p.get('id') == 'C' and isinstance(p.get('sequence'), str):
                                    p['sequence'] = p['sequence'].replace('H', 'R')
                                    modified = True
                        # also handle AF3-style 'sequences' list where each item may be {'protein': {...}}
                        if 'sequences' in data and isinstance(data['sequences'], list):
                            for item in data['sequences']:
                                if isinstance(item, dict) and 'protein' in item and isinstance(item['protein'], dict):
                                    p = item['protein']
                                    if p.get('id') == 'C' and isinstance(p.get('sequence'), str):
                                        p['sequence'] = p['sequence'].replace('H', 'R')
                                        modified = True
                        # fallback: try top-level keys that look like protein entries
                        # write back only if modified
                        if modified:
                            wt_json_out = os.path.join(wt_jobdir, f'{jobname}_wt.json')
                            with open(wt_json_out, 'w') as jf:
                                json.dump(data, jf, indent=2)
                            print('Wrote modified WT JSON to', wt_json_out)
                        else:
                            print('Warning: did not find modifiable protein C entry in', wt_json)
                except Exception as e:
                    print('Failed to modify WT JSON:', e)
            else:
                print('No json found in', wt_jobdir, 'to modify for WT')

        # Optionally run alphafold on WT jobdir json
        if args.run_alphafold and not args.dry_run:
            # find the new wt json we just wrote
            candidate_json = None
            if os.path.isfile(os.path.join(wt_jobdir, f'{jobname}_wt.json')):
                candidate_json = os.path.join(wt_jobdir, f'{jobname}_wt.json')
            else:
                # fallback any json in wt_jobdir
                if os.path.isdir(wt_jobdir):
                    cand = [os.path.join(wt_jobdir, f) for f in os.listdir(wt_jobdir) if f.endswith('.json')]
                    if cand:
                        candidate_json = cand[0]
            if candidate_json:
                af_script = os.path.abspath(args.alphafold_script)
                if not os.path.isfile(af_script):
                    print('alphafold script not found at', af_script, 'skip running WT')
                else:
                    logs_dir = os.path.join(os.getcwd(), 'logs')
                    if args.sbatch:
                        full = f'sbatch --output="{logs_dir}/{wt_jobname}.out" --error="{logs_dir}/{wt_jobname}.err" {af_script} "{candidate_json}"'
                        print('Submitting WT with sbatch:', full)
                        os.makedirs(logs_dir, exist_ok=True)
                        os.system(full)
                    else:
                        print('Running alphafold script for WT (blocking):', candidate_json)
                        os.system(f'bash {af_script} "{candidate_json}"')
            else:
                print('No WT json available to run alphafold for', wt_jobdir)

    print('\nAll jobs prepared under', os.path.abspath(out_base))

if __name__ == '__main__':
    main()
