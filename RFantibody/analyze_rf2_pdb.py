#!/usr/bin/env python3
"""
Analyze RF2 scores directly from SCORE lines in RFantibody RF2 output PDB files.
Adds contact statistics using HIS sidechain heavy atoms and backbone-only heavy atoms.
"""

from count_contacts_with_hotspot import (
    DEFAULT_CUTOFFS,
    count_contacts_multi_cutoffs,
    parse_pdb_structure,
)

import argparse
import glob
import os
import re

PRIMARY_CUTOFF = 6


def extract_scores_from_pdb_score_lines(pdb_file):
    """Extract scores from SCORE lines and CDR lengths from REMARK lines in PDB file."""
    scores = {
        'file': os.path.basename(pdb_file),
        'mean_plddt': None,
        'interaction_pae': None,
        'mean_pae': None,
        'target_aligned_antibody_rmsd': None,
        'target_aligned_cdr_rmsd': None,
        'framework_aligned_antibody_rmsd': None,
        'framework_aligned_cdr_rmsd': None,
        'framework_aligned_H1_rmsd': None,
        'framework_aligned_H2_rmsd': None,
        'framework_aligned_H3_rmsd': None,
        'framework_aligned_L1_rmsd': None,
        'framework_aligned_L2_rmsd': None,
        'framework_aligned_L3_rmsd': None,
        'H1_length': None,
        'H2_length': None,
        'H3_length': None,
        'L1_length': None,
        'L2_length': None,
        'L3_length': None,
    }
    
    # CDR region tracking
    cdr_regions = {'H1': [], 'H2': [], 'H3': [], 'L1': [], 'L2': [], 'L3': []}

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                # Extract CDR lengths from REMARK PDBinfo-LABEL lines
                if line.startswith('REMARK PDBinfo-LABEL:'):
                    match = re.match(r'REMARK PDBinfo-LABEL:\s+(\d+)\s+(H[123]|L[123])', line)
                    if match:
                        resnum = int(match.group(1))
                        cdr_name = match.group(2)
                        if cdr_name in cdr_regions:
                            cdr_regions[cdr_name].append(resnum)
                
                # Extract SCORE lines
                if line.startswith('SCORE'):
                    # Use simpler pattern: SCORE key: value
                    match = re.match(r'SCORE\s+(\S+):\s+([0-9.-]+)', line)
                    if match:
                        key = match.group(1)
                        value = float(match.group(2))
                        
                        if key == 'pred_lddt':
                            scores['mean_plddt'] = value * 100 if value <= 1.0 else value
                        elif key == 'interaction_pae':
                            scores['interaction_pae'] = value
                        elif key == 'pae':
                            scores['mean_pae'] = value
                        elif key == 'target_aligned_antibody_rmsd':
                            scores['target_aligned_antibody_rmsd'] = value
                        elif key == 'target_aligned_cdr_rmsd':
                            scores['target_aligned_cdr_rmsd'] = value
                        elif key == 'framework_aligned_antibody_rmsd':
                            scores['framework_aligned_antibody_rmsd'] = value
                        elif key == 'framework_aligned_cdr_rmsd':
                            scores['framework_aligned_cdr_rmsd'] = value
                        elif key == 'framework_aligned_H1_rmsd':
                            scores['framework_aligned_H1_rmsd'] = value
                        elif key == 'framework_aligned_H2_rmsd':
                            scores['framework_aligned_H2_rmsd'] = value
                        elif key == 'framework_aligned_H3_rmsd':
                            scores['framework_aligned_H3_rmsd'] = value
                        elif key == 'framework_aligned_L1_rmsd':
                            scores['framework_aligned_L1_rmsd'] = value
                        elif key == 'framework_aligned_L2_rmsd':
                            scores['framework_aligned_L2_rmsd'] = value
                        elif key == 'framework_aligned_L3_rmsd':
                            scores['framework_aligned_L3_rmsd'] = value
        
        # Calculate CDR lengths
        for cdr_name, residues in cdr_regions.items():
            if residues:
                scores[f'{cdr_name}_length'] = len(residues)

    except Exception as exc:
        print(f"Warning: Error reading PDB file {pdb_file}: {exc}")
        return None

    return scores





def analyze_rf2_pdb_results(rf2_dir):
    """Analyze PDB files in RF2 directory."""
    print("=== RF2 PDB Results Analysis ===")
    print(f"Analysis directory: {rf2_dir}")

    pdb_files = glob.glob(os.path.join(rf2_dir, "*_best.pdb"))

    if not pdb_files:
        print("Error: No *_best.pdb files found")
        return None

    results = []
    cutoffs = sorted(set(DEFAULT_CUTOFFS + [PRIMARY_CUTOFF]))
    for pdb_file in pdb_files:
        scores = extract_scores_from_pdb_score_lines(pdb_file)
        if not scores:
            continue

        basename = os.path.basename(pdb_file).replace("_best.pdb", "")
        match = re.match(r'design_(\d+)_dldesign_(\d+)', basename)
        if match:
            scores['design_id'] = int(match.group(1))
            scores['sequence_id'] = int(match.group(2))
            scores['design_name'] = f"design_{match.group(1)}"
            scores['sequence_name'] = f"dldesign_{match.group(2)}"

        structure = parse_pdb_structure(pdb_file)
        hotspot_sidechain, contacts_sidechain = count_contacts_multi_cutoffs(
            pdb_file, cutoffs, mode='sidechain', structure=structure, verbose=False
        )
        hotspot_backbone, contacts_backbone = count_contacts_multi_cutoffs(
            pdb_file, cutoffs, mode='backbone', structure=structure, verbose=False
        )

        for c in cutoffs:
            key_side = f'sc_<={c}A'
            key_bb = f'bb_<={c}A'
            if hotspot_sidechain is not None:
                side_contacts = set(contacts_sidechain.get(c, []))
                scores[key_side] = str(len(side_contacts))
            else:
                scores[key_side] = ''

            if hotspot_backbone is not None:
                bb_contacts = set(contacts_backbone.get(c, []))
                scores[key_bb] = str(len(bb_contacts))
            else:
                scores[key_bb] = ''

        sidechain_detail_key = f'sc_hits_<={PRIMARY_CUTOFF}A'
        backbone_detail_key = f'bb_hits_<={PRIMARY_CUTOFF}A'

        if hotspot_sidechain is not None:
            detail_contacts = set(contacts_sidechain.get(PRIMARY_CUTOFF, []))
            scores[sidechain_detail_key] = '; '.join(
                f"{resname} {chain} {resseq}" for chain, resseq, resname in sorted(detail_contacts)
            )
        else:
            scores[sidechain_detail_key] = ''

        if hotspot_backbone is not None:
            detail_bb = set(contacts_backbone.get(PRIMARY_CUTOFF, []))
            scores[backbone_detail_key] = '; '.join(
                f"{resname} {chain} {resseq}" for chain, resseq, resname in sorted(detail_bb)
            )
        else:
            scores[backbone_detail_key] = ''

        results.append(scores)

    return results

def print_analysis_results(results):
    """Print analysis results"""
    if not results:
        print("No valid results found")
        return
    
    print(f"\n=== Analysis Results ({len(results)} designs) ===")
    
    # Statistics
    valid_plddt = [r['mean_plddt'] for r in results if r['mean_plddt'] is not None]
    valid_int_pae = [r['interaction_pae'] for r in results if r['interaction_pae'] is not None]
    
    if valid_plddt:
        print(f"\npLDDT Statistics:")
        print(f"  Average: {sum(valid_plddt)/len(valid_plddt):.1f}")
        print(f"  Range: {min(valid_plddt):.1f} - {max(valid_plddt):.1f}")
    
    if valid_int_pae:
        print(f"\nInteraction PAE Statistics:")
        print(f"  Average: {sum(valid_int_pae)/len(valid_int_pae):.1f} Å")
        print(f"  Range: {min(valid_int_pae):.1f} - {max(valid_int_pae):.1f} Å")
    
    # Top 10 results (original order)
    print(f"\n=== First 10 Designs ===")
    print(f"{'#':<3} {'File':<25} {'pLDDT':<6} {'Int_PAE':<8}")
    print("-" * 44)
    
    for i, result in enumerate(results[:10], 1):
        filename = result['file'].replace('_best.pdb', '')
        plddt = f"{result['mean_plddt']:.1f}" if result['mean_plddt'] is not None else "N/A"
        int_pae = f"{result['interaction_pae']:.1f}" if result['interaction_pae'] is not None else "N/A"
        
        print(f"{i:<3} {filename:<25} {plddt:<6} {int_pae:<8}")
    
    # Group by backbone design
    print(f"\n=== Grouped by Backbone Design ===")
    backbone_groups = {}
    for result in results:
        if 'design_name' in result:
            backbone = result['design_name']
            if backbone not in backbone_groups:
                backbone_groups[backbone] = []
            backbone_groups[backbone].append(result)
    
    for backbone in sorted(backbone_groups.keys()):
        group = backbone_groups[backbone]
        # Show statistics for each backbone group
        valid_plddt = [r['mean_plddt'] for r in group if r['mean_plddt'] is not None]
        valid_int_pae = [r['interaction_pae'] for r in group if r['interaction_pae'] is not None]
        
        avg_plddt = sum(valid_plddt)/len(valid_plddt) if valid_plddt else 0
        avg_int_pae = sum(valid_int_pae)/len(valid_int_pae) if valid_int_pae else 0
        
        print(f"  {backbone}: {len(group)} designs, Avg pLDDT {avg_plddt:.1f}, Avg Int_PAE {avg_int_pae:.1f}")

def save_results_csv(results, output_file):
    """Save results to CSV file"""
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            csv_cutoffs = sorted(set(DEFAULT_CUTOFFS + [PRIMARY_CUTOFF]))
            side_headers = [f'sc_<={c}A' for c in csv_cutoffs]
            backbone_headers = [f'bb_<={c}A' for c in csv_cutoffs]
            
            # CDR-specific RMSD headers
            cdr_rmsd_headers = [
                'H1_RMSD', 'H2_RMSD', 'H3_RMSD',
                'L1_RMSD', 'L2_RMSD', 'L3_RMSD'
            ]
            
            # CDR length headers
            cdr_length_headers = [
                'H1_Length', 'H2_Length', 'H3_Length',
                'L1_Length', 'L2_Length', 'L3_Length'
            ]

            header_parts = [
                "Filename",
                "Design",
                "Sequence",
                "pLDDT",
                "Interaction_PAE",
                "Mean_PAE",
                "Target_Ab_RMSD",
                "Target_CDR_RMSD",
                "Framework_Ab_RMSD",
                "Framework_CDR_RMSD",
            ] + cdr_rmsd_headers + cdr_length_headers + side_headers + backbone_headers
            f.write(','.join(header_parts) + "\n")

            # Write data rows (without rank number)
            for result in results:
                filename = result['file'].replace('_best.pdb', '')
                design = result.get('design_name', 'N/A')
                sequence = result.get('sequence_name', 'N/A')
                plddt = result['mean_plddt'] if result['mean_plddt'] is not None else ''
                int_pae = result['interaction_pae'] if result['interaction_pae'] is not None else ''
                mean_pae = result['mean_pae'] if result['mean_pae'] is not None else ''
                target_ab_rmsd = result['target_aligned_antibody_rmsd'] if result['target_aligned_antibody_rmsd'] is not None else ''
                target_cdr_rmsd = result['target_aligned_cdr_rmsd'] if result['target_aligned_cdr_rmsd'] is not None else ''
                framework_ab_rmsd = result['framework_aligned_antibody_rmsd'] if result['framework_aligned_antibody_rmsd'] is not None else ''
                framework_cdr_rmsd = result['framework_aligned_cdr_rmsd'] if result['framework_aligned_cdr_rmsd'] is not None else ''
                
                # CDR-specific RMSDs
                cdr_rmsd_values = [
                    result.get('framework_aligned_H1_rmsd', ''),
                    result.get('framework_aligned_H2_rmsd', ''),
                    result.get('framework_aligned_H3_rmsd', ''),
                    result.get('framework_aligned_L1_rmsd', ''),
                    result.get('framework_aligned_L2_rmsd', ''),
                    result.get('framework_aligned_L3_rmsd', '')
                ]
                
                # CDR lengths
                cdr_length_values = [
                    result.get('H1_length', ''),
                    result.get('H2_length', ''),
                    result.get('H3_length', ''),
                    result.get('L1_length', ''),
                    result.get('L2_length', ''),
                    result.get('L3_length', '')
                ]
                
                side_counts = [result.get(key, '') for key in side_headers]
                bb_counts = [result.get(key, '') for key in backbone_headers]

                row_values = [
                    filename,
                    design,
                    sequence,
                    plddt,
                    int_pae,
                    mean_pae,
                    target_ab_rmsd,
                    target_cdr_rmsd,
                    framework_ab_rmsd,
                    framework_cdr_rmsd,
                ] + cdr_rmsd_values + cdr_length_values + side_counts + bb_counts

                f.write(','.join(str(v) for v in row_values) + "\n")
        
        print(f"\nResults saved to: {output_file}")
        return True
    except Exception as e:
        print(f"Failed to save CSV file: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Analyze RF2 original scores from PDB file SCORE lines')
    parser.add_argument('--rf2-dir', default='/home/outputs/rf2', 
                       help='RF2 output directory path')
    parser.add_argument('--output', default='rf2_pdb_analysis.csv',
                       help='Output CSV filename')
    parser.add_argument('--top-n', type=int, default=10,
                       help='Show top N results')
    
    args = parser.parse_args()
    
    # Analyze results
    results = analyze_rf2_pdb_results(args.rf2_dir)
    
    if results:
        # Print analysis results
        print_analysis_results(results)
        
        # Save to CSV
        save_results_csv(results, args.output)
        
        print(f"\nAnalysis complete! Found {len(results)} designs.")
        print(f"Results saved to: {args.output}")
    else:
        print("Analysis failed. Please check RF2 output files")

if __name__ == "__main__":
    main()
