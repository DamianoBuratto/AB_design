#!/usr/bin/env python3
"""Build AlphaFold3 input JSON from FASTA files."""

import json
import sys
import os
import random

def read_fasta(fasta_file):
    """Read FASTA file and return sequence string."""
    if not os.path.isfile(fasta_file):
        print(f"Error: FASTA file not found: {fasta_file}", file=sys.stderr)
        sys.exit(1)
    
    try:
        with open(fasta_file) as f:
            lines = f.readlines()
        seq = ''.join([l.strip() for l in lines if not l.startswith('>')])
        
        if not seq:
            print(f"Error: No sequence found in {fasta_file}", file=sys.stderr)
            sys.exit(1)
        
        return seq
    except Exception as e:
        print(f"Error reading {fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)

chains = [
    {"id": "H", "file": "antibody_H.fasta"},
    {"id": "L", "file": "antibody_L.fasta"},
    {"id": "A", "file": "pHLA_A.fasta"},
    {"id": "C", "file": "pHLA_C.fasta"},
]

sequences = []
for c in chains:
    seq = read_fasta(c["file"])
    sequences.append({"protein": {"id": c["id"], "sequence": seq}})
    print(f"Loaded chain {c['id']}: {len(seq)} residues")

# Generate 10 random seeds for model diversity
random.seed()  # Use system time for randomness
model_seeds = [random.randint(1000, 9999) for _ in range(10)]

# Log seeds for reproducibility
try:
    with open("af3_seeds.log", "w") as log_f:
        log_f.write(f"Generated seeds: {model_seeds}\n")
    print(f"Generated 10 random seeds: {model_seeds}")
except Exception as e:
    print(f"Warning: Could not write seed log: {e}", file=sys.stderr)

# Build AF3 input JSON
af3_input = {
    "name": "Antibody-pHLA Complex",
    "modelSeeds": model_seeds,
    "sequences": sequences,
    "dialect": "alphafold3",
    "version": 1
}

# Write JSON file
try:
    with open("example.json", "w") as f:
        json.dump(af3_input, f, indent=2)
    print("Successfully created example.json")
except Exception as e:
    print(f"Error writing JSON: {e}", file=sys.stderr)
    sys.exit(1)
