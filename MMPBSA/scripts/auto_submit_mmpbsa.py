#!/usr/bin/env python3
"""
Batch MMPBSA Auto-Submit Manager for pHLA MD results (CHARMM36m)

Workflow:
  1. Scan MD results directory for completed production replicas
     (requires DONE sentinel in both setup/ and replica_K/ sub-dirs).
  2. For each completed MD replica that has no MMPBSA result yet,
     generate and submit a SLURM job (up to MAX_JOBS total in queue).
  3. Each SLURM job:
       a. Copy topology files
       b. Pure Python: assign chain IDs via topol.top + .itp → {design}_chAB.pdb
       c. gmx make_ndx: MMPBSA index (prot = HLA chain A, pep = peptide chain B)
       d. gmx convert-tpr: protein-only TPR (default "Protein" group)
       e. Python heredoc: strip solvent from topol.top → topol_prot.top
       f. gmx trjconv: protein-only PBC-corrected trajectory
       g. gmx_MMPBSA -nogui (Linear PB method)
  4. Write {base_dir}/mmpbsa_status.md after every cycle.

Usage:
  python3 auto_submit_mmpbsa.py \\
      --md-results  /public/home/xuziyi/MD/results \\
      --output      /public/home/xuziyi/MMPBSA/results \\
      --mmpbsa-in   /public/home/xuziyi/MMPBSA/inputs/mmpbsa_LinearPB.in \\
      --script-dir  /public/home/xuziyi/MMPBSA/scripts \\
    [--max-jobs 16] [--interval 300] [--partition multi] [--dry-run]
"""

import argparse
import atexit
import fcntl
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# ── Tunables ─────────────────────────────────────────────────────────────────
DEFAULT_MAX_JOBS  = 16         # max SLURM jobs in 'multi' partition
DEFAULT_INTERVAL  = 300        # seconds between queue polls (5 min)
DEFAULT_PARTITION = "multi"
DEFAULT_REPLICAS  = 3          # number of MD replicas per design
MAX_RETRIES       = 2          # max resubmit attempts for failed jobs
# ─────────────────────────────────────────────────────────────────────────────

# HPC environment settings (adjust if your cluster differs)
HPC_CONDA_BASE   = "/public/home/xuziyi/miniconda"
HPC_CONDA_ENV    = "gmxMMPBSA"
HPC_GMX_MODULE   = "gromacs/2024.2"
HPC_GMXLIB       = "/public/home/xuziyi/FEP/force_fields/mutff"
# SLURM resources for MMPBSA (CPU-only; do not request GPU)
HPC_SLURM_CPUS   = 16
HPC_SLURM_MEM    = "16G"
HPC_SLURM_TIME   = "1-00:00:00"
HPC_NODE_EXCLUDE = (
    "node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,"
    "node11,node12,node13,node14,node15,node20,node21,node24,node26"
)


# ── Logging ───────────────────────────────────────────────────────────────────
def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


# ── SLURM helpers ─────────────────────────────────────────────────────────────
def count_user_jobs(user: str, partition: str = DEFAULT_PARTITION) -> int:
    """Count active SLURM jobs for *user* in *partition*."""
    try:
        result = subprocess.run(
            ["squeue", "-u", user, "-p", partition, "-h"],
            capture_output=True, text=True, timeout=15
        )
        lines = [l for l in result.stdout.splitlines() if l.strip()]
        return len(lines)
    except Exception as e:
        log(f"WARNING: squeue failed ({e}); assuming 0 jobs")
        return 0


def submit_slurm(script_path: Path, dry_run: bool = False) -> str:
    """Submit a SLURM script; return job-id string or '' on dry-run/error."""
    if dry_run:
        log(f"  [dry-run] would submit: {script_path}")
        return "DRY_RUN"
    try:
        result = subprocess.run(
            ["sbatch", str(script_path)],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode == 0:
            m = re.search(r"Submitted batch job (\d+)", result.stdout)
            return m.group(1) if m else "UNKNOWN"
        else:
            log(f"  ERROR sbatch: {result.stderr.strip()}")
            return ""
    except Exception as e:
        log(f"  ERROR submitting {script_path}: {e}")
        return ""


def get_current_user() -> str:
    return os.environ.get("USER", subprocess.run(
        ["whoami"], capture_output=True, text=True).stdout.strip())


# ── MD results scanning ───────────────────────────────────────────────────────
def find_completed_replicas(md_results: Path, num_replicas: int = DEFAULT_REPLICAS):
    """
    Yield (source_group, design_name, replica_num, paths_dict) for every
    MD replica that has BOTH setup/DONE and replica_K/DONE.

    paths_dict keys:
        setup_dir      : Path to setup/ directory
        replica_dir    : Path to replica_K/ directory
        tpr            : Path to md_rK.tpr
        xtc            : Path to md_rK.xtc
        gro            : Path to md_rK.gro
        topol          : Path to setup/topol.top
        npt_gro        : Path to setup/npt.gro
    """
    for source_group in sorted(md_results.iterdir()):
        if not source_group.is_dir():
            continue
        for design_dir in sorted(source_group.iterdir()):
            if not design_dir.is_dir():
                continue

            setup_dir = design_dir / "setup"
            if not (setup_dir / "DONE").exists():
                continue
            if not (setup_dir / "topol.top").exists():
                continue
            if not (setup_dir / "npt.gro").exists():
                continue

            for rep in range(1, num_replicas + 1):
                replica_dir = design_dir / f"replica_{rep}"
                if not (replica_dir / "DONE").exists():
                    continue

                tpr = replica_dir / f"md_r{rep}.tpr"
                xtc = replica_dir / f"md_r{rep}.xtc"
                gro = replica_dir / f"md_r{rep}.gro"

                if not tpr.exists() or not xtc.exists():
                    continue

                yield (
                    source_group.name,
                    design_dir.name,
                    rep,
                    {
                        "setup_dir":   setup_dir,
                        "replica_dir": replica_dir,
                        "tpr":         tpr,
                        "xtc":         xtc,
                        "gro":         gro,
                        "topol":       setup_dir / "topol.top",
                        "npt_gro":     setup_dir / "npt.gro",
                    }
                )


def mmpbsa_output_dir(output_root: Path, source_group: str,
                      design_name: str, rep: int) -> Path:
    return output_root / source_group / design_name / f"mmpbsa_r{rep}"


def is_mmpbsa_done(out_dir: Path) -> bool:
    return (out_dir / "DONE").exists()


def is_mmpbsa_submitted(out_dir: Path) -> bool:
    """Check if a SLURM script for this MMPBSA job has been generated."""
    return out_dir.exists() and any(out_dir.glob("submit_mmpbsa_*.slurm"))


def is_mmpbsa_failed(out_dir: Path) -> bool:
    """
    Detect failure in a previously submitted MMPBSA job.
    Checks both .err and .out SLURM logs for common error indicators.
    Returns True if we are confident the job has finished with an error.
    """
    if not out_dir.exists():
        return False
    # Keywords that indicate a definite failure
    fail_patterns = ["ERROR", "error:", "Traceback", "CANCELLED", "FAILED",
                     "Segmentation fault", "Bus error", "DUE TO TIME LIMIT",
                     "slurmstepd: error", "exit code"]
    for pattern in ["slurm_mmpbsa_*.err", "slurm_mmpbsa_*.out"]:
        for f in out_dir.glob(pattern):
            if f.stat().st_size == 0:
                continue
            try:
                text = f.read_text(errors="replace")
                if any(kw in text for kw in fail_patterns):
                    return True
            except OSError:
                pass
    return False


# ── Topology helper ───────────────────────────────────────────────────────────
SOLVENT_NAMES = frozenset([
    "SOL", "TIP3", "TIP3P", "WAT", "HOH",   # water
    "NA", "NA+", "SOD",                       # sodium
    "CL", "CL-", "CLA",                       # chloride
    "K", "K+", "MG", "CA", "ZN",             # other ions
])

TOPOLOGY_STRIP_SCRIPT = r"""
import sys, re

with open('topol.top', 'r') as f:
    lines = f.readlines()

in_molecules = False
out = []
for line in lines:
    stripped = line.strip()
    if re.match(r'\[\s*molecules\s*\]', stripped):
        in_molecules = True
        out.append(line)
        continue
    if in_molecules:
        if stripped == '' or stripped.startswith(';') or stripped.startswith('['):
            out.append(line)
        else:
            mol_name = stripped.split()[0] if stripped.split() else ''
            if mol_name.upper() in {
                'SOL','TIP3','TIP3P','WAT','HOH',
                'NA','NA+','SOD','CL','CL-','CLA',
                'K','K+','MG','CA','ZN'
            }:
                sys.stderr.write(f"  Removed solvent entry: {stripped}\n")
            else:
                out.append(line)
    else:
        out.append(line)

with open('topol_prot.top', 'w') as f:
    f.writelines(out)
print("topol_prot.top written")
"""

GRO_TO_PDB_CHAIN_SCRIPT = r"""
import re
from pathlib import Path


def parse_itp_atom_counts(itp_path: Path):
    counts = {}
    current_mol = None
    in_atoms = False
    atom_count = 0
    for raw in itp_path.read_text(errors='replace').splitlines():
        line = raw.split(';', 1)[0].strip()
        if not line:
            continue
        if line.startswith('[') and line.endswith(']'):
            section = line[1:-1].strip().lower()
            if section == 'moleculetype':
                if current_mol and atom_count > 0:
                    counts[current_mol] = atom_count
                current_mol = None
                in_atoms = False
                atom_count = 0
            elif section == 'atoms':
                in_atoms = True
            else:
                in_atoms = False
            continue
        if current_mol is None:
            parts = line.split()
            if parts:
                current_mol = parts[0]
            continue
        if in_atoms:
            parts = line.split()
            if parts and parts[0].isdigit():
                atom_count += 1
    if current_mol and atom_count > 0:
        counts[current_mol] = atom_count
    return counts


def parse_topol(topol_path: Path):
    include_files = []
    molecules = []
    in_molecules = False

    for raw in topol_path.read_text(errors='replace').splitlines():
        line = raw.split(';', 1)[0].strip()
        if not line:
            continue
        m = re.match(r'#include\s+"([^"]+)"', line)
        if m:
            include_files.append(Path(m.group(1)).name)
            continue
        if re.match(r'\[\s*molecules\s*\]', line):
            in_molecules = True
            continue
        if in_molecules:
            if line.startswith('['):
                break
            parts = line.split()
            if len(parts) >= 2 and parts[1].isdigit():
                molecules.append((parts[0], int(parts[1])))
    return include_files, molecules


def build_mol_atom_counts(topol_path: Path):
    include_files, molecules = parse_topol(topol_path)
    mol_atom = {}
    for fn in include_files:
        fp = topol_path.parent / fn
        if not fp.is_file() or fp.suffix.lower() != '.itp':
            continue
        mol_atom.update(parse_itp_atom_counts(fp))

    expanded = []
    for mol, n in molecules:
        if mol not in mol_atom:
            continue
        for _ in range(n):
            expanded.append((mol, mol_atom[mol]))
    return expanded


def parse_gro_atoms(gro_path: Path):
    lines = gro_path.read_text(errors='replace').splitlines()
    natom = int(lines[1].strip())
    atom_lines = lines[2:2 + natom]
    atoms = []
    for line in atom_lines:
        resnr = int(line[0:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnr = int(line[15:20])
        x = float(line[20:28]) * 10.0
        y = float(line[28:36]) * 10.0
        z = float(line[36:44]) * 10.0
        atoms.append((resnr, resname, atomname, atomnr, x, y, z))
    return atoms


SOLVENT = {
    'SOL', 'TIP3', 'TIP3P', 'WAT', 'HOH',
    'NA', 'NA+', 'SOD', 'CL', 'CL-', 'CLA',
    'K', 'K+', 'MG', 'CA', 'ZN'
}


def write_pdb_chains(gro_path: Path, topol_path: Path, out_pdb: Path):
    atoms = parse_gro_atoms(gro_path)
    mols = build_mol_atom_counts(topol_path)
    if len(mols) < 2:
        raise RuntimeError('Could not identify first two molecules from topol/itp')

    pep_count = mols[0][1]
    hla_count = mols[1][1]
    total = pep_count + hla_count
    if total > len(atoms):
        raise RuntimeError('Protein atom count exceeds atoms in GRO')

    pdb_lines = []
    serial = 1
    last_resid = {'A': None, 'B': None}
    resseq = {'A': 0, 'B': 0}

    for i, (resnr, resname, atomname, _atomnr, x, y, z) in enumerate(atoms[:total], start=1):
        chain = 'B' if i <= pep_count else 'A'
        if resname.upper() in SOLVENT:
            continue
        if last_resid[chain] != resnr:
            resseq[chain] += 1
            last_resid[chain] = resnr
        element = (atomname.strip()[0] if atomname.strip() else 'C').upper()
        pdb_lines.append(
            f"ATOM  {serial:5d} {atomname:<4s} {resname:>3s} {chain}{resseq[chain]:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element:>2s}"
        )
        serial += 1

    pdb_lines.append('TER')
    pdb_lines.append('END')
    out_pdb.write_text('\n'.join(pdb_lines) + '\n')
    print(f'Wrote {out_pdb} with chain B=peptide, chain A=HLA; atoms={serial-1}')


write_pdb_chains(Path('npt.gro'), Path('topol.top'), Path(f"{DESIGN}_chAB.pdb"))
"""


# ── SLURM script generation ───────────────────────────────────────────────────
def generate_mmpbsa_slurm(
    source_group: str,
    design_name: str,
    rep: int,
    paths: dict,
    out_dir: Path,
    mmpbsa_in: Path,
    script_dir: Path,
    gmx: str = "gmx",
    partition: str = DEFAULT_PARTITION,
) -> Path:
    """
    Generate a self-contained SLURM script for one MMPBSA calculation.

    The script:
      1. Loads GROMACS and activates gmxMMPBSA conda env
      2. Copies topology + trajectory files
      3. Pure Python: assign chain IDs (topol.top + .itp) → {design}_chAB.pdb
      4. gmx make_ndx: MMPBSA index (chain A=prot, chain B=pep)
      5. gmx make_ndx + convert-tpr: protein-only TPR
      6. Python: strips solvent from topol.top → topol_prot.top
      7. gmx trjconv: protein-only PBC-corrected trajectory
      8. gmx_MMPBSA (nogui, LinearPB)
      9. Touches DONE on success
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    job_name   = f"mmpbsa_{design_name}_r{rep}"
    script_path = out_dir / f"submit_mmpbsa_{design_name}_r{rep}.slurm"

    setup_dir   = paths["setup_dir"]
    replica_dir = paths["replica_dir"]
    tpr         = paths["tpr"]
    xtc         = paths["xtc"]
    npt_gro     = paths["npt_gro"]

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --exclude={HPC_NODE_EXCLUDE}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={HPC_SLURM_CPUS}
#SBATCH --mem={HPC_SLURM_MEM}
#SBATCH --time={HPC_SLURM_TIME}
#SBATCH --output={out_dir}/slurm_mmpbsa_%j.out
#SBATCH --error={out_dir}/slurm_mmpbsa_%j.err

# ===========================================================================
# gmx_MMPBSA calculation: {design_name} replica {rep}
#   Source group : {source_group}
#   MD setup dir : {setup_dir}
#   MD replica   : {replica_dir}
#   MMPBSA out   : {out_dir}
# ===========================================================================

# ── Load GROMACS ──────────────────────────────────────────────────────────
module load {HPC_GMX_MODULE}
export GMXLIB={HPC_GMXLIB}
export GMX_MAXCONSTRWARN=-1

GMX="{gmx}"
DESIGN="{design_name}"
REP={rep}

# ── Activate gmxMMPBSA conda environment ─────────────────────────────────
source {HPC_CONDA_BASE}/etc/profile.d/conda.sh
conda activate {HPC_CONDA_ENV} || {{ echo "ERROR: conda activate {HPC_CONDA_ENV} failed"; exit 1; }}
set -euo pipefail
export OMP_NUM_THREADS={HPC_SLURM_CPUS}
echo "Python: $(python3 --version)"
gmx_MMPBSA --version 2>&1 | head -1 || true
echo ""

# ── Enter output directory ────────────────────────────────────────────────
mkdir -p {out_dir}
cd {out_dir} || {{ echo "ERROR: cannot cd to {out_dir}"; exit 1; }}

echo "======================================================"
echo " MMPBSA: ${{DESIGN}} r${{REP}}  on $(hostname)  job $SLURM_JOB_ID"
echo " $(date)"
echo "======================================================"

# ── Step 1: Copy topology files ───────────────────────────────────────────
echo ""
echo "[1/7] Copying topology files from setup/..."
cp {setup_dir}/topol.top ./                        || {{ echo "ERROR: missing topol.top"; exit 1; }}
cp {setup_dir}/*.itp ./          2>/dev/null       || true
cp {setup_dir}/posre*.itp ./     2>/dev/null       || true
cp {npt_gro} ./npt.gro                             || {{ echo "ERROR: missing npt.gro";  exit 1; }}
cp {tpr} ./md_r${{REP}}.tpr                        || {{ echo "ERROR: missing tpr";      exit 1; }}
cp {xtc} ./md_r${{REP}}.xtc                        || {{ echo "ERROR: missing xtc";      exit 1; }}
echo "    Done."

# ── Step 2: Build chain-assigned PDB via Python (no VMD dependency) ──────
# Uses topol.top + .itp molecule atom counts to split first two molecules:
# molecule #1 -> chain B (peptide), molecule #2 -> chain A (HLA)
# Output: ${{DESIGN}}_chAB.pdb
echo ""
echo "[2/7] Building chain-assigned PDB without VMD..."
python3 - <<PYEOF || {{ echo "ERROR: Python chain assignment failed"; exit 1; }}
DESIGN = "${{DESIGN}}"
{GRO_TO_PDB_CHAIN_SCRIPT}
PYEOF
echo "    Written: ${{DESIGN}}_chAB.pdb  (chain A=HLA, chain B=peptide)"

# ── Step 3: Create MMPBSA index file ─────────────────────────────────────
# Uses chAB.pdb (with proper chains) so group numbers are deterministic:
#   Group 0 = all atoms  (from 'keep 0')
#   Group 1 = chain A = HLA   → named 'prot'
#   Group 2 = chain B = peptide → named 'pep'
echo ""
echo "[3/7] Creating MMPBSA index (prot=HLA, pep=peptide)..."
printf "keep 0\\nsplitch 0\\nname 1 prot\\nname 2 pep\\n\\nq\\n" | \\
    $GMX make_ndx -f ${{DESIGN}}_chAB.pdb -o index_mmpbsa.ndx > index_mmpbsa.log 2>&1 \\
    || {{ echo "ERROR: make_ndx for MMPBSA failed."; cat index_mmpbsa.log; exit 1; }}
echo "    Done: index_mmpbsa.ndx"

# ── Step 4: Create protein-only TPR ──────────────────────────────────────
echo ""
echo "[4/7] Converting TPR to protein-only..."
# Create index selecting only protein (no water, no ions)
printf '"Protein"\\nq\\n' | \\
    $GMX make_ndx -f md_r${{REP}}.tpr -o index_protein.ndx > index_protein.log 2>&1 \\
    || {{ echo "ERROR: make_ndx for protein failed."; cat index_protein.log; exit 1; }}

echo "Protein" | \\
    $GMX convert-tpr -s md_r${{REP}}.tpr -n index_protein.ndx \\
    -o ${{DESIGN}}_prot.tpr > convert_tpr.log 2>&1 \\
    || {{ echo "ERROR: convert-tpr failed."; cat convert_tpr.log; exit 1; }}
echo "    Done: ${{DESIGN}}_prot.tpr"

# ── Step 5: Strip solvent from topology ──────────────────────────────────
echo ""
echo "[5/7] Creating protein-only topology (topol_prot.top)..."
python3 - <<'PYEOF'
{TOPOLOGY_STRIP_SCRIPT}
PYEOF
echo "    Done: topol_prot.top"

# ── Step 6: Extract protein-only PBC-corrected trajectory ────────────────
echo ""
echo "[6/7] Extracting protein trajectory (protein-only, PBC corrected)..."
# Center on Protein, output protein atoms only
printf "Protein\\nProtein\\n" | \\
    $GMX trjconv \\
    -f md_r${{REP}}.xtc \\
    -s md_r${{REP}}.tpr \\
    -n index_protein.ndx \\
    -o ${{DESIGN}}_protein_r${{REP}}.xtc \\
    -pbc mol -center > trjconv.log 2>&1 \\
    || {{ echo "ERROR: trjconv failed."; cat trjconv.log; exit 1; }}
NFRAMES=$(grep "^fr" trjconv.log 2>/dev/null | wc -l || echo "?")
echo "    Done: ${{DESIGN}}_protein_r${{REP}}.xtc  (frames: ~$NFRAMES)"

# ── Step 7: Run gmx_MMPBSA ────────────────────────────────────────────────
echo ""
echo "[7/7] Running gmx_MMPBSA (Linear PB)..."
echo "    Input     : {mmpbsa_in}"
echo "    Topology  : topol_prot.top"
echo "    TPR       : ${{DESIGN}}_prot.tpr"
echo "    Index     : index_mmpbsa.ndx  (-cg 1 2)"
echo "    Trajectory: ${{DESIGN}}_protein_r${{REP}}.xtc"
echo "    Reference : ${{DESIGN}}_chAB.pdb"
echo ""

gmx_MMPBSA -O -nogui \\
    -i {mmpbsa_in} \\
    -cp topol_prot.top \\
    -cs ${{DESIGN}}_prot.tpr \\
    -ci index_mmpbsa.ndx \\
    -cg 1 2 \\
    -ct ${{DESIGN}}_protein_r${{REP}}.xtc \\
    -cr ${{DESIGN}}_chAB.pdb \\
    -eo gmx_MMPBSA_plot.csv \\
    > mmpbsa.log 2>&1 \\
    || {{ echo "ERROR: gmx_MMPBSA failed. Check mmpbsa.log"; tail -30 mmpbsa.log; exit 1; }}

# ── Verify output ─────────────────────────────────────────────────────────
if [[ ! -f "FINAL_RESULTS_MMPBSA.dat" ]]; then
    echo "ERROR: FINAL_RESULTS_MMPBSA.dat not found — gmx_MMPBSA may have failed."
    tail -30 mmpbsa.log
    exit 1
fi

echo ""
echo "======================================================"
echo " MMPBSA COMPLETED: ${{DESIGN}} r${{REP}}"
echo " $(date)"
echo "======================================================"
echo ""
grep -A5 "DELTA TOTAL" FINAL_RESULTS_MMPBSA.dat | head -10 || true

touch DONE
echo "DONE sentinel written."
"""

    script_path.write_text(script)
    script_path.chmod(0o755)
    return script_path


# ── Status report ──────────────────────────────────────────────────────────────
def write_status_report(
    md_results: Path,
    output_root: Path,
    base_dir: Path,
    num_replicas: int,
) -> None:
    """Write a Markdown status table to {base_dir}/mmpbsa_status.md."""
    report_path = base_dir / "mmpbsa_status.md"
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        f"# MMPBSA Status Report",
        f"",
        f"Generated: {now}",
        f"MD results : `{md_results}`",
        f"MMPBSA out : `{output_root}`",
        f"",
    ]

    cols = ["Design", "Group"] + [f"r{r}" for r in range(1, num_replicas + 1)]
    header = "| " + " | ".join(cols) + " |"
    sep    = "| " + " | ".join(["---"] * len(cols)) + " |"
    lines += [header, sep]

    for grp in sorted(md_results.iterdir()):
        if not grp.is_dir():
            continue
        for design_dir in sorted(grp.iterdir()):
            if not design_dir.is_dir():
                continue
            row = [design_dir.name, grp.name]
            for rep in range(1, num_replicas + 1):
                out_dir = mmpbsa_output_dir(output_root, grp.name, design_dir.name, rep)
                if is_mmpbsa_done(out_dir):
                    row.append("✓ DONE")
                elif out_dir.exists():
                    if is_mmpbsa_failed(out_dir):
                        row.append("✗ failed")
                    else:
                        row.append("⚙ running")
                else:
                    md_rep_dir = design_dir / f"replica_{rep}"
                    if not (md_rep_dir / "DONE").exists():
                        row.append("— MD pend")
                    else:
                        row.append("— queued")
            lines.append("| " + " | ".join(row) + " |")

    lines.append("")
    with open(report_path, "w") as f:
        f.write("\n".join(lines))
    log(f"Status report written: {report_path}")


# ── State persistence ─────────────────────────────────────────────────────────
STATE_FILE = ".mmpbsa_job_state.json"
LOCK_FILE = ".mmpbsa_manager.lock"

def load_state(base_dir: Path) -> dict:
    p = base_dir / STATE_FILE
    if p.exists():
        return json.loads(p.read_text())
    return {}


def save_state(base_dir: Path, state: dict) -> None:
    (base_dir / STATE_FILE).write_text(json.dumps(state, indent=2))


def acquire_single_instance_lock(base_dir: Path):
    """Ensure only one auto-submit manager instance runs per base_dir."""
    lock_path = base_dir / LOCK_FILE
    lock_handle = open(lock_path, "a+")
    try:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock_handle.seek(0)
        owner = lock_handle.read().strip() or "UNKNOWN"
        log(f"ERROR: another manager instance is already running (owner={owner}).")
        log(f"       lock file: {lock_path}")
        log("       stop the old process or remove stale lock after confirming it is dead.")
        sys.exit(2)

    lock_handle.seek(0)
    lock_handle.truncate()
    lock_handle.write(f"pid={os.getpid()} host={os.uname().nodename} started={datetime.now().isoformat()}\n")
    lock_handle.flush()

    def _release_lock():
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
            lock_handle.close()
        except OSError:
            pass

    atexit.register(_release_lock)
    return lock_handle


# ── Main loop ─────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Batch MMPBSA auto-submit manager for pHLA MD results"
    )
    parser.add_argument("--md-results", required=True,
                        help="MD results root directory (contains replicaN/ subdirs)")
    parser.add_argument("--output", required=True,
                        help="MMPBSA outputs root directory")
    parser.add_argument("--mmpbsa-in", required=True,
                        help="Path to mmpbsa .in input file")
    parser.add_argument("--script-dir", required=True,
                        help="Path to MMPBSA scripts directory")
    parser.add_argument("--base-dir", default=None,
                        help="Base dir for status report and state file (default: --output)")
    parser.add_argument("--max-jobs", type=int, default=DEFAULT_MAX_JOBS,
                        help=f"Max SLURM jobs in partition (default: {DEFAULT_MAX_JOBS})")
    parser.add_argument("--interval", type=int, default=DEFAULT_INTERVAL,
                        help=f"Polling interval in seconds (default: {DEFAULT_INTERVAL})")
    parser.add_argument("--replicas", type=int, default=DEFAULT_REPLICAS,
                        help=f"Number of MD replicas to look for (default: {DEFAULT_REPLICAS})")
    parser.add_argument("--gmx", default="gmx",
                        help="GROMACS executable (default: gmx)")
    parser.add_argument("--partition", default=DEFAULT_PARTITION,
                        help=f"SLURM partition (default: {DEFAULT_PARTITION}, recommended: multi)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Generate scripts but do not submit SLURM jobs")
    parser.add_argument("--once", action="store_true",
                        help="Run one submission cycle then exit (no polling)")
    parser.add_argument("--report", action="store_true",
                        help="Print status report and exit")
    args = parser.parse_args()

    md_results  = Path(args.md_results).resolve()
    output_root = Path(args.output).resolve()
    mmpbsa_in   = Path(args.mmpbsa_in).resolve()
    script_dir  = Path(args.script_dir).resolve()
    base_dir    = Path(args.base_dir).resolve() if args.base_dir else output_root

    if not md_results.is_dir():
        log(f"ERROR: --md-results not found: {md_results}")
        sys.exit(1)
    if not mmpbsa_in.is_file():
        log(f"ERROR: --mmpbsa-in not found: {mmpbsa_in}")
        sys.exit(1)

    output_root.mkdir(parents=True, exist_ok=True)
    base_dir.mkdir(parents=True, exist_ok=True)
    _lock_handle = acquire_single_instance_lock(base_dir)

    user  = get_current_user()
    state = load_state(base_dir)

    log("=" * 60)
    log(f"  MMPBSA Auto-Submit Manager")
    log(f"  MD results : {md_results}")
    log(f"  Output     : {output_root}")
    log(f"  MMPBSA in  : {mmpbsa_in}")
    log(f"  Script dir : {script_dir}")
    log(f"  Partition  : {args.partition}  (mode: SLURM)")
    log(f"  Max jobs   : {args.max_jobs}")
    log(f"  User       : {user}")
    log(f"  Dry-run    : {args.dry_run}")
    log("=" * 60)

    if args.report:
        write_status_report(md_results, output_root, base_dir, args.replicas)
        return

    cycle = 0
    while True:
        cycle += 1
        log(f"\n── Cycle {cycle}  ({datetime.now().strftime('%H:%M:%S')}) ──")

        # ── Count current jobs ─────────────────────────────────────────────
        current_jobs = count_user_jobs(user, partition=args.partition)
        free_slots   = max(0, args.max_jobs - current_jobs)
        log(f"  Queue [{args.partition}]: {current_jobs}/{args.max_jobs} jobs  ({free_slots} free slots)")

        submitted_this_cycle = 0

        # ── Scan for pending MMPBSA work ───────────────────────────────────
        for src_grp, design, rep, paths in find_completed_replicas(md_results, args.replicas):
            out_dir = mmpbsa_output_dir(output_root, src_grp, design, rep)

            # Already done?
            if is_mmpbsa_done(out_dir):
                continue

            # Already submitted (script exists) and not failed?
            key = f"{src_grp}/{design}/r{rep}"
            if is_mmpbsa_submitted(out_dir):
                if not is_mmpbsa_failed(out_dir):
                    continue   # job running or queued (or completed without DONE — rare)
                else:
                    log(f"  !! Failed (resubmitting): {key}")
                    retries = state.get(key, {}).get("retries", 0)
                    if retries >= MAX_RETRIES:
                        log(f"  !! Gave up after {retries} retries: {key}")
                        continue
                    state.setdefault(key, {})["retries"] = retries + 1

            # Submit if slots available
            if free_slots <= 0:
                log(f"  Queue full, deferring: {key}")
                continue

            script = generate_mmpbsa_slurm(
                src_grp, design, rep, paths, out_dir,
                mmpbsa_in, script_dir, gmx=args.gmx,
                partition=args.partition,
            )

            log(f"  → Submitting [{args.partition}]: {key}")
            job_id = submit_slurm(script, dry_run=args.dry_run)
            if job_id:
                state.setdefault(key, {})["job_id"] = job_id
                state.setdefault(key, {})["submitted_at"] = datetime.now().isoformat()
                free_slots -= 1
                submitted_this_cycle += 1
                log(f"    Submitted job {job_id}: {key}")

        save_state(base_dir, state)
        write_status_report(md_results, output_root, base_dir, args.replicas)

        log(f"  Submitted {submitted_this_cycle} job(s) this cycle.")

        if args.once or args.dry_run:
            log("Exiting after one cycle (--once / --dry-run).")
            break

        log(f"  Sleeping {args.interval}s (Ctrl+C to stop)...")
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
