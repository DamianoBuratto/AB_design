#!/usr/bin/env python3
"""
Batch MMPBSA Auto-Submit Manager for Ab-pHLA binding (CHARMM36m)

Computes the binding free energy of ANTIBODY (H+L chains) to the
pHLA complex (peptide P + HLA A), NOT the peptide-HLA binding.

Workflow:
  1. Scan MD results for completed production replicas
     (requires DONE sentinel in both setup/ and replica_K/).
  2. For each completed MD replica without an MMPBSA result,
     generate and submit a SLURM job (up to MAX_JOBS total in queue).
  3. Each SLURM job executes a 7-step pipeline:
         1/7  Copy topology files from MD setup/
         2/7  Pure Python chain-ID assignment (topol.top + .itp) ->
                {design}_chPAHL.pdb  (P=peptide, A=HLA, H=Ab-heavy, L=Ab-light)
                index_protein.ndx     (group Protein, all 4 chain atoms)
                index_mmpbsa.ndx      (group receptor=pHLA P+A, group ligand=Ab H+L)
         3/7  gmx trjconv: extract full protein (4-chain) reference structure from npt.gro
         4/7  Python: trim topol.top to first 4 molecules, remove posre includes
                -> topol_prot.top (4 protein chains, no position restraints)
         5/7  gmx grompp: build protein-only TPR from stripped topology + Protein.gro
                (define= empty, no position restraints, -maxwarn 2)
         6/7  gmx trjconv: protein-only PBC-corrected trajectory
                (center=Protein, output=Protein, TPR=full-system md_rN.tpr)
         7/7  gmx_MMPBSA -O -nogui  (Linear PB, -cg 0 1)
             receptor = pHLA (chains P+A), ligand = Ab (chains H+L)
  4. Write {base_dir}/mmpbsa_status.md after every cycle.

Usage:
  python3 auto_submit_mmpbsa.py \\
      --md-results  /public/home/xuziyi/MD/results \\
      --output      /public/home/xuziyi/MMPBSA/results \\
      --mmpbsa-in   /public/home/xuziyi/MMPBSA/inputs/mmpbsa_LinearPB.in \\
      [--targets-file /path/to/targets.tsv] \\
    [--max-jobs 16] [--interval 300] [--partition multi] [--dry-run] [--once] [--report]
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

# ── Tunables ──────────────────────────────────────────────────────────────────
DEFAULT_MAX_JOBS  = 16
DEFAULT_INTERVAL  = 300        # seconds between queue polls
DEFAULT_PARTITION = "multi"
DEFAULT_REPLICAS  = 3
MAX_RETRIES       = 2
# ─────────────────────────────────────────────────────────────────────────────

# HPC environment -- adjust if cluster settings change
HPC_CONDA_BASE   = "/public/home/xuziyi/miniconda"
HPC_CONDA_ENV    = "gmxMMPBSA"
HPC_GMX_MODULE   = "gromacs/2024.2"
HPC_GMXLIB       = "/public/home/xuziyi/FEP/force_fields/mutff"
HPC_SLURM_CPUS   = 16
HPC_SLURM_MEM    = "32G"         # ~6400 atoms (2x larger than pHLA-only)
HPC_SLURM_TIME   = "2-00:00:00" # more time for larger system
HPC_NODE_EXCLUDE = (
    "node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,"
    "node11,node12,node13,node14,node15,node20,node21,node24,node26"
)

STATE_FILE = ".mmpbsa_job_state.json"
LOCK_FILE  = ".mmpbsa_manager.lock"


# ── Logging ───────────────────────────────────────────────────────────────────
def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


# ── SLURM helpers ─────────────────────────────────────────────────────────────
def get_current_user() -> str:
    return os.environ.get("USER", subprocess.run(
        ["whoami"], capture_output=True, text=True).stdout.strip())


def count_user_jobs(user: str, partition: str = DEFAULT_PARTITION) -> int:
    """Count active/pending SLURM jobs for *user* in *partition*."""
    try:
        result = subprocess.run(
            ["squeue", "-u", user, "-p", partition, "-h"],
            capture_output=True, text=True, timeout=15,
        )
        return len([l for l in result.stdout.splitlines() if l.strip()])
    except Exception as e:
        log(f"WARNING: squeue failed ({e}); assuming 0 jobs")
        return 0


def submit_slurm(script_path: Path, dry_run: bool = False) -> str:
    """Submit a SLURM script; return job-id string or '' on error."""
    if dry_run:
        log(f"  [dry-run] would submit: {script_path}")
        return "DRY_RUN"
    try:
        result = subprocess.run(
            ["sbatch", str(script_path)],
            capture_output=True, text=True, timeout=30,
        )
        if result.returncode == 0:
            m = re.search(r"Submitted batch job (\d+)", result.stdout)
            return m.group(1) if m else "UNKNOWN"
        log(f"  ERROR sbatch: {result.stderr.strip()}")
        return ""
    except Exception as e:
        log(f"  ERROR submitting {script_path}: {e}")
        return ""


def get_active_jobs(user: str, partition: str = DEFAULT_PARTITION) -> dict[str, str]:
    """Return {job_name: job_id} for active/pending jobs for *user* in *partition*."""
    try:
        result = subprocess.run(
            ["squeue", "-u", user, "-p", partition, "-h", "-o", "%i|%j"],
            capture_output=True, text=True, timeout=15,
        )
        jobs = {}
        for line in result.stdout.splitlines():
            line = line.strip()
            if not line or "|" not in line:
                continue
            job_id, job_name = line.split("|", 1)
            jobs[job_name.strip()] = job_id.strip()
        return jobs
    except Exception as e:
        log(f"WARNING: squeue(job names) failed ({e}); assuming no active MMPBSA jobs")
        return {}


def slurm_job_name(source_group: str, design_name: str, rep: int) -> str:
    return f"mmpbsa_{source_group}_{design_name}_r{rep}"


# ── Single-instance lock ──────────────────────────────────────────────────────
def acquire_single_instance_lock(base_dir: Path):
    lock_path   = base_dir / LOCK_FILE
    lock_handle = open(lock_path, "a+")
    try:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock_handle.seek(0)
        owner = lock_handle.read().strip() or "UNKNOWN"
        log(f"ERROR: another manager instance is running (owner={owner}).")
        log(f"       lock: {lock_path}")
        log("       Stop the old process, or remove stale lock if process is dead.")
        sys.exit(2)

    lock_handle.seek(0)
    lock_handle.truncate()
    lock_handle.write(
        f"pid={os.getpid()} host={os.uname().nodename} "
        f"started={datetime.now().isoformat()}\n"
    )
    lock_handle.flush()

    def _release():
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
            lock_handle.close()
        except OSError:
            pass

    atexit.register(_release)
    return lock_handle


# ── MD results scanning ───────────────────────────────────────────────────────
def find_completed_replicas(
    md_results: Path,
    num_replicas: int = DEFAULT_REPLICAS,
    targets: set[tuple[str, str]] | None = None,
):
    """
    Yield (source_group, design_name, replica_num, paths_dict) for every
    MD replica with DONE in both setup/ and replica_K/.

    If *targets* is provided, only yield designs matching
    (source_group, design_name) in the set.
    """
    for source_group in sorted(md_results.iterdir()):
        if not source_group.is_dir():
            continue
        for design_dir in sorted(source_group.iterdir()):
            if not design_dir.is_dir():
                continue
            if targets is not None and (source_group.name, design_dir.name) not in targets:
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
                        "npt_gro":     setup_dir / "npt.gro",
                        "topol":       setup_dir / "topol.top",
                    },
                )


# ── MMPBSA output helpers ─────────────────────────────────────────────────────
def mmpbsa_output_dir(output_root: Path, source_group: str,
                      design_name: str, rep: int) -> Path:
    return output_root / source_group / design_name / f"mmpbsa_r{rep}"


def is_mmpbsa_done(out_dir: Path) -> bool:
    return (out_dir / "DONE").exists() or (out_dir / "FINAL_RESULTS_MMPBSA.dat").exists()


def has_submission_artifacts(out_dir: Path) -> bool:
    """True if output dir has .slurm script or SLURM log files (may be stale)."""
    return out_dir.exists() and (
        any(out_dir.glob("submit_mmpbsa_*.slurm"))
        or any(out_dir.glob("slurm_mmpbsa_*.out"))
        or any(out_dir.glob("slurm_mmpbsa_*.err"))
    )


def is_mmpbsa_failed(out_dir: Path) -> bool:
    """
    Return True only when unambiguous failure keywords appear in SLURM logs.
    Avoids false positives from normal GROMACS informational messages.
    """
    if not out_dir.exists():
        return False
    fail_keywords = [
        "Traceback (most recent call last)",
        "CANCELLED",
        "DUE TO TIME LIMIT",
        "Segmentation fault",
        "Bus error",
        "slurmstepd: error",
        "exit code 1",
        "ERROR: ",      # uppercase: only our explicit echo "ERROR: ..." guards
        "FATAL ERROR",  # GROMACS fatal
    ]
    for pattern in ("slurm_mmpbsa_*.err", "slurm_mmpbsa_*.out"):
        for f in out_dir.glob(pattern):
            if f.stat().st_size == 0:
                continue
            try:
                text = f.read_text(errors="replace")
                if any(kw in text for kw in fail_keywords):
                    return True
            except OSError:
                pass
    return False


# ── State persistence ─────────────────────────────────────────────────────────
def load_state(base_dir: Path) -> dict:
    p = base_dir / STATE_FILE
    if p.exists():
        try:
            return json.loads(p.read_text())
        except Exception:
            pass
    return {}


def save_state(base_dir: Path, state: dict) -> None:
    (base_dir / STATE_FILE).write_text(json.dumps(state, indent=2))


def parse_targets_file(targets_file: Path) -> set[tuple[str, str]]:
    """
    Parse whitespace-separated target lines: <outer_replica> <design_name>.
    Blank lines and lines starting with '#' are ignored.
    """
    targets = set()
    for lineno, raw in enumerate(targets_file.read_text(errors="replace").splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(
                f"Invalid targets entry at {targets_file}:{lineno}: "
                f"expected 2 columns (<outer_replica> <design_name>), got: {raw!r}"
            )
        targets.add((parts[0], parts[1]))
    if not targets:
        raise ValueError(f"No targets found in {targets_file}")
    return targets


# ── Embedded Python scripts for the SLURM heredocs ───────────────────────────

# Inserted into an UNQUOTED heredoc (<<PYEOF) so bash expands ${DESIGN}.
# Reads npt.gro + topol.top + .itp to identify all 4 protein molecules.
# molecule #1 (peptide)  -> chain P (receptor part 1)
# molecule #2 (HLA)     -> chain A (receptor part 2)
# molecule #3 (Ab heavy) -> chain H (ligand part 1)
# molecule #4 (Ab light) -> chain L (ligand part 2)
# Outputs: {design}_chPAHL.pdb, index_protein.ndx, index_mmpbsa.ndx
_CHAIN_ASSIGN_PY = r"""
import re
from pathlib import Path


def parse_itp_atom_counts(itp_path):
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
        if in_atoms and line.split()[0].isdigit():
            atom_count += 1
    if current_mol and atom_count > 0:
        counts[current_mol] = atom_count
    return counts


def parse_topol(topol_path):
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


def build_mol_atoms(topol_path):
    includes, molecules = parse_topol(topol_path)
    mol_atom = {}
    for fn in includes:
        fp = topol_path.parent / fn
        if fp.is_file() and fp.suffix.lower() == '.itp':
            mol_atom.update(parse_itp_atom_counts(fp))
    expanded = []
    for mol, n in molecules:
        if mol not in mol_atom:
            continue
        for _ in range(n):
            expanded.append((mol, mol_atom[mol]))
    return expanded


def parse_gro_atoms(gro_path):
    lines = gro_path.read_text(errors='replace').splitlines()
    natom = int(lines[1].strip())
    atoms = []
    for line in lines[2:2 + natom]:
        resnr    = int(line[0:5])
        resname  = line[5:10].strip()
        atomname = line[10:15].strip()
        x = float(line[20:28]) * 10.0
        y = float(line[28:36]) * 10.0
        z = float(line[36:44]) * 10.0
        atoms.append((resnr, resname, atomname, x, y, z))
    return atoms


design   = "${DESIGN}"
gro_path = Path('npt.gro')
top_path = Path('topol.top')
out_pdb  = Path(f'{design}_chPAHL.pdb')

mols  = build_mol_atoms(top_path)
if len(mols) < 4:
    raise RuntimeError(f'Need >= 4 protein molecules for Ab-pHLA MMPBSA; found {len(mols)} in topol+itp. '
                       f'Check that the MD setup includes antibody chains (H, L).')

pep_count  = mols[0][1]   # molecule 1: peptide
hla_count  = mols[1][1]   # molecule 2: HLA alpha
abh_count  = mols[2][1]   # molecule 3: antibody heavy chain
abl_count  = mols[3][1]   # molecule 4: antibody light chain
total_prot = pep_count + hla_count + abh_count + abl_count

atoms = parse_gro_atoms(gro_path)
if total_prot > len(atoms):
    raise RuntimeError(f'Protein atoms ({total_prot}) exceed GRO atoms ({len(atoms)})')

pdb_lines = []
serial    = 1

mol_ends = [0]
for i, (_, count) in enumerate(mols[:4]):
    mol_ends.append(mol_ends[-1] + count)

for i, (resnr, resname, atomname, x, y, z) in enumerate(atoms[:total_prot], start=1):
    if   i <= mol_ends[1]: chain = 'P'  # peptide
    elif i <= mol_ends[2]: chain = 'A'  # HLA
    elif i <= mol_ends[3]: chain = 'H'  # Ab heavy
    elif i <= mol_ends[4]: chain = 'L'  # Ab light
    else: chain = 'X'
    element = (atomname.lstrip('0123456789')[0] if atomname.strip() else 'C').upper()
    pdb_lines.append(
        f"ATOM  {serial:5d} {atomname:<4s} {resname:>3s} {chain}{resnr:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element:>2s}"
    )
    serial += 1

pdb_lines += ['TER', 'END']
out_pdb.write_text('\n'.join(pdb_lines) + '\n')
print(f'chain_assign: {out_pdb}  chain_P(pep)={pep_count}  chain_A(HLA)={hla_count}  '
      f'chain_H(AbH)={abh_count}  chain_L(AbL)={abl_count}  total={serial-1}')


def _write_ndx_group(fh, name, atoms):
    fh.write(f'[ {name} ]\n')
    for start in range(0, len(atoms), 15):
        fh.write(' '.join(str(a) for a in atoms[start:start + 15]) + '\n')


# index_protein.ndx: all 4 protein chain atoms for grompp and trjconv
with open('index_protein.ndx', 'w') as fh:
    _write_ndx_group(fh, 'Protein', range(1, total_prot + 1))
print(f'index_protein.ndx: {total_prot} atoms (all 4 protein chains)')

# index_mmpbsa.ndx for Ab-pHLA binding:
#   receptor = pHLA (peptide + HLA, chains P+A): atoms 1..pHLA_end
#   ligand   = antibody (H+L, chains):         atoms pHLA_end+1..total_prot
# gmx_MMPBSA -cg uses 0-based group indexing:
#   group 0 = receptor (pHLA), group 1 = ligand (Ab)
pHLA_end = pep_count + hla_count

with open('index_mmpbsa.ndx', 'w') as fh:
    _write_ndx_group(fh, 'receptor', range(1, pHLA_end + 1))
    _write_ndx_group(fh, 'ligand', range(pHLA_end + 1, total_prot + 1))
print(f'index_mmpbsa.ndx: receptor(pHLA)={pHLA_end} atoms [P+A], ligand(Ab)={abh_count+abl_count} atoms [H+L]')
print(f'  gmx_MMPBSA -cg 0 1  =>  group 0=receptor(pHLA), group 1=ligand(antibody)')
"""

# Inserted into a QUOTED heredoc (<<'PYEOF') -- no bash expansion needed.
# Keep first 4 molecules (peptide, HLA, Ab-heavy, Ab-light), remove solvent/ions + posre.
_TOPOL_STRIP_PY = r"""
import re, sys
from pathlib import Path

def strip_posre_includes(filepath):
    lines = filepath.read_text(errors='replace').splitlines(keepends=True)
    dropped = 0
    out = []
    for line in lines:
        s = line.strip()
        if re.search(r'#include\s+"posre', s):
            dropped += 1
            sys.stderr.write(f'  drop posre: {filepath.name}: {s}\n')
            continue
        out.append(line)
    if dropped > 0:
        filepath.write_text(''.join(out))
    return dropped

N_KEEP = 4   # peptide, HLA, Ab-heavy, Ab-light

lines = Path('topol.top').read_text().splitlines(keepends=True)
in_mol = False
mol_count = 0
out = []
for line in lines:
    s = line.strip()
    if re.match(r'\[\s*molecules\s*\]', s):
        in_mol = True
        out.append(line)
        continue
    if in_mol:
        if not s or s.startswith(';') or s.startswith('['):
            out.append(line)
            if s.startswith('['):
                in_mol = False
        else:
            mol_count += 1
            if mol_count <= N_KEEP:
                out.append(line)
                sys.stderr.write(f'  keep: {s}\n')
            else:
                sys.stderr.write(f'  drop: {s}\n')
    else:
        out.append(line)

Path('topol_prot.top').write_text(''.join(out))

itp_posre = 0
for itp in sorted(Path('.').glob('topol_Protein_chain_*.itp')):
    itp_posre += strip_posre_includes(itp)

print(f'topol_prot.top: kept {min(mol_count, N_KEEP)}/{mol_count} molecules, removed {itp_posre} posre includes from ITP files')
"""


# \u2500\u2500 SLURM script generation \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500

# Template uses @PLACEHOLDER@ markers \u2014 no f-string escaping issues.
_SLURM_TEMPLATE = r"""#!/bin/bash
#SBATCH --job-name=@JOB_NAME@
#SBATCH --partition=@PARTITION@
#SBATCH --exclude=@NODE_EXCLUDE@
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=@CPUS@
#SBATCH --mem=@MEM@
#SBATCH --time=@TIME@
#SBATCH --output=@OUT_DIR@/slurm_mmpbsa_%j.out
#SBATCH --error=@OUT_DIR@/slurm_mmpbsa_%j.err

# ===========================================================================
# gmx_MMPBSA:  @DESIGN@  replica @REP@
#   MD source group : @SOURCE_GROUP@
#   MD setup dir    : @SETUP_DIR@
#   MMPBSA out dir  : @OUT_DIR@
# ===========================================================================

# \u2500\u2500 Environment \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
module load @GMX_MODULE@
export GMXLIB=@GMXLIB@
export GMX_MAXCONSTRWARN=-1
export OMP_NUM_THREADS=@CPUS@

source @CONDA_BASE@/etc/profile.d/conda.sh
conda activate @CONDA_ENV@ || { echo "ERROR: conda activate @CONDA_ENV@ failed"; exit 1; }
set -euo pipefail

GMX="@GMX@"
DESIGN="@DESIGN@"
REP=@REP@

echo "======================================================================"
echo " MMPBSA (Ab-pHLA): ${DESIGN} r${REP}  node=$(hostname)  job=$SLURM_JOB_ID"
echo " $(date)"
echo "======================================================================"
echo "Python: $(python3 --version)"
gmx_MMPBSA --version 2>&1 | head -1 || true
echo ""

mkdir -p @OUT_DIR@
cd @OUT_DIR@ || { echo "ERROR: cannot cd to @OUT_DIR@"; exit 1; }

# \u2500\u2500 Step 1/7: Copy input files \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
echo "[1/7] Copying topology and trajectory files..."
cp @SETUP_DIR@/topol.top ./                    || { echo "ERROR: missing topol.top"; exit 1; }
cp @SETUP_DIR@/*.itp ./          2>/dev/null   || true
cp @NPT_GRO@         ./npt.gro                 || { echo "ERROR: missing npt.gro";  exit 1; }
cp @TPR@             ./md_r${REP}.tpr        || { echo "ERROR: missing TPR";      exit 1; }
cp @XTC@             ./md_r${REP}.xtc        || { echo "ERROR: missing XTC";      exit 1; }
echo "    OK"

# \u2500\u2500 Step 2/7: Chain-ID assignment (pure Python, no VMD) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
# Reads npt.gro + topol.top + .itp; assigns 4-chain IDs:
#   molecule #1 -> chain P (peptide), #2 -> chain A (HLA),
#   molecule #3 -> chain H (Ab-heavy), #4 -> chain L (Ab-light).
# Outputs: ${DESIGN}_chPAHL.pdb, index_protein.ndx, index_mmpbsa.ndx
echo ""
echo "[2/7] Assigning chain IDs via Python (topol.top + itp)..."
python3 - <<PYEOF || { echo "ERROR: chain assignment failed"; exit 1; }
@CHAIN_ASSIGN_PY@
PYEOF
echo "    OK: ${DESIGN}_chPAHL.pdb  index_protein.ndx  index_mmpbsa.ndx"

# ── Step 3/7: Extract full protein reference structure from npt.gro ────
# Extract all 4 protein chain atoms from npt.gro using index_protein.ndx.
# This GRO file provides coordinates and box vectors for grompp.
echo ""
echo "[3/7] Extracting protein reference structure (Protein.gro)..."
printf "Protein\n" | \
    $GMX trjconv -s md_r${REP}.tpr -f npt.gro -n index_protein.ndx \
    -o Protein.gro > extract_gro.log 2>&1 \
    || { echo "ERROR: failed to extract protein GRO"; cat extract_gro.log; exit 1; }
GRO_ATOMS=$(sed -n '2p' Protein.gro 2>/dev/null | tr -d ' ' || echo "?")
echo "    OK: Protein.gro  ($GRO_ATOMS atoms)"

# ── Step 4/7: Trim topology (4 molecules, no position restraints) ──────
# Removes molecules 5+ (solvent, ions) AND all position-restraint
# includes (posre_*.itp) from ITP files.  Although grompp with define=
# (empty) disables #ifdef-guarded posre blocks, some posre includes are
# unconditional and must be removed from the files.  It also prevents
# grompp from looking for posre coordinate sections that may not exist
# in the protein-only structure.
echo ""
echo "[4/7] Trimming topology -> topol_prot.top (no posre)..."
python3 - <<'PYEOF'
@TOPOL_STRIP_PY@
PYEOF
echo "    OK: topol_prot.top"

# ── Step 5/7: Build protein-only TPR via grompp (no position restraints) ─
# We must use grompp (NOT convert-tpr) because convert-tpr inherits
# position-restraint definitions whose reference coordinates are incomplete
# in the protein subset, causing "Position restraint coordinates are missing"
# inside gmx_MMPBSA.  grompp with define= (empty) + posre-stripped ITP files
# creates a clean TPR with no position restraints.
# -maxwarn 2: the 4-chain protein system may have net charge, triggering
# harmless Ewald warnings; some designs also have +1e on 1 chain.
echo ""
echo "[5/7] Building protein-only TPR via grompp (no posre)..."
cat > grompp_mmpbsa.mdp << 'MDPEOF'
integrator = steep
nsteps = 0
emtol = 1000
define =
cutoff-scheme = Verlet
nstlist = 10
rlist = 1.2
coulombtype = PME
rcoulomb = 1.2
rvdw = 1.2
pbc = xyz
MDPEOF
$GMX grompp -f grompp_mmpbsa.mdp -c Protein.gro -p topol_prot.top \
    -o ${DESIGN}_prot.tpr -maxwarn 2 > grompp_mmpbsa.log 2>&1 \
    || { echo "ERROR: grompp failed"; cat grompp_mmpbsa.log; exit 1; }
echo "    OK: ${DESIGN}_prot.tpr"

# ── Step 6/7: Protein-only PBC-corrected trajectory ──────────────────────
# Use the full-system TPR (md_r${REP}.tpr) as -s for proper centering with
# solvent; the protein-only TPR (design_prot.tpr) lacks solvent molecules.
# trjconv selects all 4 protein chain atoms via -n index_protein.ndx.
echo ""
echo "[6/7] Extracting protein-only trajectory (PBC corrected)..."
printf "Protein\nProtein\n" | \
    $GMX trjconv \
    -f md_r${REP}.xtc \
    -s md_r${REP}.tpr \
    -n index_protein.ndx \
    -o ${DESIGN}_prot_r${REP}.xtc \
    -pbc mol -center > trjconv.log 2>&1 \
    || { echo "ERROR: trjconv failed"; cat trjconv.log; exit 1; }

# GROMACS 2024 trjconv outputs progress sporadically; parse "Last frame" for
# real frame count.  Fall back to file-size heuristic if unavailable.
LAST_FRAME=$(grep -oP 'Last frame\s+\K\d+' trjconv.log 2>/dev/null || echo "")
if [ -n "$LAST_FRAME" ]; then
    NFRAMES=$((LAST_FRAME + 1))
else
    # Approximate frame count from XTC file size.
    # Each compressed XTC frame ≈ 6400 atoms * 4.2 bytes/atom ≈ 26.9 KB.
    XTC_BYTES=$(stat --printf="%s" "${DESIGN}_prot_r${REP}.xtc" 2>/dev/null \
                || stat -f%z "${DESIGN}_prot_r${REP}.xtc" 2>/dev/null || echo "0")
    NFRAMES=$((XTC_BYTES / 27000))
    if [ "$NFRAMES" -lt 5 ]; then NFRAMES="?"; fi
fi
echo "    OK: ${DESIGN}_prot_r${REP}.xtc  (frames: $NFRAMES)"

# ── Step 6b: Validate extracted trajectory ──────────────────────────────────
echo ""
echo "[6b/7] Validating protein trajectory..."
TPR_NATOM=$($GMX dump -s ${DESIGN}_prot.tpr 2>/dev/null | grep -m1 'natoms' | grep -oP '\d+' || echo "")
echo "    TPR atoms: ${TPR_NATOM:-?}"
# XTC atom count: read first frame header via gmx trjconv -dump 0 and
# number of atoms from the PDB or just trust the trjconv step.
# gmx dump on XTC is too slow and produces garbled output, so skip it.
# Instead verify XTC file size is reasonable (> 100 KB).
XTC_SIZE=$(stat --printf="%s" "${DESIGN}_prot_r${REP}.xtc" 2>/dev/null \
           || stat -f%z "${DESIGN}_prot_r${REP}.xtc" 2>/dev/null || echo "0")
echo "    XTC size: ${XTC_SIZE} bytes  (frames: $NFRAMES)"
if [ "$XTC_SIZE" -lt 102400 ] 2>/dev/null; then
    echo "ERROR: ${DESIGN}_prot_r${REP}.xtc is only $XTC_SIZE bytes; trajectory extraction likely failed"
    echo "       Check trjconv.log and verify the source MD trajectory exists and is valid."
    cat trjconv.log
    exit 1
fi

# ── Step 7/7: gmx_MMPBSA ──────────────────────────────────────────────────────
# index_mmpbsa.ndx written by step 2:
#   group 0 = receptor (pHLA: peptide P + HLA A)
#   group 1 = ligand   (antibody: heavy H + light L)
# -cg 0 1: receptor = group 0, ligand = group 1  (0-based NDX indexing)
echo ""
echo "[7/7] Running gmx_MMPBSA (LinearPB)..."
echo "  Files:"
echo "    topol_prot.top  = $(wc -l < topol_prot.top 2>/dev/null || echo '?') lines"
echo "    TPR atoms       = ${TPR_NATOM:-?}"
echo "    NDX groups      = $(grep -c '^\[' index_mmpbsa.ndx 2>/dev/null || echo '?')"
echo "    PDB written     = $(grep -c '^ATOM' ${DESIGN}_chPAHL.pdb 2>/dev/null || echo '?') atoms"
gmx_MMPBSA -O -nogui \
    -i @MMPBSA_IN@ \
    -cp topol_prot.top \
    -cs ${DESIGN}_prot.tpr \
    -ci index_mmpbsa.ndx \
    -cg 0 1 \
    -ct ${DESIGN}_prot_r${REP}.xtc \
    -cr ${DESIGN}_chPAHL.pdb \
    -eo gmx_MMPBSA_plot.csv \
    > mmpbsa.log 2>&1 \
    || { echo "ERROR: gmx_MMPBSA failed. See mmpbsa.log"; tail -60 mmpbsa.log; exit 1; }

[[ ! -f "FINAL_RESULTS_MMPBSA.dat" ]] && \
    { echo "ERROR: FINAL_RESULTS_MMPBSA.dat missing after gmx_MMPBSA"; \
       tail -40 mmpbsa.log; exit 1; }

echo ""
echo "======================================================================"
echo " MMPBSA COMPLETED: ${DESIGN} r${REP}  $(date)"
echo "======================================================================"
grep -E -A5 'ΔTOTAL|DELTA TOTAL' FINAL_RESULTS_MMPBSA.dat | head -10 || true

touch DONE
echo "DONE sentinel written."
"""


def generate_mmpbsa_slurm(
    source_group: str,
    design_name: str,
    rep: int,
    paths: dict,
    out_dir: Path,
    mmpbsa_in: Path,
    gmx: str = "gmx",
    partition: str = DEFAULT_PARTITION,
) -> Path:
    """Generate a self-contained SLURM script for one MMPBSA run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    script_path = out_dir / f"submit_mmpbsa_{design_name}_r{rep}.slurm"

    replacements = {
        "@DESIGN@":          design_name,
        "@REP@":             str(rep),
        "@JOB_NAME@":        slurm_job_name(source_group, design_name, rep),
        "@PARTITION@":       partition,
        "@NODE_EXCLUDE@":    HPC_NODE_EXCLUDE,
        "@CPUS@":            str(HPC_SLURM_CPUS),
        "@MEM@":             HPC_SLURM_MEM,
        "@TIME@":            HPC_SLURM_TIME,
        "@OUT_DIR@":         str(out_dir),
        "@SOURCE_GROUP@":    source_group,
        "@SETUP_DIR@":       str(paths["setup_dir"]),
        "@NPT_GRO@":        str(paths["npt_gro"]),
        "@TPR@":             str(paths["tpr"]),
        "@XTC@":             str(paths["xtc"]),
        "@GMX_MODULE@":      HPC_GMX_MODULE,
        "@GMXLIB@":          HPC_GMXLIB,
        "@CONDA_BASE@":      HPC_CONDA_BASE,
        "@CONDA_ENV@":       HPC_CONDA_ENV,
        "@GMX@":             gmx,
        "@MMPBSA_IN@":       str(mmpbsa_in),
        "@CHAIN_ASSIGN_PY@": _CHAIN_ASSIGN_PY,
        "@TOPOL_STRIP_PY@":  _TOPOL_STRIP_PY,
    }

    script = _SLURM_TEMPLATE
    for marker, value in replacements.items():
        script = script.replace(marker, value)

    script_path.write_text(script)
    script_path.chmod(0o755)
    return script_path

# ── Status report ─────────────────────────────────────────────────────────────
def write_status_report(
    md_results: Path,
    output_root: Path,
    base_dir: Path,
    num_replicas: int,
    targets: set[tuple[str, str]] | None = None,
    targets_file: Path | None = None,
) -> None:
    """Write a Markdown status table to {base_dir}/mmpbsa_status.md."""
    report_path = base_dir / "mmpbsa_status.md"
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    header = ("| Design | Group | "
              + " | ".join(f"r{r}" for r in range(1, num_replicas + 1))
              + " |")
    sep    = "|--------|-------|" + "--------|" * num_replicas

    lines = [
        "# MMPBSA Status Report", "",
        f"Generated  : {now}",
        f"MD results : `{md_results}`",
        f"MMPBSA out : `{output_root}`",
        f"Targets    : `{targets_file}`" if targets_file else "Targets    : all completed MD designs",
        "", header, sep,
    ]

    for grp in sorted(md_results.iterdir()):
        if not grp.is_dir():
            continue
        for design_dir in sorted(grp.iterdir()):
            if not design_dir.is_dir():
                continue
            if targets is not None and (grp.name, design_dir.name) not in targets:
                continue
            row = [design_dir.name, grp.name]
            for rep in range(1, num_replicas + 1):
                out_dir = mmpbsa_output_dir(output_root, grp.name,
                                            design_dir.name, rep)
                if is_mmpbsa_done(out_dir):
                    row.append("OK done")
                elif out_dir.exists():
                    row.append("FAIL" if is_mmpbsa_failed(out_dir) else "run")
                else:
                    md_rep = design_dir / f"replica_{rep}"
                    row.append("pend MD" if not (md_rep / "DONE").exists()
                                else "queued")
            lines.append("| " + " | ".join(row) + " |")

    lines.append("")
    report_path.write_text("\n".join(lines))
    log(f"Status report -> {report_path}")


# ── Main loop ─────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Batch MMPBSA auto-submit manager for Ab-pHLA binding",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--md-results", required=True,
                        help="MD results root (contains replica1/ replica2/ ... subdirs)")
    parser.add_argument("--output",     required=True,
                        help="MMPBSA outputs root directory")
    parser.add_argument("--mmpbsa-in",  required=True,
                        help="Path to gmx_MMPBSA .in input file")
    parser.add_argument("--base-dir",   default=None,
                        help="Base dir for status report and lock (default: --output)")
    parser.add_argument("--gmx",        default="gmx",
                        help="GROMACS executable (default: gmx)")
    parser.add_argument("--targets-file", default=None,
                        help="Optional whitespace-separated file: <outer_replica> <design_name> (one pair per line)")
    parser.add_argument("--max-jobs",   type=int, default=DEFAULT_MAX_JOBS,
                        help=f"Max SLURM jobs in partition (default: {DEFAULT_MAX_JOBS})")
    parser.add_argument("--interval",   type=int, default=DEFAULT_INTERVAL,
                        help=f"Queue poll interval in seconds (default: {DEFAULT_INTERVAL})")
    parser.add_argument("--replicas",   type=int, default=DEFAULT_REPLICAS,
                        help=f"MD replicas per design (default: {DEFAULT_REPLICAS})")
    parser.add_argument("--partition",  default=DEFAULT_PARTITION,
                        help=f"SLURM partition (default: {DEFAULT_PARTITION})")
    parser.add_argument("--dry-run",    action="store_true",
                        help="Generate scripts but do not submit SLURM jobs")
    parser.add_argument("--once",       action="store_true",
                        help="Run one submission cycle then exit (no polling)")
    parser.add_argument("--report",     action="store_true",
                        help="Print status report and exit")
    args = parser.parse_args()

    md_results  = Path(args.md_results).resolve()
    output_root = Path(args.output).resolve()
    mmpbsa_in   = Path(args.mmpbsa_in).resolve()
    base_dir    = Path(args.base_dir).resolve() if args.base_dir else output_root
    targets_file = Path(args.targets_file).resolve() if args.targets_file else None
    targets: set[tuple[str, str]] | None = None

    if not md_results.is_dir():
        log(f"ERROR: --md-results not found: {md_results}"); sys.exit(1)
    if not mmpbsa_in.is_file():
        log(f"ERROR: --mmpbsa-in not found: {mmpbsa_in}");   sys.exit(1)
    if targets_file and not targets_file.is_file():
        log(f"ERROR: --targets-file not found: {targets_file}"); sys.exit(1)

    if targets_file:
        try:
            targets = parse_targets_file(targets_file)
        except ValueError as e:
            log(f"ERROR: {e}"); sys.exit(1)
        log(f"  Targets    : {targets_file}  ({len(targets)} design groups)")

    output_root.mkdir(parents=True, exist_ok=True)
    base_dir.mkdir(parents=True, exist_ok=True)

    if args.report:
        write_status_report(
            md_results, output_root, base_dir, args.replicas,
            targets=targets, targets_file=targets_file,
        )
        return

    _lock = acquire_single_instance_lock(base_dir)  # noqa: F841

    user  = get_current_user()
    state = load_state(base_dir)

    log("=" * 64)
    log("  MMPBSA Auto-Submit Manager")
    log(f"  MD results : {md_results}")
    log(f"  Output     : {output_root}")
    log(f"  MMPBSA in  : {mmpbsa_in}")
    if targets_file:
        log(f"  Targets    : {targets_file}  ({len(targets)} design groups)")
    else:
        log("  Targets    : all completed MD designs")
    log(f"  Partition  : {args.partition}   max_jobs={args.max_jobs}")
    log(f"  Replicas   : {args.replicas}    interval={args.interval}s")
    log(f"  User       : {user}   dry_run={args.dry_run}")
    log("=" * 64)

    cycle = 0
    while True:
        cycle += 1
        log(f"\n-- Cycle {cycle}  ({datetime.now().strftime('%H:%M:%S')}) --")

        current_jobs = count_user_jobs(user, partition=args.partition)
        active_jobs  = get_active_jobs(user, partition=args.partition)
        free_slots   = max(0, args.max_jobs - current_jobs)
        log(f"  Queue [{args.partition}]: {current_jobs}/{args.max_jobs} jobs  "
            f"({free_slots} free)")

        submitted_this_cycle = 0

        for src_grp, design, rep, paths in find_completed_replicas(
                md_results, args.replicas, targets=targets):
            out_dir = mmpbsa_output_dir(output_root, src_grp, design, rep)
            key     = f"{src_grp}/{design}/r{rep}"
            job_name = slurm_job_name(src_grp, design, rep)

            if is_mmpbsa_done(out_dir):
                continue

            if job_name in active_jobs:
                state.setdefault(key, {}).update({
                    "job_id": active_jobs[job_name],
                    "seen_active_at": datetime.now().isoformat(),
                })
                continue

            if has_submission_artifacts(out_dir):
                if not is_mmpbsa_failed(out_dir):
                    log(f"  STALE sub: {key} has artifacts but no active job; resubmitting")
                retries = state.get(key, {}).get("retries", 0)
                if retries >= MAX_RETRIES:
                    log(f"  GAVE UP ({retries} retries): {key}")
                    continue
                log(f"  RETRY ({retries + 1}/{MAX_RETRIES}): {key}")
                state.setdefault(key, {})["retries"] = retries + 1

            if free_slots <= 0:
                log(f"  Queue full, deferring: {key}")
                continue

            script = generate_mmpbsa_slurm(
                src_grp, design, rep, paths, out_dir,
                mmpbsa_in, gmx=args.gmx, partition=args.partition,
            )
            log(f"  -> Submitting: {key}")
            job_id = submit_slurm(script, dry_run=args.dry_run)
            if job_id:
                state.setdefault(key, {}).update({
                    "job_id":        job_id,
                    "submitted_at":  datetime.now().isoformat(),
                    "job_name":      job_name,
                })
                free_slots -= 1
                submitted_this_cycle += 1
                log(f"     job {job_id}: {key}")

        save_state(base_dir, state)
        write_status_report(
            md_results, output_root, base_dir, args.replicas,
            targets=targets, targets_file=targets_file,
        )
        log(f"  Submitted {submitted_this_cycle} job(s) this cycle.")

        if args.once or args.dry_run:
            log("Exiting after one cycle (--once / --dry-run).")
            break

        log(f"  Sleeping {args.interval}s  (Ctrl+C to stop)...")
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
