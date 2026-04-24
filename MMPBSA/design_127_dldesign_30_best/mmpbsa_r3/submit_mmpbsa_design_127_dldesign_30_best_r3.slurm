#!/bin/bash
#SBATCH --job-name=mmpbsa_replica3_design_127_dldesign_30_best_r3
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node15,node20,node21,node24,node26
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=/public/home/xuziyi/MMPBSA/results_top15/replica3/design_127_dldesign_30_best/mmpbsa_r3/slurm_mmpbsa_%j.out
#SBATCH --error=/public/home/xuziyi/MMPBSA/results_top15/replica3/design_127_dldesign_30_best/mmpbsa_r3/slurm_mmpbsa_%j.err

# ===========================================================================
# gmx_MMPBSA:  design_127_dldesign_30_best  replica 3
#   MD source group : replica3
#   MD setup dir    : /public/home/xuziyi/MD/results/replica3/design_127_dldesign_30_best/setup
#   MMPBSA out dir  : /public/home/xuziyi/MMPBSA/results_top15/replica3/design_127_dldesign_30_best/mmpbsa_r3
# ===========================================================================

# \u2500\u2500 Environment \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
module load gromacs/2024.2
export GMXLIB=/public/home/xuziyi/FEP/force_fields/mutff
export GMX_MAXCONSTRWARN=-1
export OMP_NUM_THREADS=16

source /public/home/xuziyi/miniconda/etc/profile.d/conda.sh
conda activate gmxMMPBSA || { echo "ERROR: conda activate gmxMMPBSA failed"; exit 1; }
set -euo pipefail

GMX="gmx"
DESIGN="design_127_dldesign_30_best"
REP=3

echo "======================================================================"
echo " MMPBSA (Ab-pHLA): ${DESIGN} r${REP}  node=$(hostname)  job=$SLURM_JOB_ID"
echo " $(date)"
echo "======================================================================"
echo "Python: $(python3 --version)"
gmx_MMPBSA --version 2>&1 | head -1 || true
echo ""

mkdir -p /public/home/xuziyi/MMPBSA/results_top15/replica3/design_127_dldesign_30_best/mmpbsa_r3
cd /public/home/xuziyi/MMPBSA/results_top15/replica3/design_127_dldesign_30_best/mmpbsa_r3 || { echo "ERROR: cannot cd to /public/home/xuziyi/MMPBSA/results_top15/replica3/design_127_dldesign_30_best/mmpbsa_r3"; exit 1; }

# \u2500\u2500 Step 1/7: Copy input files \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
echo "[1/7] Copying topology and trajectory files..."
cp /public/home/xuziyi/MD/results/replica3/design_127_dldesign_30_best/setup/topol.top ./                    || { echo "ERROR: missing topol.top"; exit 1; }
cp /public/home/xuziyi/MD/results/replica3/design_127_dldesign_30_best/setup/*.itp ./          2>/dev/null   || true
cp /public/home/xuziyi/MD/results/replica3/design_127_dldesign_30_best/setup/npt.gro         ./npt.gro                 || { echo "ERROR: missing npt.gro";  exit 1; }
cp /public/home/xuziyi/MD/results/replica3/design_127_dldesign_30_best/replica_3/md_r3.tpr             ./md_r${REP}.tpr        || { echo "ERROR: missing TPR";      exit 1; }
cp /public/home/xuziyi/MD/results/replica3/design_127_dldesign_30_best/replica_3/md_r3.xtc             ./md_r${REP}.xtc        || { echo "ERROR: missing XTC";      exit 1; }
echo "    OK"

# \u2500\u2500 Step 2/7: Chain-ID assignment (pure Python, no VMD) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
# Reads npt.gro + topol.top + .itp; assigns 4-chain IDs:
#   molecule #1 -> chain P (peptide), #2 -> chain A (HLA),
#   molecule #3 -> chain H (Ab-heavy), #4 -> chain L (Ab-light).
# Outputs: ${DESIGN}_chPAHL.pdb, index_protein.ndx, index_mmpbsa.ndx
echo ""
echo "[2/7] Assigning chain IDs via Python (topol.top + itp)..."
python3 - <<PYEOF || { echo "ERROR: chain assignment failed"; exit 1; }

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


SOLVENT = frozenset([
    'SOL', 'TIP3', 'TIP3P', 'WAT', 'HOH',
    'NA', 'NA+', 'SOD', 'CL', 'CL-', 'CLA',
    'K', 'K+', 'MG', 'CA', 'ZN',
])


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
    if resname.upper() in SOLVENT:
        continue
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
    -i /public/home/xuziyi/MMPBSA/inputs/mmpbsa_LinearPB.in \
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
