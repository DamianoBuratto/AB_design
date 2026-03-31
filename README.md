# AB_design

A set of scripts to perform Ab initio antibody design using RFantibody integrated with AlphaFold3 and MD simulations using Gromacs.

---

## Pipeline Workflow

```
Stage 1: DESIGN (RFantibody)         
  RFdiffusion → Filter1 → ProteinMPNN → RF2 → Analysis
  Output: 20-100 high-quality antibody designs

Stage 2: VALIDATION (AF3)          
  Filter RF2 → Prepare AF3 → Run AF3 → RMSD Check
  Output: Validated Ab-pHLA complex structures

Stage 3: DYNAMICS (MD)              
  Setup → Equilibration → Production (3×200ns)
  Output: Equilibrated trajectories

Stage 4: AFFINITY (FEP)         
  Setup → Lambda Windows (32*4) → Analysis
  Output: ΔΔG_binding
```


---
## Module 1: RFantibody - Antibody Design

### Core Scripts

| Script | Function |
|--------|----------|
| `run_hotspot_batch_functional_split.sh` | Main pipeline launcher (2-stage SLURM,Dependency) |
| `run_config.txt` | Configuration file |
| `count_contacts_with_hotspot.py` | Filter designs by hotspot contacts |
| `analyze_rf2_pdb.py` | Quality assessment (pLDDT, RMSD, contacts) |

### Command Flow

```bash
cd RFantibody

# 1. Configure
vim run_config.txt

# 2. Launch pipeline
bash run_hotspot_batch_functional_split.sh
```
---

## Module 2: AF3 - Structure Validation

### Core Scripts

| Script | Function |
|--------|----------|
| `filter_and_copy_rf2_best.py` | Select top RF2 designs |
| `batch_prepare_af3.py` | Generate AF3 input JSON files |
| `smart_queue_manager.sh` | Intelligent job submission |
| `alphafold3.sh` | AF3 execution wrapper |
| `compute_rmsd.py` | RMSD calculation |

### Command Flow

```bash
cd AF3

# 1. Filter best designs
python3 filter_and_copy_rf2_best.py \
    --rf-output-dir ../RFantibody/outputs/[replica]

# 2. Prepare AF3 jobs
python3 batch_prepare_af3.py \
    --input-dir RFantibody_best \
    --templates-dir templates \
    --output-base jobs \
    --mut-only

# 3. Submit jobs
screen -S af3
bash smart_queue_manager.sh jobs ./alphafold3.sh

# 4. Compute RMSD (after completion)
bash batch_compute_rmsd.sh
```
---

## Module 3: MD - Molecular Dynamics

### Core Scripts

| Script | Function |
|--------|----------|
| `scripts_charmm/auto_submit_charmm.py` | Main submit tool: auto-manages SLURM queue |
| `scripts_charmm/setup_system_charmm.py` | System setup: chain split, renumber, pdb2gmx, solvate, minimize |
| `scripts_charmm/run_equilibration_charmm.py` | NVT + NPT equilibration |
| `scripts_charmm/run_production_charmm.py` | Production MD (3 replicas × 150 ns each) |
| `scripts_charmm/batch_submit_charmm.py` | Generate SLURM scripts only (no submission) |

### Command Flow

```bash
cd MD

# 1. Place PDB in data/

# 2. Generate workflow & auto_submit
screen -S md_run

python3 scripts_charmm/auto_submit_charmm.py \
    --output ./results \
    --mdp    ./mdp_files_charmm \
    --base-dir /public/home/xuziyi/MD \
    >> auto_submit.log 2>&1

Ctrl+A, D 
```

---

## Module 4: FEP - Free Energy Perturbation (TESTING)

### Core Scripts

| Script | Function |
|--------|----------|
| `scripts/setup_fep.sh` | Full setup: GRO→PDB, pmx mutate, pdb2gmx, solvate, genion, replica copy, reverse dirs |
| `scripts/submit_fep.sh` | Batch-generate and submit SLURM jobs|
| `scripts/analyze_fep.py` | MBAR post-processing → ΔΔG_binding with overlap matrices and convergence plots |
| `mdps/em_fep.mdp` | Energy minimization MDP (steep→CG, pmx hybrid topology) |
| `mdps/nvt_fep.mdp` | NVT equilibration MDP (300 ps, V-rescale, 310 K) |
| `mdps/npt_fep.mdp` | NPT equilibration MDP (500 ps, Berendsen) |
| `mdps/prod_fep.mdp` | Production HREX MDP (1.5 ns, 32 λ windows, C-rescale) |

### Command Flow

```bash
cd FEP

# --- FORWARD ---
# 1. Setup (forward only, ~5-10 min on HPC)
#    Runs: GRO->PDB, pmx mutate, pdb2gmx, solvate, genion
#    Produces: outputs/<DESIGN>/forward/bound_R1~R3 + unbound_R1~R3
bash scripts/setup_fep.sh

# 2. Submit forward jobs (6 SLURM jobs: bound/unbound x R1-R3)
bash scripts/submit_fep.sh

# --- REVERSE (after forward PROD is complete) ---
# 3. Manually extract last frame of lambda31 (ARG end-state) from each forward run.
#    BOUND and UNBOUND have different atom counts - must be extracted separately.
DESIGN=design_373_dldesign_13_best
echo "0" | gmx trjconv \
    -f outputs/${DESIGN}/forward/bound_R1/lambda31/PROD/prod.xtc \
    -s outputs/${DESIGN}/forward/bound_R1/lambda31/PROD/prod.tpr \
    -o rev_start_bound.gro -dump -1

echo "0" | gmx trjconv \
    -f outputs/${DESIGN}/forward/unbound_R1/lambda31/PROD/prod.xtc \
    -s outputs/${DESIGN}/forward/unbound_R1/lambda31/PROD/prod.tpr \
    -o rev_start_unbound.gro -dump -1

# 4. Setup reverse: reuses forward hybrid topology; only starting GRO differs.
#    --rev-gro and --rev-gro-unbound are required to place starting structures.
bash scripts/setup_fep.sh --direction reverse \
    --rev-gro  rev_start_bound.gro \
    --rev-gro-unbound rev_start_unbound.gro

# 5. Submit reverse jobs (6 SLURM jobs)
bash scripts/submit_fep.sh reverse

# --- ANALYSIS (after all PROD complete: forward + reverse) ---
# 6. MBAR analysis -> ddG_results.txt + overlap/convergence plots
conda activate fep_env
python scripts/analyze_fep.py outputs/${DESIGN}
```

