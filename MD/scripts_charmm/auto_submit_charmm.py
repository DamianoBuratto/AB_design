#!/usr/bin/env python3
"""
Auto-submit manager for pHLA MD (CHARMM36m).

Workflow:
  1. Scan data/replica*/ for all PDBs to run.
  2. Submit up to MAX_JOBS / JOBS_PER_SYSTEM systems at a time.
  3. Each cycle: detect failures in already-submitted designs, retry once if possible.
  4. Poll the SLURM queue every POLL_INTERVAL seconds; submit more as slots free up.
  5. Write base_dir/md_status.md after every cycle.

Success detection (sentinel-file approach):
  - SLURM scripts write  setup/DONE  and  replica_N/DONE  on success.
  - Python checks: DONE present + log ends with success marker + XTC exists (replicas).
  - If DONE is absent after the job has left the queue → stage is "failed".

Retry logic:
  - Setup fail   → re-submit setup + all replicas  (max MAX_RETRIES times)
  - Replica fail → re-submit that replica only      (max MAX_RETRIES times)
  - After MAX_RETRIES: logged as "gave up"; node name recorded in .job_state.json
    and shown in md_status.md for diagnosis.

Directory conventions (all relative to --base-dir, default = script's parent):
  data/               All PDBs to run, organised into replica subdirs:
                        data/replica1/  data/replica2/  data/replica3/ …
  data_pre/           Archived PDBs — do not submit, do not touch their results.
  results/            Output root; mirrors the data/ subdir structure:
                        data/replica1/design.pdb → results/replica1/design/
  md_status.md        Real-time Markdown status report (written to base_dir/)

Usage:
  python3 auto_submit_charmm.py \\
      --output   ./results \\
      --mdp      ./mdp_files_charmm \\
      [--base-dir  /public/home/xuziyi/MD] \\
      [--max-jobs  16] \\
      [--interval  120] \\
      [--dry-run]
"""

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# ── tunables ────────────────────────────────────────────────────────────────
JOBS_PER_SYSTEM  = 4          # 1 setup + 3 replicas
DEFAULT_MAX_JOBS = 16         # QOS limit
DEFAULT_INTERVAL = 120        # seconds between queue polls
DEFAULT_REPLICAS = 3
# ────────────────────────────────────────────────────────────────────────────


def log(msg: str):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def count_user_jobs(user: str) -> int:
    """Return how many jobs (any state) the user has in the 'multi' partition.

    Only the multi partition is counted so that FEP jobs (quick partition)
    do not consume MD slot budget.
    """
    try:
        result = subprocess.run(
            ["squeue", "-u", user, "-p", "multi", "-h"],
            capture_output=True, text=True, timeout=15
        )
        lines = [l for l in result.stdout.splitlines() if l.strip()]
        return len(lines)
    except Exception as e:
        log(f"WARNING: squeue failed ({e}); assuming 0 jobs")
        return 0


def get_current_user() -> str:
    return os.environ.get("USER", subprocess.run(
        ["whoami"], capture_output=True, text=True).stdout.strip())


def find_pending_pdbs(data_dir: Path) -> list[Path]:
    """All PDB files in data/ (any depth) sorted by relative path."""
    return sorted(data_dir.rglob("*.pdb"))


def output_dir_for(pdb: Path, data_dir: Path, output_root: Path) -> Path:
    """Mirror the PDB's subdir position within data/ into the output root."""
    rel    = pdb.relative_to(data_dir)
    subdir = rel.parent
    return output_root / subdir / pdb.stem


def is_already_submitted(pdb: Path, data_dir: Path, output_root: Path) -> bool:
    """A system is considered submitted if its output directory exists."""
    return output_dir_for(pdb, data_dir, output_root).exists()


# ── SLURM script generation (mirrors batch_submit_charmm.py) ─────────────

def generate_setup_slurm(pdb_name: str, setup_dir: Path,
                         pdb_file: Path, mdp_dir: Path,
                         script_dir: Path, gmx: str) -> Path:
    script = f"""#!/bin/bash
#SBATCH --job-name=setup_{pdb_name}
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node20,node21,node24,node26
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=3-00:00:00
#SBATCH --output={setup_dir}/slurm_setup_%j.out
#SBATCH --error={setup_dir}/slurm_setup_%j.err

CONDA_BASE=/public/home/xuziyi/miniconda
CONDA_ENV=pHLA_MD_env
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate $CONDA_ENV || {{ echo "ERROR: conda activate failed"; exit 1; }}
echo "Python: $(python3 --version)"

module load gromacs/2024.2
module load cuda
export GMXLIB=/public/home/xuziyi/FEP/force_fields/mutff
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_MAXCONSTRWARN=-1

GMX="{gmx}"
MDP_DIR="{mdp_dir}"
SCRIPT_DIR="{script_dir}"
PDB_FILE="{pdb_file}"

# Record execution context
NODE=$(hostname)
JOBID=$SLURM_JOB_ID
echo "node=$NODE" >> {setup_dir}/run_info.txt
echo "setup_jobid=$JOBID" >> {setup_dir}/run_info.txt
echo "setup_started=$(date -Iseconds)" >> {setup_dir}/run_info.txt

echo "======================================"
echo "Setup: {pdb_name} (CHARMM36m) on $NODE"
echo "======================================"

cd {setup_dir}
python3 $SCRIPT_DIR/setup_system_charmm.py \\
    -i "$PDB_FILE" -w {setup_dir} \\
    -ions "$MDP_DIR/ions.mdp" -em "$MDP_DIR/minim.mdp" \\
    -gmx "$GMX" -ff charmm36m-mut -water tip3p \\
    -box triclinic -d 1.2 -conc 0.10 || {{ echo "ERROR: setup failed"; exit 1; }}

python3 $SCRIPT_DIR/run_equilibration_charmm.py \\
    -w {setup_dir} -nvt "$MDP_DIR/nvt.mdp" -npt "$MDP_DIR/npt.mdp" \\
    -gmx "$GMX" || {{ echo "ERROR: equilibration failed"; exit 1; }}

echo "setup_finished=$(date -Iseconds)" >> {setup_dir}/run_info.txt
touch {setup_dir}/DONE
echo "Setup + equilibration done for {pdb_name}"
"""
    out = setup_dir / f"submit_charmm_{pdb_name}_setup.slurm"
    out.write_text(script)
    out.chmod(0o755)
    return out


def generate_replica_slurm(pdb_name: str, setup_dir: Path, replica_dir: Path,
                            rep: int, mdp_dir: Path, script_dir: Path,
                            gmx: str) -> Path:
    script = f"""#!/bin/bash
#SBATCH --job-name=prod_{pdb_name}_r{rep}
#SBATCH --partition=multi
#SBATCH --exclude=node1,node2,node3,node4,node5,node6,node7,node8,node9,node10,node11,node12,node13,node14,node20,node21,node24,node26
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output={replica_dir}/slurm_r{rep}_%j.out
#SBATCH --error={replica_dir}/slurm_r{rep}_%j.err

CONDA_BASE=/public/home/xuziyi/miniconda
CONDA_ENV=pHLA_MD_env
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate $CONDA_ENV || {{ echo "ERROR: conda activate failed"; exit 1; }}

module load gromacs/2024.2
module load cuda
export GMXLIB=/public/home/xuziyi/FEP/force_fields/mutff
export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_MAXCONSTRWARN=-1

GMX="{gmx}"
MDP_DIR="{mdp_dir}"
SCRIPT_DIR="{script_dir}"

# Record execution context
NODE=$(hostname)
JOBID=$SLURM_JOB_ID
echo "r{rep}_node=$NODE" >> {setup_dir}/run_info.txt
echo "r{rep}_jobid=$JOBID" >> {setup_dir}/run_info.txt
echo "r{rep}_started=$(date -Iseconds)" >> {setup_dir}/run_info.txt

[ ! -f "{setup_dir}/npt.gro" ] && \\
    {{ echo "ERROR: equilibration not complete (npt.gro missing)"; exit 1; }}

cd {replica_dir}
python3 $SCRIPT_DIR/run_production_charmm.py \\
    -w {replica_dir} -md "$MDP_DIR/production.mdp" \\
    -r {rep} -setup {setup_dir} -time 150 \\
    -gmx "$GMX" --no-analysis || {{ echo "ERROR: production failed"; exit 1; }}

echo "r{rep}_finished=$(date -Iseconds)" >> {setup_dir}/run_info.txt
touch {replica_dir}/DONE
echo "Replica {rep} done for {pdb_name}"
"""
    out = setup_dir.parent / f"submit_charmm_{pdb_name}_r{rep}.slurm"
    out.write_text(script)
    out.chmod(0o755)
    return out


def submit_system(pdb: Path, pdb_dir: Path, data_dir: Path,
                  mdp_dir: Path, script_dir: Path, gmx: str,
                  num_replicas: int, dry_run: bool,
                  retry_setup: bool = False,
                  retry_replicas: list[int] | None = None) -> bool:
    """
    Generate scripts + submit jobs for one PDB.

    Normal first-time submission: retry_setup=False, retry_replicas=None
      → submits setup + all replicas (replicas depend on setup via afterok)

    Retry setup: retry_setup=True
      → clears DONE, resubmits setup; replicas re-queued with afterok dependency

    Retry individual replicas: retry_replicas=[1,3]
      → only resubmits those replicas; setup already done, no dependency needed

    pdb_dir: output directory for this design (output_root/replica_N/design_name)
    Returns True on success.
    """
    name = pdb.stem

    if dry_run:
        action = "setup+replicas"
        if retry_setup:
            action = "RETRY setup+replicas"
        elif retry_replicas:
            action = f"RETRY replicas {retry_replicas}"
        log(f"[DRY-RUN] Would submit {name} ({action})")
        return True

    setup_dir = pdb_dir / "setup"
    setup_dir.mkdir(parents=True, exist_ok=True)

    replica_dirs = []
    for i in range(1, num_replicas + 1):
        rd = pdb_dir / f"replica_{i}"
        rd.mkdir(parents=True, exist_ok=True)
        replica_dirs.append(rd)

    setup_slurm = generate_setup_slurm(name, setup_dir, pdb,
                                       mdp_dir, script_dir, gmx)
    replica_slurms = [
        generate_replica_slurm(name, setup_dir, replica_dirs[i - 1],
                               i, mdp_dir, script_dir, gmx)
        for i in range(1, num_replicas + 1)
    ]

    # ── retry setup: clear old DONE so sentinel logic works correctly ──────
    if retry_setup:
        done_flag = setup_dir / "DONE"
        if done_flag.exists():
            done_flag.unlink()
        log(f"  Retrying setup for {name}")

    # ── submit setup (unless we're only retrying specific replicas) ────────
    if not retry_replicas:
        try:
            result = subprocess.run(
                ["sbatch", "--parsable", str(setup_slurm)],
                capture_output=True, text=True, check=True, timeout=30
            )
            setup_job_id = result.stdout.strip()
            log(f"  Setup job submitted: {setup_job_id}  ({name})")
        except subprocess.CalledProcessError as e:
            log(f"  ERROR: Setup submission failed for {name}: {e.stderr.strip()}")
            return False
    else:
        setup_job_id = None  # replicas run without setup dependency

    # ── decide which replicas to submit ───────────────────────────────────
    reps_to_submit = retry_replicas if retry_replicas else list(range(1, num_replicas + 1))

    for i in reps_to_submit:
        rslurm = replica_slurms[i - 1]
        # Clear old DONE for retried replicas
        done_flag = pdb_dir / f"replica_{i}" / "DONE"
        if done_flag.exists():
            done_flag.unlink()
        try:
            sbatch_cmd = ["sbatch", "--parsable"]
            if setup_job_id:
                sbatch_cmd.append(f"--dependency=afterok:{setup_job_id}")
            sbatch_cmd.append(str(rslurm))
            result = subprocess.run(
                sbatch_cmd,
                capture_output=True, text=True, check=True, timeout=30
            )
            log(f"  Replica {i} job submitted: {result.stdout.strip()}  ({name})")
        except subprocess.CalledProcessError as e:
            log(f"  ERROR: Replica {i} submission failed for {name}: {e.stderr.strip()}")

    return True


# ── status & completion helpers ──────────────────────────────────────────────

MAX_RETRIES = 1   # how many times a failed stage can be re-submitted

def _log_contains(path: Path, marker: str) -> bool:
    """Return True if *path* exists and any non-empty line contains *marker*."""
    if not path.exists():
        return False
    try:
        lines = [l for l in path.read_text(errors="replace").splitlines() if l.strip()]
        return any(marker in line for line in lines)
    except OSError:
        return False


def check_setup_complete(setup_dir: Path) -> str:
    """
    Inspect the setup directory.

    Returns one of:
      'complete'  DONE sentinel present (implying both setup+equil logs succeeded)
      'running'   directory exists, no DONE, no error marker in logs
      'failed'    explicit error marker in logs (✗) OR npt.gro missing but
                  equilibration.log exists and doesn't end in success
      'missing'   directory does not exist

    Note: we cannot distinguish "still running" from "silently crashed" without
    sacct.  We conservatively call it 'running' unless there is a positive error
    signal from the log.
    """
    if not setup_dir.exists():
        return "missing"
    if (setup_dir / "DONE").exists():
        return "complete"
    # Check for explicit failure markers in logs
    for logname, fail_marker in [
        ("setup_system.log",  "✗"),
        ("equilibration.log", "✗"),
        ("setup_system.log",  "SYSTEM SETUP FAILED"),
        ("equilibration.log", "EQUILIBRATION FAILED"),
    ]:
        if _log_contains(setup_dir / logname, fail_marker):
            return "failed"
    return "running"


def check_replica_complete(replica_dir: Path, rep: int) -> str:
    """
    Inspect one production replica directory.

    Returns one of:
      'complete'  DONE sentinel present AND at least one XTC file exists
      'running'   directory exists, no DONE yet, no error marker in log
      'failed'    explicit error marker in production log (✗ or FAILED)
                  OR DONE present but no XTC (interrupted after sentinel write)
      'missing'   directory does not exist
    """
    if not replica_dir.exists():
        return "missing"
    if (replica_dir / "DONE").exists():
        # Belt-and-suspenders: also verify XTC file exists
        if list(replica_dir.glob("*.xtc")):
            return "complete"
        # DONE written but no XTC — treat as failed
        return "failed"
    prod_log = replica_dir / f"production_r{rep}.log"
    for marker in ["✗", "PRODUCTION MD FAILED", "ERROR: production failed"]:
        if _log_contains(prod_log, marker):
            return "failed"
    return "running"


# ── per-design state file (.job_state.json) ──────────────────────────────────

def _state_path(pdb_dir: Path) -> Path:
    return pdb_dir / ".job_state.json"


def load_state(pdb_dir: Path) -> dict:
    """
    Load or initialise the per-design state file.

    Schema:
      setup_retries: int
      replica_retries: {1: int, 2: int, 3: int}
      failed_nodes: list[str]   # nodes recorded when a job failed
    """
    p = _state_path(pdb_dir)
    if p.exists():
        try:
            return json.loads(p.read_text())
        except Exception:
            pass
    return {"setup_retries": 0, "replica_retries": {}, "failed_nodes": []}


def save_state(pdb_dir: Path, state: dict) -> None:
    _state_path(pdb_dir).write_text(json.dumps(state, indent=2))


def read_node_from_run_info(setup_dir: Path, key: str) -> str | None:
    """Read a key=value pair from setup_dir/run_info.txt."""
    info = setup_dir / "run_info.txt"
    if not info.exists():
        return None
    for line in info.read_text().splitlines():
        if line.startswith(f"{key}="):
            return line.split("=", 1)[1].strip()
    return None


def generate_status_report(
    output_root: Path,
    base_dir: Path,
    num_replicas: int,
) -> None:
    """
    Scan data/ and output_root; write a Markdown status table to
    base_dir/md_status.md (human-readable, real-time updated).

    Columns per design:
      Source | Design | Setup | R1 | R2 | R3 | Failed nodes | Notes

    data_pre/ entries are listed as archived (read-only, no result check).
    """
    report_file = base_dir / "md_status.md"
    data_dir    = base_dir / "data"
    datapre_dir = base_dir / "data_pre"

    # ---------- Collect all active designs ----------
    # key: (source_label, design_name)  value: pdb_dir Path or None
    designs: dict[tuple[str, str], Path | None] = {}

    def _register_data(pdb: Path) -> None:
        rel    = pdb.relative_to(data_dir)
        subdir = str(rel.parent) if rel.parent != Path(".") else ""
        source = f"data/{subdir}" if subdir else "data"
        outdir = (output_root / subdir / pdb.stem) if subdir else (output_root / pdb.stem)
        designs.setdefault((source, pdb.stem), outdir)

    if data_dir.exists():
        for pdb in sorted(data_dir.rglob("*.pdb")):
            _register_data(pdb)

    # Catch designs whose output dir exists but PDB may have been moved
    if output_root.exists():
        for setup_dir in sorted(output_root.rglob("setup")):
            if not setup_dir.is_dir():
                continue
            design_dir = setup_dir.parent
            rel        = design_dir.relative_to(output_root)
            subdir     = str(rel.parent) if rel.parent != Path(".") else ""
            source     = f"data/{subdir}" if subdir else "data"
            designs.setdefault((source, design_dir.name), design_dir)

    # ---------- data_pre (archived) ----------
    pre_entries: list[tuple[str, str]] = []
    if datapre_dir.exists():
        for pdb in sorted(datapre_dir.rglob("*.pdb")):
            rel    = pdb.relative_to(datapre_dir)
            subdir = str(rel.parent) if rel.parent != Path(".") else ""
            source = f"data_pre/{subdir}" if subdir else "data_pre"
            pre_entries.append((source, pdb.stem))

    # ---------- Status symbols ----------
    SYM = {
        "complete": "✅ done",
        "running":  "🔄 run",
        "failed":   "❌ FAIL",
        "missing":  "—",
    }

    now   = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines = [
        f"# pHLA MD Status  —  {now}",
        f"",
        f"Base dir : `{base_dir}`  ",
        f"Results  : `{output_root}`  ",
        f"",
    ]

    n_sub = n_done = n_fail = 0
    prev_source = None

    rep_headers = " | ".join(f"R{i}" for i in range(1, num_replicas + 1))
    header = f"| Source | Design | Setup | {rep_headers} | Failed nodes | Notes |"
    sep    = "|--------|--------|-------|" + "-------|" * num_replicas + "--------------|-------|"

    def _md_row(source, name, setup_sym, rep_syms, nodes, notes):
        rep_cols = " | ".join(rep_syms)
        return f"| {source} | {name} | {setup_sym} | {rep_cols} | {nodes} | {notes} |"

    for (source, name) in sorted(designs):
        if source != prev_source:
            lines += ["", f"## {source}", "", header, sep]
            prev_source = source

        pdb_dir = designs[(source, name)]
        submitted = pdb_dir is not None and (pdb_dir / "setup").is_dir()

        if not submitted:
            lines.append(_md_row(source, name, "—",
                                 ["—"] * num_replicas, "", "not submitted"))
            continue

        n_sub += 1
        setup_dir    = pdb_dir / "setup"
        setup_status = check_setup_complete(setup_dir)
        setup_sym    = SYM[setup_status]

        # Load retry state
        state       = load_state(pdb_dir)
        failed_nodes = state.get("failed_nodes", [])
        nodes_str    = ", ".join(sorted(set(failed_nodes))) if failed_nodes else ""

        rep_syms  = []
        all_ok    = (setup_status == "complete")
        any_fail  = (setup_status == "failed")
        notes_parts = []

        if setup_status == "failed":
            s_retries = state.get("setup_retries", 0)
            if s_retries >= MAX_RETRIES:
                notes_parts.append(f"setup gave up (retried {s_retries}x)")

        for i in range(1, num_replicas + 1):
            rdir   = pdb_dir / f"replica_{i}"
            status = check_replica_complete(rdir, i)
            rep_syms.append(SYM[status])
            if status != "complete":
                all_ok = False
            if status == "failed":
                any_fail = True
                r_retries = state.get("replica_retries", {}).get(str(i), 0)
                if r_retries >= MAX_RETRIES:
                    notes_parts.append(f"r{i} gave up (retried {r_retries}x)")

        if all_ok:
            n_done += 1
        if any_fail:
            n_fail += 1

        lines.append(_md_row(source, name, setup_sym, rep_syms,
                             nodes_str, "; ".join(notes_parts)))

    # ---------- Archived section ----------
    if pre_entries:
        lines += ["", "## Archived (data_pre)", "", header, sep]
        for (source, name) in sorted(pre_entries):
            lines.append(_md_row(source, name, "archived",
                                 ["n/a"] * num_replicas, "", ""))

    # ---------- Summary ----------
    total = len(designs)
    lines += [
        "",
        "---",
        f"**Active designs: {total}** | Submitted: {n_sub} | "
        f"All replicas done: {n_done} | Has failures: {n_fail} | "
        f"Archived: {len(pre_entries)}",
        "",
        "_Status key: ✅ done = log+XTC confirmed | "
        "🔄 run = in-progress | ❌ FAIL = log error / DONE missing after job | "
        "— = not started_",
    ]

    report_file.write_text("\n".join(lines) + "\n")
    log(f"Status report → {report_file}")


# ── main loop ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Auto-submit manager for pHLA MD (CHARMM36m)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--output",   required=True,
                        help="Output root directory for simulation results")
    parser.add_argument("--mdp",      required=True,
                        help="Directory containing MDP files")
    parser.add_argument("--base-dir", default=None,
                        help="Base directory (contains data/, data_pre/). "
                             "Default: parent of this script.")
    parser.add_argument("--gmx",      default="gmx",
                        help="GROMACS executable (default: gmx)")
    parser.add_argument("--max-jobs", type=int, default=DEFAULT_MAX_JOBS,
                        help=f"QOS job limit (default: {DEFAULT_MAX_JOBS})")
    parser.add_argument("--replicas", type=int, default=DEFAULT_REPLICAS,
                        help=f"Replicas per system (default: {DEFAULT_REPLICAS})")
    parser.add_argument("--interval", type=int, default=DEFAULT_INTERVAL,
                        help=f"Polling interval in seconds (default: {DEFAULT_INTERVAL})")
    parser.add_argument("--once",     action="store_true",
                        help="Submit one batch then exit (no polling loop)")
    parser.add_argument("--dry-run",  action="store_true",
                        help="Print what would be submitted without actually submitting")
    parser.add_argument("--report",   action="store_true",
                        help="Generate a status report and exit (no submission)")
    args = parser.parse_args()

    script_dir  = Path(__file__).parent.resolve()
    base_dir    = Path(args.base_dir).resolve() if args.base_dir else script_dir.parent
    data_dir    = base_dir / "data"
    output_root = Path(args.output).resolve()
    mdp_dir     = Path(args.mdp).resolve()
    user        = get_current_user()

    jobs_per = 1 + args.replicas          # setup + replicas
    max_systems = args.max_jobs // jobs_per

    # ── report-only mode ────────────────────────────────────────────────────
    if args.report:
        generate_status_report(output_root, base_dir, args.replicas)
        return

    log("=" * 60)
    log("pHLA MD Auto-submit Manager")
    log(f"  base_dir    : {base_dir}")
    log(f"  data/       : {data_dir}")
    log(f"  output      : {output_root}")
    log(f"  max_jobs    : {args.max_jobs}  →  max {max_systems} systems at once")
    log(f"  poll every  : {args.interval}s")
    log("=" * 60)

    for d in [data_dir, output_root, mdp_dir]:
        if not d.exists():
            log(f"ERROR: Required directory not found: {d}")
            sys.exit(1)

    while True:
        current_jobs = count_user_jobs(user)
        available_slots = max(0, args.max_jobs - current_jobs)
        can_submit = available_slots // jobs_per

        log(f"Queue: {current_jobs} jobs active/pending  →  {can_submit} system slot(s) free")

        submitted_this_cycle = 0

        # ── Step 1: check already-submitted designs for failures & retry ───
        # Iterate all output dirs that have a setup subdir
        if output_root.exists():
            for setup_dir in sorted(output_root.rglob("setup")):
                if not setup_dir.is_dir():
                    continue
                pdb_dir    = setup_dir.parent
                design_name = pdb_dir.name

                state = load_state(pdb_dir)

                # Find the source PDB (need it for re-submission)
                rel    = pdb_dir.relative_to(output_root)
                subdir = rel.parent             # e.g. replica1
                pdb_candidate = data_dir / subdir / f"{design_name}.pdb"
                if not pdb_candidate.exists():
                    # PDB no longer in data/ — can't retry
                    continue

                # ── check setup ──────────────────────────────────────────
                setup_status = check_setup_complete(setup_dir)
                if setup_status == "failed":
                    retries = state.get("setup_retries", 0)
                    if retries < MAX_RETRIES and can_submit > 0:
                        # Record failed node before retry
                        node = read_node_from_run_info(setup_dir, "node")
                        if node:
                            state.setdefault("failed_nodes", []).append(node)
                        state["setup_retries"] = retries + 1
                        save_state(pdb_dir, state)
                        log(f"  RETRY setup ({retries + 1}/{MAX_RETRIES}): {design_name}"
                            f"  [failed node: {node or 'unknown'}]")
                        submit_system(pdb_candidate, pdb_dir, data_dir,
                                      mdp_dir, script_dir, args.gmx,
                                      args.replicas, args.dry_run,
                                      retry_setup=True)
                        can_submit    -= 1
                        submitted_this_cycle += 1
                    elif retries >= MAX_RETRIES:
                        log(f"  GAVE UP on setup for {design_name} "
                            f"(failed {retries} time(s))")
                    continue  # don't check replicas if setup isn't done

                if setup_status != "complete":
                    continue  # still running

                # ── check individual replicas ─────────────────────────────
                failed_reps: list[int] = []
                for i in range(1, args.replicas + 1):
                    rdir   = pdb_dir / f"replica_{i}"
                    status = check_replica_complete(rdir, i)
                    if status == "failed":
                        retries = state.get("replica_retries", {}).get(str(i), 0)
                        if retries < MAX_RETRIES:
                            failed_reps.append(i)
                        else:
                            log(f"  GAVE UP on {design_name} r{i} "
                                f"(failed {retries} time(s))")

                if failed_reps and can_submit > 0:
                    # Record nodes for failed replicas
                    for i in failed_reps:
                        node = read_node_from_run_info(setup_dir, f"r{i}_node")
                        if node:
                            state.setdefault("failed_nodes", []).append(node)
                        state.setdefault("replica_retries", {})[str(i)] = (
                            state.get("replica_retries", {}).get(str(i), 0) + 1
                        )
                    save_state(pdb_dir, state)
                    log(f"  RETRY replicas {failed_reps} for {design_name}")
                    submit_system(pdb_candidate, pdb_dir, data_dir,
                                  mdp_dir, script_dir, args.gmx,
                                  args.replicas, args.dry_run,
                                  retry_replicas=failed_reps)
                    # Retrying replicas only consumes 1 slot (no setup job)
                    can_submit    -= 1
                    submitted_this_cycle += 1

        # ── Step 2: submit new (not yet submitted) designs ─────────────────
        for pdb in find_pending_pdbs(data_dir):
            if can_submit <= 0:
                break
            if is_already_submitted(pdb, data_dir, output_root):
                continue

            pdb_dir = output_dir_for(pdb, data_dir, output_root)
            log(f"Submitting: {pdb.name}  →  {pdb_dir.relative_to(output_root)}")
            ok = submit_system(
                pdb, pdb_dir, data_dir, mdp_dir, script_dir,
                args.gmx, args.replicas, args.dry_run
            )
            if ok:
                can_submit           -= 1
                submitted_this_cycle += 1

        remaining_data = len([p for p in find_pending_pdbs(data_dir)
                               if not is_already_submitted(p, data_dir, output_root)])

        log(f"Submitted this cycle: {submitted_this_cycle} | "
            f"New designs remaining: {remaining_data}")

        # ── write Markdown status report every cycle ───────────────────────
        generate_status_report(output_root, base_dir, args.replicas)

        if args.once:
            log("--once flag set, exiting.")
            break

        if remaining_data == 0:
            # Check if all submitted designs are in a terminal state (complete or gave-up)
            all_terminal = True
            if output_root.exists():
                for setup_dir in output_root.rglob("setup"):
                    if not setup_dir.is_dir():
                        continue
                    pdb_dir = setup_dir.parent
                    setup_status = check_setup_complete(setup_dir)
                    if setup_status == "running":
                        all_terminal = False
                        break
                    if setup_status == "complete":
                        state = load_state(pdb_dir)
                        for i in range(1, args.replicas + 1):
                            rdir = pdb_dir / f"replica_{i}"
                            if check_replica_complete(rdir, i) == "running":
                                all_terminal = False
                                break
                    if not all_terminal:
                        break
            if all_terminal:
                log("All designs processed and queue is empty. Done.")
                break
            else:
                log(f"All PDBs submitted; waiting for running jobs to finish...")

        log(f"Sleeping {args.interval}s ...")
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
