#!/usr/bin/env python3
"""
Auto-submit manager for FEP ΔΔG calculations (HREX-FEP, CHARMM36m-mut).

Reads MD results from ~/MD/results/replica1/<design>/replica_{1..3}/ and processes
all 3 inner MD replicas per design as independent FEP calculation slots.
Never touches MD directory contents; copies inputs to FEP/inputs/replica{N}/<design>/.

Directory layout
----------------
MD source (read-only):
    ~/MD/results/replica1/<design>/       ← outer "replica1" groups all designs
        replica_1/                        ← inner replica WITH underscore (3 repeats)
            md_r1.gro                     ← final MD frame for replica 1
            md_r1.tpr                     ← matching TPR
        replica_2/
            md_r2.gro
            md_r2.tpr
        replica_3/
            md_r3.gro
            md_r3.tpr
        setup/
            system.pdb                    ← pdb2gmx output WITH hydrogens (shared)

FEP input structure (created by this script, mirrors MD source):
    ~/FEP/inputs/replica{N}/<design>/
        replica_{N}/
            md_r{N}.gro                   ← INPUT_STRUCTURE
            md_r{N}.tpr                   ← INPUT_TPR
        setup/
            system.pdb                    ← INPUT_PDB

FEP output structure (created by setup_fep.sh):
    ~/FEP/outputs/replica{N}/<design>/
        forward/
            bound_R{1..3}/lambda{0..31}/PROD/prod.xvg   ← HREX replicas
            unbound_R{1..3}/lambda{0..31}/PROD/prod.xvg
        reverse/  (set up manually after forward completes)

Usage
-----
    python scripts/auto_submit_fep.py --dry-run    # preview only
    python scripts/auto_submit_fep.py --once       # one-shot then exit
    python scripts/auto_submit_fep.py              # continuous (in screen)
    python scripts/auto_submit_fep.py --report     # status report only

    Live log:  tail -f ~/FEP/logs/auto_fep.log
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION  —  edit these paths before running on HPC
# ═══════════════════════════════════════════════════════════════════════════════

# Root directories (absolute paths on HPC)
FEP_ROOT    = Path("/public/home/xuziyi/FEP")
MD_ROOT     = Path("/public/home/xuziyi/MD")

# FEP subdirectories (derived from FEP_ROOT)
FEP_INPUTS_DIR  = FEP_ROOT / "inputs"
FEP_OUTPUTS_DIR = FEP_ROOT / "outputs"
FEP_SCRIPTS_DIR = FEP_ROOT / "scripts"
FEP_STATUS_FILE = FEP_ROOT / "fep_status.md"

# MD source: outer replica dirs (no underscore) contain all designs
# ~/MD/results/replica{1..N_MD_REPLICAS}/<design>/...
MD_RESULTS_DIR = MD_ROOT / "results"

# Number of outer MD replica directories to scan (replica1, replica2, replica3)
N_MD_REPLICAS = 3

# Inner replica used as FEP starting structure (always replica_1 with underscore)
# ~/MD/results/replica{N}/<design>/replica_1/md_r1.gro
INNER_REP_DIR = "replica_1"    # subdirectory under each design dir
INNER_GRO     = "md_r1.gro"   # GRO filename inside INNER_REP_DIR
INNER_TPR     = "md_r1.tpr"   # TPR filename inside INNER_REP_DIR
MD_SRC_PDB    = Path("setup/system.pdb")   # design-level PDB (no inner replica prefix)

# SLURM limits
QUICK_JOB_LIMIT  = 30  # max jobs in quick partition at once
JOBS_PER_SLOT    = 6   # jobs per (design, md_replica): bound_R{1..3} + unbound_R{1..3}

# Polling
DEFAULT_INTERVAL = 120  # seconds

# Retry
MAX_RETRIES = 1         # retry attempts per (design, md_replica) slot

# Sentinel file written by submit_fep.sh SLURM jobs on success
FEP_N_LAMBDA  = 32
LAST_WIN      = FEP_N_LAMBDA - 1   # 31

# Log file for real-time progress (appended across runs)
FEP_LOG_FILE  = FEP_ROOT / "logs" / "auto_fep.log"

# ═══════════════════════════════════════════════════════════════════════════════

# Module-level log file handle; set in main() before the first log() call
_log_fh = None


def log(msg: str) -> None:
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    if _log_fh is not None:
        _log_fh.write(line + "\n")
        _log_fh.flush()


# ── SLURM helpers ─────────────────────────────────────────────────────────────

def get_current_user() -> str:
    return os.environ.get("USER", subprocess.run(
        ["whoami"], capture_output=True, text=True).stdout.strip())


def count_quick_jobs(user: str) -> int:
    """Count jobs in the 'quick' partition for this user (any state)."""
    try:
        result = subprocess.run(
            ["squeue", "-u", user, "-p", "quick", "-h"],
            capture_output=True, text=True, timeout=15,
        )
        lines = [l for l in result.stdout.splitlines() if l.strip()]
        return len(lines)
    except Exception as e:
        log(f"WARNING: squeue failed ({e}); assuming 0 jobs")
        return 0


def is_job_in_queue(user: str, job_id: str) -> bool:
    """Return True if the given SLURM job ID is still in the queue."""
    try:
        result = subprocess.run(
            ["squeue", "-u", user, "-j", job_id, "-h"],
            capture_output=True, text=True, timeout=15,
        )
        return bool(result.stdout.strip())
    except Exception:
        return False


# ── MD input discovery ────────────────────────────────────────────────────────

def find_available_work() -> list[tuple[str, int]]:
    """Return sorted list of (design_name, md_replica) pairs where all required
    MD input files exist.

    Scans MD_RESULTS_DIR/replica{1..N}/<design>/INNER_REP_DIR/ for GRO+TPR,
    and MD_RESULTS_DIR/replica{N}/<design>/setup/system.pdb for the PDB."""
    work = []
    for rep_n in range(1, N_MD_REPLICAS + 1):
        outer_dir = MD_RESULTS_DIR / f"replica{rep_n}"
        if not outer_dir.exists():
            log(f"  NOTE: outer MD replica dir not found: {outer_dir}")
            continue
        for d in sorted(outer_dir.iterdir()):
            if not d.is_dir():
                continue
            gro = d / INNER_REP_DIR / INNER_GRO
            tpr = d / INNER_REP_DIR / INNER_TPR
            pdb = d / MD_SRC_PDB
            if gro.is_file() and tpr.is_file() and pdb.is_file():
                work.append((d.name, rep_n))
            else:
                missing = [name for name, p in [
                    (f"{INNER_REP_DIR}/{INNER_GRO}", gro),
                    (f"{INNER_REP_DIR}/{INNER_TPR}", tpr),
                    (str(MD_SRC_PDB), pdb),
                ] if not p.is_file()]
                log(f"  SKIP {d.name} (outer replica{rep_n}): missing {', '.join(missing)}")
    return work


# ── FEP input preparation ─────────────────────────────────────────────────────

def fep_replica_inputs_dir(design: str, md_replica: int) -> Path:
    """Input directory for a specific (design, md_replica) pair."""
    return FEP_INPUTS_DIR / f"replica{md_replica}" / design


def prepare_fep_inputs(design: str, md_replica: int, dry_run: bool) -> bool:
    """Copy GRO/TPR/PDB from MD results into FEP/inputs/replica{N}/<design>/.

    Source  (MD, outer replica N, inner always replica_1):
        MD/results/replica{N}/<design>/replica_1/md_r1.gro  ← INPUT_STRUCTURE
        MD/results/replica{N}/<design>/replica_1/md_r1.tpr  ← INPUT_TPR
        MD/results/replica{N}/<design>/setup/system.pdb     ← INPUT_PDB

    Destination (FEP inputs, same relative layout):
        FEP/inputs/replica{N}/<design>/replica_1/md_r1.gro
        FEP/inputs/replica{N}/<design>/replica_1/md_r1.tpr
        FEP/inputs/replica{N}/<design>/setup/system.pdb

    Idempotent: already-copied files are not overwritten.
    Returns True on success."""
    src_dir = MD_RESULTS_DIR / f"replica{md_replica}" / design
    dst_dir = fep_replica_inputs_dir(design, md_replica)

    files = {
        src_dir / INNER_REP_DIR / INNER_GRO: dst_dir / INNER_REP_DIR / INNER_GRO,
        src_dir / INNER_REP_DIR / INNER_TPR: dst_dir / INNER_REP_DIR / INNER_TPR,
        src_dir / MD_SRC_PDB:                dst_dir / MD_SRC_PDB,
    }

    if dry_run:
        log(f"  [DRY-RUN] Would copy inputs for {design}/replica{md_replica}")
        return True

    try:
        for src, dst in files.items():
            if dst.exists():
                continue   # already copied — do not overwrite
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        return True
    except Exception as e:
        log(f"  ERROR: copy inputs for {design}/replica{md_replica}: {e}")
        return False


# ── FEP status helpers ────────────────────────────────────────────────────────

def fep_replica_dir(design: str, md_replica: int) -> Path:
    """Output directory for a specific (design, md_replica) pair."""
    return FEP_OUTPUTS_DIR / f"replica{md_replica}" / design


def is_setup_done(design: str, md_replica: int) -> bool:
    """Return True when forward setup is present (any *ions.gro exists in bound_R1)."""
    marker = fep_replica_dir(design, md_replica) / "forward" / "bound_R1"
    return any(marker.glob("*ions.gro")) if marker.is_dir() else False


def hrex_rep_complete(design: str, md_replica: int,
                      direction: str, state: str, hrex_rep: int) -> bool:
    """Return True if lambda31/PROD/prod.gro exists for an HREX replica."""
    f = (fep_replica_dir(design, md_replica) / direction
         / f"{state}_R{hrex_rep}" / f"lambda{LAST_WIN}" / "PROD" / "prod.gro")
    return f.is_file() and f.stat().st_size > 0


def count_complete_hrex(design: str, md_replica: int, direction: str) -> dict:
    """Return {state: n_complete} for bound and unbound HREX replicas."""
    return {
        state: sum(1 for r in range(1, 4)
                   if hrex_rep_complete(design, md_replica, direction, state, r))
        for state in ("bound", "unbound")
    }


def slot_forward_complete(design: str, md_replica: int) -> bool:
    """Return True when all 6 forward PROD jobs are done for this (design, replica) slot."""
    c = count_complete_hrex(design, md_replica, "forward")
    return c["bound"] == 3 and c["unbound"] == 3


# ── Per-design state file ─────────────────────────────────────────────────────

def _state_path(design: str, md_replica: int) -> Path:
    return fep_replica_dir(design, md_replica) / ".fep_state.json"


def load_state(design: str, md_replica: int) -> dict:
    p = _state_path(design, md_replica)
    if p.exists():
        try:
            return json.loads(p.read_text())
        except Exception:
            pass
    return {"retries": 0, "slurm_job_ids": []}


def save_state(design: str, md_replica: int, state: dict) -> None:
    p = _state_path(design, md_replica)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(state, indent=2))


# ── Setup + submission ────────────────────────────────────────────────────────

def run_setup_fep(design: str, md_replica: int, dry_run: bool) -> bool:
    """Run setup_fep.sh --direction forward --md-replica N for this (design, replica).
    Patches DESIGN_NAME and MD_REPLICA in a temp copy so the master script is unchanged.
    Returns True on success."""
    import re
    setup_script = FEP_SCRIPTS_DIR / "setup_fep.sh"
    if not setup_script.exists():
        log(f"  ERROR: setup_fep.sh not found at {setup_script}")
        return False

    if dry_run:
        log(f"  [DRY-RUN] Would run setup_fep.sh for {design}/replica{md_replica}")
        return True

    tmp_dir = fep_replica_inputs_dir(design, md_replica)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp_script = tmp_dir / "setup_fep_run.sh"

    script_text = setup_script.read_text()
    script_text = re.sub(r'^(DESIGN_NAME=).*$',
                         f'DESIGN_NAME="{design}"', script_text, flags=re.MULTILINE)
    script_text = re.sub(r'^(MD_REPLICA=).*$',
                         f'MD_REPLICA={md_replica}', script_text, flags=re.MULTILINE)
    tmp_script.write_text(script_text)
    tmp_script.chmod(0o755)

    log(f"  Running setup_fep.sh for {design}/replica{md_replica}...")
    setup_log = fep_replica_dir(design, md_replica) / "setup_fep.log"
    setup_log.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(setup_log, "w") as flog:
            result = subprocess.run(
                ["bash", str(tmp_script), "--direction", "forward"],
                stdout=flog, stderr=subprocess.STDOUT, timeout=900,
            )
        if result.returncode != 0:
            log(f"  ERROR: setup_fep.sh failed for {design}/replica{md_replica}. See {setup_log}")
            return False
        log(f"  Setup done for {design}/replica{md_replica}")
        return True
    except subprocess.TimeoutExpired:
        log(f"  ERROR: setup_fep.sh timed out for {design}/replica{md_replica}")
        return False
    except Exception as e:
        log(f"  ERROR: setup_fep.sh exception for {design}/replica{md_replica}: {e}")
        return False


def run_submit_fep(design: str, md_replica: int, direction: str,
                   dry_run: bool, force: bool = False) -> list[str]:
    """Run submit_fep.sh for this (design, md_replica) and direction.
    Returns list of submitted SLURM job IDs (empty on failure or dry-run)."""
    import re
    submit_script = FEP_SCRIPTS_DIR / "submit_fep.sh"
    if not submit_script.exists():
        log(f"  ERROR: submit_fep.sh not found at {submit_script}")
        return []

    if dry_run:
        log(f"  [DRY-RUN] Would submit {direction} for {design}/replica{md_replica}")
        return []

    tmp_dir = fep_replica_inputs_dir(design, md_replica)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp_script = tmp_dir / "submit_fep_run.sh"

    script_text = submit_script.read_text()
    script_text = re.sub(r'^(DESIGN_NAME=).*$',
                         f'DESIGN_NAME="{design}"', script_text, flags=re.MULTILINE)
    script_text = re.sub(r'^(MD_REPLICA=).*$',
                         f'MD_REPLICA={md_replica}', script_text, flags=re.MULTILINE)
    tmp_script.write_text(script_text)
    tmp_script.chmod(0o755)

    cmd = ["bash", str(tmp_script), direction]
    if force:
        cmd.append("--force")

    log(f"  Submitting {direction} for {design}/replica{md_replica}...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        job_ids = []
        for line in result.stdout.splitlines():
            print(line)
            if line.startswith("SUBMIT:"):
                parts = line.split()
                if parts and parts[-1].isdigit():
                    job_ids.append(parts[-1])

        if result.returncode != 0:
            log(f"  WARNING: submit_fep.sh returned {result.returncode} "
                f"for {design}/replica{md_replica}")
        return job_ids
    except subprocess.TimeoutExpired:
        log(f"  ERROR: submit_fep.sh timed out for {design}/replica{md_replica}")
        return []
    except Exception as e:
        log(f"  ERROR: submit_fep.sh exception for {design}/replica{md_replica}: {e}")
        return []


# ── Status report ─────────────────────────────────────────────────────────────

def generate_status_report(all_work: list[tuple[str, int]]) -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Group by design for display
    design_replicas: dict[str, list[int]] = {}
    for design, rep_n in all_work:
        design_replicas.setdefault(design, []).append(rep_n)

    lines = [
        f"# FEP Status  —  {now}",
        "",
        f"MD source : `{MD_RESULTS_DIR}/replica{{1..{N_MD_REPLICAS}}}/<design>/{INNER_REP_DIR}/`",
        f"FEP root  : `{FEP_ROOT}`",
        "",
        "| Design | MD rep | Inputs | Setup | fwd bound (h1/h2/h3) | fwd unbound (h1/h2/h3) | Done | Notes |",
        "|--------|--------|--------|-------|----------------------|------------------------|------|-------|",
    ]

    n_total = n_done = n_fail = 0

    for design in sorted(design_replicas):
        for rep_n in sorted(design_replicas[design]):
            n_total += 1
            inp_dir   = fep_replica_inputs_dir(design, rep_n)
            inputs_ok = (inp_dir / INNER_REP_DIR / INNER_GRO).is_file()
            setup_ok  = is_setup_done(design, rep_n)
            st        = load_state(design, rep_n)

            def hrex_sym(state: str, r: int, _rep_n: int = rep_n) -> str:
                if hrex_rep_complete(design, _rep_n, "forward", state, r):
                    return "✅"
                d = fep_replica_dir(design, _rep_n) / "forward" / f"{state}_R{r}"
                return "🔄" if d.exists() else "—"

            b_syms = " / ".join(hrex_sym("bound",   r) for r in range(1, 4))
            u_syms = " / ".join(hrex_sym("unbound", r) for r in range(1, 4))
            done   = "✅" if slot_forward_complete(design, rep_n) else "—"

            retries = st.get("retries", 0)
            notes = []
            if retries > 0:
                notes.append(f"retry {retries}×")
            if retries >= MAX_RETRIES and not slot_forward_complete(design, rep_n) and setup_ok:
                notes.append("gave up")
                n_fail += 1
            if slot_forward_complete(design, rep_n):
                n_done += 1

            inp_sym   = "✅" if inputs_ok else "❌"
            setup_sym = "✅" if setup_ok else (
                "🔄" if (fep_replica_dir(design, rep_n) / "forward").exists() else "—")
            lines.append(
                f"| {design} | r{rep_n} | {inp_sym} | {setup_sym} "
                f"| {b_syms} | {u_syms} | {done} | {'; '.join(notes)} |"
            )

    lines += [
        "",
        "---",
        f"**Total slots: {n_total}** | Complete: {n_done} | Failed/gave-up: {n_fail}",
        "",
        "_Key: ✅ done | 🔄 submitted/running | — not started | ❌ missing_",
        "_Columns: h1/h2/h3 = HREX replicas within each FEP calculation_",
    ]

    FEP_STATUS_FILE.write_text("\n".join(lines) + "\n")
    log(f"Status report → {FEP_STATUS_FILE}")


# ── Main loop ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="FEP auto-submit manager (HREX-FEP, all MD replicas)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--dry-run",  action="store_true",
                        help="Print what would happen without submitting")
    parser.add_argument("--once",     action="store_true",
                        help="Submit one batch then exit (no polling)")
    parser.add_argument("--report",   action="store_true",
                        help="Generate status report and exit")
    parser.add_argument("--interval", type=int, default=DEFAULT_INTERVAL,
                        help=f"Poll interval in seconds (default: {DEFAULT_INTERVAL})")
    args = parser.parse_args()

    global _log_fh
    FEP_LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    _log_fh = open(FEP_LOG_FILE, "a", encoding="utf-8")

    try:
        _run(args)
    finally:
        _log_fh.close()
        _log_fh = None


def _run(args) -> None:
    user = get_current_user()

    if args.report:
        all_work = find_available_work()
        generate_status_report(all_work)
        return

    max_slots = QUICK_JOB_LIMIT // JOBS_PER_SLOT
    log("=" * 64)
    log("FEP Auto-submit Manager")
    log(f"  MD source    : {MD_RESULTS_DIR}/replica{{1..{N_MD_REPLICAS}}}/<design>/{INNER_REP_DIR}/")
    log(f"  FEP root     : {FEP_ROOT}")
    log(f"  Job limit    : {QUICK_JOB_LIMIT} (quick partition)")
    log(f"  Jobs/slot    : {JOBS_PER_SLOT}  → max {max_slots} (design, outer_replica) slots at once")
    log(f"  Outer reps   : {N_MD_REPLICAS} (replica1..replica{N_MD_REPLICAS}, no underscore)")
    log(f"  Inner GRO    : {INNER_REP_DIR}/{INNER_GRO}  (always, WITH underscore)")
    log(f"  Log file     : {FEP_LOG_FILE}")
    log(f"  Dry-run      : {args.dry_run}")
    log("=" * 64)

    while True:
        all_work = find_available_work()
        if not all_work:
            log("No (design, md_replica) pairs found in MD results. Exiting.")
            break

        current_jobs = count_quick_jobs(user)
        available    = max(0, QUICK_JOB_LIMIT - current_jobs)
        can_submit   = available // JOBS_PER_SLOT

        log(f"Quick-partition jobs: {current_jobs}/{QUICK_JOB_LIMIT}  "
            f"→ {can_submit} slot(s) free")

        submitted_this_cycle = 0

        # ── Phase 1: retry failed (design, md_replica) slots ────────────────────
        for design, md_replica in all_work:
            if can_submit <= 0:
                break
            if slot_forward_complete(design, md_replica):
                continue

            st      = load_state(design, md_replica)
            retries = st.get("retries", 0)
            if retries == 0 or retries >= MAX_RETRIES:
                continue   # Phase 2 handles first-time; gave up otherwise

            if is_setup_done(design, md_replica):
                c = count_complete_hrex(design, md_replica, "forward")
                if c["bound"] + c["unbound"] < 6:
                    log(f"  RETRY ({retries + 1}/{MAX_RETRIES}): {design}/replica{md_replica}")
                    job_ids = run_submit_fep(design, md_replica, "forward",
                                            args.dry_run, force=True)
                    st["retries"] = retries + 1
                    st.setdefault("slurm_job_ids", []).extend(job_ids)
                    save_state(design, md_replica, st)
                    if job_ids or args.dry_run:
                        can_submit           -= 1
                        submitted_this_cycle += 1

        # ── Phase 2: new (design, md_replica) slots ───────────────────────────
        for design, md_replica in all_work:
            if can_submit <= 0:
                break
            if slot_forward_complete(design, md_replica):
                continue

            st = load_state(design, md_replica)
            if st.get("retries", 0) > 0:
                continue   # handled in Phase 1

            if is_setup_done(design, md_replica):
                log(f"  Setup done, submitting: {design}/replica{md_replica}")
                job_ids = run_submit_fep(design, md_replica, "forward", args.dry_run)
            else:
                log(f"  New slot: {design}/replica{md_replica}")
                ok = prepare_fep_inputs(design, md_replica, args.dry_run)
                if not ok:
                    continue
                ok = run_setup_fep(design, md_replica, args.dry_run)
                if not ok:
                    log(f"  Setup failed for {design}/replica{md_replica}, skipping")
                    continue
                job_ids = run_submit_fep(design, md_replica, "forward", args.dry_run)

            if job_ids or args.dry_run:
                st["retries"] = 0
                st.setdefault("slurm_job_ids", []).extend(job_ids or [])
                save_state(design, md_replica, st)
                can_submit           -= 1
                submitted_this_cycle += 1

        log(f"Submitted this cycle: {submitted_this_cycle}")

        generate_status_report(all_work)

        if args.once:
            log("--once flag set, exiting.")
            break

        all_done = all(
            slot_forward_complete(d, r) or load_state(d, r).get("retries", 0) >= MAX_RETRIES
            for d, r in all_work
        )
        if all_done:
            log("All (design, replica) slots complete or gave up. Done.")
            break

        remaining = sum(1 for d, r in all_work
                        if not slot_forward_complete(d, r)
                        and load_state(d, r).get("retries", 0) < MAX_RETRIES)
        log(f"Slots not yet complete: {remaining}. Sleeping {args.interval}s ...")
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
