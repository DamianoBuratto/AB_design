#!/usr/bin/env python3
"""
batch_analysis.py — Batch MD RMSD analysis for all designs.

Scans  BASE_DIR/{replica1,replica2,replica3}/design_NAME/
for every design folder, runs antibody + CDR backbone RMSD analysis for each
simulation replica (1-3), and outputs:

  Per-design (matching Analysis.ipynb format exactly):
    analysis_results/RMSD_antibody_replica_N.csv
    analysis_results/RMSD_antibody_statistics_summary.csv
    analysis_results/RMSD_antibody_all_replicas.png   ← same as notebook §7 plot
    analysis_results/CDR_RMSD_replica_N.csv
    analysis_results/CDR_RMSD_statistics_summary.csv
    analysis_results/CDR_RMSD_loops.png               ← same as notebook §8 plot

  Master summary CSV:
    BASE_DIR/batch_summary/master_rmsd_summary.csv    ← one row per design×replica

Usage:
    python3 batch_analysis.py
    python3 batch_analysis.py --base-dir ~/MD/results
    python3 batch_analysis.py --base-dir ~/MD/results --force
    python3 batch_analysis.py --base-dir ~/MD/results --stats-skip-ns 50
    nohup python3 batch_analysis.py --base-dir ~/MD/results > batch_analysis.log 2>&1 &
"""

import argparse
import re
import sys
import traceback
from collections import defaultdict
from pathlib import Path

# Non-interactive Agg backend MUST be set before importing pyplot.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis import transformations

# ── Default configuration ─────────────────────────────────────────────────────
_DEFAULT_BASE_DIR      = Path.home() / 'MD' / 'results'
_DESIGN_GROUPS         = ['replica1', 'replica2', 'replica3']
_ACTIVE_REPLICAS       = [1, 2, 3]
_DEFAULT_STATS_SKIP_NS = 100
_REF_FRAME             = 0
_RMSD_STEP             = 1
_PHLA_CHAINS           = ['A', 'P']
_ANTIBODY_CHAINS       = ['H', 'L']
_CDR_LABEL_TO_NAME     = {'H1': 'H-CDR1', 'H2': 'H-CDR2', 'H3': 'H-CDR3',
                          'L1': 'L-CDR1', 'L2': 'L-CDR2', 'L3': 'L-CDR3'}
_CDR_LABEL_TO_CHAIN    = {'H1': 'H', 'H2': 'H', 'H3': 'H',
                          'L1': 'L', 'L2': 'L', 'L3': 'L'}
_CDR_ORDER             = ['H-CDR1', 'H-CDR2', 'H-CDR3', 'L-CDR1', 'L-CDR2', 'L-CDR3']
_REPLICA_COLORS        = {1: '#2E86AB', 2: '#F18F01', 3: '#66BB6A'}

# ── Plot style (identical to Analysis.ipynb §1) ───────────────────────────────
def _setup_plot_style():
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'figure.dpi':         300,
        'font.size':          10,
        'axes.titlesize':     11,
        'axes.labelsize':     10,
        'xtick.labelsize':    8,
        'ytick.labelsize':    8,
        'legend.fontsize':    9,
        'lines.linewidth':    1.0,
        'axes.linewidth':     0.8,
        'grid.alpha':         0.3,
        'savefig.bbox':       'tight',
        'savefig.dpi':        300,
        'font.family':        'sans-serif',
        'font.sans-serif':    ['Arial', 'DejaVu Sans', 'Helvetica'],
    })


# ═══════════════════════════════════════════════════════════════════════════════
# §4 — Topology / trajectory helpers  (copied verbatim from Analysis.ipynb)
# ═══════════════════════════════════════════════════════════════════════════════

def find_topology(replica_dir, replica_num):
    """Return (path, type_str) for the best available topology file."""
    for name in [f'md_r{replica_num}.tpr', f'production_r{replica_num}.tpr', 'topol.tpr']:
        p = replica_dir / name
        if p.exists():
            return p, 'tpr'
    processed = replica_dir.parent / 'setup' / 'processed.pdb'
    if processed.exists():
        return processed, 'pdb'
    for name in ['npt.gro', f'md_r{replica_num}.gro']:
        p = replica_dir / name
        if p.exists():
            return p, 'gro_npt' if name == 'npt.gro' else 'gro_md'
    return None, None


def find_trajectory(replica_dir, replica_num):
    """Return path to XTC trajectory, or None."""
    for name in [f'md_r{replica_num}.xtc', f'production_r{replica_num}.xtc']:
        p = replica_dir / name
        if p.exists():
            return p
    return None


def build_chain_selection(universe, chain_letters, atom_sel='backbone'):
    """Return MDAnalysis selection string for chains (chainID or segid fallback)."""
    protein = universe.select_atoms('protein')
    try:
        chain_ids = sorted(set(c for c in protein.chainIDs if c.strip()))
        if chain_ids:
            return atom_sel + ' and (' + ' or '.join(f'chainID {c}' for c in chain_letters) + ')'
    except (AttributeError, mda.NoDataError):
        pass
    segids = sorted(set(protein.segids))
    if any('Protein_chain_' in s for s in segids):
        parts = []
        for c in chain_letters:
            match = next((s for s in segids if f'Protein_chain_{c}' in s), None)
            if match:
                parts.append(f'segid {match}')
        if parts:
            return atom_sel + ' and (' + ' or '.join(parts) + ')'
    raise ValueError(f'No chain info found (segids: {segids[:6]}). Use TPR or processed.pdb.')


def detect_cdr_from_pdb(pdb_path):
    """Parse CDR residue ranges from REMARK PDBinfo-LABEL lines in processed.pdb."""
    text = Path(pdb_path).read_text()
    cdr_pdb = {}
    for m in re.finditer(r'REMARK PDBinfo-LABEL:\s+(\d+)\s+(H[123]|L[123])', text):
        resid, label = int(m.group(1)), m.group(2)
        cdr_pdb.setdefault(label, []).append(resid)
    return {k: sorted(v) for k, v in cdr_pdb.items()}


def remap_cdr_to_universe(universe, pdb_path):
    """Convert PDB-local CDR resids to universe (TPR) resids using rank within chain."""
    cdr_pdb = detect_cdr_from_pdb(pdb_path)
    if not cdr_pdb:
        return {}

    text = Path(pdb_path).read_text()
    pdb_chains = defaultdict(set)
    for m in re.finditer(r'^ATOM.{13}.{4}(\w)\s+(\d+)', text, re.MULTILINE):
        pdb_chains[m.group(1)].add(int(m.group(2)))
    pdb_chains = {ch: sorted(rids) for ch, rids in pdb_chains.items()}

    protein = universe.select_atoms('protein')
    try:
        u_chains = {}
        for ch in set(protein.chainIDs):
            if ch.strip():
                u_chains[ch] = sorted(set(protein.select_atoms(f'chainID {ch}').resids))
    except (AttributeError, mda.NoDataError):
        return {}

    result = {}
    for label, pdb_rids in cdr_pdb.items():
        chain = _CDR_LABEL_TO_CHAIN.get(label)
        if not chain or chain not in pdb_chains or chain not in u_chains:
            print(f'    ⚠ Could not map {label}: chain {chain} missing from PDB or universe',
                  flush=True)
            continue
        try:
            rank_start = pdb_chains[chain].index(pdb_rids[0])
            rank_end   = pdb_chains[chain].index(pdb_rids[-1])
        except ValueError:
            print(f'    ⚠ {label}: resid {pdb_rids[0]} or {pdb_rids[-1]} not in PDB chain {chain}',
                  flush=True)
            continue
        u_rids = u_chains[chain]
        if rank_end >= len(u_rids):
            print(f'    ⚠ {label}: rank {rank_end} out of range for chain {chain}', flush=True)
            continue
        tpr_start, tpr_end = u_rids[rank_start], u_rids[rank_end]
        result[_CDR_LABEL_TO_NAME[label]] = (chain, (tpr_start, tpr_end))

    return result


# ═══════════════════════════════════════════════════════════════════════════════
# §5 — PBC removal  (copied verbatim from Analysis.ipynb)
# ═══════════════════════════════════════════════════════════════════════════════

def remove_pbc(universe, protein_selection='protein'):
    """Apply NoJump → Unwrap → Center PBC corrections."""
    universe.trajectory._transformations = []
    protein = universe.select_atoms(protein_selection)

    has_bonds = False
    try:
        n = len(universe.atoms.bonds)
        if n > 0:
            has_bonds = True
            print(f'    {n:,} bonds — Unwrap enabled', flush=True)
    except (AttributeError, mda.NoDataError):
        pass

    if not has_bonds:
        try:
            protein.guess_bonds()
            has_bonds = True
            print('    Bonds guessed — Unwrap enabled', flush=True)
        except Exception as e:
            print(f'    ⚠ Bond guessing failed ({e}) — NoJump + Center only', flush=True)

    if has_bonds:
        workflow = [
            transformations.nojump.NoJump(),
            transformations.unwrap(protein),
            transformations.center_in_box(protein, center='geometry'),
        ]
        print('    ✓ NoJump + Unwrap + Center', flush=True)
    else:
        workflow = [
            transformations.nojump.NoJump(),
            transformations.center_in_box(protein, center='geometry'),
        ]
        print('    ✓ NoJump + Center (Unwrap skipped)', flush=True)

    universe.trajectory.add_transformations(*workflow)
    return universe


# ═══════════════════════════════════════════════════════════════════════════════
# §7 — Antibody RMSD  (copied verbatim from Analysis.ipynb)
# ═══════════════════════════════════════════════════════════════════════════════

def compute_rmsd_antibody(universe, fit_selection, group_selection,
                          reference_frame=0, skip_frames=0, step=1, stats_skip_ns=100):
    """Fit on pHLA, compute antibody backbone RMSD (no re-fitting)."""
    R = rms.RMSD(
        universe, universe,
        select=fit_selection,
        groupselections=[group_selection],
        ref_frame=reference_frame,
        verbose=False,
    )
    R.run(start=skip_frames, step=step, verbose=False)

    arr     = R.results.rmsd
    time_ns = arr[:, 1] / 1000
    rmsd_ab = arr[:, 3]

    mask  = time_ns > stats_skip_ns
    stats = rmsd_ab[mask] if np.any(mask) else rmsd_ab
    n_stats = int(np.sum(mask)) if np.any(mask) else len(rmsd_ab)

    if not np.any(mask):
        print(f'    ⚠ trajectory < {stats_skip_ns} ns — using all frames for stats', flush=True)

    result = {
        'frames':    arr[:, 0].astype(int),
        'time_ns':   time_ns,
        'rmsd_fit':  arr[:, 2],
        'rmsd':      rmsd_ab,
        'mean':      float(np.mean(stats)),
        'std':       float(np.std(stats, ddof=1)),
        'max':       float(np.max(stats)),
        'min':       float(np.min(stats)),
        'n_stats':   n_stats,
        'duration_ns': float(time_ns[-1]) if len(time_ns) else 0.0,
    }
    print(f'    ✓ Ab RMSD (t>{stats_skip_ns}ns): {result["mean"]:.2f} ± {result["std"]:.2f} Å'
          f'  max={result["max"]:.2f} Å', flush=True)
    if result['max'] > 15.0:
        print(f'    ⚠ High RMSD — check for PBC artifacts or antibody dissociation.', flush=True)
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# §8 — CDR RMSD  (copied verbatim from Analysis.ipynb)
# ═══════════════════════════════════════════════════════════════════════════════

def build_cdr_sel(universe, chain, resid_start, resid_end, atom_sel='backbone'):
    """Build MDAnalysis selection for a CDR region (chainID or segid fallback)."""
    protein = universe.select_atoms('protein')
    try:
        chain_ids = sorted(set(c for c in protein.chainIDs if c.strip()))
        if chain_ids:
            return f'{atom_sel} and chainID {chain} and resid {resid_start}:{resid_end}'
    except (AttributeError, mda.NoDataError):
        pass
    segids = sorted(set(protein.segids))
    match = next((s for s in segids if f'Protein_chain_{chain}' in s), None)
    if match:
        return f'{atom_sel} and segid {match} and resid {resid_start}:{resid_end}'
    raise ValueError(f'Chain {chain} not found in universe')


def compute_cdr_rmsd(universe, fit_selection, cdr_residues,
                     reference_frame=0, skip_frames=0, step=1, stats_skip_ns=100):
    """Fit on pHLA, compute per-CDR backbone RMSD in a single trajectory pass."""
    cdr_names = list(cdr_residues.keys())
    cdr_sels  = []
    for name, (chain, (s, e)) in cdr_residues.items():
        sel = build_cdr_sel(universe, chain, s, e)
        n   = len(universe.select_atoms(sel))
        flag = '✓' if n > 0 else '⚠ EMPTY'
        print(f'    {name}: chain {chain}  resid {s}-{e}  {n} backbone atoms  {flag}', flush=True)
        cdr_sels.append(sel)

    empty = [n for n, (ch, (s, e)) in cdr_residues.items()
             if len(universe.select_atoms(build_cdr_sel(universe, ch, s, e))) == 0]
    if empty:
        raise ValueError(f'Empty CDR selections: {empty}. Check CDR_RESIDUES resids.')

    R = rms.RMSD(
        universe, universe,
        select=fit_selection,
        groupselections=cdr_sels,
        ref_frame=reference_frame,
        verbose=False,
    )
    R.run(start=skip_frames, step=step, verbose=False)

    arr     = R.results.rmsd
    time_ns = arr[:, 1] / 1000
    mask    = time_ns > stats_skip_ns
    out     = {'time_ns': time_ns}

    for i, name in enumerate(cdr_names):
        rmsd  = arr[:, 3 + i]
        stats = rmsd[mask] if np.any(mask) else rmsd
        out[name] = {
            'rmsd': rmsd,
            'mean': float(np.mean(stats)),
            'std':  float(np.std(stats, ddof=1)),
            'max':  float(np.max(stats)),
            'min':  float(np.min(stats)),
        }
        print(f'    {name}: {out[name]["mean"]:.2f} ± {out[name]["std"]:.2f} Å'
              f'  max={out[name]["max"]:.2f}', flush=True)
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Plot functions — produce figures identical to notebook §7 and §8
# ═══════════════════════════════════════════════════════════════════════════════

def plot_rmsd_antibody(rmsd_results, system_name, output_dir):
    """Antibody RMSD time series — all replicas on one figure (= notebook §7 plot)."""
    fig, ax = plt.subplots(figsize=(10, 5))
    for key, data in rmsd_results.items():
        if data is None:
            continue
        rep_num = int(key.split('_')[-1])
        color   = _REPLICA_COLORS.get(rep_num)
        line    = ax.plot(data['time_ns'], data['rmsd'], alpha=0.7, color=color,
                          label=f"{key}  ({data['mean']:.2f}±{data['std']:.2f} Å)")
        ax.axhline(data['mean'], linestyle='--', alpha=0.4,
                   color=color or line[0].get_color())

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title(f'Antibody backbone RMSD  [fit: pHLA | ref: frame 0]\n{system_name}')
    ax.legend(fontsize='small')
    plt.tight_layout()
    out_path = output_dir / 'RMSD_antibody_all_replicas.png'
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out_path


def plot_cdr_rmsd(cdr_rmsd_results, cdr_names, system_name, output_dir):
    """CDR RMSD 2×3 grid — all replicas per subplot (= notebook §8 CDR plot)."""
    n_cdrs = len(cdr_names)
    ncols  = 3
    nrows  = (n_cdrs + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4 * nrows), sharex=True)
    axes = np.array(axes).flatten()

    for ax_idx, cdr in enumerate(cdr_names):
        ax = axes[ax_idx]
        for rep_key, data in cdr_rmsd_results.items():
            if data is None or cdr not in data:
                continue
            rep_num = int(rep_key.split('_')[-1])
            color   = _REPLICA_COLORS.get(rep_num)
            d       = data[cdr]
            line    = ax.plot(data['time_ns'], d['rmsd'], alpha=0.7, color=color,
                              label=f"R{rep_num} ({d['mean']:.1f}±{d['std']:.1f} Å)")
            ax.axhline(d['mean'], linestyle='--', alpha=0.4,
                       color=color or line[0].get_color())
        ax.set_title(cdr, fontweight='bold')
        ax.set_ylabel('RMSD (Å)')
        ax.legend(fontsize=7)

    for ax in axes[n_cdrs:]:
        ax.set_visible(False)
    last_row_start = n_cdrs - ncols if n_cdrs > ncols else 0
    for ax in axes[last_row_start:n_cdrs]:
        ax.set_xlabel('Time (ns)')

    fig.suptitle(f'CDR Loop Backbone RMSD  [fit: pHLA | ref: frame 0]\n{system_name}',
                 fontsize=11)
    plt.tight_layout()
    out_path = output_dir / 'CDR_RMSD_loops.png'
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out_path


# ═══════════════════════════════════════════════════════════════════════════════
# Save functions — identical column names to notebook §7 and §8 save cells
# ═══════════════════════════════════════════════════════════════════════════════

def save_rmsd_antibody(rmsd_results, output_dir, stats_skip_ns):
    """Save per-replica RMSD CSVs and summary — identical to notebook §7 save cell."""
    for replica_key, data in rmsd_results.items():
        if data is None:
            continue
        df = pd.DataFrame({
            'Frame':           data['frames'],
            'Time_ns':         data['time_ns'],
            'RMSD_pHLA_fit_A': data['rmsd_fit'],
            'RMSD_antibody_A': data['rmsd'],
        })
        out = output_dir / f'RMSD_antibody_{replica_key}.csv'
        df.to_csv(out, index=False, float_format='%.3f')

    rows = [{
        'Replica':        int(k.split('_')[-1]),
        'Mean_Ab_RMSD_A': v['mean'],
        'Std_Ab_RMSD_A':  v['std'],
        'Max_Ab_RMSD_A':  v['max'],
        'Min_Ab_RMSD_A':  v['min'],
    } for k, v in rmsd_results.items() if v is not None]
    pd.DataFrame(rows).to_csv(
        output_dir / 'RMSD_antibody_statistics_summary.csv', index=False, float_format='%.3f'
    )


def save_cdr_rmsd(cdr_rmsd_results, cdr_names, output_dir):
    """Save per-replica CDR CSVs and summary — identical to notebook §8 save cell."""
    for rep_key, data in cdr_rmsd_results.items():
        if data is None:
            continue
        df = pd.DataFrame({
            'Time_ns': data['time_ns'],
            **{cdr: data[cdr]['rmsd'] for cdr in cdr_names if cdr in data},
        })
        df.to_csv(output_dir / f'CDR_RMSD_{rep_key}.csv', index=False, float_format='%.3f')

    rows = []
    for rep_key, data in cdr_rmsd_results.items():
        if data is None:
            continue
        row = {'Replica': int(rep_key.split('_')[-1])}
        for cdr in cdr_names:
            if cdr in data:
                row[f'{cdr}_mean_A'] = data[cdr]['mean']
                row[f'{cdr}_std_A']  = data[cdr]['std']
        rows.append(row)
    pd.DataFrame(rows).to_csv(
        output_dir / 'CDR_RMSD_statistics_summary.csv', index=False, float_format='%.3f'
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Per-design processing
# ═══════════════════════════════════════════════════════════════════════════════

def _make_failed_row(design_group, system_name, sim_replica, status):
    row = {
        'design_group': design_group, 'system_name': system_name,
        'sim_replica': sim_replica, 'status': status,
    }
    for col in ['duration_ns', 'n_frames_stats', 'Ab_RMSD_mean_A', 'Ab_RMSD_std_A',
                'Ab_RMSD_max_A', 'Ab_RMSD_min_A']:
        row[col] = float('nan')
    for cdr in _CDR_ORDER:
        col = cdr.replace('-', '_')
        row[f'{col}_mean_A'] = float('nan')
        row[f'{col}_std_A']  = float('nan')
    return row


def process_design(design_group, system_name, base_dir,
                   stats_skip_ns=_DEFAULT_STATS_SKIP_NS,
                   force=False):
    """
    Process all simulation replicas for one design.

    Returns
    -------
    list[dict]  — one master-CSV row per simulation replica
    """
    system_dir = base_dir / design_group / system_name
    output_dir = system_dir / 'analysis_results'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Skip if already done (unless --force)
    sentinel = output_dir / 'RMSD_antibody_statistics_summary.csv'
    if sentinel.exists() and not force:
        print(f'  [SKIP] {design_group}/{system_name} — already analysed (use --force to redo)',
              flush=True)
        # Read existing stats to fill master CSV
        rows = []
        try:
            df = pd.read_csv(sentinel)
            cdr_summary_path = output_dir / 'CDR_RMSD_statistics_summary.csv'
            cdr_df = pd.read_csv(cdr_summary_path) if cdr_summary_path.exists() else None
            for _, ab_row in df.iterrows():
                rep = int(ab_row['Replica'])
                row = {
                    'design_group': design_group, 'system_name': system_name,
                    'sim_replica':  rep,
                    'duration_ns':  float('nan'),
                    'n_frames_stats': int((pd.read_csv(
                        output_dir / f'RMSD_antibody_replica_{rep}.csv'
                    ).query(f'Time_ns > {stats_skip_ns}').shape[0])
                    ) if (output_dir / f'RMSD_antibody_replica_{rep}.csv').exists() else float('nan'),
                    'Ab_RMSD_mean_A': float(ab_row['Mean_Ab_RMSD_A']),
                    'Ab_RMSD_std_A':  float(ab_row['Std_Ab_RMSD_A']),
                    'Ab_RMSD_max_A':  float(ab_row['Max_Ab_RMSD_A']),
                    'Ab_RMSD_min_A':  float(ab_row['Min_Ab_RMSD_A']),
                    'status': 'ok (cached)',
                }
                # Read CDR stats if available
                for cdr in _CDR_ORDER:
                    col_base = cdr.replace('-', '_')
                    col_mean = f'{cdr}_mean_A'
                    col_std  = f'{cdr}_std_A'
                    if cdr_df is not None:
                        cdr_r = cdr_df[cdr_df['Replica'] == rep]
                        if not cdr_r.empty and col_mean in cdr_r.columns:
                            row[f'{col_base}_mean_A'] = float(cdr_r.iloc[0][col_mean])
                            row[f'{col_base}_std_A']  = float(cdr_r.iloc[0][col_std])
                        else:
                            row[f'{col_base}_mean_A'] = float('nan')
                            row[f'{col_base}_std_A']  = float('nan')
                    else:
                        row[f'{col_base}_mean_A'] = float('nan')
                        row[f'{col_base}_std_A']  = float('nan')
                rows.append(row)
        except Exception as e:
            print(f'    ⚠ Could not read cached results: {e}', flush=True)
        return rows

    # ── Load all simulation replicas ──────────────────────────────────────────
    universes   = {}
    load_status = {}

    for rep_num in _ACTIVE_REPLICAS:
        key       = f'replica_{rep_num}'
        rep_dir   = system_dir / f'replica_{rep_num}'
        topo, topo_type = find_topology(rep_dir, rep_num) if rep_dir.exists() else (None, None)
        traj            = find_trajectory(rep_dir, rep_num) if rep_dir.exists() else None

        if not rep_dir.exists():
            load_status[key] = 'missing_dir'
            continue
        if topo is None:
            load_status[key] = 'missing_tpr'
            continue
        if traj is None:
            load_status[key] = 'missing_xtc'
            continue

        try:
            u = mda.Universe(str(topo), str(traj))
            universes[key]   = u
            load_status[key] = 'ok'
            labels = {'tpr': 'TPR', 'pdb': 'PDB', 'gro_npt': 'GRO/npt', 'gro_md': 'GRO/md'}
            print(f'    ✓ {key}: {topo.name} [{labels.get(topo_type)}] + {traj.name}'
                  f'  {len(u.trajectory):,} frames'
                  f'  {u.trajectory.totaltime/1000:.1f} ns', flush=True)
        except Exception as e:
            load_status[key] = 'load_error'
            print(f'    ✗ {key}: load failed — {e}', flush=True)

    if not universes:
        print(f'  ✗ No replicas loaded — skipping.', flush=True)
        return [_make_failed_row(design_group, system_name, r, load_status.get(f'replica_{r}', 'missing_dir'))
                for r in _ACTIVE_REPLICAS]

    # ── Build shared selection strings ────────────────────────────────────────
    ref_u = next(iter(universes.values()))
    try:
        fit_sel   = build_chain_selection(ref_u, _PHLA_CHAINS)
        group_sel = build_chain_selection(ref_u, _ANTIBODY_CHAINS)
    except ValueError as e:
        print(f'  ✗ Chain selection failed: {e}', flush=True)
        return [_make_failed_row(design_group, system_name, r, 'chain_sel_error')
                for r in _ACTIVE_REPLICAS]

    # ── Auto-detect CDR regions ───────────────────────────────────────────────
    pdb_path    = system_dir / 'setup' / 'processed.pdb'
    cdr_residues = {}
    if pdb_path.exists():
        try:
            cdr_residues = remap_cdr_to_universe(ref_u, pdb_path)
            if cdr_residues:
                print(f'    CDR regions detected: {list(cdr_residues.keys())}', flush=True)
            else:
                print(f'    ⚠ No REMARK PDBinfo-LABEL in processed.pdb', flush=True)
        except Exception as e:
            print(f'    ⚠ CDR detection failed: {e}', flush=True)
    else:
        print(f'    ⚠ setup/processed.pdb not found — CDR analysis skipped', flush=True)

    cdr_names = [c for c in _CDR_ORDER if c in cdr_residues]

    # ── Apply PBC corrections ─────────────────────────────────────────────────
    for key, u in universes.items():
        print(f'    PBC: {key}', flush=True)
        try:
            universes[key] = remove_pbc(u)
        except Exception as e:
            print(f'    ⚠ PBC failed for {key}: {e}', flush=True)

    # ── Compute RMSD for each replica ─────────────────────────────────────────
    rmsd_results     = {}
    cdr_rmsd_results = {}

    for key, u in universes.items():
        print(f'    RMSD: {key}', flush=True)

        # Antibody RMSD
        try:
            rmsd_results[key] = compute_rmsd_antibody(
                u, fit_sel, group_sel,
                _REF_FRAME, 0, _RMSD_STEP, stats_skip_ns,
            )
        except Exception as e:
            print(f'    ✗ Ab RMSD failed: {e}', flush=True)
            rmsd_results[key] = None

        # CDR RMSD
        if cdr_residues:
            try:
                cdr_rmsd_results[key] = compute_cdr_rmsd(
                    u, fit_sel, cdr_residues,
                    _REF_FRAME, 0, _RMSD_STEP, stats_skip_ns,
                )
            except Exception as e:
                print(f'    ✗ CDR RMSD failed: {e}', flush=True)
                cdr_rmsd_results[key] = None

    # ── Save per-design CSVs ──────────────────────────────────────────────────
    # Guard: only write if at least one replica succeeded; otherwise no sentinel
    # is created and the design will be re-attempted on the next run.
    if any(v is not None for v in rmsd_results.values()):
        save_rmsd_antibody(rmsd_results, output_dir, stats_skip_ns)
    else:
        print('    ⚠ All RMSD computations failed — sentinel not written, will retry on next run.',
              flush=True)
    if cdr_rmsd_results and any(v is not None for v in cdr_rmsd_results.values()):
        save_cdr_rmsd(cdr_rmsd_results, cdr_names, output_dir)

    # ── Save plots (identical to notebook) ───────────────────────────────────
    if any(v is not None for v in rmsd_results.values()):
        png = plot_rmsd_antibody(rmsd_results, system_name, output_dir)
        print(f'    Saved: {png.name}', flush=True)

    if cdr_rmsd_results and cdr_names and any(v is not None for v in cdr_rmsd_results.values()):
        png = plot_cdr_rmsd(cdr_rmsd_results, cdr_names, system_name, output_dir)
        print(f'    Saved: {png.name}', flush=True)

    # ── Build master CSV rows ─────────────────────────────────────────────────
    rows = []
    for rep_num in _ACTIVE_REPLICAS:
        key = f'replica_{rep_num}'
        ab  = rmsd_results.get(key)
        cdr = cdr_rmsd_results.get(key)

        if ab is None:
            st = load_status.get(key, 'failed')
            rows.append(_make_failed_row(design_group, system_name, rep_num, st))
            continue

        row = {
            'design_group':   design_group,
            'system_name':    system_name,
            'sim_replica':    rep_num,
            'duration_ns':    ab['duration_ns'],
            'n_frames_stats': ab['n_stats'],
            'Ab_RMSD_mean_A': ab['mean'],
            'Ab_RMSD_std_A':  ab['std'],
            'Ab_RMSD_max_A':  ab['max'],
            'Ab_RMSD_min_A':  ab['min'],
            'status':         'ok',
        }
        for cdr_name in _CDR_ORDER:
            col_base = cdr_name.replace('-', '_')
            if cdr is not None and cdr_name in cdr:
                row[f'{col_base}_mean_A'] = cdr[cdr_name]['mean']
                row[f'{col_base}_std_A']  = cdr[cdr_name]['std']
            else:
                row[f'{col_base}_mean_A'] = float('nan')
                row[f'{col_base}_std_A']  = float('nan')
        rows.append(row)

    return rows


# ═══════════════════════════════════════════════════════════════════════════════
# Discovery
# ═══════════════════════════════════════════════════════════════════════════════

def discover_designs(base_dir, design_group_dirs=None):
    """Return list of (design_group, system_name) tuples found under base_dir."""
    if design_group_dirs is None:
        design_group_dirs = _DESIGN_GROUPS
    designs = []
    for grp in design_group_dirs:
        grp_path = base_dir / grp
        if not grp_path.exists():
            print(f'[WARN] Design group dir not found: {grp_path}', flush=True)
            continue
        found = sorted(d.name for d in grp_path.iterdir() if d.is_dir())
        print(f'  {grp}: {len(found)} designs', flush=True)
        for name in found:
            designs.append((grp, name))
    return designs


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Batch RMSD analysis for all MD designs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '--base-dir', type=Path, default=_DEFAULT_BASE_DIR,
        help='Root directory containing replica1/, replica2/, replica3/ sub-folders.',
    )
    parser.add_argument(
        '--stats-skip-ns', type=float, default=_DEFAULT_STATS_SKIP_NS,
        help='Equilibration cutoff: statistics are computed for t > this value (ns).',
    )
    parser.add_argument(
        '--design-groups', type=str, default=','.join(_DESIGN_GROUPS),
        help='Comma-separated list of design group dirs to process.',
    )
    parser.add_argument(
        '--force', action='store_true',
        help='Re-analyse even if output already exists.',
    )
    args = parser.parse_args()

    base_dir      = args.base_dir.expanduser().resolve()
    design_groups = [g.strip() for g in args.design_groups.split(',')]

    _setup_plot_style()

    print('=' * 72, flush=True)
    print('batch_analysis.py', flush=True)
    print(f'  MDAnalysis : {mda.__version__}', flush=True)
    print(f'  Base dir   : {base_dir}', flush=True)
    print(f'  Groups     : {design_groups}', flush=True)
    print(f'  Stats skip : t > {args.stats_skip_ns} ns', flush=True)
    print(f'  Force redo : {args.force}', flush=True)
    print('=' * 72, flush=True)

    if not base_dir.exists():
        print(f'ERROR: base_dir does not exist: {base_dir}', flush=True)
        sys.exit(1)

    # Discover all designs
    print('\nDiscovering designs:', flush=True)
    all_designs = discover_designs(base_dir, design_groups)
    print(f'Total: {len(all_designs)} designs × {len(_ACTIVE_REPLICAS)} replicas'
          f' = up to {len(all_designs)*len(_ACTIVE_REPLICAS)} rows\n', flush=True)

    # Process each design and collect master CSV rows
    all_rows   = []
    n_ok       = 0
    n_skip     = 0
    n_fail     = 0
    n_total    = len(all_designs)

    for idx, (design_group, system_name) in enumerate(all_designs, 1):
        print(f'[{idx:3d}/{n_total}] {design_group}/{system_name}', flush=True)
        try:
            rows = process_design(
                design_group, system_name, base_dir,
                stats_skip_ns=args.stats_skip_ns,
                force=args.force,
            )
            all_rows.extend(rows)
            for r in rows:
                st = str(r.get('status', ''))
                if 'cached' in st:
                    n_skip += 1
                elif 'ok' in st:
                    n_ok += 1
                else:
                    n_fail += 1
        except Exception:
            print(f'  ✗ Unhandled exception for {design_group}/{system_name}:', flush=True)
            traceback.print_exc()
            for rep in _ACTIVE_REPLICAS:
                all_rows.append(_make_failed_row(design_group, system_name, rep, 'exception'))
            n_fail += len(_ACTIVE_REPLICAS)
        print('', flush=True)  # blank line between designs

    # Save master CSV
    summary_dir = base_dir / 'batch_summary'
    summary_dir.mkdir(parents=True, exist_ok=True)
    master_csv = summary_dir / 'master_rmsd_summary.csv'

    # Build a tidy column order
    col_order = [
        'design_group', 'system_name', 'sim_replica',
        'duration_ns', 'n_frames_stats',
        'Ab_RMSD_mean_A', 'Ab_RMSD_std_A', 'Ab_RMSD_max_A', 'Ab_RMSD_min_A',
    ]
    for cdr in _CDR_ORDER:
        col_base = cdr.replace('-', '_')
        col_order += [f'{col_base}_mean_A', f'{col_base}_std_A']
    col_order.append('status')

    master_df = pd.DataFrame(all_rows)
    # Ensure all expected columns exist (fill missing with NaN)
    for col in col_order:
        if col not in master_df.columns:
            master_df[col] = float('nan')
    master_df = master_df[col_order]
    master_df.to_csv(master_csv, index=False, float_format='%.3f')

    print('=' * 72, flush=True)
    print(f'DONE', flush=True)
    print(f'  Rows written : {len(master_df)}', flush=True)
    print(f'  Status ok    : {n_ok}', flush=True)
    print(f'  Status fail  : {n_fail}', flush=True)
    print(f'  Status skip  : {n_skip}', flush=True)
    print(f'  Master CSV   : {master_csv}', flush=True)
    print('=' * 72, flush=True)


if __name__ == '__main__':
    main()
