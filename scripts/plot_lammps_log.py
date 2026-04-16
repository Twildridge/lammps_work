#!/usr/bin/env python3
"""
plot_lammps_log.py  —  Convergence and shear-diagnostics plots for LAMMPS runs.

Usage (general / slab_with_flow / slab_with_support):
    python plot_lammps_log.py <folder> <dataname>

Usage (shear_slab — output files include interaction + nsteps in suffix):
    python plot_lammps_log.py <folder> <dataname> \\
        --run-id <dataname>_<interaction>_<nsteps>

<folder>   : run directory containing log.lammps
<dataname> : base name used as plot title and for legacy file lookups
--run-id   : full output-file suffix for shear_slab runs
             (dataname_interaction_nsteps).  Defaults to <dataname>.

Output (saved inside <folder>/output_plots/convergence_plots/):
    {run_id}_convergence.png       — T, P, box/gel volume
    {run_id}_shear_diagnostics.png — shear-specific panels (auto-generated
                                     when shear output files or Pxz/Xz thermo
                                     columns are detected)
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')          # headless — safe on HPC clusters
import matplotlib.pyplot as plt
import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# FILE READERS
# ─────────────────────────────────────────────────────────────────────────────

def read_volume_file(filepath):
    """Single-column volume data (legacy format)."""
    data = []
    with open(filepath) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    data.append(float(line.split()[0]))
                except (ValueError, IndexError):
                    pass
    return np.array(data)


def read_timestep_volume_file(filepath):
    """Two-column step + value file (legacy format). Returns (steps, values)."""
    steps, vals = [], []
    with open(filepath) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    p = line.split()
                    steps.append(float(p[0]))
                    vals.append(float(p[1]))
                except (ValueError, IndexError):
                    pass
    return np.array(steps), np.array(vals)


def read_fix_print(filepath):
    """N-column fix print output (no header). Returns 2-D array (nrows x ncols)."""
    rows = []
    with open(filepath) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    rows.append([float(x) for x in line.split()])
                except ValueError:
                    pass
    return np.array(rows) if rows else np.empty((0, 0))


def read_ave_time(filepath):
    """fix ave/time output (2-line # header, then data).
    Returns (timesteps_1d, data_2d) where data_2d has shape (nrows, ncols-1).
    Column order matches the fix definition: xx yy zz xy xz yz."""
    rows = []
    with open(filepath) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    rows.append([float(x) for x in line.split()])
                except ValueError:
                    pass
    if not rows:
        return np.array([]), np.empty((0, 0))
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1:]


# ─────────────────────────────────────────────────────────────────────────────
# LOG PARSER
# ─────────────────────────────────────────────────────────────────────────────

def parse_lammps_log(filepath):
    """Parse all thermo blocks in log.lammps. Returns dict of np.arrays.

    Handles multiple run blocks (Phase 1 / Phase 2 / Phase 3) by appending
    rows from each block under the same keys."""
    data = {}
    reading = False
    headers = []

    with open(filepath, encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Step'):
                headers = line.split()
                reading = True
                for h in headers:
                    if h not in data:
                        data[h] = []
                continue
            if reading and ('Loop time' in line or line.startswith('WARNING')):
                reading = False
                continue
            if reading and line and not line.startswith('#'):
                try:
                    vals = line.split()
                    if len(vals) == len(headers):
                        for h, v in zip(headers, vals):
                            data[h].append(float(v))
                except ValueError:
                    pass

    return {k: np.array(v) for k, v in data.items()}


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _annotate_last30(ax, y, fmt='.3f', color='r'):
    """Dashed mean ± std line over the last 30% of y."""
    n = max(int(len(y) * 0.3), 1)
    if n < 5:
        return
    m, s = y[-n:].mean(), y[-n:].std()
    ax.axhline(m, color=color, ls='--', lw=1.2,
               label=f'Last 30%: {m:{fmt}} ± {s:{fmt}}')
    ax.legend(fontsize=8)


def _first_key(d, *keys):
    """Return the first key from *keys that exists in dict d, or None."""
    return next((k for k in keys if k in d), None)


# ─────────────────────────────────────────────────────────────────────────────
# CONVERGENCE PLOT  (T, P, box/gel volumes)
# ─────────────────────────────────────────────────────────────────────────────

def plot_convergence(data, folder, run_id, output):
    """Standard convergence figure for any sim type.

    Reads output files using run_id as the filename suffix (= dataname for
    non-shear sims; = dataname_interaction_nsteps for shear_slab)."""
    vd = os.path.join(folder, 'output_files', 'volume_data')

    box_file      = os.path.join(vd, f'box_dimensions_{run_id}.dat')
    gel_bb_file   = os.path.join(vd, f'gel_volume_bb_{run_id}.dat')
    gel_rg_file   = os.path.join(vd, f'gel_volume_rg_{run_id}.dat')
    num_dens_file = os.path.join(vd, f'num_density_{run_id}.dat')

    has_box      = os.path.exists(box_file)
    has_gel_bb   = os.path.exists(gel_bb_file)
    has_gel_rg   = os.path.exists(gel_rg_file)
    has_num_dens = os.path.exists(num_dens_file)

    n = 2 + has_box + has_gel_bb + has_gel_rg + has_num_dens
    fig, axes = plt.subplots(n, 1, figsize=(10, 3 * n), sharex=False)
    if n == 1:
        axes = [axes]
    fig.suptitle(run_id, fontsize=13, fontweight='bold')
    idx = 0
    steps = data.get('Step', np.array([]))

    # ── Temperature (supports both standard 'Temp' and shear_slab 'c_mobile_temp')
    temp_key = _first_key(data, 'Temp', 'c_mobile_temp')
    if temp_key:
        ax = axes[idx]; idx += 1
        ax.plot(steps, data[temp_key], 'b-', lw=1.5)
        ax.set_ylabel('Temperature')
        ax.grid(alpha=0.3)
        _annotate_last30(ax, data[temp_key])

    # ── Pressure
    press_key = _first_key(data, 'Press', 'c_mobile_press')
    if press_key:
        ax = axes[idx]; idx += 1
        ax.plot(steps, data[press_key], 'g-', lw=1.5)
        ax.set_ylabel('Pressure')
        ax.grid(alpha=0.3)
        _annotate_last30(ax, data[press_key])

    # ── Box volume (reads 4- or 7-column fix print box_dimensions file)
    if has_box:
        arr = read_fix_print(box_file)
        if arr.size and arr.shape[1] >= 4:
            t    = arr[:, 0]
            vols = arr[:, 1] * arr[:, 2] * arr[:, 3]   # lx * ly * lz
            is_pure = 'pure' in run_id.lower()
            y   = vols if is_pure else vols / vols[0]
            lbl = 'Box Volume (σ³)' if is_pure else 'Box Volume / Initial'
            ax = axes[idx]; idx += 1
            ax.plot(t, y, 'm-', lw=1.5)
            ax.set_ylabel(lbl)
            ax.grid(alpha=0.3)
            _annotate_last30(ax, y, fmt='.3g' if is_pure else '.4f')

    # ── Number density
    if has_num_dens:
        ts, nd = read_timestep_volume_file(num_dens_file)
        if ts.size:
            ax = axes[idx]; idx += 1
            ax.plot(ts, nd, color='darkorange', lw=1.5)
            ax.set_ylabel('Number Density (σ⁻³)')
            ax.grid(alpha=0.3)
            _annotate_last30(ax, nd, fmt='.4f')

    # ── Gel volume — bounding box (normalized)
    if has_gel_bb:
        vols = read_volume_file(gel_bb_file)
        if vols.size:
            y = vols / vols[0]
            ax = axes[idx]; idx += 1
            ax.plot(np.arange(len(y)), y, color='orange', lw=1.5)
            ax.set_ylabel('Gel Vol (BB) / Initial')
            ax.grid(alpha=0.3)
            _annotate_last30(ax, y)

    # ── Gel volume — Rg-based (normalized)
    if has_gel_rg:
        ts, vols = read_timestep_volume_file(gel_rg_file)
        if ts.size:
            y = vols / vols[0]
            ax = axes[idx]; idx += 1
            ax.plot(ts, y, color='cyan', lw=1.5)
            ax.set_ylabel('Gel Vol (Rg³) / Initial')
            ax.grid(alpha=0.3)
            _annotate_last30(ax, y)

    axes[-1].set_xlabel('Step')
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Convergence plot    → {output}")


# ─────────────────────────────────────────────────────────────────────────────
# SHEAR DIAGNOSTICS PLOT
# ─────────────────────────────────────────────────────────────────────────────

def plot_shear_diagnostics(data, folder, run_id, output):
    """Shear-specific diagnostic figure.

    Panels (shown when data is available):
      1. Temperature          — c_mobile_temp from thermo
      2. Bulk Pxz             — total xz stress component (all 3 phases)
      3. Box xz tilt          — Xz column from thermo; flat/ramp/flat pattern
      4. Polymer σ_p_xz       — time-averaged partial xz stress (Phase 3)
      5. Gel Rg dimensions    — lx_rg, ly_rg, lz_rg vs step

    These collectively answer: Did shear apply cleanly? Is the polymer stress
    response well-converged? Is the gel stable throughout?"""

    sd = os.path.join(folder, 'output_files', 'stress_data')
    vd = os.path.join(folder, 'output_files', 'volume_data')

    stress_p_file = os.path.join(sd, f'stress_tensor_polymer_{run_id}.dat')
    rg_file       = os.path.join(vd, f'gel_dimensions_rg_{run_id}.dat')

    # Decide which panels to draw
    temp_key = _first_key(data, 'c_mobile_temp', 'Temp')
    panels = []
    if temp_key:                          panels.append('temp')
    if 'Pxz' in data:                     panels.append('pxz')
    if 'Xz'  in data:                     panels.append('xz_tilt')
    if os.path.exists(stress_p_file):     panels.append('sigma_p_xz')
    if os.path.exists(rg_file):           panels.append('rg')

    if not panels:
        print("No shear diagnostic data found — skipping shear plot.")
        return

    fig, axes = plt.subplots(len(panels), 1, figsize=(10, 3 * len(panels)))
    if len(panels) == 1:
        axes = [axes]
    fig.suptitle(f'{run_id}  —  shear diagnostics', fontsize=12, fontweight='bold')
    steps = data.get('Step', np.array([]))

    for i, panel in enumerate(panels):
        ax = axes[i]
        ax.grid(alpha=0.3)

        # ── Temperature ──────────────────────────────────────────────────────
        if panel == 'temp':
            ax.plot(steps, data[temp_key], 'b-', lw=1.2)
            ax.set_ylabel('Temperature')
            _annotate_last30(ax, data[temp_key])

        # ── Bulk Pxz (all phases) ─────────────────────────────────────────
        elif panel == 'pxz':
            ax.plot(steps, data['Pxz'], color='steelblue', lw=1.2)
            ax.axhline(0, color='k', ls=':', lw=0.8)
            ax.set_ylabel('Pxz (bulk, thermo)')
            # Last 30% = Phase 3 production → annotate to preview G
            _annotate_last30(ax, data['Pxz'])

        # ── Box xz tilt (all phases) ──────────────────────────────────────
        elif panel == 'xz_tilt':
            ax.plot(steps, data['Xz'], color='darkgreen', lw=1.5)
            ax.set_ylabel('xz tilt (box)')
            # Expect: ~0 in Phase 1, linear ramp in Phase 2, flat in Phase 3
            _annotate_last30(ax, data['Xz'], fmt='.4f')

        # ── Polymer σ_p_xz (Phase 3 time-averaged, from fix ave/time) ────
        elif panel == 'sigma_p_xz':
            ts, arr = read_ave_time(stress_p_file)
            # ave/time columns (matching fix definition): xx yy zz xy xz yz
            #   index:                                     0  1  2  3  4  5
            if ts.size and arr.shape[1] >= 5:
                sigma_xz = arr[:, 4]
                ax.plot(ts, sigma_xz, color='crimson', lw=1.5)
                ax.axhline(0, color='k', ls=':', lw=0.8)
                ax.set_ylabel('σ_p_xz (polymer partial)')
                _annotate_last30(ax, sigma_xz)
                # Annotate plateau mean prominently — this feeds directly into G
                n = max(int(len(sigma_xz) * 0.3), 1)
                if n >= 5:
                    plateau = sigma_xz[-n:].mean()
                    ax.text(0.02, 0.92,
                            f'Plateau mean: {plateau:.5f}  →  G = σ_p_xz / γ',
                            transform=ax.transAxes, fontsize=9,
                            va='top', color='crimson',
                            bbox=dict(boxstyle='round,pad=0.3',
                                      facecolor='white', alpha=0.7))

        # ── Gel Rg dimensions ─────────────────────────────────────────────
        elif panel == 'rg':
            arr = read_fix_print(rg_file)
            # fix print columns: step  lx_rg  ly_rg  lz_rg
            if arr.size and arr.shape[1] >= 4:
                t = arr[:, 0]
                ax.plot(t, arr[:, 1], label='lx_Rg', color='royalblue',  lw=1.5)
                ax.plot(t, arr[:, 2], label='ly_Rg', color='darkorange', lw=1.5)
                ax.plot(t, arr[:, 3], label='lz_Rg', color='forestgreen', lw=1.5)
                ax.set_ylabel('Gel Rg dims (σ)')
                ax.legend(fontsize=8)
                # Stable Rg → gel not melting or grossly deforming under shear

    axes[-1].set_xlabel('Step')
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Shear diagnostics   → {output}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot LAMMPS log convergence and (optionally) shear diagnostics.')
    parser.add_argument('folder',
                        help='Run directory containing log.lammps')
    parser.add_argument('dataname',
                        help='Base data name (plot title; file lookup key for '
                             'non-shear sims)')
    parser.add_argument('--run-id', dest='run_id', default=None,
                        help='Full output-file suffix for shear_slab runs: '
                             'dataname_interaction_nsteps.  '
                             'Defaults to <dataname> if omitted.')
    args = parser.parse_args()

    run_id  = args.run_id or args.dataname
    logfile = os.path.join(args.folder, 'log.lammps')

    data = parse_lammps_log(logfile)
    if not data:
        print(f'No thermo data found in {logfile}')
        sys.exit(1)

    out_dir = os.path.join(args.folder, 'output_plots', 'convergence_plots')
    os.makedirs(out_dir, exist_ok=True)

    # Convergence plot (always)
    plot_convergence(data, args.folder, run_id,
                     os.path.join(out_dir, f'{run_id}_convergence.png'))

    # Shear diagnostics (auto-triggered by presence of shear output files or
    # shear-specific thermo columns Pxz / Xz)
    sd = os.path.join(args.folder, 'output_files', 'stress_data')
    shear_files_present = (
        os.path.exists(os.path.join(sd, f'stress_tensor_polymer_{run_id}.dat')) or
        os.path.exists(os.path.join(sd, f'shear_strain_{run_id}.dat'))
    )
    if shear_files_present or 'Pxz' in data or 'Xz' in data:
        plot_shear_diagnostics(data, args.folder, run_id,
                               os.path.join(out_dir, f'{run_id}_shear_diagnostics.png'))
