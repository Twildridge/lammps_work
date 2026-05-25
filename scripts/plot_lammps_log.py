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
    {run_id}_flow_diagnostics.png  — piston force + pressure panels (auto-generated
                                     when piston_force_*.dat is present;
                                     mode detected from which files exist)
    {run_id}_chempot_diagnostics.png — Widom μ_ex(z) profiles + solvent density
                                       (auto-generated when
                                       output_files/chemical_potential/ files present;
                                       produced by slab_with_support runs)
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


def read_ave_chunk(filepath):
    """Parse LAMMPS fix ave/chunk output.

    The file format written by LAMMPS is:
        # comment lines (start with #, skipped)
        timestep  n_chunks  total_count       ← 3 numeric tokens (block header)
        chunk_id  coord  Ncount  val1 ...
        ...
        timestep  n_chunks  total_count
        ...

    Note: the block header has THREE fields (timestep, n_chunks, total_count),
    not two — a previous version of this reader incorrectly checked len==2 and
    therefore never matched any header, producing an empty frame list.

    Returns a list of (timestep, coords, values) tuples where:
        timestep : float — the output step
        coords   : 1-D array of length n_chunks — bin-centre coordinates
        values   : 2-D array (n_chunks, n_cols) — Ncount, then each property
    """
    frames = []
    with open(filepath) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]

    i = 0
    while i < len(lines):
        parts = lines[i].split()
        # Block header: exactly THREE tokens (timestep  n_chunks  total_count)
        if len(parts) == 3:
            try:
                ts    = float(parts[0])
                nrows = int(parts[1])
                # parts[2] is total_count — read but not used
                i += 1
                coords, values = [], []
                for _ in range(nrows):
                    if i < len(lines):
                        row = [float(x) for x in lines[i].split()]
                        # row[0]=chunk_id  row[1]=coord  row[2]=Ncount  row[3+]=props
                        if len(row) >= 3:
                            coords.append(row[1])
                            values.append(row[2:])  # Ncount + computed properties
                        i += 1
                if coords:
                    frames.append((ts, np.array(coords), np.array(values)))
                continue
            except (ValueError, IndexError):
                pass
        i += 1
    return frames


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

    box_file         = os.path.join(vd, f'box_dimensions_{run_id}.dat')
    gel_bb_file      = os.path.join(vd, f'gel_volume_bb_{run_id}.dat')
    gel_rg_file      = os.path.join(vd, f'gel_volume_rg_{run_id}.dat')
    gel_dims_rg_file = os.path.join(vd, f'gel_dimensions_rg_{run_id}.dat')
    num_dens_file    = os.path.join(vd, f'num_density_{run_id}.dat')

    has_box         = os.path.exists(box_file)
    has_gel_bb      = os.path.exists(gel_bb_file)
    has_gel_rg      = os.path.exists(gel_rg_file)
    has_gel_dims_rg = os.path.exists(gel_dims_rg_file)
    has_num_dens    = os.path.exists(num_dens_file)

    n = 2 + has_box + has_gel_bb + has_gel_rg + has_gel_dims_rg + has_num_dens
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
        ax.plot(steps, data[temp_key], 'b-', lw=1.5, marker='o', markersize=3)
        ax.set_ylabel('Temperature')
        ax.grid(alpha=0.3)
        _annotate_last30(ax, data[temp_key])

    # ── Pressure
    press_key = _first_key(data, 'Press', 'c_mobile_press')
    if press_key:
        ax = axes[idx]; idx += 1
        ax.plot(steps, data[press_key], 'g-', lw=1.5, marker='o', markersize=3)
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
            ax.plot(t, y, 'm-', lw=1.5, marker='o', markersize=3)
            ax.set_ylabel(lbl)
            ax.grid(alpha=0.3)
            _annotate_last30(ax, y, fmt='.3g' if is_pure else '.4f')

    # ── Number density
    if has_num_dens:
        ts, nd = read_timestep_volume_file(num_dens_file)
        if ts.size:
            ax = axes[idx]; idx += 1
            ax.plot(ts, nd, color='darkorange', lw=1.5, marker='o', markersize=3)
            ax.set_ylabel('Number Density (σ⁻³)')
            ax.grid(alpha=0.3)
            _annotate_last30(ax, nd, fmt='.4f')

    # ── Gel volume — bounding box (normalized)
    if has_gel_bb:
        vols = read_volume_file(gel_bb_file)
        if vols.size:
            y = vols / vols[0]
            ax = axes[idx]; idx += 1
            ax.plot(np.arange(len(y)), y, color='orange', lw=1.5, marker='o', markersize=3)
            ax.set_ylabel('Gel Vol (BB) / Initial')
            ax.grid(alpha=0.3)
            _annotate_last30(ax, y)

    # ── Gel volume — Rg-based single-value file (non-shear sims)
    if has_gel_rg:
        ts, vols = read_timestep_volume_file(gel_rg_file)
        if ts.size:
            y = vols / vols[0]
            ax = axes[idx]; idx += 1
            ax.plot(ts, y, color='cyan', lw=1.5, marker='o', markersize=3)
            ax.set_ylabel('Gel Vol (Rg³) / Initial')
            ax.grid(alpha=0.3)
            _annotate_last30(ax, y)

    # ── Gel volume — Rg dimensions file (shear_slab: gel_dimensions_rg_*.dat)
    # Columns: step  lx_rg  ly_rg  lz_rg  →  vol = product of last three
    if has_gel_dims_rg:
        arr = read_fix_print(gel_dims_rg_file)
        if arr.size and arr.shape[1] >= 4:
            t    = arr[:, 0]
            vols = arr[:, 1] * arr[:, 2] * arr[:, 3]
            y    = vols / vols[0]
            ax   = axes[idx]; idx += 1
            ax.plot(t, y, color='cyan', lw=1.5, marker='o', markersize=3)
            ax.axhline(1.0, color='k', ls=':', lw=0.8)
            ax.set_ylabel('Gel Vol (Rg) / Initial')
            ax.grid(alpha=0.3)
            _annotate_last30(ax, y)
            # Flag deswell: if final vol < 0.97 × initial, annotate warning
            if y[-1] < 0.97:
                ax.text(0.02, 0.05, f'⚠ Final vol/vol₀ = {y[-1]:.3f}  (deswell detected)',
                        transform=ax.transAxes, fontsize=9, color='red',
                        va='bottom',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

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
      3. Box xz tilt          — Xz column from thermo; step jump at Phase 2 start then flat
      4. Polymer σ_p_xz       — time-averaged partial xz stress (Phase 3)
      5. Gel Rg dimensions    — lx_rg, ly_rg, lz_rg vs step
      6. Gel strain (CM)      — gel_strain_cm from polymer COM (4-col format);
                                starts at ~target_strain_xz after change_box affine remap,
                                then partially relaxes during Phase 2 NVT (surface layers).
                                Interior network held at target_strain_xz by periodic topology.
                                G = <sigma_p_xz> / target_strain_xz (box strain, not CM).
                                or gel_strain_box vs gel_strain_cm comparison (5-col legacy)

    These collectively answer: Did shear apply cleanly? Is the polymer stress
    response well-converged? Is the gel stable throughout? Is the gel itself
    sheared by the target strain (not slipping relative to the box)?"""

    sd = os.path.join(folder, 'output_files', 'stress_data')
    vd = os.path.join(folder, 'output_files', 'volume_data')

    stress_p_file   = os.path.join(sd, f'stress_tensor_polymer_{run_id}.dat')
    rg_file         = os.path.join(vd, f'gel_dimensions_rg_{run_id}.dat')
    shear_str_file  = os.path.join(sd, f'shear_strain_{run_id}.dat')

    # Decide which panels to draw
    temp_key = _first_key(data, 'c_mobile_temp', 'Temp')
    panels = []
    if temp_key:                              panels.append('temp')
    if 'Pxz' in data:                         panels.append('pxz')
    if 'Xz'  in data:                         panels.append('xz_tilt')
    if os.path.exists(stress_p_file):         panels.append('sigma_p_xz')
    if os.path.exists(rg_file):               panels.append('rg')
    # Gel strain panel — col 4 = gel_strain_cm (4-col rheometer format)
    # or gel_strain_compare panel if old 5-col format (box + CM)
    if os.path.exists(shear_str_file):
        _ss = read_fix_print(shear_str_file)
        if _ss.size:
            if _ss.shape[1] >= 5:
                panels.append('gel_strain_compare')
            elif _ss.shape[1] >= 4:
                panels.append('gel_strain_cm')

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
            ax.plot(steps, data[temp_key], 'b-', lw=1.2, marker='o', markersize=3)
            ax.set_ylabel('Temperature')
            _annotate_last30(ax, data[temp_key])

        # ── Bulk Pxz (all phases) ─────────────────────────────────────────
        elif panel == 'pxz':
            ax.plot(steps, data['Pxz'], color='steelblue', lw=1.2, marker='o', markersize=3)
            ax.axhline(0, color='k', ls=':', lw=0.8)
            ax.set_ylabel('Pxz (bulk, thermo)')
            # Last 30% = Phase 3 production → annotate to preview G
            _annotate_last30(ax, data['Pxz'])

        # ── Box xz tilt (all phases) ──────────────────────────────────────
        elif panel == 'xz_tilt':
            ax.plot(steps, data['Xz'], color='darkgreen', lw=1.5, marker='o', markersize=3)
            ax.set_ylabel('xz tilt (box)')
            # Expect: ~0 in Phase 1, linear ramp in Phase 2, flat in Phase 3

        # ── Polymer σ_p_xz (Phase 3 time-averaged, from fix ave/time) ────
        elif panel == 'sigma_p_xz':
            ts, arr = read_ave_time(stress_p_file)
            # ave/time columns (matching fix definition): xx yy zz xy xz yz
            #   index:                                     0  1  2  3  4  5
            if ts.size and arr.shape[1] >= 5:
                sigma_xz = arr[:, 4]
                ax.plot(ts, sigma_xz, color='crimson', lw=1.5, marker='o', markersize=3)
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
                ax.plot(t, arr[:, 1], label='lx_Rg', color='royalblue',  lw=1.5, marker='o', markersize=3)
                ax.plot(t, arr[:, 2], label='ly_Rg', color='darkorange', lw=1.5, marker='o', markersize=3)
                ax.plot(t, arr[:, 3], label='lz_Rg', color='forestgreen', lw=1.5, marker='o', markersize=3)
                ax.set_ylabel('Gel Rg dims (σ)')
                ax.legend(fontsize=8)
                # Stable Rg → gel not melting or grossly deforming under shear

        # ── Gel strain (polymer COM, change_box) — 4-column shear_strain file ─
        elif panel == 'gel_strain_cm':
            arr = read_fix_print(shear_str_file)
            # cols: step  gel_lz_initial  gel_thick  gel_strain_cm
            # With change_box approach: gel_strain_cm starts near target_strain_xz
            # (immediately after affine remap) then partially relaxes during Phase 2 NVT.
            # The plateau is expected to be < target_strain_xz — that is physically correct
            # (surface layers relax; interior network held by periodic topology).
            # G = <sigma_p_xz> / target_strain_xz  (NOT gel_strain_cm).
            if arr.size and arr.shape[1] >= 4:
                t       = arr[:, 0]
                gs_cm   = arr[:, 3]
                # Infer target strain from first value (should be ≈ target after change_box remap)
                target_gamma = gs_cm[0] if len(gs_cm) else 0.10
                ax.plot(t, gs_cm, color='darkorange', lw=1.5, marker='s', markersize=3,
                        label='gel_strain_cm  (surface COM, partially relaxes)')
                ax.axhline(target_gamma, color='steelblue', ls='--', lw=1.0,
                           label=f'γ_box = {target_gamma:.3f}  (interior, fixed)')
                ax.axhline(0, color='k', ls=':', lw=0.8)
                ax.set_ylabel('Gel Shear Strain γ')
                ax.legend(fontsize=8)
                ax.text(0.02, 0.92,
                        f'Surface plateau: γ_CM = {gs_cm[-1]:.4f}\n'
                        f'Interior (box):  γ_box = {target_gamma:.4f}',
                        transform=ax.transAxes, fontsize=9, va='top',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

        # ── Gel strain comparison: box (xz/lz) vs polymer COM method ─────
        elif panel == 'gel_strain_compare':
            arr = read_fix_print(shear_str_file)
            # shear_strain_*.dat cols: step  gel_lz_initial  xz_tilt  gel_strain_box  gel_strain_cm
            if arr.size and arr.shape[1] >= 5:
                t        = arr[:, 0]
                gs_box   = arr[:, 3]   # xz / lz
                gs_cm    = arr[:, 4]   # polymer COM method
                ax.plot(t, gs_box, color='steelblue',  lw=1.5, marker='o', markersize=3,
                        label='gel_strain_box  (xz/lz)')
                ax.plot(t, gs_cm,  color='darkorange', lw=1.5, marker='s', markersize=3,
                        label='gel_strain_cm  (polymer COM)')
                ax.set_ylabel('Gel Shear Strain γ')
                ax.legend(fontsize=8)
                # Final values
                ax.text(0.02, 0.92,
                        f'Final:  box={gs_box[-1]:.4f}   CM={gs_cm[-1]:.4f}',
                        transform=ax.transAxes, fontsize=9, va='top',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
                # Flag discrepancy (|box − CM| > 0.01 = 1 percentage point)
                if abs(gs_box[-1] - gs_cm[-1]) > 0.01:
                    ax.text(0.02, 0.06,
                            f'⚠ box/CM mismatch: Δγ={abs(gs_box[-1]-gs_cm[-1]):.4f}',
                            transform=ax.transAxes, fontsize=9, color='red', va='bottom',
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    axes[-1].set_xlabel('Step')
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Shear diagnostics   → {output}")


# ─────────────────────────────────────────────────────────────────────────────
# FLOW / COMPRESSION DIAGNOSTICS PLOT
# ─────────────────────────────────────────────────────────────────────────────

def plot_flow_diagnostics(folder, run_id, output):
    """Piston force and reservoir-pressure diagnostics for slab_with_flow runs.

    Mode is auto-detected from which output files are present:
      Compression  — piston_force_*.dat present, no piston_force_pressure_*.dat
                     Panels: piston force | gel strain | piston position
      Permeation   — piston_force_pressure_*.dat also present
                     Panels: piston force (with Phase-0 baseline) |
                             F/A vs P_res_feed vs P_res_perm (instantaneous) |
                             time-averaged reservoir pressures

    File column reference
    ─────────────────────
    piston_force_*.dat              : step  F_piston_z
    piston_force_pressure_*.dat     : step  F_z  F/A  P_res_feed  P_res_perm  P_feed×A
    strain_zz_*.dat                 : step  L0   L_current
    pressure_feed_*.dat             : fix ave/time scalar → (ts, [[P_feed], ...])
    pressure_permeate_*.dat         : fix ave/time scalar → (ts, [[P_perm], ...])
    piston_position_*.dat           : step  piston_z
    permeate_count_*.dat            : step  N_permeate  (solvent atoms in permeate region)
    """
    piston_dir  = os.path.join(folder, 'output_files', 'piston_data')
    stress_dir  = os.path.join(folder, 'output_files', 'stress_data')
    perm_dir    = os.path.join(folder, 'output_files', 'permeation_data')

    force_file      = os.path.join(piston_dir, f'piston_force_{run_id}.dat')
    force_pres_file = os.path.join(piston_dir, f'piston_force_pressure_{run_id}.dat')
    pos_file        = os.path.join(piston_dir, f'piston_position_{run_id}.dat')
    strain_file     = os.path.join(stress_dir,  f'strain_zz_{run_id}.dat')
    p_feed_file     = os.path.join(stress_dir,  f'pressure_feed_{run_id}.dat')
    p_perm_file     = os.path.join(stress_dir,  f'pressure_permeate_{run_id}.dat')
    flux_file       = os.path.join(perm_dir,    f'permeate_count_{run_id}.dat')

    has_force      = os.path.exists(force_file)
    has_force_pres = os.path.exists(force_pres_file)
    has_pos        = os.path.exists(pos_file)
    has_strain     = os.path.exists(strain_file)
    has_p_feed     = os.path.exists(p_feed_file)
    has_p_perm     = os.path.exists(p_perm_file)
    has_flux       = os.path.exists(flux_file)

    if not has_force and not has_force_pres:
        return  # nothing to do

    is_permeation = has_force_pres

    # Build panel list
    panels = []
    if has_force:
        panels.append('force')
    if is_permeation:
        panels.append('force_vs_pres')       # F/A vs P_feed vs P_perm (instantaneous)
        if has_p_feed or has_p_perm:
            panels.append('pres_timeseries') # time-averaged reservoir pressures
        if has_flux:
            panels.append('solvent_flux')    # N_permeate vs step
    else:
        # compression-specific panels
        if has_strain:
            panels.append('strain')
        if has_pos:
            panels.append('piston_pos')

    if not panels:
        return

    mode_label = 'permeation' if is_permeation else 'compression'
    fig, axes = plt.subplots(len(panels), 1,
                             figsize=(10, 3.5 * len(panels)), sharex=False)
    if len(panels) == 1:
        axes = [axes]
    fig.suptitle(f'{run_id}  —  flow diagnostics ({mode_label})',
                 fontsize=12, fontweight='bold')

    for i, panel in enumerate(panels):
        ax = axes[i]
        ax.grid(alpha=0.3)

        # ── Piston force time series (both modes) ─────────────────────────
        if panel == 'force':
            arr = read_fix_print(force_file)
            if arr.size and arr.shape[1] >= 2:
                t, fz = arr[:, 0], arr[:, 1]
                ax.plot(t, fz, color='steelblue', lw=1.5, marker='o', markersize=2)
                ax.axhline(0, color='k', ls=':', lw=0.8)
                ax.set_ylabel('F_piston_z  (ε/σ)')
                if is_permeation:
                    ax.set_title(
                        'Piston z-force  (positive = upward).  '
                        'Phase-0 region is the zero-pressure baseline.',
                        fontsize=9)
                else:
                    ax.set_title(
                        'Gel reaction force on piston  '
                        '(rises during compression → plateaus at ~10 % strain → relaxes)',
                        fontsize=9)
                    _annotate_last30(ax, fz, fmt='.4f', color='darkblue')
                    # Mark peak force
                    pk_idx = np.argmax(fz)
                    ax.annotate(f'peak = {fz[pk_idx]:.4f}',
                                xy=(t[pk_idx], fz[pk_idx]),
                                xytext=(0.6, 0.92), textcoords='axes fraction',
                                fontsize=8, color='steelblue',
                                arrowprops=dict(arrowstyle='->', color='steelblue', lw=1),
                                bbox=dict(boxstyle='round,pad=0.2',
                                          facecolor='white', alpha=0.8))

        # ── F/A vs P_feed vs P_perm (instantaneous, permeation) ──────────
        elif panel == 'force_vs_pres':
            arr = read_fix_print(force_pres_file)
            # cols: step | F_z | F/A | P_res_feed | P_res_perm | P_feed×A
            if arr.size and arr.shape[1] >= 6:
                t        = arr[:, 0]
                f_over_a = arr[:, 2]   # F / (lx·ly)  — piston pressure
                p_feed   = arr[:, 3]   # feed reservoir pressure
                p_perm   = arr[:, 4]   # permeate reservoir pressure

                ax.plot(t, p_feed,   color='tomato',         lw=1.5, marker='o',
                        markersize=2, label='P_res_feed  (below piston)')
                ax.plot(t, p_perm,   color='cornflowerblue', lw=1.2, marker='o',
                        markersize=2, alpha=0.8, label='P_res_perm  (above piston)')
                ax.plot(t, f_over_a, color='steelblue',      lw=1.5, marker='s',
                        markersize=2, ls='--',
                        label='F_piston / A  (≈ P_feed − P_perm at quasi-equil)')
                ax.axhline(0, color='k', ls=':', lw=0.8)
                ax.set_ylabel('Pressure  (ε/σ³)')
                ax.set_title(
                    'Piston pressure vs reservoir pressures\n'
                    'F/A ≈ P_feed − P_perm validates that pressure drop = piston force / area',
                    fontsize=9)
                ax.legend(fontsize=8)

        # ── Time-averaged reservoir pressures (permeation, from fix ave/time) ──
        elif panel == 'pres_timeseries':
            plotted = False
            if has_p_feed:
                ts, dat = read_ave_time(p_feed_file)
                if ts.size and dat.shape[1] >= 1:
                    ax.plot(ts, dat[:, 0], color='tomato', lw=1.5, marker='o',
                            markersize=2, label='P_res_feed  (time-averaged)')
                    _annotate_last30(ax, dat[:, 0], fmt='.4f', color='tomato')
                    plotted = True
            if has_p_perm:
                ts, dat = read_ave_time(p_perm_file)
                if ts.size and dat.shape[1] >= 1:
                    ax.plot(ts, dat[:, 0], color='cornflowerblue', lw=1.5,
                            marker='o', markersize=2,
                            label='P_res_perm  (time-averaged)')
                    plotted = True
            if plotted:
                ax.set_ylabel('Reservoir Pressure  (ε/σ³)')
                ax.set_title(
                    'Time-averaged feed and permeate reservoir pressures\n'
                    'ΔP = P_feed − P_perm drives steady-state solvent flux',
                    fontsize=9)
                ax.legend(fontsize=8)

        # ── Gel compressive strain vs step (compression) ──────────────────
        elif panel == 'strain':
            arr = read_fix_print(strain_file)
            # cols: step | L0 | L_current
            if arr.size and arr.shape[1] >= 3:
                t   = arr[:, 0]
                L0  = arr[0, 1]          # initial Rg-based thickness (constant)
                L   = arr[:, 2]          # current thickness
                eps = (L0 - L) / L0 * 100.0  # strain in %
                ax.plot(t, eps, color='darkorange', lw=1.5, marker='o', markersize=2)
                ax.axhline(10.0, color='r', ls='--', lw=1.0, label='10 % target strain')
                ax.set_ylabel('Gel Strain ε  (%)')
                ax.set_title('Gel compressive strain (Rg-based thickness)', fontsize=9)
                ax.legend(fontsize=8)
                peak  = eps.max()
                final = eps[-1]
                ax.text(0.98, 0.92,
                        f'Peak: {peak:.2f}%   Final (post-relax): {final:.2f}%',
                        transform=ax.transAxes, fontsize=9, ha='right', va='top',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

        # ── Piston z-position vs step (compression) ───────────────────────
        elif panel == 'piston_pos':
            arr = read_fix_print(pos_file)
            # cols: step | piston_z
            if arr.size and arr.shape[1] >= 2:
                t, z = arr[:, 0], arr[:, 1]
                ax.plot(t, z, color='purple', lw=1.5, marker='o', markersize=2)
                ax.set_ylabel('Piston z  (σ)')
                ax.set_title('Piston centre-of-mass z position', fontsize=9)

        # ── Solvent flux: N_permeate vs step (permeation) ─────────────────
        elif panel == 'solvent_flux':
            arr = read_fix_print(flux_file)
            # cols: step | N_permeate
            if arr.size and arr.shape[1] >= 2:
                t, n_perm = arr[:, 0], arr[:, 1]
                ax.plot(t, n_perm, color='teal', lw=1.5, marker='o', markersize=2)
                ax.set_ylabel('N_permeate  (atoms)')
                ax.set_title(
                    'Solvent atoms in permeate region vs time\n'
                    'Slope = steady-state solvent flux  (atoms / τ)',
                    fontsize=9)
                # Annotate final count
                ax.text(0.98, 0.05,
                        f'Final count: {int(n_perm[-1])}',
                        transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    axes[-1].set_xlabel('Step')
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Flow diagnostics    → {output}")


# ─────────────────────────────────────────────────────────────────────────────
# CHEMICAL POTENTIAL DIAGNOSTICS PLOT
# ─────────────────────────────────────────────────────────────────────────────

def _read_cavity_frames(filepath):
    """Read mu_z_cavity_{stem}.dat (per-frame cavity Widom output).

    Columns: step  z_center  z_lo  z_hi  p_cav  n_cavity  n_trial  mu_ex

    Returns OrderedDict  {timestep (int): {'z', 'mu_ex', 'p_cav'} arrays}
    in chronological order.
    """
    from collections import OrderedDict
    frames = OrderedDict()
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                step  = int(float(parts[0]))
                z     = float(parts[1])
                p_cav = float(parts[4]) if parts[4] != 'nan' else np.nan
                mu    = float(parts[7]) if parts[7] != 'nan' else np.nan
            except ValueError:
                continue
            if step not in frames:
                frames[step] = {'z': [], 'mu_ex': [], 'p_cav': []}
            frames[step]['z'].append(z)
            frames[step]['mu_ex'].append(mu)
            frames[step]['p_cav'].append(p_cav)
    return {s: {k: np.array(v) for k, v in d.items()} for s, d in frames.items()}


def plot_chempot_diagnostics(folder, run_id, output, thermo_data=None):
    """Cavity-biased Widom chemical potential and pore pressure diagnostics.

    Standard Widom (fix widom) is dropped — it fails for WCA fluids at liquid
    densities.  All chemical-potential information comes from cavity_widom.py.

    Panels
    ──────
    1. Solvent number density ρ_s(z) — one curve per output frame.
       Reservoir ρ_s → V̄_s = 1/ρ_s,reservoir used in panels 2 & 3.

    2. Cavity-biased Widom μ_ex(z) — per-frame evolution curves (light,
       connected) overlaid with time-averaged ± stderr (thick, p_cav coloured).

    3. Pore pressure profile p_p(z)/P_ext — evolution curves.
       p_p(z) = μ_ex(z) / V̄_s = μ_ex(z) × ρ_s,reservoir
       Normalized by the external (barostat) pressure P_ext so that the
       reservoir value equals 1.  P_ext is taken from the mean thermo pressure
       of the latter half of the production run (or 1.5 ε/σ³ if unavailable).

    File column reference
    ─────────────────────
    solvent_density_z_*.dat  (fix ave/chunk) : timestep  n_chunks  total_count /
                                               chunk  coord  Ncount  density
    mu_z_cavity_*.dat        (cavity_widom)  : step z_center z_lo z_hi p_cav
                                               n_cavity n_trial mu_ex
    mu_z_cavity_summary_*.dat               : z_center z_lo z_hi mu_ex_mean
                                              mu_ex_stderr p_cav_mean n_frames
    """
    chem_dir         = os.path.join(folder, 'output_files', 'chemical_potential')
    dens_file        = os.path.join(chem_dir, f'solvent_density_z_{run_id}.dat')
    cavity_frame_file= os.path.join(chem_dir, f'mu_z_cavity_{run_id}.dat')
    cavity_sum_file  = os.path.join(chem_dir, f'mu_z_cavity_summary_{run_id}.dat')

    has_dens         = os.path.exists(dens_file)
    has_cavity_frame = os.path.exists(cavity_frame_file)
    has_cavity_sum   = os.path.exists(cavity_sum_file)

    if not has_dens and not has_cavity_sum:
        return

    # ── Infer P_ext ──────────────────────────────────────────────────────────
    # Use latter-half mean of thermo pressure; fall back to 1.5 LJ units.
    press_ext = 1.5
    if thermo_data is not None and 'Press' in thermo_data:
        p_arr = np.asarray(thermo_data['Press'], dtype=float)
        if len(p_arr) > 10:
            press_ext = float(p_arr[len(p_arr) // 2:].mean())

    # ── Determine which panels to draw ───────────────────────────────────────
    show_density   = has_dens
    show_cavity    = has_cavity_frame or has_cavity_sum
    show_pore_pres = has_dens and has_cavity_frame   # need both for p_p

    n_panels = int(show_density) + int(show_cavity) + int(show_pore_pres)
    if n_panels == 0:
        return

    fig, axes = plt.subplots(n_panels, 1, figsize=(10, 4.5 * n_panels),
                             sharex=False)
    if n_panels == 1:
        axes = [axes]
    fig.suptitle(f'{run_id}  —  chemical potential diagnostics',
                 fontsize=12, fontweight='bold')
    panel_idx = 0

    # ── Shared state: reservoir ρ_s ──────────────────────────────────────────
    rho_res  = None   # set in density panel; used by cavity & pore panels
    vbar_res = None

    # ── Panel 1: Solvent density ─────────────────────────────────────────────
    if show_density:
        dens_frames = read_ave_chunk(dens_file)
        ax          = axes[panel_idx]; panel_idx += 1
        ax.grid(alpha=0.3)

        if dens_frames:
            n_dens   = len(dens_frames)
            colors_d = plt.cm.plasma(np.linspace(0, 1, max(n_dens, 1)))
            rho_last = z_last = None

            for i, (ts, coords, values) in enumerate(dens_frames):
                if values.shape[1] < 2:
                    continue
                rho = values[:, 1]
                ax.plot(coords, rho, color=colors_d[i], lw=1.5,
                        marker='o', markersize=2,
                        label=f'step {int(ts):,}')
                if i == n_dens - 1:
                    rho_last = rho
                    z_last   = coords

            ax.set_ylabel('ρ_s  (σ⁻³)')
            ax.set_title(
                'Solvent number density profile  ρ_s(z)\n'
                'V̄_s = 1/ρ_s,reservoir  (partial molar volume, LJ units).\n'
                'Dip inside gel = polymer exclusion, not solvent compression.',
                fontsize=9)
            ax.legend(fontsize=7, ncol=min(n_dens, 5),
                      loc='upper right', framealpha=0.7)

            RHO_THRESH_FRAC = 0.85
            if rho_last is not None and z_last is not None and len(z_last) > 4:
                rho_max  = float(rho_last.max())
                res_mask = rho_last >= RHO_THRESH_FRAC * rho_max
                if res_mask.sum() >= 2 and rho_max > 1e-6:
                    rho_res  = float(rho_last[res_mask].mean())
                    vbar_res = 1.0 / rho_res
                    res_z    = z_last[res_mask]
                    ax.axvspan(float(res_z.min()), float(res_z.max()),
                               alpha=0.10, color='steelblue',
                               label=f'reservoir bins (ρ_s ≥ {RHO_THRESH_FRAC:.0%}×max)')
                    ax.text(0.98, 0.92,
                            f'Reservoir  ρ_s = {rho_res:.4f} σ⁻³'
                            f'  ({res_mask.sum()} bins)\n'
                            f'→  V̄_s = {vbar_res:.3f} σ³/atom'
                            f'   P_ext = {press_ext:.3f} ε/σ³',
                            transform=ax.transAxes, fontsize=9,
                            ha='right', va='top',
                            bbox=dict(boxstyle='round,pad=0.3',
                                      facecolor='lightyellow', alpha=0.9))
                    ax.legend(fontsize=7, ncol=min(n_dens, 5),
                              loc='upper left', framealpha=0.7)
                else:
                    ax.text(0.98, 0.92,
                            f'Could not isolate reservoir bins\n'
                            f'(< 2 bins above {RHO_THRESH_FRAC:.0%}×max ρ_s)',
                            transform=ax.transAxes, fontsize=9,
                            ha='right', va='top', color='red',
                            bbox=dict(boxstyle='round,pad=0.3',
                                      facecolor='white', alpha=0.85))
        else:
            ax.text(0.5, 0.5, 'No data in solvent_density_z file',
                    transform=ax.transAxes, ha='center', va='center')

    # ── Panel 2: Cavity-biased Widom μ_ex(z) evolution + summary ─────────────
    if show_cavity:
        ax = axes[panel_idx]; panel_idx += 1
        ax.grid(alpha=0.3)

        # --- Per-frame evolution curves (light, connected lines) -------------
        cav_frames = {}
        if has_cavity_frame:
            try:
                cav_frames = _read_cavity_frames(cavity_frame_file)
            except Exception:
                cav_frames = {}

        if cav_frames:
            steps       = sorted(cav_frames.keys())
            n_cf        = len(steps)
            colors_cf   = plt.cm.viridis(np.linspace(0, 1, max(n_cf, 1)))
            for i, step in enumerate(steps):
                d   = cav_frames[step]
                z   = d['z']
                mu  = d['mu_ex']
                ok  = ~np.isnan(mu)
                ax.plot(z[ok], mu[ok], color=colors_cf[i], lw=1.2,
                        alpha=0.55, marker='o', markersize=3,
                        label=f'step {step:,}')

        # --- Time-averaged ± stderr (thick, p_cav coloured) ------------------
        if has_cavity_sum:
            try:
                cav = np.genfromtxt(cavity_sum_file, comments='#')
            except Exception:
                cav = np.empty((0, 0))

            if cav.ndim == 2 and cav.shape[0] >= 2 and cav.shape[1] >= 6:
                z_c   = cav[:, 0]
                mu    = cav[:, 3]
                se    = cav[:, 4]
                p_cav = cav[:, 5]
                valid = ~np.isnan(mu)

                cmap   = plt.cm.viridis
                p_norm = plt.Normalize(
                    vmin=0.0,
                    vmax=max(float(p_cav[valid].max()) if valid.any() else 1.0, 0.01))

                # Connected line through valid time-averaged points
                ax.plot(z_c[valid], mu[valid], color='k', lw=1.0,
                        alpha=0.4, zorder=3)

                for j in np.where(valid)[0]:
                    color = cmap(p_norm(p_cav[j]))
                    se_j  = se[j] if not np.isnan(se[j]) else 0.0
                    ax.errorbar(z_c[j], mu[j], yerr=se_j,
                                fmt='o', color=color, ecolor=color,
                                capsize=3, markersize=6, lw=2.0, zorder=4)

                nan_idx = np.where(~valid)[0]
                if len(nan_idx):
                    ax.scatter(z_c[nan_idx],
                               np.full(len(nan_idx), float(mu[valid].min()) if valid.any() else 0),
                               marker='x', color='grey', s=50, zorder=5,
                               label=f'{len(nan_idx)} skipped bins')

                # Colorbar
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=p_norm)
                sm.set_array([])
                plt.colorbar(sm, ax=ax,
                             label='p_cav  (void fraction of trial points)',
                             fraction=0.03, pad=0.02)

                if valid.sum() >= 2:
                    mu_v  = mu[valid]
                    dmu   = float(mu_v.max() - mu_v.min())
                    ax.text(0.02, 0.05,
                            f'Valid bins: {valid.sum()}/{len(mu)}'
                            f'   mean μ_ex = {float(mu_v.mean()):.4f} ε'
                            f'   Δμ_ex = {dmu:.4f} ε',
                            transform=ax.transAxes, fontsize=9, va='bottom',
                            bbox=dict(boxstyle='round,pad=0.3',
                                      facecolor='white', alpha=0.85))
                    if dmu > 0.05:
                        ax.text(0.02, 0.92,
                                f'⚠ Δμ_ex = {dmu:.4f} ε  (> 0.05 ε) — '
                                'chemical potential gradient present',
                                transform=ax.transAxes, fontsize=9,
                                color='red', va='top',
                                bbox=dict(boxstyle='round,pad=0.3',
                                          facecolor='white', alpha=0.85))

        if cav_frames:
            ax.legend(fontsize=7, ncol=min(len(cav_frames), 5),
                      loc='upper right', framealpha=0.7)

        ax.set_ylabel('μ_ex  (ε)')
        ax.set_title(
            'Cavity-biased Widom  μ_ex(z)\n'
            'Light curves = per-frame snapshots;  thick markers = time-avg ± stderr'
            '  (colour = p_cav).\n'
            'Flat profile → chemical equilibrium;  gradient → driving force for solvent flow.',
            fontsize=9)

    # ── Panel 3: Pore pressure p_p(z)/P_ext evolution ────────────────────────
    if show_pore_pres:
        ax = axes[panel_idx]; panel_idx += 1
        ax.grid(alpha=0.3)

        if cav_frames and rho_res is not None and press_ext > 0:
            steps     = sorted(cav_frames.keys())
            n_cf      = len(steps)
            colors_pp = plt.cm.viridis(np.linspace(0, 1, max(n_cf, 1)))

            for i, step in enumerate(steps):
                d  = cav_frames[step]
                z  = d['z']
                mu = d['mu_ex']
                ok = ~np.isnan(mu)
                # p_p(z) = μ_ex(z) / V̄_s = μ_ex(z) × ρ_s,reservoir
                # normalized by P_ext so reservoir → 1
                pp = mu[ok] * rho_res / press_ext
                ax.plot(z[ok], pp, color=colors_pp[i], lw=1.2,
                        alpha=0.7, marker='o', markersize=3,
                        label=f'step {step:,}')

            # Time-averaged pore pressure from summary file
            if has_cavity_sum:
                try:
                    cav = np.genfromtxt(cavity_sum_file, comments='#')
                    if cav.ndim == 2 and cav.shape[0] >= 2 and cav.shape[1] >= 4:
                        z_c   = cav[:, 0]
                        mu_m  = cav[:, 3]
                        valid = ~np.isnan(mu_m)
                        pp_m  = mu_m[valid] * rho_res / press_ext
                        ax.plot(z_c[valid], pp_m, color='k', lw=2.0,
                                linestyle='--', marker='D', markersize=4,
                                zorder=5, label='time-avg')
                except Exception:
                    pass

            ax.axhline(1.0, color='gray', lw=1.0, ls=':', label='P_ext')
            ax.set_ylabel(f'p_p(z) / P_ext')
            ax.set_title(
                f'Pore pressure profile  p_p(z) / P_ext\n'
                f'p_p = μ_ex × ρ_s,reservoir = μ_ex / V̄_s;   '
                f'P_ext = {press_ext:.3f} ε/σ³  (barostat).\n'
                f'Dotted line = 1 (reservoir reference).  '
                f'Deviation → osmotic pressure imbalance.',
                fontsize=9)
            ax.legend(fontsize=7, ncol=min(len(cav_frames) + 2, 6),
                      loc='upper right', framealpha=0.7)

            if has_cavity_sum:
                try:
                    pp_all = mu_m[valid] * rho_res / press_ext
                    dp     = float(pp_all.max() - pp_all.min())
                    ax.text(0.02, 0.05,
                            f'Δ(p_p/P_ext) = {dp:.4f}  '
                            f'(= {dp * press_ext:.4f} ε/σ³)',
                            transform=ax.transAxes, fontsize=9, va='bottom',
                            bbox=dict(boxstyle='round,pad=0.3',
                                      facecolor='white', alpha=0.85))
                except Exception:
                    pass
        elif rho_res is None:
            ax.text(0.5, 0.5,
                    'Reservoir ρ_s not available — cannot compute V̄_s.\n'
                    'Check solvent_density_z file.',
                    transform=ax.transAxes, ha='center', va='center', fontsize=9)
        else:
            ax.text(0.5, 0.5, 'No per-frame cavity data available.',
                    transform=ax.transAxes, ha='center', va='center', fontsize=9)

    # Shared x-label on bottom panel
    axes[-1].set_xlabel('z  (σ)')
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Chempot diagnostics → {output}")


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

    # Flow/compression diagnostics (auto-triggered when piston force files are present).
    # Mode (compression vs permeation) is detected inside plot_flow_diagnostics from
    # which output files exist — no flag needed.
    piston_dir = os.path.join(args.folder, 'output_files', 'piston_data')
    force_files_present = (
        os.path.exists(os.path.join(piston_dir, f'piston_force_{run_id}.dat')) or
        os.path.exists(os.path.join(piston_dir, f'piston_force_pressure_{run_id}.dat'))
    )
    if force_files_present:
        plot_flow_diagnostics(args.folder, run_id,
                              os.path.join(out_dir, f'{run_id}_flow_diagnostics.png'))

    # Chemical potential diagnostics (auto-triggered when solvent_density_z or
    # cavity Widom files are present in output_files/chemical_potential/).
    # Standard fix-widom output (mu_z_*.dat) is intentionally excluded —
    # it fails for WCA fluids at liquid density.  Cavity Widom only.
    chem_dir = os.path.join(args.folder, 'output_files', 'chemical_potential')
    chempot_files_present = (
        os.path.exists(os.path.join(chem_dir, f'solvent_density_z_{run_id}.dat')) or
        os.path.exists(os.path.join(chem_dir, f'mu_z_cavity_summary_{run_id}.dat'))
    )
    if chempot_files_present:
        plot_chempot_diagnostics(args.folder, run_id,
                                 os.path.join(out_dir,
                                              f'{run_id}_chempot_diagnostics.png'),
                                 thermo_data=data)
