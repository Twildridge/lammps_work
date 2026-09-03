#!/usr/bin/env python3
"""Cluster-side quick-look for the triaxial_compression stress-strain sweep.

Reads the per-level output files written by triaxial_compression.lmp (each level
tagged _c<level>) and CONSOLIDATES every quantity into a single figure per
quantity, with one curve per compression level -- instead of one PNG per level.
This is the compression analog of plot_shear_strain_sweep.py.

Produces, in <folder>/output_plots/:

  1. <stem>_piston_sweep.png
       Piston position, velocity, and force vs timestep -- all levels overlaid
       (one color per level).  Replaces the old per-level *_c<lvl>_piston.png.

  2. <stem>_stress_profiles_sweep.png
       Final-snapshot z-binned total stress profiles (isotropic-z, sigma_zz,
       sigma_xx, sigma_yy) -- one curve per level.  A compact consolidated
       stand-in for the old per-level 6x5 *_c<lvl>_stress.png grids.

The notebooks (triaxial_compression.ipynb) remain the source of truth for the
detailed analysis; this is only an automatic first look so a run's sweep can be
eyeballed without opening the notebook.  Only numpy + matplotlib are required.

Usage:
  python plot_compression_strain_sweep.py <folder> <stem> "<lvl1 lvl2 ...>" [oldsteps]
    folder : run dir (".")   stem : DATANAME_INTERACTION_TOTSTEPS
    levels : space-separated list, exactly the literals used as _c tags
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def read_xy(path):
    """fix print file with 'step value' rows -> (steps[N], values[N]). None if absent/empty."""
    if not os.path.exists(path):
        return None, None
    s, v = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    s.append(int(float(parts[0])))
                    v.append(float(parts[1]))
                except ValueError:
                    continue
    if not s:
        return None, None
    return np.array(s), np.array(v)


def read_ave_time_last(path):
    """fix ave/time (mode vector) file -> (rows[M], values[M]) for the LAST snapshot.
    Per-timestep header has 2 fields: timestep nrows.  None if absent/empty."""
    if not os.path.exists(path):
        return None, None
    with open(path) as f:
        lines = [l for l in f if l.strip() and not l.startswith('#')]
    last_vals = None
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) == 2:
            try:
                nrows = int(parts[1])
            except ValueError:
                i += 1
                continue
            vals = []
            for j in range(1, nrows + 1):
                if i + j < len(lines):
                    p = lines[i + j].split()
                    if len(p) >= 2:
                        vals.append(float(p[1]))
            if vals:
                last_vals = np.array(vals)
            i += nrows + 1
        else:
            i += 1
    if last_vals is None:
        return None, None
    rows = np.arange(1, len(last_vals) + 1)
    return rows, last_vals


def _level_colors(levels):
    vals = []
    for s in levels:
        try:
            vals.append(float(s))
        except ValueError:
            vals.append(np.nan)
    finite = [v for v in vals if np.isfinite(v)]
    if len(set(finite)) > 1:
        norm = Normalize(vmin=min(finite), vmax=max(finite))
    else:
        norm = Normalize(0, 1)
    cmap = plt.cm.viridis
    return [cmap(norm(v)) if np.isfinite(v) else cmap(0.5) for v in vals]



def find_file(path):
    """Exact path if it exists; else the same name with the <steps> tag wild-carded
    (production files are tagged with each level's auto-sized hold length, so the
    stem's TOTSTEPS need not match). Returns the newest match, or the original path."""
    if os.path.exists(path):
        return path
    import glob, re
    pat = re.sub(r'_(\d+)(_c[\d.]+\.dat)$', r'_*\2', path)
    hits = sorted(glob.glob(pat), key=os.path.getmtime)
    return hits[-1] if hits else path


def plot_piston_sweep(folder, stem, levels, colors, oldsteps=0):
    pd = os.path.join(folder, 'output_files', 'piston_data')
    plot_dir = os.path.join(folder, 'output_plots')
    os.makedirs(plot_dir, exist_ok=True)

    panels = [
        ('piston_position', 'Piston position (z)', 'Position vs time'),
        ('piston_velocity', 'Piston velocity (vz)', 'Velocity vs time'),
        ('piston_force',    'Piston force (fz)',    'Force vs time'),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    any_data = False
    for ax, (prefix, ylabel, title) in zip(axes, panels):
        for s, c in zip(levels, colors):
            steps, vals = read_xy(find_file(os.path.join(pd, f"{prefix}_{stem}_c{s}.dat")))
            if steps is None:
                continue
            ax.plot(steps, vals, lw=1.5, color=c, alpha=0.9, label=fr'$\varepsilon={s}$')
            any_data = True
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontweight='bold')
        ax.grid(alpha=0.3)
        if prefix != 'piston_position':
            ax.axhline(0, color='gray', ls='--', alpha=0.5)
    axes[-1].set_xlabel('Timestep', fontsize=11)
    if not any_data:
        plt.close(fig)
        print("  No piston files found for any level — skipping piston sweep plot.")
        return
    handles, lbls = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(loc='best', fontsize=9, title='cumulative level')
    suptitle = stem + (f'  (continuing from {oldsteps} steps)' if oldsteps > 0 else '')
    fig.suptitle(suptitle, fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = os.path.join(plot_dir, f"{stem}_piston_sweep.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_stress_profiles_sweep(folder, stem, levels, colors):
    sd = os.path.join(folder, 'output_files', 'stress_data')
    plot_dir = os.path.join(folder, 'output_plots')
    os.makedirs(plot_dir, exist_ok=True)

    # (title, polymer-file-prefix, solvent-file-prefix, extra-prefixes-for-total)
    specs = [
        (r'Isotropic $p(z)$', 'stress_z_polymer', 'stress_z_solvent',
         ['stress_z_piston', 'stress_z_support']),
        (r'$\sigma_{zz}(z)$', 'sigmazz_polymer', 'sigmazz_solvent', []),
        (r'$\sigma_{xx}(z)$', 'sigmaxx_polymer', 'sigmaxx_solvent', []),
        (r'$\sigma_{yy}(z)$', 'sigmayy_polymer', 'sigmayy_solvent', []),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    axes = axes.ravel()
    any_data = False
    for ax, (title, poly_p, solv_p, extras) in zip(axes, specs):
        for s, c in zip(levels, colors):
            comp_prefixes = [poly_p, solv_p] + extras
            total = None
            zfrac = None
            for pref in comp_prefixes:
                rows, vals = read_ave_time_last(find_file(os.path.join(sd, f"{pref}_{stem}_c{s}.dat")))
                if vals is None:
                    continue
                if total is None:
                    total = np.zeros_like(vals)
                    zfrac = (rows - 0.5) / len(rows)
                n = min(len(total), len(vals))
                total[:n] += vals[:n]
            if total is None:
                continue
            ax.plot(zfrac, total, '-o', ms=3, lw=1.5, color=c, alpha=0.9,
                    label=fr'$\varepsilon={s}$')
            any_data = True
        ax.axhline(0, color='k', ls='--', lw=0.8, alpha=0.4)
        ax.set_title(title, fontweight='bold')
        ax.set_xlabel('z (bin fraction)', fontsize=11)
        ax.set_ylabel('total stress (LJ)', fontsize=11)
        ax.grid(alpha=0.3)
    if not any_data:
        plt.close(fig)
        print("  No stress-profile files found for any level — skipping stress sweep plot.")
        return
    handles, lbls = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(loc='best', fontsize=10, title='cumulative level')
    fig.suptitle(f'{stem}\nFinal-snapshot total stress profiles (all levels)',
                 fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(plot_dir, f"{stem}_stress_profiles_sweep.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print('Usage: python plot_compression_strain_sweep.py <folder> <stem> "<lvl1 lvl2 ...>" [oldsteps]')
        sys.exit(1)
    folder = sys.argv[1]
    stem = sys.argv[2]
    levels = sys.argv[3].split()
    oldsteps = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    print(f"Compression sweep post-processing: stem={stem}  levels={levels}")
    colors = _level_colors(levels)
    try:
        plot_piston_sweep(folder, stem, levels, colors, oldsteps)
        plot_stress_profiles_sweep(folder, stem, levels, colors)
    except Exception as e:
        # Never let a quick-look plot fail the whole post-processing chain.
        print(f"  plot_compression_strain_sweep.py warning: {e}")
