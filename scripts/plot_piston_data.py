#!/usr/bin/env python3
"""
Plot piston position and velocity versus timesteps.
Reads piston position and velocity from piston_data folder.

Multi-piston aware (2026-09-16): the two-piston decks (triaxial_*_two_pist)
write ONE COLUMN SET PER PISTON, e.g.
    # step z_feed z_perm            (permeation)
    # step z_dry z_feed z_perm      (compression)
and a per-piston pressure log piston_pressure_<stem>.dat
    # step P_dry_meas P_feed_meas P_feed_app P_perm_meas P_perm_app
plus permeation_data/permeation_<stem>.dat (permeation: Q_perm vs time;
compression: solvent expelled dV_total).  When more than one piston column is
present, one panel row is drawn per quantity with one line per piston, and the
pressure / permeation panels are appended.  ONE-PISTON FILES PLOT EXACTLY AS
BEFORE (same panels, labels and output file name).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os


def read_piston_file(filepath):
    """Read LAMMPS fix print output file with format: timestep value
    (the FIRST value column; multi-column files are handled by read_table)."""
    timesteps = []
    values = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    timesteps.append(int(float(parts[0])))
                    values.append(float(parts[1]))
                except ValueError:
                    continue

    return np.array(timesteps), np.array(values)


def read_table(filepath):
    """Multi-column fix print / fix ave/time file -> (names, array).
    names come from a '# col col ...' header whose token count matches the data
    columns (fix print `title`, fix ave/time `title2`); None when absent."""
    names, rows = None, []
    with open(filepath) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                toks = s.lstrip('#').split()
                # keep the last header line whose leading tokens look like names
                if toks and not toks[0].replace('.', '').replace('-', '').isdigit():
                    names = toks
                continue
            try:
                rows.append([float(x) for x in s.split()])
            except ValueError:
                continue
    if not rows:
        return None, np.empty((0, 0))
    ncol = min(len(r) for r in rows)
    arr = np.array([r[:ncol] for r in rows])
    if names is not None:
        # fix ave/time headers may carry a trailing comment after the names
        names = names[:ncol] if len(names) >= ncol else None
    return names, arr


def piston_labels(names, ncol, prefix):
    """Piston labels from header names like z_feed / vz_perm / z_dry."""
    n = ncol - 1
    if names is not None and len(names) == ncol:
        return [nm[len(prefix):] if nm.startswith(prefix) else nm for nm in names[1:]]
    return ['piston'] if n == 1 else [f'piston{i + 1}' for i in range(n)]


def _plot_one_piston(folder, dataname, oldsteps, pos_file, vel_file):
    """The original single-piston figure, unchanged."""
    timesteps_pos, positions = read_piston_file(pos_file)
    timesteps_vel, velocities = read_piston_file(vel_file)

    if len(timesteps_pos) == 0:
        print(f"Error: No data found in {pos_file}")
        return False

    print(f"Read {len(timesteps_pos)} position points, {len(timesteps_vel)} velocity points")

    output_plot_dir = os.path.join(folder, 'output_plots')
    os.makedirs(output_plot_dir, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    if oldsteps > 0:
        fig.suptitle(f'{dataname}\n(continuing from {oldsteps} steps)',
                     fontsize=12, fontweight='bold')
    else:
        fig.suptitle(f'{dataname}', fontsize=12, fontweight='bold')

    color_pos = plt.cm.viridis(0.3)
    color_vel = plt.cm.viridis(0.7)

    axes[0].plot(timesteps_pos, positions, linewidth=1.5, color=color_pos, alpha=0.9)
    axes[0].set_ylabel('Piston Position (z)', fontsize=11)
    axes[0].grid(alpha=0.3)
    axes[0].set_title('Position vs Time', fontweight='bold')
    if len(positions) > 0:
        axes[0].axhline(y=positions[0], color='gray', linestyle='--',
                        alpha=0.5, label=f'Initial: {positions[0]:.2f}')
        axes[0].legend(loc='best', fontsize=9)

    axes[1].plot(timesteps_vel, velocities, linewidth=1.5, color=color_vel, alpha=0.9)
    axes[1].set_ylabel('Piston Velocity (vz)', fontsize=11)
    axes[1].set_xlabel('Timestep', fontsize=11)
    axes[1].grid(alpha=0.3)
    axes[1].set_title('Velocity vs Time', fontweight='bold')
    axes[1].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    if len(velocities) > 10:
        steady_start = len(velocities) // 2
        mean_vel = np.mean(velocities[steady_start:])
        std_vel = np.std(velocities[steady_start:])
        axes[1].axhline(y=mean_vel, color='red', linestyle='-', alpha=0.7,
                        label=f'Steady-state mean: {mean_vel:.2e} ± {std_vel:.2e}')
        axes[1].legend(loc='best', fontsize=9)

    plt.tight_layout()
    output_plot_file = os.path.join(output_plot_dir, f'{dataname}_piston.png')
    plt.savefig(output_plot_file, dpi=150)
    print(f"Piston plot saved to {output_plot_file}")
    plt.close()

    print("\n=== Piston Summary ===")
    print(f"Initial position: {positions[0]:.4f}")
    print(f"Final position:   {positions[-1]:.4f}")
    print(f"Total displacement: {positions[-1] - positions[0]:.4f}")
    if len(velocities) > 10:
        steady_start = len(velocities) // 2
        print(f"Steady-state velocity: {np.mean(velocities[steady_start:]):.2e} ± {np.std(velocities[steady_start:]):.2e}")
    return True


def _plot_multi_piston(folder, dataname, oldsteps, pos_file, vel_file, pos_names, pos, vel_names, vel):
    """Two-piston layout: position + velocity rows (one line per piston), then the
    per-piston pressure panel (measured solid vs applied dashed) and the
    permeation panel (Q_perm and N_permeate, or solvent expelled) when present."""
    piston_dir = os.path.join(folder, 'output_files', 'piston_data')
    perm_dir = os.path.join(folder, 'output_files', 'permeation_data')
    labels = piston_labels(pos_names, pos.shape[1], 'z_')
    vlabels = piston_labels(vel_names, vel.shape[1], 'vz_') if vel.size else labels
    print(f"Read {len(pos)} position rows for {len(labels)} pistons ({', '.join(labels)}), {len(vel)} velocity rows")

    pres_file = os.path.join(piston_dir, f'piston_pressure_{dataname}.dat')
    perm_file = os.path.join(perm_dir, f'permeation_{dataname}.dat')
    pres_names, pres = read_table(pres_file) if os.path.exists(pres_file) else (None, np.empty((0, 0)))
    perm_names, perm = read_table(perm_file) if os.path.exists(perm_file) else (None, np.empty((0, 0)))
    has_pres = pres.size > 0 and pres_names is not None
    has_perm = perm.size > 0 and perm_names is not None

    n_pan = 2 + int(has_pres) + int(has_perm)
    fig, axes = plt.subplots(n_pan, 1, figsize=(10, 3.6 * n_pan), sharex=True)
    axes = np.atleast_1d(axes)
    title = f'{dataname}' + (f'\n(continuing from {oldsteps} steps)' if oldsteps > 0 else '')
    fig.suptitle(title + '  —  two-piston', fontsize=12, fontweight='bold')
    cmap = plt.cm.viridis
    cols = {lab: cmap(0.15 + 0.7 * i / max(len(labels) - 1, 1)) for i, lab in enumerate(labels)}

    ax = axes[0]
    for i, lab in enumerate(labels):
        ax.plot(pos[:, 0], pos[:, i + 1], lw=1.5, color=cols.get(lab, 'k'), alpha=0.9,
                label=f'{lab}  (initial {pos[0, i + 1]:.2f}, final {pos[-1, i + 1]:.2f})')
    ax.set_ylabel('Piston Position (z)', fontsize=11)
    ax.set_title('Position vs Time (one line per piston)', fontweight='bold')
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=9)

    ax = axes[1]
    if vel.size:
        for i, lab in enumerate(vlabels):
            v = vel[:, i + 1]
            lbl = lab
            if len(v) > 10:
                h = len(v) // 2
                lbl += f'  (steady mean {np.mean(v[h:]):.2e} ± {np.std(v[h:]):.2e})'
            ax.plot(vel[:, 0], v, lw=1.2, color=cols.get(lab, 'k'), alpha=0.85, label=lbl)
    ax.axhline(0, color='gray', ls='--', alpha=0.5)
    ax.set_ylabel('Piston Velocity (vz)', fontsize=11)
    ax.set_title('Velocity vs Time' + ('' if vel.size else '  (no velocity file)'), fontweight='bold')
    ax.grid(alpha=0.3)
    if vel.size:
        ax.legend(loc='best', fontsize=9)

    k = 2
    if has_pres:
        ax = axes[k]
        k += 1
        t = pres[:, 0]
        for j, nm in enumerate(pres_names[1:], 1):
            if nm.endswith('_meas'):
                base = nm[len('P_'):-len('_meas')]
                ax.plot(t, pres[:, j], lw=1.5, color=cols.get(base, 'k'), alpha=0.9, label=f'P_{base} measured')
            elif nm.endswith('_app'):
                base = nm[len('P_'):-len('_app')]
                ax.plot(t, pres[:, j], lw=1.2, ls='--', color=cols.get(base, 'k'), alpha=0.8, label=f'P_{base} applied')
        ax.set_ylabel('Pressure  (ε/σ³)', fontsize=11)
        ax.set_title('Piston pressure: F_fluid/(lx·ly) measured (solid) vs applied (dashed)', fontweight='bold')
        ax.grid(alpha=0.3)
        ax.legend(loc='best', fontsize=9, ncol=2)

    if has_perm:
        ax = axes[k]
        t = perm[:, 0]
        cols_lc = [c.lower() for c in perm_names]
        qi = next((j for j, c in enumerate(cols_lc) if c.startswith('q_perm')), None)
        ni = next((j for j, c in enumerate(cols_lc) if c.startswith('n_permeate')), None)
        dvi = next((j for j, c in enumerate(cols_lc) if c == 'dv_total'), None)
        if qi is not None:
            q = perm[:, qi]
            ax.plot(t, q, lw=1.5, color='teal', alpha=0.9, label='Q_perm = A·dz_perm/dt (block-averaged)')
            if len(q) > 10:
                h = len(q) // 2
                ax.axhline(np.mean(q[h:]), color='red', ls='-', alpha=0.7,
                           label=f'steady mean: {np.mean(q[h:]):.3e} ± {np.std(q[h:]):.2e}')
            ax.set_ylabel('Q_perm  (σ³/τ)', fontsize=11)
            ax.set_title('Permeate flux from the permeate-piston displacement', fontweight='bold')
            if ni is not None:
                ax2 = ax.twinx()
                ax2.plot(t, perm[:, ni] - perm[0, ni], lw=1.0, color='0.4', alpha=0.7, label='ΔN_permeate (bead count)')
                ax2.set_ylabel('ΔN_permeate', fontsize=10)
                ax2.legend(loc='lower right', fontsize=8)
        elif dvi is not None:
            ax.plot(t, perm[:, dvi], lw=1.5, color='teal', alpha=0.9, label='dV_total = solvent expelled (A·Δz of the wet pistons)')
            ax.set_ylabel('expelled volume  (σ³)', fontsize=11)
            ax.set_title('Solvent expelled into the reservoirs', fontweight='bold')
        ax.grid(alpha=0.3)
        ax.legend(loc='best', fontsize=9)

    axes[-1].set_xlabel('Timestep', fontsize=11)
    plt.tight_layout()
    output_plot_dir = os.path.join(folder, 'output_plots')
    os.makedirs(output_plot_dir, exist_ok=True)
    output_plot_file = os.path.join(output_plot_dir, f'{dataname}_piston.png')
    plt.savefig(output_plot_file, dpi=150)
    print(f"Piston plot saved to {output_plot_file}")
    plt.close()

    print("\n=== Piston Summary (two-piston) ===")
    for i, lab in enumerate(labels):
        line = f"{lab:>10s}: z {pos[0, i + 1]:.4f} -> {pos[-1, i + 1]:.4f}  (displacement {pos[-1, i + 1] - pos[0, i + 1]:+.4f})"
        if vel.size and len(vel) > 10 and i + 1 < vel.shape[1]:
            h = len(vel) // 2
            line += f"   steady v {np.mean(vel[h:, i + 1]):.2e} ± {np.std(vel[h:, i + 1]):.2e}"
        print(line)
    if has_pres:
        h = len(pres) // 2
        for j, nm in enumerate(pres_names[1:], 1):
            print(f"{nm:>14s}: last-half mean {np.mean(pres[h:, j]):.4f} ± {np.std(pres[h:, j]):.4f}")
    return True


def plot_piston_data(folder, dataname, oldsteps=0):
    """Plot piston position and velocity versus timesteps."""

    piston_data_dir = os.path.join(folder, 'output_files', 'piston_data')
    pos_file = os.path.join(piston_data_dir, f'piston_position_{dataname}.dat')
    vel_file = os.path.join(piston_data_dir, f'piston_velocity_{dataname}.dat')

    # Check files exist. Some geometries (e.g. slab_with_support, compress_slab)
    # never write piston position/velocity files, so a missing file is expected,
    # not an error — skip cleanly instead of printing a scary "Error:" line.
    # Two-piston compression writes its whole-run positions as piston_position_run_*.
    if not os.path.exists(pos_file):
        alt = os.path.join(piston_data_dir, f'piston_position_run_{dataname}.dat')
        if os.path.exists(alt):
            pos_file = alt
        else:
            print(f"No piston position file ({os.path.basename(pos_file)}) — skipping piston plots.")
            return True

    pos_names, pos = read_table(pos_file)
    if pos.size == 0:
        print(f"Error: No data found in {pos_file}")
        return False
    n_pist = pos.shape[1] - 1
    if n_pist <= 1:
        if not os.path.exists(vel_file):
            print(f"No piston velocity file ({os.path.basename(vel_file)}) — skipping piston plots.")
            return True
        return _plot_one_piston(folder, dataname, oldsteps, pos_file, vel_file)
    vel_names, vel = read_table(vel_file) if os.path.exists(vel_file) else (None, np.empty((0, 0)))
    return _plot_multi_piston(folder, dataname, oldsteps, pos_file, vel_file, pos_names, pos, vel_names, vel)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python plot_piston_data.py <folder> <dataname> [oldsteps]")
        print("Example: python plot_piston_data.py . mydata_1.5_1.4_100000 0")
        sys.exit(1)

    folder = sys.argv[1]
    dataname = sys.argv[2]
    oldsteps = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    success = plot_piston_data(folder, dataname, oldsteps)
    if not success:
        sys.exit(1)
