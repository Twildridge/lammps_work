#!/usr/bin/env python3
"""
plot_eos.py  —  Plot P* vs ρ* from a pure_solvent_psweep.lmp or pure_polymer.lmp EOS run.

Usage:
    python plot_eos.py <folder> <dataname> <interaction>

Example:
    python plot_eos.py . pure_solvent_1000 1p0
    python plot_eos.py . pure_polymer_5beads_8x8x8 1p0

Expects files at:
    <folder>/eos_<dataname>_<interaction>.dat
    <folder>/log.lammps                                                        (for volume panel)
    <folder>/output_files/volume_data/rho_trace_<dataname>_<interaction>.dat  (optional)

Output:
    <folder>/output_plots/convergence_plots/<dataname>_<interaction>_eos.png
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os


# ── Fonts (match notebook style) ─────────────────────────────────────────────
plt.rcParams.update({
    'font.family':      'serif',
    'font.serif':       ['CMU Serif', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'axes.labelsize':   13,
    'xtick.labelsize':  11,
    'ytick.labelsize':  11,
    'legend.fontsize':  11,
})


def read_eos(filepath):
    """Read P_target, rho_mean from EOS file (skips # comment lines).
    Returns arrays in file order (= P-sweep order, NOT sorted by P)."""
    P, rho = [], []
    with open(filepath) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 2:
                    P.append(float(parts[0]))
                    rho.append(float(parts[1]))
    return np.array(P), np.array(rho)


def read_rho_trace(filepath):
    """
    Read density trace file.  Columns: label  step  rho
    Returns dict: label -> (steps, rhos)
    """
    traces = {}
    with open(filepath) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                if len(parts) == 3:
                    label, step, rho = parts[0], float(parts[1]), float(parts[2])
                    if label not in traces:
                        traces[label] = ([], [])
                    traces[label][0].append(step)
                    traces[label][1].append(rho)
    return {k: (np.array(v[0]), np.array(v[1])) for k, v in traces.items()}


def read_log_volumes(filepath, is_polymer=False):
    """
    Parse log.lammps and return the mean box volume for each NPT state point.

    Strategy:
      - Find thermo header lines (contain 'Step' and 'Volume').
      - Read subsequent numeric rows until a non-numeric line is reached.
      - Each uninterrupted block of thermo data = one run segment.
      - State points are delineated by 'reset_timestep 0' in the log, which
        causes the Step counter to restart at 0.  Eq and prod are separate
        runs within the same state point — prod continues from where eq left
        off, so its first step > 0.  We group consecutive segments: a new
        state point starts whenever a segment's first step == 0.
      - Mean volume is computed over all segments in the group (eq + prod).

    Returns: list of (mean_vol, std_vol) in log/sweep order, one per state
             point.  Returns empty list if log file not found or has no
             Volume column.
    """
    if not os.path.exists(filepath):
        return []

    # ── 1. Parse all thermo blocks ──────────────────────────────────────────
    segments = []
    vol_col  = None
    in_block = False
    cur_steps, cur_vols = [], []

    with open(filepath) as f:
        for line in f:
            stripped = line.strip()

            # Detect thermo header line
            if stripped.lower().startswith('step'):
                tokens = stripped.split()
                vol_col = None
                for ci, tok in enumerate(tokens):
                    if tok.lower() in ('volume', 'vol'):
                        vol_col = ci
                        break
                if vol_col is not None:
                    if cur_steps:
                        segments.append({'steps': cur_steps, 'vols': cur_vols})
                        cur_steps, cur_vols = [], []
                    in_block = True
                continue

            if in_block and vol_col is not None:
                tokens = stripped.split()
                if not tokens:
                    continue
                try:
                    vals = [float(t) for t in tokens]
                    if len(vals) > vol_col:
                        cur_steps.append(int(vals[0]))
                        cur_vols.append(vals[vol_col])
                except ValueError:
                    if cur_steps:
                        segments.append({'steps': cur_steps, 'vols': cur_vols})
                        cur_steps, cur_vols = [], []
                    in_block = False

    if cur_steps:
        segments.append({'steps': cur_steps, 'vols': cur_vols})

    if not segments:
        return []

    # ── 2. Group segments into state points by step-counter resets ──────────
    # A new state point begins when its eq run starts with reset_timestep 0.
    # The first step-0 group is always the gentle NVT warmup and is discarded
    # via result[1:] at the end.  For polymer, vols[1:] drops the step-0
    # carry-over snapshot from each NPT eq run (not the gentle NVT, since it
    # is discarded anyway — slicing it would make it empty and cause it to be
    # silently skipped, throwing off the result[1:] index).
    state_points = []
    group_vols   = []
    step0_count  = 0    # how many step-0 segments seen so far

    for seg in segments:
        if not seg['steps']:
            continue
        if seg['steps'][0] == 0:
            if group_vols:
                state_points.append(group_vols)
            step0_count += 1
            vols = list(seg['vols'])
            # Apply polymer carry-over fix only to NPT state points (step0_count > 1)
            group_vols = vols[1:] if (is_polymer and step0_count > 1) else vols
        else:
            group_vols.extend(seg['vols'])

    if group_vols:
        state_points.append(group_vols)

    # ── 3. Compute mean ± std per state point ───────────────────────────────
    result = []
    for vols in state_points:
        arr = np.array(vols)
        result.append((arr.mean(), arr.std()))

    return result[1:]   # drop the initial gentle NVT warmup block


def plot_eos(foldername, dataname, interaction):

    eos_file   = os.path.join(foldername, f'eos_{dataname}_{interaction}.dat')
    trace_file = os.path.join(foldername, 'output_files/volume_data',
                              f'rho_trace_{dataname}_{interaction}.dat')
    log_file   = os.path.join(foldername, 'log.lammps')

    if not os.path.exists(eos_file):
        print(f"Error: EOS file not found: {eos_file}")
        sys.exit(1)

    # Detect whether this is a polymer run (suppress VLE annotations)
    is_polymer = 'polymer' in dataname.lower()

    # Read EOS in sweep order (before sorting) so indices align with log order
    P_sweep, rho_sweep = read_eos(eos_file)

    # Read volumes from log in the same sweep order
    vol_data = read_log_volumes(log_file, is_polymer=is_polymer)
    has_vol  = len(vol_data) == len(P_sweep)
    if vol_data and not has_vol:
        print(f"Warning: log.lammps has {len(vol_data)} state points "
              f"but EOS file has {len(P_sweep)} — skipping volume panel.")

    # Sort everything by P ascending for plots
    order     = np.argsort(P_sweep)
    P         = P_sweep[order]
    rho       = rho_sweep[order]
    if has_vol:
        vol_means = np.array([vol_data[i][0] for i in order])
        vol_stds  = np.array([vol_data[i][1] for i in order])

    # Detect VLE transition (used only for solvent)
    drho     = np.diff(rho)
    idx_jump = np.argmin(drho)
    P_trans  = 0.5 * (P[idx_jump] + P[idx_jump + 1])

    # ── Bulk modulus (always computed from EOS data) ───────────────────────
    drho_dP = np.gradient(rho, P)
    beta_T  = drho_dP / rho
    K_T     = 1.0 / beta_T

    # ── Panel count ───────────────────────────────────────────────────────
    has_trace = os.path.exists(trace_file)
    nrows     = 2 + int(has_vol) + int(has_trace)   # +1 for K_T panel
    fig, axes = plt.subplots(nrows, 1, figsize=(8, 5 * nrows))
    if nrows == 1:
        axes = [axes]

    title_str = f'{dataname}_{interaction}'
    fig.suptitle(title_str, fontsize=13, fontweight='bold')

    ax_idx = 0

    # ── Panel 1: P* vs rho* ──────────────────────────────────────────────
    ax = axes[ax_idx]; ax_idx += 1
    ax.plot(rho, P, 'o-', color='steelblue', linewidth=1.8,
            markersize=5, label='NPT equilibrated')

    if not is_polymer:
        ax.axvline(x=rho[idx_jump],     color='r', linestyle='--', alpha=0.6)
        ax.axvline(x=rho[idx_jump + 1], color='r', linestyle='--', alpha=0.6)
        ax.axhspan(P[idx_jump], P[idx_jump + 1], alpha=0.12, color='red',
                   label=f'VLE region  P* ~ {P_trans:.3f}')

    ax.set_xlabel(r'Number density $\rho^*$ ($\sigma^{-3}$)')
    ax.set_ylabel(r'Pressure $P^*$')
    title_suffix = 'Polymer Melt' if is_polymer else 'Solvent'
    ax.set_title(rf'LJ Equation of State — {title_suffix}  ($T^* = 1.0$,  $r_c = 2.5\sigma$)')
    ax.grid(alpha=0.3)
    ax.legend()

    # Terminal summary
    print(f"\n{'P*':>10}  {'rho*':>12}")
    print('-' * 25)
    for p_val, r_val in zip(P, rho):
        flag = ' <- VLE' if (not is_polymer and abs(p_val - P_trans) < 0.01) else ''
        print(f"{p_val:>10.4f}  {r_val:>12.5f}{flag}")
    if not is_polymer:
        print(f"\nDetected VLE transition:  P* ~ {P_trans:.4f}")

    # ── Panel 2: P* vs box volume ─────────────────────────────────────────
    if has_vol:
        ax2 = axes[ax_idx]; ax_idx += 1
        ax2.errorbar(vol_means, P, xerr=vol_stds,
                     fmt='o-', color='darkorange', linewidth=1.8,
                     markersize=5, capsize=3, elinewidth=1.0,
                     label='Mean ± std (eq + prod)')

        if not is_polymer:
            ax2.axhspan(P[idx_jump], P[idx_jump + 1], alpha=0.12, color='red',
                        label=f'VLE region  P* ~ {P_trans:.3f}')

        ax2.set_xlabel(r'Box volume $V^*$ ($\sigma^{3}$)')
        ax2.set_ylabel(r'Pressure $P^*$')
        ax2.set_title(r'$P^*$ vs Box Volume  (averaged over eq + prod per state point)')
        ax2.grid(alpha=0.3)
        ax2.legend()

        print(f"\n{'P*':>10}  {'V_mean':>14}  {'V_std':>10}")
        print('-' * 38)
        for p_val, v_mean, v_std in zip(P, vol_means, vol_stds):
            print(f"{p_val:>10.4f}  {v_mean:>14.2f}  {v_std:>10.4f}")

    # ── Panel 3: Density traces (convergence check) ───────────────────────
    if has_trace:
        ax3    = axes[ax_idx]; ax_idx += 1
        traces = read_rho_trace(trace_file)
        cmap   = plt.cm.coolwarm
        labels = sorted(traces.keys())
        n      = len(labels)
        for k, label in enumerate(labels):
            steps, rhos = traces[label]
            color = cmap(k / max(n - 1, 1))
            ax3.plot(np.arange(len(rhos)), rhos, '-', color=color,
                     linewidth=1.0, alpha=0.8, label=label)

        ax3.set_xlabel('Sample index (within each window)')
        ax3.set_ylabel(r'$\rho^*$ ($\sigma^{-3}$)')
        ax3.set_title('Density trace per state point  (coolwarm: high→low P*)')
        ax3.grid(alpha=0.3)
        if n <= 20:
            ax3.legend(ncol=2, fontsize=8)

    # ── Panel 4: Isothermal bulk modulus K_T vs P* ───────────────────────
    ax4 = axes[ax_idx]; ax_idx += 1
    ax4.plot(P, K_T, 's-', color='mediumseagreen', linewidth=1.8, markersize=5)
    ax4.set_xlabel(r'Pressure $P^*$')
    ax4.set_ylabel(r'$K_T^*$')
    ax4.set_title(r'Isothermal Bulk Modulus  $K_T = \rho \left(\partial P / \partial \rho\right)_T$')
    ax4.grid(alpha=0.3)

    # Print K_T table
    print(f"\n{'P*':>10}  {'K_T*':>12}")
    print('-' * 25)
    for p_val, k_val in zip(P, K_T):
        print(f"{p_val:>10.4f}  {k_val:>12.4f}")

    plt.tight_layout(pad=3.0)

    out_dir = os.path.join(foldername, 'output_plots', 'convergence_plots')
    os.makedirs(out_dir, exist_ok=True)
    outname = os.path.join(out_dir, f'{dataname}_{interaction}_eos.png')
    plt.savefig(outname, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to {outname}")


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    plot_eos(
        foldername  = sys.argv[1],
        dataname    = sys.argv[2],
        interaction = sys.argv[3],
    )