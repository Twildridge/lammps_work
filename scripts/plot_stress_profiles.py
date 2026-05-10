#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import sys
import os
import re

# HOW TO RUN BY ITSELF
# cd ~/Documents/lammps_runs/slab_with_flow_walled.....
#
# python ~/Documents/lammps_work/scripts/plot_stress_profiles.py . \
# walled_slab_support_5beads_tall_4_1.0_1.05_15000000_1.0_1.05_1000000 0



def read_ave_time_file(filepath):
    """Read LAMMPS ave/time output file with format: timestep nrows, then row pressure."""
    data_by_time = []

    with open(filepath, 'r') as f:
        lines = [line for line in f if not line.startswith('#') and line.strip()]

        i = 0
        while i < len(lines):
            parts = lines[i].split()
            if len(parts) == 2:  # Timestep line
                timestep = int(parts[0])
                nrows = int(parts[1])
                pressures = []

                for j in range(1, nrows + 1):
                    if i + j < len(lines):
                        p_parts = lines[i + j].split()
                        if len(p_parts) == 2:
                            pressures.append(float(p_parts[1]))

                if pressures:
                    rows = np.arange(1, len(pressures) + 1)
                    data_by_time.append((timestep, rows, np.array(pressures)))

                i += nrows + 1
            else:
                i += 1

    return data_by_time

def get_box_dims(folder, dataname):
    """Extract box dimensions from data file in the working directory."""
    parts = dataname.split('_')

    base_parts = []
    for part in parts:
        if re.match(r'\d+\.\d+', part):
            break
        base_parts.append(part)

    base_name = '_'.join(base_parts)

    possible_paths = [
        os.path.join(folder, 'data_files', f'{base_name}.data'),
        os.path.join(folder, f'final_config_{dataname}.data'),
        os.path.join(folder, f'final_flow_{dataname}.data'),
        os.path.join(folder, 'data_files', f'{dataname}.data'),
    ]

    data_file = None
    for path in possible_paths:
        if os.path.exists(path):
            data_file = path
            print(f"Found data file: {path}")
            break

    if not data_file:
        print(f"Warning: Could not find data file, using default box dimensions")
        return {'x': 100.0, 'y': 100.0, 'z': 50.0}

    box_dims = {}
    with open(data_file, 'r') as f:
        for line in f:
            if 'xlo xhi' in line:
                vals = line.split()
                box_dims['x'] = float(vals[1]) - float(vals[0])
            elif 'ylo yhi' in line:
                vals = line.split()
                box_dims['y'] = float(vals[1]) - float(vals[0])
            elif 'zlo zhi' in line:
                vals = line.split()
                box_dims['z'] = float(vals[1]) - float(vals[0])
                break
    return box_dims

def check_stress_data_exists(folder, dataname):
    """Check if any stress data files exist."""
    data_dir = os.path.join(folder, 'output_files', 'stress_data')
    dims = ['x', 'y', 'z']
    components = ['polymer', 'solvent', 'piston', 'support']

    for dim in dims:
        for comp in components:
            f = os.path.join(data_dir, f'stress_{dim}_{comp}_{dataname}.dat')
            if os.path.exists(f):
                return True
    return False

def plot_stress_profiles(folder, dataname, oldsteps):
    """Plot pressure profiles.

    Rows 0-2  (existing): isotropic partial pressure binned by x, y, z.
               Columns: Polymer | Solvent | Piston | Support | Total
    Rows 3-5  (new):     individual diagonal stress components, all z-binned.
               Row 3: σ_zz  — polymer, solvent, [n/a], [n/a], total
               Row 4: σ_xx  — polymer, solvent, [n/a], [n/a], total
               Row 5: σ_yy  — polymer, solvent, [n/a], [n/a], total
    """
    box_dims = get_box_dims(folder, dataname)

    # 6 rows x 5 cols
    fig, axes = plt.subplots(6, 5, figsize=(30, 20))
    if oldsteps > 0:
        fig.suptitle(f'{dataname} (continuing from {oldsteps} steps)', fontsize=14, fontweight='bold')
    else:
        fig.suptitle(f'{dataname} (fresh run)', fontsize=14, fontweight='bold')

    labels      = ['X', 'Y', 'Z']
    dims        = ['x', 'y', 'z']
    binWidth    = 2
    col_titles  = ['Polymer', 'Solvent', 'Piston', 'Support', 'Total']
    components  = ['polymer', 'solvent', 'piston', 'support']

    data_dir = os.path.join(folder, 'output_files', 'stress_data')

    # Per-column shared y-limits for the isotropic rows (rows 0-2)
    ylims = [[float('inf'), float('-inf')] for _ in range(5)]

    # Per-(component, col) y-limits for the diagonal rows (rows 3-5)
    # Keys: 'zz', 'xx', 'yy'; values: list of [lo, hi] per column
    diag_ylims = {k: [[float('inf'), float('-inf')] for _ in range(5)]
                  for k in ('zz', 'xx', 'yy')}

    # Colormap: built from all timesteps in the polymer-x isotropic file
    all_timesteps = []
    poly_file = os.path.join(data_dir, f'stress_x_polymer_{dataname}.dat')
    if os.path.exists(poly_file):
        poly_data = read_ave_time_file(poly_file)
        plot_interval = max(1, len(poly_data) // 10)
        all_timesteps = [t for i, (t, _, _) in enumerate(poly_data) if i % plot_interval == 0]

    cmap = plt.cm.viridis
    if all_timesteps:
        norm = Normalize(vmin=min(all_timesteps), vmax=max(all_timesteps))
    else:
        norm = Normalize(vmin=0, vmax=1)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    def _alpha(t):
        if len(all_timesteps) < 2:
            return 0.8
        return 0.6 + 0.4 * (t - min(all_timesteps)) / (max(all_timesteps) - min(all_timesteps))

    def _annotate(ax, data):
        """Add first/last timestep labels to an axis."""
        if not data:
            return
        ax.text(0.02, 0.98, f'First: t={data[0][0]}',
                transform=ax.transAxes, fontsize=9, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        ax.text(0.98, 0.98, f'Last: t={data[-1][0]}',
                transform=ax.transAxes, fontsize=9, va='top', ha='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    def _not_measured(ax):
        ax.text(0.5, 0.5, 'not measured', transform=ax.transAxes,
                fontsize=11, color='gray', ha='center', va='center', style='italic')
        ax.set_facecolor('#f5f5f5')
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # ──────────────────────────────────────────────────────────────────────
    # ROWS 0-2: isotropic partial pressure, binned by x / y / z
    # ──────────────────────────────────────────────────────────────────────
    for row, (label, dim) in enumerate(zip(labels, dims)):
        comp_data = {}
        for comp in components:
            fpath = os.path.join(data_dir, f'stress_{dim}_{comp}_{dataname}.dat')
            comp_data[comp] = read_ave_time_file(fpath) if os.path.exists(fpath) else []

        max_frames    = max((len(v) for v in comp_data.values()), default=1)
        plot_interval = max(1, max_frames // 10)

        # Columns 0-3: individual components
        for col_idx, comp in enumerate(components):
            data = comp_data[comp]
            ax   = axes[row, col_idx]

            if comp in ('piston', 'support') and dim in ('x', 'y'):
                _not_measured(ax)
                if row == 0:
                    ax.set_title(col_titles[col_idx], fontweight='bold', fontsize=12)
                continue

            for i, (t, rows, P) in enumerate(data):
                if i % plot_interval != 0:
                    continue
                coords_norm = (rows * binWidth - binWidth / 2) / box_dims[dim]
                ylims[col_idx][0] = min(ylims[col_idx][0], P.min())
                ylims[col_idx][1] = max(ylims[col_idx][1], P.max())
                ax.plot(coords_norm, P, linewidth=2.0, alpha=_alpha(t), color=cmap(norm(t)))

            _annotate(ax, data)
            ax.set_ylabel(f'Partial stress ({label})', fontsize=11)
            ax.set_xlabel(f'{label}/L{label}', fontsize=11)
            ax.set_xlim(0, 1)
            ax.grid(alpha=0.3)
            if row == 0:
                ax.set_title(col_titles[col_idx], fontweight='bold', fontsize=12)

        # Column 4: total
        ax_tot    = axes[row, 4]
        available = {comp: comp_data[comp] for comp in components if comp_data[comp]}
        if available:
            ts_maps = {comp: {t: (r, P) for t, r, P in data}
                       for comp, data in available.items()}
            ref_comp      = max(ts_maps, key=lambda c: len(ts_maps[c]))
            ref_timesteps = sorted(ts_maps[ref_comp].keys())
            ref_plotted   = [t for i, t in enumerate(ref_timesteps) if i % plot_interval == 0]
            first_t = last_t = None

            for t in ref_plotted:
                coords_common = np.linspace(0, 1, 200)
                P_total       = np.zeros(200)
                for comp, ts_map in ts_maps.items():
                    nearest_t    = min(ts_map.keys(), key=lambda tt: abs(tt - t))
                    rows_c, P_c  = ts_map[nearest_t]
                    coords_c     = (rows_c * binWidth - binWidth / 2) / box_dims[dim]
                    P_total     += np.interp(coords_common, coords_c, P_c, left=0, right=0)
                ylims[4][0] = min(ylims[4][0], P_total.min())
                ylims[4][1] = max(ylims[4][1], P_total.max())
                ax_tot.plot(coords_common, P_total, linewidth=2.0,
                            alpha=_alpha(t), color=cmap(norm(t)))
                if first_t is None:
                    first_t = t
                last_t = t

            if first_t is not None:
                ax_tot.text(0.02, 0.98, f'First: t={first_t}',
                            transform=ax_tot.transAxes, fontsize=9, va='top', ha='left',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                ax_tot.text(0.98, 0.98, f'Last: t={last_t}',
                            transform=ax_tot.transAxes, fontsize=9, va='top', ha='right',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        ax_tot.set_ylabel(f'Total stress ({label})', fontsize=11)
        ax_tot.set_xlabel(f'{label}/L{label}', fontsize=11)
        ax_tot.set_xlim(0, 1)
        ax_tot.grid(alpha=0.3)
        if row == 0:
            ax_tot.set_title('Total', fontweight='bold', fontsize=12)

    # ──────────────────────────────────────────────────────────────────────
    # ROWS 3-5: individual diagonal components (σ_zz, σ_xx, σ_yy), z-binned
    # Columns 0, 1 = polymer, solvent.  Columns 2, 3 = not measured.
    # Column 4 = polymer + solvent total for that component.
    # ──────────────────────────────────────────────────────────────────────
    diag_specs = [
        # (row_index, comp_key, poly_file_prefix, solv_file_prefix, y_label, col4_label)
        (3, 'zz', 'sigmazz_polymer', 'sigmazz_solvent',
         r'$\sigma_{zz}$ (z-binned)', r'Total $\sigma_{zz}$'),
        (4, 'xx', 'sigmaxx_polymer', 'sigmaxx_solvent',
         r'$\sigma_{xx}$ (z-binned)', r'Total $\sigma_{xx}$'),
        (5, 'yy', 'sigmayy_polymer', 'sigmayy_solvent',
         r'$\sigma_{yy}$ (z-binned)', r'Total $\sigma_{yy}$'),
    ]

    for row, comp_key, poly_prefix, solv_prefix, ylabel, tot_ylabel in diag_specs:
        poly_path = os.path.join(data_dir, f'{poly_prefix}_{dataname}.dat')
        solv_path = os.path.join(data_dir, f'{solv_prefix}_{dataname}.dat')
        poly_data = read_ave_time_file(poly_path) if os.path.exists(poly_path) else []
        solv_data = read_ave_time_file(solv_path) if os.path.exists(solv_path) else []

        comp_datasets = {'polymer': poly_data, 'solvent': solv_data}
        max_frames    = max(len(poly_data), len(solv_data), 1)
        plot_interval = max(1, max_frames // 10)

        # Columns 0 (polymer) and 1 (solvent)
        for col_idx, (comp_label, data) in enumerate(
                [('Polymer', poly_data), ('Solvent', solv_data)]):
            ax = axes[row, col_idx]

            for i, (t, rows, P) in enumerate(data):
                if i % plot_interval != 0:
                    continue
                coords_norm = (rows * binWidth - binWidth / 2) / box_dims['z']
                diag_ylims[comp_key][col_idx][0] = min(diag_ylims[comp_key][col_idx][0], P.min())
                diag_ylims[comp_key][col_idx][1] = max(diag_ylims[comp_key][col_idx][1], P.max())
                ax.plot(coords_norm, P, linewidth=2.0,
                        alpha=_alpha(t), color=cmap(norm(t)))

            _annotate(ax, data)
            ax.set_ylabel(ylabel, fontsize=11)
            ax.set_xlabel('z/Lz', fontsize=11)
            ax.set_xlim(0, 1)
            ax.grid(alpha=0.3)
            if row == 3:
                ax.set_title(col_titles[col_idx], fontweight='bold', fontsize=12)

        # Columns 2, 3: not measured
        for col_idx in (2, 3):
            ax = axes[row, col_idx]
            _not_measured(ax)
            if row == 3:
                ax.set_title(col_titles[col_idx], fontweight='bold', fontsize=12)

        # Column 4: total = polymer + solvent for this component
        ax_tot = axes[row, 4]
        avail  = {k: v for k, v in comp_datasets.items() if v}
        if avail:
            ts_maps = {k: {t: (r, P) for t, r, P in v} for k, v in avail.items()}
            ref_key  = max(ts_maps, key=lambda k: len(ts_maps[k]))
            ref_ts   = sorted(ts_maps[ref_key].keys())
            plotted  = [t for i, t in enumerate(ref_ts) if i % plot_interval == 0]
            first_t  = last_t = None

            for t in plotted:
                coords_common = np.linspace(0, 1, 200)
                P_total       = np.zeros(200)
                for k, ts_map in ts_maps.items():
                    nearest_t   = min(ts_map.keys(), key=lambda tt: abs(tt - t))
                    rows_c, P_c = ts_map[nearest_t]
                    coords_c    = (rows_c * binWidth - binWidth / 2) / box_dims['z']
                    P_total    += np.interp(coords_common, coords_c, P_c, left=0, right=0)
                diag_ylims[comp_key][4][0] = min(diag_ylims[comp_key][4][0], P_total.min())
                diag_ylims[comp_key][4][1] = max(diag_ylims[comp_key][4][1], P_total.max())
                ax_tot.plot(coords_common, P_total, linewidth=2.0,
                            alpha=_alpha(t), color=cmap(norm(t)))
                if first_t is None:
                    first_t = t
                last_t = t

            if first_t is not None:
                ax_tot.text(0.02, 0.98, f'First: t={first_t}',
                            transform=ax_tot.transAxes, fontsize=9, va='top', ha='left',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                ax_tot.text(0.98, 0.98, f'Last: t={last_t}',
                            transform=ax_tot.transAxes, fontsize=9, va='top', ha='right',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        ax_tot.set_ylabel(tot_ylabel, fontsize=11)
        ax_tot.set_xlabel('z/Lz', fontsize=11)
        ax_tot.set_xlim(0, 1)
        ax_tot.grid(alpha=0.3)
        if row == 3:
            ax_tot.set_title('Total', fontweight='bold', fontsize=12)

    # ── Apply shared y-limits ──────────────────────────────────────────────
    # Isotropic rows (0-2): one shared limit per column
    for col_idx in range(5):
        if ylims[col_idx][0] != float('inf'):
            for row in range(3):
                axes[row, col_idx].set_ylim(ylims[col_idx])

    # Diagonal rows (3-5): each component gets its own shared limit per column
    for (row, comp_key, *_) in diag_specs:
        for col_idx in range(5):
            lo, hi = diag_ylims[comp_key][col_idx]
            if lo != float('inf'):
                axes[row, col_idx].set_ylim(lo, hi)

    # ── Row-section labels in the left margin ─────────────────────────────
    section_labels = {
        0: 'Isotropic\n(x-binned)',
        1: 'Isotropic\n(y-binned)',
        2: 'Isotropic\n(z-binned)',
        3: r'$\sigma_{zz}$' + '\n(z-binned)',
        4: r'$\sigma_{xx}$' + '\n(z-binned)',
        5: r'$\sigma_{yy}$' + '\n(z-binned)',
    }
    for row_idx, section_label in section_labels.items():
        axes[row_idx, 0].set_ylabel(
            section_label + '\n' + axes[row_idx, 0].get_ylabel(),
            fontsize=11)

    # ── Colorbar ──────────────────────────────────────────────────────────
    if all_timesteps:
        cbar_ax = fig.add_axes([0.93, 0.15, 0.012, 0.7])
        cbar    = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Timestep', rotation=270, labelpad=20, fontsize=12)

    plt.tight_layout(rect=[0, 0, 0.92, 1])
    output_dir = os.path.join(folder, 'output_plots')
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f'{dataname}_stress.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Stress profile saved to {out_path}")
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python plot_stress_profiles.py <folder> <dataname> <oldsteps>")
        sys.exit(1)

    folder   = sys.argv[1]
    dataname = sys.argv[2]
    oldsteps = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    if check_stress_data_exists(folder, dataname):
        plot_stress_profiles(folder, dataname, oldsteps)
    else:
        print(f"No stress data found for {dataname}, skipping stress plots")
