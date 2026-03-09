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
    """Plot pressure profiles: polymer, solvent, piston, support (cols 0-3), and total (col 4)."""
    box_dims = get_box_dims(folder, dataname)
    
    # 3 rows (x, y, z) x 5 cols (polymer, solvent, piston, support, total)
    fig, axes = plt.subplots(3, 5, figsize=(30, 10))
    if oldsteps > 0:
        fig.suptitle(f'{dataname} (continuing from {oldsteps} steps)', fontsize=14, fontweight='bold')
    else:
        fig.suptitle(f'{dataname} (fresh run)', fontsize=14, fontweight='bold')
    
    labels = ['X', 'Y', 'Z']
    dims = ['x', 'y', 'z']
    binWidth = 2
    col_titles = ['Polymer', 'Solvent', 'Piston', 'Support', 'Total']
    components = ['polymer', 'solvent', 'piston', 'support']
    
    data_dir = os.path.join(folder, 'output_files', 'stress_data')
    
    # Per-column shared y-limits (one per component + total)
    ylims = [[float('inf'), float('-inf')] for _ in range(5)]
    
    # First pass: collect all timesteps from polymer x file for colormap
    all_timesteps = []
    poly_file = os.path.join(data_dir, f'stress_x_polymer_{dataname}.dat')
    if os.path.exists(poly_file):
        poly_data = read_ave_time_file(poly_file)
        plot_interval = max(1, len(poly_data) // 10)
        all_timesteps = [t for i, (t, _, _) in enumerate(poly_data) if i % plot_interval == 0]
    
    if all_timesteps:
        cmap = plt.cm.viridis
        norm = Normalize(vmin=min(all_timesteps), vmax=max(all_timesteps))
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
    else:
        cmap = plt.cm.viridis
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

    for row, (label, dim) in enumerate(zip(labels, dims)):
        # Load data for all four components
        comp_data = {}
        for comp in components:
            fpath = os.path.join(data_dir, f'stress_{dim}_{comp}_{dataname}.dat')
            comp_data[comp] = read_ave_time_file(fpath) if os.path.exists(fpath) else []

        # Determine plot_interval from whichever component has the most frames
        max_frames = max((len(v) for v in comp_data.values()), default=1)
        plot_interval = max(1, max_frames // 10)

        # ----- Columns 0-3: individual components -----
        for col_idx, comp in enumerate(components):
            data = comp_data[comp]
            ax = axes[row, col_idx]

            # Piston and support are z-only; mark x/y panels as not measured
            if comp in ('piston', 'support') and dim in ('x', 'y'):
                ax.text(0.5, 0.5, 'not measured', transform=ax.transAxes,
                        fontsize=11, color='gray', ha='center', va='center',
                        style='italic')
                ax.set_facecolor('#f5f5f5')
                ax.tick_params(left=False, bottom=False,
                               labelleft=False, labelbottom=False)
                if row == 0:
                    ax.set_title(col_titles[col_idx], fontweight='bold', fontsize=12)
                continue

            for i, (t, rows, P) in enumerate(data):
                if i % plot_interval != 0:
                    continue
                coords_norm = (rows * binWidth - binWidth / 2) / box_dims[dim]
                ylims[col_idx][0] = min(ylims[col_idx][0], P.min())
                ylims[col_idx][1] = max(ylims[col_idx][1], P.max())
                color = cmap(norm(t))
                ax.plot(coords_norm, P, linewidth=2.0, alpha=_alpha(t), color=color)

            _annotate(ax, data)
            ax.set_ylabel(f'Partial stress ({label})', fontsize=11)
            ax.set_xlabel(f'{label}/L{label}', fontsize=11)
            ax.set_xlim(0, 1)
            ax.grid(alpha=0.3)
            if row == 0:
                ax.set_title(col_titles[col_idx], fontweight='bold', fontsize=12)

        # ----- Column 4: total stress (sum of all four components) -----
        ax_tot = axes[row, 4]

        # Find common timesteps across all components that have data
        available = {comp: comp_data[comp] for comp in components if comp_data[comp]}
        if available:
            # Build dict of {timestep: (rows, P)} for each component
            ts_maps = {}
            for comp, data in available.items():
                ts_maps[comp] = {t: (rows, P) for t, rows, P in data}

            # Use timesteps from the component with the most frames
            ref_comp = max(ts_maps, key=lambda c: len(ts_maps[c]))
            ref_timesteps = sorted(ts_maps[ref_comp].keys())
            ref_plotted = [t for i, t in enumerate(ref_timesteps) if i % plot_interval == 0]

            first_total_data = None
            last_total_data = None

            for t in ref_plotted:
                coords_common = np.linspace(0, 1, 200)
                P_total = np.zeros(200)

                for comp, ts_map in ts_maps.items():
                    # Find nearest available timestep for this component
                    if not ts_map:
                        continue
                    nearest_t = min(ts_map.keys(), key=lambda tt: abs(tt - t))
                    rows_c, P_c = ts_map[nearest_t]
                    coords_c = (rows_c * binWidth - binWidth / 2) / box_dims[dim]
                    P_total += np.interp(coords_common, coords_c, P_c, left=0, right=0)

                ylims[4][0] = min(ylims[4][0], P_total.min())
                ylims[4][1] = max(ylims[4][1], P_total.max())
                color = cmap(norm(t))
                ax_tot.plot(coords_common, P_total, linewidth=2.0, alpha=_alpha(t), color=color)

                if first_total_data is None:
                    first_total_data = t
                last_total_data = t

            if first_total_data is not None:
                ax_tot.text(0.02, 0.98, f'First: t={first_total_data}',
                            transform=ax_tot.transAxes, fontsize=9, va='top', ha='left',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                ax_tot.text(0.98, 0.98, f'Last: t={last_total_data}',
                            transform=ax_tot.transAxes, fontsize=9, va='top', ha='right',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        ax_tot.set_ylabel(f'Total stress ({label})', fontsize=11)
        ax_tot.set_xlabel(f'{label}/L{label}', fontsize=11)
        ax_tot.set_xlim(0, 1)
        ax_tot.grid(alpha=0.3)
        if row == 0:
            ax_tot.set_title('Total', fontweight='bold', fontsize=12)

    # Apply shared y-limits per column
    for col_idx in range(5):
        if ylims[col_idx][0] != float('inf'):
            for row in range(3):
                axes[row, col_idx].set_ylim(ylims[col_idx])

    # Colorbar
    if all_timesteps:
        cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
        cbar = fig.colorbar(sm, cax=cbar_ax)
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
    
    folder = sys.argv[1]
    dataname = sys.argv[2]
    oldsteps = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    
    if check_stress_data_exists(folder, dataname):
        plot_stress_profiles(folder, dataname, oldsteps)
    else:
        print(f"No stress data found for {dataname}, skipping stress plots")