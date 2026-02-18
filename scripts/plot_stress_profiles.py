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
    # Extract base dataname without interaction and timesteps
    parts = dataname.split('_')
    
    # Find where interaction starts (format: number.number)
    base_parts = []
    for part in parts:
        if re.match(r'\d+\.\d+', part):  # Found interaction parameter
            break
        base_parts.append(part)
    
    base_name = '_'.join(base_parts)
    
    # Try multiple possible locations for the data file
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
    
    for dim in dims:
        poly_file = os.path.join(data_dir, f'stress_{dim}_polymer_{dataname}.dat')
        solv_file = os.path.join(data_dir, f'stress_{dim}_solvent_{dataname}.dat')
        if os.path.exists(poly_file) or os.path.exists(solv_file):
            return True
    return False

def plot_stress_profiles(folder, dataname, oldsteps):
    """Plot pressure profiles for polymer (left), solvent (middle), and total (right)."""
    box_dims = get_box_dims(folder, dataname)
    
    fig, axes = plt.subplots(3, 3, figsize=(20, 10))
    if oldsteps > 0:
        fig.suptitle(f'{dataname} (continuing from {oldsteps} steps)', fontsize=14, fontweight='bold')
    else:
        fig.suptitle(f'{dataname} (fresh run)', fontsize=14, fontweight='bold')
    
    labels = ['X', 'Y', 'Z']
    dims = ['x', 'y', 'z']
    binWidth = 2
    
    data_dir = os.path.join(folder, 'output_files', 'stress_data')
    
    polymer_ylims = [float('inf'), float('-inf')]
    solvent_ylims = [float('inf'), float('-inf')]
    total_ylims = [float('inf'), float('-inf')]
    
    # First pass: collect all timesteps to set up colormap
    all_timesteps = []
    poly_file = os.path.join(data_dir, f'stress_x_polymer_{dataname}.dat')
    if os.path.exists(poly_file):
        poly_data = read_ave_time_file(poly_file)
        plot_interval = max(1, len(poly_data) // 10)
        all_timesteps = [t for i, (t, _, _) in enumerate(poly_data) if i % plot_interval == 0]
    
    # Set up colormap normalized to timestep range
    if all_timesteps:
        cmap = plt.cm.viridis
        norm = Normalize(vmin=min(all_timesteps), vmax=max(all_timesteps))
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
    
    for row, (label, dim) in enumerate(zip(labels, dims)):
        poly_file = os.path.join(data_dir, f'stress_{dim}_polymer_{dataname}.dat')
        solv_file = os.path.join(data_dir, f'stress_{dim}_solvent_{dataname}.dat')
        
        poly_data = read_ave_time_file(poly_file) if os.path.exists(poly_file) else []
        solv_data = read_ave_time_file(solv_file) if os.path.exists(solv_file) else []
        
        if poly_data:
            plot_interval = max(1, len(poly_data) // 10)
            
            for i, (t, rows, P) in enumerate(poly_data):
                if i % plot_interval != 0:
                    continue
                coords_norm = (rows * binWidth - binWidth/2) / box_dims[dim]
                polymer_ylims[0] = min(polymer_ylims[0], P.min())
                polymer_ylims[1] = max(polymer_ylims[1], P.max())
                
                # Use continuous colormap based on timestep
                color = cmap(norm(t))
                alpha = 0.6 + 0.4 * (t - min(all_timesteps)) / (max(all_timesteps) - min(all_timesteps))
                
                axes[row, 0].plot(coords_norm, P, linewidth=2.0, alpha=alpha,
                                 color=color, label=f't={t}' if i // plot_interval < 3 else None)
            
            # Annotate first and last curves
            if len(poly_data) > 0:
                t_first, rows_first, P_first = poly_data[0]
                coords_first = (rows_first * binWidth - binWidth/2) / box_dims[dim]
                axes[row, 0].text(0.02, 0.98, f'First: t={t_first}', 
                                 transform=axes[row, 0].transAxes,
                                 fontsize=9, va='top', ha='left',
                                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                
                t_last = poly_data[-1][0]
                axes[row, 0].text(0.98, 0.98, f'Last: t={t_last}', 
                                 transform=axes[row, 0].transAxes,
                                 fontsize=9, va='top', ha='right',
                                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        
        axes[row, 0].set_ylabel(f'Partial stress ({label})', fontsize=11)
        axes[row, 0].set_xlabel(f'{label}/L{label}', fontsize=11)
        axes[row, 0].set_xlim(0, 1)
        axes[row, 0].grid(alpha=0.3)
        if row == 0:
            axes[row, 0].set_title('Polymer', fontweight='bold', fontsize=12)
        
        if solv_data:
            plot_interval = max(1, len(solv_data) // 10)
            
            for i, (t, rows, P) in enumerate(solv_data):
                if i % plot_interval != 0:
                    continue
                coords_norm = (rows * binWidth - binWidth/2) / box_dims[dim]
                solvent_ylims[0] = min(solvent_ylims[0], P.min())
                solvent_ylims[1] = max(solvent_ylims[1], P.max())
                
                color = cmap(norm(t))
                alpha = 0.6 + 0.4 * (t - min(all_timesteps)) / (max(all_timesteps) - min(all_timesteps))
                
                axes[row, 1].plot(coords_norm, P, linewidth=2.0, alpha=alpha, color=color)
            
            # Annotate first and last
            if len(solv_data) > 0:
                t_first = solv_data[0][0]
                t_last = solv_data[-1][0]
                axes[row, 1].text(0.02, 0.98, f'First: t={t_first}', 
                                 transform=axes[row, 1].transAxes,
                                 fontsize=9, va='top', ha='left',
                                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                axes[row, 1].text(0.98, 0.98, f'Last: t={t_last}', 
                                 transform=axes[row, 1].transAxes,
                                 fontsize=9, va='top', ha='right',
                                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        
        axes[row, 1].set_ylabel(f'Partial stress ({label})', fontsize=11)
        axes[row, 1].set_xlabel(f'{label}/L{label}', fontsize=11)
        axes[row, 1].set_xlim(0, 1)
        axes[row, 1].grid(alpha=0.3)
        if row == 0:
            axes[row, 1].set_title('Solvent', fontweight='bold', fontsize=12)
        
        # Total stress (interpolate and sum)
        if poly_data and solv_data:
            plot_interval = max(1, len(poly_data) // 10)
            
            for i, ((t_p, rows_p, P_p), (t_s, rows_s, P_s)) in enumerate(zip(poly_data, solv_data)):
                if i % plot_interval != 0:
                    continue
                
                coords_p = (rows_p * binWidth - binWidth/2) / box_dims[dim]
                coords_s = (rows_s * binWidth - binWidth/2) / box_dims[dim]
                
                # Create common grid
                coords_common = np.linspace(0, 1, 200)
                P_p_interp = np.interp(coords_common, coords_p, P_p, left=0, right=0)
                P_s_interp = np.interp(coords_common, coords_s, P_s, left=0, right=0)
                P_total = P_p_interp + P_s_interp
                
                total_ylims[0] = min(total_ylims[0], P_total.min())
                total_ylims[1] = max(total_ylims[1], P_total.max())
                
                color = cmap(norm(t_p))
                alpha = 0.6 + 0.4 * (t_p - min(all_timesteps)) / (max(all_timesteps) - min(all_timesteps))
                
                axes[row, 2].plot(coords_common, P_total, linewidth=2.0, alpha=alpha, color=color)
            
            # Annotate first and last
            if len(poly_data) > 0:
                t_first = poly_data[0][0]
                t_last = poly_data[-1][0]
                axes[row, 2].text(0.02, 0.98, f'First: t={t_first}', 
                                 transform=axes[row, 2].transAxes,
                                 fontsize=9, va='top', ha='left',
                                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                axes[row, 2].text(0.98, 0.98, f'Last: t={t_last}', 
                                 transform=axes[row, 2].transAxes,
                                 fontsize=9, va='top', ha='right',
                                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        
        axes[row, 2].set_ylabel(f'Total stress ({label})', fontsize=11)
        axes[row, 2].set_xlabel(f'{label}/L{label}', fontsize=11)
        axes[row, 2].set_xlim(0, 1)
        axes[row, 2].grid(alpha=0.3)
        if row == 0:
            axes[row, 2].set_title('Total', fontweight='bold', fontsize=12)
    
    # Set shared y-limits
    if polymer_ylims[0] != float('inf'):
        for row in range(3):
            axes[row, 0].set_ylim(polymer_ylims)
    if solvent_ylims[0] != float('inf'):
        for row in range(3):
            axes[row, 1].set_ylim(solvent_ylims)
    if total_ylims[0] != float('inf'):
        for row in range(3):
            axes[row, 2].set_ylim(total_ylims)
    
    # Add colorbar for timestep progression
    if all_timesteps:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Timestep', rotation=270, labelpad=20, fontsize=12)
    
    plt.tight_layout(rect=[0, 0, 0.9, 1])
    output_dir = os.path.join(folder, 'output_plots')
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f'{dataname}_stress.png'), dpi=150, bbox_inches='tight')
    print(f"Stress profile saved to {os.path.join(output_dir, f'{dataname}_stress.png')}")
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