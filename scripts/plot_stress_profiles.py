#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

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
    """Extract box dimensions from data file."""
    data_file = f'slab_with_flow/data_files/equil_{dataname}.data'

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
    data_dir = f'{folder}/output_files/stress_data'
    dims = ['x', 'y', 'z']
    
    for dim in dims:
        poly_file = f'{data_dir}/stress_{dim}_polymer_{dataname}.dat'
        solv_file = f'{data_dir}/stress_{dim}_solvent_{dataname}.dat'
        if os.path.exists(poly_file) or os.path.exists(solv_file):
            return True
    return False

def check_volume_data_exists(folder, dataname):
    """Check if any volume data files exist."""
    data_dir = f'{folder}/output_files/volume_data'
    dims = ['x', 'y', 'z']
    
    for dim in dims:
        poly_file = f'{data_dir}/vol_{dim}_polymer_{dataname}.dat'
        solv_file = f'{data_dir}/vol_{dim}_solvent_{dataname}.dat'
        if os.path.exists(poly_file) or os.path.exists(solv_file):
            return True
    return False

def plot_stress_profiles(folder, dataname, datasteps):
    """Plot pressure profiles for polymer (left), solvent (middle), and total (right)."""
    box_dims = get_box_dims(folder, dataname)
    
    fig, axes = plt.subplots(3, 3, figsize=(18, 10))
    fig.suptitle(f'{dataname} (original {datasteps} steps)', fontsize=14, fontweight='bold')
    
    labels = ['X', 'Y', 'Z']
    dims = ['x', 'y', 'z']
    colors = plt.cm.viridis(np.linspace(0, 1, 10))
    binWidth = 2
    
    data_dir = f'{folder}/output_files/stress_data'
    
    polymer_ylims = [float('inf'), float('-inf')]
    solvent_ylims = [float('inf'), float('-inf')]
    total_ylims = [float('inf'), float('-inf')]
    
    for row, (label, dim) in enumerate(zip(labels, dims)):
        poly_file = f'{data_dir}/stress_{dim}_polymer_{dataname}.dat'
        solv_file = f'{data_dir}/stress_{dim}_solvent_{dataname}.dat'
        
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
                
                axes[row, 0].plot(coords_norm, P, linewidth=1.5, alpha=0.7,
                                 color=colors[(i // plot_interval) % len(colors)],
                                 label=f't={t}' if row == 0 else None)
        
        axes[row, 0].set_ylabel(f'Partial stress ({label})')
        axes[row, 0].set_xlabel(f'{label}/L{label}')
        axes[row, 0].grid(alpha=0.3)
        if row == 0:
            axes[row, 0].legend(loc='best', fontsize=7, ncol=2)
            axes[row, 0].set_title('Polymer', fontweight='bold')
        
        if solv_data:
            plot_interval = max(1, len(solv_data) // 10)
            
            for i, (t, rows, P) in enumerate(solv_data):
                if i % plot_interval != 0:
                    continue
                coords_norm = (rows * binWidth - binWidth/2) / box_dims[dim]
                solvent_ylims[0] = min(solvent_ylims[0], P.min())
                solvent_ylims[1] = max(solvent_ylims[1], P.max())
                
                axes[row, 1].plot(coords_norm, P, linewidth=1.5, alpha=0.7,
                                 color=colors[(i // plot_interval) % len(colors)])
        
        axes[row, 1].set_ylabel(f'Partial stress ({label})')
        axes[row, 1].set_xlabel(f'{label}/L{label}')
        axes[row, 1].set_xlim(0, 1)
        axes[row, 1].grid(alpha=0.3)
        if row == 0:
            axes[row, 1].set_title('Solvent', fontweight='bold')
        
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
                
                axes[row, 2].plot(coords_common, P_total, linewidth=1.5, alpha=0.7,
                                 color=colors[(i // plot_interval) % len(colors)])
        
        axes[row, 2].set_ylabel(f'Total stress ({label})')
        axes[row, 2].set_xlabel(f'{label}/L{label}')
        axes[row, 2].set_xlim(0, 1)
        axes[row, 2].grid(alpha=0.3)
        if row == 0:
            axes[row, 2].set_title('Total', fontweight='bold')
    
    if polymer_ylims[0] != float('inf'):
        for row in range(3):
            axes[row, 0].set_ylim(polymer_ylims)
    if solvent_ylims[0] != float('inf'):
        for row in range(3):
            axes[row, 1].set_ylim(solvent_ylims)
    if total_ylims[0] != float('inf'):
        for row in range(3):
            axes[row, 2].set_ylim(total_ylims)
    
    plt.tight_layout()
    plt.savefig(f'{folder}/output_plots/stress_plots/{dataname}_stress.png', dpi=150)
    print(f"Stress profile saved to {folder}/output_plots/stress_plots/{dataname}_stress.png")
    plt.close()

def plot_volume_fraction_profiles(folder, dataname, datasteps):
    """Plot volume fraction profiles for polymer (left), solvent (middle), and total (right)."""
    box_dims = get_box_dims(folder, dataname)
    
    fig, axes = plt.subplots(3, 3, figsize=(18, 10))
    fig.suptitle(f'{dataname} Volume Fractions (original {datasteps} steps)', fontsize=14, fontweight='bold')
    
    labels = ['X', 'Y', 'Z']
    dims = ['x', 'y', 'z']
    colors = plt.cm.viridis(np.linspace(0, 1, 10))
    binWidth = 0.5
    
    data_dir = f'{folder}/output_files/volume_data'
    
    polymer_ylims = [float('inf'), float('-inf')]
    solvent_ylims = [float('inf'), float('-inf')]
    total_ylims = [float('inf'), float('-inf')]
    
    for row, (label, dim) in enumerate(zip(labels, dims)):
        # Calculate bin volume for this direction
        if dim == 'x':
            bin_volume = binWidth * box_dims['y'] * box_dims['z']
        elif dim == 'y':
            bin_volume = box_dims['x'] * binWidth * box_dims['z']
        else:  # z
            bin_volume = box_dims['x'] * box_dims['y'] * binWidth
        
        poly_file = f'{data_dir}/vol_{dim}_polymer_{dataname}.dat'
        solv_file = f'{data_dir}/vol_{dim}_solvent_{dataname}.dat'
        
        poly_data = read_ave_time_file(poly_file) if os.path.exists(poly_file) else []
        solv_data = read_ave_time_file(solv_file) if os.path.exists(solv_file) else []
        
        if poly_data:
            plot_interval = max(1, len(poly_data) // 10)
            
            for i, (t, rows, V) in enumerate(poly_data):
                if i % plot_interval != 0:
                    continue
                coords_norm = (rows * binWidth - binWidth/2) / box_dims[dim]
                phi = V / bin_volume
                polymer_ylims[0] = min(polymer_ylims[0], phi.min())
                polymer_ylims[1] = max(polymer_ylims[1], phi.max())
                
                axes[row, 0].plot(coords_norm, phi, linewidth=1.5, alpha=0.7,
                                 color=colors[(i // plot_interval) % len(colors)],
                                 label=f't={t}' if row == 0 else None)
        
        axes[row, 0].set_ylabel(f'Volume fraction ({label})')
        axes[row, 0].set_xlabel(f'{label}/L{label}')
        axes[row, 0].grid(alpha=0.3)
        if row == 0:
            axes[row, 0].legend(loc='best', fontsize=7, ncol=2)
            axes[row, 0].set_title('Polymer', fontweight='bold')
        
        if solv_data:
            plot_interval = max(1, len(solv_data) // 10)
            
            for i, (t, rows, V) in enumerate(solv_data):
                if i % plot_interval != 0:
                    continue
                coords_norm = (rows * binWidth - binWidth/2) / box_dims[dim]
                phi = V / bin_volume
                solvent_ylims[0] = min(solvent_ylims[0], phi.min())
                solvent_ylims[1] = max(solvent_ylims[1], phi.max())
                
                axes[row, 1].plot(coords_norm, phi, linewidth=1.5, alpha=0.7,
                                 color=colors[(i // plot_interval) % len(colors)])
        
        axes[row, 1].set_ylabel(f'Volume fraction ({label})')
        axes[row, 1].set_xlabel(f'{label}/L{label}')
        axes[row, 1].set_xlim(0, 1)
        axes[row, 1].grid(alpha=0.3)
        if row == 0:
            axes[row, 1].set_title('Solvent', fontweight='bold')
        
        # Total volume fraction (interpolate and sum)
        if poly_data and solv_data:
            plot_interval = max(1, len(poly_data) // 10)
            
            for i, ((t_p, rows_p, V_p), (t_s, rows_s, V_s)) in enumerate(zip(poly_data, solv_data)):
                if i % plot_interval != 0:
                    continue
                
                coords_p = (rows_p * binWidth - binWidth/2) / box_dims[dim]
                coords_s = (rows_s * binWidth - binWidth/2) / box_dims[dim]
                phi_p = V_p / bin_volume
                phi_s = V_s / bin_volume
                
                # Create common grid
                coords_common = np.linspace(0, 1, 200)
                phi_p_interp = np.interp(coords_common, coords_p, phi_p, left=0, right=0)
                phi_s_interp = np.interp(coords_common, coords_s, phi_s, left=0, right=0)
                phi_total = phi_p_interp + phi_s_interp
                
                total_ylims[0] = min(total_ylims[0], phi_total.min())
                total_ylims[1] = max(total_ylims[1], phi_total.max())
                
                axes[row, 2].plot(coords_common, phi_total, linewidth=1.5, alpha=0.7,
                                 color=colors[(i // plot_interval) % len(colors)])
        
        axes[row, 2].set_ylabel(f'Total volume fraction ({label})')
        axes[row, 2].set_xlabel(f'{label}/L{label}')
        axes[row, 2].set_xlim(0, 1)
        axes[row, 2].grid(alpha=0.3)
        if row == 0:
            axes[row, 2].set_title('Total', fontweight='bold')
    
    if polymer_ylims[0] != float('inf'):
        for row in range(3):
            axes[row, 0].set_ylim(polymer_ylims)
    if solvent_ylims[0] != float('inf'):
        for row in range(3):
            axes[row, 1].set_ylim(solvent_ylims)
    if total_ylims[0] != float('inf'):
        for row in range(3):
            axes[row, 2].set_ylim(total_ylims)
    
    plt.tight_layout()
    plt.savefig(f'{folder}/output_plots/volfrac_plots/{dataname}_volume.png', dpi=150)
    print(f"Volume fraction profile saved to {folder}/output_plots/volfrac_plots/{dataname}_volume.png")
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python plot_stress_profiles.py <folder> <dataname> <datasteps>")
        sys.exit(1)
    
    folder = sys.argv[1]
    dataname = sys.argv[2]
    datasteps = sys.argv[3] if len(sys.argv) > 3 else 'N/A'
    
    # Check if data exists before creating plots
    stress_exists = check_stress_data_exists(folder, dataname)
    volume_exists = check_volume_data_exists(folder, dataname)
    
    if stress_exists:
        os.makedirs(f'{folder}/output_plots/stress_plots', exist_ok=True)
        plot_stress_profiles(folder, dataname, datasteps)
    else:
        print(f"No stress data found for {dataname}, skipping stress plots")
    
    if volume_exists:
        os.makedirs(f'{folder}/output_plots/volfrac_plots', exist_ok=True)
        plot_volume_fraction_profiles(folder, dataname, datasteps)
    else:
        print(f"No volume data found for {dataname}, skipping volume fraction plots")
