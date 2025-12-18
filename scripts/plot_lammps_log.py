#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def read_volume_file(filepath):
    """Read single-column volume data."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    data.append(float(line.split()[0]))
                except (ValueError, IndexError):
                    continue
    return np.array(data)

def read_timestep_volume_file(filepath):
    """Read two-column timestep + volume data."""
    timesteps = []
    volumes = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    parts = line.split()
                    timesteps.append(float(parts[0]))
                    volumes.append(float(parts[1]))
                except (ValueError, IndexError):
                    continue
    return np.array(timesteps), np.array(volumes)

def parse_lammps_log(filepath='log.lammps'):
    """Parse LAMMPS log file and extract thermo data."""
    data = {}
    reading_data = False
    headers = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('Step'):
                headers = line.split()
                reading_data = True
                for h in headers:
                    data[h] = []
                continue
            
            if reading_data and ('Loop time' in line or line.startswith('WARNING')):
                reading_data = False
                continue
            
            if reading_data and line and not line.startswith('#'):
                try:
                    values = line.split()
                    if len(values) == len(headers):
                        for h, v in zip(headers, values):
                            data[h].append(float(v))
                except ValueError:
                    continue
    
    for key in data:
        data[key] = np.array(data[key])
    
    return data

def plot_convergence(data, foldername, dataname, output='convergence.png'):
    """Plot temperature, pressure, normalized box volume, and gel volumes."""
    
    # Try to read volume files
    box_vol_file = os.path.join(foldername, 'output_files/volume_data', 
                                 f'box_dimensions_{dataname}.dat')
    gel_bb_file = os.path.join(foldername, 'output_files/volume_data', 
                                f'gel_volume_bb_{dataname}.dat')
    gel_rg_file = os.path.join(foldername, 'output_files/volume_data', 
                                f'gel_volume_rg_{dataname}.dat')
    
    has_box = os.path.exists(box_vol_file)
    has_gel_bb = os.path.exists(gel_bb_file)
    has_gel_rg = os.path.exists(gel_rg_file)
    
    num_plots = 2  # temp + pressure always
    if has_box:
        num_plots += 1
    if has_gel_bb:
        num_plots += 1
    if has_gel_rg:
        num_plots += 1
    
    fig, axes = plt.subplots(num_plots, 1, figsize=(10, 3*num_plots))
    if num_plots == 1:
        axes = [axes]
    
    fig.suptitle(dataname, fontsize=14, fontweight='bold')
    
    plot_idx = 0
    
    # Temperature
    if 'Temp' in data:
        axes[plot_idx].plot(data['Step'], data['Temp'], 'b-', linewidth=2.0)
        axes[plot_idx].set_ylabel('Temperature')
        axes[plot_idx].grid(alpha=0.3)
        n_last = int(len(data['Temp']) * 0.3)
        if n_last > 10:
            last_mean = data['Temp'][-n_last:].mean()
            last_std = data['Temp'][-n_last:].std()
            axes[plot_idx].axhline(y=last_mean, color='r', linestyle='--', 
                                   label=f'Last 30%: {last_mean:.3f} ± {last_std:.3f}')
            axes[plot_idx].legend()
        plot_idx += 1
    
    # Pressure
    if 'Press' in data:
        axes[plot_idx].plot(data['Step'], data['Press'], 'g-', linewidth=2.0)
        axes[plot_idx].set_ylabel('Pressure')
        axes[plot_idx].grid(alpha=0.3)
        n_last = int(len(data['Press']) * 0.3)
        if n_last > 10:
            last_mean = data['Press'][-n_last:].mean()
            last_std = data['Press'][-n_last:].std()
            axes[plot_idx].axhline(y=last_mean, color='r', linestyle='--',
                                   label=f'Last 30%: {last_mean:.3f} ± {last_std:.3f}')
            axes[plot_idx].legend()
        plot_idx += 1
    
    # Box Volume (normalized)
    if has_box:
        timesteps, box_dims = read_timestep_volume_file(box_vol_file)
        # Compute volume from dimensions (assuming format: timestep Lx Ly Lz)
        # Actually file has: timestep Lx Ly Lz
        box_vols = []
        with open(box_vol_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) == 4:
                        lx, ly, lz = float(parts[1]), float(parts[2]), float(parts[3])
                        box_vols.append(lx * ly * lz)
        box_vols = np.array(box_vols)
        
        if len(box_vols) > 0:
            vol_normalized = box_vols / box_vols[0]
            axes[plot_idx].plot(timesteps, vol_normalized, 'm-', linewidth=2.0)
            axes[plot_idx].set_ylabel('Box Volume / Initial')
            axes[plot_idx].grid(alpha=0.3)
            
            n_last = int(len(vol_normalized) * 0.3)
            if n_last > 10:
                last_mean = vol_normalized[-n_last:].mean()
                last_std = vol_normalized[-n_last:].std()
                axes[plot_idx].axhline(y=last_mean, color='r', linestyle='--',
                                       label=f'Last 30%: {last_mean:.3f} ± {last_std:.3f}')
                axes[plot_idx].legend()
        plot_idx += 1
    
    # Gel Volume - Bounding Box (normalized)
    if has_gel_bb:
        gel_bb_vols = read_volume_file(gel_bb_file)
        if len(gel_bb_vols) > 0:
            gel_bb_normalized = gel_bb_vols / gel_bb_vols[0]
            timesteps_gel = np.arange(len(gel_bb_vols))
            
            axes[plot_idx].plot(timesteps_gel, gel_bb_normalized, 'orange', linewidth=2.0)
            axes[plot_idx].set_ylabel('Gel Volume (BB) / Initial')
            axes[plot_idx].grid(alpha=0.3)
            
            n_last = int(len(gel_bb_normalized) * 0.3)
            if n_last > 10:
                last_mean = gel_bb_normalized[-n_last:].mean()
                last_std = gel_bb_normalized[-n_last:].std()
                axes[plot_idx].axhline(y=last_mean, color='r', linestyle='--',
                                       label=f'Last 30%: {last_mean:.3f} ± {last_std:.3f}')
                axes[plot_idx].legend()
        plot_idx += 1
    
    # Gel Volume - Radius of Gyration (normalized)
    if has_gel_rg:
        timesteps_rg, gel_rg_vols = read_timestep_volume_file(gel_rg_file)
        if len(gel_rg_vols) > 0:
            gel_rg_normalized = gel_rg_vols / gel_rg_vols[0]
            
            axes[plot_idx].plot(timesteps_rg, gel_rg_normalized, 'cyan', linewidth=2.0)
            axes[plot_idx].set_ylabel('Gel Volume (Rg³) / Initial')
            axes[plot_idx].grid(alpha=0.3)
            
            n_last = int(len(gel_rg_normalized) * 0.3)
            if n_last > 10:
                last_mean = gel_rg_normalized[-n_last:].mean()
                last_std = gel_rg_normalized[-n_last:].std()
                axes[plot_idx].axhline(y=last_mean, color='r', linestyle='--',
                                       label=f'Last 30%: {last_mean:.3f} ± {last_std:.3f}')
                axes[plot_idx].legend()
        plot_idx += 1
    
    axes[-1].set_xlabel('Step')
    
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Plot saved to {output}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python plot_lammps_log.py <folder> <dataname>")
        sys.exit(1)
    
    foldername = sys.argv[1]
    dataname = sys.argv[2]
    
    filepath = os.path.join(foldername, 'log.lammps')
    data = parse_lammps_log(filepath)
    
    if not data:
        print(f"No thermo data found in {filepath}")
        sys.exit(1)
    
    output = os.path.join(foldername, 'output_plots/convergence_plots', f'{dataname}_convergence.png')
    os.makedirs(os.path.join(foldername, 'output_plots/convergence_plots'), exist_ok=True)
    
    plot_convergence(data, foldername, dataname, output)


