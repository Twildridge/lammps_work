#!/usr/bin/env python3
import sys
import os
import re
import glob
import shutil
import numpy as np
import matplotlib.pyplot as plt

# HOW TO RUN
# cd ~/Documents/lammps_runs/slab_with_support_*_<latest_timestamp>
# module load anaconda3
# python ~/Documents/lammps_work/scripts/write_tracking.py . "slab_support_5beads_10x10x5_rho6_extra_padding43_1.5_1.4_40000" "1"

def get_tracking_file_path():
    """Get path to central tracking file in lammps_work."""
    # Look for lammps_work directory
    home = os.path.expanduser('~')
    tracking_file = os.path.join(home, 'Documents', 'lammps_work', 'tracking.txt')
    return tracking_file

def parse_data_file(foldername, dataname, suffix=""):
    """Extract box dimensions and atom count from LAMMPS data file (excluding support atoms)."""
    # Extract base dataname without interaction and timesteps
    # Format: slab_support_5beads_10x10x5_rho6_extra_padding43_1.5_1.4_20000
    parts = dataname.split('_')
    
    # Find where interaction starts (format: number.number)
    base_parts = []
    for part in parts:
        if re.match(r'\d+\.\d+', part):  # Found interaction parameter
            break
        base_parts.append(part)
    
    base_name = '_'.join(base_parts)
    
    # Remove suffix digit if present
    if suffix and base_name.endswith(suffix):
        base_name = base_name[:-len(suffix)]
    
    # Look for data file in working directory
    data_file = os.path.join(foldername, 'data_files', f'{base_name}.data')
    
    if not os.path.exists(data_file):
        print(f"Data file not found: {data_file}")
        return None, None
    
    box_dims = {}
    natoms = 0
    num_support = 0
    
    reading_atoms = False
    with open(data_file, 'r') as f:
        for line in f:
            if 'atoms' in line and not reading_atoms:
                natoms = int(line.split()[0])
            elif 'xlo xhi' in line:
                vals = line.split()
                box_dims['x'] = float(vals[1]) - float(vals[0])
            elif 'ylo yhi' in line:
                vals = line.split()
                box_dims['y'] = float(vals[1]) - float(vals[0])
            elif 'zlo zhi' in line:
                vals = line.split()
                box_dims['z'] = float(vals[1]) - float(vals[0])
            elif line.strip() == 'Atoms':
                reading_atoms = True
            elif reading_atoms and line.strip() and not line.startswith('Bonds'):
                parts = line.split()
                if len(parts) >= 3:
                    atom_type = int(parts[2])
                    if atom_type in [4, 5]:  # Support and piston atoms
                        num_support += 1
    
    natoms_mobile = natoms - num_support
    print(f"Mobile atoms: {natoms_mobile}")
    return box_dims, natoms_mobile

def parse_lammps_log(filepath):
    """Extract wall time from LAMMPS log file."""
    wall_time = None
    with open(filepath, 'r') as f:
        for line in f:
            if 'Loop time of' in line:
                try:
                    wall_time = float(line.split()[3])
                except:
                    pass
    return wall_time

def parse_tracking_file(tracking_file):
    """Parse tracking.txt and extract all simulation data."""
    data = []
    
    if not os.path.exists(tracking_file):
        return data
    
    with open(tracking_file, 'r') as f:
        lines = f.readlines()[2:]  # Skip header
        
        for line in lines:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            
            sim_name = parts[0]
            natoms = int(parts[4])
            time_str = parts[5]
            
            # Parse time to seconds
            if ':' in time_str:
                time_parts = time_str.split(':')
                if len(time_parts) == 2:
                    time_sec = int(time_parts[0]) * 60 + int(time_parts[1])
                else:
                    time_sec = int(time_parts[0]) * 3600 + int(time_parts[1]) * 60 + int(time_parts[2])
            else:
                continue
            
            # Extract beads, padding iteration, and nsteps
            beads_match = re.search(r'(\d+)beads', sim_name)
            padding_match = re.search(r'extra_padding(\d+)', sim_name)
            nsteps_match = re.search(r'_(\d+)$', sim_name)
            
            beads = int(beads_match.group(1)) if beads_match else 0
            padding = int(padding_match.group(1)) if padding_match else 1
            nsteps = int(nsteps_match.group(1)) if nsteps_match else 0
            
            data.append({
                'name': sim_name,
                'natoms': natoms,
                'time_sec': time_sec,
                'beads': beads,
                'padding': padding,
                'nsteps': nsteps
            })
    
    return data

def write_tracking_file(dataname, box_dims, natoms, wall_time):
    """Write or append to central tracking file in lammps_work."""
    tracking_file = get_tracking_file_path()
    
    # Create directory if needed
    os.makedirs(os.path.dirname(tracking_file), exist_ok=True)
    
    # Backup existing tracking file
    if os.path.exists(tracking_file):
        shutil.copy(tracking_file, tracking_file.replace('.txt', '_backup.txt'))

    # Convert wall time to min:sec format
    if wall_time:
        minutes = int(wall_time // 60)
        seconds = int(wall_time % 60)
        time_str = f"{minutes}:{seconds:02d}"
    else:
        time_str = "N/A"
        wall_time = float('inf')
    
    # Read existing entries
    entries = []
    if os.path.exists(tracking_file):
        with open(tracking_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 2:
                for line in lines[2:]:
                    if line.strip():
                        parts = line.strip().split()
                        if len(parts) >= 5:
                            time_field = parts[-1]
                            if ':' in time_field:
                                time_parts = time_field.split(':')
                                if len(time_parts) == 2:
                                    time_seconds = int(time_parts[0]) * 60 + int(time_parts[1])
                                else:
                                    time_seconds = int(time_parts[0]) * 3600 + int(time_parts[1]) * 60 + int(time_parts[2])
                            elif time_field == 'N/A':
                                time_seconds = float('inf')
                            else:
                                time_seconds = float(time_field)
                            entries.append((line.strip(), time_seconds))
    
    # Check if entry already exists
    entry_exists = False
    for existing_entry, _ in entries:
        existing_name = existing_entry.split()[0]
        if existing_name == dataname:
            entry_exists = True
            print(f"Entry '{dataname}' already exists in tracking file, skipping.")
            break
    
    if not entry_exists:
        # Add new entry
        box_str = f"{box_dims.get('x', 0):<10.2f} {box_dims.get('y', 0):<10.2f} {box_dims.get('z', 0):<10.2f}"
        new_entry = f"{dataname:<70} {box_str} {natoms:<10} {time_str:<15}"
        entries.append((new_entry, wall_time if wall_time != float('inf') else float('inf')))
    
    # Sort by wall time
    entries.sort(key=lambda x: x[1])
    
    # Write sorted entries
    with open(tracking_file, 'w') as f:
        f.write(f"{'Simulation':<70} {'Box X':<10} {'Box Y':<10} {'Box Z':<10} {'Atoms':<10} {'Simulation Time':<15}\n")
        f.write("-" * 125 + "\n")
        for entry, _ in entries:
            f.write(entry + "\n")
    
    print(f"Tracking info written to {tracking_file}")

def plot_performance(data, output_dir):
    """Create performance plots."""
    if not data:
        return
    
    # Separate data by third digit of padding AND whether padding starts with 4 (piston)
    data_sim_only = []
    data_stress = []
    data_volume = []
    data_stressvol = []
    data_sim_only_piston = []
    data_stress_piston = []
    data_volume_piston = []
    data_stressvol_piston = []
    
    for d in data:
        padding_str = str(d['padding'])
        has_piston = padding_str.startswith('4')
        
        if len(padding_str) == 2:  # e.g., 31, 32, 41, 42
            if has_piston:
                data_sim_only_piston.append(d)
            else:
                data_sim_only.append(d)
        elif len(padding_str) >= 3:
            third_digit = int(padding_str[2])
            if third_digit == 1:
                if has_piston:
                    data_stress_piston.append(d)
                else:
                    data_stress.append(d)
            elif third_digit == 2:
                if has_piston:
                    data_volume_piston.append(d)
                else:
                    data_volume.append(d)
            elif third_digit == 3:
                if has_piston:
                    data_stressvol_piston.append(d)
                else:
                    data_stressvol.append(d)
    
    if not any([data_sim_only, data_stress, data_volume, data_stressvol,
                data_sim_only_piston, data_stress_piston, data_volume_piston, data_stressvol_piston]):
        print("No entries found, skipping plots")
        return
    
    # Plot 1: Time/Timestep vs Atoms (log-log)
    fig, ax = plt.subplots(figsize=(10, 8))

    types = [
        ('sim_only', data_sim_only, 'o', 'k', 'Simulation only'),
        ('stress', data_stress, 's', 'r', 'Sim. + stress'),
        ('volume', data_volume, '^', 'b', 'Sim. + volume'),
        ('stressvol', data_stressvol, 'D', 'purple', 'Sim. + stress/vol'),
        ('sim_only_piston', data_sim_only_piston, 'o', 'green', 'Simulation only (piston)'),
        ('stress_piston', data_stress_piston, 's', 'green', 'Sim. + stress (piston)'),
        ('volume_piston', data_volume_piston, '^', 'green', 'Sim. + volume (piston)'),
        ('stressvol_piston', data_stressvol_piston, 'D', 'green', 'Sim. + stress/vol (piston)')
    ]
    
    from matplotlib.lines import Line2D
    legend_elements = []
    
    for type_name, type_data, marker, color, label in types:
        if not type_data:
            continue
            
        all_atoms = []
        all_time = []
        
        unique_beads = sorted(set(d['beads'] for d in type_data))
        for beads in unique_beads:
            subset = [d for d in type_data if d['beads'] == beads]
            if subset:
                atoms = [d['natoms'] for d in subset]
                time_per_step = [d['time_sec'] / d['nsteps'] for d in subset]
                all_atoms.extend(atoms)
                all_time.extend(time_per_step)
                ax.scatter(atoms, time_per_step, color=color, marker=marker,
                          s=100, alpha=0.7)
        
        # Fit line - only if at least 3 unique atom counts
        unique_atom_counts = len(set(all_atoms))
        if unique_atom_counts >= 3:
            log_atoms = np.log10(all_atoms)
            log_time = np.log10(all_time)
            coeffs = np.polyfit(log_atoms, log_time, 1)
            slope = coeffs[0]
            
            atoms_fit = np.logspace(np.log10(min(all_atoms)), np.log10(max(all_atoms)), 100)
            time_fit = 10**(coeffs[0] * np.log10(atoms_fit) + coeffs[1])
            ax.plot(atoms_fit, time_fit, '--', color=color, linewidth=2, alpha=0.5)
            
            # Triangle
            x_tri = (min(all_atoms) * max(all_atoms)) ** 0.5
            y_tri = 10**(coeffs[0] * np.log10(x_tri) + coeffs[1])
            dx = x_tri * 0.2
            dy = y_tri * (10**(slope * np.log10(1.2)) - 1)
            
            ax.plot([x_tri, x_tri + dx], [y_tri, y_tri], '-', color=color, linewidth=1.5)
            ax.plot([x_tri + dx, x_tri + dx], [y_tri, y_tri + dy], '-', color=color, linewidth=1.5)
            ax.plot([x_tri, x_tri + dx], [y_tri, y_tri + dy], '-', color=color, linewidth=1.5)
            ax.text(x_tri + dx * 0.5, y_tri + dy * 1.3, f'{slope:.2f}', 
                    ha='center', fontsize=11, fontweight='bold', color=color)
        
        legend_elements.append(Line2D([0], [0], marker=marker, color='w', 
                                     markerfacecolor=color, markersize=10, label=label))
    
    ax.legend(handles=legend_elements, loc='upper left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Number of Atoms')
    ax.set_ylabel('Computation Time per Timestep (s/step)')
    ax.set_title('Computation Time per Timestep vs Number of Atoms (Log-Log)')
    ax.grid(alpha=0.3, which='both')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'time_vs_atoms.png'), dpi=150, bbox_inches='tight')
    print(f"Saved {os.path.join(output_dir, 'time_vs_atoms.png')}")
    plt.close()

    # Plot 2: Time vs Timesteps
    fig, ax = plt.subplots(figsize=(10, 8))
    
    unique_beads_all = sorted(set(d['beads'] for d in data))
    unique_padding_all = sorted(set(d['padding'] for d in data))
    colors_all = plt.cm.tab10(np.linspace(0, 1, len(unique_padding_all)))
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
    
    for i, padding in enumerate(unique_padding_all):
        for j, beads in enumerate(unique_beads_all):
            subset = [d for d in data if d['padding'] == padding and d['beads'] == beads]
            if subset:
                nsteps = [d['nsteps'] for d in subset]
                times = [d['time_sec'] / 60 for d in subset]
                ax.scatter(nsteps, times, color=colors_all[i], marker=markers[j % len(markers)],
                         s=100, alpha=0.7, label=f'{beads}beads, padding{padding}')
    
    ax.set_xlabel('Number of Timesteps')
    ax.set_ylabel('Computation Time (minutes)')
    ax.set_title('Computation Time vs Number of Timesteps')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'time_vs_timesteps.png'), dpi=150, bbox_inches='tight')
    print(f"Saved {os.path.join(output_dir, 'time_vs_timesteps.png')}")
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python write_tracking.py <folder> <dataname> [suffix]")
        sys.exit(1)
    
    foldername = sys.argv[1]
    dataname = sys.argv[2]
    suffix = sys.argv[3] if len(sys.argv) > 3 else ""
    
    # Parse data file
    box_dims, natoms = parse_data_file(foldername, dataname, suffix)
    
    # Parse log file
    logfile = os.path.join(foldername, 'log.lammps')
    wall_time = parse_lammps_log(logfile)
    
    if box_dims and natoms:
        write_tracking_file(dataname, box_dims, natoms, wall_time)
        
        # Generate performance plots in lammps_work directory
        tracking_file = get_tracking_file_path()
        output_dir = os.path.dirname(tracking_file)
        data = parse_tracking_file(tracking_file)
        plot_performance(data, output_dir)
    else:
        print("Error: Could not parse data file")