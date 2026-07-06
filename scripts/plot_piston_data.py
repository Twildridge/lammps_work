#!/usr/bin/env python3
"""
Plot piston position and velocity versus timesteps.
Reads piston position and velocity from piston_data folder.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


def read_piston_file(filepath):
    """Read LAMMPS fix print output file with format: timestep value"""
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
                    timesteps.append(int(parts[0]))
                    values.append(float(parts[1]))
                except ValueError:
                    continue
    
    return np.array(timesteps), np.array(values)


def plot_piston_data(folder, dataname, oldsteps=0):
    """Plot piston position and velocity versus timesteps."""
    
    # Input file paths
    piston_data_dir = os.path.join(folder, 'output_files', 'piston_data')
    pos_file = os.path.join(piston_data_dir, f'piston_position_{dataname}.dat')
    vel_file = os.path.join(piston_data_dir, f'piston_velocity_{dataname}.dat')
    
    # Check files exist. Some geometries (e.g. slab_with_support, compress_slab)
    # never write piston position/velocity files, so a missing file is expected,
    # not an error — skip cleanly instead of printing a scary "Error:" line.
    if not os.path.exists(pos_file):
        print(f"No piston position file ({os.path.basename(pos_file)}) — skipping piston plots.")
        return True
    if not os.path.exists(vel_file):
        print(f"No piston velocity file ({os.path.basename(vel_file)}) — skipping piston plots.")
        return True
    
    # Read data
    timesteps_pos, positions = read_piston_file(pos_file)
    timesteps_vel, velocities = read_piston_file(vel_file)
    
    if len(timesteps_pos) == 0:
        print(f"Error: No data found in {pos_file}")
        return False
    
    print(f"Read {len(timesteps_pos)} position points, {len(timesteps_vel)} velocity points")
    
    # Create output directory for plots
    output_plot_dir = os.path.join(folder, 'output_plots')
    os.makedirs(output_plot_dir, exist_ok=True)
    
    # Create figure with two subplots
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Title
    if oldsteps > 0:
        fig.suptitle(f'{dataname}\n(continuing from {oldsteps} steps)', 
                     fontsize=12, fontweight='bold')
    else:
        fig.suptitle(f'{dataname}', fontsize=12, fontweight='bold')
    
    # Color scheme matching other scripts
    color_pos = plt.cm.viridis(0.3)
    color_vel = plt.cm.viridis(0.7)
    
    # Plot position
    axes[0].plot(timesteps_pos, positions, linewidth=1.5, color=color_pos, alpha=0.9)
    axes[0].set_ylabel('Piston Position (z)', fontsize=11)
    axes[0].grid(alpha=0.3)
    axes[0].set_title('Position vs Time', fontweight='bold')
    
    # Add horizontal line at initial position for reference
    if len(positions) > 0:
        axes[0].axhline(y=positions[0], color='gray', linestyle='--', 
                        alpha=0.5, label=f'Initial: {positions[0]:.2f}')
        axes[0].legend(loc='best', fontsize=9)
    
    # Plot velocity
    axes[1].plot(timesteps_vel, velocities, linewidth=1.5, color=color_vel, alpha=0.9)
    axes[1].set_ylabel('Piston Velocity (vz)', fontsize=11)
    axes[1].set_xlabel('Timestep', fontsize=11)
    axes[1].grid(alpha=0.3)
    axes[1].set_title('Velocity vs Time', fontweight='bold')
    
    # Add horizontal line at zero velocity
    axes[1].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    # Add horizontal line at mean velocity (excluding initial transient)
    if len(velocities) > 10:
        # Use last 50% of data for steady-state estimate
        steady_start = len(velocities) // 2
        mean_vel = np.mean(velocities[steady_start:])
        std_vel = np.std(velocities[steady_start:])
        axes[1].axhline(y=mean_vel, color='red', linestyle='-', alpha=0.7,
                        label=f'Steady-state mean: {mean_vel:.2e} ± {std_vel:.2e}')
        axes[1].legend(loc='best', fontsize=9)
    
    plt.tight_layout()
    
    # Save plot
    output_plot_file = os.path.join(output_plot_dir, f'{dataname}_piston.png')
    plt.savefig(output_plot_file, dpi=150)
    print(f"Piston plot saved to {output_plot_file}")
    plt.close()
    
    # Print summary statistics
    print("\n=== Piston Summary ===")
    print(f"Initial position: {positions[0]:.4f}")
    print(f"Final position:   {positions[-1]:.4f}")
    print(f"Total displacement: {positions[-1] - positions[0]:.4f}")
    if len(velocities) > 10:
        steady_start = len(velocities) // 2
        print(f"Steady-state velocity: {np.mean(velocities[steady_start:]):.2e} ± {np.std(velocities[steady_start:]):.2e}")
    
    return True


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


