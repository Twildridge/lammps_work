# LAMMPS Simulation Workflow on Bridges-2

This is the setup for running coarse-grained hydrogel simulations on Bridges-2. The system handles ~400k-1M beads and stores trajectory files in scratch to avoid filling up the home directory quota.

## Folder Organization

```
~/Documents/
├── lammps_work/                          # Git repository - sync this across machines
│   ├── scripts/
│   │   ├── run_lammps.sh                # Main runner - handles LAMMPS execution
│   │   ├── plot_lammps_log.py           # Plots convergence (T, P, volumes)
│   │   ├── plot_stress_profiles.py      # Stress and volume fraction profiles
│   │   └── write_tracking.py            # Performance tracking across runs
│   ├── simulations/
│   │   └── slab_with_support/
│   │       ├── slab_with_support.lmp    # LAMMPS input script
│   │       └── slab_with_support.batch  # SLURM batch script
│   ├── tracking.txt                     # Central performance log
│   └── README.md                        # This file
│
├── lammps_data/
│   └── input_data/                      # Initial configurations (.data files)
│
└── lammps_runs/                         # Timestamped working directories
    └── slab_*_20251219_143022/          # One per job submission
        ├── data_files/                  # Symlink to input data
        ├── output_files/                # Stress, volume, thermo data
        ├── output_plots/                # Convergence and profile plots
        ├── traj_files/                  # Symlink → scratch (see below)
        ├── log.lammps
        └── restart_*.restart            # For continuing runs

/ocean/projects/chm250028p/dpollard/
└── lammps_trajectories/                 # Multi-GB trajectory files live here
    └── slab_*_20251219_143022/
        └── polymer_*.lammpstrj
```

**Key point**: Trajectory files go to scratch (`/ocean/projects/...`) automatically. Everything else stays in home directory. This keeps you under the 10 GB quota.

## Running Simulations

### Fresh Run (starting from .data file)

Edit `slab_with_support.batch`:
```bash
NSTEPS=4000000
OLDSTEPS=0        # Zero means fresh run
TYPE=""           # Optional: stress, volume, or stressvol
```

Submit:
```bash
cd ~/Documents/lammps_work/simulations/slab_with_support
sbatch slab_with_support.batch
```

Check status:
```bash
squeue -u $USER
ssh r164    # Replace with your node number
htop        # Watch those cores work
```

### Continuation Run (restart from previous)

After your first run completes, check the restart file name:
```bash
ls ~/Documents/lammps_runs/slab_*_latest/restart*.restart
# Example output: restart_slab_support_5beads_..._1.5_1.4_4000000.restart
#                                                              ^^^^^^^ This is your OLDSTEPS
```

Edit batch script:
```bash
NSTEPS=4000000
OLDSTEPS=4000000  # Total timesteps from previous run
TYPE=""
```

Submit the same way. LAMMPS will automatically find the restart file and continue from there. The new run will output files named with total timesteps (8000000 in this case).

## What the Scripts Do

**run_lammps.sh**: Orchestrates everything. Creates working directory, symlinks data files, runs LAMMPS, calls post-processing scripts. You shouldn't need to edit this unless changing the workflow.

**plot_lammps_log.py**: Reads `log.lammps` and plots temperature, pressure, and volume convergence. Shows mean ± std for the last 30% of data so you can tell if things actually equilibrated.

**plot_stress_profiles.py**: Generates spatial profiles of stress and volume fraction along x, y, z axes. Needs the box dimensions from your .data file, so make sure it's accessible.

**write_tracking.py**: Logs performance data (atoms, runtime, timesteps) to a central tracking file and generates scaling plots. Useful for optimizing resource requests.

## Performance Notes

Optimal configuration on Bridges-2:
- **120 tasks per node** (not 128 - communication overhead kills you at 128)
- **1 CPU per task**
- **Pure MPI** (no hybrid MPI/OpenMP for this system size)

Expected runtime for 4M timesteps with 400k beads: ~20 hours on 120 cores.

## Common Issues

**"Disk quota exceeded"**: Trajectory files filled up home directory. Check that `traj_files/` is a symlink to scratch:
```bash
ls -la ~/Documents/lammps_runs/slab_*/traj_files
# Should show: traj_files -> /ocean/projects/...
```

**"Data file not found"**: Make sure your .data file is in `~/Documents/lammps_data/input_data/` with the exact name (no `_1.5_1.4_20000` suffix).

**Jobs stuck in CG state**: Contact Bridges-2 support. This is a system issue, not your fault.

**Simulation diverged**: Your gel probably exploded. Check that your interaction parameters make sense and you're not using too large a timestep. Also verify you're actually reading the correct data file.

## Notes

- Post-processing runs automatically after LAMMPS finishes (plots and tracking)
- Each run gets a timestamped directory so you don't accidentally overwrite data
- The tracking file aggregates performance across all runs for easy comparison
- Don't delete restart files if you plan to continue the run later

Now go equilibrate some hydrogels. May your simulations converge and your atoms stay bonded.