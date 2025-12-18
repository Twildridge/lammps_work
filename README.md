# LAMMPS Workflow Setup Guide

## Folder Structure

### On MacBook (in Documents/)
```
Documents/
├── lammps_work/                    # Git repository
│   ├── scripts/
│   │   ├── run_lammps.sh
│   │   ├── plot_lammps_log.py
│   │   ├── plot_stress_profiles.py
│   │   └── write_tracking.py
│   ├── simulations/
│   │   ├── slab_with_support/
│   │   │   ├── slab_with_support.lmp
│   │   │   └── slab_with_support.batch
│   │   └── slab_with_flow/
│   │       └── slab_with_flow.lmp
│   └── README.md
│
└── lammps_data/                    # NOT in Git
    └── input_data/
        └── (all your .data files)
```

### On Pod/Bridges-2 (in ~/Documents/)
```
Documents/
├── lammps_work/                    # Cloned from GitHub
│   ├── scripts/
│   └── simulations/
│
├── lammps_data/                    # Upload specific .data files here
│   └── input_data/
│
└── lammps_runs/                    # Auto-created working directories
    └── (timestamped run directories)
```

## Setup Instructions

### 1. Local Setup (MacBook)


# Make folder easier to get to—call it "docs"
nano ~/.zshrc
# add this line below to this zsh file
ln -s "$HOME/Library/Mobile Documents/com~apple~CloudDocs/Documents - Dylan's MacBook Air" ~/docs
# exit file and run below line
source ~/.zshrc


```bash
cd ~/docs
mkdir -p lammps_work/scripts lammps_work/simulations
mkdir -p lammps_data/input_data

# Initialize Git repo
cd lammps_work
git init
git remote add origin <your-github-repo-url>


# Make scripts executable
chmod +x scripts/run_lammps.sh

# Commit and push
git add .
git commit -m "Initial commit"
git push -u origin main
```

### 2. Cluster Setup (Pod/Bridges-2)
```bash
# SSH to cluster
ssh pod-login1.cnsi.ucsb.edu

cd ~/Documents
mkdir -p lammps_data/input_data

# Clone repository
git clone <your-github-repo-url>
cd lammps_work
chmod +x scripts/run_lammps.sh

# Upload specific .data files you need
# Use scp or rsync from your Mac:
# scp ~/Documents/lammps_data/input_data/slab_*.data user@pod:~/Documents/lammps_data/input_data/
```

### 3. Running Simulations

#### Submit a job:
```bash
cd ~/Documents/lammps_work/simulations/slab_with_support
sbatch slab_with_support.batch
```

#### Monitor:
```bash
squeue -u $USER
tail -f slurm_<jobid>.out
```

## Key Changes from Old Setup

1. **Data files separated**: Input .data files live in `lammps_data/`, not in Git
2. **Working directories**: Each run creates a timestamped directory in `lammps_runs/`
3. **Output organized**: All outputs (stress, volume, trajectories, plots) go to working directory
4. **Scripts version controlled**: All .lmp, .sh, .py, .batch files tracked in Git
5. **Synced across machines**: Pull latest scripts on clusters with `git pull`

## Workflow Example

1. Edit scripts locally on Mac (in VS Code)
2. Commit and push changes: `git commit -am "Updated parameters" && git push`
3. SSH to cluster via VS Code
4. Pull changes: `cd ~/Documents/lammps_work && git pull`
5. Submit job: `cd simulations/slab_with_support && sbatch slab_with_support.batch`
6. Results appear in `~/Documents/lammps_runs/<timestamped_directory>/`

## Batch Script Arguments

Edit the variables in `slab_with_support.batch`:
- `FOLDER`: simulation folder name (e.g., "slab_with_support")
- `DATANAME`: .data file name without extension
- `INTERACTION`: "epsSS_epsSP" format (e.g., "1.5_1.6")
- `NSTEPS`: number of timesteps to run
- `OLDSTEPS`: previous timesteps (for restart, default: NSTEPS)
- `DATASTEPS`: timesteps in data file (default: NSTEPS)
- `TYPE`: "stress", "volume", "stressvol", or empty
