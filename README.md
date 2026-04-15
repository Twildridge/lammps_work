# LAMMPS Hydrogel Simulation — Project Guide

This repository contains everything needed to run, analyze, and reproduce coarse-grained molecular dynamics simulations of tetrahedral hydrogel slabs. The goals are to extract poroelastic observables — pore pressure, volume fraction profiles, cooperative diffusivity, and elastic moduli (M and G) — directly from MD trajectories.

> **New to this project?** Read sections 1–4 in order. Come back to the later sections as you need them.
> **Returning user?** Use the section headers to jump to what you need.

---

## Table of Contents

1. [What the simulations model](#1-what-the-simulations-model)
2. [Repository layout](#2-repository-layout)
3. [First-time setup](#3-first-time-setup)
4. [Building input data files](#4-building-input-data-files)
5. [Running simulations](#5-running-simulations)
6. [Downloading output from clusters](#6-downloading-output-from-clusters)
7. [Analysis](#7-analysis)
8. [Syncing with GitHub](#8-syncing-with-github)
9. [Cluster reference](#9-cluster-reference)
10. [Common issues](#10-common-issues)
11. [LJ units quick reference](#11-lj-units-quick-reference)

---

## 1. What the simulations model

The system is a **coarse-grained tetrahedral hydrogel** — a crosslinked polymer network swollen with solvent — represented using the Kremer–Grest bead-spring model in LAMMPS. Each bead represents roughly one Kuhn segment (~0.76 nm, ~93 g/mol for PEG). The simulation runs in reduced Lennard-Jones (LJ) units where length, energy, and mass are all set to 1.

**Physical analog:** PEG-based hydrogel membranes used in filtration, drug delivery, and soft robotics.

**Research goals (poroelasticity observables):** pore pressure and volume fraction profiles at equilibrium, longitudinal modulus M from polymer partial stress under axial compression, shear modulus G from box shearing, cooperative diffusivity Dc and solvent flux from permeation runs, and verification of the osmotic pressure balance ΔP_th = ΔP + ΔΠ.

The LJ units cheat sheet (`lj_units_cheat_sheet.md`) maps simulation values to physical quantities (MPa, nm, ps).

---

## 2. Repository layout

```
lammps_work/                    ← This git repository
│
├── simulations/                ← One folder per simulation type
│   ├── slab_with_support/      ← Gel equilibration and compression (main workhorse)
│   ├── slab_with_flow/         ← Permeation and compression with piston forcing
│   ├── shear_slab/             ← Shear modulus measurement (xz shear of isolated gel)
│   ├── solvent_phase/          ← Pure solvent equation of state sweep
│   ├── solvent_pure/           ← Single-state pure solvent run
│   ├── polymer_phase/          ← Pure polymer equation of state sweep
│   └── polymer_pure/           ← Single-state pure polymer run
│
├── scripts/                    ← Shell and Python scripts
│   ├── run_lammps.sh           ← Job runner for Expanse
│   ├── run_lammps_bridges.sh   ← Job runner for Bridges-2
│   ├── run_lammps_pod.sh       ← Job runner for Pod
│   ├── git_sync.sh             ← One-command GitHub sync (MacBook)
│   ├── add_walls_to_slab.ipynb         ← Build slab data file (no angles)
│   ├── add_walls_with_angles.ipynb     ← Build slab data file (with angles)
│   ├── slab_with_support.ipynb         ← Build basic slab geometry
│   ├── slab_with_support_angled.ipynb  ← Build angled-chain slab geometry
│   ├── isolate_gel.ipynb               ← Extract just the swollen gel from a run
│   ├── split_gel_slab.ipynb            ← Split a gel slab into pieces
│   ├── pure_polymer.ipynb              ← Build pure polymer data file
│   ├── pure_solvent_1.ipynb            ← Build pure solvent data file
│   ├── longitudinal_modulus_analysis.ipynb  ← Calculate M from compression run
│   ├── flow_poroelasticity_analysis.ipynb   ← Calculate Dc, flux, pore pressure
│   ├── plot_lammps_log.py      ← Plot T, P, volume convergence from log.lammps
│   ├── plot_stress_profiles.py ← Plot stress and volume fraction profiles
│   ├── plot_piston_data.py     ← Plot piston position and velocity
│   └── write_tracking.py       ← Log performance data to tracking.txt
│
├── lammps_data/                ← Small .data files committed to git
│   └── (input_data/ lives outside the repo — see §4)
│
├── README.md                   ← This file
├── CLAUDE_CONTEXT.md           ← Context file for Cowork AI sessions
├── lj_units_cheat_sheet.md     ← Unit conversions and parameter reference
├── expanse_lammps_guide.md     ← Expanse-specific setup and GitHub guide
├── slurm_commands_and_compiling.md  ← Cluster commands and LAMMPS build notes
├── slab_data_file_info.md      ← Log of every .data file that has been generated
├── documenting_pod_runs.md     ← HPC performance benchmarks
└── tracking.txt                ← Legacy performance log (not actively maintained)
```

**What lives outside the repo** (on each machine):

```
~/Documents/
├── lammps_work/        ← This repo (git-tracked)
├── lammps_data/
│   └── input_data/     ← .data files (large; NOT in git — copy manually)
└── lammps_runs/        ← Timestamped output from each job (NOT in git)
    └── slab_*_20250101_120000/
        ├── data_files/     ← symlink to input .data file
        ├── output_files/   ← stress, volume, piston .dat files
        ├── output_plots/   ← auto-generated plots
        ├── traj_files/     ← symlink to scratch (large trajectory files)
        └── log.lammps
```

---

## 3. First-time setup

### 3a. MacBook

**Clone the repository:**
```bash
cd ~/Documents
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
chmod +x scripts/*.sh
```

**Make the `lsync` shortcut available** (only needed once):
```bash
echo "alias lsync='bash ~/docs/grad_research/lammps/lammps_work/scripts/git_sync.sh'" >> ~/.zshrc
source ~/.zshrc
```

After this, typing `lsync` from any terminal will pull, commit, and push in one step. See [§8](#8-syncing-with-github) for details.

**Create the data file directory** (if it doesn't exist):
```bash
mkdir -p ~/Documents/lammps_data/input_data
```

---

### 3b. On a cluster (first time)

> **What is a cluster?** A cluster is a collection of powerful computers ("nodes") you access remotely. You write a script describing your job (how many CPUs, how long), submit it to a queue manager called SLURM, and SLURM runs it when resources are available. You don't interact with the job while it runs.

SSH into the cluster, then clone the repo:

```bash
# Example for Expanse:
ssh <username>@login.expanse.sdsc.edu

cd ~/Documents
mkdir -p lammps_data/input_data    # for .data files you'll scp over
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
chmod +x scripts/*.sh
```

**Configure git credentials** (so `git pull` works without a password prompt):
```bash
git config --global user.name "Twildridge"
git config --global user.email "semiinfiniteslab@icloud.com"
git config --global credential.helper store
# Then do: git pull  — it'll ask for your GitHub PAT once, then store it
```

> A **Personal Access Token (PAT)** replaces your GitHub password for command-line use. Generate one at github.com → Settings → Developer settings → Personal access tokens → Tokens (classic). Check the `repo` scope, copy it, and paste it when prompted by `git pull`.

**Copy your .data files to the cluster:**
```bash
# From MacBook terminal:
scp ~/Documents/lammps_data/input_data/your_file.data \
    <username>@data.expanse.sdsc.edu:~/Documents/lammps_data/input_data/
```

Full details for each cluster are in `expanse_lammps_guide.md` and `slurm_commands_and_compiling.md`.

---

## 4. Building input data files

Before running a gel simulation you need a `.data` file — a text file describing every atom's position, bond topology, and atom types. These are generated on your MacBook using Jupyter notebooks in `scripts/`.

> **What is a Jupyter notebook?** An interactive Python document where you run code cells one at a time and see results immediately. Open JupyterLab from your terminal with `jupyter lab`, then navigate to the file.

### Which notebook to use

| Notebook | What it builds | When to use |
|----------|---------------|-------------|
| `add_walls_to_slab.ipynb` | Gel slab + flat support + piston (no chain angles) | Standard compression/flow runs |
| `add_walls_with_angles.ipynb` | Same but with harmonic and cosine angle terms | When running angle-restrained chains |
| `slab_with_support.ipynb` | Basic slab geometry builder | Reference / older geometry |
| `slab_with_support_angled.ipynb` | Angled-chain slab geometry | Angled geometry variants |
| `isolate_gel.ipynb` | Extracts the swollen polymer (+solvent) from a finished slab run | Input for `shear_slab` and modulus analysis |
| `split_gel_slab.ipynb` | Splits a slab into polymer-only and solvent-only files | Isolated component analysis |
| `pure_polymer.ipynb` | Pure polymer box (no solvent) | EOS and baseline runs |
| `pure_solvent_1.ipynb` | Pure solvent box | EOS and solvent calibration |

### Typical workflow for a new slab

1. Open `add_walls_to_slab.ipynb`
2. Set the unit cell dimensions (e.g. `10×10×8`), beads per chain, solvent density
3. Run all cells → generates a `.data` file in `../lammps_data_files_local/`
4. The output filename encodes all key parameters (e.g. `slab_support_5beads_tall_3.data`)
5. Copy to `~/Documents/lammps_data/input_data/` and then to the cluster

Generated file specs are logged in `slab_data_file_info.md` so you can always look up what was built and when.

---

## 5. Running simulations

### 5a. How it works

When you submit a job, the flow is:

```
You edit .batch file → sbatch → SLURM queues job → node allocated →
run_lammps[_bridges].sh runs → creates timestamped lammps_runs/ dir →
pulls latest scripts from GitHub → runs LAMMPS → post-processing plots
```

The `run_lammps` scripts handle all the bookkeeping automatically: creating timestamped working directories, symlinking the data file, routing trajectory files to scratch storage, and running post-processing after LAMMPS finishes. You only need to edit the `.batch` file.

### 5b. Simulation types

| Folder | What it simulates | Typical use |
|--------|------------------|-------------|
| `slab_with_support/` | Gel equilibration or axial compression with piston | Equilibrate gel; measure M (longitudinal modulus) |
| `slab_with_flow/` | Piston-forced permeation or compression of a pre-equilibrated walled gel | Measure flux, Dc, pore pressure profiles |
| `shear_slab/` | xz shear of an isolated swollen gel (from `isolate_gel.ipynb`) | Measure G (shear modulus) |
| `solvent_phase/` | Pure solvent pressure sweep across many state points | Build solvent EOS |
| `solvent_pure/` | Single pure solvent run | Baseline pressure/density check |
| `polymer_phase/` | Pure polymer pressure sweep | Build polymer EOS |
| `polymer_pure/` | Single pure polymer run | Baseline |

Each folder contains:
- `<name>.lmp` — the LAMMPS input script (controls the physics; rarely needs editing)
- `<name>.batch` — the SLURM job script (controls cluster resources; **you edit this**)
- `<name>_bridges.batch` — Bridges-2 version (if applicable)
- `<name>_pod.batch` — Pod version (if applicable)

### 5c. Editing the batch file

Open the `.batch` file for your target cluster. The section you edit is at the bottom:

```bash
FOLDER="slab_with_flow"
DATANAME="walled_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000"
INTERACTION="1.0_1.0"   # epsSS_epsSP
NSTEPS=2000000
OLDSTEPS=0              # 0 = fresh run; set to previous total steps for continuation
TYPE=""                 # stress, volume, stressvol, or leave empty
```

- **`DATANAME`**: the filename of your `.data` file in `lammps_data/input_data/`, without the `.data` extension
- **`INTERACTION`**: `epsSS_epsSP` — solvent–solvent and polymer–solvent LJ well depths (see cheat sheet)
- **`NSTEPS`**: number of timesteps to run (1 million steps ≈ 4100 τ ≈ 19 ns for PEG)
- **`OLDSTEPS`**: set to the total timesteps already completed if continuing a previous run; 0 otherwise
- **`TYPE`**: optional suffix that gets appended to the data file name in output (for labelling stress/volume variants)

For `shear_slab`, `NSTEPS` controls only Phase 3 (production). Phase 1 (NPT equilibration, 50k steps) and Phase 2 (shear application, 200k steps) are hardcoded in `shear_slab.lmp`.

### 5d. Submitting a job

```bash
# SSH into the cluster, then:
cd ~/Documents/lammps_work/simulations/slab_with_flow
sbatch slab_with_flow_bridges.batch    # Bridges-2
sbatch slab_with_flow.batch            # Expanse
```

> **What is `sbatch`?** It submits a job script to SLURM's queue. SLURM schedules it when the requested nodes are free. You get a job ID back immediately; the simulation runs in the background.

**Check your job status:**
```bash
squeue -u $USER                        # see all your queued/running jobs
watch -n 5 squeue -u $USER            # auto-refresh every 5 seconds (Ctrl-C to stop)
```

**Monitor a running job:**
```bash
ssh <node_name>     # e.g. ssh r164 (Bridges-2) — node name is in squeue output
htop                # live CPU/memory usage per core
```

**Read the LAMMPS log while it runs:**
```bash
# From the timestamped working directory (check run_lammps output for the exact path):
tail -f ~/Documents/lammps_runs/slab_with_flow_*/log.lammps
```

### 5e. Continuing a run from a restart file

LAMMPS writes `.restart` files at regular intervals so you can pick up where you left off if a job times out or you need more steps.

1. Find the restart file in the run's working directory:
   ```bash
   ls ~/Documents/lammps_runs/slab_*_latest/restart*.restart
   # The number at the end is the timestep: e.g. restart_..._4000000.restart → OLDSTEPS=4000000
   ```

2. Edit the batch file:
   ```bash
   NSTEPS=4000000     # additional steps you want to run
   OLDSTEPS=4000000   # total steps completed so far
   ```

3. Submit again with `sbatch`. LAMMPS will find the restart file automatically.

### 5f. SLURM resource guidelines

| Cluster | Partition (CPU) | Cores/node | Optimal ntasks/node | Max walltime |
|---------|----------------|------------|---------------------|--------------|
| Bridges-2 | `RM` | 128 | **120** (128 causes comm overhead) | 5 days |
| Expanse | `compute` | 128 | 128 | 2 days |
| Pod | varies | 40–80 | 40 | check queue |

**Expected runtimes** (4M timesteps, 400k–600k atoms):
- 1 node (120 cores, Bridges-2): ~10–11 hours
- 4 nodes (480 cores, Bridges-2): ~3 hours
- 2 nodes (256 cores, Expanse): ~5–6 hours

Benchmarks across system sizes are in `documenting_pod_runs.md`.

---

## 6. Downloading output from clusters

After a job finishes, the output `.dat` files are in the working directory on the cluster. Trajectory files (`.lammpstrj`) are stored separately in scratch storage (the `traj_files/` symlink points there).

**Download output data files to MacBook** (into `flow_data_local/`):
```bash
# From MacBook terminal — adjust path to match your run
rsync -avP <username>@login.expanse.sdsc.edu:\
  ~/Documents/lammps_runs/slab_with_flow_<dataname>_<interaction>_<timestamp>/output_files/ \
  ~/Documents/lammps/flow_data_local/
```

**Download a trajectory file** (large — use the data transfer node):
```bash
rsync -avPz <username>@data.bridges2.psc.edu:\
  /ocean/projects/chm250028p/dpollard/lammps_trajectories/<run_folder>/<file>.lammpstrj.gz \
  ~/Downloads/
# Drop the 'z' flag if the file is already gzip-compressed
```

> **Trajectory files are large** (often 1–50 GB). Only download if you need to analyze atom positions directly (e.g. volume fraction profiles). The `.dat` output files are much smaller and contain pre-computed quantities.

---

## 7. Analysis

### 7a. Python scripts (auto-run after each job)

These scripts are called automatically by `run_lammps.sh` after LAMMPS finishes. You can also run them manually from any `lammps_runs/` working directory.

**`plot_lammps_log.py`**
Reads `log.lammps` and plots temperature, pressure, and volume vs. timestep. Shows mean ± std for the last 30% of the run to verify equilibration.
```bash
python ~/Documents/lammps_work/scripts/plot_lammps_log.py . "<dataname>_<interaction>_<totsteps>"
```

**`plot_stress_profiles.py`**
Reads the binned stress `.dat` files and plots partial stress and volume fraction profiles along x, y, and z.
```bash
python ~/Documents/lammps_work/scripts/plot_stress_profiles.py . "<dataname>_<interaction>_<totsteps>" <oldsteps>
```

**`plot_piston_data.py`**
Plots piston position and velocity vs. time. Useful for checking that the piston is moving at the right rate during compression/permeation runs.

**`write_tracking.py`**
Appends a row to `tracking.txt` with atoms, runtime, and timesteps, and regenerates the performance scaling plot. This was useful for benchmarking but is no longer actively maintained.

---

### 7b. Jupyter analysis notebooks

Open these on your MacBook in JupyterLab (`jupyter lab`), pointing them at data files in `flow_data_local/`.

**`longitudinal_modulus_analysis.ipynb`** ← *Current: Milestone 3 (done)*
Reads polymer partial stress vs. strain from a `slab_with_flow` compression run and extracts the longitudinal modulus M = Δσ_zz / Δε_zz. Inputs: `stress_tensor_polymer_*.dat`, `shear_strain_*.dat`.

**`flow_poroelasticity_analysis.ipynb`** ← *Current: Milestones 5–6*
Reads volume fraction and stress profiles from permeation runs. Extracts pore pressure profiles, cooperative diffusivity Dc (from mode-decay fitting), solvent flux, and osmotic/elastic pressure contributions. Inputs: `stress_profile_z_*.dat`, traj files.

**`shear_modulus_analysis.ipynb`** ← *Upcoming: Milestone 4*
Will read the box-integrated polymer stress tensor from `shear_slab` runs and extract shear modulus G = σ_xz / γ_xz. Inputs: `stress_tensor_polymer_*.dat`, `shear_strain_*.dat`.

---

## 8. Syncing with GitHub

The repo lives at `https://github.com/Twildridge/lammps_work`. Changes flow like this:

```
MacBook (edit scripts/notebooks)
    ↓  lsync
GitHub
    ↓  git pull (automatic at start of each cluster job)
Expanse / Bridges-2 / Pod
```

### MacBook → GitHub

```bash
lsync                        # auto-generates commit message from changed files
lsync "my custom message"    # use a specific commit message
```

`lsync` does: `git pull --rebase` → `git add -A` → `git commit` → `git push`. If there are no changes, it exits cleanly.

### Clusters → GitHub

If you edit scripts directly on a cluster, you need `lsync` there too. Add the alias once to each cluster's shell config:

```bash
# Run once on each cluster (Expanse, Bridges-2, Pod):
echo "alias lsync='bash ~/Documents/lammps_work/scripts/git_sync.sh'" >> ~/.bashrc
source ~/.bashrc
```

After that, `lsync` works identically on the cluster — pull, commit, push. Without this setup, edits made on a cluster will not be pushed to GitHub automatically.

Clusters also pull automatically at the start of every job (built into `run_lammps.sh` and `run_lammps_bridges.sh`). To pull manually at any time:
```bash
cd ~/Documents/lammps_work && git pull
```

### Handling merge conflicts

If you edited the same file on two machines without syncing between them, `git pull --rebase` will stop with a conflict error. Fix it like this:

```bash
git status                  # identify the conflicting file(s)
# Open each conflicting file — look for markers like:
#   <<<<<<< HEAD
#   (your local version)
#   =======
#   (incoming version from GitHub)
#   >>>>>>> origin/main
# Edit the file to keep what you want, removing all the marker lines.
git add <filename>
git rebase --continue
git push
```

The simplest way to avoid conflicts: always run `lsync` (or `git pull`) before you start editing on any machine.

---

### SSH key setup (if `git push` gives a login error)

This is needed once per cluster the first time you try to push. The error usually looks like `Permission denied (publickey)` or `fatal: Authentication failed`.

**1. Generate an SSH key on the cluster:**
```bash
ssh-keygen -t ed25519 -C "pollard@ucsb.edu"
# Press Enter for all prompts (default location, no passphrase needed)
```

**2. Print your public key:**
```bash
cat ~/.ssh/id_ed25519.pub
# Copy the entire line of output
```

**3. Add the key to GitHub:**
- Go to [github.com/settings/ssh/new](https://github.com/settings/ssh/new)
- Title: something descriptive like `"Bridges-2 cluster"` or `"Expanse cluster"`
- Key: paste the full line from step 2
- Click **Add SSH key**

**4. Test the connection:**
```bash
ssh -T git@github.com
# Expected: "Hi Twildridge! You've successfully authenticated..."
```

**5. Switch your repo remote to SSH:**
```bash
git remote set-url origin git@github.com:Twildridge/lammps_work.git
```

**6. Try pushing again:**
```bash
git push
```

Repeat steps 1–6 for each cluster (Expanse, Bridges-2, Pod) — each needs its own key added to GitHub. On your MacBook, this is usually handled automatically by the macOS keychain, so you shouldn't need to do it there.

---

### What is and isn't tracked by git

The `.gitignore` excludes large files so the repo stays fast:

| Excluded | Why |
|----------|-----|
| `*.lammpstrj`, `*.restart*` | Trajectory/restart files can be hundreds of GB |
| `output_files/`, `traj_files/` | Per-run outputs — too large and numerous |
| `log.lammps`, `*.log` | Redundant with run directories |
| `final_*.data` | LAMMPS-written final configs |
| `slurm_*.out`, `*.o*.*`, `*.e*.*` | SLURM stdout/stderr logs |
| `.DS_Store`, `__pycache__/`, `.ipynb_checkpoints/` | OS and IDE clutter |

Everything in `simulations/`, `scripts/`, and the docs is tracked. `.data` input files are **not** tracked (too large) — manage them manually with `scp`/`rsync`.

---

## 9. Cluster reference

### Expanse (SDSC) — `run_lammps.sh`

| Setting | Value |
|---------|-------|
| Login node | `login.expanse.sdsc.edu` |
| Data transfer node | `data.expanse.sdsc.edu` |
| CPU partition | `compute` |
| GPU partition | `gpu` |
| Account flag | Required: `--account=csb197` |
| Cores/node | 128 |
| Max walltime (compute) | 48 hours |
| Scratch | `/expanse/lustre/scratch/$USER/temp_project` |
| LAMMPS binary | loaded via `module load LAMMPS/...` then `lmp` |

Load modules (in this order):
```bash
module reset
module load gcc/10.2.0
module load openmpi/4.1.3
module load python/3.8.12
```

---

### Bridges-2 (PSC) — `run_lammps_bridges.sh`

| Setting | Value |
|---------|-------|
| Login node | `bridges2.psc.edu` |
| Data transfer node | `data.bridges2.psc.edu` |
| CPU partition | `RM` |
| GPU partition | `GPU` |
| Cores/node | 128 (use 120 for optimal MPI performance) |
| Max walltime | 5 days |
| Scratch | `/ocean/projects/chm250028p/dpollard/` |
| Home quota | 10 GB — trajectory files must go to scratch |
| LAMMPS binary | `/opt/packages/LAMMPS/lammps-22Jul2025/build-RM-gcc13.3.1/lmp` |

Load modules (in this order):
```bash
module purge
module load python/3.8.6
module load openmpi/5.0.8-gcc13.3.1
module load cuda/12.6.1
module load intel-mkl/2023.2.0
module load LAMMPS/22Jul25-gcc
```

> **Bridges-2 quota:** If you see "Disk quota exceeded", trajectory files have filled your 10 GB home directory. Check that `traj_files/` is a symlink pointing to scratch: `ls -la ~/Documents/lammps_runs/*/traj_files`

---

### Pod (CNSI, UCSB) — `run_lammps_pod.sh`

| Setting | Value |
|---------|-------|
| Login node | `pod-login1.cnsi.ucsb.edu` |
| Cores/node | 40 (optimal) |
| Internode bandwidth | 100 Gb/s (half of Bridges-2's 200 Gb/s — avoid multi-node for large systems) |
| Note | GPUs available (L40S); check `documenting_pod_runs.md` for benchmarks |

> **Pod requires the UCSB campus VPN.** You must be connected before you can SSH in. Download and install Ivanti Secure Access from [it.ucsb.edu/ivanti-secure-access-campus-vpn/get-connected-campus-vpn](https://it.ucsb.edu/ivanti-secure-access-campus-vpn/get-connected-campus-vpn), then connect to the UCSB VPN before running `ssh pod-login1.cnsi.ucsb.edu`. Expanse and Bridges-2 do not require a VPN.

---

## 10. Common issues

**"Data file not found"**
The `.data` file must be in `~/Documents/lammps_data/input_data/` with the exact name matching `DATANAME` in the batch file (no extra suffixes).

**"Disk quota exceeded" on Bridges-2**
Trajectory files ended up in your home directory. Verify the symlink: `ls -la ~/Documents/lammps_runs/*/traj_files` should show `→ /ocean/projects/...`. If not, move the `.lammpstrj` files to scratch manually.

**Job stuck in `CG` state on Bridges-2**
This is a cluster-side issue. Contact PSC support.

**Simulation diverged ("gel evaporated")**
Check that `epsSP` and `epsSS` are physically reasonable. Verify you're reading the correct data file. Ensure `timestep = 0.005` (larger timesteps can cause FENE bond divergence).

**`xztilt` variable parse error in `shear_slab.lmp`**
If LAMMPS reports an unknown keyword, replace `xztilt` with `xz` in the `variable shear_strain_current` definition (older LAMMPS versions use `xz` rather than `xztilt` in variable expressions).

**`fix ave/chunk norm none` error**
`norm none` requires LAMMPS ≥ March 2020. On older builds, replace the `prof_z_polymer/solvent` fixes with the reduce/chunk + ave/time approach used in `slab_with_flow.lmp`.

**`git pull` fails on cluster ("merge conflict")**
Run `git stash` to set aside local changes, then `git pull`, then `git stash pop` to restore them. If conflicts persist, resolve them manually or ask for help.

---

## 11. LJ units quick reference

Full details are in `lj_units_cheat_sheet.md`. Key conversions for PEG/water:

| Quantity | Multiply LJ value by | To get |
|----------|---------------------|--------|
| Length | 0.76 | nm |
| Time (τ) | 4.6 | ps |
| Pressure (P*) | 9.4 | MPa |
| Force (F*) | 5.4 | pN |
| Temperature (T*) | 300 | K (at T*=1) |
| 10⁶ steps at dt=0.005 | → | ~23 ns |

**Current simulation parameters:**
- Temperature: T* = 1.0 (= 300 K)
- Pressure (NPT target): P* = 1.5
- Timestep: dt = 0.005 τ
- FENE bonds: K=30, R₀=1.5, ε=1, σ=1 (Kremer–Grest standard)
- Pair interactions: WCA (rc = 1.122σ, purely repulsive) — consistent across all active simulation types

---

*Last updated: 2026-04-15. For questions, contact Dylan Pollard (pollard@ucsb.edu).*
