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
│   ├── slab_with_flow/         ← Permeation OR compression with piston forcing (compression_mode switch)
│   ├── shear_slab/             ← Shear modulus measurement (plate-driven xz shear; current G workflow)
│   ├── slab_elongation/        ← LEGACY uniaxial elongation prototype — superseded by shear_slab
│   ├── solvent_phase/          ← Pure solvent equation of state sweep
│   ├── solvent_pure/           ← Single-state pure solvent run
│   ├── polymer_phase/          ← Pure polymer equation of state sweep
│   └── polymer_pure/           ← Single-state pure polymer run
│
├── scripts/                    ← Shell and Python scripts
│   ├── run_lammps.sh           ← Job runner for Expanse
│   ├── run_lammps_bridges.sh   ← Job runner for Bridges-2
│   ├── run_lammps_pod.sh       ← Job runner for Pod
│   ├── run_lammps_pod_old.sh   ← Legacy Pod runner (kept for reference)
│   ├── build_lammps.sh         ← One-shot LAMMPS build script (Expanse-style cmake)
│   ├── git_sync.sh             ← One-command GitHub sync (MacBook + clusters)
│   ├── add_walls_to_slab.ipynb         ← Build slab data file (no angles)
│   ├── add_walls_with_angles.ipynb     ← Build slab data file (with angles)
│   ├── slab_with_support.ipynb         ← Build basic slab geometry
│   ├── slab_with_support_angled.ipynb  ← Build angled-chain slab geometry
│   ├── slab_with_support_old.ipynb     ← Legacy slab builder (kept for reference)
│   ├── slab_with_support_old_2.ipynb   ← Legacy slab builder (kept for reference)
│   ├── isolate_gel.ipynb               ← Extract just the swollen gel from a run
│   ├── split_gel_slab.ipynb            ← Split a gel slab into pieces
│   ├── add_plates_to_gel.ipynb         ← Attach shear plates to isolated gel (input for shear_slab)
│   ├── add_more_plates_to_gel.ipynb    ← Variant: plates on all six faces (experimental)
│   ├── pure_polymer.ipynb              ← Build pure polymer data file
│   ├── pure_solvent_1.ipynb            ← Build pure solvent data file
│   ├── compression_analysis.ipynb           ← Current: M, Dc, pore pressure, volume fractions (compression runs)
│   ├── longitudinal_modulus_analysis.ipynb  ← Older M-only notebook (superseded by compression_analysis)
│   ├── permeation_analysis.ipynb            ← Flow profiles & pore pressure evolution (compression_mode=0)
│   ├── flow_poroelasticity_analysis.ipynb   ← Older flow notebook (superseded by permeation_analysis)
│   ├── shear_analysis.ipynb                 ← G, N1/N2, stress profiles (shear_slab output)
│   ├── compression_analysis_backup.ipynb    ← Frozen snapshot of compression notebook
│   ├── plot_lammps_log.py      ← Plot T, P, volume convergence from log.lammps (+ shear diagnostics)
│   ├── plot_stress_profiles.py ← Plot stress and volume fraction profiles
│   ├── plot_piston_data.py     ← Plot piston position and velocity
│   ├── plot_eos.py             ← Plot P* vs ρ* for EOS sweeps (solvent_phase / polymer_phase)
│   └── write_tracking.py       ← Log performance data to tracking.txt
│
├── lammps_data/                ← Reserved for small committed .data files (currently empty;
│                                  input_data/ lives outside the repo — see §4)
│
├── README.md                   ← This file
├── CLAUDE_CONTEXT.md           ← Context file for Cowork AI sessions (gitignored)
├── .gitattributes              ← Runs nbstripout on every notebook commit (see §8)
├── lj_units_cheat_sheet.md     ← Unit conversions and parameter reference
├── expanse_lammps_guide.md     ← Expanse-specific setup and GitHub guide
├── slurm_commands_and_compiling.md  ← Cluster commands and LAMMPS build notes
├── slab_data_file_info.md      ← Log of every .data file that has been generated
├── documenting_pod_runs.md     ← HPC performance benchmarks
├── time_vs_atoms.png           ← Benchmark figure (scaling vs system size)
├── time_vs_timesteps.png       ← Benchmark figure (wall-clock vs timestep count)
├── tracking.txt                ← Legacy performance log (not actively maintained)
└── tracking_backup.txt         ← Snapshot of tracking.txt
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
| `isolate_gel.ipynb` | Extracts the swollen polymer (+solvent) from a finished slab run | Pre-step for `shear_slab`; modulus analysis |
| `add_plates_to_gel.ipynb` | Attaches rigid shear plates on the x-faces of an isolated gel (atom type 4, harmonic-bonded to surface polymer) | **Required input for `shear_slab.lmp`** |
| `add_more_plates_to_gel.ipynb` | Variant that adds plates on all six faces | Experimental six-face confinement runs |
| `split_gel_slab.ipynb` | Splits a slab into polymer-only and solvent-only files | Isolated component analysis |
| `pure_polymer.ipynb` | Pure polymer box (no solvent) | EOS and baseline runs |
| `pure_solvent_1.ipynb` | Pure solvent box | EOS and solvent calibration |

Typical shear-modulus pipeline: run `slab_with_support` to equilibrate → `isolate_gel.ipynb` to strip the support/piston → `add_plates_to_gel.ipynb` to attach plates → submit `shear_slab.lmp` with the `*_with_plates.data` file.

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
| `slab_with_flow/` | Piston-forced permeation **or** compression of a pre-equilibrated walled gel — toggled by `variable compression_mode` (0 = permeation, 1 = compression). Phase 0.5 piston pre-compression brings the gel to P* = 1.5 before production. | Measure flux, Dc, pore pressure profiles (mode 0); load gel to measure M (mode 1) |
| `shear_slab/` | Plate-driven xz shear of an isolated swollen gel with attached plates (input from `add_plates_to_gel.ipynb`). Phase 1a NPT (50k) + Phase 1b NVT (100k) + Phase 2 shear with `fix halt` at γ = 10% + Phase 3 NVT production. | Measure G (shear modulus) from ⟨σ_p,xz⟩ / γ |
| `slab_elongation/` | Legacy uniaxial elongation prototype (NPT → fix deform z → NVT) | Superseded by `shear_slab/`; folder kept temporarily |
| `solvent_phase/` | Pure solvent pressure sweep across many state points | Build solvent EOS |
| `solvent_pure/` | Single pure solvent run | Baseline pressure/density check |
| `polymer_phase/` | Pure polymer pressure sweep | Build polymer EOS |
| `polymer_pure/` | Single pure polymer run | Baseline |

> **Legacy `.lmp` files inside simulation folders** (`slab_with_flow_old1.lmp`, `slab_with_flow_old2.lmp`, `solvent_pure_old_no_piston_support.lmp`) are kept for reference but are no longer the active scripts — `run_lammps*.sh` always picks `<folder>.lmp`.

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

For `shear_slab`, `NSTEPS` controls only Phase 3 (production). Phase 1a (NPT equilibration, 50k steps), Phase 1b (NVT lock, 100k steps), and Phase 2 (shear up to 570k steps, halted automatically at γ = 10%) are hardcoded inside `shear_slab.lmp`.

For `slab_with_flow`, the operating mode is set **inside the `.lmp` file**, not the batch file:

```lammps
variable compression_mode equal 1    # 0 = permeation, 1 = compression
```

The same `.lmp` script runs both modes — change the variable, push, and resubmit.

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

These scripts run automatically at the end of every job. If you update a script and want to replot without rerunning the simulation, run them manually from inside the timestamped run directory.

#### Finding your run directory and run_id

Run directories live at `~/Documents/lammps_runs/` and are named:
```
<folder>_<DATANAME>_<INTERACTION>_<TIMESTAMP>/
```

The **run_id** used in all output filenames is `<DATANAME>_<INTERACTION>_<TOTSTEPS>`. You can always find it from the output files themselves:
```bash
ls ~/Documents/lammps_runs/<run_dir>/output_files/stress_data/
# e.g.: stress_tensor_polymer_isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000.dat
#                              └─────────────────────── run_id ───────────────────────────────────────────┘
```

`TOTSTEPS = OLDSTEPS + NSTEPS` from the batch file (for a fresh run, `TOTSTEPS = NSTEPS`).

---

#### Re-running the scripts manually

Navigate into the run directory first — all scripts use `.` as the folder argument:

```bash
cd ~/Documents/lammps_runs/<run_dir>
```

---

**`plot_lammps_log.py`** — T, P, volume convergence (all sim types) + shear diagnostics (shear_slab only, auto-detected)

For `slab_with_flow`, `slab_with_support`, etc.:
```bash
python ~/Documents/lammps_work/scripts/plot_lammps_log.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000
```

For `shear_slab` — pass only the base `DATANAME` as the title, and the full run_id via `--run-id` for file lookups:
```bash
python ~/Documents/lammps_work/scripts/plot_lammps_log.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000 \
    --run-id isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000
```

Or equivalently (same result as the auto-run call from `run_lammps.sh`):
```bash
python ~/Documents/lammps_work/scripts/plot_lammps_log.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000
```

Output plots are saved to `./output_plots/convergence_plots/`.

---

**`plot_stress_profiles.py`** — partial stress and volume fraction profiles (slab_with_flow / slab_with_support)
```bash
python ~/Documents/lammps_work/scripts/plot_stress_profiles.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000 \
    0
# third argument is OLDSTEPS (0 for a fresh run)
```

---

**`plot_piston_data.py`** — piston position and velocity (slab_with_flow only)
```bash
python ~/Documents/lammps_work/scripts/plot_piston_data.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000 \
    0
```

---

**`write_tracking.py`** — performance log (all sim types; not actively maintained)
```bash
python ~/Documents/lammps_work/scripts/write_tracking.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000 \
    ""
# third argument is the TYPE suffix (stress/volume/stressvol or empty string)
```

---

### 7b. Jupyter analysis notebooks

Open these on your MacBook in JupyterLab (`jupyter lab`), pointing them at data files in `flow_data_local/<sim_type>/<RUN_ID>/`. Each notebook has a config cell near the top — only `RUN_ID` and `sim_name` change between runs; all paths derive from those.

**`compression_analysis.ipynb`** ← *Current: Milestone 3 (M) + Milestone 5 (Dc)*
Reads partial stress, volume fraction, and pore pressure data from a `slab_with_flow` compression run (`compression_mode = 1`). Extracts: longitudinal modulus M (both from network-stress integration and Voronoi-tessellated φ_p/φ_s decomposition), cooperative diffusivity Dc from φ_p(z,t) relaxation, and pore-pressure profiles. Inputs: `stress_tensor_polymer_*.dat`, `stress_profile_z_*.dat`, trajectory file. Supersedes `longitudinal_modulus_analysis.ipynb` (still present for reference).

**`permeation_analysis.ipynb`** ← *Current: Milestones 5–6*
Reads volume fraction and stress profiles from `slab_with_flow` permeation runs (`compression_mode = 0`). Extracts pore-pressure profiles φ_p(z), φ_s(z), p_p(z), network stress σ'(z), and solvent flux. Plots each trajectory dump as a separate curve to show temporal evolution. Supersedes `flow_poroelasticity_analysis.ipynb` (still present for reference).

**`shear_analysis.ipynb`** ← *Current: Milestone 4 (G)*
Reads the bulk-region polymer stress tensor from a `shear_slab` Phase 3 production run and extracts G = ⟨σ_p,xz⟩ / γ_cm, plus normal stress differences N1 / N2, x-profile stress plots, and a polymer/solvent poroelastic decomposition. Inputs: `stress_tensor_polymer_*.dat`, `stress_profile_x_polymer_*.dat`, `shear_strain_*.dat`. Atoms within 3σ of either plate are excluded from all stress computes.

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

> **`nbstripout` (notebook output stripping):** the `.gitattributes` file at the repo root runs `nbstripout` as a clean filter on every `*.ipynb` commit, so notebook *outputs* (figures, large data printouts) are never staged — only code and markdown. One-time install on each machine: `pip install nbstripout && (cd ~/Documents/lammps_work && nbstripout --install)`.

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

**`Invalid syntax in variable formula` on first run of `shear_slab.lmp`**
This was caused by `xztilt`/`yztilt` not being valid LAMMPS thermo keywords — the correct keywords are `xz` and `yz`. Already fixed in the current script.

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

*Last updated: 2026-05-15. For questions, contact Dylan Pollard (pollard@ucsb.edu).*
