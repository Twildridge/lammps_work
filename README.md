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
│   ├── triaxial_compression/   ← CURRENT axial-compression workflow (periodic slab, force-piston) — split from slab_with_flow
│   ├── triaxial_permeation/    ← CURRENT permeation workflow (periodic slab, force-piston) — split from slab_with_flow
│   ├── shear_slab/             ← Shear modulus measurement (plate-driven xz shear; current G workflow)
│   ├── compress_slab/          ← Bulk modulus K measurement (plate-driven isotropic six-face compression; undergrad project, in development)
│   ├── solvent_phase/          ← Pure solvent equation of state sweep
│   ├── solvent_pure/           ← Single-state pure solvent run
│   ├── polymer_phase/          ← Pure polymer equation of state sweep
│   └── polymer_pure/           ← Single-state pure polymer run
│
├── scripts/                    ← Shell and Python scripts
│   ├── run_lammps.sh           ← Job runner for Expanse
│   ├── run_lammps_bridges.sh   ← Job runner for Bridges-2
│   ├── run_lammps_pod.sh       ← Job runner for Pod
│   ├── build_lammps.sh         ← One-shot LAMMPS build script (Expanse-style cmake)
│   ├── git_sync.sh             ← One-command GitHub sync (MacBook + clusters)
│   ├── isolate_gel.py          ← CLI: strip bath/walls and define interior control volume (used by volmix_sweep)
│   ├── add_walls_to_slab.ipynb         ← Build slab data file (no angles)
│   ├── add_walls_with_angles.ipynb     ← Build slab data file (with angles)
│   ├── slab_with_support.ipynb         ← Build basic slab geometry
│   ├── slab_with_support_periodic.ipynb ← CURRENT: xy-periodic crosslinked slab (bonds wrap x,y only; finite-z; p p p; one support+piston per z-period; no side padding). Input for triaxial_* runs
│   ├── slab_with_support_angled.ipynb  ← Build angled-chain slab geometry
│   ├── isolate_gel.ipynb               ← Extract just the swollen gel from a run
│   ├── split_gel_slab.ipynb            ← Split a gel slab into pieces
│   ├── add_plates_to_gel.ipynb         ← Attach shear plates to isolated gel (input for shear_slab)
│   ├── add_more_plates_to_gel.ipynb    ← Variant: plates on all six faces — input for compress_slab
│   ├── pure_polymer.ipynb              ← Build pure polymer data file
│   ├── pure_solvent_1.ipynb            ← Build pure solvent data file
│   ├── triaxial_compression.ipynb           ← CURRENT: M, network stress, pore pressure vs compression-level sweep (triaxial_compression runs)
│   ├── triaxial_permeation.ipynb            ← CURRENT: piston/thickness/stress/density/permeate + partial-vs-ss, with Phase 1.5 reference overlays (triaxial_permeation runs)
│   ├── bulk_modulus_analysis.ipynb          ← Drained vs osmotic bulk modulus K
│   ├── shear_analysis.ipynb                 ← G, N1/N2, stress profiles (shear_slab output)
│   ├── volume_of_mixing.ipynb               ← ΔV_mix(P*) and φ(P*) across pressure sweep (syncs from Expanse via paramiko; requires isolated_* data files)
│   ├── plot_lammps_log.py      ← Plot T, P, volume convergence from log.lammps (+ shear diagnostics)
│   ├── plot_compression_strain_sweep.py ← Plot stress-strain / M across a triaxial_compression sweep
│   ├── plot_shear_strain_sweep.py ← Plot stress-strain / G across a shear_slab sweep
│   ├── split_gel.py            ← CLI: split isolated gel into polymer-only / solvent-only (used by volmix_sweep)
│   ├── plot_stress_profiles.py ← Plot stress and volume fraction profiles
│   ├── plot_piston_data.py     ← Plot piston position and velocity
│   └── plot_eos.py             ← Plot P* vs ρ* for EOS sweeps (solvent_phase / polymer_phase)
│
├── lammps_data/                ← Reserved for small committed .data files (currently empty;
│                                  input_data/ lives outside the repo — see §4)
│
├── README.md                   ← This file
├── shear_slab_notes.md         ← Append-only engineering decision log (gitignored) — NOT current-state; this README is
├── .gitattributes              ← Runs nbstripout on every notebook commit (see §8)
├── lj_units_cheat_sheet.md     ← Unit conversions and parameter reference
├── expanse_lammps_guide.md     ← Expanse-specific setup and GitHub guide
├── slurm_commands_and_compiling.md  ← Cluster commands and LAMMPS build notes
├── slab_data_file_info.md      ← Log of every .data file that has been generated
└── documenting_pod_runs.md     ← HPC performance benchmarks
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

**Activate the Python environment** (required for post-processing scripts):

The run scripts load Python differently depending on the cluster. No extra steps are needed on Expanse or Bridges-2 — `module load anaconda3/...` is called automatically inside `run_lammps.sh`. On **Pod**, you need to activate a conda environment once before submitting jobs:

```bash
# Pod only — one-time setup (and any time you open a new terminal):
module load miniconda
conda activate lammps_analysis
```

If `lammps_analysis` doesn't exist yet on Pod, create it:
```bash
module load miniconda
conda create -n lammps_analysis python=3.11 numpy scipy matplotlib pandas -y
conda activate lammps_analysis
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
| `add_more_plates_to_gel.ipynb` | Variant that adds plates on all six faces | **Required input for `compress_slab.lmp`** (bulk modulus K) |
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
| `slab_with_support/` | Gel equilibration (free-swelling) or axial compression with piston. NPT uses **`aniso`** so x, y, z each relax independently to P* — the gel reaches its true equilibrium swelling instead of being locked to the data file's aspect ratio (see **Barostat choice** note at the end of this guide). A `pre_swell` knob scales the lattice constant `a` so the gel starts near its swollen equilibrium (faster convergence; pdamp raised 1→5). | Equilibrate gel; measure M (longitudinal modulus) |
| `triaxial_compression/` | **Current** axial-compression workflow. Periodic (`p p p`) laterally-unconfined slab from `slab_with_support_periodic.ipynb`, no side walls. A **force-controlled piston** (`setforce`+`aveforce`+`nve`) loads the gel in z; the box is not barostatted during loading. Supports a **cumulative pressure sweep** (`COMPRESSIONS` array in the `.batch`, `_c<level>` output tags) — each level is a piston pressure and the run measures the resulting strain, so **M = pressure / strain**. | Measure M vs compression level |
| `triaxial_permeation/` | **Current** permeation workflow. Same periodic slab and force-piston machinery, driving solvent through the network at constant piston force. | Measure flux, Dc, pore-pressure profiles |
| `shear_slab/` | Plate-driven xz shear of an isolated swollen gel with attached plates (input from `add_plates_to_gel.ipynb`). Phase 1a NPT (50k) + Phase 1b NVT (100k) + Phase 2 shear with `fix halt` at γ = 10% + Phase 3 NVT production. | Measure G (shear modulus) from ⟨σ_p,xz⟩ / γ |
| `compress_slab/` | **In development (undergrad project).** Isotropic bulk-modulus analogue of `shear_slab`: an isolated gel with plates on all six faces (input from `add_more_plates_to_gel.ipynb`) is compressed simultaneously along x, y, and z by driving all six plates inward. Steps through a 3-point cumulative volumetric-strain ladder (ε_vol = 0.015, 0.030, 0.045), holding + measuring the equilibrated network stress at each stage, then fits ΔP'_net vs ε_vol (slope = K) — a 3-point linear fit rather than the single-point `K_single` estimate in `bulk_modulus_analysis.ipynb`, which its own header notes carries a biasing assumption. Writes `bulk_modulus_plot_data_*.dat` for `bulk_modulus_analysis.ipynb` to read. No dedicated analysis notebook yet. | Measure drained bulk modulus K |
| `solvent_phase/` | Pure solvent pressure sweep across many state points | Build solvent EOS |
| `solvent_pure/` | Single pure solvent run | Baseline pressure/density check |
| `polymer_phase/` | Pure polymer pressure sweep | Build polymer EOS |
| `polymer_pure/` | Single pure polymer run | Baseline |

> **Legacy `.lmp` files inside simulation folders** (`solvent_pure_old_no_piston_support.lmp`) are kept for reference but are no longer the active scripts — `run_lammps*.sh` always picks `<folder>.lmp`.

Each folder contains:
- `<name>.lmp` — the LAMMPS input script (controls the physics; rarely needs editing)
- `<name>.batch` — the SLURM job script (controls cluster resources; **you edit this**)
- `<name>_bridges.batch` — Bridges-2 version (if applicable)
- `<name>_pod.batch` — Pod version (if applicable)

### 5c. Editing the batch file

Open the `.batch` file for your target cluster. The section you edit is at the bottom:

```bash
FOLDER="triaxial_compression"
DATANAME="final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000000"
INTERACTION="1.0_1.0"   # epsSS_epsSP
NSTEPS=4000000
TYPE=""                 # stress, volume, stressvol, or leave empty
```

- **`DATANAME`**: the filename of your `.data` file in `lammps_data/input_data/`, without the `.data` extension
- **`INTERACTION`**: `epsSS_epsSP` — solvent–solvent and polymer–solvent LJ well depths (see cheat sheet)
- **`NSTEPS`**: number of timesteps to run (1 million steps ≈ 4100 τ ≈ 19 ns for PEG)
- **`TYPE`**: optional suffix that gets appended to the data file name in output (for labelling stress/volume variants)

To extend a finished run, use `continue_sim.sh` (§5e) rather than resubmitting this batch file — there's no more `OLDSTEPS`-based resubmission path (removed 2026-08-06; it always reran full setup from the prior run's final `.data` file anyway, `continue_sim.sh` is the real restart).

For `shear_slab`, `NSTEPS` controls only Phase 3 (production). Phase 1a (NPT equilibration, 50k steps), Phase 1b (NVT lock, 100k steps), and Phase 2 (shear up to 570k steps, halted automatically at γ = 10%) are hardcoded inside `shear_slab.lmp`.

### 5d. Submitting a job

```bash
# SSH into the cluster, then:
cd ~/Documents/lammps_work/simulations/triaxial_compression
sbatch triaxial_compression_bridges.batch    # Bridges-2
sbatch triaxial_compression.batch            # Expanse
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
tail -f ~/Documents/lammps_runs/triaxial_compression_*/log.lammps
```

### 5e. Continuing a run with `continue_sim.sh`

`continue_sim.sh` picks up from where a finished run left off — no restart files, no editing batch scripts. It reads the SLURM output file to find the original working directory and auto-detects all run parameters from there. This is a **real restart** (skip setup, keep going) — contrast with editing `NSTEPS` in a `.batch` file and resubmitting, which is a fresh job that reruns all setup from scratch (see §5c).

**When to use it:** you want more steps from a completed run. As of 2026-08-06, supported for `slab_with_support`, `solvent_pure`, `polymer_pure`, `triaxial_compression`, `triaxial_permeation`, and `shear_slab`. Not supported: `solvent_phase`/`polymer_phase` (their internal P-sweeps complete in one invocation — "continuing" isn't a meaningful operation) or the `volmix_sweep` pipeline (its own SLURM-chained orchestration). `compress_slab` is a separate project — not wired up here.

#### What it does per folder

The script always skips setup and always resumes whatever the folder's production behavior actually is — for driven simulations that means *continuing the drive*, not freezing into a hold:

| Folder | Skipped | Runs |
|--------|---------|------|
| `slab_with_support` | Soft push-off, minimize, NVT ramp, NPT (`aniso`) warm-up | `aniso` NPT production with volume/dimension outputs |
| `solvent_pure`, `polymer_pure` | Box rescale/harmonic pre-relax, gentle Langevin ramp | More NPT production steps |
| `triaxial_compression` (sweep) | Phase 0/1.25/1.5 setup, **and** the non-equilibrium piston drive-to-target | Equilibration-only measurement hold, extended — at whichever `_c<level>` was last reached (auto-detected, never re-sweeps) |
| `shear_slab` (sweep) | Phase 1a/1b equilibration, **and** the non-equilibrium plate shear drive | Equilibration-only production hold, extended — at whichever `_g<strain>` was last reached (auto-detected, never re-sweeps) |
| `triaxial_permeation` (not a sweep) | Phase 0/0.5/1.5 setup, **and** the piston reposition/WCA-relax/force-ramp | The constant-pressure forcing drive itself, extended — continuation here means **keep forcing solvent through the gel**, never a passive hold |

For the two sweep folders, `continue_sim.sh` scans the original run's `output_files/stress_data/` for the highest `_c<level>`/`_g<level>` tag present and passes just that one value back — the `.lmp` script's sweep variable becomes a one-element list, so it runs exactly once at that level and exits, instead of re-driving through the whole ladder.

All the same output files are produced (stress profiles, chemical potential, piston/permeate data, trajectories, `log.lammps`) into a fresh `continuation_{timestamp}/` subfolder, so nothing from the original run is overwritten. The only intentional omissions are the setup trajectory and any ε=0/zero-flux reference files — those already exist from the original run.

#### How to run it

1. **Navigate to the simulation folder** where you submitted the original job:
   ```bash
   cd ~/Documents/lammps_work/simulations/slab_with_support
   ```

2. **Run the continuation** using the SLURM job ID from the output file name and your desired extra steps:
   ```bash
   ~/Documents/lammps_work/scripts/continue_sim.sh <job_id> <nsteps>
   ```
   For example, if your output file is `slab_support.o49772594.exp-14-05`:
   ```bash
   ~/Documents/lammps_work/scripts/continue_sim.sh 49772594 500000
   ```
   That's it — no other arguments needed. The folder you're in determines everything else (data-file naming, whether it's a sweep, which value to auto-detect).

   > **Tip:** Add `~/Documents/lammps_work/scripts` to your `$PATH` in `~/.bashrc` so you can just type `continue_sim.sh 49772594 500000` directly.

3. **Results appear** in a `continuation_{timestamp}/` subfolder inside the original run's working directory:
   ```
   ~/Documents/lammps_runs/slab_with_support_{dataname}_{interaction}_{timestamp}/
   ├── output_files/           ← original run outputs
   ├── log.lammps              ← original log
   └── continuation_20250602_143012/
       ├── output_files/       ← continuation outputs (stress, chempot, etc.)
       ├── traj_files/         ← symlink to scratch
       └── log.lammps          ← continuation log
   ```

#### What the script actually does under the hood

1. Finds `*.o{job_id}.*` in the current directory (the SLURM output file).
2. Reads the line `Working directory: /path/...` that `run_lammps.sh` printed when the job ran — this gives the original output folder.
3. Looks up the current folder's output-file prefix (`final_config`, `final_tricomp`, `final_triperm`, `final_shear`, `puresolv`, or `purepol`) and finds `<prefix>_*.data` inside that folder.
4. Parses `dataname`, `epsSS`, `epsSP` from the filename (last three `_`-delimited tokens).
5. For `triaxial_compression`/`shear_slab` only: scans `output_files/stress_data/` for the highest `_c<level>`/`_g<level>` tag and passes that single value back as the sweep variable.
6. Creates `continuation_{timestamp}/` inside the original folder and symlinks the data file in.
7. Passes `-var cont 1` (plus the sweep variable, if applicable) to LAMMPS, which triggers the `jump`/`if` logic in each `.lmp` script that bypasses setup — and, for the driven sims, bypasses the drive/ramp itself too, so continuation always means "keep going," never "add an artificial hold."

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

### 5g. Volume-of-mixing sweep (`volmix_sweep.sh`)

Computes ΔV_mix(P*) = V_gel_mixed − V_polymer_pure − V_solvent_pure across 11 pressures (P* = 1.0–2.0, step 0.1).

**Pipeline per pressure** (SLURM chain):

1. **`slab_with_support`** — 600k-step NPT equilibration (3 nodes, 128–64 tasks/node). Copies `final_config_*.data` to `~/Documents/lammps_data/slab_with_support/`.
2. **`isolate_gel.py`** — strips support/piston and bath solvent; defines a control volume as the inner 88% of the polymer distribution (cv_percentile = 6%) in all three axes; writes `isolated_*.data` with box = CV bounds.
3. **`gel_mixed`** and **`split_gel.py`** — run concurrently after isolate. `split_gel.py` splits the isolated gel into `_polymer_only.data` and `_solvent_only.data`.
4. **`solvent_pure`** and **`polymer_pure`** — 100k-step NPT runs of the split components; run concurrently.

The sweep processes one pressure at a time (WINDOW = 1). The last `polymer_pure` job auto-submits the next pressure via a lightweight launcher job.

**Submit commands:**

```bash
cd ~/Documents/lammps_work/simulations/volmix_sweep

# Full sweep (slab runs not yet done):
bash volmix_sweep.sh

# Skip slab step if final_configs already exist in lammps_data/slab_with_support/:
bash volmix_sweep.sh --skip-slab

# Resume from a specific pressure index (e.g. P*=1.6 = index 6):
bash volmix_sweep.sh --from 6 --skip-slab
```

**Control volume note:** `isolate_gel.py --cv-percentile` overrides the 6% default at runtime if you want to test a different inset without editing the script.

**Analysis:** `volume_of_mixing.ipynb` (§7b) syncs `box_dimensions_*.dat` files from Expanse and computes ΔV_mix(P*).

---

## 6. Downloading output from clusters

After a job finishes, the output `.dat` files are in the working directory on the cluster. Trajectory files (`.lammpstrj`) are stored separately in scratch storage (the `traj_files/` symlink points there).

**Download output data files to MacBook** (into `flow_data_local/`):
```bash
# From MacBook terminal — adjust path to match your run
rsync -avP <username>@login.expanse.sdsc.edu:\
  ~/Documents/lammps_runs/triaxial_compression_<dataname>_<interaction>_<timestamp>/output_files/ \
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

`TOTSTEPS = NSTEPS` from the batch file (there's no more `OLDSTEPS`-based cumulative counting — see §5c).

---

#### Re-running the scripts manually

Navigate into the run directory first — all scripts use `.` as the folder argument:

```bash
cd ~/Documents/lammps_runs/<run_dir>
```

---

**`plot_lammps_log.py`** — T, P, volume convergence (all sim types) + shear diagnostics (shear_slab only, auto-detected)

For `triaxial_compression`, `triaxial_permeation`, `slab_with_support`, etc.:
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

**`plot_stress_profiles.py`** — partial stress and volume fraction profiles (`slab_with_support`, `triaxial_compression` single-level runs, `triaxial_permeation`, etc.)
```bash
python ~/Documents/lammps_work/scripts/plot_stress_profiles.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000 \
    0
# third argument is OLDSTEPS — always 0 now (the old batch-file continuation path was removed 2026-08-06)
```

---

**`plot_piston_data.py`** — piston position and velocity (`triaxial_compression` single-level runs, `triaxial_permeation`, etc. — compression *sweeps* and `shear_slab` use their own consolidated plotters instead, see §5b)
```bash
python ~/Documents/lammps_work/scripts/plot_piston_data.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000 \
    0
```

---

### 7b. Jupyter analysis notebooks

Open these on your MacBook in JupyterLab (`jupyter lab`), pointing them at data files in `flow_data_local/<sim_type>/<RUN_ID>/`. Each notebook has a config cell near the top — only `RUN_ID` and `sim_name` change between runs; all paths derive from those.

**`triaxial_compression.ipynb`** (current)
Reads a `triaxial_compression` run (or cumulative pressure sweep). Extracts the longitudinal modulus M from the piston stress vs. strain relationship (M = pressure / strain), network stress σ′ = total − reservoir baseline, and flat pore-pressure profiles. Sweep conventions: `COMP_LEVELS`/`DETAIL_LEVEL` select which `_c<level>` tags to load; per-level reservoir normalisation; block-bootstrap piston CI; z-grid origin includes box zlo (the stress z-coordinates are plotted as z/Lz with support/piston lines). Supersedes `compression_analysis.ipynb`.

**`triaxial_permeation.ipynb`** (current)
Reads a `triaxial_permeation` run. Six panels — piston, thickness, stress, density, permeate, and partial-vs-ss stress — with Phase 1.5 reference overlays. Supersedes `permeation_analysis.ipynb`.

**`bulk_modulus_analysis.ipynb`**
Drained vs. osmotic bulk modulus K. The osmotic K_osm = Π − dW/dV carries the absolute swelling pressure (large); the *drained* K should be computed like M (network stress response), not from the osmotic branch.

**`volume_of_mixing.ipynb`**
Computes ΔV_mix(P*) = V_mixed − V_pure_solvent − V_pure_polymer across the pressure sweep (P* = 1.0–2.0). Cell 2 syncs `box_dimensions_*.dat` files directly from Expanse via `paramiko` SFTP — no SSH keys required; prompts for password and TOTP code in the notebook. Subsequent cells parse the box dimension files, time-average volumes over the last 50% of each run, and plot both ΔV_mix and the individual component volumes vs P*. Requires `paramiko` (`pip install paramiko`). Data lands in `flow_data_local/volmix_sweep/p{P}/`.

**`shear_analysis.ipynb`**
Reads the bulk-region polymer stress tensor from a `shear_slab` Phase 3 production run and extracts G = ⟨σ_p,xz⟩ / γ_cm, plus normal stress differences N1 / N2, x-profile stress plots, and a polymer/solvent poroelastic decomposition. Inputs: `stress_tensor_polymer_*.dat`, `stress_profile_x_polymer_*.dat`, `shear_strain_*.dat`. Atoms within 3σ of either plate are excluded from all stress computes.

---

### 7c. Running Python scripts manually on a cluster

You may want to rerun post-processing after a job without relaunching LAMMPS — for example, after updating an analysis script, or to run `cavity_widom.py` which is not called automatically on Bridges-2. **Note:** `cavity_widom.py`'s excess-chemical-potential workflow (including the `--p-ext`/`--exclusion-buffer`/`--piston-eps` flags below) was built specifically for the now-removed `slab_with_flow`; it hasn't been ported to `triaxial_compression`/`triaxial_permeation`, which have no equivalent postprocess.sh hook for it yet. `plot_stress_profiles.py` and `plot_piston_data.py` further down are unaffected — those work for the current folders.

All scripts below assume you are **inside the run's working directory** on the cluster:

```bash
cd ~/Documents/lammps_runs/<run_dir>
# e.g. cd ~/Documents/lammps_runs/slab_with_support_slab_support_5beads_tall_rho04_1.0_1.0_20260705_124556
```

Set these variables once at the top of your shell session — everything else is derived from them:

```bash
DATANAME="slab_support_5beads_tall_rho04"
INTERACTION="1.0_1.0"
TOTSTEPS=3000000
EPSSS="${INTERACTION%%_*}"   # first part:  e.g. 1.0
EPSSP="${INTERACTION##*_}"   # second part: e.g. 1.0
SCRIPTS=~/Documents/lammps_work/scripts
```

#### Load Python — Expanse

```bash
module load anaconda3/2021.05/q4munrg
```

#### Load Python — Bridges-2

```bash
module load anaconda3/2024.10-1
```

---

#### `cavity_widom.py` — excess chemical potential μ_ex(z)

This is the main script to run manually, especially on Bridges-2 where it is not called automatically by `run_lammps_bridges.sh`.

**Legacy `slab_with_flow` (compression mode — piston-eps 0, exclusion buffer 2σ); kept as a worked example, folder removed 2026-08-06:**

```bash
python "$SCRIPTS/cavity_widom.py" \
  --traj "traj_files/widom_${DATANAME}_${INTERACTION}_${TOTSTEPS}.lammpstrj" \
  --out-dir "output_files/chemical_potential" \
  --out-stem "${DATANAME}_${INTERACTION}_${TOTSTEPS}" \
  --eps-sp "$EPSSP" --eps-ss "$EPSSS" \
  --n-bins 40 --n-trial 50000 --r-cavity 0.5 \
  --exclusion-buffer 2.0 --piston-eps 0.0 \
  --p-ext 1.8 --temperature 1.0
```

**`slab_with_support` (no piston forcing — piston-eps 1, no exclusion buffer):**

```bash
python "$SCRIPTS/cavity_widom.py" \
  --traj "traj_files/widom_${DATANAME}_${INTERACTION}_${TOTSTEPS}.lammpstrj" \
  --out-dir "output_files/chemical_potential" \
  --out-stem "${DATANAME}_${INTERACTION}_${TOTSTEPS}" \
  --eps-sp "$EPSSP" --eps-ss "$EPSSS" \
  --n-bins 40 --n-trial 50000 --r-cavity 0.5 \
  --piston-eps 1.0 \
  --p-ext 1.5 --temperature 1.0
```

The trajectory file is in `traj_files/` (symlink to scratch). If scratch has been purged, you will need to re-download it from the cluster or rerun the simulation.

---

#### `plot_lammps_log.py` — T, P, volume convergence + μ_ex diagnostics

```bash
python "$SCRIPTS/plot_lammps_log.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}"
```
No `--p-ext` flag needed for any current folder — it was a `slab_with_flow`-only option (now removed).

Output saved to `./output_plots/`.

---

#### `plot_stress_profiles.py` — partial stress and volume fraction profiles

```bash
python "$SCRIPTS/plot_stress_profiles.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0
```

---

#### `plot_piston_data.py` — piston position, velocity, force (`triaxial_compression` single-level runs, `triaxial_permeation`, etc.)

```bash
python "$SCRIPTS/plot_piston_data.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0
```

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
`norm none` requires LAMMPS ≥ March 2020. On older builds, replace the `prof_z_polymer/solvent` fixes with an equivalent `compute reduce` + `fix ave/time` approach (see git history for `slab_with_flow.lmp` pre-2026-08-06 removal for a worked example).

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

## 12. Barostat choice (equilibration)

Match the barostat to the geometry:

- **Free-swelling gel in a solvent bath** (`slab_with_support` equilibration): use `fix npt … aniso P P pdamp`. Each of x, y, z is barostatted independently to the target pressure, so the box adopts whatever aspect ratio balances σxx = σyy = σzz = P and the gel relaxes to its own equilibrium shape — the analogue of a hydrogel free to swell in all directions. **Do not use `iso` here:** `iso` controls only the mean (hydrostatic) pressure and freezes the box aspect ratio, so any anisotropic stress or z-padding baked into the data file is never relaxed. Switching `slab_with_support` from `iso` → `aniso` (2026-06-24) fixed exactly this: the gel had been stuck artificially swollen along z, and with `aniso` it reaches a noticeably taller, true equilibrium swelling. Use `couple xy` only if you must enforce lateral isotropy (a free gel reaches it anyway); use `tri` only to relax shear stress (lets the box tilt).
- **Piston-driven runs** (`triaxial_compression`, `triaxial_permeation`): production is `fix nvt` with the box fixed and the piston as the sole z-actuator — **no box barostat in production** (a barostat would double-control z and fight the piston). The only barostat is the Phase-0.5 pre-equilibration, run with the piston/support temporarily on `nve`+`setforce` so they scale with the box.

  As of **2026-07-29** this Phase-0.5 barostat is **`fix nph z`** (z-only), replacing the earlier `aniso`/`iso` forms. Rationale: transverse (xx, yy) stresses build up in the polymer network during compression or permeation; a scalar (`iso`) or per-axis (`aniso`) barostat would let those transverse stresses perturb the box and drift the reservoir pressure. Barostatting **z only** targets the zz stress component directly, holding the solvent reservoir at Pzz = P* = 1.5 while x, y box dimensions stay fixed at the periodic slab's equilibrium extent. Because the target is now Pzz (not the full scalar `Press`), the old `+0.41` kinetic offset used by `slab_with_flow`'s `iso` convention is dropped — all three scripts target `npt_P05_target = P_target`.

  These three scripts run on the **periodic** slab geometry (`slab_with_support_periodic.ipynb`): laterally periodic (`p p p`), no side walls, one support+piston sheet per z-period. Because the slab fills the x-y plane, x and y are free to relax and the piston stress is simply `c_piston_fz / (lx*ly)` with no multi-period normalisation.

When adding lateral walls to an equilibrated config, `add_walls_to_slab.ipynb` unwraps the gel via image flags before measuring its extent, so a gel that has drifted across a periodic face (a little "image pollution" in the visualizer) does not corrupt the wall placement. Under `boundary p p p` the wrapped sliver itself is harmless to the run (bonds use the minimum image; `compute com` unwraps). With `RECENTER_LATERAL = True` (default) it also shifts the mobile group (gel + solvent) in x/y so the polymer COM lands at the lateral box center — the support/piston plates and the z-axis are left untouched, which keeps gel↔support/piston contact along the loading axis and relies on the plates being laterally larger than the gel (the notebook checks this and warns if the gel would exceed the support footprint).

---

*Last updated: 2026-07-29. For questions, contact Dylan Pollard (pollard@ucsb.edu).*
