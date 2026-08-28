# Running Simulations

[← back to README](../README.md)

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

To extend a finished run, use `continue_sim.sh` ([§5e below](#5e-continuing-a-run-with-continue_simsh)) rather than resubmitting this batch file — there's no more `OLDSTEPS`-based resubmission path (removed 2026-08-06; it always reran full setup from the prior run's final `.data` file anyway, `continue_sim.sh` is the real restart).

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

`continue_sim.sh` picks up from where a finished run left off — no restart files, no editing batch scripts. It reads the SLURM output file to find the original working directory and auto-detects all run parameters from there. This is a **real restart** (skip setup, keep going) — contrast with editing `NSTEPS` in a `.batch` file and resubmitting, which is a fresh job that reruns all setup from scratch (see [§5c above](#5c-editing-the-batch-file)).

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

Benchmarks across system sizes are in [`documenting_pod_runs.md`](../documenting_pod_runs.md).

### 5g. Volume-of-mixing sweep (`volmix_sweep.sh`)

> **Archived 2026-08** — superseded by the calibration sweep (partial-molar-volume NPT sweeps; see `lammps_work/archive/README.md`). The script now lives in `lammps_work/archive/volmix_sweep.sh`; this section is kept because the calibration driver inherits its SLURM machinery.

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

**Analysis:** `volume_of_mixing.ipynb` ([§7b, Jupyter analysis notebooks](analysis.md#7b-jupyter-analysis-notebooks)) syncs `box_dimensions_*.dat` files from Expanse and computes ΔV_mix(P*).

---
