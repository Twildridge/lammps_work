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
| `triaxial_compression_two_pist/` | **Standard compression workflow since 2026-09-22** (two-piston, 2026-09-16) — drained consolidation at **constant bath pressure**: both solvent reservoirs are closed by NPT-pistons (Marioni et al. 2026, Eq. 3: `setforce 0 0 NULL` + `aveforce NULL NULL v_fz` + `nve`, `fz = ∓P·lx·ly/N − C·v_com`, unthermostatted) held at `P_target`, and a solvent-transparent **load piston** (type 7) loads the network through the **same cumulative strain sweep** as the archived one-piston `triaxial_compression` (`STRAIN_TARGETS` → `COMPRESSIONS`, `_c<strain>` tags, auto-sized holds, symmetric drive: load piston down + support up, `drive_split`). Phases: minimize + 50k Langevin NVT (pistons frozen) → **NPT-piston settle** (`NPT_PISTON_STEPS`, replaces the 2M aniso NPH Phase 0) → symmetric seating (load piston down, support up) → ε = 0 reference → sweep. Input = the converter file (`slab_two_pistons.ipynb`). Piston files carry `[load \| feed \| perm]` columns; `piston_pressure_*` logs the bath check and `permeation_*_c<lvl>` the solvent expelled. Solvent expelled raises the feed piston by ~strain·L₀: the converter's `margin_feed` must cover the deepest level (box-face `fix halt` otherwise). | M, G, D_c, κ vs strain at constant bath pressure (`triaxial_compression_{single,sweep}_two_pist.ipynb`) |
| `triaxial_permeation_two_pist/` | **Standard permeation workflow since 2026-09-22** (two-piston, 2026-09-16) — the two-reservoir replica of Marioni et al. Fig. 2B (rotated 90°, support below the gel): feed piston (type 5, top) at `P_target + DP_PISTON`, permeate piston (type 6, bottom) at `P_target`, both NPT-pistons; solvent–piston WCA always on, so the one-piston reposition / WCA-switch / overlap-relax stages and the z-only NPH are gone. Phases: minimize + Langevin NVT → NPT-piston settle → zero-flux reference (`_ref`) → `P_feed` ramp in 5 stages → one continuous constant-dP drive for `NSTEPS` (**no sweep**). Primary flux = permeate-piston displacement (`permeation_data/permeation_*`: `Q_perm = A·dz_perm/dt`), bead count kept as a cross-check. Halts when the feed reservoir thins to 2 σ or a piston nears a box face. | Q_perm, permeability k = Q L/(A dP), reservoir pressures (`triaxial_permeation_single_two_pist.ipynb`) |
| `triaxial_compression/`, `triaxial_permeation/` | **Archived 2026-09-22 → `archive/simulations/`.** The one-piston decks (single velocity/force-controlled piston + support, periodic slab, aniso-NPH Phase 0). Superseded by the two-piston decks above, which keep their strain sweep, hold sizing, observables and file conventions. Existing runs remain analysable (`triaxial_compression_{single,sweep}.ipynb`, `triaxial_permeation.ipynb`); to run one again `git mv` its folder back under `simulations/` — everything keys on the folder name. | — |
| `shear_slab/` | Plate-driven xz shear of a gel with attached plates. Current input (2026-09-22): the 14000002 periodic slab converted by `scripts/slab_shear_plates.py` (isolate → plates normal to z on the slab's faces; the deck drives them ±x, so z = gap and x = shear like the compression decks; plate separation ~120 σ, laterally periodic), i.e. the same gel state as the two-piston compression sweep. The older isolated-cube input (`isolate_gel.ipynb` → `add_plates_to_gel.ipynb`, plates on the x-faces driven ±z) needs the pre-2026-09-22 deck. Phase 1a **skipped** (`phase1a_mode = 0`, default since 2026-09-22): the plated slab arrives at P* = 1.5 (bulk gel pressure 1.50 at step 0) and the converter deletes the solvent the plate atoms displace (`--solvent-delete 0.87`, calibrated so the bulk gel pressure in the hold is 1.50; without it the plates' excluded volume in a box the gel cannot drain from read 1.53). Do not use the old box NPT (`phase1a_mode = 1`, the isolated-cube protocol) on this input: it targets the thermo `press`, which in this box (vacuum behind the position-controlled plates, plate–gel contact virial) is not the gel pressure. Then Phase 1b NVT (`nsteps_equil` = 100k, plates frozen), then a **cumulative strain sweep** (`STRAINS_LIST=(...)` in the `.batch`, `_g<strain>` tags): each level drives the plates until `fix halt` fires at the target γ (plate-based since 2026-09-22: γ = relative plate z-displacement / plate separation, both prescribed by `fix move`; the old bounding-box denominator overstated the sheared thickness by ~5 %), freezes them and holds for `NSTEPS` while the shear stress relaxes. Since 2026-09-22 the deck follows the triaxial decks' layout and conventions: index-style phase knobs for smoke tests, stem-suffixed whole-run trackers for `plot_lammps_log.py`, per-level loggers spanning drive + hold, epoch-aligned holds (chunk/atom lock), and the same `mobile` group / `c_mobile_temp` thermo names. One physics-affecting addition the same day: the mobile atoms get Maxwell-Boltzmann velocities at T = 1 before Phase 1a (`velocity create`, seed `vel_seed` from `run_lammps.sh`) — the `*_with_plates.data` file has no Velocities section, and the Nose-Hoover thermostat started from rest could not heat the minimised gel (T crept at 0.03–0.5 through Phase 1, then overshot: FENE bonds > R0, lost atoms or a crash in local smoke tests). | Measure G (shear modulus) from ⟨σ_p,xz⟩_bulk / γ; compare with G = (σ′_zz − σ′_xx)/2ε from the two-piston compression sweep |
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

### 5b′. The two-piston sequence (2026-09-16) — the standard triaxial workflow since 2026-09-22

```
slab_with_support (aniso NPH, PISTON_TRANSPARENT=1)  →  final_config_…_14000002.data
    → scripts/slab_two_pistons.ipynb  (converter, on your Mac)  →  …_14000002_two_pist.data
    → copy to ~/Documents/lammps_data/input_data/ on Expanse
    → sbatch triaxial_permeation_two_pist.batch   and/or   sbatch triaxial_compression_two_pist.batch
```

Knobs exported by the `*_two_pist.batch` files and forwarded by `run_lammps.sh` as `-var` (same pattern as
`PISTON_TRANSPARENT`; other engines ignore them):

| env var | `-var` | default | meaning |
|---|---|---|---|
| `PRESS_TARGET` (positional) | `press_target` | 1.5 | bath pressure on the wet pistons |
| `DP_PISTON` | `dp_piston` | 0.1 | permeation only: `P_feed = P_target + dP` (a single value, never a list) |
| `PISTON_MASS` | `piston_mass` | 1000 | mass of every piston bead (types 5/6/7), applied with `mass` after `read_data` |
| `C_PIST_FRAC` | `c_pist_frac` | 1.0 | damping / critical damping; the deck prints `C_crit`, ω and the piston period |
| `NPT_PISTON_STEPS` | `npt_piston_steps` | 1 000 000 | Phase-1 settle length (aim for ≥ 3–5 printed periods) |
| `SETTLE_HALT` | `settle_halt` | 0 | 1 = end the settle early once both pistons are at rest |
| `COMPRESSIONS` | `compressions` | 0.10 | compression only: the cumulative strain sweep (inherited from the archived one-piston deck) |

Deck-only knobs (`-var` override, index-style): `K_solv` (bulk modulus for the damping estimate, 10),
`min_iter`, `phase0_steps`, `ref_avg_steps`, `ref_nfreq`, `ramp_steps`, `t_seat`, `nsteps_settle`,
`v_piston_prod`, `volume_freq`, `thermo_freq`, `strain_freq`, `flux_freq` — handy for short smoke tests.
`DRIVE_SPLIT` (env, default 0.5; forwarded as `-var drive_split`) sets the piston's share of each compression level's
gap closure in both compression decks: 0.5 = symmetric drive (piston down + support up), 1.0 = the old top-only drive.
The archived one-piston `archive/simulations/triaxial_compression/triaxial_compression.lmp` exposes the same kind of knobs since 2026-09-18 (`drive_split`,
`v_piston_prod`, `npt_p05_steps`, `t_seat`, `nsteps_settle`, `ref_avg_steps`, `ref_nfreq`, `min_iter`), e.g. a local
smoke test: `-var npt_p05_steps 2000 -var t_seat 10000 -var nsteps_settle 2000 -var ref_avg_steps 2000 -var ref_nfreq 1000 -var v_piston_prod 0.4 -var min_iter 100 -var hold_auto 0 -var nsteps 4000`.
The two-piston decks have no continuation path (`continue_sim.sh` does not apply): every run is a fresh start
from the converter's data file.

### 5c. Editing the batch file

Open the `.batch` file for your target cluster. The section you edit is at the bottom:

```bash
FOLDER="triaxial_compression_two_pist"
DATANAME="final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000002_two_pist"
INTERACTION="1.0_1.0"   # epsSS_epsSP
NSTEPS=9000000          # flat-hold fallback + setup-file tag (the holds themselves are auto-sized per level)
TYPE=""                 # stress, volume, stressvol, or leave empty
PRESS_TARGET=1.5        # bath pressure on the wet pistons
```

- **`DATANAME`**: the filename of your `.data` file in `lammps_data/input_data/`, without the `.data` extension
- **`INTERACTION`**: `epsSS_epsSP` — solvent–solvent and polymer–solvent LJ well depths (see cheat sheet)
- **`NSTEPS`**: number of timesteps to run (1 million steps ≈ 4100 τ ≈ 19 ns for PEG)
- **`TYPE`**: optional suffix that gets appended to the data file name in output (for labelling stress/volume variants)

To extend a finished run, use `continue_sim.sh` ([§5e below](#5e-continuing-a-run-with-continue_simsh)) rather than resubmitting this batch file — there's no more `OLDSTEPS`-based resubmission path (removed 2026-08-06; it always reran full setup from the prior run's final `.data` file anyway, `continue_sim.sh` is the real restart).

For `shear_slab`, `NSTEPS` is the Phase 3 production hold per strain level (and the output-file tag; the hold is padded up to the next stress-averaging epoch). Phase 1a (`phase1a_mode`: 0 = skip, default; 1 = the old box NPT, `nsteps_npt` = 50k), Phase 1b (NVT, `nsteps_equil` = 100k) and the Phase 2 drive cap (`nsteps_shear` = 570k; every drive halts at its strain target long before) are index-style knobs in `shear_slab.lmp` since 2026-09-22, so a local smoke test is `-var nsteps_equil 4000 -var min_iter 100 -var nsteps_shear_base 30000 -var thermo_p1 500 -var stress_freq 200 -var strains 0.005 0.01 -var nsteps 4000` (verified 2026-09-22 with the LAMMPS GUI `lmp -sf omp` on the Mac; ~4 min). Check the gel pressure with the bulk partial-stress sum (notebook Step 7 or the `stress_tensor_*` files), never with the thermo `press`. The deck wraps its `bond_style hybrid` in `suffix off`/`suffix on`: see the OpenMP note in [common_issues.md](common_issues.md).

### 5d. Submitting a job

```bash
# SSH into the cluster, then:
cd ~/Documents/lammps_work/simulations/triaxial_compression_two_pist
sbatch triaxial_compression_two_pist.batch   # Expanse (the two-piston decks have Expanse batch files only)
# Bridges-2 / Pod variants exist for shear_slab, slab_with_support, … (<name>_bridges.batch / <name>_pod.batch)
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
tail -f ~/Documents/lammps_runs/triaxial_compression_two_pist/triaxial_compression_two_pist_*/log.lammps
```

### 5e. Continuing a run with `continue_sim.sh`

`continue_sim.sh` picks up from where a finished run left off — no restart files, no editing batch scripts. It reads the SLURM output file to find the original working directory and auto-detects all run parameters from there. This is a **real restart** (skip setup, keep going) — contrast with editing `NSTEPS` in a `.batch` file and resubmitting, which is a fresh job that reruns all setup from scratch (see [§5c above](#5c-editing-the-batch-file)).

**When to use it:** you want more steps from a completed run. Not for the two-piston decks
(`triaxial_*_two_pist`, the standard triaxial workflow, which are always fresh starts). As of 2026-08-06, supported for `slab_with_support`, `solvent_pure`, `polymer_pure`, `shear_slab`, and the archived one-piston `triaxial_compression` / `triaxial_permeation` (run it from a folder that holds the original SLURM output file next to the archived `.lmp` / `.batch`; the script keys on the folder name). Not supported: `solvent_phase`/`polymer_phase` (their internal P-sweeps complete in one invocation — "continuing" isn't a meaningful operation) or the `volmix_sweep` pipeline (its own SLURM-chained orchestration). `compress_slab` is a separate project — not wired up here.

#### What it does per folder

The script always skips setup and always resumes whatever the folder's production behavior actually is — for driven simulations that means *continuing the drive*, not freezing into a hold:

| Folder | Skipped | Runs |
|--------|---------|------|
| `slab_with_support` | Soft push-off, minimize, NVT ramp, NPT (`aniso`) warm-up | `aniso` NPT production with volume/dimension outputs |
| `solvent_pure`, `polymer_pure` | Box rescale/harmonic pre-relax, gentle Langevin ramp | More NPT production steps |
| `triaxial_compression` (sweep; archived deck) | Phase 0/1.25/1.5 setup, **and** the non-equilibrium piston drive-to-target | Equilibration-only measurement hold, extended — at whichever `_c<level>` was last reached (auto-detected, never re-sweeps) |
| `shear_slab` (sweep) | Phase 1a/1b equilibration, **and** the non-equilibrium plate shear drive | Equilibration-only production hold, extended — at whichever `_g<strain>` was last reached (auto-detected, never re-sweeps) |
| `triaxial_permeation` (not a sweep; archived deck) | Phase 0/0.5/1.5 setup, **and** the piston reposition/WCA-relax/force-ramp | The constant-pressure forcing drive itself, extended — continuation here means **keep forcing solvent through the gel**, never a passive hold |

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
