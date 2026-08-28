# Analysis

[← back to README](../README.md)

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

`TOTSTEPS = NSTEPS` from the batch file (there's no more `OLDSTEPS`-based cumulative counting — see [§5c, editing the batch file](running_simulations.md#5c-editing-the-batch-file)).

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

**`plot_piston_data.py`** — piston position and velocity (`triaxial_compression` single-level runs, `triaxial_permeation`, etc. — compression *sweeps* and `shear_slab` use their own consolidated plotters instead, see [§5b, simulation types](running_simulations.md#5b-simulation-types))
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

**`volume_of_mixing.ipynb`** *(archived 2026-08 → `lammps_work/archive/`; superseded by the calibration-sweep analysis — see `archive/README.md`)*
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

> **Archived 2026-08** — moved to `lammps_work/archive/`; the widom_traj dumps that feed it are now off by default (`SKIP_WIDOM=1`). See `archive/README.md` for how to re-enable.

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
