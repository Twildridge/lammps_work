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

For `triaxial_*_two_pist`, `slab_with_support`, the archived one-piston triaxial runs, etc.:
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

**`plot_piston_data.py`** — piston position and velocity (`triaxial_compression` single-level runs, `triaxial_permeation`, etc. — compression *sweeps* and `shear_slab` use their own consolidated plotters instead, see [§5b, simulation types](running_simulations.md#5b-simulation-types)). Multi-piston aware since 2026-09-16: the two-piston files (`# step z_feed z_perm` headers) get one line per piston plus a measured-vs-applied pressure panel and a `Q_perm` / solvent-expelled panel; one-piston files plot exactly as before. `plot_stress_profiles.py` overlays every `stress_z_piston_<sheet>` file in the Piston column and `plot_lammps_log.py` adds the bath-check and `Q_perm` panels to the flow diagnostics.
```bash
python ~/Documents/lammps_work/scripts/plot_piston_data.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000 \
    0
```

---

### 7b. Jupyter analysis notebooks

Open these on your MacBook in JupyterLab (`jupyter lab`), pointing them at data files in `flow_data_local/<sim_type>/<RUN_ID>/`. Each notebook has a config cell near the top — only `RUN_ID` and `sim_name` change between runs; all paths derive from those. For the triaxial-compression notebooks leave `NSTEPS = None`: the `<steps>` tag in the file names is each level's auto-sized hold length and is resolved from the files.

**`triaxial_compression_single.ipynb`** (one-piston runs — one strain level; the one-piston deck was archived on 2026-09-22, existing runs stay analysable)
Analysis of **one** applied-strain level of a `triaxial_compression` run: a single-level run, or one level picked out of a sweep (`LEVEL = "0.15"`). Eleven figures in a fixed order: strain diagnostic; solvent volume fraction (mass fraction / Voronoi / λ-calibrated Voronoi, reference vs compressed); total σ_zz, σ_xx, σ_yy evolution normalised by the bath pressure `P_BARO`; solvent + polymer partial stress with the total superimposed; network stress σ′_zz, σ′_xx, σ′_yy (Terzaghi split) — **reference and final equilibrated state only** with 95 % bands (the evolution curves were dropped 2026-09-17: p_pore is not uniform while the gel consolidates, so subtracting one reservoir value is only valid once equilibrated); piston pressure (linear + log); M (network vs piston); the anisotropy σ′_zz/σ′_xx, σ′_zz/σ′_yy vs step with propagated error bars; G = (σ′_zz − σ′_ii)/2ε from xx and from yy (wall-trimmed membrane interior); the D_c consolidation fit; and κ = D_c/M; then (12) the thermodynamic pressure P_th = −⅓ tr(σ^t) evolution and (13) the osmotic pressure Π = −⅓ tr(σ′) in the reference and final states (both positive under compression in the profiles' sign convention). Figure conventions: every printed number ≤ 3 significant figures, legends placed by `tri.smart_legend` so they never cover data. Layout: **Config → sync → load/compute** (one section), **figures** (one call per figure), and all method notes in a **Notes** markdown at the end. Requires the `sigmaxx_*` / `sigmayy_*` profiles the `.lmp` writes (synced automatically); the pair/bond dumps are *not* needed.

**`triaxial_compression_sweep.ipynb`** (one-piston runs — whole sweep)
Same figures for **every** level of a strain sweep, overlaid (colour = level, reference dashed): profiles per level (network stress and osmotic pressure as final states only), then M, G, D_c and κ **vs applied strain**, plus the P_th evolution and final-state Π overlays. `COMP_LEVELS` must match `STRAIN_TARGETS=(...)` in `triaxial_compression.batch`. Per-level numbers come from the same `load_level` as the single-level notebook, so the two never disagree.

**`lib/triaxial.py`** — the analysis code behind both notebooks
Readers, the Terzaghi network/pore split, the drift-tested block-bootstrap plateau window behind M_piston, the G estimate, the D_c consolidation fit + hold-adequacy check, the λ-calibrated Voronoi φ_s (via `lib/volfrac.py`), the Expanse sync and every figure. Knobs live in `tri.Config` (defaults documented there; the ones that matter for M/G/D_c/κ are spelled out in each notebook's Config cell). Edit the module, re-run the first notebook cell (`importlib.reload`) — no kernel restart needed.

**`triaxial_compression.ipynb`** (long-form original, kept)
The full diagnostic notebook the two above were distilled from (2026-09-02). Still the place for the solvent-phase stress W_s,zz/V_solv diagnostics, the ss/pp pair-virial reconstruction and cross-virial check, the piston–gel contact analysis and the Widom-insertion appendix. Sweep conventions: `COMP_LEVELS`/`DETAIL_LEVEL` select which `_c<level>` tags to load; block-bootstrap piston CI; z-grid origin includes box zlo. Supersedes `compression_analysis.ipynb`.

**`triaxial_permeation.ipynb`** (one-piston runs; deck archived 2026-09-22)
Reads a one-piston `triaxial_permeation` run. Eight sections — piston, thickness, total stress, density, permeate, partial-vs-ss stress, the thermodynamic pressure P_th = −⅓ tr(σ^t) evolution and the osmotic pressure Π = −⅓ tr(σ′) (zero-flux reference, and the final state with the caveat that p_pore is not uniform under flow) — with Phase 1.5 reference overlays. Supersedes `permeation_analysis.ipynb`.

**`triaxial_compression_single_two_pist.ipynb`** / **`triaxial_compression_sweep_two_pist.ipynb`** (two-piston compression, 2026-09-16 — the standard decks since 2026-09-22)
Same Config → sync → load → figures pattern and the same eleven (single) / twelve (sweep) figures as the one-piston
notebooks — the load piston is "the piston" (its force is the network load; the piston files carry it in the first
value column, `[load | feed | perm]`) — plus the **wet-piston bath check** (`fig_wet_pistons[_sweep]`: P_feed, P_perm
measured on the sheets vs `P_target` over the hold) and the **solvent expelled** (`fig_solvent_expelled`, feed-piston
rise + permeate-piston descent). `Config(mode="compression", two_pist=True)` points the sync at
`lammps_runs/triaxial_compression_two_pist` and pulls `piston_pressure`, `permeation`, `pressure_reservoirs` too.
The pore-pressure baseline is the feed-reservoir interior (the box top is vacuum in this geometry), rebuilt per
stress snapshot from the measured feed-piston plane (`piston_position`, column `z_feed`) with every bin kept entirely
≥ `Config.res_wall_margin` (3 σ) clear of the piston sheet — the bin touching a wet piston under-reads σ_zz by ~0.07
(depletion layer + wall virial booked on the piston atoms), which before 2026-09-22 shifted every σ′ by +0.016.
The summary prints *plates over the level* (support, load piston, feed and permeate pistons: start → end).

**`triaxial_permeation_single_two_pist.ipynb`** (two-piston permeation, 2026-09-16, standard since 2026-09-22 — no sweep notebook)
`Config(mode="permeation", two_pist=True)`; `tri.load_permeation` builds one dict `P`: total / partial / network
stress and density evolutions (cividis, bold final, zero-flux `_ref` baseline dashed), both wet pistons
(displacement; block-averaged `F_fluid/(lx ly)` from `piston_force_avg` measured vs applied; reservoir virial
pressures dotted, pf cadence since 2026-09-23), the flux `Q_perm` as the **slope of the permeate bead count** over
independent `q_win_steps` windows with its standard error and a slope drift test (`plateau_window_slopes`,
2026-09-23; the `z_perm` slope and the legacy block-bootstrap mean of the piston-velocity trace are drawn beside
it — the velocity blocks are dominated by the sheet's thermal jitter and their scatter is not an error bar), and
the permeability `k = Q_perm L/(A ΔP)` with the applied and the measured ΔP (CIs in quadrature), and — since
2026-09-23 — the three solvent volume-fraction estimators of the compression notebooks (`tri.add_perm_volume_fractions`:
mass fraction from every density snapshot, Voronoi and λ-calibrated Voronoi on `VOR_EVO_FRAMES` frames over the drive
plus `VOR_MAX_FRAMES` inside the steady window). Because the pore pressure falls feed → permeate across the membrane
under flow, the calibration pressure handed to λ(φ_p, P) is chosen per bin by `Config.P_CAL_MODE`: `'pore'` (the
permeation default — the ramp between the measured feed and permeate reservoir baselines across the membrane),
`'const'` (P_CAL everywhere — the compression default, correct for a drained equilibrium) or `'thermo'` (the local
P_th = −⅓ tr σ^t); the same knob applies to `add_volume_fractions` in the compression notebooks. Eleven figures
(`fig_perm_pistons`, `fig_total_stress`, `fig_partial_stress`, `fig_network_stress`, `fig_perm_density`,
`fig_perm_volfrac` — reference vs steady state, the P_local(z) used and the resulting λ(z) —, `fig_perm_volfrac_evolution`,
`fig_perm_flux`, `fig_perm_permeability`, `fig_thermo_pressure`, `fig_osmotic_pressure`).

**Tests** (`scripts/tests/`, 2026-09-16): `run_plot_tests.sh` builds synthetic one-piston and two-piston run
trees (`make_fixtures.py`) and runs the three plotters on both formats; `run_notebook_tests.sh` executes the
two-piston notebooks headlessly on the same tree; `lint_lmp.py` statically checks a deck (definitions before use,
unfix/undump pairing, labels, forwarded `-var`s).

**`bulk_modulus_analysis.ipynb`**
Drained vs. osmotic bulk modulus K. The osmotic K_osm = Π − dW/dV carries the absolute swelling pressure (large); the *drained* K should be computed like M (network stress response), not from the osmotic branch.

**`volume_of_mixing.ipynb`** *(archived 2026-08 → `lammps_work/archive/`; superseded by the calibration-sweep analysis — see `archive/README.md`)*
Computes ΔV_mix(P*) = V_mixed − V_pure_solvent − V_pure_polymer across the pressure sweep (P* = 1.0–2.0). Cell 2 syncs `box_dimensions_*.dat` files directly from Expanse via `paramiko` SFTP — no SSH keys required; prompts for password and TOTP code in the notebook. Subsequent cells parse the box dimension files, time-average volumes over the last 50% of each run, and plot both ΔV_mix and the individual component volumes vs P*. Requires `paramiko` (`pip install paramiko`). Data lands in `flow_data_local/volmix_sweep/p{P}/`.

**`shear_analysis_single.ipynb`** / **`shear_analysis_sweep.ipynb`** (2026-09-23; replace the 13-step `shear_analysis.ipynb`, git history)
The shear counterparts of `triaxial_compression_{single,sweep}.ipynb`, same **Config → sync → load → figures → Notes** layout; all analysis code lives in **`lib/shear.py`**, which imports the statistics, readers, palette and legend placement from `lib/triaxial.py` so the shear numbers are formed with the same code as the compression numbers. `sh.Config` points at `flow_data_local/shear/<RUN_ID>` and Expanse `lammps_runs/shear_slab`. Per level: the plate-based held γ; **G_network** = plateau-averaged interior ⟨σ_p,xz⟩ increment from the γ = 0 reference over γ (t-interval over bins ⊕ reference CI, exactly `M_network`); **G_series** = the fine `stress_series` ⟨σ_p,xz⟩(t) over the drift-tested block-bootstrap plateau window, increment over γ (exactly `M_piston`); **G_total** from σ^t_xz as the poroelastic check; N1/N2; the bulk P_th = −⅓ tr σ^t (the P* = 1.5 check); the transverse D_c fit of u_x(z,t) (even sine modes) and **κ = D_c/G**, the same hydraulic permeability as D_c/M in compression. Ten figures each (single: strain, total stress, partial σ_xz, network σ_p,xz reference + final, stress history, G comparison, N1/N2, D_c fit, κ, P_th; sweep: the overlays plus G vs γ, the stress–strain curve with least-squares and through-origin slopes, N1/N2, D_c, κ, P_th vs level). Needs the `_ref` files and `stress_series` files the deck writes since 2026-09-23 (older runs are reported absolute, flagged). Since the 2026-09-22 re-layout of `shear_slab.lmp` the strain is **plate-based**: γ = (x_plate,top − x_plate,bot − x₀)/plate_sep from the plate COMs (both prescribed by `fix move`), replacing the old bounding-box-gap denominator and outermost-1 %-bead numerator, which overstated the sheared thickness by ~5 % because `add_plates_to_gel` puts each plate ~1–3 σ inside the outermost beads (the bulk region is now cut `plate_excl` inside the plate planes for the same reason). The `shear_strain_*` files keep their four columns (step, gel_lz_initial, **plate_sep**, **γ**); the measured surface-COM strain goes to `shear_strain_surface_*_g<γ>.dat` as a slip diagnostic. The per-level `shear_strain_*_g<γ>.dat` spans the drive **and** the hold (the notebook and `plot_shear_strain_sweep.py` take its last row as the level's γ, i.e. the held strain), the hold is padded to the next stress-averaging epoch, and stem-suffixed whole-run trackers (`box_dimensions`, `gel_dimensions_rg`, `gel_volume_{bb,rg}`, `shear_strain`) feed `plot_lammps_log.py`'s convergence and shear-diagnostics figures; the files and columns the notebook reads are unchanged. To compare with the two-piston compression G, use the same bath pressure (P* = 1.5: the converter's calibrated solvent deletion around the plates puts the gel there, and Step 7's bulk `<p_total>` is the check — the thermo `press` of this box is not the gel pressure) and the G = (σ′_zz − σ′_xx)/2ε values from `triaxial_compression_sweep_two_pist.ipynb`.

---

### 7c. Running Python scripts manually on a cluster

You may want to rerun post-processing after a job without relaunching LAMMPS — for example, after updating an analysis script, or to run `cavity_widom.py` which is not called automatically on Bridges-2. **Note:** `cavity_widom.py`'s excess-chemical-potential workflow (including the `--p-ext`/`--exclusion-buffer`/`--piston-eps` flags below) was built specifically for the now-removed `slab_with_flow`; it hasn't been ported to the triaxial decks (two-piston, or the archived one-piston ones), which have no equivalent postprocess.sh hook for it yet. `plot_stress_profiles.py` and `plot_piston_data.py` further down are unaffected — those work for the current folders.

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
