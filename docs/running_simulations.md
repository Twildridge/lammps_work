# Running Simulations

[← back to README](../README.md)

### 5a. How it works

When you submit a job, the flow is:

```
edit the .batch file → sbatch → SLURM queues the job → nodes allocated →
run_lammps.sh creates a timestamped lammps_runs/ directory →
pulls the latest scripts from GitHub → runs LAMMPS → post-processing plots
```

`run_lammps.sh` does the bookkeeping: the timestamped working directory, the symlink to the data file, the trajectory files routed to scratch, and the post-processing after LAMMPS finishes. The `.batch` file is the only thing to edit.

### 5b. Simulation types

| Folder | What it simulates | Typical use |
|--------|------------------|-------------|
| `slab_with_support/` | Gel equilibration (free-swelling) or axial compression with piston. NPT uses **`aniso`** so x, y, z each relax independently to P* — the gel reaches its true equilibrium swelling instead of being locked to the data file's aspect ratio (see **Barostat choice** note at the end of this guide). A `pre_swell` knob scales the lattice constant `a` so the gel starts near its swollen equilibrium (faster convergence; pdamp raised 1→5). **Isolated gel** (a finite block in its bath, no sheets; `includeWalls = False` in the builder): use `slab_with_support_isolated.batch` (the 12-cell cube, `square_rho10_nowalls`: rhoSolv 1.0, renamed from `rho04_…_03_nowalls` on 2026-10-11; the Expanse outputs of its 2026-10 run carry the old stem) or `slab_with_support_isolated_tall.batch` (2026-10-11: the 9 × 9 × 21 tall gel, `tall_rho04_nowalls`, 435 k atoms, ~15 h for 8 M steps on 4 nodes). `slab_with_support_isolated_tall_rot45.batch` (2026-10-25) runs the same protocol on `tall_rho04_nowalls_rot45`, the tall gel's network with the diamond lattice rotated 45° about y (448 k atoms, same cost): its final config, converted with `add_more_plates_to_gel.py --layout x_plates`, is the shear-anisotropy test input (see `building_data_files.md` and the batch header). `slab_with_support_single.batch` (2026-10-25) does the same for `tall_rho04_nowalls_random`, the rot45 block with its strands re-wired into a random f = 4 network by bond switching (no lattice directions: its shear G must agree with its compression G whatever the orientation if the discrepancy is anisotropy; smoke-tested locally on a switched 4 × 3 × 8 mini gel, 200 k steps, T 1.00, P 1.50, no FENE warnings). Both set `PIN_GEL=1`; the tall batch sets `BARO_MODE=iso` because nothing elastic spans the box of an isolated gel, so under `aniso` the box shape has no restoring force (the July-8 run of this gel with sheets collapsed 62 × 62 × 122 → 48 × 47 × 212 and its frozen support was carried 26 σ off by the barostat's `dilate all`). Smoke-test a few hundred thousand steps first: `box_dimensions` keeps the data file's aspect ratio (iso), `gel_dimensions` (bounding box) sits at ~52 × 52 × 115–116 (lattice 47.5 × 47.5 × 110.9 at a = 5.28, plus ~4–5 σ of outermost-bead excursions per axis; the deck path was smoke-tested locally on a 3 × 3 × 7-cell version 2026-10-11), `c_mobile_temp` 1.00. | Equilibrate gel; measure M (longitudinal modulus) |
| `triaxial_compression_two_pist/` | **Standard compression workflow since 2026-09-22** (two-piston, 2026-09-16) — drained consolidation at **constant bath pressure**: both solvent reservoirs are closed by NPT-pistons (Marioni et al. 2026, Eq. 3: `setforce 0 0 NULL` + `aveforce NULL NULL v_fz` + `nve`, `fz = ∓P·lx·ly/N − C·v_com`, unthermostatted) held at `P_target`, and a solvent-transparent **load piston** (type 7) loads the network through the **same cumulative strain sweep** as the archived one-piston `triaxial_compression` (`STRAIN_TARGETS` → `COMPRESSIONS`, `_c<strain>` tags, auto-sized holds, symmetric drive: load piston down + support up, `drive_split`). Phases: minimize + 50k Langevin NVT (pistons frozen) → **NPT-piston settle** (`NPT_PISTON_STEPS`, replaces the 2M aniso NPH Phase 0) → symmetric seating (load piston down, support up) → ε = 0 reference → sweep. Input = the converter file (`slab_two_pistons_periodic.ipynb`). Piston files carry `[load \| feed \| perm]` columns; `piston_pressure_*` logs the bath check and `permeation_*_c<lvl>` the solvent expelled. Solvent expelled raises the feed piston by ~strain·L₀: the converter's `margin_feed` must cover the deepest level (box-face `fix halt` otherwise). | M, G, D_c, κ vs strain at constant bath pressure (`triaxial_compression_{single,sweep}_two_pist.ipynb`) |
| `triaxial_permeation_two_pist/` | **Standard permeation workflow since 2026-09-22** (two-piston, 2026-09-16) — the two-reservoir replica of Marioni et al. Fig. 2B (rotated 90°, support below the gel): feed piston (type 5, top) at `P_target + DP_PISTON`, permeate piston (type 6, bottom) at `P_target`, both NPT-pistons; solvent–piston WCA always on, so the one-piston reposition / WCA-switch / overlap-relax stages and the z-only NPH are gone. Phases: minimize + Langevin NVT → NPT-piston settle → zero-flux reference (`_ref`) → `P_feed` ramp in `n_ramp` × `ramp_steps` = 100 × 10 k = 1 M steps (since 2026-09-28; was 5 × 1 k, which rang the undamped sheets for ~2 M steps; both are `index` variables, `-var n_ramp N -var ramp_steps S`) with the polymer displacement reference reset at the ramp start and its profile (`disp_z_polymer`, every `disp_nfreq` = 50 k steps, block-averaged) live through it → one continuous constant-dP drive for `NSTEPS` (**no sweep**). Primary flux = permeate-piston displacement (`permeation_data/permeation_*`: `Q_perm = A·dz_perm/dt`), bead count kept as a cross-check. Halts when the feed reservoir thins to 2 σ or a piston nears a box face. | Q_perm, permeability k = Q L/(A dP), reservoir pressures, and D_c / M from the polymer displacement with a q(t) flux check (`triaxial_permeation_single_two_pist.ipynb`) |
| `triaxial_compression/`, `triaxial_permeation/` | **Archived 2026-09-22 → `archive/simulations/`.** The one-piston decks (single velocity/force-controlled piston + support, periodic slab, aniso-NPH Phase 0). Superseded by the two-piston decks above, which keep their strain sweep, hold sizing, observables and file conventions. Existing runs remain analysable (`triaxial_compression_{single,sweep}.ipynb`, `triaxial_permeation.ipynb`); to run one again `git mv` its folder back under `simulations/` — everything keys on the folder name. | — |
| `shear_slab/` (**archived 2026-10-28** → `archive/simulations/shear_slab/`) | Plate-driven xz shear of a laterally periodic slab with plates on its z-faces (input from `scripts/lib/slab_shear_plates_periodic.py`; Phase 1a skipped, Phase 1b NVT, a cumulative `STRAINS_LIST` sweep with `_g<strain>` tags, an optional Phase 4 unload). Superseded by `shear_slab_two_pist`; `shear_analysis_{single,sweep}.ipynb` with `DECK = "shear_slab"` still read its runs. The archive README has the details. | (archived) |
| `shear_slab_two_pist/` | **Shear on the two-piston system — the current shear workflow.** One deck since 2026-10-28: **`shear_slab_two_pist.lmp` = the ISOLATED-gel deck**, run by `shear_slab_two_pist_single.batch` (one level; since 2026-10-28 the random-network gel at γ = 0.10) or `shear_slab_two_pist_sweep.batch` (0.01 → 0.03 → 0.05 → 0.1). The laterally periodic deck of 2026-10-03/05 (`shear_slab_two_pist_periodic.lmp` + its batch, `DECK_FILE=shear_slab_two_pist_periodic`) and the fine-cadence batches are in `archive/simulations/shear_slab_two_pist/`. **Isolated-gel deck (2026-10-13; grip plates 2026-10-14):** input = the converter file of an equilibrated isolated gel with **layout `x_plates`** (`lib/add_more_plates_to_gel.py` on the `slab_with_support_single.batch` final config, `..._x_plates.data`: gel, bath, the two wet NPT-pistons that hold the one connected bath at `P_target`, and two solvent-transparent grip plates, types 8/9, parked 3 σ outside the gel's x-faces; a `six_plates` file works too, the deck deletes the other four plates after `read_data`, IDs kept). The gel is periodic in no direction. After the piston settle the deck **seats the two plates onto the gel's left and right x-faces** (its long sides; Phase 1.25, `compress_slab`'s seating: common travel time `t_seat` = 100k to `contact_gap` = 1.12 σ short of the Rg face (`SEAT_MODE=rg`) or the outermost bead (`bb`), then `seat_settle_steps` = 200k fixed). Each **grip = plate + the polymer within `layer_thick` = 3 σ of its seat plane** (plus everything beyond), **rigid in all three directions** (`setforce 0 0 0` + `velocity set`, plain `nve`) and driven −z (left) / +z (right): **x is the gap, z the shear direction** (the transpose of the periodic deck). Fully rigid because a finite block gripped on two faces has a rigid-rotation mode about y that reproduces the prescribed ±dz with zero strain if the grips could slide in x; the periodic network cannot rotate, this one can. Before the grips exist the free gel is pinned by `fix momentum` on the polymer (linear + angular, `pin_gel_every` = 100). Stresses are measured in a **core**: `plate_excl` = 5 σ inside the gripped faces, `core_margin_y` = 6 σ inside the free y faces and `core_margin_z` = 20 σ inside the free z faces (σ_xz must vanish on a free surface, so the top and bottom carry Saint-Venant end zones ~W/2 deep); profiles along x (`stress_profile_x_*`, `disp_z_polymer{,_cum}_*`, bins over the whole box inside the core y,z window, reduced x/lx) are the gap profiles, and `stress_profile_z_*` (along z inside the core x,y window) is the end-zone diagnostic for choosing `core_margin_z`. Strain: γ = relative z-displacement of the two plates / `W_drv` = the x-distance between the inner edges of the grip layers; tracking slabs (`track_thick` = 2 σ) give the interior strain; `shear_strain_rg_*` holds Rg_xx², Rg_xz. `plate_pressure_*`: the plates' seat planes (live from step 0, so the seating travel is in the file), `S_r`, `S_l` (constraint z-force / Rg face area = σ_xz at each grip), `N_r`, `N_l` (net normal traction the grips carry) and `Np_r`, `Np_l` (polymer-only x-force on each plate / face area = the seated preload). Same phases, sweep, `_g<strain>` tags and knobs as the periodic deck (plus `core_margin_y`, `core_margin_z`, `pin_gel_every`, `seat_mode`, `t_seat`, `seat_settle_steps`). **Periodic deck:** the `shear_slab` protocol without plates on the compression deck's `*_two_pist.data` file read unchanged; the polymer beads within 3 σ of each Rg z-face are rigid in x/y (prescribed ±x translation, free in z), the support and load piston stay frozen and unseated, tracking slabs and `S_top`/`S_bot` as above. (2026-10-05: replaced the first version's COM-steering force, whose thermally free layer beads swamped the D_c profiles.) **Unload (Phase 2c, 2026-10-28, `UNLOAD=1` in both batches):** after every level's hold the grips are driven back at `vshear` until γ ≤ 0, frozen, and the gel relaxes for `UNLOAD_STEPS` (0 = `NSTEPS`), the level's file set again tagged `_u<strain>`; the next level shears from γ = 0. The notebooks read the residual stresses and the transverse relaxation **D_c** there (`lib/shear.py fit_Dc_unload`, diffusive vs one-time model); the hold-onset `disp_z_polymer` and the fine-cadence files are no longer written. No `continue_sim.sh`; plain `bond_style fene` so `-sf omp` is safe. | G, N1/N2, D_c, κ at constant bath pressure from the same notebooks as `shear_slab`: `DECK = "shear_slab_two_pist_iso"` (isolated gel; the lib relabels x ↔ z at read time so the figures stay in the usual gap-"z" frame) or `DECK = "shear_slab_two_pist"` (periodic); the isolated run is the test of whether the laterally periodic sheared slab overstates G |
| `compress_slab/` | Bulk modulus K (rebuilt 2026-09-29 on the template of `triaxial_compression_two_pist`; the 2026-08 fixed-volume deck is in `archive/`). An **isolated** gel (input from `lib/add_more_plates_to_gel.py`) in a bath held at P* by two wet NPT-pistons; six solvent-transparent plates are seated onto its faces (`SEAT_MODE=bb` since 2026-10-19: the outermost beads, zero-load reference) and closed together through `STRAIN_TARGETS` of cumulative **volumetric** strain (one level, 0.10, in `compress_slab_single.batch`), with a measured ε_vol = 0 reference window and auto-sized holds. Files per level `_<hold>_c<level>`; profiles in core columns along z, x, y; since 2026-10-05 also the cumulative polymer displacement since the reference along each axis (`disp_{z,alongx,alongy}_polymer_cum`), whose steady slopes give the network's own ε_vol (`lib/bulk.py` `K_STRAIN='disp'`, the bulk analogue of the compression notebooks' `M_STRAIN`). **Unload (Phase 2c, 2026-10-28, `UNLOAD=1` in both batches):** after a level's hold the six plates open at the drive speed to the seated gaps + 2 `UNLOAD_MARGIN` (3 σ per plate, clear of the re-swollen shell) and freeze, and the gel re-swells freely for `UNLOAD_STEPS` (0 = the hold length), files `_<hold>_u<level>` with the level's full file set; the next level compresses from the reference again. The free cube is the 3-D unconstrained problem: the notebooks fit its re-swelling for `D_m = (K/M) D_c` and `κ = D_m/K` (`lib/bulk.py fit_Dm_unload`); the compression hold gives `K` only (the hold-onset `displace/atom` files and the held-cube `D_c` fit were dropped the same day). Two batches: `compress_slab_single.batch` = one level (`STRAIN_TARGETS=(0.10)`) + its unload, ~8 M steps; `compress_slab_sweep.batch` = `(0.05 0.075 0.10 0.13 0.16)`, each level followed by its unload, ~31 M steps, 48 h. Same env knobs as the two-piston decks plus `SEAT_MODE`, `HOLD_AUTO`, `UNLOAD` / `UNLOAD_STEPS` / `UNLOAD_MARGIN`. Always a fresh start (no `continue_sim.sh`). |
| `solvent_phase/` | Pure solvent pressure sweep across many state points | Build solvent EOS |
| `solvent_pure/` | Single pure solvent run | Baseline pressure/density check |
| `polymer_phase/` | Pure polymer pressure sweep | Build polymer EOS |
| `polymer_pure/` | Single pure polymer run | Baseline |

> **Legacy `.lmp` files inside simulation folders** (`solvent_pure_old_no_piston_support.lmp`) are kept for reference but are no longer the active scripts — `run_lammps*.sh` always picks `<folder>.lmp`.

Each folder contains:
- `<name>.lmp` — the LAMMPS input script (controls the physics; rarely needs editing)
- `<name>.batch` — the SLURM job script (controls cluster resources; **you edit this**)
- (`<name>_bridges.batch` / `<name>_pod.batch` — the Bridges-2 and Pod versions — were archived to `archive/pod_bridges/` on 2026-10-11; Expanse is the only cluster in use)

### 5b′. The two-piston sequence (2026-09-16) — the standard triaxial workflow since 2026-09-22

```
slab_with_support (aniso NPH, PISTON_TRANSPARENT=1)  →  final_config_…_14000002.data
    → scripts/slab_two_pistons_periodic.ipynb  (converter, on your Mac)  →  …_14000002_two_pist.data
    → copy to ~/Documents/lammps_data/input_data/ on Expanse
    → sbatch triaxial_permeation_two_pist_single.batch   and/or   sbatch triaxial_compression_two_pist_sweep.batch
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
`v_piston_prod`, `volume_freq`, `thermo_freq`, `strain_freq`, `flux_freq`, used for short smoke tests.
Local smoke test of `triaxial_compression_two_pist` including the Phase 2c unload (verified 2026-10-09 with the LAMMPS GUI `lmp -sf omp -pk omp 8` on the 290k-atom file, 28 min; the first attempt caught a formula error in the never-executed unload block, so run it after ANY edit to a block the production runs have not reached yet): `-var phase0_steps 200 -var npt_piston_steps 400 -var min_iter 50 -var t_seat 1000 -var nsteps_settle 400 -var ref_avg_steps 400 -var ref_nfreq 200 -var nsteps 2000 -var compressions 0.02 -var hold_auto 0 -var unload 1 -var unload_steps 2000 -var disp_fine_nfreq 500` (plus the usual `dataname` / `interaction` / `epsSS` / `epsSP` / `press_target` / `vel_seed` / `skip_widom` / `compressions_list` and the `piston_mass` / `c_pist_frac` / `settle_halt` / `drive_split` / `disp_fine_bin` / `unload_margin` vars `run_lammps.sh` forwards; work dir with `data_files/<dataname>.data`, `traj_files/` and the `output_files/*` subfolders). Pass: 15 `_<hold>_u0.02` files (all the level loggers plus `disp_z_polymer{,_fine}`), the log's `[unload 0.02] retraction done` line with the gap back at `L0_gap` and the fine profile files with one block per `disp_fine_nfreq` steps of hold. LAMMPS gotcha the error came from: `${...}` is NOT expanded inside an immediate `$(...)`, so write such lines as plain `equal` variables or use `v_` references inside the `$()`.
`DRIVE_SPLIT` (env, default 0.5; forwarded as `-var drive_split`) sets the piston's share of each compression level's
gap closure in both compression decks: 0.5 = symmetric drive (piston down + support up), 1.0 = the old top-only drive.
The archived one-piston `archive/simulations/triaxial_compression/triaxial_compression.lmp` exposes the same kind of knobs since 2026-09-18 (`drive_split`,
`v_piston_prod`, `npt_p05_steps`, `t_seat`, `nsteps_settle`, `ref_avg_steps`, `ref_nfreq`, `min_iter`), e.g. a local
smoke test: `-var npt_p05_steps 2000 -var t_seat 10000 -var nsteps_settle 2000 -var ref_avg_steps 2000 -var ref_nfreq 1000 -var v_piston_prod 0.4 -var min_iter 100 -var hold_auto 0 -var nsteps 4000`.

Local smoke test of `compress_slab` including the Phase 2c unload (verified 2026-10-28 with the LAMMPS GUI `lmp -sf omp -pk omp 8` on the 726k-atom `..._six_plates` file, 25 min): `-var phase0_steps 200 -var npt_piston_steps 400 -var min_iter 50 -var t_seat 1000 -var nsteps_settle 400 -var ref_avg_steps 400 -var ref_nfreq 200 -var nsteps 2000 -var compressions 0.02 -var compressions_list "0.02" -var hold_auto 0 -var unload 1 -var unload_steps 2000 -var unload_margin 1.0 -var v_piston_prod 0.4 -var thermo_freq 500 -var volume_freq 200 -var strain_freq 100 -var aniso_nevery 100 -var aniso_nrepeat 2 -var aniso_nfreq 200` (plus `dataname` / `interaction` / `epsSS` / `epsSP` / `press_target` / `seat_mode` / `drive_split` / `piston_mass` / `c_pist_frac` / `settle_halt`; work dir with `data_files/<dataname>.data`, `traj_files/` and the `output_files/*` subfolders). Pass: the log's `[unload 0.02] retraction done` line with every gap at its target (seated gap + 2 margins), no `hit its cap` warning, 46 `_<hold>_u0.02` files (the level's loggers, the 18 stress profiles, 3 densities, the bath pressures and the 3 cumulative displacement profiles), `strain sweep complete (unload = 1)`. `scripts/tests/make_fixtures.flatten_to_flow_data` turns the work dir into the notebooks' layout, on which `bulk.load_level` must report the unload and run `fig_Dm`.
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

For the archived `shear_slab`, `NSTEPS` is the Phase 3 production hold per strain level (and the output-file tag; the hold is padded up to the next stress-averaging epoch). Phase 1a (`phase1a_mode`: 0 = skip, default; 1 = the old box NPT, `nsteps_npt` = 50k), Phase 1b (NVT, `nsteps_equil` = 100k) and the Phase 2 drive cap (`nsteps_shear` = 570k; every drive halts at its strain target long before) are index-style knobs in `shear_slab.lmp` since 2026-09-22, so a local smoke test is `-var nsteps_equil 4000 -var min_iter 100 -var nsteps_shear_base 30000 -var thermo_p1 500 -var stress_freq 200 -var strains 0.005 0.01 -var nsteps 4000` (verified 2026-09-22 with the LAMMPS GUI `lmp -sf omp` on the Mac; ~4 min). Check the gel pressure with the bulk partial-stress sum (notebook Step 7 or the `stress_tensor_*` files), never with the thermo `press`. The deck wraps its `bond_style hybrid` in `suffix off`/`suffix on`: see the OpenMP note in [common_issues.md](common_issues.md).

For `shear_slab_two_pist`, `NSTEPS` is likewise the hold per strain level and the file tag; the other lengths are index knobs (`phase0_steps` = 50k, `npt_piston_steps` = `NPT_PISTON_STEPS`, `ref_avg_steps` = 150k, drive cap `nsteps_shear` = 570k). Local smoke test of the **isolated-gel deck** (verified 2026-10-14 with the LAMMPS GUI `lmp -sf omp -pk omp 4`, 24 s on a 3 × 3 × 7-cell mini gel converted with `add_more_plates_to_gel.py` layout `x_plates`; the generator's random solvent must be pushed off first — soft pair + `nve/limit`, then LJ — or the converter's overlap check and the deck's minimize both fail): `-var phase0_steps 200 -var npt_piston_steps 400 -var min_iter 50 -var ref_avg_steps 400 -var ref_nfreq 200 -var series_freq 100 -var stress_freq 100 -var nsteps_shear_base 3000 -var nsteps_shear 6000 -var num_dump_frames 1 -var thermo_freq 500 -var core_margin_y 3 -var core_margin_z 6 -var t_seat 400 -var seat_settle_steps 400 -var strains 0.005 0.01 -var nsteps 2000` (the two core margins must leave a core thicker than 2 bins on the mini gel; the deck quits with a message otherwise). Pass: the plates seat at the Rg faces (`plate_pressure_*` columns 2–3 move from the parked to the seat planes), γ held exactly through each hold (`shear_strain_*` column 4 constant), zero FENE warnings, all `stress_profile_x_*` / `stress_profile_z_*` / `disp_z_polymer*` files written, and `lib/shear.py` (`DECK = "shear_slab_two_pist_iso"`) finds the exact grip planes and runs. **With the unload** (verified 2026-10-28 with the LAMMPS GUI `lmp -sf omp -pk omp 8` on the 542k-atom random-gel `_x_plates` file): add `-var unload 1 -var unload_steps 2000` (and `-var strains 0.02 -var strains_list "0.02" -var nsteps_shear_base 2000 -var nsteps_shear 6000 -var seat_mode rg` on the full-size file); pass = the log's `[unload 0.02] retraction done: gamma = ...` line near 0, the level's file set written again tagged `_2000_u0.02`, and `lib/shear.py` reading the level with `L['U']` and a `D_c` from `fit_Dc_unload` every figure. The periodic deck's smoke test (verified 2026-10-03, ~4 min on the 290k-atom file) is the same line without the core margins and with `-var nsteps_shear_base 30000`, run as `-in shear_slab_two_pist_periodic.lmp`. Submit with `sbatch shear_slab_two_pist_single.batch` (isolated: `DATANAME` = the `..._x_plates` stem, copy it to `input_data/` first) or `sbatch shear_slab_two_pist_periodic.batch` (periodic: the data file is the one already in `input_data/` for the compression runs); defaults `STRAINS_LIST=(0.1)`, `NSTEPS=9000000`, 4 nodes × 32 tasks × 4 threads. **Strain sweep (2026-10-19): `sbatch shear_slab_two_pist_sweep.batch`** — the isolated-gel deck through the cumulative levels `STRAINS_LIST=(0.01 0.03 0.05 0.1)` in one job (the question: is the G ≈ 0.19 of the single 0.1 level strain-dependent?), `NSTEPS=3000000` per hold (the 0.01 level's core σ_xz is read to ~6 % from the shear_4 noise; the whole run is ~14.6 M steps, shear_4's 10.8 M took under ~30 h) and `REF_AVG_STEPS=1000000` (env, new in `run_lammps.sh`: `-var ref_avg_steps`, forwarded only when set; the γ = 0 reference is subtracted from every level and its 150 k default was 20 % of the 0.01 signal). The deck honours **`-var wall_timeout`** since the same date (`timer timeout`, as `slab_with_support.lmp`): a sweep that outruns the job ends cleanly — the level in flight is cut short, the loop exits (later levels are not started; drop the cut level from the notebook's `STRAINS`) and the final `write_data` still runs. **Loop-exit fix (2026-10-21):** sweep 55082229 finished all four levels and then died on `Substitution for illegal variable strains` (no final `write_data`): `next` deletes the index variable once the list is exhausted and skips the following command, so the `sweep_cut` print after the loop read an undefined `${strains}`. The deck now saves each level's tag in `last_strain` before `next` and the normal exit jumps over the budget-cut message; nothing after the loop may reference `${strains}`. Every level's files were intact, and `postprocess.sh` (run by hand) now also draws the stress-strain curve for this deck. LAMMPS checks the clock every `wall_check_every` (1000) steps counted within each `run`, so a run shorter than that is never checked (every production run here is ≥ 26 k steps; a smoke test sets `-var wall_check_every 100`). Analyse with `shear_analysis_sweep.ipynb`, `DECK = "shear_slab_two_pist_iso"`, `STRAINS = ["0.01", "0.03", "0.05", "0.1"]`, `NSTEPS = None`. **Fine-cadence run (2026-10-22): `sbatch shear_slab_two_pist_fine.batch`** — the same four levels with three changes aimed at the transverse relaxation at the start of each hold, which sweep 55082229 could not see (its hold-referenced `disp_z_polymer` is a static pattern from the first 300k-step block on, pinned to the hold-onset `displace/atom` snapshot, and the 5000-step series are at plateau in their first block). `DISP_FINE_NFREQ=250` / `DISP_FINE_SPAN=200000` (`-var disp_fine_nfreq` / `disp_fine_span`, deck defaults 0 / 100000): the cumulative profile `v_uz_cum` on the `disp_z_polymer_cum` bins every 250 steps from the drive start through the first 200k steps of the hold → `disp_z_polymer_fine_*_g<strain>`; the hold is run as two segments (the first padded to a `disp_fine_nf` epoch and capped at the hold, the fine fix unfixed between them; all other hold fixes span both). `NSTEPS_SHEAR_BASE=10000` (`-var nsteps_shear_base`): the drive ~28× faster (0.046 σ/τ per grip, a 0.10 increment in ~9k steps) so the ramp does not hide the relaxation; `STRESS_FREQ=100` (`-var stress_freq`): the halt check every 100 steps keeps the overshoot at ~0.001 (it is also the fix print cadence). `run_lammps.sh` forwards the three new knobs only when set. `lib/shear.py` `fit_Dc_fine` reads the fine file against the late plateau of the coarse `_cum` file (same quantity, same bins, so the onset snapshot cancels) and the notebooks' figures 9b / 7b draw it. Smoke-tested locally 2026-10-22 (LAMMPS GUI `lmp -sf omp -pk omp 6`, the full 435k-atom `_x_plates` file, 9.5 min) with the isolated deck's smoke line plus `-var disp_fine_nfreq 100 -var disp_fine_span 700 -var wall_check_every 100`: fine files with one block per 100 steps from the drive start to hold + 700, both hold segments, the coarse files unchanged (10 blocks), no FENE warnings, `lib/shear.py` runs every figure including the two new ones. Not submitted. **Single-level 5 % run (2026-10-25): `sbatch shear_slab_two_pist_fine_g005.batch`** — the fine batch with `STRAINS_LIST=(0.05)` (one ~4.5 k-step drive from the γ = 0 reference to 0.05, the largest one-step transient the protocol allows, then the 3 M-step hold; ~5.1 M steps, `--time 24:00:00`). Analyse with `shear_analysis_single.ipynb`, `STRAINS = ["0.05"]`. The same batch with `DATANAME` set to the rotated gel's `_x_plates` stem is the anisotropy comparison (same level, same cadence, orientation the only difference).

**`DECK_FILE`** (env, optional, 2026-10-13): run `simulations/$FOLDER/$DECK_FILE.lmp` instead of `$FOLDER.lmp` — for a folder that keeps more than one deck (the archived `shear_slab_two_pist_periodic.batch` exported `DECK_FILE=shear_slab_two_pist_periodic`; no current batch uses it). Output folders and file names still key on `FOLDER`.

### 5d. Submitting a job

```bash
# SSH into the cluster, then:
cd ~/Documents/lammps_work/simulations/triaxial_compression_two_pist
sbatch triaxial_compression_two_pist_sweep.batch   # Expanse (the two-piston decks have Expanse batch files only)
# every folder keeps TWO batch files (2026-10-28): <folder>_single.batch (one level / one run) and <folder>_sweep.batch
# (the strain sweep); here _single = one level 0.10 with the fine-cadence displacement files and the UNLOAD (plates back
# to the seated gap + a 25M-step free re-swelling hold tagged _u0.10, the trusted D_c; knobs UNLOAD / UNLOAD_STEPS /
# UNLOAD_MARGIN), _sweep = the cumulative strain sweep.  The _half / _quarter gels and the permeation _quarter /
# _quarter_smalldp batches are in archive/simulations/ (see archive/README.md).
```

`sbatch` submits the job script to SLURM's queue and returns a job ID at once; SLURM starts the job when the requested nodes are free.

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

`continue_sim.sh` continues a finished run without restart files or batch edits. It reads the SLURM output file to find the original working directory and takes every run parameter from there. It is a restart in the proper sense (setup skipped, production continued), unlike editing `NSTEPS` in a `.batch` file and resubmitting, which is a fresh job that redoes the setup ([§5c above](#5c-editing-the-batch-file)).

**When to use it:** to add steps to a completed run. Not for the two-piston decks
(`triaxial_*_two_pist`, the standard triaxial workflow, which are always fresh starts). As of 2026-08-06, supported for `slab_with_support`, `solvent_pure`, `polymer_pure`, `shear_slab`, and the archived one-piston `triaxial_compression` / `triaxial_permeation` (run it from a folder that holds the original SLURM output file next to the archived `.lmp` / `.batch`; the script keys on the folder name). Not supported: `solvent_phase`/`polymer_phase` (their internal P-sweeps complete in one invocation — "continuing" isn't a meaningful operation) or the `volmix_sweep` pipeline (its own SLURM-chained orchestration). `compress_slab` is a separate project — not wired up here.

#### What it does per folder

The script skips the setup and resumes the folder's production behaviour. For a driven simulation that means continuing the drive, not freezing into a hold:

| Folder | Skipped | Runs |
|--------|---------|------|
| `slab_with_support` | Soft push-off, minimize, NVT ramp, NPT (`aniso`) warm-up | `aniso` NPT production with volume/dimension outputs |
| `solvent_pure`, `polymer_pure` | Box rescale/harmonic pre-relax, gentle Langevin ramp | More NPT production steps |
| `triaxial_compression` (sweep; archived deck) | Phase 0/1.25/1.5 setup, **and** the non-equilibrium piston drive-to-target | Equilibration-only measurement hold, extended — at whichever `_c<level>` was last reached (auto-detected, never re-sweeps) |
| `shear_slab` (sweep; archived deck) | Phase 1a/1b equilibration, **and** the non-equilibrium plate shear drive | Equilibration-only production hold, extended — at whichever `_g<strain>` was last reached (auto-detected, never re-sweeps) |
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
   No other arguments are needed: the current folder determines the data-file naming, whether the run is a sweep and which level to pick up.

   With `~/Documents/lammps_work/scripts` on `$PATH` in `~/.bashrc`, the command is `continue_sim.sh 49772594 500000`.

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

#### What the script does

1. Finds `*.o{job_id}.*` in the current directory (the SLURM output file).
2. Reads the line `Working directory: /path/...` that `run_lammps.sh` printed when the job ran — this gives the original output folder.
3. Looks up the current folder's output-file prefix (`final_config`, `final_tricomp`, `final_triperm`, `final_shear`, `puresolv`, or `purepol`) and finds `<prefix>_*.data` inside that folder.
4. Parses `dataname`, `epsSS`, `epsSP` from the filename (last three `_`-delimited tokens).
5. For `triaxial_compression`/`shear_slab` only: scans `output_files/stress_data/` for the highest `_c<level>`/`_g<level>` tag and passes that single value back as the sweep variable.
6. Creates `continuation_{timestamp}/` inside the original folder and symlinks the data file in.
7. Passes `-var cont 1` (plus the sweep variable where there is one) to LAMMPS; the `jump`/`if` logic in each `.lmp` script then skips the setup and, for the driven decks, the drive or ramp as well, so the continuation carries on the production block rather than adding a hold.

#### Runs longer than the wall limit: automatic chaining (2026-09-30)

Expanse's `compute` partition allows 48 h per job, and a 14 M-step `slab_with_support` run of the 1.13 M-atom `_nowalls` file needs 55–60 h on 4 nodes (the 238 k-atom periodic slab ran at 78 Matom-step/s there). Such a job chains itself instead of being split by hand:

1. `run_lammps.sh` passes the deck a wall-clock budget, `-var wall_timeout` = seconds left in the SLURM job minus `WALL_MARGIN` (default 2400 s; both are env knobs set in the `.batch`). `slab_with_support.lmp` arms it with `timer timeout … every 1000` right at start-up, so setup time counts too.
2. When the budget runs out, the production block in flight returns early (`WARNING: Wall time limit reached`), its snapshot `final_config_*_snap<N>.data` is written as usual, the loop exits, and the deck logs one bookkeeping line — `>>> PRODUCTION_STEPS_DONE <done> OF <nsteps> OLDSTEPS <o> TOTSTEPS <t> STATUS INCOMPLETE` — **without** writing the un-suffixed `final_config_*.data`, which stays the completion marker.
3. With `CHAIN=1` exported by the `.batch`, `run_lammps.sh` then runs `continue_sim.sh --chain <job_id>` *before* post-processing. That reads the bookkeeping line, resumes from the newest snapshot with `cont=1` (no push-off/minimize/warm-up), asks for exactly `totsteps − (oldsteps + done)` more steps and keeps the original `totsteps` tag on every output file. The new job (`cont_<job_id>`, borrowing the `#SBATCH` lines and `export`s of the `.batch` whose `--job-name` matches the output file) chains again if it, too, runs out of time.
4. The final `final_config_<dataname>_<interaction>_<totsteps>.data` is in the **last** `continuation_<timestamp>/` directory under the original run directory (continuations of continuations are siblings there, never nested). Each job's own `continuation_*/log.lammps` and `output_files/` cover its share of the steps.

The margin covers the final `write_data` of a ~1 M-atom system, the MPI teardown grace, the chain submission and the plots. If the chain breaks (node failure, `sbatch` refused), re-launch it from the simulation folder with `continue_sim.sh --chain <last job id>`; `--dry-run` prints the batch script it would submit, `--batch <file>` picks another `.batch` to borrow SLURM settings from. `--chain` on a job whose log says `STATUS COMPLETE` submits nothing. Decks without `wall_timeout` (every other folder) run exactly as before: they ignore the variable, and without the bookkeeping line nothing is chained.

Note the classic `continue_sim.sh <job_id> <nsteps>` *adds* `<nsteps>` on top of a finished run (outputs tagged `<nsteps>`, `oldsteps = 0`); when it finds no final file it falls back to the newest snapshot and says so. Also fixed the same day: `slab_with_support.lmp` now declares `pair_style`/`bond_style` before `read_data` when `cont = 1` (a `write_data` snapshot carries `Pair Coeffs`/`Bond Coeffs` sections, which `read_data` refuses otherwise — the continuation path had never actually run for this deck), and the continuation job now forwards `BARO_MODE`, `PIN_GEL`, `PISTON_TRANSPARENT` and `PRESS_TARGET` to LAMMPS and post-processes through `postprocess.sh` like `run_lammps.sh`.

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

Benchmarks across system sizes are in [`archive/pod_bridges/documenting_pod_runs.md`](../archive/pod_bridges/documenting_pod_runs.md).

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
