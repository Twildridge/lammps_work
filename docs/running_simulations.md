# Running Simulations

[← back to README](../README.md)

### 5a. How it works

```
edit the .batch file → sbatch → batch_prologue.sh pulls the repo and re-runs the fresh batch →
run_lammps.sh creates a timestamped lammps_runs/ directory → mpirun lmp → continuation (CHAIN=1) → postprocess.sh
```

A `.batch` file has six parts: the `#SBATCH` header (with a `#POD` line under it), the line `source "$HOME/Documents/lammps_work/scripts/batch_prologue.sh"`, a short description, the four run parameters (`FOLDER`, `DATANAME`, `INTERACTION`, `NSTEPS`), the exported deck knobs, and the call `bash "$LAMMPS_WORK_DIR/scripts/run_lammps.sh" "$FOLDER" "$DATANAME" "$INTERACTION" "$NSTEPS"`.

- `scripts/batch_prologue.sh` does the git pull and re-exec (so the knobs are read from the current file, not the copy `sbatch` took), sets the Open MPI / UCX environment and loads the modules.
- `scripts/lib/site.sh` holds the site settings: LAMMPS binary, data / runs / scratch paths, `MPI_MAP_BY`, `WALL_MARGIN`, modules.
- `scripts/lib/run_common.sh` holds the knob table and the launch / completion logic shared by `run_lammps.sh` and `continue_sim.sh`.
- `scripts/run_lammps.sh` does the bookkeeping: the timestamped working directory, the deck copied into it, the symlink to the data file, the trajectory files routed to scratch, the wall-clock budget, and the post-processing after LAMMPS finishes.

Inside the working directory every `.data` / `.restart` file is in `data_files/`, every `.dat` series in `output_files/`, and every figure in `output_plots/`. `README.txt` there explains the layout and the file names, and `run_info.txt` records what was run ([repository_layout.md](repository_layout.md#inside-a-run-folder)).

`DRY_RUN=1 bash <file>.batch` (outside SLURM) prints the resolved run and the full `mpirun` line without creating anything.

**Resources.** Every production batch requests 6 nodes x 128 tasks x 1 thread, and `run_lammps.sh` maps ranks with `--map-by ppr:128:node:pe=1` (`MPI_MAP_BY=core`; `node` restores the old round-robin placement). §5f has the benchmark behind this choice.

### 5b. Simulation types

| Folder | What it simulates | Typical use |
|--------|------------------|-------------|
| `slab_with_support/` | Free-swelling equilibration of a gel in its solvent bath | Equilibrated gel for every other deck |
| `chain_length_testing/` | The same deck on periodic slabs of different strand length, ended by an Rg-plateau halt | Equilibrated slabs for the chain-length comparison and the rot45 standard gel |
| `uniaxial_compression/` | Two-piston drained compression at constant bath pressure, strain sweep with optional unload | M, G, D_c, κ vs strain |
| `uniaxial_permeation/` | Two-reservoir permeation at constant dP, then velocity-controlled and frozen-network stages | Q_perm, permeability, reservoir pressures, D_c, M |
| `simple_shear/` | Shear on the two-piston system, laterally periodic gel (default) or isolated gel with grip plates | G, N1/N2 at constant bath pressure |
| `isotropic_compression/` | Six-plate volumetric compression of an isolated gel, with unload | K, D_m, κ |
| `pure_solvent/` | Pure solvent, NPT at one pressure or a pressure sweep | Baseline pressure/density check; solvent EOS |
| `pure_polymer/` | Pure polymer (or the mixed gel), NPT at one pressure or a pressure sweep | Mixing-volume reference; polymer EOS |

Each folder contains:
- `<name>.lmp`, the LAMMPS input script. It controls the physics and rarely needs editing. `run_lammps.sh` always picks `<folder>.lmp` unless `DECK_FILE` names another deck (§5c).
- `<name>_single.batch` (one level or one run) and, where there is a sweep, `<name>_sweep.batch`. These are the SLURM job scripts and the files to edit.

Expanse is the main cluster and Pod the backup; the `#POD` line of each batch carries the Pod allocation (§5d). The older per-cluster `<name>_bridges.batch` / `<name>_pod.batch` files are in `archive/pod_bridges/` and `archive/simulations/`.

The folders were renamed on 2026-11-18 (`triaxial` → `uniaxial`, `shear_slab_two_pist` → `simple_shear`, `compress_slab` → `isotropic_compression`). Data files, run IDs and the Expanse run folders made before then keep the old names.

#### `slab_with_support`

Gel equilibration by free swelling. The barostat default is `aniso` (`BARO_MODE`), so x, y and z each relax independently to P* and the gel reaches its true equilibrium swelling instead of being locked to the data file's aspect ratio ([physics_reference.md](physics_reference.md#barostat-choice-equilibration)). The barostat damping `pdamp` is 5, looser than 1 to avoid ringing during swelling. The data-file generators have a `pre_swell` knob that scales the lattice constant `a` so the gel starts near its swollen equilibrium and converges faster.

Batches:
- `slab_with_support_single.batch` equilibrates an **isolated gel** (a finite block in its bath, no sheets; `includeWalls = False` in the builder) at P* = 1.5. `DATANAME` and `NSTEPS` come from the environment when set: `DATANAME=<generator output stem in input_data/> NSTEPS=8000000 sbatch slab_with_support_single.batch`. Defaults: `slab_support_5beads_tall_rho04_nowalls_random`, 8 M steps, `--time=10:00:00`, `BARO_MODE=iso PIN_GEL=1 PISTON_TRANSPARENT=1 CHAIN=1`.
- `slab_with_support_periodic_random.batch` equilibrates the x/y-periodic random slab (14 M steps, 24 h, `PISTON_TRANSPARENT=1 BARO_MODE=aniso PIN_GEL=0 CHAIN=1`); see [Pipelines](#pipelines).

Settings for an isolated gel:
- `BARO_MODE=iso`, because nothing elastic spans the box of an isolated gel, so under `aniso` the box shape has no restoring force. The July-8 run of the tall gel with sheets collapsed from 62 × 62 × 122 to 48 × 47 × 212 and its frozen support was carried 26 σ off by the barostat's `dilate all`.
- `PIN_GEL=1` zeroes the polymer's linear and angular momentum every `pin_gel_every` = 100 steps, so the free block neither drifts nor turns.

Isolated gels in use:

| Stem | Gel | Atoms |
|---|---|---|
| `square_rho10_nowalls` | 12-cell cube, rhoSolv 1.0 (named `rho04_…_03_nowalls` until 2026-10-11; the Expanse outputs of its 2026-10 run carry the old stem) | |
| `tall_rho04_nowalls` | 9 × 9 × 21 tall gel | 435 k |
| `tall_rho04_nowalls_rot45` | the tall gel with the diamond lattice rotated 45° about y | 448 k |
| `tall_rho04_nowalls_random` | the rot45 block re-wired into a random f = 4 network by bond switching (no lattice directions) | |

An 8 M-step run of the tall gel took about 15 h on 4 nodes at the pre-2026-11-08 layout. The final config, converted with `add_more_plates_to_gel.py` (`--layout x_plates` for `simple_shear`, `six_plates` for `isotropic_compression`), is the input of those decks ([building_data_files.md](building_data_files.md)). The dedicated batches `slab_with_support_isolated{,_tall,_tall_rot45}.batch` are in `archive/simulations/slab_with_support/`.

Smoke-test a few hundred thousand steps before a long run. For the tall gel: `box_dimensions` keeps the data file's aspect ratio (iso), `gel_dimensions` (bounding box) sits at about 52 × 52 × 115–116 (lattice 47.5 × 47.5 × 110.9 at a = 5.28, plus 4–5 σ of outermost-bead excursions per axis), and `c_mobile_temp` is 1.00.

#### `chain_length_testing`

Beads-per-strand sweep of the periodic slab. The deck `chain_length_testing.lmp` is a symlink to `slab_with_support.lmp`, so there is one deck and `run_lammps.sh` keys the run directories on this folder. It equilibrates `slab_support_periodic_{N}beads_tall_rho04_{random,rot45}` (`PISTON_TRANSPARENT=1 BARO_MODE=aniso PIN_GEL=0 CHAIN=1`, 36 h per job) with the Rg-plateau halt. Submission, budgets and the halt are described under [Pipelines](#chain-length-sweep). Analysis: `scripts/chain_length_comparison.ipynb`.

#### `uniaxial_compression`

The standard compression workflow: drained consolidation at constant bath pressure.

- **System.** Both solvent reservoirs are closed by NPT-pistons (Marioni et al. 2026, Eq. 3: `setforce 0 0 NULL` + `aveforce NULL NULL v_fz` + `nve`, `fz = ∓P·lx·ly/N − C·v_com`, unthermostatted) held at `P_target`. A solvent-transparent load piston (type 7) loads the network.
- **Input.** The converter file `…_two_pist.data` (`scripts/slab_two_pistons_periodic.ipynb` or `lib/slab_two_pistons_periodic.py`).
- **Phases.** Minimize and 50k steps of Langevin NVT with the pistons frozen; NPT-piston settle (`NPT_PISTON_STEPS`); symmetric seating (load piston down, support up); ε = 0 reference; the cumulative strain sweep (`COMPRESSIONS`), each level with an auto-sized hold (`HOLD_AUTO=1`, τ₁ = L²/(π² D_c)) or a flat hold of `NSTEPS` (`HOLD_AUTO=0`); optionally the Phase 2c unload after a level (`UNLOAD=1`: plates back to the seated gap, then free re-swelling for `UNLOAD_STEPS`).
- **Drive.** Symmetric: load piston down and support up (`DRIVE_SPLIT`, §5b′).
- **Files.** Level files carry `_<hold>_c<strain>`, unload files `_<hold>_u<strain>`. Piston files carry the columns load, feed, perm. `piston_pressure_*` logs the bath check and `permeation_*_c<lvl>` the solvent expelled.
- **Feed margin.** The solvent expelled raises the feed piston by about strain·L₀, so the converter's `margin_feed` must cover the deepest level; otherwise the box-face `fix halt` ends the run.
- **Lateral momentum pin.** `velocity mobile zero linear` and `fix momentum 100 linear 1 1 0` act on the mobile group from the production thermostat onward (z free). Without it the gel slides along the laterally pinned pistons at the thermal COM speed left by the Langevin phases (about 2.5e-3 σ/τ: 107 σ over the random run's hold, about 190 σ over the 2026-10-10 lattice run's). In runs made before 2026-11-05 the x/y stress-profile files `stress_x_*` / `stress_y_*` are therefore all zeros from mid-hold, because their regions follow the unwrapped polymer COM out of the box. The notebooks never read those files, and every z-binned quantity is region-free and unaffected.
- **No `wall_timeout`.** A run that outgrows the time limit loses the end of its unload (post-processed as PARTIAL).

Batches:
- `uniaxial_compression_single.batch`: one level (ε = 0.10), flat hold of 8 M steps (`NSTEPS_OVERRIDE` changes it; 50 k τ > 3 τ₁ of the hold's slow mode), `DISP_FINE_NFREQ=5000` (displacement profiles every 25 τ on hold and unload), `UNLOAD=1 UNLOAD_STEPS=14000000` (70 k τ = 4.1 τ₁ at D_c 0.12, L 142), `--time=28:00:00`. `DATANAME` comes from the environment; the default is the 5-bead rot45 `_two_pist_pm40` stem.
- `uniaxial_compression_sweep.batch`: `COMPRESSIONS="0.10 0.15 0.20 0.3 0.4 0.5 0.7"`, auto-sized holds, 48 h. For a short trial use `COMPRESSIONS="0.05"`, `HOLD_AUTO=0` and a small `NSTEPS`.
- `uniaxial_compression_random.batch`: the single protocol on the random periodic slab (9 M hold, 25 M-step unload, 48 h).
- `chain_random_periodic.batch`: the bridge job that converts an equilibrated slab and submits the two-piston runs ([Pipelines](#pipelines)).

Analysis: `uniaxial_compression_{single,sweep}.ipynb`.

#### `uniaxial_permeation`

The standard permeation workflow: the two-reservoir replica of Marioni et al. Fig. 2B, rotated 90° with the support below the gel. The feed piston (type 5, top) is at `P_target + DP_PISTON`, the permeate piston (type 6, bottom) at `P_target`; both are NPT-pistons and the solvent–piston WCA is always on. There is no strain sweep: one continuous drive.

Phases, in order:

1. Minimize and Langevin NVT, then the NPT-piston settle.
2. Zero-flux reference (files `_ref`).
3. `P_feed` ramp in `n_ramp` × `ramp_steps` = 100 × 10 k = 1 M steps. Both are `index` variables (`-var n_ramp N -var ramp_steps S`). A 5 × 1 k ramp rang the undamped sheets for about 2 M steps. The polymer displacement reference is reset at the ramp start and its profile (`disp_z_polymer`, every `disp_nfreq` = 50 k steps, block-averaged) is live through the ramp.
4. **Drive** at constant dP for `NSTEPS`. D_c comes from the displacement transient and M from the steady displacement profile.
5. **Phase 2v** (`VCTRL=1`, the deck default): the feed piston is taken off its pressure law and moved at a constant velocity for `VCTRL_STEPS` (batch default 4 M). The velocity is the piston's mean velocity over the last `VCTRL_AVG_STEPS` of the drive (displacement over time; 0 = 40 % of `NSTEPS`), or `VFEED=<speed>` (σ/τ, toward the gel) prescribes it. The permeate piston stays at `PRESS_TARGET`. The flux is prescribed and the pressure drop is measured: κ_v = q L/ΔP_meas on the live network. Files `_<NSTEPS>_v`, plus `vctrl_info_…_v.dat` with the velocity and the window it came from.
6. **Phase 2d, frozen network** (`FROZEN=1`, the deck default): every polymer bead is stopped in its steady state (velocities zeroed, polymer out of the integrator with its forces zeroed, thermostat re-created on the solvent alone) and the same feed velocity is kept for `FROZEN_STEPS` (0 = `NSTEPS`; batch default 4 M): κ_f = q L/ΔP_meas on the rigid network, at a flux identical to stage 5. Files `_<NSTEPS>_f`. The polymer is released afterwards and the feed piston returns to pressure control.
7. **Phase 2c, unload** (`UNLOAD=1`, off by default because D_c comes from the drive transient): every production fix is closed, `P_feed` is ramped back to `P_target` in the same `n_ramp` × `ramp_steps` staircase, and the gel re-swells at zero dP for `UNLOAD_STEPS` (0 = `NSTEPS`; batch default 8 M). Every logger is re-opened on files tagged `_<NSTEPS>_u` (pistons, permeation / flux, gel dimensions, reservoir pressures, σ profiles and solvent density at the stress cadence, `disp_z_polymer` reset at the ramp-down start). The network stress drops to zero at both faces at once, so the re-swelling has the drive's modes (τ₁ = L²/(π² D_c)); `uni.add_perm_unload` fits D_c and M to it (notebook cell 12u). It is the permeation analogue of the compression deck's unload.

The primary flux is the permeate-piston displacement (`permeation_data/permeation_*`: `Q_perm = A·dz_perm/dt`); the bead count is kept as a cross-check. Permeability is k = Q L/(A dP). The deck halts when the feed reservoir thins to 2 σ or a piston nears a box face. The reservoir is measured from the gel's Rg top, so the frozen stage has a third halt (2026-11-24): the feed piston within `feed_frozen_gap` (1.5 σ, `FEED_FROZEN_GAP`) of the outermost frozen bead, which the piston would otherwise run into (job 55987856: the bead is 8 σ above the Rg top, `P_feed` reached 10¹⁴). The log's `[frozen] feed piston …` line gives the room and the step at which the limit will be reached, and `[frozen] ENDED EARLY` says which halt fired. The feed reservoir has to last through stages 5 and 6: at q ≈ 2e-4 σ/τ each 4 M steps takes about 4 σ, and the deck halts a stage cleanly when 2 σ remain.

**Why velocity control.** The stage exists because of where the membrane's thermal noise sits. At fixed dP it is in the flux (σ_q² = 2 k_BT κ/(A L T)); at fixed flux it is in dP. κ is equally precise either way, but q_s(t) is flat under velocity control. The velocity is set once and not steered toward a target pressure (a loop on the pressure would put the noise back into the piston), so dP settles near `DP_PISTON` rather than on it, within the precision of the drive's flux, and κ uses the measured value. The drive is written as two `run`s so the deck can note the feed-piston position at the start of the velocity window; the production files are unaffected apart from one extra frame in the `gel_flow` movie. `VCTRL=0` runs the frozen stage at the applied dP instead (the sequence before 2026-11-16).

**What the frozen stage measures (`FROZEN_AT=end`, the default).** One job gives both permeabilities. Both stages share the flow-compacted geometry, so the frozen κ is what a rigid porous medium of that shape gives, and the gap between it and the live κ is what the network's own motion adds. `uni.add_perm_frozen` reads it, and the notebook reports the ratio κ_f/κ_v (the reconciliation table and figure 4); for weak flow it should be close to 1.

**`FROZEN_AT=ref`** (`-var frozen_at ref`) is a different quantity and a separate job: the network is frozen before the ramp, in the zero-flux reference state, and the whole production runs on the undeformed rigid network. The normal production files are then the frozen flux, and Phases 2v/2d/2c are skipped. Its Darcy κ is the undeformed gel's, the value to hold against the compression side's D_c/M at small strain; the live drive compacts the support side, and its κ is a harmonic mean over that profile. Submit with `FROZEN_AT=ref NSTEPS_OVERRIDE=4000000 DATANAME=<stem> sbatch uniaxial_permeation_single.batch` (about 6 M steps).

Batch: `uniaxial_permeation_single.batch` (48 h). `DATANAME`, `NSTEPS_OVERRIDE` and every knob named above are read from the environment.

The batch plans the drive from the compression run of the same data file. `uniaxial_compression_single.ipynb` (`save_M_record`) writes `simulations/uniaxial_permeation/compression_records/<DATANAME>.json` whenever the level has an unload D_c; after an `lsync` the batch reads it through `scripts/lib/perm_plan.py` and sets `DP_PISTON` = 0.24 M (feed-face strain dP/(2M) = 0.12), `NSTEPS` = 4.1 τ₁ with τ₁ = L0²/(π² D_c) from the unload D_c, and `VCTRL_STEPS` = `FROZEN_STEPS` = 4 M, shortened when the feed piston's travel (flux q = dP D_c/(M L0), plus 10 %) would not fit in the feed reservoir. `python3 scripts/lib/perm_plan.py <DATANAME>` prints the plan without submitting. A value set in the environment wins over the plan. A data file with no record is refused unless `DP_PISTON` and `NSTEPS_OVERRIDE` are both given; the standard 5-bead rot45 gel falls back to dP 0.1 and 10.8 M steps (4.1 τ₁ at the drive D_c 0.147, L 138), and `FROZEN_AT=ref` runs keep dP 0.1. The rules are the constants at the top of `perm_plan.py`. The full sequence has not yet been run on Expanse. Analysis: `uniaxial_permeation_single.ipynb` (D_c / M from the polymer displacement with a q(t) flux check).

#### `simple_shear`

Shear on the two-piston system. The folder holds two decks, and every batch names its deck with `DECK_FILE`. Both systems have the two external NPT pistons that hold the bath at `PRESS_TARGET`.

| Gel | Deck | Batches | Input file | Notebook `DECK` |
|---|---|---|---|---|
| periodic in x and y (the default) | `simple_shear_periodic.lmp` | `simple_shear_periodic_single.batch`, `simple_shear_periodic_sweep.batch` | `<stem>_two_pist_pm40.data`, the file the compression and permeation decks read | `"simple_shear"` |
| isolated block with grip plates | `simple_shear_isolated.lmp` | `simple_shear_isolated_single.batch`, `simple_shear_isolated_sweep.batch` | `<stem>_x_plates.data` from `lib/add_more_plates_to_gel.py` | `"simple_shear_iso"` |

`simple_shear.lmp`, `simple_shear_single.batch` and `simple_shear_sweep.batch` are symlinks to the periodic files. Before 2026-11-14 the unsuffixed names were the isolated deck; `simulations/simple_shear/README.md` lists the old names, and runs made before 2026-11-18 are in `lammps_runs/shear_slab_two_pist/` on the cluster. Both decks use plain `bond_style fene`, so `-sf omp` is safe. Analysis: `simple_shear_{single,sweep}.ipynb`; for the isolated deck the lib relabels x ↔ z at read time so the figures stay in the usual gap-"z" frame. The isolated run is the test of whether the laterally periodic sheared slab overstates G.

Common to both decks: a cumulative strain sweep (`STRAINS`, space-separated) with files tagged `_g<strain>`; `NSTEPS` is the hold per level; the γ = 0 reference window `REF_AVG_STEPS` (deck default 150 k) is subtracted from every level; driven layers of `layer_thick` = 3 σ; tracking slabs of `track_thick` = 2 σ just inside them give the interior strain.

**Periodic deck.** The gel is periodic in x and y and reads the compression deck's two-piston file unchanged, so compression, permeation and shear are measured on one gel state (the larger vacuum below the permeate piston is free). The polymer beads within 3 σ of each Rg z-face are rigid in x/y (prescribed ±x translation, free in z) and are driven until γ reaches the target, then held; G = σ_xz/γ, with `S_top` / `S_bot` logged. The support and load piston stay frozen and unseated. The deck has no unload phase and no `wall_timeout`: a job that outgrows `--time` is killed mid-level, and finished levels stay on disk.

- `simple_shear_periodic_single.batch`: `STRAINS=0.1`, `NSTEPS` 9 M (`NSTEPS_OVERRIDE`), `--time=12:00:00`. Measured: 10.8 M steps took 6.8 h for the 5-bead rot45 gel on 6 nodes (job 55663498) and 7.2 h for the 15-bead gel on 10 nodes (55663503). More levels or a longer hold need `sbatch --time=…` (about 6.5 h per extra 9 M-step level).
- `simple_shear_periodic_sweep.batch`: `STRAINS="0.01 0.03 0.05 0.1"`, `NSTEPS=3000000`, `REF_AVG_STEPS=1000000`, `--time=20:00:00` (about 14.5 M steps, about 9 h at the 5-bead rate). Bigger gels: `sbatch --nodes=10 --time=…`.
- Both read `DATANAME` and `STRAINS` from the environment (default: the 5-bead rot45 pm40 stem): `DATANAME=<stem> sbatch simple_shear_periodic_single.batch`.

**Isolated-gel deck.** The gel is periodic in no direction.

- *Input.* The converter file of an equilibrated isolated gel with layout `x_plates` (`lib/add_more_plates_to_gel.py` on the `slab_with_support_single.batch` final config, `…_x_plates.data`): gel, bath, the two wet NPT-pistons that hold the one connected bath at `P_target`, and two solvent-transparent grip plates (types 8/9) parked 3 σ outside the gel's x-faces. A `six_plates` file works too; the deck deletes the other four plates after `read_data` and keeps the IDs.
- *Free-gel pin.* Before the grips exist the gel is pinned by `fix momentum` on the polymer (linear + angular, `pin_gel_every` = 100).
- *Seating (Phase 1.25).* After the piston settle the two plates are seated onto the gel's left and right x-faces (its long sides): common travel time `t_seat` = 100k to `contact_gap` = 1.12 σ short of the Rg face (`SEAT_MODE=rg`, the deck default) or of the outermost bead (`bb`), then `seat_settle_steps` = 200k fixed.
- *Grips.* Each grip is the plate plus the polymer within `layer_thick` = 3 σ of its seat plane (and everything beyond), rigid in all three directions (`setforce 0 0 0` + `velocity set`, plain `nve`), driven −z (left) / +z (right). x is the gap and z the shear direction, the transpose of the periodic deck. The grips are fully rigid because a finite block gripped on two faces has a rigid-rotation mode about y that reproduces the prescribed ±dz with zero strain if the grips could slide in x; the periodic network cannot rotate, this one can.
- *Core.* Stresses are measured `plate_excl` = 5 σ inside the gripped faces, `core_margin_y` = 6 σ inside the free y faces and `core_margin_z` = 20 σ inside the free z faces. σ_xz must vanish on a free surface, so the top and bottom carry Saint-Venant end zones about W/2 deep.
- *Profiles.* `stress_profile_x_*` and `disp_z_polymer{,_cum}_*` (bins over the whole box inside the core y,z window, reduced x/lx) are the gap profiles. `stress_profile_z_*` (along z inside the core x,y window) is the end-zone diagnostic for choosing `core_margin_z`.
- *Strain.* γ = relative z-displacement of the two plates / `W_drv`, the x-distance between the inner edges of the grip layers. `shear_strain_rg_*` holds Rg_xx², Rg_xz.
- *`plate_pressure_*`.* The plates' seat planes (live from step 0, so the seating travel is in the file); `S_r`, `S_l` (constraint z-force / Rg face area = σ_xz at each grip); `N_r`, `N_l` (net normal traction the grips carry); `Np_r`, `Np_l` (polymer-only x-force on each plate / face area = the seated preload).
- *Unload (Phase 2c, `UNLOAD=1`).* After a level's hold the grips are driven back at `vshear` until γ ≤ 0 and frozen, and the gel relaxes for `UNLOAD_STEPS` (0 = `NSTEPS`). The level's file set is written again tagged `_u<strain>`, and the next level shears from γ = 0. The notebooks read the residual stresses and the transverse relaxation D_c there (`lib/shear.py fit_Dc_unload`, diffusive vs one-time model). The deck skips the unload when the wall budget ran out inside the hold.
- *Wall-clock guard.* The deck honours `-var wall_timeout` (`timer timeout`, as `slab_with_support.lmp`): a sweep that outruns the job ends cleanly. The level in flight is cut short, the loop exits (later levels are not started; drop the cut level from the notebook's `STRAINS`) and the final `write_data` still runs. LAMMPS checks the clock every `wall_check_every` (1000) steps counted within each `run`, so a run shorter than that is never checked (every production run here is ≥ 26 k steps; a smoke test sets `-var wall_check_every 100`).
- *Loop rule.* `next` deletes the index variable once the list is exhausted and skips the following command, so nothing after the sweep loop may reference `${strains}`. The deck saves each level's tag in `last_strain` before `next`. Sweep 55082229 finished all four levels and then died on `Substitution for illegal variable strains` without a final `write_data`; its level files are intact. `postprocess.sh` draws the stress-strain curve for this deck.

Isolated batches:
- `simple_shear_isolated_single.batch`: a plain load-only run. `DATANAME` = the random-gel `…_x_plates` stem (copy it to `input_data/` first), `STRAINS="0.1"`, `NSTEPS=9000000`, `REF_AVG_STEPS=1000000`, `UNLOAD=0`, `--time=16:00:00` for about 11.4 M steps (needs ≥ 0.75 M steps/h; a short budget cuts the hold, whose files stay intact).
- `simple_shear_isolated_sweep.batch`: `STRAINS="0.01 0.03 0.05 0.1"` in one job (the question is whether the G ≈ 0.19 of the single 0.1 level is strain-dependent), `NSTEPS=3000000` per hold (the 0.01 level's core σ_xz is read to about 6 % from the shear_4 noise), `REF_AVG_STEPS=1000000` (the 150 k default was 20 % of the 0.01 signal), `UNLOAD=1 UNLOAD_STEPS=0`, `--time=36:00:00`. The run is about 14.6 M steps plus the unloads. Analyse with `simple_shear_sweep.ipynb`, `DECK = "simple_shear_iso"`, `STRAINS = ["0.01", "0.03", "0.05", "0.1"]` (as written in the batch: "0.1", not "0.10"), `NSTEPS = None`.

**D_c in shear (parked).** Two routes are wired in the isolated deck and switched off.

- *(A) Hold transient.* A fast drive (`NSTEPS_SHEAR_BASE`, steps per 0.10 increment, deck default 285000; `HALT_EVERY`, the drive's halt check alone, so the fix print cadence `STRESS_FREQ` can stay at 1000) plus the fine displacement file (`DISP_FINE_NFREQ`, `DISP_FINE_SPAN`, deck defaults 0 / 100000). The file `disp_z_polymer_fine_*_g<γ>` holds the cumulative u_z(x) on the `disp_z_polymer_cum` bins every N steps from the drive start through the first span of the hold; it has its own chunk compute and the hold is run in two segments. `lib/shear.py fit_Dc_fine` fits the lag to the coarse `_cum` file's late plateau (same quantity, same bins, so the onset snapshot cancels) in sin(2πkẑ) modes, because the drive lag is antisymmetric about the mid-plane, and gives the diffusive reading per mode beside the one-time (G/η_s) reading (single notebook figure 9b, sweep notebook 7b; κ falls back to it).
- *Drive speed limits.* A local sweep on the random gel (γ 0.10) found no FENE or neighbour-list trouble at `NSTEPS_SHEAR_BASE` 1000 / 300 / 100 (grip 0.43 / 1.45 / 4.34 σ/τ). At 100 the free group heats to T 1.23 and the global Nosé–Hoover takes about 5 τ to settle it, so a faster drive wants a streaming-unbiased thermostat temperature (`compute temp/profile` binned along the gap on `nvt_prod`), never Langevin: a fixed-frame drag immobilizes the solvent and manufactures the diffusive mode.
- *(B) Phase 2c unload.* `UNLOAD=1`, `UNLOAD_STEPS`; `fit_Dc_unload`, figures 9/11. Its profile is block-averaged at NSTEPS/10, too coarse unless the unload gets its own logger (not built).
- *Why parked.* The same sweep showed the interior reaching 80 % of the applied strain within 4.5 τ: the relaxation is the viscous G/η_s mode of the entrained solvent, on which D_c k² is a few-percent correction. The project's D_c stays the compression unload's free re-swelling.

The fine-cadence batches that used route A (`simple_shear_fine.batch`, `_fine_g005.batch`, `_fine_g005_rot45.batch`: `DISP_FINE_NFREQ=250`, `DISP_FINE_SPAN=200000`, `NSTEPS_SHEAR_BASE=10000`, `STRESS_FREQ=100`) are in `archive/simulations/simple_shear/`.

#### `isotropic_compression`

Bulk modulus K, on the template of `uniaxial_compression`.

- **System.** An isolated gel (input from `lib/add_more_plates_to_gel.py --layout six_plates`) in a bath held at P* by two wet NPT-pistons. Six solvent-transparent plates are seated onto its faces and closed together through `COMPRESSIONS` of cumulative **volumetric** strain, with a measured ε_vol = 0 reference window and auto-sized holds.
- **Seating.** Both batches set `SEAT_MODE=rg`; the deck default is `bb` (§5b′, seat modes and guards).
- **Files.** Per level `_<hold>_c<level>`. Profiles are in core columns along z, x, y, and include the cumulative polymer displacement since the reference along each axis (`disp_{z,alongx,alongy}_polymer_cum`), whose steady slopes give the network's own ε_vol (`lib/isotropic.py` `K_STRAIN='disp'`, the bulk analogue of the compression notebooks' `M_STRAIN`).
- **Unload (Phase 2c, `UNLOAD=1` in both batches).** After a level's hold the six plates open at the drive speed to the seated gaps + 2 `UNLOAD_MARGIN` (3 σ per plate, clear of the re-swollen shell) and freeze, and the gel re-swells freely for `UNLOAD_STEPS` (0 = the hold length). Files `_<hold>_u<level>` carry the level's full file set. The next level compresses from the reference again.
- **What is fitted.** The compression hold gives K only. The free cube is the 3-D unconstrained problem: the notebooks fit its re-swelling for `D_m = (K/M) D_c` and `κ = D_m/K` (`lib/isotropic.py fit_Dm_unload`). A run submitted with `UNLOAD=0` has no D_m.
- **Knobs.** The two-piston knobs plus `SEAT_MODE`, `HOLD_AUTO`, `UNLOAD` / `UNLOAD_STEPS` / `UNLOAD_MARGIN`. The converter's `margin_feed + margin_perm` must cover the deepest level or the box-face halt ends the run.

Batches, both reading `DATANAME` from the environment (default: the square rho-1.0 cube's `…_six_plates` stem): `DATANAME=<six_plates stem> sbatch isotropic_compression_single.batch`.
- `isotropic_compression_single.batch`: one level (0.10) and its unload, about 8 M steps, `--time=27:00:00`. For a short trial use `HOLD_AUTO=0`, a small `NSTEPS`, `UNLOAD_STEPS=NSTEPS`, `NPT_PISTON_STEPS=100000`.
- `isotropic_compression_sweep.batch`: `COMPRESSIONS="0.05 0.075 0.10 0.13 0.16"`, each level followed by its unload, about 31 M steps, 48 h.
- `chain_cube_six_plates.batch`: the bridge job of the cube pipeline ([Pipelines](#isolated-cube)).

Analysis: `isotropic_compression_{single,sweep}.ipynb`. The 2026-08 fixed-volume deck is in `archive/simulations/isotropic_compression/`.

#### `pure_solvent` and `pure_polymer`

Each folder has one deck with two modes.

- `pure_solvent_single.batch`: NPT at one pressure (`PRESS_TARGET`), block-averaged volume.
- `pure_solvent_sweep.batch` (`PHASE_SWEEP=1`): the deck loops over its own pressure list (P* = 1.0 → 0.005), `NSTEPS_EQ` + `NSTEPS_PROD` steps per state point (200 k + 100 k; `NSTEPS_PROD` must be divisible by 100), and writes `eos_*.dat`. `NSTEPS` is unused.
- `pure_polymer_single.batch`: WCA, aniso NPT at one pressure.
- `pure_polymer_sweep.batch` (`PHASE_SWEEP=1`): the pressure loop, with the attractive LJ potential (cutoff 2.5) and an iso barostat, on a relaxed FENE configuration.

#### Archived decks

`archive/README.md` has the details of everything below. Everything keys on the folder name, so to run an archived deck again `git mv` its folder back under `simulations/`.

- `uniaxial_compression_one_pist/`, `uniaxial_permeation_one_pist/` (archived 2026-09-22): the one-piston decks (single velocity/force-controlled piston + support, periodic slab, aniso-NPH Phase 0). The two-piston decks keep their strain sweep, hold sizing, observables and file conventions. Existing runs remain analysable (`uniaxial_compression_one_pist_{single,sweep}.ipynb`, `uniaxial_permeation_one_pist.ipynb`).
- `simple_shear_plates/` (archived 2026-10-28): plate-driven xz shear of a laterally periodic slab with plates on its z-faces (input from `scripts/lib/slab_shear_plates_periodic.py`). `simple_shear_{single,sweep}.ipynb` with `DECK = "simple_shear_plates"` still read its runs. For those runs, check the gel pressure with the bulk partial-stress sum (notebook Step 7 or the `stress_tensor_*` files), never with the thermo `press`; the deck wraps its `bond_style hybrid` in `suffix off`/`suffix on` (OpenMP note in [common_issues.md](common_issues.md)).
- `archive/simulations/solvent_pure/solvent_pure_old_no_piston_support.lmp`: a legacy deck kept for reference.

### 5b′. The two-piston sequence (2026-09-16) — the standard uniaxial workflow since 2026-09-22

```
slab_with_support (aniso NPH, PISTON_TRANSPARENT=1)  →  final_config_…_14000002.data
    → scripts/slab_two_pistons_periodic.ipynb  (converter, on your Mac)  →  …_14000002_two_pist.data
    → copy to ~/Documents/lammps_data/input_data/ on Expanse
    → sbatch uniaxial_permeation_single.batch   and/or   sbatch uniaxial_compression_sweep.batch
```

The two-piston decks (`uniaxial_*`, `isotropic_compression`, `simple_shear`) have no continuation path: every run is a fresh start from the converter's data file, and `continue_sim.sh` does not apply.

#### Knobs

Exported by the batch files and forwarded by `run_lammps.sh` as `-var` (decks that do not define a knob ignore it):

| env var | `-var` | default | meaning |
|---|---|---|---|
| `PRESS_TARGET` | `press_target` | 1.5 | bath pressure on the wet pistons |
| `DP_PISTON` | `dp_piston` | 0.1 | permeation only: `P_feed = P_target + dP` (a single value, never a list) |
| `PISTON_MASS` | `piston_mass` | 1000 | mass of every piston bead (types 5/6/7), applied with `mass` after `read_data` |
| `C_PIST_FRAC` | `c_pist_frac` | 1.0 | damping / critical damping; the deck prints `C_crit`, ω and the piston period |
| `NPT_PISTON_STEPS` | `npt_piston_steps` | 1 000 000 | Phase-1 settle length (aim for ≥ 3–5 printed periods) |
| `SETTLE_HALT` | `settle_halt` | 0 | 1 = end the settle early once both pistons are at rest |
| `COMPRESSIONS` | `compressions` | 0.10 | compression only: the cumulative strain sweep |
| `DRIVE_SPLIT` | `drive_split` | 0.5 | both compression decks: the piston's share of each level's gap closure. 0.5 = symmetric drive (piston down + support up), 1.0 = top-only drive |
| `L0_AVG_STEPS` | `l0_avg_steps` | 500 000 | window of the reference-thickness means (below) |
| `SAVE_STRAINED` | `save_strained` | 1 | strained-state restart at the end of each hold (below) |

Deck-only knobs (`-var` override, index-style), used for short smoke tests: `K_solv` (bulk modulus for the damping estimate, 10), `min_iter`, `phase0_steps`, `ref_avg_steps`, `ref_nfreq`, `ramp_steps`, `t_seat`, `nsteps_settle`, `v_piston_prod`, `volume_freq`, `thermo_freq`, `strain_freq`, `flux_freq`.

#### Reference thickness: `L0_AVG_STEPS`

`L0_AVG_STEPS` is the window of the reference-thickness means in `uniaxial_compression`, `isotropic_compression` and `uniaxial_permeation`.

- **Why a window.** The free gel's thickness breathes thermally (rms ≈ √(kT·L/(M·A)), correlation time ≈ τ₁: 0.4 σ over about 2 M steps on the 5-bead rot45 gel, 0.5–0.9 σ on the 15-bead one), and the outermost bead adds 0.4–0.5 σ of single-step scatter. The decks therefore do not take single-step readings.
- **What is averaged.** The faces the plates seat on are moving-window means over the last `l0_avg_steps` of the NPT-piston settle (mean offset of the outermost bead from the polymer COM, added to the COM at the seating start). `L0_bb` / `L0_rg` (and `V0_bb` / `V0_rg` in `isotropic_compression`) are means over the end of the seated settle. The permeation deck averages its post-settle `L0` and prints a steady `L_ss` over the end of the drive. Each log line prints the single-step value beside the mean.
- **Limits.** The window removes the fast scatter. It averages the slow breathing only as far as it spans τ₁, so to beat that down as well raise `NPT_PISTON_STEPS` (the free-gel phase the seating faces come from) together with `L0_AVG_STEPS`, at about 40 min of wall time per extra 1 M steps. `L0_AVG_STEPS=0` restores the single-step readings.
- **Analysis side.** The permeation `PERM_L0_STEPS` (bounding-box `L_0` before the ramp, `L_ss` before the unload ramp-down) is 500 k steps.

#### Seat modes and guards (`SEAT_MODE`, `SEAT_GAP_RG`, `FRINGE_MAX`, `PRESTRAIN_MAX`, `RG_STRAIN_MIN`)

These apply to `uniaxial_compression` and `isotropic_compression`; all are forwarded as `-var`. The shear and permeation decks have no plates seated on the outermost bead and do not use them.

The strain bookkeeping is the bounding box in both seat modes: `L0_bb` is the bounding-box thickness of the seated gel, the drive closes the gap by strain × `L0_bb`, and the plates are the bounding box from then on. The convergence plot therefore shows the applied compression (BB) and the network's response (Rg).

- `SEAT_MODE=bb` (deck default) seats the plates on the outermost polymer bead: a zero-load reference, right for a sharp surface (rot45 lattices: `L_bb − L_rg` = 4 and 8 σ).
- `SEAT_MODE=rg` seats on the Rg face (COM ± `L_rg`/2, the half-density face within about 1 σ on both gel types), so the plates fold the fringe onto the dense gel and the reference is the gel in contact with the plates. Its preload and pre-strain are printed (`>>> reference check`) and the moduli are increments from that reference.
- `SEAT_GAP_RG` (deck default 2.5 σ since 2026-11-24, both decks; 1.12 before) is how far outside the Rg face an `rg`-seated plate stops. At 1.12 the plates compress the whole network, because the polymer is depleted next to a plate: the Rg edge ends 1.5–1.9 σ per face inside where the plate lets the beads go. On the wide random 10-bead slab that was a pre-strain of 0.049 in thickness (job 55793779, the run behind M = 0.197); on the random 10-bead cube 0.115 in volume, with the density of the whole interior up by 9–15 % and a plate preload of 0.0128 (job 56000507, stopped by `PRESTRAIN_MAX`). The free Rg size is kept at about 2.7–2.8 σ (cube) and 3.1 σ (slab); 2.5 keeps the plates in contact with a predicted pre-strain of about 2 % (cube, volume) and 1.5 % (slab). `SEAT_GAP_RG=1.12` repeats an older run. `isotropic_compression` writes the value in use into `plate_position_ref` and `lib/isotropic.py` reads it from there (it enters the reference edges `L0` = plate gap − 2 × gap); `bb` seating keeps 1.12.

`bb` fails on a gel with a wide fringe. The wide random 10-bead gel has 13–16 σ of fringe (attached loops, no dangling ends) against a 9 σ closure at strain 0.10: in job 55745206 the plates stopped 3–5 σ outside the dense faces, the Rg thickness fell 2.5 % for a nominal 10 %, and the stress (0.0045 ± 0.0007) never levelled (7,600 SU).

Three checks stop a run that is not what was specified, each with a `>>> GUARD …` line and exit status 1 (job FAILED, files so far kept; a value of 0 turns a check off):

| Knob | Default | When | Quits when |
|---|---|---|---|
| `FRINGE_MAX` | 1.0 | `bb` only, about 2 min in | the fringe exceeds that multiple of the first level's closure (wide random gel 1.7; rot45 5-bead 0.28, 15-bead 0.54) |
| `PRESTRAIN_MAX` | 0.10 | after the seated settle, about 1.5 h in | the seating alone has compressed the Rg thickness by more than that |
| `RG_STRAIN_MIN` | 0.5 | end of every hold, before the unload and the later levels | the Rg strain is below that fraction of the applied strain (rot45: 0.83–0.86; job 55745206: 0.27) |

`isotropic_compression` has the same three guards (per-axis fringe, volumetric Rg pre-strain, Rg volume strain over the applied strain). With the defaults they stop the 12-cell cube in `bb` mode (fringe ≈ 5 σ per axis against 2.4 σ of closure at ε_vol 0.10; the `bb` sweep of 2026-10-19 reached 4–14 % of the applied strain in Rg volume), so that cube needs `SEAT_MODE=rg` (2026-10-11 run: 77 %, pre-strain 0.079). Both `isotropic_compression` batch files set `SEAT_MODE=rg`; the environment overrides it.

Unload margin:
- `uniaxial_compression` (deck default `unload_margin` 0): with `SEAT_MODE=rg` and the unload, add `UNLOAD_MARGIN` ≈ half the fringe so the plates clear the re-swelling gel.
- `isotropic_compression` in `rg` mode opens each plate by `UNLOAD_MARGIN` plus half the largest free-gel fringe (`>>> unload opening` line), so the margin counts from the free gel's outermost-bead face in both modes and 3.0 serves any fringe.

Laptop check of the `rg` seating on the rot45 12-cell cube (fringe 5.2 / 5.9 / 7.3 σ, shortened phases): seating pre-strain 0.018 in Rg volume, preload 0.060; plates at 0.100 gave bounding box 0.099 and Rg volume 0.057 after a 20k-step hold and rising (the `bb` run of 2026-11-13 on the same file: 0.069 and 0.002–0.004 after 3 M steps). The fringe-aware unload opening parsed and set its targets in a second short run, which was stopped during the retraction; it has not yet been run to the end. Status of the wide-gel measurements and the compression rerun: [random10_gel_plan.md](random10_gel_plan.md).

#### Strained-state restart

`uniaxial_compression`, `isotropic_compression`, both `simple_shear` decks and `uniaxial_permeation` write `data_files/strained_<deck>_<stem>_<stage>.restart` at the end of every level's hold, before the unload (permeation: at the end of the constant-dP drive), so the unload can be repeated from the same equilibrated strained state. The file is about 150 MB per level for the 770 k-atom gel; `SAVE_STRAINED=0` turns it off. No deck reads it back yet. Details: [repository_layout.md](repository_layout.md#inside-a-run-folder).

#### Local smoke tests

Run the smoke test after any edit to a block the production runs have not reached yet: the first `uniaxial_compression` unload test caught a formula error in the never-executed unload block. All tests below were run with the LAMMPS GUI build on the Mac. Each needs a work dir with `data_files/<dataname>.data`, `traj_files/` and the `output_files/*` subfolders, and the usual `dataname` / `interaction` / `epsSS` / `epsSP` / `press_target` vars plus the others `run_lammps.sh` forwards.

LAMMPS gotchas these tests caught:
- `${...}` is not expanded inside an immediate `$(...)`. Write such lines as plain `equal` variables or use `v_` references inside the `$()`.
- An apostrophe inside a double-quoted `print` line stops LAMMPS with "Unmatched single quote".
- A `fix ave/chunk` (`ave one`) removed in mid-window leaves its lock on the chunk/atom compute, and the next fix on that compute stops LAMMPS with "Two fix commands using same compute chunk/atom command in incompatible ways". The fine and drive-stage displacement fixes of the two compression decks therefore make their compute with the fix and delete it with it (2026-11-21). Likewise the Phase 2v displacement logger has its own compute (`chunk_zv_polymer`), deleted after the stage.

**`uniaxial_compression` with the Phase 2c unload** (`lmp -sf omp -pk omp 8`, 290k-atom file, 28 min):

```
-var phase0_steps 200 -var npt_piston_steps 400 -var min_iter 50 -var t_seat 1000 -var nsteps_settle 400 -var ref_avg_steps 400 -var ref_nfreq 200 -var nsteps 2000 -var compressions 0.02 -var hold_auto 0 -var unload 1 -var unload_steps 2000 -var disp_fine_nfreq 500
```

Also pass `vel_seed` / `skip_widom` / `compressions_list` and the `piston_mass` / `c_pist_frac` / `settle_halt` / `drive_split` / `disp_fine_bin` / `unload_margin` vars. Pass: 15 `_<hold>_u0.02` files (all the level loggers plus `disp_z_polymer{,_fine}`), the log's `[unload 0.02] retraction done` line with the gap back at `L0_gap`, and the fine profile files with one block per `disp_fine_nfreq` steps of hold.

**`isotropic_compression` with the Phase 2c unload** (`lmp -sf omp -pk omp 8`, 726k-atom `…_six_plates` file, 25 min):

```
-var phase0_steps 200 -var npt_piston_steps 400 -var min_iter 50 -var t_seat 1000 -var nsteps_settle 400 -var ref_avg_steps 400 -var ref_nfreq 200 -var nsteps 2000 -var compressions 0.02 -var compressions_list "0.02" -var hold_auto 0 -var unload 1 -var unload_steps 2000 -var unload_margin 1.0 -var v_piston_prod 0.4 -var thermo_freq 500 -var volume_freq 200 -var strain_freq 100 -var aniso_nevery 100 -var aniso_nrepeat 2 -var aniso_nfreq 200
```

Also pass `seat_mode` / `drive_split` / `piston_mass` / `c_pist_frac` / `settle_halt`, and `disp_fine_nfreq 300` for the fine unload profiles (three more `_u` files, `disp_*_polymer_cum_fine`). Pass: the log's `[unload 0.02] retraction done` line with every gap at its target (seated gap + 2 margins), no `hit its cap` warning, 46 `_<hold>_u0.02` files (the level's loggers, the 18 stress profiles, 3 densities, the bath pressures and the 3 cumulative displacement profiles), and `strain sweep complete (unload = 1)`. `scripts/tests/make_fixtures.flatten_to_flow_data` turns the work dir into the notebooks' layout, on which `iso.load_level` must report the unload and run `fig_Dm`.

**`uniaxial_permeation`, base line** (290k-atom `_two_pist_pm40` file):

```
-var min_iter 100 -var phase0_steps 1000 -var npt_piston_steps 3000 -var ref_avg_steps 2000 -var ref_nfreq 1000 -var n_ramp 3 -var ramp_steps 1000 -var disp_nfreq 2000 -var volume_freq 500 -var thermo_freq 500 -var flux_freq 500 -var nsteps 6000
```

- *Full sequence* (`lmp -sf omp -pk omp 8`, 156k-atom `_two_pist_quarter_perm` file, about 35 min): add `-var vctrl_steps 4000 -var frozen_steps 4000 -var unload 1 -var unload_steps 4000 -var perm_stress_curves 3`. Pass: the log's `[vctrl] done`, `[frozen] done`, `feed piston back to pressure control`, `[unload] done` and `Simulation complete` lines, exit 0; 17 `_6000_v`, 15 `_6000_f` and 23 `_6000_u` files; `z_feed` in the `_v` and `_f` files linear in time at `vfeed_set` to the files' print precision; the files load through `add_perm_vctrl` / `add_perm_frozen` / `add_perm_unload`. The pressures and κ of such a run mean nothing, because the pistons are still ringing after a 3000-step settle.
- *Frozen stage and unload without velocity control:* add `-var frozen 1 -var frozen_steps 4000 -var unload 1 -var unload_steps 6000 -var perm_stress_curves 3`. Pass: the log's `[frozen] polymer frozen` / `[frozen] done` / `[frozen] polymer released` lines and 15 `_6000_f` files; the `[unload] ramp-down done` and `[unload] done` lines and 23 `_6000_u` files (all the production loggers plus `disp_z_polymer`, σ profiles and `solvent_density_z`); `Simulation complete`, exit 0.
- *`FROZEN_AT=ref`* (`-var frozen_at ref`, 12 min): the log's `FROZEN REFERENCE NETWORK` line, every `disp_z_polymer` row exactly zero through the production, no `_f` / `_u` files, `Simulation complete`, exit 0.

**`simple_shear`, isolated deck** (`lmp -sf omp -pk omp 4`, 24 s on a 3 × 3 × 7-cell mini gel converted with `add_more_plates_to_gel.py` layout `x_plates`). The generator's random solvent must be pushed off first (soft pair + `nve/limit`, then LJ), or the converter's overlap check and the deck's minimize both fail.

```
-var phase0_steps 200 -var npt_piston_steps 400 -var min_iter 50 -var ref_avg_steps 400 -var ref_nfreq 200 -var series_freq 100 -var stress_freq 100 -var nsteps_shear_base 3000 -var nsteps_shear 6000 -var num_dump_frames 1 -var thermo_freq 500 -var core_margin_y 3 -var core_margin_z 6 -var t_seat 400 -var seat_settle_steps 400 -var strains 0.005 0.01 -var nsteps 2000
```

The two core margins must leave a core thicker than 2 bins on the mini gel; the deck quits with a message otherwise. Pass: the plates seat at the Rg faces (`plate_pressure_*` columns 2–3 move from the parked to the seat planes), γ is held exactly through each hold (`shear_strain_*` column 4 constant), zero FENE warnings, all `stress_profile_x_*` / `stress_profile_z_*` / `disp_z_polymer*` files written, and `lib/shear.py` (`DECK = "simple_shear_iso"`) finds the exact grip planes and runs.

- *With the unload* (`lmp -sf omp -pk omp 8`, 542k-atom random-gel `_x_plates` file): add `-var unload 1 -var unload_steps 2000`, and on the full-size file `-var strains 0.02 -var strains_list "0.02" -var nsteps_shear_base 2000 -var nsteps_shear 6000 -var seat_mode rg`. Pass: the log's `[unload 0.02] retraction done: gamma = ...` line near 0, the level's file set written again tagged `_2000_u0.02`, and `lib/shear.py` reading the level with `L['U']` and a `D_c` from `fit_Dc_unload` in every figure.
- *With the fine displacement file* (`lmp -sf omp -pk omp 6`, the full 435k-atom `_x_plates` file, 9.5 min): add `-var disp_fine_nfreq 100 -var disp_fine_span 700 -var wall_check_every 100`. Pass: fine files with one block per 100 steps from the drive start to hold + 700, both hold segments, the coarse files unchanged (10 blocks), no FENE warnings.

**`simple_shear`, periodic deck** (about 4 min on the 290k-atom file): the isolated line without the core margins and with `-var nsteps_shear_base 30000`, run as `-in simple_shear_periodic.lmp`.

### 5c. Editing the batch file

The section to edit is below the `source` line:

```bash
FOLDER=uniaxial_compression
DATANAME=final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000002_two_pist
INTERACTION=1.0_1.0
NSTEPS=9000000          # fallback hold with HOLD_AUTO=0, and the tag of the setup-phase files
export PRESS_TARGET=1.5
export COMPRESSIONS="0.10 0.15 0.20 0.3 0.4 0.5 0.7"
```

- **`DATANAME`**: the filename of the `.data` file in `lammps_data/input_data/`, without the `.data` extension. A driver can point `run_lammps.sh` at another folder with `LAMMPS_INPUT_DATA_OVERRIDE`; `calibration_sweep.sh` uses `input_data/calibration_sweep/` for its ~370 prep files and run symlinks. Several batches read `DATANAME` from the environment, so one batch serves every gel: `DATANAME=<stem> sbatch <batch>`.
- **`INTERACTION`**: `epsSS_epsSP`, the solvent–solvent and polymer–solvent LJ well depths (see cheat sheet).
- **`NSTEPS`**: number of timesteps (1 million steps ≈ 4100 τ ≈ 19 ns for PEG). It is also the output-file tag. In the compression decks it is the flat-hold length (with `HOLD_AUTO=1` the holds are auto-sized per level and `NSTEPS` is the fallback), in permeation the drive length, in `simple_shear` the hold per strain level, and under the Rg halt a budget.
- **Deck knobs** are exported environment variables named like the deck's LAMMPS variable in upper case (`PRESS_TARGET` → `-var press_target`). `run_lammps.sh` forwards a knob only when it is set, so a knob left out keeps the deck's own default. The full list is `KNOBS` in `scripts/lib/run_common.sh`, and the decks' `variable … index` lines give the defaults.
- **Lists**: `STRAINS` and `COMPRESSIONS` are space-separated lists.
- **Run control**: `CHAIN=1` chains continuation jobs (§5e), `WALL_MARGIN` sets the deck's stopping margin, `VEL_SEED` sets the RNG seed (default 12345), and `COMM_CUTOFF=5.0` restores the pre-2026-11-07 communication cutoff.
- **`DECK_FILE`** (optional): run `simulations/$FOLDER/$DECK_FILE.lmp` instead of `$FOLDER.lmp`, for a folder that keeps more than one deck (the `simple_shear` batches export `DECK_FILE=simple_shear_periodic` or `simple_shear_isolated`). Output folders and file names still key on `FOLDER`.

For `simple_shear` the other lengths are index knobs: `phase0_steps` = 50k, `npt_piston_steps` = `NPT_PISTON_STEPS`, `ref_avg_steps` = 150k, and the drive cap `nsteps_shear` = 570k (every drive halts at its strain target long before).

To extend a finished run, use `continue_sim.sh` ([§5e below](#5e-continuing-a-run-with-continue_simsh)). Resubmitting the batch with a larger `NSTEPS` is a fresh job that redoes the setup.

### 5d. Submitting a job

```bash
# SSH into the cluster, then:
cd ~/Documents/lammps_work/simulations/uniaxial_compression
sbatch uniaxial_compression_sweep.batch
# every folder keeps <folder>_single.batch (one level / one run) and, where there is a sweep, <folder>_sweep.batch;
# a knob can be set for one submission on the command line: DATANAME=<stem> sbatch uniaxial_compression_single.batch
#
# From any folder, on either cluster: lsbatch finds the batch file by name and submits it from its own folder.
# On Pod it is the only way in (it reads the batch file's #POD line); on Expanse it is the same as the sbatch above.
lsbatch uniaxial_compression_sweep.batch
DATANAME=<stem> lsbatch --nodes=4 uniaxial_compression_single.batch
```

`lsbatch` is a shell alias for `scripts/submit.sh`, installed with `lsync` by `bash install.sh` (once per machine, then a new terminal). Each batch file holds two allocations: the `#SBATCH` lines are Expanse's and the `#POD` line under them is Pod's. `lsbatch` submits with the one of the cluster it runs on. The settings behind each cluster (binary, modules, MPI, scratch) are the two blocks of `scripts/lib/site.sh`. See [cluster_reference.md](cluster_reference.md), Pod.

`sbatch` submits the job script to SLURM's queue and returns a job ID at once; SLURM starts the job when the requested nodes are free.

**Check the job status:**
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
tail -f ~/Documents/lammps_runs/uniaxial_compression/uniaxial_compression_*/log.lammps
```

### 5e. Continuing a run with `continue_sim.sh`

`continue_sim.sh` continues a finished run without restart files or batch edits. It reads the SLURM output file to find the original working directory and takes every run parameter from there. It is a restart in the proper sense (setup skipped, production continued), unlike editing `NSTEPS` in a `.batch` file and resubmitting, which is a fresh job that redoes the setup ([§5c above](#5c-editing-the-batch-file)).

**Output files from before 2026-11-18.** The run folders under `lammps_runs/` were moved into the renamed parents that day ([repository_layout.md](repository_layout.md#moved-on-the-cluster-2026-11-18)), and an output file written before it still names the old parent on its `Working directory:` line. When that path is gone, the script looks the run up under the new parent (`RENAMED_RUN_PARENT` in the script) and prints a `NOTE:` with the path it uses. This only finds the folder. The runs moved into `pure_solvent/` and `pure_polymer/` wrote `final_config_*.data`, not the `puresolv_*` / `purepol_*` files the current decks write and the script resumes from, so they stop at "no `puresolv_*.data`"; the other renamed folders have no continuation path at all. `slab_with_support` and `chain_length_testing` were not renamed and are unaffected.

**When to use it:** to add steps to a completed run, or to finish one that hit its wall-clock budget (`--chain`, below).

- Supported folders (the ones whose decks have a `cont=1` path): `slab_with_support`, `chain_length_testing`, `pure_solvent`, `pure_polymer`.
- Not for the two-piston decks (`uniaxial_*`, `isotropic_compression`, `simple_shear`: always fresh starts) or the phase sweeps (their internal P-sweeps complete in one invocation).
- The archived one-piston and `simple_shear_plates` decks are not supported; the 2026-11-07 version of the script in git history handles them if one is ever revived.

#### What it does per folder

The script skips the setup and resumes the folder's production behaviour.

| Folder | Skipped | Runs |
|--------|---------|------|
| `slab_with_support`, `chain_length_testing` | Soft push-off, minimize, NVT ramp, NPT (`aniso`) warm-up | `aniso` NPT production with volume/dimension outputs |
| `pure_solvent`, `pure_polymer` | Box rescale/harmonic pre-relax, gentle Langevin ramp | More NPT production steps |

The same output files are produced (stress profiles, chemical potential, piston/permeate data, trajectories, `log.lammps`) in a fresh `continuation_{timestamp}/` subfolder, so nothing from the original run is overwritten. The only omissions are the setup trajectory and any ε=0/zero-flux reference files, which already exist from the original run. The continuation job forwards every deck knob (`BARO_MODE`, `PIN_GEL`, `PISTON_TRANSPARENT`, `PRESS_TARGET`, the Rg-halt knobs and the rest) and post-processes through `postprocess.sh` like `run_lammps.sh`.

#### How to run it

1. **Go to the simulation folder** where the original job was submitted:
   ```bash
   cd ~/Documents/lammps_work/simulations/slab_with_support
   ```

2. **Run the continuation** with the SLURM job ID from the output file name and the number of extra steps:
   ```bash
   ~/Documents/lammps_work/scripts/continue_sim.sh <job_id> <nsteps>
   ```
   For example, if the output file is `slab_support.o49772594.exp-14-05`:
   ```bash
   ~/Documents/lammps_work/scripts/continue_sim.sh 49772594 500000
   ```
   No other arguments are needed: the current folder determines the data-file naming. With `~/Documents/lammps_work/scripts` on `$PATH` in `~/.bashrc`, the command is `continue_sim.sh 49772594 500000`.

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

`continue_sim.sh <job_id> <nsteps>` adds `<nsteps>` on top of a finished run (outputs tagged `<nsteps>`, `oldsteps = 0`). When it finds no final file it falls back to the newest snapshot and says so.

#### What the script does

1. Finds `*.o{job_id}.*` in the current directory (the SLURM output file).
2. Reads the line `Working directory: /path/...` that `run_lammps.sh` printed when the job ran. This gives the original output folder.
3. Looks up the current folder's output-file prefix (`final_config`, `puresolv` or `purepol`) and finds `<prefix>_*.data` in that folder's `data_files/` (runs before 2026-11-17: in the folder itself).
4. Parses `dataname`, `epsSS`, `epsSP` from the filename (last three `_`-delimited tokens).
5. Borrows the `#SBATCH` lines of the `.batch` whose `--job-name` matches the output file (else the batch a chained job inherited, else `<folder>_single.batch`, else the folder's only batch; `--batch <file>` overrides) and every deck knob: the environment first, then the batch's assignment lines.
6. Creates `continuation_{timestamp}/` inside the original folder, copies the deck and symlinks the data file in.
7. Passes `-var cont 1` to LAMMPS; the `jump`/`if` logic in the deck then skips the setup and carries on the production block. With `cont = 1` the deck declares `pair_style`/`bond_style` before `read_data`, because a `write_data` snapshot carries `Pair Coeffs`/`Bond Coeffs` sections that `read_data` refuses otherwise.

#### Runs longer than the wall limit: automatic chaining

Expanse's `compute` partition allows 48 h per job. A 14 M-step `slab_with_support` run of the 1.13 M-atom `_nowalls` file needed 55–60 h on 4 nodes at the pre-2026-11-08 layout (the 238 k-atom periodic slab ran at 78 Matom-step/s there). Such a job chains itself instead of being split by hand:

1. `run_lammps.sh` passes the deck a wall-clock budget, `-var wall_timeout` = seconds left in the SLURM job minus `WALL_MARGIN` (default 2400 s; an env knob). `slab_with_support.lmp` arms it with `timer timeout … every 1000` right at start-up, so setup time counts too.
2. When the budget runs out, the production block in flight returns early (`WARNING: Wall time limit reached`), its snapshot `data_files/final_config_*_snap<N>.data` is written as usual, the loop exits, and the deck logs one bookkeeping line, `>>> PRODUCTION_STEPS_DONE <done> OF <nsteps> OLDSTEPS <o> TOTSTEPS <t> STATUS INCOMPLETE`. It does **not** write the un-suffixed `final_config_*.data`, which stays the completion marker.
3. With `CHAIN=1` exported by the `.batch`, `run_lammps.sh` then runs `continue_sim.sh --chain <job_id>` before post-processing. That reads the bookkeeping line, resumes from the newest snapshot with `cont=1` (no push-off/minimize/warm-up), asks for exactly `totsteps − (oldsteps + done)` more steps and keeps the original `totsteps` tag on every output file. The new job (`cont_<job_id>`, borrowing the `#SBATCH` lines and `export`s of the `.batch` whose `--job-name` matches the output file) chains again if it also runs out of time.
4. The final `final_config_<dataname>_<interaction>_<totsteps>.data` is in `data_files/` of the **last** `continuation_<timestamp>/` directory under the original run directory. Continuations of continuations are siblings there, never nested. Each job's own `continuation_*/log.lammps` and `output_files/` cover its share of the steps.

The margin covers the final `write_data` of a ~1 M-atom system, the MPI teardown grace, the chain submission and the plots.

- If the chain breaks (node failure, `sbatch` refused), re-launch it from the simulation folder with `continue_sim.sh --chain <last job id>`.
- `--dry-run` prints the batch script it would submit; `--batch <file>` picks another `.batch` to borrow SLURM settings from.
- `--chain` on a job whose log says `STATUS COMPLETE` submits nothing.
- Decks without the chaining path ignore the variable, and without the bookkeeping line nothing is chained. (`simple_shear_isolated.lmp` uses `wall_timeout` only to end its sweep cleanly, §5b.)

### 5f. SLURM resource guidelines

| Cluster | Partition (CPU) | Cores/node | Optimal ntasks/node | Max walltime |
|---------|----------------|------------|---------------------|--------------|
| Bridges-2 | `RM` | 128 | **120** (128 causes comm overhead) | 5 days |
| Expanse | `compute` | 128 | 128 | 2 days |
| Pod | default | 40 | 40 (not re-measured since 2026-10) | `sinfo -o "%P %l"` |

#### Pod per-user limits (queried 2026-11-19)

Pod's cap is the QOS, not the account: `pod-users` carries no limits of its own, so each partition's QOS is the whole story. The `batch` partition (the default; 62 nodes x 40 cores = 2480 CPUs, default walltime 3 days, max 58 days) runs under the `normal` QOS:

| QOS | Partitions | Max CPUs per user | Max running jobs | Max queued jobs | Partition walltime |
|-----|------------|-------------------|------------------|-----------------|--------------------|
| `normal` | `batch`, `normal` | **320** (8 full nodes) | 32 | 600 | 58 days (default 3) |
| `short` | `short` (1 node) | none | none | none | 2 h |
| `largemem` | `largemem` (4 nodes) | none | none | none | 37.5 days |
| `gpu` | `gpu` (14 nodes) | none | 8 | none | 10 days |

So at most eight 40-core jobs run at once per user; the rest sit with reason `QOSMaxCpuPerUserLimit` (a `Priority` reason means waiting on other users, not the cap). Size a sweep or node screen so no more than 8 x 40-core jobs need to be concurrent. Re-check with:

```bash
sacctmgr show qos format=Name,MaxTRESPU,MaxJobsPU,MaxSubmitPU
sacctmgr show assoc user=$USER format=Account,Partition,QOS,GrpTRES,MaxTRESPU
```

#### Pod node health (2026-11-20)

On Pod `lsbatch` adds an `#SBATCH --exclude` line for the nodes the node screen marks as suspect (throttling, slow) or still to be re-measured, and prints it with the header. The list is refreshed by a new screen when the last one is older than 7 days (1-core jobs) or 30 days (whole-node runs); `lsbatch` starts that screen itself, in the background, and the screen's own summary job updates the list when it is done. Nothing changes for Expanse or in the batch files. Details, overrides and how to read the table: [cluster_reference.md](cluster_reference.md), "Node health".

**Expected runtimes** (4M timesteps, 400k–600k atoms):
- 1 node (120 cores, Bridges-2): ~10–11 hours
- 4 nodes (480 cores, Bridges-2): ~3 hours
- 2 nodes (256 cores, Expanse): ~5–6 hours

Benchmarks across system sizes are in [`archive/pod_bridges/documenting_pod_runs.md`](../archive/pod_bridges/documenting_pod_runs.md).

#### Expanse allocation benchmark (2026-11-07)

Folder `simulations/allocation_benchmark/`, rot45 slab, N = 5 and 15. The debug-partition layout scan found, for both N:

- 128 ranks x 1 thread per node beats every hybrid split (64x2 by 7-30 %, 32x4 and 16x8 by 2-3x).
- `--map-by ppr:128:node:pe=1` (ranks fill a node before the next one, "core") instead of `--map-by node` gives 1.4-1.6x on 2 nodes.
- A comm cutoff of 2.5 instead of 5.0 gives another 1.4x with identical step-0 forces.

The compute-partition node scan (1-16 nodes, 128 x 1, core, cutoff 2.5) then gave, against the 2026-11-06 production setup (4 nodes, 128 x 1, map-by node, cutoff 5.0):

| N | nodes | steps/s | SU per M steps | speed vs baseline | SU vs baseline |
|---|---|---|---|---|---|
| 5 | 1 | 362 | 98 | 1.5x | 0.16x |
| 5 | 4 | 710 | 200 | 3.0x | 0.33x |
| 5 | 6 | 893 | 239 | 3.8x | 0.40x |
| 15 | 1 | 189 | 188 | 1.2x | 0.21x |
| 15 | 4 | 368 | 387 | 2.3x | 0.43x |
| 15 | 6 | 460 | 463 | 2.9x | 0.52x |
| 15 | 8 | 501 | 567 | 3.2x | 0.63x |

Repeat spread at 3-6 nodes is 20-30 % (node placement), so 4 and 6 nodes tie in speed at N = 5.

**Adopted 2026-11-08:** 6 nodes x 128 x 1 in every production batch, `MPI_MAP_BY=core` in `scripts/lib/site.sh`, cutoff 2.5 in the decks. This is 2.9-3.8x the 2026-11-06 production speed at 40-52 % of its SU per step. 4 nodes is about 20 % cheaper per step and 25 % slower, for when SUs matter more than wall time (`sbatch --nodes=4 …`); 1 node is the cheapest per step. Data: `~/Documents/lammps_runs/allocation_benchmark/` (`nodescan_rot45.{png,txt}`, `layouts_rot45.{png,txt}`), README in the benchmark folder. The two-piston batches ran 32 x 4 until 2026-11-08 (0.55-0.75 M steps/h); `--time` budgets that were set for that layout now have slack.

### 5g. Volume-of-mixing sweep (`volmix_sweep.sh`)

Archived in 2026-08 and superseded by the calibration sweep (partial-molar-volume NPT sweeps; see `lammps_work/archive/README.md` and `simulations/calibration_sweep/README.md`). The script is `lammps_work/archive/volmix_sweep.sh`; the calibration driver inherits its SLURM chaining machinery. Its data files were deleted on 2026-11-07: `lammps_data/slab_with_support/` is gone, as are the old `pstar*`/`snap*`/`rep*` `_polymer_only` / `_solvent_only` / `_volmix_manifest` files in `input_data/`. The current calibration inputs are in `input_data/calibration_sweep/`.

The sweep computed ΔV_mix(P*) = V_gel_mixed − V_pure_polymer − V_pure_solvent across 11 pressures (P* = 1.0–2.0, step 0.1), one pressure at a time, as a SLURM chain: a 600k-step `slab_with_support` equilibration, `isolate_gel.py` (control volume = the inner 88 % of the polymer distribution, `--cv-percentile` default 6 %), `split_gel.py`, then 100k-step `pure_solvent` and `pure_polymer` runs of the split components. Analysis: `volume_of_mixing.ipynb` ([§7b, Jupyter analysis notebooks](analysis.md#7b-jupyter-analysis-notebooks)).

---

### Pipelines

#### Chain-length sweep

`cd simulations/chain_length_testing && bash submit_chain_length_sweep.sh [N ...]` submits one job per strand length, each running `chain_length_testing_sweep.batch` with `NBEADS` and `NETWORK` exported and named `chain_len_<N>b[_<network>]`. The script checks that every input file exists in `input_data/` before submitting. The input files come from `scripts/lib/build_beads_sweep.py` ([building_data_files.md](building_data_files.md)).

| `NETWORK` | Default N | Step budget (`NSTEPS`) | Rg block |
|---|---|---|---|
| `random` (default) | 10 15 20 30 | 24, 30, 35, 45 M (the 5-bead slab's 14 M scaled by the mesh ratio) | 1.5, 2, 2.5, 3 M |
| `rot45` | 5 15 | 20 M, 60 M | 1 M, 2 M |

The random budgets were 21–36 h on 4 nodes at the pre-2026-11-08 layout. Each job has 36 h and `CHAIN=1`.

**Rg-plateau halt** (`slab_with_support.lmp`, off unless `RG_HALT=1`). `NSTEPS` is only the budget; production ends on its own once the gel's block-mean R_g has stopped moving.

- Production runs in blocks of `RG_BLOCK` steps. `fix ave/time` averages the polymer `compute gyration` over each block (one sample every `RG_EVERY` = 10000 steps), and after each block the relative change of the block mean from the previous block's is compared with `RG_TOL`.
- `RG_NCONSEC` consecutive blocks under the tolerance print `RG_PLATEAU_DETECTED`. The deck then runs a tail of `RG_EXTRA` × the production steps done so far (whole blocks, at least one), prints `RG_PLATEAU_CONVERGED`, leaves the loop and reports `STATUS COMPLETE`, so the `final_config` is written and `CHAIN=1` stops.
- The settings (`RG_TOL` 0.15 %, `RG_NCONSEC` 3, `RG_EXTRA` 10 %) were set on the 5-bead random run's R_g: after 10 M steps its block-to-block change was 0.02–0.14 %, the plateau is detected at 13 M and the tail ends the run at 15 M, 5–10 % past where its convergence plot goes flat.
- Longer strands get longer blocks rather than a looser tolerance, because a slow drift inside a short block would pass as a plateau.
- A continuation job (wall limit first) restarts the block count, so a chain hop costs up to `RG_NCONSEC + 1` extra blocks. `continue_sim.sh` carries the knobs through the job environment.
- The `<totsteps>` tag of every output file stays the budget. The log's `RG_BLOCK` lines (step, block mean, previous, relative change) and `PRODUCTION_STEPS_DONE` say how far the run went; `scripts/chain_length_comparison.ipynb` reads them.

Then, per N: pull `final_config_<DATANAME>_1.0_1.0_<budget>.data` to `input_data/`, convert it with `lib/slab_two_pistons_periodic.py`, run `uniaxial_compression` on it (one level, no unload for now) and point the notebook's `CFG[N]` at the run.

#### Random periodic slab to two-piston compression

`bash scripts/submit_random_periodic_pipeline.sh` (once, from the login node) does three things:

1. Builds `slab_support_periodic_5beads_tall_rho04_random.data` into `input_data/` if it is missing (`lib/periodic_data_file_generator.py --network random`).
2. Submits `simulations/slab_with_support/slab_with_support_periodic_random.batch`.
3. Submits `simulations/uniaxial_compression/chain_random_periodic.batch` with `--dependency=afterany` on it.

The bridge runs on one `shared` core when the slab job ends. It finds the `final_config_…_14000000.data` (continuation subdirs included; if the slab chained a continuation it re-submits itself after that job), copies it to `input_data/`, converts it with `lib/slab_two_pistons_periodic.py` (`…_two_pist.data`) and submits `uniaxial_compression_random.batch`. Outputs: `slab_periodic_random.o<job>.*`, `chain_random_periodic.o<job>.*`, `uniax_comp_random.o<job>.*` in their simulation folders.

#### The rot45 standard gel

The periodic rot45 lattice slab ([building_data_files.md](building_data_files.md): 5 and 15 beads per strand, 15 % taller than the 14000002 slab) is the standard gel of the compression + permeation programme since 2026-11-07. The files are in `lammps_data_files_local/` (`slab_support_periodic_{5,15}beads_tall_rho04_rot45.data`). Copy them to Expanse `~/Documents/lammps_data/input_data/`, then from the login node:

```bash
cd ~/Documents/lammps_work/simulations/chain_length_testing
NETWORK=rot45 bash submit_chain_length_sweep.sh 5      # -> job J5  (chain_len_5b_rot45;  budget 20 M, Rg halt)
NETWORK=rot45 bash submit_chain_length_sweep.sh 15     # -> job J15 (chain_len_15b_rot45; budget 60 M)
cd ../uniaxial_compression
# 5 beads: convert and submit BOTH two-piston runs when the equilibration ends
sbatch --dependency=afterany:J5 --export=ALL,SLAB_JOB=J5,SLAB_FOLDER=chain_length_testing,SLAB_DATANAME=slab_support_periodic_5beads_tall_rho04_rot45,SLAB_NSTEPS=20000000,MARGIN_PERM=40,MARGIN_FEED=30,SUBMIT="compression permeation" chain_random_periodic.batch
# 15 beads: the compression first -- its M sets the permeation dP (header of the permeation batch)
sbatch --dependency=afterany:J15 --export=ALL,SLAB_JOB=J15,SLAB_FOLDER=chain_length_testing,SLAB_DATANAME=slab_support_periodic_15beads_tall_rho04_rot45,SLAB_NSTEPS=60000000,MARGIN_PERM=40,MARGIN_FEED=30,SUBMIT=compression chain_random_periodic.batch
```

`SLAB_NSTEPS` is the equilibration's budget, the tag of its `final_config` file.

**The bridge** (`chain_random_periodic.batch`, knobs `SLAB_*` / `MARGIN_PERM` / `MARGIN_FEED` / `SUBMIT`; it follows continuation jobs) copies `final_config_<DATANAME>_1.0_1.0_<SLAB_NSTEPS>.data` to `input_data/`, converts it once per gel (`slab_two_pistons_periodic.py --margin-perm 40 --margin-feed 30` → `<stem>_two_pist_pm40.data`) and submits `uniaxial_compression_single.batch` and/or `uniaxial_permeation_single.batch` with `DATANAME` exported. By hand it is the same three commands (the converter line is in the chain-length batch header; then `DATANAME=<stem> sbatch <batch>`).

**One file for three decks.** The compression, permeation and periodic shear decks all read the `_two_pist_pm40` file. The permeate margin is vacuum and costs the compression nothing, and the shear deck only needs the 7-type two-piston layout.

**`--margin-feed 30` is required for the rot45 gels.** The converter's default 15 σ feed-side vacuum leaves 13 σ of headroom above the feed piston, and a 0.10 strain on these taller gels (L_bb 148 / 155 σ) expels about 15 σ of solvent into the feed reservoir. The deck's `box_guard` (feed piston within 2 σ of zhi) would halt the drive at ε ≈ 0.08–0.09; the converter prints the cap (`absorbs the solvent expelled up to a compression strain of ~0.08`). 30 σ gives a cap of about 0.18.

**Current files.** The 5-bead slab ended at 20 M steps and the 15-bead at 60 M, so the stems are `…_5beads_…_rot45_1.0_1.0_20000000` and `…_15beads_…_rot45_1.0_1.0_60000000`. The converted files have 425 k and 831 k atoms (conversion on the login node: 21 s and 1 min).

**The 15-bead gel** has 1.9× the atoms (its cross-section matches the polymer beads per z-bin). Submit its two-piston runs with `sbatch --nodes=10 …` (the command-line option overrides the script's `#SBATCH --nodes`) or halve its hold / unload / drive budgets; the full length on 6 nodes is untimed, so also add `--time=48:00:00`. For its permeation, set `DP_PISTON` from the compression M so that dP/(2M) ≈ 0.12 (dP ≈ 0.24 M): `DATANAME=<15-bead stem> DP_PISTON=<value> sbatch uniaxial_permeation_single.batch`.

**Measured wall times (2026-11-11).** Compression with unload (25.2 M steps): 15.8 h for the 5-bead gel on 6 nodes (job 55663497) and 18.4 h for the 15-bead gel on 10 nodes (55663502). The periodic shear times are in §5b. Measured 6-node equilibration rates: 2.1–3.1 M steps/h for the 5-bead slab (354 k atoms), 1.75 M/h for the 15-bead (660 k). The batch limits (28 h compression, 12 h periodic shear) are ≥ 1.5× the slower run; neither deck has a `wall_timeout`, so longer holds or extra strain levels need `sbatch --time=…`.

**Permeation budgets.** The permeation batch's 48 h was sized at the old layout's ~0.55 M steps/h (about 31 h, about 47 h with the optional unload; the frozen reference run about 11 h) and has not been re-timed. The random slab ran at 0.75 M/h on 288 k atoms there and lost the end of its 25 M unload to the 48 h limit.

**Analysis.** `uniaxial_compression_single.ipynb` (LEVEL 0.10, the unload figures) and `uniaxial_permeation_single.ipynb` with the new `DATANAME` and a `RUN_ID` of your choice, then `manuscript_figures.ipynb`.

#### Isolated cube

The isolated-cube pipeline on Expanse, in order; every step reads `DATANAME` from the environment.

1. Build the cube with `scripts/data_file_generator.ipynb` (`NETWORK = "rot45"` for the rotated lattice; `gen.build` writes an isolated gel, no sheets) and put the file in `~/Documents/lammps_data/input_data/<gel>.data`. For a cube give `units_x`, `units_y`, `units_z` of equal edge length: the rotated builds count BCT cells of a/√2 along x and z, so x and z take about 1.41× the y count.
2. `DATANAME=<gel> sbatch slab_with_support_single.batch` (free swelling at P* = 1.5, `BARO_MODE=iso PIN_GEL=1`). The run dir holds `final_config_<gel>_1.0_1.0_<NSTEPS>.data`.
3. `python scripts/lib/add_more_plates_to_gel.py --layout six_plates --input <that final_config>`, and copy the `…_six_plates.data` it writes to `input_data/`.
4. `DATANAME=<final_config stem>_six_plates sbatch isotropic_compression_single.batch`.

`bash scripts/submit_cube_pipeline.sh [NETWORK] [N] [BEADS]` (defaults `random 10 10` since 2026-11-22: the cube twin of the wide random 10-bead slab, about 78 σ per edge and 1.2 M atoms; `rot45 12 5` is the earlier lattice cube) runs all four steps in one submission, with `chain_cube_six_plates.batch` as the bridge job. Strands longer than 5 beads are laid as random walks (`chain_placement='walk'`, now also in `data_file_generator.py`), the bath starts at 0.44 with wider pads for the surface fringe, and the equilibration has a 16 M-step budget with the Rg-plateau halt. The compression is the single batch's defaults: Rg seating, one 0.10 level, the unload with fine displacement profiles. K comes from the hold and D_m only from the Phase 2c re-swelling (`UNLOAD=1`, set in both batches; `UNLOAD_STEPS=0` makes it as long as the hold). Analyse with `isotropic_compression_single.ipynb`, `DATANAME` = the six-plate stem, `LEVEL "0.10"`.
