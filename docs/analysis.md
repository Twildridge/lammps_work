# Analysis

[← back to README](../README.md)

### 7a. Python scripts (auto-run after each job)

`scripts/postprocess.sh` runs these plotters at the end of every job. To replot after changing a script, run them by hand from inside the timestamped run directory.

#### Run directory and run_id

Run directories live at `~/Documents/lammps_runs/` and are named:
```
<folder>_<DATANAME>_<INTERACTION>_<TIMESTAMP>/
```

The run_id in every output file name is `<DATANAME>_<INTERACTION>_<TOTSTEPS>`, where `TOTSTEPS = NSTEPS` from the batch file ([§5c, editing the batch file](running_simulations.md#5c-editing-the-batch-file)). It can be read off the output files:
```bash
ls ~/Documents/lammps_runs/<run_dir>/output_files/stress_data/
# e.g.: stress_tensor_polymer_isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000.dat
#                              └─────────────────────── run_id ───────────────────────────────────────────┘
```

Run directories made since 2026-11-17 carry a `README.txt` (what every subfolder and file-name tag means) and a `run_info.txt` (job, deck, input file, every LAMMPS variable). The layout is described in [repository_layout.md](repository_layout.md#inside-a-run-folder).

#### Replotting by hand

Change into the run directory first. Every script takes `.` as the folder argument and the run_id as the second argument.

```bash
cd ~/Documents/lammps_runs/<run_dir>
```

**`plot_lammps_log.py`** plots T, P and volume convergence for every simulation type, plus the shear diagnostics for `simple_shear_plates` runs (detected automatically). Plots go to `./output_plots/convergence_plots/`.

```bash
python ~/Documents/lammps_work/scripts/lib/plot_lammps_log.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000
```

This is the call `run_lammps.sh` makes. For `simple_shear_plates` the base `DATANAME` can be passed as the title and the full run_id given through `--run-id` for the file lookups:

```bash
python ~/Documents/lammps_work/scripts/lib/plot_lammps_log.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000 \
    --run-id isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000
```

**`plot_stress_profiles.py`** plots the partial stress and volume fraction profiles (`slab_with_support`, `uniaxial_permeation` and similar folders). The third argument is `OLDSTEPS`, which is always 0. It overlays every `stress_z_piston_<sheet>` file in the Piston column.

```bash
python ~/Documents/lammps_work/scripts/lib/plot_stress_profiles.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000 \
    0
```

**`plot_piston_data.py`** plots piston position and velocity (`uniaxial_compression`, `uniaxial_permeation` and similar folders). Two-piston files (`# step z_feed z_perm` headers) get one line per piston, a measured-vs-applied pressure panel and a `Q_perm` / solvent-expelled panel; one-piston files are plotted with a single line. For two-piston runs `plot_lammps_log.py` also adds the bath-check and `Q_perm` panels to its flow diagnostics.

```bash
python ~/Documents/lammps_work/scripts/lib/plot_piston_data.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000 \
    0
```

Compression sweeps and `simple_shear` runs use their own consolidated plotters (`plot_compression_strain_sweep.py`, `plot_shear_strain_sweep.py`); see [§5b, simulation types](running_simulations.md#5b-simulation-types).

---

### 7b. Jupyter analysis notebooks

Open the notebooks on the MacBook in JupyterLab (`jupyter lab`). They read data files in `flow_data_local/raw_data/<sim_type>/<RUN_ID>/`. Each notebook has a config cell near the top; only `RUN_ID` and `sim_name` change between runs, and all paths derive from those. In the uniaxial-compression notebooks leave `NSTEPS = None`: the `<steps>` tag in the file names is each level's hold length and is resolved from the files.

| Notebook (`scripts/`) | Library (`scripts/lib/`) | What it produces |
|---|---|---|
| `uniaxial_compression_single.ipynb`, `uniaxial_compression_sweep.ipynb` | `uniaxial.py` (`volfrac.py`) | Two-piston compression, one level or a whole sweep: M, G, the unload D_c, κ = D_c/M, stresses, volume fractions |
| `uniaxial_permeation_single.ipynb` | `uniaxial.py` (`volfrac.py`, `psd.py`) | Two-piston permeation: D_c and M from the polymer displacement, κ from the velocity-controlled and frozen stages, their reconciliation |
| `isotropic_compression_single.ipynb`, `isotropic_compression_sweep.ipynb` | `isotropic.py` | Bulk modulus K from the hold, D_m = (K/M) D_c and κ = D_m/K from the free re-swelling unload |
| `simple_shear_single.ipynb`, `simple_shear_sweep.ipynb` | `shear.py` | Shear modulus G, N1/N2, unload D_c and κ = D_c/G |
| `manuscript_figures.ipynb` | `manuscript.py`, `compare.py`, `chain_length.py` | The paper's figures, each editable in its cell; then the figures that combine decks: D_c(φ_s), κ(φ_s), the equation of state σ′(φ_s) |
| `chain_length_comparison.ipynb`, `calibration_analysis.ipynb` | `chain_length.py` (chain-length notebook) | Beads-per-strand sweep and the φ_f calibration; described in `scripts/README.md` |

There is no permeation sweep notebook.

#### Conventions shared by the notebooks

Every notebook has the layout Config → sync → load/compute → figures (one call per figure) → Notes (a markdown cell with the method notes). All analysis code is in the library module. After editing a module, re-run the first notebook cell (`importlib.reload`); no kernel restart is needed.

`lib/uniaxial.py` holds the readers, the Terzaghi network/pore split, the drift-tested block-bootstrap plateau window behind `M_piston`, the G estimate, the unload D_c fit, the λ-calibrated Voronoi φ_s (through `lib/volfrac.py`), the cluster sync and every figure. Its knobs live in `uni.Config`, with the defaults documented there; the ones that matter for M, G, D_c and κ are spelled out in each notebook's Config cell. `lib/shear.py` and `lib/isotropic.py` import the statistics, readers, palette and legend placement from it, so the shear and bulk numbers are formed with the same code as the compression numbers.

**Sync from Expanse or Pod.** The sync cell of the uniaxial, shear and isotropic notebooks (`sync_from_expanse` → `uni.sync_pull`) pulls from the cluster named by `CLUSTER` in the Config: `'expanse'` (password + TOTP), `'pod'` (password; the campus VPN has to be connected) or `'auto'`, the default. `'auto'` asks the cluster this run was last pulled from (remembered in `<DATA_DIR>/.sync_cluster`), else Expanse, and then the other one for whatever is still missing; Pod is asked only when its login node answers, so with the VPN down it is skipped with a note and the missing files are not recorded as absent. Hosts, users and paths are in `CLUSTERS` in `lib/uniaxial.py`; `RUNS_ROOT` / `TRAJ_ROOT` stay the Expanse paths and are translated for Pod (`$HOME/Documents/lammps_runs`, `/scratch/$USER/lammps_trajectories`). The Pod path has been tested with a mocked login only, not yet against Pod. `chain_length_comparison.ipynb` and `calibration_analysis.ipynb` still pull from Expanse only.

Every printed number has at most 3 significant figures, and legends are placed by `uni.smart_legend` so that they do not cover data.

Two `uni.Config` knobs control the figure style of the three two-piston notebooks. Both are on by default. `iso.Config` sets both off, and the shear and cross-deck notebooks do not use these figures.

**`FLIP_Z = True`** draws every z-profile figure with the feed / load piston on the left and the permeate / support on the right, as in Marioni et al. The `z/L` axes become `1 − z/L`, running 0 → 1 with the two ends named in the label (feed → permeate, or load piston → support); the σ- and ζ-unit axes are inverted and say so in their labels. It affects plotting only (`zn`, `finish_axes`, `_zlim`, `flip_z_axis`); every mask, window and fit uses z itself.

**`PARTIAL_NORM = 'share'`** makes `fig_total_stress`, `fig_partial_stress`, `fig_thermo_pressure` and their `_sweep` versions draw a second figure under the original (same cell, same plot number, file stem `…_norm`). `'raw'` draws the originals only. The second figure removes each species' share of the back pressure and divides by the driving pressure (`uni.partial_norm`):

- `w_s = σ_s/σᵗ` per z-bin and snapshot (for the trace, `w_s = tr σ_s/tr σᵗ`)
- `σ_s* = (σ_s − w_s P_ref)/ΔP`
- `σ_p* = (σ_p − (1 − w_s) P_ref)/ΔP`
- `σᵗ* = (σᵗ − P_ref)/ΔP`

Under permeation `P_ref` and `ΔP` are the applied permeate pressure and the feed − permeate difference. In compression they are `P_BARO` and the level's plateau load-piston increment `dP_pist`. The shares come from the partial stresses alone, with no mass or volume fractions. This is the form of Marioni's partial-`P_zz` figure but it does not reproduce his end points (perm_3: solvent* 0.50 → 0.38 and polymer* 0.49 → 0.61 from the feed face to the support face); see the Notes cells.

The normalised trace figure (`fig_thermo_pressure_norm`) also carries the guide `ΔP_th = (4/3)(G/M) ΔP_ext`. This is the drop of the total `P_th` from the feed reservoir's value at the feed face to the permeate side when the lateral stresses carry `(M − 2G)/M` of the axial load; in a compression hold it is a flat `1 − 4G/3M`. One guide is drawn per G estimate: `G_REF` (the shear notebooks' 0.2) and `G_COMP_REF` (the compression holds' lateral-stress value, about 0.05). Both default to `None`. `M` is the run's own primary estimate, or `M_REF` if set. `calibration_analysis.ipynb` ends with this figure drawn from perm_3 through its own `uni.Config`, as the input to the `P_CAL_MODE` question.

#### `uniaxial_compression_single.ipynb` and `uniaxial_compression_sweep.ipynb`

The single notebook analyses one applied-strain level: a single-level run, or one level picked out of a sweep (`LEVEL = "0.15"`). The sweep notebook overlays every level (colour = level, reference dashed) and adds M, G, D_c and κ against applied strain and the stress–strain curve. `COMP_LEVELS` must match `STRAIN_TARGETS=(...)` in the batch file. Per-level numbers come from the same `load_level` in both notebooks, so the two cannot disagree.

**Inputs and sync.** `Config(mode="compression")` points the sync at `lammps_runs/uniaxial_compression` (`two_pist=True` is the default) and pulls the stress and density profiles, including the `sigmaxx_*` / `sigmayy_*` profiles, and the `piston_pressure`, `permeation` and `pressure_reservoirs` files. The pair/bond dumps are not needed. The load piston is "the piston": its force is the network load, and the piston files carry it in the first value column (`[load | feed | perm]`).

**Figures.**

- Strain diagnostic and the displacement profile (`fig_strain`, `fig_disp_profile[_sweep]`).
- Solvent volume fraction: mass fraction, Voronoi and λ-calibrated Voronoi, reference against compressed (`fig_volfrac`).
- Total σ_zz, σ_xx, σ_yy evolution normalised by the bath pressure `P_BARO` (1.5).
- Solvent and polymer partial stress with the total superimposed.
- Network stress σ′_zz, σ′_xx, σ′_yy from the Terzaghi split, in the reference and final equilibrated states only, with 95 % bands. No evolution curves are drawn because p_pore is not uniform while the gel consolidates, so subtracting one reservoir value is valid only once equilibrated.
- Reservoir normal-stress check (figure 5b, below).
- Piston pressure (linear and log).
- M, network against piston (`fig_M`; the sweep adds `fig_stress_strain_sweep`).
- The anisotropy σ′_zz/σ′_xx and σ′_zz/σ′_yy against step with propagated error bars (`fig_ratio`).
- G = (σ′_zz − σ′_ii)/2ε from xx and from yy, on the wall-trimmed membrane interior (`fig_G`).
- The unload D_c fit (`fig_Dc_reswelling`, `fig_Dc_sweep`) and κ = D_c/M (`fig_kappa[_sweep]`).
- The same profile fit against time (`fig_Dc_reswelling_time`, single notebook figure 10b): each snapshot projected on the fit's cosine modes, drawn as the first mode's share of the thickness change against the fitted decay, residuals underneath. It is a second view of the figure 10 fit, not a second estimate. The thickness-trace fits (bounding-box and Rg thickness against the uniform-strain series) were removed on 2026-11-21: the bounding box follows the surface fringe, the fixed asymptote fails when the gel does not return to the seated reference, and the free-asymptote Rg trace is a weaker fit of the same re-swelling event.
- Wet-piston bath check (`fig_wet_pistons[_sweep]`): P_feed and P_perm measured on the sheets against `P_target` over the hold.
- Solvent rate (`fig_solvent_rate`; `fig_strain_solvent` stacks it under the strain diagnostic): the rate at which solvent enters or leaves the gel over the level's whole load and unload, one curve. The gel is the slab between the planes below which the fractions `SOLV_Q` and 1 − `SOLV_Q` of the polymer mass lie (default 0.025); its solvent content changes by ρ_s,0·A times the change of their separation, read from the per-bin polymer counts of the displacement files. The Rg thickness is not used because it weights the mid-plane 1.5 and the faces 0 during a transient. No profiles are written during the ramp or the retraction, so those stretches are interpolated and drawn dashed. For one level, `fig_strain` draws the load and unload on one time axis with the phases shaded. Decks since 2026-11-21 also write `disp_z_polymer_drive` during the ramp and the retraction when `DISP_FINE_NFREQ` > 0 (every `disp_drive_nfreq` = 1000 steps), and the curve is then measured throughout; on older runs those two stretches are interpolated and drawn dashed. `fig_strain_solvent(..., layout='twin')` puts the strain (left axis) and the rate (right axis) on one panel with their zero lines aligned.
- The thermodynamic pressure P_th = −⅓ tr(σᵗ) evolution, and the osmotic pressure Π = −⅓ tr(σ′) in the reference and final states. Both are positive under compression in the profiles' sign convention.

The summary prints the plates over the level (support, load piston, feed and permeate pistons: start → end).

**Pore-pressure baseline.** The baseline is the feed-reservoir interior, because the box top is vacuum in this geometry. It is rebuilt per stress snapshot from the measured feed-piston plane (`piston_position`, column `z_feed`), keeping only bins that lie entirely at least `Config.res_wall_margin` (3 σ) clear of the piston sheet. The bin touching a wet piston under-reads σ_zz by about 0.07 (depletion layer plus the wall virial booked on the piston atoms), which shifts every σ′ by +0.016 if it is included. Results computed before 2026-09-22 carry that shift.

**Moduli.** `M_network` is the plateau-averaged interior network-stress increment from the ε = 0 reference divided by the strain. `M_piston` is the load-piston pressure over the drift-tested block-bootstrap plateau window, as an increment, divided by the strain. The strain is the network's own (`M_STRAIN = 'disp'`, the slope of the steady displacement profile u_z(Z)), and its CI is added in quadrature. `RELAX_SYS = True` fits an exponential tail to the hold and reports the unrelaxed stress as a systematic δ_sys on M.

**Reservoir normal-stress check (figure 5b, `fig_reservoir_normal_stress[_sweep]`).** The figure shows σᵗ_xx, σᵗ_yy and σᵗ_zz over the hold in the gel interior and in the two solvent reservoirs (per-snapshot piston-tracking masks, the pore-baseline windows), as plateau increments from the ε = 0 reference. The reservoirs are the in-situ control: a fluid must stay at P_bath in every component. A lateral rise confined to the gel is therefore network stress (σ′_xx; the box quotes G_comp/M = ½(1 − σ′_xx/σ′_zz)), while a lateral drift in the reservoirs would be a box or normalisation artefact.

The sweep version normalises the increments by the gel's Δσᵗ_zz per level and plots the implied G_comp/M against `G_REF`/M_net. On run 7 the reservoirs read 0 within ±0.003 at every level while the gel's lateral increment is 0.77–0.97 of the axial one, so the compression lateral stress is real network stress.

##### D_c from the unload (free re-swelling)

The compression D_c is the unload profile fit (`load_unload`, `fit_Dc_unload`; `DC_PRIMARY = 'unload'`). κ = D_c/M uses that D_c only, and a level without an unload has no D_c or κ. The hold contributes its geometry only (`hold_geometry`).

**Deck side.** With `UNLOAD=1` (`-var unload 1`, deck Phase 2c), after a level's hold both plates are driven back at the drive speeds to the seated reference gap (`UNLOAD_MARGIN` σ beyond it if set) and frozen, the displacement reference is reset, and the gel re-swells for `UNLOAD_STEPS` (0 = the hold length). The unload writes its own files tagged `_<tag>_u<lvl>`: fine and coarse profiles, bounding box, plate positions, load-piston force, reservoir pressures, solvent uptake, and σ_zz of polymer and solvent. The next level compresses from the reference again, so levels are matched compression / re-swelling pairs. `simulations/uniaxial_compression/uniaxial_compression_single.batch` is the one-level job on the full gel: level 0.10, an 8 M-step hold, a 14 M-step unload, and fine displacement profiles every 25 τ (`DISP_FINE_NFREQ=5000`).

**Boundary-value problem.** The unload starts from the uniform equilibrated compression, so there is no ramp-made skin. The plates no longer touch the gel, the network stress at both faces is zero, the strain there is pinned at the reference value (Dirichlet), and solvent crosses the faces. These are the permeation transient's boundary conditions, not the hold's. The modes are sin kπζ at k²π²D_c/L², so **τ₁ = L²/(π²D_c)**, four times the hold's τ₁ = L²/(4π²D_c) at the same D_c; the hold only redistributes solvent between the compressed skins and the interior, with no flux through the plates. The faces detach from the retracting plates within a few τ, so t = 0 of the trace is the retraction start.

**Two readings.**

- Profile fit (the trusted estimate, `L['unload']['prof']`): the displacement profiles since the reset, in the compressed gel's material coordinate with the COM drift removed, fitted with modes cos kπζ and free amplitudes.
- Thickness trace (a check): L_bb(t) − L_hold against the exact uniform-initial-condition series ΔL/ΔL∞ = 1 − Σ_odd 8/(kπ)² e^{−k²π²D_c t/L²}, with ΔL∞ = L_ref − L_hold fixed, which leaves one free parameter. The fit with ΔL∞ free is a further check.

`UNLOAD_L` picks the L of τ₁: `'hold'` (default) is the compressed thickness; `'ref'` reads (L_ref/L_hold)² larger. `load_level` prints both readings and `print_summary` has a `D_c_unload` column. `tests/dc_convention_test.py` integrates this boundary-value problem (u′ = s at both faces, mean u fixed) by finite differences, and both unload readings return the solver's D_c within 0.1 %.

##### Why the hold D_c is not used

`lib/uniaxial.py` no longer fits the compression hold. `fit_Dc`, `fit_Dc_stress`, the D_c(t) and D_c(φ_s) fits, `fig_Dc` and `print_hold_check` are in `archive/scripts/lib/uniaxial_load_Dc.py`, which runs against the current library (`ldc.add_load_Dc(cfg, R, L)`; the single notebook has an optional cell for it). `tests/dc_convention_test.py` imports `fit_Dc` from there.

The archived fit uses the modes of a slab held between two drained plates at the same bath pressure, cos/sin(2mπζ) at 4m²π²D_c/L², τ₁ = L²/(4π²D_c). (Terzaghi's constant-load series at k²π²D_c/L² is the wrong set for this problem and returns about 5× the true D_c.) `fit_Dc` fits the first `DC_FIT_TAU1` (6) τ₁ of the hold and flags a window or mode-count disagreement above `DC_STAB_FLAG` (20 %). `fit_Dc_stress` fits the load-piston stress trace F_load/A, which is the network stress at the held face because the load piston is solvent-transparent, from `DC_STRESS_TMIN` (250 τ).

The reasons the hold value is not trusted:

- A hold is not one poroelastic transient. On the quarter gel (comp_3 level 0.30, network strain 0.15, fine file) the first ~6 τ₁ follow D_c = 0.109 (κ = D_c/M = 0.22, near the permeation Darcy κ of 0.24–0.32). A slower process then adds +0.012 of interior strain (20 % of the transient) and takes 8 % off the load-piston stress between ~1000 and ~4000 τ. A whole-hold fit with poroelastic modes only absorbs it by halving D_c to 0.054, and that fit is unstable in the number of modes (N = 8 gives 0.021) and in the trim.
- The hold D_c depends on what is assumed about the onset state. Fitting the run 7 hold with the mode amplitudes a constant-rate ramp predicts, A_m ∝ (1 − e^{−λ_m t_r})/m³, the decay rates give 0.051, the spectrum shape 0.065–0.085 and a fit with no free amplitudes 0.063–0.067.
- The interior is already strained at the hold onset. The D_c-independent mode-1 amplitude is half the prediction: the hold redistributes 0.058 of strain against a closure of 0.111 and a plateau profile strain of 0.082, so the interior was strained by about 0.02 during the ramp (Rg extent −4.5 σ during the ramp, load-piston spike ~0.3 against a drained 0.03, interior strain slightly negative over the first ~1500 τ of the hold). This is an undrained response of the ρ = 0.47 solvent, a storage term the 1-D incompressible model lacks.
- On run 7 (coarse file: 2120-τ blocks, the first ending 2775 τ = 0.6 τ₁ into the hold, so the fast phase sits inside one block) the profile fit gives 0.051 / 0.046 at levels 0.10 / 0.15 against the stress trace 0.073 ± 0.008 / 0.069 ± 0.013 (τ₁ 4500 / 4200 τ; the 0.15 tail bends to a slower component beyond ~5000 τ), that is κ 0.13 → 0.19 against Darcy 0.24–0.32. Levels above ε ≈ 0.2 are held for fewer than 5 τ₁.
- Other explanations were tested and excluded. Projecting out the per-snapshot rigid translation (the barycentric velocity v⁰ = φ_s v_s + φ_p v_p is uniform in z, so it drops out of the strain equation but not of the displacement one the fits use) changes D_c by less than 1 % (run 7 level 0.10: 0.0507 against 0.0506). The finite-difference convention test passes for the hold, permeation and unload fitters. Excluding the faces from the hold fit does not raise its D_c. Windows starting after the ramp-end skin has relaxed still read 0.05. A global strain-dependent D(ε) fit gives the hold half the re-swelling D at equal local strain.

The hold's consolidation and the free re-swelling are different processes on the same gel. The unload starts from a uniform state with free faces and shows none of the above; its profile fit gives 0.12.

#### `uniaxial_permeation_single.ipynb`

`Config(mode="permeation")`. `uni.load_permeation` builds one dict `P` holding the total, partial and network stress and density evolutions (cividis, final curve bold, the zero-flux `_ref` baseline dashed), both wet pistons, the flux and the permeability. The notebook has four sections:

1. Setup.
2. Results: figure 1 the run as a timeline (`fig_perm_stages`), figure 2 D_c (`fig_perm_Dc`), figure 3 M (`fig_perm_M`), figure 4 the reconciliation (`fig_perm_reconcile`).
3. Checks: the constant-dP flux (`fig_perm_flux`), the flux the fits imply (`fig_perm_flux_check`), M against the compression deck (`fig_perm_vs_compression`, using the M record the compression notebooks save with `save_M_record`), the unload (`fig_perm_Dc_unload`) and `fig_perm_v0_check`.
4. The state of the gel during the drive: pistons, stresses, volume fractions, pore sizes.

`fig_network_stress` and `fig_osmotic_pressure` are not called under permeation, because p_pore(z) is not a reservoir baseline under flow and the deck does not measure the solvent chemical potential. `fig_perm_permeability`, `fig_perm_frozen` and `fig_perm_density` are also not called; all of these functions remain in the library. Figure references inside `lib/uniaxial.py` use function names, not numbers.

The sync logs in when a trajectory is missing. An unresolved `<steps>` tag matches digits only, so a permeation pattern cannot stage a compression level's `_c<lvl>` dump.

##### Pistons, flux and the constant-dP κ

- Pistons: displacement, and the block-averaged `F_fluid/(lx ly)` from `piston_force_avg` measured against applied, with the reservoir virial pressures dotted.
- Flux: `Q_perm` is the slope of the permeate bead count over independent `q_win_steps` (150 000) windows, with its standard error and a slope drift test (`plateau_window_slopes`). The `z_perm` slope is drawn beside it. The block-bootstrap mean of the piston-velocity trace is only printed: the velocity blocks are dominated by the sheet's thermal jitter and their scatter is not an error bar.
- Permeability: `κ = Q_perm L/(A ΔP_ext)` with the measured ΔP_ext (primary) and the applied one, CIs in quadrature.

##### D_c and M from the polymer displacement

`uni.add_perm_displacement` reads `disp_z_polymer` (the deck's `fix ave/chunk` of the per-atom z displacement since the reset, every `disp_nfreq` = 50 k steps, block-averaged) and the bounding-box thickness trace.

**Model.** 1-D poroelasticity, `∂u_z/∂t = q(t) + D_c ∂²u_z/∂z²`, with `u_z = 0` at the support and a free feed face. The pistons prescribe the pressure drop and the total stress is uniform, so the network stress at the support equals ΔP from t = 0⁺. The strain then obeys the diffusion equation with Dirichlet ends, and the displacement modes are 1 − cos(kπζ) at k²π²D_c/L₀², τ₁ = L₀²/(π²D_c). The drive is the applied-ΔP history read from `piston_pressure`, so any ramp shape is treated exactly.

**D_c (`fig_perm_Dc`).** The profile snapshots are fitted with free modal amplitudes (window-averaged modes). Since 2026-11-24 the model is driven by the polymer-support pair force / A (`support_force`), which is the network stress at the pinned face, instead of the applied ΔP ramp. The two differ whenever the gel does not rest on the support at zero flux: in random_10b_wide_perm_1 the gel floated 8 σ above the plate, fell through the whole ramp carrying no load, landed at the ramp end and took 2.3 times the steady load from the pistons' momentum before settling 1.8 M steps later. With the applied ramp the fit gave R² 0.87; with the support force 0.99. The fit domain is the dense gel only (`PERM_DENSE_FRAC`): displacement binned by current z is biased where the density has a gradient, and L₀ ends at the half-density top edge of the first snapshot's count profile. The loader refits three variants (theory amplitudes A_k = 2e₀/(kπ)², and the first and second half of the snapshots after the load arrives); their range is the interval D_c carries into `perm_reconcile`. The bounding-box trace fit is off by default (`PERM_TRACE_FIT`): on a gel with a fringe the bounding box follows the fringe.

**Skin test (`perm_skin_test`, panel d).** A layer at the support that is less permeable than the bulk acts as a hydraulic resistance in series at the drainage face. For a thin layer this is the same as applying the support load a distance ℓ below the pinned plane, where ℓ is the skin's excess equivalent thickness (bulk permeability times the skin's resistance, minus its real thickness). The linear solution is refitted (scale and D_c) on a slab of length L₀ + ℓ for a scan of ℓ. A skin shows as a residual minimum at ℓ > 0, and D_c(ℓ = 0)/D_c(ℓ*) is the factor by which it lowers the whole-gel D_c. On synthetic profiles with ℓ = 8 σ and D_c = 0.2 the scan returns 7.9 σ and 0.199 (the no-skin fit of the same profiles gives 0.119). random_10b_wide_perm_1: the residual is smallest at ℓ = 0 and passes the 25 % tolerance (`PERM_SKIN_TOL`) at 5.5 σ, so a skin lowers D_c by less than 30 %.

**M (`fig_perm_M`).** M = ΔP_ext L_ss/(2|u_F|) with L_ss = L₀(1 + c), the steady thickness of the dense gel; u_F from the per-atom pairing over the dense reference bins (primary), from the steady-profile parabola, and from the bounding box (kept for reference; it counts the fringe). The CIs include the measured ΔP_ext in quadrature. The loader prints κ = D_c/M next to the Darcy κ.

**Flux check (`fig_perm_flux_check`).** The flux the fits imply, Q(t) = A[(D_c/M) ΔP(t)/L₀ − d⟨u_z⟩/dt], is drawn against the measured piston-velocity trace, the N_perm-slope Q and the cumulative bead count. Nothing in it is fitted to the flux.

perm_2 (dP = 0.1) gives D_c 0.18 (profile) / 0.21 ± 0.03 (trace), M 0.48, κ = D_c/M 0.4 against Darcy 0.29. The profile fit is under-determined there (six stress-cadence snapshots, the transient unrecorded), which is why the deck uses the 50 k cadence and a 1 M-step ramp.

Knobs:

| Knob | Default | Meaning |
|---|---|---|
| `PERM_GAP` | 1.0 σ | The pinned polymer face is z_support + `PERM_GAP` |
| `PERM_L0_MODE` | `'dense'` | L₀ of the fit: half-density top edge of the first snapshot's polymer count profile minus the pinned face; `'face'` is the reference polymer-stress edge (in the fringe on the random gel), `'bb'` the bounding-box thickness |
| `PERM_DENSE_FRAC` | 0.8 | Bins of the fit domain, of the contact layer and of the per-atom profile hold at least this fraction of the bulk polymer count; 0 restores the `Ncount_min` floor |
| `PERM_L0_STEPS` | 500 000 | Window of the bounding-box L₀ before the ramp and of the unload's L_ss before the ramp-down, ending at the last sample before the ramp (100 000 before 2026-11-11) |
| `PERM_FORCING` | `'support'` | Load history driving the model: the polymer-support force / A (falls back to `'applied'` without a `support_force` file); `'applied'` is the ΔP ramp, `'measured'` the reservoir pressure difference |
| `PERM_DC_AMPS` | `'free'` | Modal amplitudes: fitted per mode, or `'theory'` (A_k = 2e₀/(kπ)², two parameters) |
| `PERM_DC_VARIANTS` | `True` | Refit the variants that give the D_c range |
| `PERM_TRACE_FIT` | `False` | Fit the bounding-box thickness trace as a second D_c |
| `PERM_DC_WINDOW_AVG` | `True` | Average the fitted modes over each `ave/chunk` window |
| `PERM_DC_BOUNDS` | (1e-6, 50.0) | D_c search interval |
| `PERM_DISP_RESET` | `'auto'` | Where the deck reset the displacement: ramp start (decks since 2026-09-28) or ramp end (older decks), detected automatically |
| `PERM_TRACE_SKIP` | 0 | Steps after the ramp start left out of the trace fit, to skip a short-ramp run's piston ringing |

##### Velocity-controlled and frozen stages

Functions: `perm_vctrl_info`, `_perm_stage`, `add_perm_vctrl`, `add_perm_frozen`, `perm_reconcile`; paths from `Config.vpath` and `Config.fpath`. Knobs: `PERM_VCTRL` and `PERM_FROZEN` (both `True`: read the stage when its files exist), `PERM_VCTRL_SKIP` and `PERM_FROZEN_SKIP` (both 0 steps). `PERM_FROZEN_END` (0 = the whole stage) is the last step of the frozen stage that is read: the stage's windows are trailing ones, so a stage whose end is not flow through the gel any more has to be cut there. With the polymer stopped, the velocity-controlled feed piston can reach the gel's outermost beads, which cannot give way (job 55987856: contact from step 12.78 M, `P_feed` 10¹⁴ by 13.3 M); `add_perm_frozen` prints a warning when |ΔP| exceeds `P_BARO` in the samples it reads.

**Velocity-controlled stage (`_v` files, deck Phase 2v).** After the drive the feed piston moves at a prescribed velocity and ΔP is measured. `add_perm_vctrl` fills `P['vctrl']`:

- q = |`vfeed_set`| from `vctrl_info`.
- ΔP = the block-bootstrap mean of `P_feed_meas − P_perm_meas` over the longest trailing window (up to the whole stage) in which neither ΔP nor the N_permeate slopes drift.
- L = the stage's own Rg thickness.
- κ_v = q L/ΔP (key `v_measured`, the primary).
- `Q_ratio` = Q_N/(q A), the N-slope flux over the same window, is a check: it is 1 up to the feed/permeate density ratio when the gel passes what the piston pushes.

**Frozen stage (`_f` files, deck Phase 2d, `FROZEN=1`).** The polymer is stopped and the solvent is still driven. When `vctrl_info` exists, `add_perm_frozen` treats the `_f` files as it does the `_v` files and gives κ_f; its `ratio` is taken against the velocity-controlled live stage when there is one, otherwise against the drive.

Without a velocity-controlled stage it uses the constant-dP form (files `_<NSTEPS>_f`, optional in the sync). Q is the N_permeate slope over independent `q_win_steps` windows of the frozen interval, with the same drift-tested plateau as the live drive; dP is the measured piston dP over that window; L is the live steady thickness, since the network no longer moves. κ_f = Q L/(A dP) with its CI, the z_perm-slope variant and the ratio κ_frozen/κ_live per estimator are reported, and `print_perm_summary` adds the FROZEN line. Flow through a rigid network and weak flow through the swollen one should agree to leading order, so the ratio measures the network's own response to the drag during the live drive.

**Reconciliation.** `perm_reconcile` collects D_c (profile, trace), M (primary), κ_v, κ_f, the drive's constant-dP κ and D_c/M into `P['recon']` with the closure ratios (D_c/M)/κ and κ_f/κ_v, and prints the table. Runs without `_v` files reconcile on the constant-dP κ (perm_3: D_c/M = 0.335 against 0.231, ratio 1.45 [1.23, 1.67]).

##### Reading the Darcy κ of the live drive

The Darcy κ is not a reference-state property. On perm_3 the consolidation-fit D_c is 0.147 (profile) / 0.159 (trace) while κ_Darcy × M_perm is 0.231 × 0.416 = 0.096, a factor of about 1.5 apart on one run; the compression unload gives 0.120.

- Under ΔP = 0.1 the support side compacts to φ_s ≈ 0.33 (feed side 0.44, zero flux 0.447). Q L/(A ΔP) is the harmonic mean of κ(φ_s) over that profile, so it reads below the κ of the undeformed gel.
- The small-ΔP quarter-gel runs show this directly: κ 0.30–0.32 at ΔP 0.02–0.05 against 0.24 at 0.1, so κ(ΔP → 0) × M ≈ 0.13, at the transient D_c and the unload D_c.
- The consolidation transient starts from the undeformed state and weighs the early, uncompressed modes, which is why its D_c sits above κ_Darcy M.
- `FROZEN_AT=ref` runs the drive on the network frozen in the reference state (no compaction, no network response). Its Darcy κ is the undeformed gel's and the one to hold against the compression side. `FROZEN_AT=end` (the end-of-drive frozen stage, the batch default) isolates the network's dynamic response at fixed compaction.
- L, A, ΔP and ρ_s,0 in the Darcy formula were checked on perm_3 (ΔP measured 0.0998, the N_permeate-slope and z_perm-slope fluxes agree, ρ_s,0 is the bulk reservoir density). The formula is not where the factor comes from.

##### Permeation unload

Functions: `add_perm_unload`, `load_perm_unload`, `perm_forcing_unload`, `fig_perm_Dc_unload`. Knobs: `PERM_UNLOAD` (`True`), `PERM_UNLOAD_TRACE_SKIP` (0), `PERM_UNLOAD_L` (`'ss'`).

Runs of `uniaxial_permeation.lmp` with `UNLOAD=1` ramp dP back to zero after the drive and record a free re-swelling hold on files tagged `_<NSTEPS>_u` (`Config.upath(name, '')` in permeation mode; optional in the sync).

- `load_perm_unload` builds the same `disp` dict as the drive's: `disp_z_polymer_…_u`, reset at the ramp-down start; the forcing from the applied-dP decrements of `piston_pressure_…_u`, with t₀ the first sample minus one sampling interval; and a bounding-box trace made of the drive's last `PERM_L0_STEPS` samples (the compressed steady thickness L_ss) followed by the unload's.
- `fit_perm_Dc` runs on it unchanged except for the L (`PERM_UNLOAD_L = 'ss'`: L₀(1 + c) of the drive fit, the state the re-swelling starts from; `'ref'` is the drive's L₀) and the "steady" window (the last `plateau_frac` of the re-swelling).
- The modes are the drive's. When dP returns to zero the network stress vanishes at both faces at the same time (support: the dP it carried; feed face: free), so the strain is Dirichlet at both ends and τ₁ = L²/(π²D_c). This is the loading transient's boundary-value problem, not the compression unload's, whose faces are freed from a uniform strain.
- The rigid shift u₀ (contact layer) is the network lifting off the support.
- `perm_M_estimates` on the re-swollen parabola, thickness and trace asymptote gives M(unload) = dP L/(2|u_F|). There is no Lagrangian variant because the unload writes no trajectory.
- κ(unload) = D_c(unload)/M is formed with the drive's primary M and with the unload's own. `print_perm_summary` adds the UNLOAD line.

The stress / Terzaghi companions of the unload are not built, although the files are written (σ_zz/xx/yy of polymer and solvent, solvent density, reservoir pressures).

##### Volume fractions and pore sizes

`uni.add_perm_volume_fractions` computes the three solvent volume-fraction estimators of the compression notebooks: the mass fraction from every density snapshot, and Voronoi and λ-calibrated Voronoi on `VOR_EVO_FRAMES` (5) frames over the drive plus `VOR_MAX_FRAMES` (4) inside the steady window.

The pore pressure falls from feed to permeate across the membrane under flow, so the calibration pressure handed to λ(φ_p, P) is chosen per bin by `Config.P_CAL_MODE`:

- `'pore'`: the ramp between the measured feed and permeate reservoir baselines across the membrane. The permeation notebook uses this.
- `'const'`: `P_CAL` (1.5) everywhere. This is the `Config` default and is correct for a drained equilibrium (compression).
- `'thermo'`: the local P_th = −⅓ tr σᵗ.

The same knob applies to `add_volume_fractions` in the compression notebooks.

`fig_perm_volfrac` shows reference against steady state, the P_local(z) used and the resulting λ(z). `fig_perm_volfrac_evolution` shows the mass fraction and the calibrated Voronoi fraction. Its box gives the steady φ_s drop across the gel (feed face − permeate face) from `perm_dphi`, a quadratic fit over the interior bins, next to the in-gel mean change from the reference. Panel (c) turns that drop into `D_c = q_s L / ln[(1 − φ_P)/(1 − φ_F)]`, the integral of `q_s = D_c/(1 − φ_s) dφ_s/dz` with `q_s = Q_perm/A`, for the whole gel and per half.

The same frames feed `lib/psd.py` (`uni.add_perm_psd`, `fig_perm_psd`): the geometric, probe-accessible porosity per bin (`PSD_R_PROBE` = 0.5 σ) and the pore-size distribution by the Gelb–Gubbins largest-included-sphere covering (cKDTree distance field on a grid, one Euclidean distance transform per radius level), reference against steady state, with `ξ = √κ` printed next to the geometric mean pore diameter.

The remaining state figures are `fig_perm_pistons`, `fig_total_stress` (normalised by P_perm), `fig_partial_stress` and `fig_thermo_pressure`.

#### `isotropic_compression_single.ipynb` and `isotropic_compression_sweep.ipynb`

These are the bulk-modulus counterparts of the uniaxial compression notebooks. `lib/isotropic.py` takes its readers, statistics, Terzaghi split and profile figures from `lib/uniaxial.py`. `iso.Config` points at `flow_data_local/raw_data/isotropic_compression/<RUN_ID>` and Expanse `lammps_runs/isotropic_compression`. Profiles are read in a core column through the gel along z, x and y (bath | gel | bath).

**Per level.**

- The measured plate ε_vol and the network's own Rg and bounding-box volumetric strains.
- **K_network**: the plateau-averaged interior ⟨Π⟩ increment from the ε_vol = 0 reference, over ε_vol, per profile axis and their mean (formed exactly as `M_network`).
- **K_plates**: the mean plate pressure over the drift-tested block-bootstrap plateau window, as an increment, over ε_vol (formed exactly as `M_piston`).
- `K_net(Rg)`: the same stress over the network's own strain.
- δ_sys from the plate tail (`RELAX_SYS`).
- The isotropy of Δσ′. There is no G under isotropic loading.
- The bath check. There is no solvent-expelled reading: the six plates are solvent-transparent, so the solvent squeezed out of the gel stays in the same bath and the wet pistons do not move.

**No diffusivity from the hold.** The hold-onset `displace/atom` files are a one-snapshot reference, and the `bb`-seated hold mixes the shell folding in with the drainage.

**D_m from the unload.** With `UNLOAD=1` (deck Phase 2c, set in both batches), after the hold the six plates open to the seated gaps + 2 `UNLOAD_MARGIN` and the cube re-swells freely, writing the level's full file set as `_<hold>_u<level>`. With all six faces unloaded this is the three-dimensional unconstrained problem, whose diffusivity is **D_m = (K/M) D_c = K/ζ**, not the confined D_c = M/ζ of the uniaxial holds.

- `load_unload` / `fit_Dm_unload` fit the free-cube modes to the cumulative displacement profiles referenced to the re-swollen plateau: along axis a the lag relaxes in odd cos(kπζ) modes at exp(−k²π² D_m t/L_a²) times the transverse slab factors `Σ_odd 8/(kπ)² exp(−k²π² D_m t/L_b²)`, with a free constant per block and a time-independent amplitude per mode. τ₁ = L²/(3π² D_m). The lengths are the compressed Rg edges at the hold end (`UNLOAD_L = 'hold'`), a plateau, where the Rg edge is exact.
- One estimator since 2026-11-21: the per-axis fits and the joint fit (the headline). The bracket next to D_m is the range of the per-axis fits, not a confidence interval. A level is flagged as not a measurement when the per-axis fits disagree by more than `DM_AXIS_SPREAD` (2), a fit sits at a `DM_BOUNDS` edge, or the profile blocks are longer than τ₁. `fig_Dm` draws the profiles and the numbers, `fig_Dm_time` the same fit against time (first-mode amplitude per block). Decks since 2026-11-21 write fine-cadence copies of the unload profiles (`disp_*_polymer_cum_fine`, `DISP_FINE_NFREQ` = 5000 steps in both batches, from the start of the retraction), which the fit prefers; the standard files have one block per unload/20, about one τ₁ on a 12-cell cube.
- The Rg-edge, Rg-volume and interior-Π estimators were removed. An Rg edge weights the normal strain by (3/2)(1 − 4X²/L²) (mid-plane 1.5, faces 0), so during a transient it lags the true edge and does not follow `1 − R(t)`; fitted with `R(t)` it reads D_m low. The interior Π relaxation was a few block means with a free offset and amplitude.
- κ = D_m/K from the same run. With `M_REF` set, the closure G = ¾(M − K) and the confined D_c = (M/K) D_m are also given (`fig_Dm`, `fig_Dm_sweep`, `fig_closure`).

**Other figures.** Network anisotropy: `fig_ratio` (σ′_ii/Π over the hold, per axis) and `fig_anisotropy_run` (the whole-run `stress_aniso` file, optionally beside the `slab_with_support` run's file). Solvent volume fractions in the core columns along x, y, z: mass fraction, Voronoi and λ-calibrated Voronoi (`add_volume_fractions`, one tessellation per frame; `fig_volfrac`, `fig_volfrac_evolution`).

**Knobs beyond `uni.Config`.** `AXES` (`('z', 'x', 'y')`), `K_STRAIN` (`'disp'`), `PLATE_AREA` (`'gap'`), `K_SUBTRACT_REF` (`True`), `EPS_MIN` (1e-3), `SS_FIT_MODE` (`'auto'`), `SS_FIT_TOL` (0.05), `SS_FIT_NMIN` (3), `M_REF` and `G_REF` (`None`), `HOST` (`'expanse'`), `DM_UNLOAD` (`True`), `DM_N_MODES` (3), `DM_TRACE_SKIP` (0 τ), `DM_BOUNDS` ((1e-4, 10.0)), `DM_AXIS_SPREAD` (2), `core_margin` (6.0 σ).

#### `simple_shear_single.ipynb` and `simple_shear_sweep.ipynb`

These are the shear counterparts of the uniaxial compression notebooks; the code is in `lib/shear.py`. `sh.Config` points at `flow_data_local/raw_data/simple_shear/<RUN_ID>`. `STRAINS` lists the levels, whose files are tagged `_g<strain>`.

##### Decks

`DECK` in the Config cell selects the run type. The same `_g<strain>` files, columns and figures apply to all three.

| `DECK` | Run type | Expanse folder |
|---|---|---|
| `"simple_shear_iso"` | The isolated gel in the two-piston bath, grip plates seated on the x-faces and driven ∓z. This is the current workflow. | `lammps_runs/simple_shear` |
| `"simple_shear_plates"` (`Config` default) | The archived plate-driven deck | `lammps_runs/simple_shear_plates` |
| `"simple_shear"` | The archived laterally periodic two-piston slab (deck file `simple_shear_periodic.lmp`) | `lammps_runs/simple_shear` |

Runs of the two archived decks still load.

**Isolated deck.** The gap axis is the box x and the shear direction is z, the transpose of every other shear deck. `lib/shear.py` maps the file names (`stress_profile_x_*` into the usual `stress_profile_z` slot, `disp_z_polymer*` into `disp_x_polymer*`) and relabels the stress axes x ↔ z and the box lx ↔ lz at read time (a module flag set by `Config`). Every number and figure is therefore in the notebooks' usual frame:

- The plotted "z" is the deck's x, "u_x" is the deck's u_z, σ_xz is itself, and N1 = σ_shear − σ_gap.
- The dash-dot planes are the grips (`plane_label` "grip"), at the plates' seat planes (`plate_pressure_*` columns 2–3).
- The solvent-only "film" bins are the bath beyond the gripped faces (p_pore = the bath).
- `S_top` / `S_bot` are the right / left grip tractions.

Two things are not mapped and must be read directly: the deck's `stress_profile_z_*` end-zone diagnostic along the shear direction, and the `N_r` / `N_l` (net) and `Np_r` / `Np_l` (polymer-only, the seated preload) grip normal tractions in `plate_pressure_*`.

**Periodic two-piston deck.** The dash-dot planes are the gel's Rg faces. The bulk starts 5 σ inside them (`plate_excl` 5, `wall_margin` 6: driven layer plus tracking slab; the `Config` defaults of 3 and 4 are for the plate deck). The solvent-only "film" used for p_pore is the two reservoirs (full bins only), `G_plate` is the face-layer force per area, and the unload D_c fit maps the bulk bins onto the driven-layer COM planes (`H_drv` apart). The residual figure (`fig_residual`) skips itself for a run without an unload, as it does for `simple_shear_plates` runs without `-var unload 1`.

##### Quantities per level

- The held strain γ.
- **G_network**: the plateau-averaged interior ⟨σ_p,xz⟩ increment from the γ = 0 reference, over γ (t-interval over bins ⊕ reference CI, formed exactly as `M_network`).
- **G_series**: the fine `stress_series` ⟨σ_p,xz⟩(t) over the drift-tested block-bootstrap plateau window, as an increment, over γ (formed exactly as `M_piston`).
- **G_total** from σᵗ_xz, as the poroelastic check.
- N1 and N2.
- The bulk P_th = −⅓ tr σᵗ (the P* = 1.5 check).
- D_c from the level's unload and **κ = D_c/G**, the same hydraulic permeability as D_c/M in compression.

The notebooks need the `_ref` files and `stress_series` files that decks since 2026-09-23 write; older runs are reported as absolute values and flagged.

**Figures.** Single: strain, the displacement profile (1b), total stress, the reference normal stresses (2b), partial σ_xz, network σ_p,xz in the reference and final states, stress history, G comparison, N1/N2, the D_c fit, κ and P_th. Sweep: the overlays, plus G against γ, the stress–strain curve with least-squares and through-origin slopes, and N1/N2, D_c, κ and P_th against level.

**Figure 2b (`fig_ref_normal_stress`).** The three total normal stresses before shearing: bulk σᵗ_xx, σᵗ_yy, σᵗ_zz and P_th over the γ = 0 reference window with the lateral pre-stress ⟨σᵗ_xx − σᵗ_zz⟩ and ⟨σᵗ_yy − σᵗ_zz⟩ (block bootstrap), and their window-mean z-profiles. The reservoirs must read P_bath in all three. A lateral component away from σᵗ_zz inside the gel would be a pre-stress of the slab: the network would be transversely isotropic, so the compression (C₃₃ − C₁₃)/2 and the shear C₄₄ would differ. On shear_3 (two-piston) and film7 (plates) both differences are 0 within ±0.001–0.0016.

##### Strain

**Prescribed strain (plate deck).** γ = (x_plate,top − x_plate,bot − x₀)/plate_sep from the plate COMs, both prescribed by `fix move`. A bounding-box-gap denominator with an outermost-1 %-bead numerator overstates the sheared thickness by about 5 %, because `add_plates_to_gel_periodic` puts each plate 1–3 σ inside the outermost beads. For the same reason the bulk region is cut `plate_excl` inside the plate planes.

**Strain files.** The `shear_strain_*` files have four columns (step, gel_lz_initial, plate_sep, γ). The measured surface-COM strain goes to `shear_strain_surface_*_g<γ>.dat` as a slip diagnostic. The per-level `shear_strain_*_g<γ>.dat` spans the drive and the hold; the notebook and `plot_shear_strain_sweep.py` take its last row as the level's held γ. The hold is padded to the next stress-averaging epoch. Stem-suffixed whole-run trackers (`box_dimensions`, `gel_dimensions_rg`, `gel_volume_{bb,rg}`, `shear_strain`) feed `plot_lammps_log.py`'s convergence and shear-diagnostics figures.

**The strain G divides by (`G_STRAIN`, the counterpart of `M_STRAIN`).** G and κ divide by the strain of the network, not the prescribed boundary strain.

- `'disp'` (default): the interior slope of the steady cumulative displacement profile u_x(z) from `disp_x_polymer_cum_*`. This is the plateau mean minus the γ = 0 reference-window profile, with a line fit over the populated bins between the planes and `DISP_TRIM_BINS` = 3 dropped per face. The CI is the larger of the least-squares and frame-to-frame intervals. z does not change under shear, so the current-z bins need no conversion.
- `'rg'`: (Rg_xz − ⟨Rg_xz⟩_ref)/Rg_zz² from `shear_strain_rg_*`.
- `'trk'`: the tracking-slab strain (two-piston deck: hold plateau minus its γ = 0 reading).
- `'applied'`: the prescribed plate / driven-layer strain.

When a source is missing the chain falls back in that order. The chosen strain's CI is added in quadrature to every G. Runs of decks older than 2026-10-03 have neither file and fall back to the prescribed strain with a printed FALLBACK note. Figure 1b (`sh.fig_disp_profile`) shows the profiles, line fit and residuals with every strain measure listed, and `print_summary` lists them per level together with the G the prescribed strain would give. `L['gamma']` is the strain used and `L['gamma_applied']` the prescribed one.

##### Shear D_c

**From the unload.** D_c is read from each level's unload (deck Phase 2c; `fit_Dc_unload`, `fig_Dc_unload`, `fig_Dc_sweep`; knobs `DC_UNLOAD` = `True`, `DC_UNLOAD_N_MODES` = 3, `DC_UNLOAD_BOUNDS` = (1e-4, 100.0), `DC_UNLOAD_SKIP` = 0 τ). The fit is to the lag of the cumulative displacement profile to its late plateau after the grips return to γ = 0, in sine modes with one D_c. The one-time Kelvin–Voigt alternative is fitted beside it and the residuals' verdict is printed.

**Why not from the hold.**

- The hold-referenced `disp_z_polymer` is pinned to the single `displace/atom` snapshot at the hold onset, so every block carries the same static 0.04–0.09 σ pattern (block 1 against block 10 correlation 0.96–1.00 in sweep 55082229). A model A_k(1 − e^{−t/τ_k}) can absorb a constant only by sending D_c to the top of its bounds.
- Nothing is left to relax by the first block. The 5000-step series (tracking-slab strain, σ_xz) are at their plateau in the first 25 τ of every hold, and the first displacement block differs from the late mean by at most 0.02 σ, against the 0.2–0.6 σ ramp lag a Tanaka–Fillmore D_c of 0.03–0.1 would leave.
- In the two-fluid picture a transverse network mode whose wavelength exceeds the screening length (≈ √κ_D ≈ 0.5 σ) carries the solvent with it and relaxes at the Kelvin–Voigt rate G/η_s ≈ 0.1 τ⁻¹, the same for every mode, not diffusively. Pure shear changes no volume and generates no pore pressure, so κ does not enter.
- The D_c figures of the periodic plate runs shear_1 / shear_2 show the same static pattern at ±2.5–2.8 σ and are an artefact. `simple_shear_plates.lmp` zeroes the angular momentum of polymer + solvent every 100 steps (`fix momentum … angular`) during the drive and the hold, which removes most of the interior shear rate each time, leaves the interior nearly unstrained at the halt (2.7 of 3.0 σ missing at ẑ = 0.25) and lets it catch up before the first 900k-step block. The two-piston decks zero linear momentum only, and shear_3 (same gel, same ramp) shows no lag.

**Fine-cadence hold file (`fit_Dc_fine`, `fig_Dc_fine[_sweep]`).** For runs that write the fine file, `fit_Dc_fine` forms the lag still to relax, [u_late(ẑ) − u(ẑ, t)]/L, where u_late is the plateau of the coarse `_cum` file so that the onset snapshot cancels. It projects the lag per snapshot on sin(2πkẑ), k = 1..`DC_FINE_N_MODES` (3), fits each amplitude with A_k e^{−t/τ_k} + c_k over the hold, and reports τ_k, the diffusive reading D_c,k = L²/(4π²k²τ_k) (equal for every k if diffusive) and the one-time reading 1/τ₁ = G/η_s. A joint one-D_c fit and a joint one-τ fit are compared by residual (`pref`). On synthetic data both readings recover their inputs to 1 % and the preference is unambiguous either way.

##### Comparing with the compression G

Use the same bath pressure and the G = (σ′_zz − σ′_xx)/2ε values from `uniaxial_compression_sweep.ipynb`. The converter's calibrated solvent deletion around the plates puts the gel at P* = 1.5, and the notebook's bulk `<p_total>` is the check. The thermo `press` of this box is not the gel pressure.

#### `manuscript_figures.ipynb`

Named `cross_deck_comparison.ipynb` until 2026-11-22. It opens with the figures chosen for the paper: (1) strain and solvent exchange over the compression cycle (`cmp.fig_strain_solvent`, twin axis), (2) displacement evolution in re-swelling and, once a live-network permeation run is loaded, in permeation (`ms.fig_displacement_evolution`), (3) the four networks with their frozen-pore fraction (`cl.fig_network_compare` with `cl.load_void_persistence`).

`lib/manuscript.py` finishes a figure after its library has drawn it. `ms.capture()` holds back the library's `plt.show()`; `ms.finish(fig, name, text, size, fonts, limits, legends, style)` applies the cell's edits and the global settings of `ms.set_style` (font, sizes per kind of text, line-width and marker scales, grid, title switches, dpi, formats, output folder) and saves to `flow_data_local/plots/manuscript/`. Every piece of text has a key (`suptitle`, `ax0.xlabel`, `ax0.legend.2`, `ax1.text.0`, `fig.text.0`); `ms.print_text(fig)` lists them. In a cell's `TEXT` dict a string replaces the text, `''` removes it and `None` keeps the library's wording, which is the setting for text that carries a fitted number.

The frozen-pore fraction Φ is the void persistence of `lib/screen_voids.py`: a grid point is a pore in a frame when no polymer bead centre lies within the probe radius, and Φ is the share of the pore pattern that stays in place over the window (0 renewed, 1 permanent). `cl.load_void_persistence` runs it on each gel's compression reference trajectory (16 frames, 750 τ) and caches the result under `raw_data/chain_length/void_persistence/`. At a probe radius of 1 σ: rot45 5-bead 0.05, rot45 15-bead 0.04, random 5-bead 0.31, random 10-bead 0.07.

The cross-deck cells that follow are skipped while no permeation deck is loaded. The frozen-reference run `random_10b_wide_perm_fref_1` is left out: a frozen network's stress is a single instant's value, so its profiles are about 50 times rougher from bin to bin than a live run's (0.077 against 0.0014 in total stress).

This notebook builds figures from more than one deck's run. One Config per deck, copied from that deck's own notebook, goes into `CFG`. `cmp.load_deck` loads each run with its own library (`lib/uniaxial.py` for the two-piston compression and permeation decks, `lib/shear.py` for `simple_shear_plates`) into `DECKS[name] = dict(kind, cfg, R, levels | P)`. The comparison functions in `lib/compare.py` take the decks they need. Plots go to `flow_data_local/plots/comparison/<NAME>`.

##### D_c(φ_s) and κ(φ_s) (`cmp.phi_dependence`)

- The compression levels are uniform states. Each has σ′ and the displacement-strain secant M = Δσ′/ε_disp; together they give the equilibrium curve σ′(φ_s) (PCHIP, Monte-Carlo CI) for the closures.
- The compression-side D_c is the unload profile fit of every level that has one (`cmp.unload_states`, over all compression decks in `CFG`; at present the `_andrelax` run, level 0.10), with κ = D_c/M from that level's secant M. Closure (ii) and `print_phi_kinds` use it too. The hold D_c is a load measurement and is not used; `compression_states` carries no D_c.
- The steady permeation membrane spans a range of φ_s at one flux. Locally D_c = q φ_p/(dφ_s/dz), from the flux and profile slope only, and κ = D_c/M_perm with the run's own M_perm = ΔP L/(2|u_F|).
- No modulus is derived from the volume fraction, because the φ_s profiles are not trusted to that level; every κ uses a displacement-strain M.
- `PHI_KIND` picks the volume fraction: `'cal'` (λ-calibrated Voronoi, default), `'vor'` or `'mf'`.

Run 7 + perm_3 + the unload of the `_andrelax` run: the two decks agree on σ′(φ_s) to about 6 %. The unload D_c (0.12 at φ_s 0.40) sits at the upper end of the steady permeation's local D_c (≈ 0.09, CI 0.06–0.15) and above κ_Darcy × M_perm (0.096).

##### Opening figures (permeation deck only)

`cmp.fig_thermo_pressure_norm` draws the permeation notebook's normalised thermodynamic-pressure figure with its two `ΔP_th = (4/3)(G/M) ΔP_ext` guides (`G_REF` = 0.2 from the shear decks, `G_COMP_REF` = 0.05 from the compression holds, both set in the notebook's deck Config) and copies the PNG into the comparison folder.

`cmp.fig_solvent_measures` recasts the permeation notebook's volume-fraction figure (7a) as concentration:

- Left axis: `c_s = φ_s^cal ρ_s,0/M_s` in beads/σ³ (`M_s` = 1).
- Right axis, scaled by `ρ_s,0`: the mass (= number) fraction `ρ_s/ρ_s,0` and the probe porosity `ε_g`.
- Reference dashed, steady solid, with 95 % bands; the individual steady Voronoi frames of `c_s` are drawn faint (`CONC_FRAMES`).
- The reservoirs are trimmed to `CONC_MARGIN` (5) σ; on the feed side they run up to the piston.
- Panel (b) shows the two other measures minus `φ_s^cal` (`CONC_DIFF_PANEL`).

The deck Config runs the PSD pass with `PSD_R_PROBE` = 0.45 (the permeation notebook uses 0.5), and `load_deck` calls `uni.add_perm_psd` when `PSD_ENABLE` is set. The reason for 0.45: `psd.py` measures the grid distance from the bead surface, `d = r_nn − σ_bead/2`, and calls a point void when `d ≥ r_probe`, so the probe's centre must clear a polymer bead centre by 0.5 + `r_probe` σ. With 0.5 that is 1.0 σ, the LJ zero-crossing. A thermal solvent bead at kT = ε reaches `U_LJ(r) = kT` at r = 0.97 σ, so 0.45 (clearance 0.95 σ) is the solvent-bead-sized probe and 0.5 slightly over-excludes. The choice shifts `ε_g` by a few hundredths, not its gap of about 0.16 below `φ_s^cal`. On perm_3 the interior `φ_s^cal` is 0.39; the mass fraction is 0.42 (+0.03: solvent packed on the polymer that the calibration assigns to the network); `ε_g` is 0.23 (−0.16: the excluded shell); and the feed-reservoir mass fraction is 1.03 (denser solvent under the feed pressure).

#### Tests

The tests are in `scripts/tests/`.

- `run_plot_tests.sh` builds synthetic one-piston, two-piston and six-plate run trees (`make_fixtures.py`) and runs the three plotters on them.
- `run_notebook_tests.sh` executes the two-piston notebooks headlessly on the same tree.
- `lint_lmp.py` statically checks a deck: definitions before use, unfix/undump pairing, labels, forwarded `-var`s.
- `lib_headless_test.py` asserts the library's stage readings on the `two_piston_perm` fixture, which carries both the velocity-controlled and the frozen stage.
- `dc_convention_test.py` checks the hold, permeation and unload D_c fitters against a finite-difference solver.
- The six-plate fixture is built from a known K and a known D_m (an unload with the free-cube series), and the test asserts that the library returns both.

#### Archived notebooks

These are described in `archive/README.md`. Existing runs stay analysable.

- `archive/scripts/uniaxial_compression_one_pist_single.ipynb`, `…_one_pist_sweep.ipynb` and `uniaxial_permeation_one_pist.ipynb` analyse runs of the one-piston decks. Open them from `scripts/`, because their paths assume that folder.
- `archive/scripts/uniaxial_compression_one_pist.ipynb` is the long-form diagnostic notebook. It is still the place for the solvent-phase stress W_s,zz/V_solv diagnostics, the ss/pp pair-virial reconstruction and cross-virial check, the piston–gel contact analysis and the Widom-insertion appendix.
- `archive/volume_of_mixing.ipynb` computes ΔV_mix(P*) = V_mixed − V_pure_solvent − V_pure_polymer across the pressure sweep (P* = 1.0–2.0) from `box_dimensions_*.dat` files synced from Expanse with `paramiko`. The calibration-sweep analysis supersedes it.
- `archive/isotropic_compression.ipynb` is the predecessor of the two isotropic notebooks.

---

### 7c. Running Python scripts manually on a cluster

Post-processing can be rerun after a job without relaunching LAMMPS, for instance after changing an analysis script. Every command below runs from inside the run's working directory on the cluster:

```bash
cd ~/Documents/lammps_runs/<run_dir>
# e.g. cd ~/Documents/lammps_runs/slab_with_support_slab_support_5beads_tall_rho04_1.0_1.0_20260705_124556
```

Set these variables once per shell session; the commands below derive everything from them:

```bash
DATANAME="slab_support_5beads_tall_rho04"
INTERACTION="1.0_1.0"
TOTSTEPS=3000000
EPSSS="${INTERACTION%%_*}"   # first part:  e.g. 1.0
EPSSP="${INTERACTION##*_}"   # second part: e.g. 1.0
SCRIPTS=~/Documents/lammps_work/scripts      # the plotters live in $SCRIPTS/lib/
```

#### Load Python

Expanse (primary):

```bash
module load anaconda3/2021.05/q4munrg
```

Pod (backup):

```bash
module load miniconda
conda activate lammps_analysis
```

Bridges-2 is no longer used.

#### Plotters

`plot_lammps_log.py` plots T, P and volume convergence. Output goes to `./output_plots/`. No `--p-ext` flag is needed for any current folder.

```bash
python "$SCRIPTS/lib/plot_lammps_log.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}"
```

`plot_stress_profiles.py` plots the partial stress and volume fraction profiles.

```bash
python "$SCRIPTS/lib/plot_stress_profiles.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0
```

`plot_piston_data.py` plots piston position, velocity and force (`uniaxial_compression`, `uniaxial_permeation` and similar folders).

```bash
python "$SCRIPTS/lib/plot_piston_data.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0
```

#### `cavity_widom.py`: excess chemical potential μ_ex(z) (archived)

The script is in `lammps_work/archive/`, and the `widom_traj` dumps that feed it are off by default (`SKIP_WIDOM=1`); `archive/README.md` explains how to re-enable them. It was written for the removed `slab_with_flow` deck and was never ported to the uniaxial decks, which have no `postprocess.sh` hook for it. It is run by hand.

For `slab_with_support` (no piston forcing: `--piston-eps 1`, no exclusion buffer):

```bash
python ~/Documents/lammps_work/archive/cavity_widom.py \
  --traj "traj_files/widom_${DATANAME}_${INTERACTION}_${TOTSTEPS}.lammpstrj" \
  --out-dir "output_files/chemical_potential" \
  --out-stem "${DATANAME}_${INTERACTION}_${TOTSTEPS}" \
  --eps-sp "$EPSSP" --eps-ss "$EPSSS" \
  --n-bins 40 --n-trial 50000 --r-cavity 0.5 \
  --piston-eps 1.0 \
  --p-ext 1.5 --temperature 1.0
```

The legacy `slab_with_flow` compression mode used the same command with `--exclusion-buffer 2.0 --piston-eps 0.0 --p-ext 1.8`.

The trajectory file is in `traj_files/`, a symlink to scratch. If scratch has been purged, the trajectory must be downloaded again from the cluster or the simulation rerun.
