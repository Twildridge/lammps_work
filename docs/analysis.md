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
python ~/Documents/lammps_work/scripts/lib/plot_lammps_log.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000
```

For `shear_slab` — pass only the base `DATANAME` as the title, and the full run_id via `--run-id` for file lookups:
```bash
python ~/Documents/lammps_work/scripts/lib/plot_lammps_log.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000 \
    --run-id isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000
```

Or equivalently (same result as the auto-run call from `run_lammps.sh`):
```bash
python ~/Documents/lammps_work/scripts/lib/plot_lammps_log.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000
```

Output plots are saved to `./output_plots/convergence_plots/`.

---

**`plot_stress_profiles.py`** — partial stress and volume fraction profiles (`slab_with_support`, `triaxial_compression` single-level runs, `triaxial_permeation`, etc.)
```bash
python ~/Documents/lammps_work/scripts/lib/plot_stress_profiles.py \
    . \
    isolated_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_1.0_1.0_500000 \
    0
# third argument is OLDSTEPS — always 0 now (the old batch-file continuation path was removed 2026-08-06)
```

---

**`plot_piston_data.py`** — piston position and velocity (`triaxial_compression` single-level runs, `triaxial_permeation`, etc. — compression *sweeps* and `shear_slab` use their own consolidated plotters instead, see [§5b, simulation types](running_simulations.md#5b-simulation-types)). Multi-piston aware since 2026-09-16: the two-piston files (`# step z_feed z_perm` headers) get one line per piston plus a measured-vs-applied pressure panel and a `Q_perm` / solvent-expelled panel; one-piston files plot exactly as before. `plot_stress_profiles.py` overlays every `stress_z_piston_<sheet>` file in the Piston column and `plot_lammps_log.py` adds the bath-check and `Q_perm` panels to the flow diagnostics.
```bash
python ~/Documents/lammps_work/scripts/lib/plot_piston_data.py \
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
**D_c fit corrected 2026-09-28** (`fit_Dc`, both compression notebooks): a slab held between two drained plates at
the same bath pressure relaxes in the periodic strain modes cos/sin(2mπζ) — `u` in sin(2mπζ) and cos(2mπζ) − 1 at
4m²π²D_c/L², **τ₁ = L²/(4π²D_c)** — not in the (2ζ − 1) + cos(kπζ) shapes at k²π²D_c/L² used before (Terzaghi's
constant-load series bent to vanish at both plates; on the exact relaxation it returned ~5× the true D_c). The decay
times never changed, the conversion to D_c did: comp_1 at ε = 0.1 now reads D_c = 0.043 and κ = D_c/M = 0.13 (was
0.49 / 1.5). The hold-adequacy line uses the corrected τ₁ for the fitted D_c and keeps the deck's own sizing formula
for the *slow* line; levels above ε ≈ 0.2 now show as held for fewer than 5 τ₁, and their D_c is correspondingly less
certain (the `RELAX_SYS` tail fit is the direct measure of the unrelaxed stress). **Figure 5b, reservoir normal-stress check (2026-10-06, `fig_reservoir_normal_stress[_sweep]`):** σᵗ_xx, σᵗ_yy, σᵗ_zz over the hold in the gel interior and in the two solvent reservoirs (per-snapshot piston-tracking masks, the pore-baseline windows), as plateau increments from the ε = 0 reference. The reservoirs are the in-situ control — a fluid must stay at P_bath in every component — so a lateral rise confined to the gel is network stress (σ'_xx; the box quotes G_comp/M = ½(1 − σ'_xx/σ'_zz)), while a lateral drift in the reservoirs would be a box / normalisation artefact. Sweep: increments normalised by the gel's Δσᵗ_zz per level + the implied G_comp/M against `G_REF`/M_net (Config `G_REF`, the shear notebooks' G). Run 7: reservoirs 0 within ±0.003 at every level while the gel's lateral increment is 0.77–0.97 of the axial one → the compression lateral stress is real network stress, not an artefact.

**D_c fit window and stress-trace D_c (2026-10-09, `fit_Dc` / `fit_Dc_stress`; both compression notebooks, no notebook change).** A hold is not one poroelastic transient. On the quarter gel (comp_3 level 0.30, network strain 0.15, fine file) the first ~6 τ₁ follow D_c = 0.109 — κ = D_c/M = 0.22, at the permeation Darcy κ of 0.24–0.32 — and a slower process then adds +0.012 of interior strain (20 % of the transient) and takes 8 % off the load-piston stress between ~1000 and ~4000 τ; the whole-hold fit with poroelastic modes only absorbs it by halving D_c to 0.054 (the "factor 2" between the hold and the permeation κ), and that fit is unstable in the number of modes (N = 8: 0.021) and the trim. `fit_Dc` therefore fits the **first `DC_FIT_TAU1` (6) τ₁ of the hold**, with τ₁ of the fit itself iterated to self-consistency; it keeps the whole-hold value as `F['Dc_all']`, scans 1 / 2 / 4 / 10 τ₁, the whole hold and 3 / 8 modes (`F['stab']`, one line in `load_level`) and flags a disagreement above `DC_STAB_FLAG` (20 %; `!` in `print_summary`). κ = D_c/M uses the windowed value; `DC_FIT_TAU1 = 0` restores the old behaviour. The independent check is `fit_Dc_stress`: the load piston is solvent-transparent, so F_load/A (`piston_force_avg`, 25-τ blocks) is the network stress at the held face, which relaxes in the same cos(2mπζ) modes — the full held-slab series (300 modes, Gaussian-skin hold-onset state) fitted from `DC_STRESS_TMIN` (250 τ) and a single exponential on the tail from `DC_STRESS_TAIL_TAU1` τ₁; the fit is refused with a reason when the transient precedes the window (the quarter gel, τ₁ ≈ 160 τ). Stored as `L['Dc_stress']`, printed next to the profile D_c with their ratio, drawn in `fig_Dc` (now 2 × 2: (c) the interior strain change vs hold time with the windowed and the whole-hold fits, (d) the stress trace with both of its fits and the profile τ₁ for comparison), in `fig_Dc_sweep` (a) and in `fig_kappa_sweep`. Run 7 (coarse file: 2120-τ blocks, the first ending 2775 τ = 0.6 τ₁ into the hold, so the fast phase sits inside one block and the window changes nothing): profile 0.051 / 0.046 at 0.10 / 0.15 against the stress trace 0.073 ± 0.008 / 0.069 ± 0.013 (τ₁ 4500 / 4200 τ; the 0.15 tail bends to a slower component beyond ~5000 τ), i.e. κ 0.13 → 0.19 against Darcy 0.24–0.32. `simulations/triaxial_compression_two_pist/triaxial_compression_two_pist_fine.batch` is the one-job test on the full gel: run 7's level 0.10 with the fine file, a 10 M-step hold and the new **unload** (below). Housekeeping in the same change: `fig_kappa_sweep` clips an error bar whose M CI crosses zero, and `add_perm_psd` no longer calls `float()` on 1-element arrays (numpy ≥ 2).

**Unload / free re-swelling (2026-10-09; deck Phase 2c `-var unload 1` = `UNLOAD=1`, `load_unload` / `fit_Dc_unload`, no notebook change).** The skin-free control of the hold's D_c: after a level's hold both plates are driven back at the drive speeds to the seated reference gap (`UNLOAD_MARGIN` σ beyond it if set), frozen, the displacement reference reset, and the gel left to re-swell for `UNLOAD_STEPS` (0 = the hold length) with its own files tagged `_<tag>_u<lvl>` (fine + coarse profiles, bounding box, plate positions, load-piston force, reservoir pressures, solvent uptake, σ_zz polymer / solvent); the next level then compresses from the reference again, so levels are matched compression / re-swelling pairs. It starts from the UNIFORM equilibrated compression, so there is no ramp-made skin, and its boundary conditions are the permeation transient's, not the hold's: the plates no longer touch the gel, the network stress at both faces is zero, the strain there is pinned at the reference value (Dirichlet), solvent crosses the faces — modes sin kπζ at k²π²D_c/L², **τ₁ = L²/(π²D_c), four times the hold's** at the same D_c (the hold only redistributes solvent between the compressed skins and the interior with no flux through the plates). The faces detach from the retracting plates within a few τ, so t = 0 of the trace is the retraction start. Two readings: the thickness trace L_bb(t) − L_hold against the exact uniform-IC series ΔL/ΔL∞ = 1 − Σ_odd 8/(kπ)² e^{−k²π²D_c t/L²} with ΔL∞ = L_ref − L_hold fixed (**one free parameter**; the free-ΔL∞ fit is the check), and the profiles since the reset in the compressed gel's material coordinate, COM drift removed, modes cos kπζ with free amplitudes and the same window / stability scan as `fit_Dc`. `UNLOAD_L` picks the L of τ₁ ('hold' = the compressed thickness, as the hold's own fit, so the two compare like for like; 'ref' reads (L_ref/L_hold)² larger). `load_level` prints both next to the hold's windowed / whole / stress-trace D_c with the ratios, `print_summary` gets a `D_c_unload` column, and `fig_Dc` grows a third row — (e) the trace with the series at every D_c reading, (f) the profiles — when the level has an unload. `tests/dc_convention_test.py` integrates this third boundary-value problem (u′ = s at both faces, mean u fixed) by finite differences and both unload readings return the solver's D_c within 0.1 %.

**`triaxial_permeation_single_two_pist.ipynb`** (two-piston permeation, 2026-09-16, standard since 2026-09-22 — no sweep notebook)
`Config(mode="permeation", two_pist=True)`; `tri.load_permeation` builds one dict `P`: total / partial / network
stress and density evolutions (cividis, bold final, zero-flux `_ref` baseline dashed), both wet pistons
(displacement; block-averaged `F_fluid/(lx ly)` from `piston_force_avg` measured vs applied; reservoir virial
pressures dotted, pf cadence since 2026-09-23), the flux `Q_perm` as the **slope of the permeate bead count** over
independent `q_win_steps` windows with its standard error and a slope drift test (`plateau_window_slopes`,
2026-09-23; the `z_perm` slope is drawn beside it; the legacy block-bootstrap mean of the piston-velocity trace is
only printed since 2026-09-24 — the velocity blocks are dominated by the sheet's thermal jitter and their scatter is
not an error bar), the permeability `κ = Q_perm L/(A ΔP_ext)` with the measured ΔP_ext (primary) and the applied one
(CIs in quadrature), and — since
2026-09-23 — the three solvent volume-fraction estimators of the compression notebooks (`tri.add_perm_volume_fractions`:
mass fraction from every density snapshot, Voronoi and λ-calibrated Voronoi on `VOR_EVO_FRAMES` frames over the drive
plus `VOR_MAX_FRAMES` inside the steady window). Because the pore pressure falls feed → permeate across the membrane
under flow, the calibration pressure handed to λ(φ_p, P) is chosen per bin by `Config.P_CAL_MODE`: `'pore'` (the
permeation default — the ramp between the measured feed and permeate reservoir baselines across the membrane),
`'const'` (P_CAL everywhere — the compression default, correct for a drained equilibrium) or `'thermo'` (the local
P_th = −⅓ tr σ^t); the same knob applies to `add_volume_fractions` in the compression notebooks. Since 2026-09-24
the same frames also feed **`lib/psd.py`** (`tri.add_perm_psd`): the geometric, probe-accessible porosity per bin and the
pore-size distribution by the Gelb–Gubbins largest-included-sphere covering (cKDTree distance field on a grid, one
Euclidean distance transform per radius level), reference vs steady state, with `ξ = √κ` printed next to the geometric
mean pore diameter. Eleven figures (`fig_perm_pistons`, `fig_total_stress` — normalised by P_perm —, `fig_partial_stress`,
`fig_network_stress` and `fig_osmotic_pressure` — drawn **blank** under permeation until the deck measures the solvent
chemical potential, since p_pore(z) is not a reservoir baseline under flow —, `fig_perm_volfrac` — reference vs steady
state, the P_local(z) used and the resulting λ(z) —, `fig_perm_volfrac_evolution` (mass fraction + calibrated Voronoi),
`fig_perm_psd`, `fig_perm_flux`, `fig_perm_permeability`, `fig_thermo_pressure`; `fig_perm_density` still exists but left the
notebook). The sync logs in when a trajectory is missing and an unresolved `<steps>` tag matches digits only, so a
permeation pattern can no longer stage a compression level's `_c<lvl>` dump (2026-09-24).
Since **2026-09-28** the notebook also measures **D_c and M from the polymer displacement** (`tri.add_perm_displacement`,
figures 12–14): `disp_z_polymer` (the deck's `fix ave/chunk` of the per-atom z displacement since the reset, every
`disp_nfreq` = 50 k steps, block-averaged) and the bounding-box thickness trace. Model: 1-D poroelasticity
`∂u_z/∂t = q(t) + D_c ∂²u_z/∂z²` with `u_z = 0` at the support and a free feed face; because the pistons prescribe the
*pressure drop* and the total stress is uniform, the network stress at the support equals ΔP from t = 0⁺, so the strain
obeys the diffusion equation with Dirichlet ends and the displacement modes are 1 − cos(kπζ) at k²π²D_c/L₀², τ₁ =
L₀²/(π²D_c); the drive is the applied-ΔP history read from `piston_pressure` (any ramp shape is exact), the reset
convention (ramp start since 2026-09-28, ramp end in older decks) is auto-detected. `fig_perm_Dc`: the profile snapshots
with free modal amplitudes (window-averaged modes, D_c) and the feed-face trace `L_bb(t) − L₀` against the exact zero-IC
series (a second D_c and the asymptote u_F(∞)). `fig_perm_M`: **M = ΔP_ext L₀/(2|u_F|)** from the steady thickness
change (primary), the steady-profile parabola and the trace asymptote, CIs with the measured ΔP_ext in quadrature; the
loader prints κ = D_c/M next to the Darcy κ. `fig_perm_flux_check` (figure 14, a check): the flux the fits imply,
Q(t) = A[(D_c/M) ΔP(t)/L₀ − d⟨u_z⟩/dt], against the measured piston-velocity trace, the N_perm-slope Q and the cumulative
bead count — nothing in it is fitted to the flux. perm_2 (dP = 0.1): D_c 0.18 (profile) / 0.21 ± 0.03 (trace), M 0.48,
κ = D_c/M 0.4 vs Darcy 0.29; the profile fit is under-determined there (six stress-cadence snapshots, the transient
unrecorded — hence the 50 k cadence and the 1 M-step ramp in the deck). Knobs: `PERM_GAP`, `PERM_L0_STEPS`,
`PERM_DC_WINDOW_AVG`, `PERM_DC_BOUNDS`, `PERM_DISP_RESET`, `PERM_TRACE_SKIP` (leave a short-ramp run's piston ringing out
of the trace fit).

**Figure conventions shared by the three two-piston notebooks (2026-10-06)** — two `tri.Config` knobs, both on by default.
`FLIP_Z = True` draws every z-profile figure with the feed / load piston on the left and the permeate / support on the
right (as in Marioni et al.): the `z/L` axes become `1 − z/L`, running 0 → 1 with the two ends named in the label
(feed → permeate, or load piston → support); the σ- and ζ-unit axes are inverted and say so in their labels. It is
plotting only (`zn` / `finish_axes` / `_zlim` / `flip_z_axis`; every mask, window and fit uses z itself), and `False`
reproduces the old figures pixel for pixel. `PARTIAL_NORM = 'share'` makes `fig_total_stress`, `fig_partial_stress`, `fig_thermo_pressure` and their `_sweep`
versions draw a second figure under the original (same cell, same plot number, file stem `…_norm`; the original PNGs are
untouched) with each species' share of the back pressure removed and the result divided by the driving pressure:
`w_s = σ_s/σᵗ` per z-bin and snapshot, `σ_s* = (σ_s − w_s P_ref)/ΔP`, `σ_p* = (σ_p − (1 − w_s) P_ref)/ΔP`,
`σᵗ* = (σᵗ − P_ref)/ΔP` (`tri.partial_norm`; `P_ref`, `ΔP` = the applied permeate pressure and feed − permeate difference
under permeation, `P_BARO` and the level's plateau load-piston increment `dP_pist` in compression; the trace is treated the
same way with `w_s = tr σ_s/tr σᵗ`). The shares come from the partial stresses alone — no mass or volume fractions. It is
the form of Marioni's partial-`P_zz` figure but does not reproduce his end points (perm_3: solvent* 0.50 → 0.38, polymer*
0.49 → 0.61 from the feed face to the support face); see the Notes cells. The normalised trace figure
(`fig_thermo_pressure_norm`, 2026-10-07) also carries the guide `ΔP_th = (4/3)(G/M) ΔP_ext` — the drop of the total `P_th`
from the feed reservoir's value at the feed face to the permeate side when the lateral stresses carry `(M − 2G)/M` of the
axial load (one guide per G estimate: `G_REF`, the shear notebooks' 0.2, and `G_COMP_REF`, the holds' lateral-stress ~0.05; `M` = the run's own primary estimate or `M_REF`; a flat `1 − 4G/3M` in a compression hold).
`calibration_analysis.ipynb` ends with that figure drawn from perm_3 (its own `tri.Config`, "for now" that run) as the
input to the `P_CAL_MODE` question. `'raw'` draws the originals only. `bulk.Config` overrides both knobs to off; the shear
and cross-deck notebooks do not use these figures.

**Tests** (`scripts/tests/`, 2026-09-16; six-plate `compress_slab` fixture added 2026-09-29): `run_plot_tests.sh` builds synthetic one-piston, two-piston and six-plate run
trees (`make_fixtures.py`) and runs the three plotters on both formats; `run_notebook_tests.sh` executes the
two-piston notebooks headlessly on the same tree; `lint_lmp.py` statically checks a deck (definitions before use,
unfix/undump pairing, labels, forwarded `-var`s).

**`bulk_modulus_analysis_single.ipynb`** / **`bulk_modulus_analysis_sweep.ipynb`** (2026-09-29; replace `bulk_modulus_analysis.ipynb`, now in `archive/`)
The bulk-modulus counterparts of `triaxial_compression_{single,sweep}_two_pist.ipynb`, same **Config → sync → load → figures → Notes** layout; all analysis code lives in **`lib/bulk.py`**, which takes its readers, statistics, Terzaghi split, consolidation fit and profile figures from `lib/triaxial.py`. `bulk.Config` points at `flow_data_local/bulk/<RUN_ID>` and Expanse `lammps_runs/compress_slab`. Profiles are read in a core column through the gel along z, x and y (bath | gel | bath). Per level: the measured plate ε_vol and the network's own Rg / BB volumetric strains; **K_network** = plateau-averaged interior ⟨Π⟩ increment from the ε_vol = 0 reference over ε_vol, per profile axis and their mean (exactly `M_network`); **K_plates** = mean plate pressure over the drift-tested block-bootstrap plateau window, increment over ε_vol (exactly `M_piston`); `K_net(Rg)` = the same stress over the network's own strain; δ_sys from the plate tail (`RELAX_SYS`); the isotropy of Δσ′ (there is no G under isotropic loading); apparent D_c per axis; bath check and solvent expelled; with `M_REF` the closure G = ¾(M − K) and κ = D_c/M. `D_c` is the **held-cube fit** (`fit_Dc_cube`: every face drained at the bath pressure, zero displacement at the centre, so `sin 2lπζ` modes only; one `D_c` for the three axes; `DC_N_T` adds the mixed transverse modes as a cross-check). Network anisotropy: `fig_ratio` (σ′_ii/Π over the hold, per axis) and `fig_anisotropy_run` (whole-run `stress_aniso`, optionally beside the `slab_with_support` run's file). Solvent volume fractions in the core columns along x, y, z: mass fraction, Voronoi and λ-calibrated Voronoi (`add_volume_fractions`, one tessellation per frame; `fig_volfrac`, `fig_volfrac_evolution`). Knobs beyond `tri.Config`: `AXES`, `K_STRAIN`, `PLATE_AREA`, `K_SUBTRACT_REF`, `M_REF`, `G_REF`, `HOST`, `DC_N_T`, `core_margin`. `scripts/tests/` builds a six-plate fixture from a known K and asserts the library returns it.

**`volume_of_mixing.ipynb`** *(archived 2026-08 → `lammps_work/archive/`; superseded by the calibration-sweep analysis — see `archive/README.md`)*
Computes ΔV_mix(P*) = V_mixed − V_pure_solvent − V_pure_polymer across the pressure sweep (P* = 1.0–2.0). Cell 2 syncs `box_dimensions_*.dat` files directly from Expanse via `paramiko` SFTP — no SSH keys required; prompts for password and TOTP code in the notebook. Subsequent cells parse the box dimension files, time-average volumes over the last 50% of each run, and plot both ΔV_mix and the individual component volumes vs P*. Requires `paramiko` (`pip install paramiko`). Data lands in `flow_data_local/volmix_sweep/p{P}/`.

**`shear_analysis_single.ipynb`** / **`shear_analysis_sweep.ipynb`** (2026-09-23; replace the 13-step `shear_analysis.ipynb`, git history)
The shear counterparts of `triaxial_compression_{single,sweep}.ipynb`, same **Config → sync → load → figures → Notes** layout; all analysis code lives in **`lib/shear.py`**, which imports the statistics, readers, palette and legend placement from `lib/triaxial.py` so the shear numbers are formed with the same code as the compression numbers. `sh.Config` points at `flow_data_local/shear/<RUN_ID>` and Expanse `lammps_runs/shear_slab` (`lammps_runs/shear_slab_two_pist` with `DECK = "shear_slab_two_pist"`, see the end of this entry). Per level: the plate-based held γ; **G_network** = plateau-averaged interior ⟨σ_p,xz⟩ increment from the γ = 0 reference over γ (t-interval over bins ⊕ reference CI, exactly `M_network`); **G_series** = the fine `stress_series` ⟨σ_p,xz⟩(t) over the drift-tested block-bootstrap plateau window, increment over γ (exactly `M_piston`); **G_total** from σ^t_xz as the poroelastic check; N1/N2; the bulk P_th = −⅓ tr σ^t (the P* = 1.5 check); the transverse D_c fit of u_x(z,t) (even sine modes) and **κ = D_c/G**, the same hydraulic permeability as D_c/M in compression. Ten figures each (single: strain, total stress, partial σ_xz, network σ_p,xz reference + final, stress history, G comparison, N1/N2, D_c fit, κ, P_th; sweep: the overlays plus G vs γ, the stress–strain curve with least-squares and through-origin slopes, N1/N2, D_c, κ, P_th vs level). Needs the `_ref` files and `stress_series` files the deck writes since 2026-09-23 (older runs are reported absolute, flagged). Since the 2026-09-22 re-layout of `shear_slab.lmp` the strain is **plate-based**: γ = (x_plate,top − x_plate,bot − x₀)/plate_sep from the plate COMs (both prescribed by `fix move`), replacing the old bounding-box-gap denominator and outermost-1 %-bead numerator, which overstated the sheared thickness by ~5 % because `add_plates_to_gel` puts each plate ~1–3 σ inside the outermost beads (the bulk region is now cut `plate_excl` inside the plate planes for the same reason). The `shear_strain_*` files keep their four columns (step, gel_lz_initial, **plate_sep**, **γ**); the measured surface-COM strain goes to `shear_strain_surface_*_g<γ>.dat` as a slip diagnostic. The per-level `shear_strain_*_g<γ>.dat` spans the drive **and** the hold (the notebook and `plot_shear_strain_sweep.py` take its last row as the level's γ, i.e. the held strain), the hold is padded to the next stress-averaging epoch, and stem-suffixed whole-run trackers (`box_dimensions`, `gel_dimensions_rg`, `gel_volume_{bb,rg}`, `shear_strain`) feed `plot_lammps_log.py`'s convergence and shear-diagnostics figures; the files and columns the notebook reads are unchanged. To compare with the two-piston compression G, use the same bath pressure (P* = 1.5: the converter's calibrated solvent deletion around the plates puts the gel there, and Step 7's bulk `<p_total>` is the check — the thermo `press` of this box is not the gel pressure) and the G = (σ′_zz − σ′_xx)/2ε values from `triaxial_compression_sweep_two_pist.ipynb`.

*Two decks, one pair of notebooks (2026-10-03).* `DECK` in the Config cell selects the run type: `"shear_slab"` (default, plate-driven) or `"shear_slab_two_pist"` (the two-piston system with rigid-in-plane face layers). Same `_g<strain>` files, columns and figures; for the two-piston deck the sync points at `lammps_runs/shear_slab_two_pist`, the dash-dot planes are the gel's Rg faces, the bulk starts 5 σ inside them (`plate_excl` 5, `wall_margin` 6: driven layer + tracking slab), the solvent-only "film" used for p_pore is the two reservoirs (full bins only), `G_plate` is the face-layer force per area, the D_c fit maps the bulk bins onto the driven-layer COM planes (`H_drv` apart), and there is no unload, so the residual figure skips itself (it also does for `shear_slab` runs without `-var unload 1`, the default since 2026-10-03). *Figure 2b (both notebooks, 2026-10-05, `fig_ref_normal_stress`):* the three total normal stresses **before shearing** — bulk σᵗ_xx, σᵗ_yy, σᵗ_zz and P_th over the γ = 0 reference window with the lateral pre-stress ⟨σᵗ_xx − σᵗ_zz⟩, ⟨σᵗ_yy − σᵗ_zz⟩ (block bootstrap), and their window-mean z-profiles. The reservoirs must read P_bath in all three; a lateral component away from σᵗ_zz inside the gel would be a pre-stress of the slab (transversely isotropic network, so the compression (C₃₃ − C₁₃)/2 and the shear C₄₄ would differ). shear_3 (two-piston) and film7 (plates): both differences are 0 within ±0.001–0.0016.

*The strain G divides by — `G_STRAIN` (2026-10-03), the shear counterpart of `M_STRAIN`.* G and κ divide by the strain of the **network**, not the prescribed boundary strain: `'disp'` (default) = the interior slope of the steady **cumulative** displacement profile u_x(z) (`disp_x_polymer_cum_*`: plateau mean minus the γ = 0 reference-window profile, line fit over the populated bins between the planes with `DISP_TRIM_BINS` = 3 dropped per face; CI = max of the least-squares and frame-to-frame intervals) — z does not change under shear, so the current-z bins need no conversion. Fallback chain when a source is missing: `'rg'` = (Rg_xz − ⟨Rg_xz⟩_ref)/Rg_zz² from `shear_strain_rg_*` → `'trk'` = the tracking-slab strain (two-piston deck: hold plateau minus its γ = 0 reading) → `'applied'` = the prescribed plate / driven-layer strain. The chosen strain's CI is added in quadrature to every G (as ε_M's is to M). Runs of decks older than 2026-10-03 have neither file and fall back to the prescribed strain with a printed FALLBACK note, so their G is unchanged. New figure **1b** (`sh.fig_disp_profile`): the profile(s), line fit and residuals, with every strain measure listed; `print_summary` lists them per level together with the G the prescribed strain would give. `L['gamma']` is the strain used, `L['gamma_applied']` the prescribed one.

---

**`cross_deck_comparison.ipynb`** + **`lib/compare.py`** (2026-10-03)

Figures built from **more than one deck's run**. One Config per deck (copied from that deck's own notebook) goes into `CFG`; `cmp.load_deck` loads each run with its own library (`lib/triaxial.py` for the two-piston compression and permeation decks, `lib/shear.py` for `shear_slab`) into `DECKS[name] = dict(kind, cfg, R, levels | P)`; comparison functions live in `lib/compare.py` and take the decks they need. Plots go to `flow_data_local/plots/comparison/<NAME>`.

First comparison, `cmp.phi_dependence`: **D_c(φ_s) and κ(φ_s)**. The compression levels are uniform states: they give σ′(φ_s), the current-frame tangent modulus `M_cur = φ_p dσ′/dφ_p`, and per level `κ = D_c/M_cur`. The steady permeation membrane spans a range of φ_s at one flux: locally `D_c = q φ_p/(dφ_s/dz)` (flux and profile slope only) and `κ = D_c/M_cur`. `PHI_KIND` picks the volume fraction (`'cal'` λ-calibrated Voronoi, `'vor'`, `'mf'`). Run 7 + perm_3: the two decks agree on σ′(φ_s) to ~6 %, but at the same φ_s the steady permeation D_c (≈ 0.09) is ~2× the hold-relaxation D_c (≈ 0.05).

Opening figures (2026-10-08, permeation deck only, before the φ_s comparison): `cmp.fig_thermo_pressure_norm` draws the permeation notebook's normalised thermodynamic-pressure figure with its two `ΔP_th = (4/3)(G/M) ΔP_ext` guides (`G_REF` = 0.2 from the shear decks, `G_COMP_REF` = 0.05 from the compression holds; the deck Config in the notebook carries them) and copies the PNG into the comparison folder; `cmp.fig_solvent_measures` is figure 7a of the permeation notebook recast as **concentration** `c_s = φ_s^cal ρ_s,0/M_s` (left axis, beads/σ³; `M_s` = 1) against the **mass (= number) fraction** `ρ_s/ρ_s,0` and the **probe porosity** `ε_g` (right axis, scaled by `ρ_s,0`), reference dashed / steady solid with 95 % bands, the individual steady Voronoi frames of `c_s` faint (`CONC_FRAMES`), the reservoirs trimmed to `CONC_MARGIN` σ (feed side: up to the piston), and a panel (b) with the two other measures minus `φ_s^cal` (`CONC_DIFF_PANEL`). The deck Config runs the PSD pass with `PSD_R_PROBE` = 0.45 (the permeation notebook uses 0.5), so `load_deck` now calls `tri.add_perm_psd` when `PSD_ENABLE` is set. perm_3: interior `φ_s^cal` 0.39, mass fraction 0.42 (+0.03: solvent packed on the polymer that the calibration assigns to the network), `ε_g` 0.23 (−0.16: the excluded shell), feed-reservoir mass fraction 1.03 (denser solvent under the feed pressure).

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
SCRIPTS=~/Documents/lammps_work/scripts      # the plotters live in $SCRIPTS/lib/ (since 2026-10-01)
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
python "$SCRIPTS/lib/plot_lammps_log.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}"
```
No `--p-ext` flag needed for any current folder — it was a `slab_with_flow`-only option (now removed).

Output saved to `./output_plots/`.

---

#### `plot_stress_profiles.py` — partial stress and volume fraction profiles

```bash
python "$SCRIPTS/lib/plot_stress_profiles.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0
```

---

#### `plot_piston_data.py` — piston position, velocity, force (`triaxial_compression` single-level runs, `triaxial_permeation`, etc.)

```bash
python "$SCRIPTS/lib/plot_piston_data.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0
```

---
