# Claude Context — LAMMPS / Hydrogel Poroelasticity Research

> **How to use this file**: At the start of any new Cowork session, tell Claude:
> *"Read CLAUDE_CONTEXT.md and let's continue."*
> Claude will read this file and pick up where we left off.
> Update this file at the end of each session (or ask Claude to update it).

---

## Researcher

- **Name**: Dylan Pollard
- **Email**: pollard@ucsb.edu
- **Institution**: UCSB

---

## Project Overview

Coarse-grained molecular dynamics simulations of **tetrahedral hydrogel slabs** using LAMMPS. The goal is to extract **poroelasticity observables** from the simulations — primarily pore pressure, volume fraction profiles, cooperative diffusivity, flux, and elastic moduli — to characterize hydrogel transport and mechanical properties.

**System**: Coarse-grained polymer gel + solvent beads in LJ reduced units. Slabs contain polymer chains, crosslinkers, solvent, support beads, and piston beads. Typical system sizes: 200k–1.5M atoms. Long runs (up to 20M timesteps) on HPC clusters.

---

## Research Milestones

1. **Minimize volume of mixing** (incompressible regime)
   - ✅ Done: Found ~1% volume change at P* = 1.5

2. **Verify uniform pore pressure profiles at thermal equilibrium**
   - ⬜ Edit LAMMPS script to ensure P* = 1.5 from the beginning
   - ⬜ Match trajectory output frequency with volume measuring frequency (reduces data processing time)

3. **Calculate M (longitudinal modulus)**
   - ✅ Done: Found from network stress (polymer partial stress)
   - Script: `scripts/compression_analysis.ipynb`

4. **Calculate G (shear modulus)** ← *CURRENT FOCUS*
   - Method: Plate-driven shear (±z on x-face plates); record gel stresses in Phase 3
   - Verify G decreases with decreasing crosslink density
   - Simulation script: `simulations/shear_slab/shear_slab.lmp` (finalized 2026-04-28)
   - Input: isolated swollen gel + attached plates from `scripts/add_plates_to_gel.ipynb`
   - G = ⟨σ_p_xz⟩ / γ_final from Phase 3 time-averaged stress tensor
   - Analysis notebook: `scripts/shear_analysis.ipynb` ← created 2026-04-28
   - Still need: run shear_slab, verify G vs crosslink density series

5. **Run permeation simulations**
   - Add reservoir pressure measurements to find ΔP
   - Add solvent bead counter to find flux
   - Find Dc from flux × L / Δφ
   - Simulation type: `slab_with_flow` with `compression_mode = 0` (permeation mode)

6. **Calculate ΔPth and verify ΔP = ΔP + ΔΠ**

---

## Folder Structure

```
lammps/                          ← workspace root (synced to MacBook via iCloud)
├── lammps_work/                 ← Git repo (synced to GitHub)
│   ├── simulations/
│   │   ├── slab_with_support/   ← equilibration/compression sims
│   │   ├── slab_with_flow/      ← permeation + compression sims (slab_with_flow.lmp)
│   │   ├── shear_slab/          ← shear modulus sims
│   │   │   ├── shear_slab.lmp   ← xz shear sim; 3 phases: NVT + shear + NVT production
│   │   │   └── shear_slab.batch ← Bridges-2 batch script
│   │   ├── polymer_phase/
│   │   ├── polymer_pure/
│   │   ├── solvent_phase/
│   │   └── solvent_pure/
│   ├── scripts/
│   │   ├── compression_analysis.ipynb    ← M, Dc, pore pressure, volume fractions
│   │   ├── permeation_analysis.ipynb     ← flow profiles, pore pressure evolution (future)
│   │   ├── shear_analysis.ipynb          ← G, N1, N2, profiles, poroelastic decomp (created 2026-04-28)
│   │   ├── add_plates_to_gel.py          ← builds *_with_plates.data for shear_slab
│   │   ├── add_plates_to_gel.ipynb
│   │   ├── add_walls_to_slab.ipynb
│   │   ├── add_walls_with_angles.ipynb
│   │   ├── slab_with_support.ipynb
│   │   ├── slab_with_support_angled.ipynb
│   │   ├── isolate_gel.ipynb
│   │   ├── split_gel_slab.ipynb
│   │   ├── pure_polymer.ipynb
│   │   ├── pure_solvent_1.ipynb
│   │   ├── plot_lammps_log.py
│   │   ├── plot_stress_profiles.py
│   │   ├── plot_piston_data.py
│   │   ├── git_sync.sh
│   │   ├── run_lammps.sh / run_lammps_pod.sh / run_lammps_bridges.sh
│   │   └── write_tracking.py
│   ├── lammps_data/             ← small data files (in repo)
│   ├── .gitattributes           ← strips notebook outputs before git commit (nbstripout)
│   ├── README.md
│   ├── expanse_lammps_guide.md
│   ├── slab_data_file_info.md
│   ├── documenting_pod_runs.md
│   ├── lj_units_cheat_sheet.md
│   └── git_setup_instructions.md
│
├── flow_data_local/             ← local copy of HPC output (NOT in git)
│   ├── compression/             ← slab_with_flow compression_mode=1 output
│   │   ├── rho04_p1.5_600k_2M/ ← canonical run (P*=1.5, current)
│   │   ├── rho04_500k_2M/      ← older run (P*=1.0)
│   │   └── rho04_p1.05_15M_10M/← oldest run
│   ├── shear/                   ← shear_slab output (future)
│   │   └── <RUN_ID>/
│   ├── permeation/              ← slab_with_flow compression_mode=0 output (future)
│   │   └── <RUN_ID>/
│   ├── plots/                   ← auto-saved analysis plots
│   │   ├── compression/<RUN_ID>/
│   │   ├── shear/<RUN_ID>/
│   │   └── permeation/<RUN_ID>/
│   └── traj_files.nosync/       ← large .lammpstrj files (iCloud nosync)
│
├── lammps_data_files_local/     ← large .data files (NOT in git, local only)
│
├── add_walls_to_slab.ipynb      ← (root-level, older versions)
└── slab_with_support_4.ipynb
```

## Analysis Notebook CONFIG Pattern

Every analysis notebook has a config cell (cell 3) at the top — **only these two lines change between runs:**

```python
RUN_ID   = "rho04_p1.5_600k_2M"   # folder name inside flow_data_local/<sim_type>/
sim_name = "walled_slab_support_5beads_tall_rho04_p1.5_..."  # LAMMPS file prefix
```

Path variables are then derived automatically:
- `DATA_DIR`  = `../../flow_data_local/<sim_type>/<RUN_ID>/`  (input .dat files)
- `PLOT_DIR`  = `../../flow_data_local/plots/<sim_type>/<RUN_ID>/`  (auto-created)
- `TRAJ_FILE` = `../../flow_data_local/traj_files.nosync/gel_flow_<sim_name>.lammpstrj`

## nbstripout (notebook output stripping)

Notebook outputs are automatically stripped before every git commit via `nbstripout`.
One-time setup on each machine:
```bash
pip install nbstripout
cd ~/docs/grad_research/lammps/lammps_work
nbstripout --install
```
The `.gitattributes` file in the repo root ensures the rule is applied everywhere.

---

## HPC Clusters

| Cluster   | Usage                        | Notes                                      |
|-----------|------------------------------|--------------------------------------------|
| Bridges-2 | Long runs (10M–20M steps)    | 120 MPI tasks/node optimal; 10 GB home quota; traj → /ocean/projects/chm250028p/dpollard/ |
| Expanse   | Medium runs; varies by queue | ~73 steps/s at 640 MPI × 5 nodes          |
| Pod       | Short/misc runs (used least) | 100 Gb/s internode (slower than Bridges-2) |

**Download trajectories**:
```bash
rsync -avP $USER@data.bridges2.psc.edu:<PATH>.lammpstrj.gz ~/Downloads
```

---

## Simulation Types

### `slab_with_flow.lmp` (primary sim for poroelasticity)
- **compression_mode = 0**: Permeation — piston pushes solvent through gel, measure flux
- **compression_mode = 1**: Compression — piston loads gel directly, solvent squeezed out
- Atom types: 1=polymer, 2=crosslinker, 3=solvent, 4=support, 5=piston, 6=walls
- Current P* target: 1.5

### `slab_with_support.lmp`
- Equilibration and compression runs
- Used for M calculation (longitudinal modulus)

### `slab_elongation.lmp`
- Will be used for G (shear modulus) calculation

---

## Key Parameters (LJ reduced units)

- `sigma = 1.0`, `temp_target = 1.0`, `P_target = 1.5`
- `timestep = 0.005` → 1e6 steps ≈ 4.1e3 τ
- Solvent density explored: ρ = 0.2–1.0
- Crosslinker configs: rho02, rho03, rho04, rho06, rho08

---

## Workflow (Day-to-Day)

1. Edit LAMMPS `.lmp` / `.batch` scripts in `lammps_work/simulations/`
2. Commit and push to GitHub (sync to clusters via `git pull`)
3. Submit jobs via `sbatch` on cluster
4. Copy output `.dat` files from cluster → `flow_data_local/` on MacBook
5. Analyze in JupyterLab using notebooks in `lammps_work/scripts/`
6. Update documentation (README.md, expanse_lammps_guide.md)

---

## shear_slab.lmp — Design Notes (updated 2026-04-28)

- **Input data file**: isolated gel WITH attached shear plates, produced by `scripts/add_plates_to_gel.ipynb`.
- **Atom types**: 1=polymer, 2=crosslinker, 3=solvent, 4=plate
- **Bond types**: 1=FENE (polymer-polymer), 2=harmonic k=30 r0=1 (plate-polymer)
- **Phases**:
  - Phase 1a (NPT, 50k steps): plates frozen (`fix move linear 0 0 0`); gel runs `fix npt iso P_target`; box x adjusts to reach P=1.5. Plate atoms scale with box (reduced coords fixed).
  - Phase 1b (NVT, 100k steps): box locked at NPT-final volume; plates frozen; gel thermalises.
  - Phase 2 (shear, up to 570k steps): plates driven ±vshear in z; `fix halt` stops run when `gel_strain_cm ≥ target_strain_xz` (10%).
  - Phase 3 (NVT production, `nsteps` from bash): plates frozen at Phase 2 final positions; gel NVT. G = ⟨σ_p_xz⟩_bulk / gel_strain_cm_final.
- **Shear method**: **`fix move linear 0 0 ±vshear`** on `plate_right` / `plate_left`. Right plate (+x face): +vshear in z. Left plate (−x face): −vshear in z. Plates translate purely in z (vx=vy=0 → no rotation). Harmonic bonds transmit shear from plate → surface polymer → interior network via FENE. Shear plane is **xz** (gap=x, shear direction=z).
- **Plate orientation**: plates are flat **yz-planes** (normal to x), placed at x = gel_xlo − 0.5σ and x = gel_xhi + 0.5σ. The gel is tall in z and short in x; plates tile the full y-z box face.
- **Plate grouping**: `plate_right` = type 4 atoms with x > xmid; `plate_left` = type 4 atoms with x < xmid. Groups defined at startup from region + intersect (static).
- **add_plates_to_gel.ipynb**: square lattice spacing 1.5σ; offset 0.5σ from polymer surface; each plate atom bonded (harmonic) to nearest polymer atom within 2.5σ cutoff (~90–93% bonding rate). ~2550 plate atoms per plate (rho04 tall dataset). Output: 4 atom types, 2 bond types.
- **vshear**: `1.1 × target_strain × gel_gap / (2 × nsteps_shear_base × dt)`. `nsteps_shear_base = 285000` anchors the speed; `nsteps_shear = 570000` is a 2× safety buffer. `fix halt` fires at gel_strain_cm ≥ 0.10 → Phase 2 always stops at exactly the target strain.
- **Strain measurement**: `gel_strain_cm = (zcm_right − zcm_left − zcm_sep_0) / gel_gap`. Tracking groups `poly_right_s` / `poly_left_s` are outermost 1% of backbone (in x) on each side — bonded directly to the plates → accurate surface displacement.
- **Bulk stress region**: atoms within 5σ of each plate face excluded from all stress calculations (`poly_bulk` / `solv_bulk` groups, x ∈ [gel_xlo+5, gel_xhi−5]). Stress normalised by **bulk volume** = `bulk_lx × ly × lz` (not full box volume).
- **Group structure**: `polymer` (1,2), `backbone` (1), `crosslinks` (2), `solvent` (3), `gel` (1,2,3), `plate` (4), `plate_right`/`plate_left` (type 4, split at xmid), `poly_right_s`/`poly_left_s` (backbone, outermost 1% in x), `poly_bulk`/`solv_bulk` (bulk interior, 5σ from walls).
- **Minimize**: plates frozen with `fix setforce 0 0 0` + `velocity plate set 0 0 0` during minimize (5000/50000 steps). Allows gel to relax around plates before dynamics.
- **Trajectory dumps**: both `traj_setup` and `traj_prod` dump **`all`** atoms (includes plate type 4 → plates visible in OVITO).
- **Bond style**: `bond_style hybrid fene harmonic`; `bond_coeff 1 fene 30.0 1.5 1.0 1.0`; `bond_coeff 2 harmonic 30.0 1.0`
- **Pair style**: WCA for all pairs including plate (lj/cut 1.122, ε=1, σ=1)
- **Gel gap (x)**: gel_gap ≈ 47.9σ (rho04 tall dataset, polymer surface-to-surface)
- **Triclinic setup**: `change_box all triclinic` before any dump (LAMMPS restriction).
- **Output frequencies**: `thermo_freq=25000`, `dump_freq=25000`, `stress_freq=10000`
- **Stress outputs**:
  - `stress_tensor_polymer/solvent_*.dat`: bulk-integrated 6-component tensor vs time (bulk_vol normalised); σ_xz is col 5 → used for G
  - `stress_profile_z_polymer/solvent_*.dat`: z-binned 6-component profile (bulk atoms only, bin volume = bulk_lx × ly × binWidth)
  - `shear_strain_*.dat`: step, gel_lz_initial, gel_gap, gel_strain_cm
  - `box_dimensions_*.dat`: step lx ly lz xy xz yz
  - `gel_dimensions_rg_*.dat`: step lx_rg ly_rg lz_rg
  - `polymer_com_*.dat`: step polymer_com_x _y _z
- **G extraction**: G = ⟨stress_p_xz⟩ / gel_strain_cm_final (from stress_tensor_polymer, col 5)
- **Optimal node count**: 1–2 nodes on Bridges-2/Expanse (~171k atoms). Rule: ≥1000 atoms/MPI task.
- **c_gel_press note**: `compute gel_press all pressure gel_temp pair bond` reports **virial-only** pressure (~0.86) because `pair bond` keywords exclude the kinetic NkT/V term. Full pressure is `press` in thermo (~1.5). Both are real; they measure different things.

## shear_analysis.ipynb — Design Notes (created 2026-04-28)

- Located at `scripts/shear_analysis.ipynb`
- Reads output from `shear_slab.lmp` Phase 3 production run
- Planned analyses: G from ⟨σ_p_xz⟩/γ, normal stress differences N1/N2, z-profile stress plots, poroelastic decomposition (polymer + solvent contributions), G vs crosslink density series

## Active Tasks / Open Questions

- [ ] Run `add_plates_to_gel.ipynb` on each isolated gel data file to generate `*_with_plates.data`
- [ ] Run shear_slab on 1–2 nodes and copy output to `flow_data_local/shear/`
- [ ] Populate `shear_analysis.ipynb` once first results arrive; verify G vs crosslink density series
- [ ] Milestone 2 cleanup: align traj/volume output frequencies
- [ ] Keep README.md and expanse_lammps_guide.md up to date
- [ ] Delete `slab_elongation/` folder (superseded by shear_slab)

---

## GitHub Sync Workflow

**Script**: `scripts/git_sync.sh` — one command to pull → stage → commit → push.

**Setup (run once on MacBook)**:
```bash
# Make the script executable (already done)
chmod +x ~/docs/grad_research/lammps/lammps_work/scripts/git_sync.sh

# Add alias to ~/.zshrc
echo "alias lsync='bash ~/docs/grad_research/lammps/lammps_work/scripts/git_sync.sh'" >> ~/.zshrc
source ~/.zshrc
```

**Daily usage on MacBook**:
```bash
lsync                     # auto-generates commit message from changed files
lsync "added shear sim"   # custom commit message
```

**On clusters (Expanse / Bridges-2 / Pod)** — pull before each job:
```bash
cd ~/Documents/lammps_work && git pull
```
Optionally prepend `git pull` to `run_lammps.sh` / `run_lammps_bridges.sh` so it happens automatically before every LAMMPS run.

**What .gitignore already excludes**: `*.lammpstrj`, `*.restart*`, `log.lammps`, `*.log`, `output_files/`, `traj_files/`, `final_config_*.data`, `final_flow_*.data`, slurm output files (`slurm_*.out`, `*.o*.*`, `*.e*.*`), `.DS_Store`, `__pycache__/`, `.ipynb_checkpoints/`.

**Remote**: `https://github.com/Twildridge/lammps_work.git`

## Automation / Efficiency Goals

- ✅ GitHub sync: `lsync` alias on MacBook; `git pull` on clusters
- Better organize `flow_data_local` and `lammps_data_files_local`
- `tracking.txt` is no longer actively maintained — consider deprecating or replacing
- Any non-thinking tasks (file management, docs updates, routine analysis) are candidates for automation via Cowork

---

## Notes / Decisions Log

- *2026-04-15*: First Cowork session. Created this context file. Migrating from Claude web "MD Work" project to Cowork for better file integration. `tracking.txt` is legacy — not actively used.
- *2026-04-15*: Created `simulations/shear_slab/shear_slab.lmp` and `shear_slab.batch`. New `shear_slab` folder replaces `slab_elongation` (old folder still present — delete manually). Shear sim uses xz box deformation + ±x forces on top/bottom 5% of polymer; outputs all 6 stress components for both polymer and solvent (box-integrated tensor + z-profiles), plus bulk pressure from thermo_style. BB (bounding box) computes removed; only minz/maxz kept for shear layer definitions. `target_strain_xz` set to 0.10; batch script set up for Expanse.
- *2026-04-16*: Debugged `shear_slab.lmp` variable formula errors. Two issues found: (1) `xztilt` is not a valid LAMMPS keyword anywhere → replaced with `xz` throughout; (2) LAMMPS cannot parse `$(c_maxz - ${shear_frac}*(c_maxz - c_minz))` — complex `$()` expressions with repeated compute references plus substituted float literals fail. Fix: use already-snapshotted constant variables (`v_gel_zhi`, `v_gel_zlo`, `v_gel_thick`) instead → `variable top_zlo equal $(v_gel_zhi - v_shear_frac*v_gel_thick)` and same pattern for `bot_zhi`. Rule: inside `$()`, keep expressions simple; use `v_` references to pre-evaluated constants rather than `c_` references with repeated arithmetic.
- *2026-04-16*: Fixed `shear_slab.lmp` triclinic box error. `fix deform xz erate` requires a triclinic box; the isolated gel data file creates an orthogonal box. `change_box all triclinic` cannot be called while a dump is active (LAMMPS restriction). Fix: moved `change_box all triclinic` to Phase 1 setup, right before `dump traj_setup` is defined. NPT iso works fine on triclinic boxes. Rule: `change_box` must precede all dump definitions.
- *2026-04-16*: Removed `addforce` on polymer top/bottom layers from shear_slab. `fix deform` with `remap x` already applies affine strain to all atoms — addforce was redundant and introduced non-uniform internal stresses. Pure strain-control (fix deform only) is the correct approach for elastic modulus measurement.
- *2026-04-16*: Fixed `compute chunk/atom bin/1d` error for triclinic boxes: requires `units reduced`. Added `variable binWidth_reduced equal $(v_binWidth/lz)` snapshotted after Phase 2 (lz frozen). Rule: inside `$()`, use `v_varname` not `${varname}`.
- *2026-04-16*: Updated `plot_lammps_log.py`: added argparse with `--run-id` flag; added shear diagnostics plot (auto-triggered by `Pxz`/`Xz` thermo columns or shear output files); fixed `UnicodeDecodeError` in log parser with `errors='ignore'` (LAMMPS echoes non-ASCII chars from print statements into log); added `marker='o', markersize=3` to all plot lines so sparse data is visible. Same encoding fix applied to `write_tracking.py`. Updated README §7a with concrete replot commands using full run_id format.
- *2026-04-16*: HPC scaling for shear_slab: ~165k atoms → optimal at 1–2 nodes (120–240 MPI tasks, ~700–1400 atoms/task). At 4 nodes (480 tasks) → ~345 atoms/task: communication-dominated, slower. slab_with_flow scales better (larger system, optimal ~5 nodes) because more atoms/task keeps compute > communication. Rule: target ≥1000 atoms/MPI task.
- *2026-04-16*: "Stale file handle" MPI error at Phase 3 start is a transient NFS filesystem issue on the cluster, not a code bug. Rerun the job.
- *2026-04-15*: Created `scripts/git_sync.sh` — one-command pull/commit/push with auto-generated commit messages. Add `lsync` alias to `~/.zshrc` on MacBook (and `~/.bashrc` on each cluster).
- *2026-04-15*: Pod VPN (Ivanti Secure Access) troubleshooting — broken install required clearing remnants (`sudo rm -rf` of app + support files) then **rebooting before reinstalling**. Reboot fixed it. Download from https://it.ucsb.edu/ivanti-secure-access-campus-vpn/get-connected-campus-vpn. If it won't open after reinstall: reboot first.
- *2026-04-17*: **Gel rigid-body rotation** — isolated gel in a snug box has no friction against walls; any shear perturbation allows rigid-body rotation. Attempted fixes in order, each failing: (1) `fix momentum linear 0 0 1` (removes z-drift but not rotation); (2) `fix setforce NULL NULL 0.0` on driven layers (caused FENE cascade — zeroes forces violate Newton's 3rd law); (3) `fix addforce ±x` equal-and-opposite (over-shearing + jagged surfaces — see below); (4) **final solution**: `fix move linear` rheometer — bottom layer frozen at v=(0,0,0), top layer at v=(vshear,0,0). Bottom physically cannot translate → rotation impossible. z=0 on both layers eliminates z-drift. FENE safe (forces computed normally).
- *2026-04-17*: **`fix addforce` failure modes** — applying a constant force per atom to top/bottom driven layers causes two problems: (1) Each atom responds independently → uneven elastic displacement → jagged surfaces instead of flat rheometer plates. (2) Elastic equilibrium strain ≈ F/(G·A) ≈ 25% for our gel, far exceeding the 10% target. Box tilt reached |xz| > lx/2, violating LAMMPS triclinic stability and causing "Bond atom missing in image check" errors. Solution: replace with `fix move linear` for exact strain control.
- *2026-04-17*: **`fix halt` timer persistence bug** — `fix halt` internally calls `timer->force_timeout()`, which persists across `run` commands even after `unfix halt_shear`. Result: Phase 3's `run ${nsteps}` exits immediately at 0 steps, producing no thermo output and no trajectory. Rule: **never use `fix halt` in multi-phase scripts with sequential `run` blocks**. Solution: removed `fix halt` entirely; Phase 2 runs full `nsteps_shear` steps.
- *2026-04-17*: **`group intersect` with dynamic groups** — LAMMPS forbids `group intersect` when either operand is a dynamic group. Also: `group intersect` takes **group IDs only** — region IDs cannot be passed directly. Pattern: always create an intermediate static group via `group tmp_name region reg_name` first, then `group final intersect groupA tmp_name`. `polymer_top` / `polymer_bot` (defined by `dynamic/chunk`) are dynamic → error "Cannot intersect groups using a dynamic group". Fix: capture static region snapshots first (`group poly_top_all region reg_driven_top`), then intersect two static groups (`group poly_top_s intersect polymer poly_top_all`). The `polymer` group is always static (defined by atom type).
- *2026-04-17*: **Nested `${}` inside `$()` is invalid LAMMPS syntax** — `$(${target_strain_xz}*lz)` causes "Invalid syntax in variable formula". Fix: use `v_varname` reference inside the `$()` expression instead: `$(v_target_strain_xz*lz)`. Rule: inside `$()`, reference variables as `v_name`, never as `${name}`.
- *2026-04-17*: **`fix print append` is not a valid keyword** — LAMMPS `fix print` has no `append` option; the keyword causes "Illegal fix print command" error. To append across phases without re-opening the file, keep the fix running continuously (do not unfix and redefine mid-script). Alternatively, accept that Phase 2 and Phase 3 write separate files.
- *2026-04-17*: **`shear_strain_*.dat` now has 4 columns** (updated from 3): step, gel_lz_initial, xz_tilt (box tilt from `$(xz)`), gel_strain (computed from top/bot CM separation via `v_gel_strain`). Written during Phase 2 only (fix unfixed at end of Phase 2; Phase 3 does not re-define it). γ = col4 directly, or col3/col2 as cross-check.
- *2026-04-17*: **`fix move` + `fix deform remap none` synchronization** — set `vshear = erate × lz` (evaluated once at Phase 2 start with `variable vshear equal $(v_erate*lz)`). Both top-surface CM displacement and box xz increase at `erate × lz` per unit time → they are exactly synchronized. After `nsteps_shear` steps both equal `target_strain_xz × lz`. No `change_box` correction needed.
- *2026-04-17*: **`fix move` rheometer → FENE bond failure at driven/interior interface** — Observed in log: bond atoms 18841–45089, length growing 1.423036 → 1.4231 → ... monotonically at ≈ vshear/step (0.01138 σ/τ ÷ 0.005 = 0.00569 σ/step × 2 ≈ measured rate). Root cause: `fix move` overrides velocity AFTER force integration; the FENE restoring force is computed correctly but the driven atom's velocity is reset to vshear regardless → bond stretches without bound. Fix: **switch to `fix deform remap x`** (affine remapping of ALL atoms). Adjacent bonded atoms (Δz ≈ 0.1–0.5σ) accumulate at most target_strain_xz × 0.5σ = 0.05σ of relative stretch total — nowhere near FENE R0 = 1.5σ. Driven/interior interface eliminated entirely.
- *2026-04-17*: **`remap x` + `fix momentum linear`** — With `remap x`, COM drifts in x at erate × z_COM × dt/step; `fix momentum linear 1 1 1` (every 100 steps) removes this. `fix momentum angular` is NOT used in Phase 2 because a shear flow has net angular momentum (L ∝ erate × M × Rz²) — zeroing it would remove the shear gradient itself. In Phase 3, `fix mom_zero all_mobile momentum 100 linear 1 1 1 angular` removes both COM drift and residual angular momentum from Phase 2 without contaminating σ_xz (virial is computed from forces, not velocities).
- *2026-04-18*: **`gel_strain` corrected to `xz/lz`** — Previous formula `xz/v_gel_lz_initial` (where gel_lz_initial is the Rg-based gel thickness ≈ gel_thick) overstates γ by factor lz/gel_thick ≈ 1.018 because lz = gel_thick + 2σ buffer. With `remap x`, xz = erate × lz × t exactly, so xz/lz = erate × t = target_strain_xz. The Rg-based gel_lz_initial is still stored as a reference but is no longer used in the strain formula.
- *2026-04-18*: **`fix deform remap x` fails for elastic gels** — gel COM strain saturates at ~45% of target (observed: 0.045 for γ_target=0.10). Root cause: FENE/crosslink restoring forces generate velocities opposing each step's remap displacement; NVT propagates these velocities. Elastic solid resists continuous rate-controlled shear. Fix: replace Phase 2 with `change_box all xz final ${target_xz} remap units box` (one instantaneous affine step) followed by `run 0` + NVT at LOCKED box geometry. The x-periodic bond topology (bonds crossing x-PBC carry image offset of xz in x) enforces gel_strain = target_strain_xz — NVT cannot un-shear the gel without an NPT barostat. FENE safety: max per-bond extension at 10% shear ≈ 0.097σ (6% of R0=1.5σ — safe). `nsteps_shear` now = NVT relaxation steps at fixed geometry (was: ramping steps with fix deform).
- *2026-04-18*: **Added polymer COM gel strain tracking** — `poly_top_s` / `poly_bot_s` groups (top/bottom 10% of gel by z, static) + `compute xcm_top/xcm_bot com` + `xcm_sep_0` snapshot → `gel_strain_cm = (xcm_top_x − xcm_bot_x − xcm_sep_0) / gel_thick`. Now 4th column of `shear_strain_*.dat`. With change_box approach, gel_strain_cm starts at ~target_strain_xz and partially relaxes (expected for surface layers); interior network held by periodic topology.
- *2026-04-18*: **`change_box` approach abandoned** — reverted to change_box + NVT after fix move failures, but Dylan noted that the gel is surrounded by a thin solvent film. Instantaneous box tilt pushes the solvent shell; the gel translates sideways rather than shearing internally. change_box works for systems without solvent film or with periodic gel networks, but not here.
- *2026-04-19*: **Phase 2 momentum fix split** — `fix mom_shear interior_mobile momentum 100 linear 1 1 1 angular` was causing visible +x COM drift because gradient atoms (28319 top vs 28530 bottom) had unequal net momentum not covered by interior_mobile group. Also, `angular` on all_mobile would strip shear-flow angular momentum from gradient zone. Fix: two separate fixes — `fix mom_lin all_mobile momentum 100 linear 1 1 1` (linear COM zeroed system-wide) + `fix mom_rot interior_mobile momentum 100 angular` (rotation removed from interior only).
- *2026-04-19*: **Exponential velocity gradient, shear_frac 0.30→0.15** — Linear 30% gradient failed at step ~330k (97% of Phase 2): Δv per bond at inner boundary ≈ 0.030×vshear was still too large. Switched to exponential profile v(ζ)=vshear×(exp(k·ζ)−1)/(exp(k)−1), k=3, with shear_frac=0.15. Inner-boundary Δv per bond ≈ 0.0094×vshear (3.2× safer). Outer-surface gradient is steep but all bonds there connect two gradient atoms (no discontinuity). gel_strain_cm will reach ~0.9×target (tracking atoms at ζ≈0.93–1.0 get 0.81–1.0×vshear); G calculation uses measured gel_strain_cm anyway. Angular momentum zeroed for interior_mobile in Phase 2 and 3 (angular keyword confirmed present). Velocity gradients are purely horizontal (x only; NULL for vy/vz).
- *2026-04-19*: **Gradient `fix move variable` — argument order bug found and fixed** — First attempt used wrong argument order: `fix move_top_grad poly_top_grad move variable v_vx_top_grad NULL NULL NULL NULL NULL`. The `fix move variable` arg order is **x y z vx vy vz** (positions first, then velocities). Passing the velocity variable as arg 1 prescribed x-position ≈ `vshear * Δz/grad_thick` ≈ 0.001σ (atoms stayed essentially frozen), not x-velocity. Result: gel_strain_cm ≈ 6e-5 (effectively zero). Also attempted 5-layer staircase with `fix move linear` to avoid the issue, but hit LAMMPS max-32-groups limit (35 groups needed). Fix: correct argument order `NULL NULL NULL v_vx_top_grad NULL NULL` — prescribes vx only, leaves positions free. Current working syntax:
  ```lammps
  variable vx_top_grad atom "v_vshear * (z - v_grad_top_zlo) / v_grad_thick"
  fix move_top_grad poly_top_grad move variable NULL NULL NULL v_vx_top_grad NULL NULL
  ```
- *2026-04-20*: **Gradient zone rotation — diagnosed and abandoned entire approach** — Even with vz=vy=0 added to `fix move variable` (to prevent rotation), the exponential gradient approach failed with FENE errors around step ~120k-200k. The gradient zone rotation was from free vz; fixing it via `v_zero_v` eliminated rotation but FENE failures persisted from a different mechanism. Decision: abandon all polymer-driven gradient approaches and switch to the **plate-driven approach** (see below).
- *2026-04-20*: **Switched to plate-driven shear** — Created `scripts/add_plates_to_gel.py` which adds rigid graphene-like plates (atom type 4, square lattice spacing 1.5σ) to the top/bottom of the isolated gel. Each plate atom is bonded to the nearest surface polymer via a harmonic bond (k=30, r0=1σ, no max extension). `shear_slab.lmp` completely rewritten: Phase 1 NVT with plates frozen; Phase 2 `fix move linear ±vshear 0 0` on plate groups; Phase 3 plates frozen + NVT production. Bond style changed to `hybrid fene harmonic`. This approach eliminates the driven/interior FENE interface entirely: plate atoms (type 4) are driven; all FENE bonds are between polymer atoms only; harmonic bonds transmit shear to surface polymer; distributed gel strain ~0.09% per bond at 10% target strain — far below FENE limit.
- *2026-04-28*: **shear_slab.lmp — plates are x-face yz-planes (not top/bottom)** — `add_plates_to_gel.ipynb` places plates on the x-faces of the gel (yz-planes, normal to x). Shear direction is z (not x). `plate_right` / `plate_left` groups split by xmid (not zmid). `fix move linear 0 0 ±vshear` drives plates in z. Gel shears in xz plane. Old notes referencing `plate_top`/`plate_bot`, z-split, or x-direction shear are superseded.
- *2026-04-28*: **Trajectory dumps now include plates** — Both `traj_setup` and `traj_prod` previously used `gel` group (type 4 plate atoms absent from trajectory). Changed to `all` → plates now visible in OVITO. The "half-plates at top and bottom" artifact was the flat z-faces of the gel slab being mistaken for plates (plates were invisible); the "left half pushed down" was correct xz-shear + PBC wraparound visualisation, not a physics bug.
- *2026-04-28*: **Phase 1 changed to NPT (50k) + NVT (100k)** — Old Phase 1 was NVT-only with frozen plates. New Phase 1a: `fix npt gel iso P_target` (50k steps, plates frozen with `fix move linear 0 0 0`); box x adjusts to reach P=1.5; plate atoms scale with box (their reduced coords stay fixed). Phase 1b: NVT at NPT-final volume (100k steps, plates remain frozen). This ensures the gel enters shear at the correct equilibrium density. Minimize improved: plates frozen with `fix setforce 0 0 0` + `velocity plate set 0 0 0` during minimize (5000/50000 steps) so only gel atoms relax. Previously minimise ran on all atoms including plates (500/5000 steps), leaving bond-stretch strains that caused high initial PE (~12).
- *2026-04-28*: **Bulk stress region** — Atoms within 5σ of each plate face excluded from stress measurements (`poly_bulk` / `solv_bulk`, x ∈ [gel_xlo+5, gel_xhi−5]). These near-wall atoms carry artifactual constraint stresses from the harmonic plate-polymer bonds. All stress computes (`stress/atom`, `reduce sum`, `chunk/atom`) and output fixes (`ave/time`, `ave/chunk`) updated to use bulk groups. Stress normalised by `bulk_vol = bulk_lx × ly × lz` (not `vol`); z-bin stress normalised by `bulk_lx × ly × binWidth`. Using `vol` was incorrect and would understate stress magnitudes.
- *2026-04-28*: **Shear speed +10%, doubled steps, `fix halt` at target strain** — `vshear` increased to `1.1 × target_strain × gel_gap / (2 × nsteps_shear_base × dt)` (base = 285000 steps). `nsteps_shear` doubled to 570000 as a safety buffer. `fix halt_shear all halt ${stress_freq} v_gel_strain_cm >= ${target_strain_xz} error continue` stops Phase 2 as soon as gel_strain_cm ≥ 0.10, then Phase 3 proceeds from the exact target strain. Pattern mirrors `slab_with_flow.lmp` compression halt.
- *2026-04-28*: **`c_gel_press` vs `press`** — `compute gel_press all pressure gel_temp pair bond` reports virial-only pressure (~0.86) because listing `pair bond` keywords excludes the kinetic NkT/V contribution. Full system pressure is `press` in thermo (~1.5). The NPT barostat targets `press`, not `c_gel_press`. Both are physically meaningful (configurational vs total pressure).
- *2026-04-28*: **`shear_analysis.ipynb` created** — New analysis notebook at `scripts/shear_analysis.ipynb` for extracting G, N1, N2, z-profile stresses, and poroelastic decomposition from shear_slab Phase 3 output.
