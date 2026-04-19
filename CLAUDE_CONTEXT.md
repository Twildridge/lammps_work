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
   - Script: `scripts/longitudinal_modulus_analysis.ipynb`

4. **Calculate G (shear modulus)** ← *CURRENT FOCUS*
   - Method: Shear the simulation box (xz plane) and record gel stresses
   - Verify G decreases with decreasing crosslink density
   - Simulation script: `simulations/shear_slab/shear_slab.lmp` (created 2026-04-15)
   - Input: isolated swollen gel from `isolate_gel.ipynb` → `lammps_data/input_data/`
   - G = σ_xz / γ_xz from Phase 3 time-averaged stress tensor
   - Still need: run shear_slab, write shear_modulus_analysis.ipynb

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
│   │   ├── shear_slab/          ← shear modulus sims (renamed from slab_elongation)
│   │   │   ├── shear_slab.lmp   ← xz shear sim; 3 phases: NPT + shear + NVT production
│   │   │   └── shear_slab.batch ← Bridges-2 batch script
│   │   ├── slab_elongation/     ← OLD (uniaxial elongation); keep or delete manually
│   │   ├── polymer_phase/
│   │   ├── polymer_pure/
│   │   ├── solvent_phase/
│   │   └── solvent_pure/
│   ├── scripts/
│   │   ├── longitudinal_modulus_analysis.ipynb   ← M calculation (done)
│   │   ├── flow_poroelasticity_analysis.ipynb    ← permeation analysis
│   │   ├── add_walls_to_slab.ipynb               ← data file builder
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
│   │   ├── run_lammps.sh / run_lammps_pod.sh / run_lammps_bridges.sh
│   │   └── write_tracking.py
│   ├── lammps_data/             ← small data files (in repo)
│   ├── README.md                ← workflow documentation
│   ├── expanse_lammps_guide.md  ← cluster-specific guide
│   ├── slab_data_file_info.md   ← data file specs log
│   ├── documenting_pod_runs.md  ← HPC performance notes
│   ├── lj_units_cheat_sheet.md
│   ├── git_setup_instructions.md
│   └── tracking.txt             ← legacy performance log (not actively used)
│
├── flow_data_local/             ← local copy of HPC output (NOT in git)
│   ├── partial_stress_data/     ← .dat files: stress_x/y/z_polymer/solvent, strain, piston
│   └── flow_plots/              ← output plots from analysis notebooks
│
├── lammps_data_files_local/     ← large .data files (NOT in git, local only)
│
├── add_walls_to_slab.ipynb      ← (root-level, older versions)
└── slab_with_support_4.ipynb
```

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

## shear_slab.lmp — Design Notes

- **3 phases**: (1) NPT iso equilibration → (2) instantaneous `change_box xz final` + NVT relaxation → (3) NVT production at fixed deformed box → G = ⟨σ_p_xz⟩ / target_strain_xz
- **Shear method** (FINAL, as of 2026-04-18): `change_box all xz final ${target_xz} remap` applies the full target_strain_xz instantaneously (affine remap). NVT then runs at the LOCKED box geometry for `nsteps_shear` steps. No driven/interior interface; no FENE failures possible.
- **Why this works for elastic gels**: Bonds crossing the x-PBC carry an image offset of xz in x — a topological constraint. NVT at fixed box shape cannot relax this interior strain without an NPT barostat on xz. Surface layers (poly_top_s/poly_bot_s) do partially relax; `gel_strain_cm` tracks this. That is EXPECTED and physically correct — not a sign of failure.
- **Strain convention**:
  - `gel_strain_box = xz / lz = target_strain_xz` — box-level, interior network, fixed throughout Phase 2 and 3
  - `gel_strain_cm = (xcm_top_x − xcm_bot_x − xcm_sep_0) / gel_thick` — surface layers, partially relaxes during NVT
  - **G uses box strain**: G = ⟨σ_p_xz⟩ / target_strain_xz (NOT gel_strain_cm)
- **Why NOT `fix deform remap x`** (tried, failed 2026-04-18): NVT propagates FENE/crosslink restoring velocities that oppose each step's affine remap. Elastic gel COM strain saturates at ~45% of target (observed: gel_strain_cm plateaued at 0.045 for γ_target=0.10). Cannot reach target strain with continuous rate-controlled shear of an elastic network.
- **Why NOT `fix move` rheometer** (tried, failed 2026-04-18): `fix move` overrides atom velocity post-integration every step → hard velocity discontinuity at driven/interior boundary → FENE bond straddling that boundary stretches monotonically. Tried: reducing shear_frac 0.1→0.05, slowing vshear 30%, excluding crosslinkers from driven groups. All iterations pushed warning onset later (step 92k→107k→149k) but never eliminated the failure. Root cause is structural — there will always be exactly one bond pair at the interface.
- **Why NOT `fix addforce`**: constant force per atom → jagged surfaces, and elastic equilibrium strain ≈ F/(G·A) is uncontrolled → massive over-shearing (hit LAMMPS triclinic limit |xz| ≤ lx/2).
- **Why NOT `fix halt`**: `fix halt` calls `timer->force_timeout()` internally; this flag persists across `run` commands even after `unfix halt_shear`. Phase 3's `run ${nsteps}` exits at 0 steps.
- **`poly_top_s` / `poly_bot_s` groups**: top/bottom 5% of gel by z-extent, backbone type 1 only. Used only for COM tracking (xcm_top/xcm_bot computes) — NOT driven in the current implementation. Two-step definition required: `group top_all region reg_top` then `group poly_top_s intersect backbone top_all` (LAMMPS forbids intersect with dynamic groups).
- **Pair style**: WCA (`lj/cut 1.122`), consistent with slab_with_flow
- **Bond coeffs**: `30.0 1.5 1.0 1.0` (Kremer-Grest FENE). Max bond extension from change_box at 10% shear ≈ target_strain_xz × Δz_typical ≈ 0.05–0.1σ (well within R0=1.5σ).
- **Parameters to tune**: `target_strain_xz` (default 0.10), `nsteps_shear` (default 285000 = NVT relaxation steps post-change_box)
- **Triclinic setup**: `change_box all triclinic` placed before `dump traj_setup` in Phase 1 (LAMMPS restriction: cannot change box topology while dumps are active).
- **z-binning**: `compute chunk/atom bin/1d` requires `units reduced` for triclinic boxes. `binWidth_reduced = $(v_binWidth/lz)` snapshotted post-Phase 2 when lz is frozen.
- **Output frequencies**: `thermo_freq=10000`, `stress_freq=5000`, `num_stress_curves=10`
- **Thermo columns**: includes `Pxz` and `Xz` — used by `plot_lammps_log.py` to auto-detect shear runs
- **`gel_strain_cm` variable**: `gel_strain_cm = (xcm_top_x − xcm_bot_x − xcm_sep_0) / gel_thick`. Starts at ~target_strain_xz immediately after change_box, then partially relaxes during Phase 2 NVT. Written to `shear_strain_*.dat` (4 columns: step, gel_lz_initial, gel_thick, gel_strain_cm).
- **Stress outputs**:
  - `stress_tensor_polymer/solvent_*.dat`: box-integrated 6-component tensor vs time (key for G); xz is col 5 (1-indexed)
  - `stress_profile_z_polymer/solvent_*.dat`: z-binned 6-component profile
  - `shear_strain_*.dat`: **4 columns**: step, gel_lz_initial, gel_thick, gel_strain_cm (surface COM)
  - `box_dimensions_*.dat`: step lx ly lz xy xz yz (7 columns; xz = target_xz = constant in Phase 3)
  - `gel_dimensions_rg_*.dat`: step lx_rg ly_rg lz_rg
- **G extraction**: G = ⟨stress_p_xz⟩ / target_strain_xz (time-averaged in Phase 3, box-integrated polymer stress divided by box-level shear strain)
- **Solvent stress**: isotropic at equilibrium (σ_s_xz ≈ 0); G comes entirely from polymer network
- **Optimal node count**: 1–2 nodes on Bridges-2/Expanse (~165k atoms → ≥700 atoms/task). Rule: target ≥1000 atoms/MPI task.
- **IMPORTANT**: `slab_elongation/` folder still exists — delete manually after verifying shear_slab works

## Active Tasks / Open Questions

- [ ] Run shear_slab on 1–2 nodes (not 4 — too small a system) and copy output to `flow_data_local/`
- [ ] Create `shear_modulus_analysis.ipynb` once first results are in
- [ ] Verify G decreases with decreasing crosslink density (rho02 → rho08 series)
- [ ] Milestone 2 cleanup: edit LAMMPS script so P* = 1.5 from t=0; align traj/volume output frequencies
- [ ] Keep README.md and expanse_lammps_guide.md up to date

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
- *2026-04-17*: **`group intersect` with dynamic groups** — LAMMPS forbids `group intersect` when either operand is a dynamic group. `polymer_top` / `polymer_bot` (defined by `dynamic/chunk`) are dynamic → error "Cannot intersect groups using a dynamic group". Fix: capture static region snapshots first (`group poly_top_all region reg_driven_top`), then intersect two static groups (`group poly_top_s intersect polymer poly_top_all`). The `polymer` group is always static (defined by atom type).
- *2026-04-17*: **Nested `${}` inside `$()` is invalid LAMMPS syntax** — `$(${target_strain_xz}*lz)` causes "Invalid syntax in variable formula". Fix: use `v_varname` reference inside the `$()` expression instead: `$(v_target_strain_xz*lz)`. Rule: inside `$()`, reference variables as `v_name`, never as `${name}`.
- *2026-04-17*: **`fix print append` is not a valid keyword** — LAMMPS `fix print` has no `append` option; the keyword causes "Illegal fix print command" error. To append across phases without re-opening the file, keep the fix running continuously (do not unfix and redefine mid-script). Alternatively, accept that Phase 2 and Phase 3 write separate files.
- *2026-04-17*: **`shear_strain_*.dat` now has 4 columns** (updated from 3): step, gel_lz_initial, xz_tilt (box tilt from `$(xz)`), gel_strain (computed from top/bot CM separation via `v_gel_strain`). Written during Phase 2 only (fix unfixed at end of Phase 2; Phase 3 does not re-define it). γ = col4 directly, or col3/col2 as cross-check.
- *2026-04-17*: **`fix move` + `fix deform remap none` synchronization** — set `vshear = erate × lz` (evaluated once at Phase 2 start with `variable vshear equal $(v_erate*lz)`). Both top-surface CM displacement and box xz increase at `erate × lz` per unit time → they are exactly synchronized. After `nsteps_shear` steps both equal `target_strain_xz × lz`. No `change_box` correction needed.
- *2026-04-17*: **`fix move` rheometer → FENE bond failure at driven/interior interface** — Observed in log: bond atoms 18841–45089, length growing 1.423036 → 1.4231 → ... monotonically at ≈ vshear/step (0.01138 σ/τ ÷ 0.005 = 0.00569 σ/step × 2 ≈ measured rate). Root cause: `fix move` overrides velocity AFTER force integration; the FENE restoring force is computed correctly but the driven atom's velocity is reset to vshear regardless → bond stretches without bound. Fix: **switch to `fix deform remap x`** (affine remapping of ALL atoms). Adjacent bonded atoms (Δz ≈ 0.1–0.5σ) accumulate at most target_strain_xz × 0.5σ = 0.05σ of relative stretch total — nowhere near FENE R0 = 1.5σ. Driven/interior interface eliminated entirely.
- *2026-04-17*: **`remap x` + `fix momentum linear`** — With `remap x`, COM drifts in x at erate × z_COM × dt/step; `fix momentum linear 1 1 1` (every 100 steps) removes this. `fix momentum angular` is NOT used in Phase 2 because a shear flow has net angular momentum (L ∝ erate × M × Rz²) — zeroing it would remove the shear gradient itself. In Phase 3, `fix mom_zero all_mobile momentum 100 linear 1 1 1 angular` removes both COM drift and residual angular momentum from Phase 2 without contaminating σ_xz (virial is computed from forces, not velocities).
- *2026-04-18*: **`gel_strain` corrected to `xz/lz`** — Previous formula `xz/v_gel_lz_initial` (where gel_lz_initial is the Rg-based gel thickness ≈ gel_thick) overstates γ by factor lz/gel_thick ≈ 1.018 because lz = gel_thick + 2σ buffer. With `remap x`, xz = erate × lz × t exactly, so xz/lz = erate × t = target_strain_xz. The Rg-based gel_lz_initial is still stored as a reference but is no longer used in the strain formula.
- *2026-04-18*: **`fix deform remap x` fails for elastic gels** — gel COM strain saturates at ~45% of target (observed: 0.045 for γ_target=0.10). Root cause: FENE/crosslink restoring forces generate velocities opposing each step's remap displacement; NVT propagates these velocities. Elastic solid resists continuous rate-controlled shear. Fix: replace Phase 2 with `change_box all xz final ${target_xz} remap units box` (one instantaneous affine step) followed by `run 0` + NVT at LOCKED box geometry. The x-periodic bond topology (bonds crossing x-PBC carry image offset of xz in x) enforces gel_strain = target_strain_xz — NVT cannot un-shear the gel without an NPT barostat. FENE safety: max per-bond extension at 10% shear ≈ 0.097σ (6% of R0=1.5σ — safe). `nsteps_shear` now = NVT relaxation steps at fixed geometry (was: ramping steps with fix deform).
- *2026-04-18*: **Added polymer COM gel strain tracking** — `poly_top_s` / `poly_bot_s` groups (top/bottom 10% of gel by z, static) + `compute xcm_top/xcm_bot com` + `xcm_sep_0` snapshot → `gel_strain_cm = (xcm_top_x − xcm_bot_x − xcm_sep_0) / gel_thick`. Now 4th column of `shear_strain_*.dat`. With change_box approach, gel_strain_cm starts at ~target_strain_xz and partially relaxes (expected for surface layers); interior network held by periodic topology.
- *2026-04-18*: **`fix move` rheometer approach ultimately abandoned** — Three rounds of incremental fixes (shear_frac 0.1→0.05, vshear 30% slower, backbone-only driven groups) pushed FENE warning onset from step 92k→107k→149k but never eliminated the failure. Root cause: `fix move` overrides velocity post-integration every step, creating a hard velocity discontinuity at the driven/interior boundary. The one bond straddling that boundary always stretches monotonically. **Final decision: revert to `change_box` + NVT approach** (as described in the current Design Notes above). Key insight: gel_strain_cm tracks surface relaxation (expected and physically correct); interior gel is held at target_strain_xz by x-periodic bond topology. G = ⟨σ_p_xz⟩ / target_strain_xz using box strain, not gel_strain_cm.
