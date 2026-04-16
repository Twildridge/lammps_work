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

- **3 phases**: (1) NPT iso equilibration → (2) xz shear via `fix deform xz erate` + `±x addforce` on top/bottom 5% polymer → (3) NVT production at fixed deformed box
- **Pair style**: WCA (`lj/cut 1.122`), consistent with slab_with_flow (NOT full LJ as in old slab_elongation)
- **Bond coeffs**: `30.0 1.5 1.0 1.0` (Kremer-Grest FENE, matches slab_with_flow)
- **Parameters to tune**: `target_strain_xz` (default 0.05), `shear_force` (default 0.05), `nsteps_shear` (default 200000)
- **Pressure output**: all 6 components in `thermo_style` (`pxx pyy pzz pxy pxz pyz`) + `c_mobile_press` for direct comparison with partial stresses
- **Stress outputs**:
  - `stress_tensor_polymer/solvent_*.dat`: box-integrated 6-component tensor vs time (key for G)
  - `stress_profile_z_polymer/solvent_*.dat`: z-binned 6-component profile via `fix ave/chunk norm none`
  - `shear_strain_*.dat`: xz tilt + γ vs step (covers Phases 2 and 3)
  - `box_dimensions_*.dat`: lx ly lz xytilt xztilt yztilt vs step
- **G extraction**: G = `stress_p_xz` / `shear_strain_current` (time-averaged in Phase 3)
- **IMPORTANT**: `slab_elongation/` folder still exists — delete manually after verifying shear_slab works

## Active Tasks / Open Questions

- [ ] Prepare isolated gel data file for shear_slab: run `isolate_gel.ipynb` on an equilibrated slab, copy .data to cluster `lammps_data/input_data/`
- [ ] Update DATANAME in `shear_slab.batch` to match actual isolated gel file
- [ ] Run shear_slab simulation and copy output to `flow_data_local/`
- [ ] Create `shear_modulus_analysis.ipynb` once first results are in
- [ ] Verify G decreases with decreasing crosslink density (rho02 → rho08 series)
- [ ] Milestone 2 cleanup: edit LAMMPS script so P* = 1.5 from t=0; align traj/volume output frequencies
- [ ] Keep README.md and expanse_lammps_guide.md up to date
- [ ] Explore GitHub automation (commit/push helpers) across local machine + clusters

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
- *2026-04-15*: Created `scripts/git_sync.sh` — one-command pull/commit/push with auto-generated commit messages. Add `lsync` alias to `~/.zshrc` on MacBook (and `~/.bashrc` on each cluster).
- *2026-04-15*: Pod VPN (Ivanti Secure Access) troubleshooting — broken install required clearing remnants (`sudo rm -rf` of app + support files) then **rebooting before reinstalling**. Reboot fixed it. Download from https://it.ucsb.edu/ivanti-secure-access-campus-vpn/get-connected-campus-vpn. If it won't open after reinstall: reboot first.
