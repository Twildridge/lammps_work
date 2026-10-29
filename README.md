> ## ⚠️ AI-generated documentation — please read first
>
> **Every document in this repository (this README, everything under `docs/`, the `*.md` notes, the method
> notes at the end of the notebooks, and the comments in the scripts and LAMMPS decks) is 100 % AI-generated
> (Claude, Anthropic), written in the course of developing and analysing the simulations.** It is working
> documentation of the code and its state, not a scientific publication, and it has not been independently
> reviewed.
>
> **No manuscript, poster, talk, thesis chapter, or any other published result from this project will contain
> AI-guided text.** The simulations, the analysis choices and every published statement are the author's own;
> the AI-generated material here documents the tooling only.

# LAMMPS Hydrogel Simulation

The build script for the Expanse binary these jobs use (`22Jul2025_update3`, MPI + OpenMP, CMake), the SLURM log of that build and the package list are in [`lammps_build/`](lammps_build/README.md).

Coarse-grained molecular dynamics simulations of tetrahedral hydrogel slabs (Kremer–Grest bead-spring model, LAMMPS), used to extract poroelastic observables directly from MD trajectories: pore pressure, volume fraction profiles, cooperative diffusivity, and elastic moduli (longitudinal M, shear G, bulk K).

**Physical analog:** PEG-based hydrogel membranes used in filtration, drug delivery, and soft robotics.

---

## Start here

[docs/getting_started.md](docs/getting_started.md) takes a new user from cloning the repo to a first finished simulation and its plots. It assumes no git experience.

The rest of this README is an index of the reference docs.

---

## Repository layout (short version)

```
lammps_work/
├── simulations/    ← One folder per simulation type — see running_simulations.md
├── scripts/        ← CURRENT notebooks (build input files, analyze output) + the job-runner / sync / post-process shell scripts — see scripts/README.md
│   └── lib/        ← ALL the Python: analysis modules (triaxial/shear/bulk/volfrac/psd), converters, CLI plotters (since 2026-10-01)
├── lammps_build/   ← How the Expanse LAMMPS binary was built (sbatch script + its log)
├── docs/           ← Everything below this point in more depth (start with getting_started.md)
├── archive/        ← Superseded, unmaintained: one-piston triaxial decks (2026-09-22), volmix/Widom tools (2026-08), one-piston notebooks + old builders, Pod / Bridges-2 runners + batches (2026-10-11) — see archive/README.md
└── *.md            ← A few standalone reference docs (unit cheat sheet, cluster guides — indexed below)
```

The full tree, one line per notebook and script, is in [docs/repository_layout.md](docs/repository_layout.md). Input `.data` files are not in the repo: they live in `lammps_data_files_local/` next to it on the Mac (catalogued in its `README.md`) and in `~/Documents/lammps_data/input_data/` on the cluster.

---

## Where to find things

| Topic | Doc |
|---|---|
| **First-time setup + your first simulation** | [docs/getting_started.md](docs/getting_started.md) |
| What each simulation type measures, editing/submitting batch jobs, `continue_sim.sh` | [docs/running_simulations.md](docs/running_simulations.md) |
| Building `.data` input files (which notebook, for which geometry) | [docs/building_data_files.md](docs/building_data_files.md) |
| Downloading output & trajectories from a cluster | [docs/downloading_output.md](docs/downloading_output.md) |
| Running analysis scripts and notebooks on output | [docs/analysis.md](docs/analysis.md) |
| **Triaxial compression / permeation — the standard decks are the two-piston ones** (`triaxial_compression_two_pist`, `triaxial_permeation_two_pist`; standard since 2026-09-22): converter `slab_two_pistons_periodic.ipynb` → decks → `triaxial_compression_{single,sweep}_two_pist.ipynb` (M, G, D_c, κ vs strain) / `triaxial_permeation_single_two_pist.ipynb` (Q_perm, Darcy κ, and since 2026-09-28 D_c and M from the polymer displacement with a q(t) flux check); code in `scripts/lib/triaxial.py`; physics note (the D_c fits of both modes were re-derived 2026-09-28 — held slab: τ₁ = L²/(4π²D_c); prescribed ΔP: τ₁ = L²/(π²D_c)) | [building_data_files.md](docs/building_data_files.md) · [running_simulations.md §5b′](docs/running_simulations.md) · [analysis.md §7b](docs/analysis.md#7b-jupyter-analysis-notebooks) · [physics_reference.md](docs/physics_reference.md#npt-piston-reservoirs-two-piston-sequence) |
| **Cross-deck comparison** (2026-10-03): `cross_deck_comparison.ipynb` + `scripts/lib/compare.py` — figures built from more than one deck's run (compression sweep + permeation: D_c(φ_s), κ(φ_s), the equation-of-state and flux closures). Opening figures (2026-10-08, permeation deck): the normalised thermodynamic pressure with the two `ΔP_th = (4/3)(G/M) ΔP_ext` guides (`G_shear` = 0.2 from the shear decks, `G_comp` = 0.05 from the compression holds), and the **three solvent measures** of the steady membrane — concentration `c_s = φ_s^cal ρ_s,0/M_s` (λ-calibrated Voronoi; `M_s` = 1, so `c_s = φ_s^cal ρ_s,0` in beads/σ³), the mass (= number) fraction `ρ_s/ρ_s,0`, and the probe porosity `ε_g` (`lib/psd.py`, PoreBlazer-style as in Marioni et al.). **Why `r_probe` = 0.45 σ there** (the permeation notebook keeps 0.5): `psd.py` measures the grid distance from the bead *surface*, `d = r_nn − σ_bead/2`, and calls a point void when `d ≥ r_probe`, so the probe's centre must clear a polymer bead centre by `0.5 + r_probe` σ. With 0.5 that is 1.0 σ, the LJ zero-crossing; a thermal solvent bead at kT = ε actually reaches `U_LJ(r) = kT` at `r = 0.97 σ`, so 0.45 (clearance 0.95 σ) is the solvent-bead-sized probe, and 0.5 slightly over-excludes. The choice shifts `ε_g` by a few hundredths, not its ~0.16 gap below `φ_s^cal` (the excluded shell around every bead), which is the point of the figure | [analysis.md §7b](docs/analysis.md#7b-jupyter-analysis-notebooks) |
| Shear modulus G, plate-driven periodic slab (**archived 2026-10-28** → `archive/simulations/shear_slab/`): `scripts/lib/slab_shear_plates_periodic.py` (periodic slab → plates on its z-faces, 2026-09-22) → `shear_slab` deck (plates driven ±x, z = gap) → `shear_analysis_{single,sweep}.ipynb` with `DECK = "shear_slab"` (its runs still load); superseded by the two-piston shear below | [running_simulations.md §5b](docs/running_simulations.md#5b-simulation-types) · [analysis.md §7b](docs/analysis.md#7b-jupyter-analysis-notebooks) |
| Shear modulus G **on the two-piston system** (the current shear workflow): `simulations/shear_slab_two_pist/` — one deck since 2026-10-28 (the laterally periodic deck and its batch are in `archive/simulations/shear_slab_two_pist/`). **Isolated gel** (`shear_slab_two_pist.lmp`; `shear_slab_two_pist_single.batch` = one level, since 2026-10-28 the random-network gel at γ = 0.10; `shear_slab_two_pist_sweep.batch` = 0.01 → 0.03 → 0.05 → 0.1): a gel equilibrated free in its bath (`slab_with_support_single.batch`) → `lib/add_more_plates_to_gel.py` with `layout = "x_plates"` (`..._x_plates.data`: two grip plates parked outside the x-faces) → the deck seats the plates onto the gel's left/right x-faces; grip = plate + the polymer within 3 σ of it, rigid in all directions, driven −z/+z (x = gap, z = shear), stresses in a core away from the grips and the free faces, bath at P* by the two wet pistons → `shear_analysis_{single,sweep}.ipynb` with `DECK = "shear_slab_two_pist_iso"` (the lib relabels x ↔ z at read time). Cumulative strain sweep (drive until the target γ, then hold; `_g<strain>` files), G divided by the network strain (`G_STRAIN`: steady displacement-profile slope), as M does in compression. **Unload (Phase 2c, 2026-10-28, `UNLOAD=1` in the sweep batch only; the single batch is load-only):** after every hold the grips go back to γ = 0 and the gel relaxes (`_u<strain>` files): residual stresses, and the transverse relaxation **D_c is read there** (`fit_Dc_unload`: lag of the cumulative displacement profile to its late plateau, diffusive vs one-time model side by side) — nothing is read from the hold's own transient any more (the hold-onset and fine-cadence displacement files were dropped from the deck) | [running_simulations.md §5b](docs/running_simulations.md#5b-simulation-types) · [analysis.md §7b](docs/analysis.md#7b-jupyter-analysis-notebooks) |
| **Bulk modulus K and the free-swelling diffusivity D_m** (`compress_slab`, rebuilt 2026-09-29; unload 2026-10-28): an isolated gel in a bath held at P* by two wet pistons, six solvent-transparent plates seated on its faces and closed to a volumetric strain (`compress_slab_single.batch` = one level 0.10, `compress_slab_sweep.batch` = 0.05 → 0.16); the hold gives **K** (network ⟨Π⟩ and plate estimators, increments from the ε_vol = 0 reference, divided by the network's own displacement-profile strain), then the plates open and the cube **re-swells freely** (Phase 2c, `_u<level>` files): the free cube is the 3-D unconstrained problem, whose diffusivity is **D_m = (K/M) D_c** — fitted with the exact free-cube series to the Rg edges, the Rg volume strain, the interior Π and the displacement profiles (`lib/bulk.py fit_Dm_unload`), **κ = D_m/K** from the same run; no diffusivity is read from the compression hold. Input: `data_file_generator.ipynb` → `slab_with_support_single.batch` → `add_more_plates_to_gel.py --layout six_plates`. Notebooks `bulk_modulus_analysis_{single,sweep}.ipynb` | [running_simulations.md §5b](docs/running_simulations.md#5b-simulation-types) · [analysis.md §7b](docs/analysis.md#7b-jupyter-analysis-notebooks) |
| **Shear-anisotropy test** (2026-10-25): the shear decks give G ≈ 0.15–0.18 and every compression-side measurement 0.04–0.06 because the gel is a diamond lattice with its cube axes along the box, i.e. cubic, not isotropic — an xz shear of the box is a cube-face shear (C₄₄) while a z compression's deviatoric response is (C₁₁ − C₁₂)/2. Three gels of the same strands decide it: the aligned tall gel, the **45°-rotated** one (`data_file_generator.ipynb`, `NETWORK = "rot45"`: cube axes along the xz box diagonals, so the same shear measures (C₁₁ − C₁₂)/2; `slab_support_5beads_tall_rho04_nowalls_rot45.data`, `slab_with_support_single.batch` with that DATANAME) and the **random f = 4 network** (`NETWORK = "random"`, `bond_switch = 3`: Wooten–Winer–Weaire bond switching, functionality and strand count kept, no lattice directions; `lib/data_file_generator.py`'s `verify_network_isotropy` check reads the file back against random directions, interior beyond 7 σ isotropic, Born Zener ratio 0.985 ± 0.013; `..._nowalls_random.data`, `slab_with_support_single.batch`). Each is equilibrated, converted with `add_more_plates_to_gel.py --layout x_plates` and sheared at one level with **`shear_slab_two_pist_single.batch`** (γ = 0.10, plain load-only hold, no unload; the aligned gel ran as shear_4 / shear_5_sweep, the rot45 gel as shear_6 at γ = 0.05; the random gel's x_plates file is the batch's DATANAME since 2026-10-28 and is the next run, shear_7). D_c in shear is parked: the fast-drive + fine-file route and the unload route are wired and switched off, see [running_simulations.md §5b](docs/running_simulations.md#5b-simulation-types). Rotated G ≈ 0.05 and random G between → lattice anisotropy; rotated G ≈ 0.17 → doubt the compression lateral stress | [building_data_files.md](docs/building_data_files.md) · [running_simulations.md §5b](docs/running_simulations.md#5b-simulation-types) |
| One-piston triaxial decks (`triaxial_compression`, `triaxial_permeation`) — **archived 2026-09-22** to `archive/simulations/`; their existing runs are still read by `archive/scripts/triaxial_compression_{single,sweep}.ipynb` / `archive/scripts/triaxial_permeation.ipynb` (archived 2026-10-11; open them from `scripts/`) | [archive/README.md](archive/README.md) · [docs/analysis.md §7b](docs/analysis.md#7b-jupyter-analysis-notebooks) |
| `lsync` / GitHub sync, credentials, SSH keys, merge conflicts | [docs/github_sync.md](docs/github_sync.md) |
| First-time cluster account setup (SSH, Python env, `.data` transfer) | [docs/cluster_setup.md](docs/cluster_setup.md) |
| Cluster login nodes, partitions, walltimes, scratch paths | [docs/cluster_reference.md](docs/cluster_reference.md) |
| LJ → physical units, current simulation parameters, barostat choices | [docs/physics_reference.md](docs/physics_reference.md) |
| Fixing a specific error message | [docs/common_issues.md](docs/common_issues.md) |
| Full annotated repo tree | [docs/repository_layout.md](docs/repository_layout.md) |

**Standalone reference docs** (repo root, unchanged):
[`lj_units_cheat_sheet.md`](lj_units_cheat_sheet.md) · [`expanse_lammps_guide.md`](expanse_lammps_guide.md) · [`slurm_commands_and_compiling.md`](slurm_commands_and_compiling.md) · [`slab_data_file_info.md`](slab_data_file_info.md) · [`archive/pod_bridges/documenting_pod_runs.md`](archive/pod_bridges/documenting_pod_runs.md)

---

*Last updated: 2026-10-28 (two batch files per simulation folder, `_single` + `_sweep`, the rest in `archive/simulations/`; `shear_slab` archived; D_c read from the unload in `compress_slab` (D_m = (K/M) D_c, free re-swelling) and in `shear_slab_two_pist` (Phase 2c); generator code → `lib/data_file_generator.py` with the `NETWORK` knob, notebook renamed `data_file_generator.ipynb`; converters write no sidecar / png; 2026-10-25: shear-anisotropy test row: rot45 and random-network gels, `bond_switch` + isotropy check in `slab_with_support.ipynb`, the fine-cadence and single-level shear batches; 2026-10-11: Pod and Bridges-2 runners, batches and benchmarks → `archive/pod_bridges/`; scripts/ reorganised: one-piston notebooks and old builders → `archive/scripts/`, LAMMPS build → `lammps_build/`, `scripts/README.md` added; 2026-10-08: cross-deck comparison row: opening figures + the `r_probe` = 0.45 σ justification; 2026-10-01: every `.py` moved from `scripts/` into `scripts/lib/`). For questions, contact Dylan Pollard (pollard@ucsb.edu).*
