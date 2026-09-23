# LAMMPS Hydrogel Simulation

Coarse-grained molecular dynamics simulations of tetrahedral hydrogel slabs (Kremer–Grest bead-spring model, LAMMPS), used to extract poroelastic observables directly from MD trajectories: pore pressure, volume fraction profiles, cooperative diffusivity, and elastic moduli (longitudinal M, shear G, bulk K).

**Physical analog:** PEG-based hydrogel membranes used in filtration, drug delivery, and soft robotics.

---

## New to this repo? Start here

**→ [docs/getting_started.md](docs/getting_started.md)** — a complete, copy-pasteable walkthrough from "I've never used this repo" to "I ran a simulation and looked at its output." No git experience assumed.

Everything else in this README is a short index into deeper reference material for once you're up and running.

---

## Repository layout (short version)

```
lammps_work/
├── simulations/    ← One folder per simulation type — see running_simulations.md
├── scripts/        ← Notebooks (build input files) + shell/Python scripts (run jobs, analyze output)
├── docs/           ← Everything below this point in more depth (start with getting_started.md)
├── lammps_data/    ← Reserved for small committed .data files (usually empty — see building_data_files.md)
├── archive/        ← Superseded, unmaintained: one-piston triaxial decks (2026-09-22), volmix/Widom tools (2026-08) — see archive/README.md
└── *.md            ← A few standalone reference docs (unit cheat sheet, cluster guides — indexed below)
```

Full annotated tree, including every notebook/script and what it's for: [docs/repository_layout.md](docs/repository_layout.md).

---

## Where to find things

| Topic | Doc |
|---|---|
| **First-time setup + your first simulation** | [docs/getting_started.md](docs/getting_started.md) |
| What each simulation type measures, editing/submitting batch jobs, `continue_sim.sh` | [docs/running_simulations.md](docs/running_simulations.md) |
| Building `.data` input files (which notebook, for which geometry) | [docs/building_data_files.md](docs/building_data_files.md) |
| Downloading output & trajectories from a cluster | [docs/downloading_output.md](docs/downloading_output.md) |
| Running analysis scripts and notebooks on output | [docs/analysis.md](docs/analysis.md) |
| **Triaxial compression / permeation — the standard decks are the two-piston ones** (`triaxial_compression_two_pist`, `triaxial_permeation_two_pist`; standard since 2026-09-22): converter `slab_two_pistons.ipynb` → decks → `triaxial_compression_{single,sweep}_two_pist.ipynb` (M, G, D_c, κ vs strain) / `triaxial_permeation_single_two_pist.ipynb`; code in `scripts/lib/triaxial.py`; physics note | [building_data_files.md](docs/building_data_files.md) · [running_simulations.md §5b′](docs/running_simulations.md) · [analysis.md §7b](docs/analysis.md#7b-jupyter-analysis-notebooks) · [physics_reference.md](docs/physics_reference.md#npt-piston-reservoirs-two-piston-sequence) |
| Shear modulus G: `slab_shear_plates.py` (periodic slab → plates on its z-faces, 2026-09-22) → `shear_slab` deck (plates driven ±x, z = gap; laid out like the triaxial decks) → `shear_analysis.ipynb`; compare with G = (σ′_zz − σ′_xx)/2ε from the two-piston compression sweep | [running_simulations.md §5b](docs/running_simulations.md#5b-simulation-types) · [analysis.md §7b](docs/analysis.md#7b-jupyter-analysis-notebooks) |
| One-piston triaxial decks (`triaxial_compression`, `triaxial_permeation`) — **archived 2026-09-22** to `archive/simulations/`; their existing runs are still read by `triaxial_compression_{single,sweep}.ipynb` / `triaxial_permeation.ipynb` | [archive/README.md](archive/README.md) · [docs/analysis.md §7b](docs/analysis.md#7b-jupyter-analysis-notebooks) |
| `lsync` / GitHub sync, credentials, SSH keys, merge conflicts | [docs/github_sync.md](docs/github_sync.md) |
| First-time cluster account setup (SSH, Python env, `.data` transfer) | [docs/cluster_setup.md](docs/cluster_setup.md) |
| Cluster login nodes, partitions, walltimes, scratch paths | [docs/cluster_reference.md](docs/cluster_reference.md) |
| LJ → physical units, current simulation parameters, barostat choices | [docs/physics_reference.md](docs/physics_reference.md) |
| Fixing a specific error message | [docs/common_issues.md](docs/common_issues.md) |
| Full annotated repo tree | [docs/repository_layout.md](docs/repository_layout.md) |

**Standalone reference docs** (repo root, unchanged):
[`lj_units_cheat_sheet.md`](lj_units_cheat_sheet.md) · [`expanse_lammps_guide.md`](expanse_lammps_guide.md) · [`slurm_commands_and_compiling.md`](slurm_commands_and_compiling.md) · [`slab_data_file_info.md`](slab_data_file_info.md) · [`documenting_pod_runs.md`](documenting_pod_runs.md)

---

*Last updated: 2026-09-22. For questions, contact Dylan Pollard (pollard@ucsb.edu).*
