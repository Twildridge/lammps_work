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
| Triaxial-compression analysis: `triaxial_compression_single.ipynb` (one strain level) / `triaxial_compression_sweep.ipynb` (M, G, D_c, κ vs strain), code in `scripts/lib/triaxial.py` | [docs/analysis.md §7b](docs/analysis.md#7b-jupyter-analysis-notebooks) |
| `lsync` / GitHub sync, credentials, SSH keys, merge conflicts | [docs/github_sync.md](docs/github_sync.md) |
| First-time cluster account setup (SSH, Python env, `.data` transfer) | [docs/cluster_setup.md](docs/cluster_setup.md) |
| Cluster login nodes, partitions, walltimes, scratch paths | [docs/cluster_reference.md](docs/cluster_reference.md) |
| LJ → physical units, current simulation parameters, barostat choices | [docs/physics_reference.md](docs/physics_reference.md) |
| Fixing a specific error message | [docs/common_issues.md](docs/common_issues.md) |
| Full annotated repo tree | [docs/repository_layout.md](docs/repository_layout.md) |

**Standalone reference docs** (repo root, unchanged):
[`lj_units_cheat_sheet.md`](lj_units_cheat_sheet.md) · [`expanse_lammps_guide.md`](expanse_lammps_guide.md) · [`slurm_commands_and_compiling.md`](slurm_commands_and_compiling.md) · [`slab_data_file_info.md`](slab_data_file_info.md) · [`documenting_pod_runs.md`](documenting_pod_runs.md)

---

*Last updated: 2026-09-02. For questions, contact Dylan Pollard (pollard@ucsb.edu).*
