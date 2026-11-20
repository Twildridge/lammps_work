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

Coarse-grained molecular dynamics of tetrahedral hydrogel slabs (Kremer–Grest bead-spring model, LAMMPS), used to extract poroelastic observables directly from MD trajectories: pore pressure, volume fraction profiles, cooperative diffusivity, permeability, and elastic moduli (longitudinal M, shear G, bulk K).

**Physical analog:** PEG-based hydrogel membranes used in filtration, drug delivery, and soft robotics.

## Recent runs

The jobs most recently started on any cluster. Each job records itself when it starts and when it ends, and `lsync` rebuilds this table ([how it works](docs/run_log.md)). The last column is the folder of that run under `flow_data_local/` on the laptop, filled in once a notebook has pulled its files.

<!-- RECENT_RUNS:BEGIN (written by scripts/lib/run_record.py; edits here are overwritten) -->
| Started | Cluster | Job | Deck | Lattice | Initial data file | Final data file(s) | Simulated | Wall time | Status | Laptop RUN_ID |
|---|---|---|---|---|---|---|---|---|---|---|
| 2026-11-17 03:30 | expanse |  | uniaxial_compression | random, 10 beads/strand | `final_config_slab_support_periodic_10beads_wide_rho04_random_1.0_1.0_16000000_two_pist_pm40.data` | `final_tricomp_final_config_slab_support_periodic_10beads_wide_rho04_random_1.0_1.0_16000000_two_pist_pm40_1.0_1.0_4000000.data` | 13.65 M steps = 68,250 τ | 10 h 01 min | complete |  |
| 2026-11-15 22:44 | expanse |  | uniaxial_compression | random, 10 beads/strand | `final_config_slab_support_periodic_10beads_wide_rho04_random_1.0_1.0_16000000_two_pist_pm40.data` | `final_tricomp_final_config_slab_support_periodic_10beads_wide_rho04_random_1.0_1.0_16000000_two_pist_pm40_1.0_1.0_4000000.data` | 13.65 M steps = 68,250 τ | 9 h 53 min | complete |  |
| 2026-11-15 22:44 | expanse |  | simple_shear | random, 10 beads/strand | `final_config_slab_support_periodic_10beads_wide_rho04_random_1.0_1.0_16000000_two_pist_pm40.data` | `final_shear_final_config_slab_support_periodic_10beads_wide_rho04_random_1.0_1.0_16000000_two_pist_pm40_1.0_1.0_3000000.data` | 4.50 M steps = 22,500 τ | 3 h 19 min | complete |  |
| 2026-11-15 22:44 | expanse |  | uniaxial_permeation | random, 10 beads/strand | `final_config_slab_support_periodic_10beads_wide_rho04_random_1.0_1.0_16000000_two_pist_pm40.data` | `final_triperm_final_config_slab_support_periodic_10beads_wide_rho04_random_1.0_1.0_16000000_two_pist_pm40_1.0_1.0_4000000.data` | 6.21 M steps = 31,050 τ | 4 h 31 min | complete |  |
| 2026-11-15 03:54 | expanse |  | chain_length_testing | random, 10 beads/strand | `slab_support_periodic_10beads_wide_rho04_random.data` | `final_config_slab_support_periodic_10beads_wide_rho04_random_1.0_1.0_16000000.data` | 14.62 M steps = 73,100 τ | 7 h 02 min | complete |  |

The 5 most recent of 21 recorded jobs; all of them: [`runs/README.md`](runs/README.md).
<!-- RECENT_RUNS:END -->

## Start here

```bash
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
bash install.sh        # once per machine, then open a new terminal
```

That is the whole setup on any machine. `install.sh` adds two commands, `lsync` (pull, commit and push in one step) and `lsbatch` (submit a job from any folder on a cluster), and on a cluster it also installs what the jobs need (on Pod: the conda environment for the plotters and, if missing, the LAMMPS binary, built on the login node in about 10 min) and ends with a check of the whole setup. [docs/getting_started.md](docs/getting_started.md) goes from there to a first finished simulation and its plots, and assumes no git experience; [docs/cluster_setup.md](docs/cluster_setup.md) is the cluster page.

## Repository layout

```
lammps_work/
├── install.sh      one-time setup of a clone
├── simulations/    one folder per simulation type: the LAMMPS deck (.lmp) and its batch files
├── scripts/        analysis notebooks, the job runner and the sync / submit scripts
│   └── lib/        all the Python (analysis modules, data-file builders, plotters) and the cluster settings
├── docs/           the guides indexed below
├── lammps_build/   how the Expanse and Pod LAMMPS binaries were built
└── archive/        superseded decks, notebooks and runners, kept for reading old runs
```

Input `.data` files are not in the repo. They live in `lammps_data_files_local/` next to the clone on the author's Mac (catalogued in its `README.md`) and in `~/Documents/lammps_data/input_data/` on a cluster. The full tree, and the table of names changed on 2026-11-18 (`triaxial` → `uniaxial`, `shear_slab` → `simple_shear`, `compress_slab` → `isotropic_compression`), is in [docs/repository_layout.md](docs/repository_layout.md).

## Simulation types

| Folder | Measures | Notebook in `scripts/` |
|---|---|---|
| `slab_with_support/` | Free-swelling equilibration of a gel; the starting point of every other run | builders: `data_file_generator.ipynb` (isolated gel), `slab_with_support_periodic.ipynb` (periodic slab) |
| `uniaxial_compression/` | Longitudinal modulus M, cooperative diffusivity D_c and permeability κ from a drained compression and its unload | `uniaxial_compression_{single,sweep}.ipynb` |
| `uniaxial_permeation/` | Darcy permeability κ, D_c and M from pressure- or velocity-driven flow through the gel | `uniaxial_permeation_single.ipynb` |
| `simple_shear/` | Shear modulus G | `simple_shear_{single,sweep}.ipynb` |
| `isotropic_compression/` | Bulk modulus K and the free-swelling diffusivity D_m | `isotropic_compression_{single,sweep}.ipynb` |
| `pure_solvent/`, `pure_polymer/` | Reference equations of state | builders: `pure_solvent.ipynb`, `pure_polymer.ipynb` |
| `chain_length_testing/` | Equilibration across strand lengths | `chain_length_comparison.ipynb` |
| `calibration_sweep/` | Solvent volume-fraction calibration across pressure and interaction strength | `calibration_analysis.ipynb` |
| `allocation_benchmark/` | Fastest nodes × ranks × threads on a cluster | its `README.md` |

The notebooks listed are analysis notebooks unless marked as builders of the input file. `slab_two_pistons_periodic.ipynb` converts an equilibrated slab into the two-piston input of the uniaxial and shear decks. `cross_deck_comparison.ipynb` combines the compression and permeation runs (D_c and κ against solvent fraction, the equation-of-state and flux closures). What each deck does and its knobs are in [docs/running_simulations.md](docs/running_simulations.md); what each notebook computes is in [docs/analysis.md](docs/analysis.md).

## Clusters

Production runs on Expanse (SDSC). Pod (CNSI, UCSB) is the backup since 2026-11-18 and uses the same runner, decks and batch files. Each batch file carries an allocation for both, and `lsbatch <name>.batch` submits with the one of the cluster it is typed on. The settings of each cluster are in `scripts/lib/site.sh` and are compared side by side in [docs/cluster_reference.md](docs/cluster_reference.md), with the node counts and layouts the allocation benchmark found fastest.

## Where to find things

| Topic | Doc |
|---|---|
| First-time setup and a first simulation | [docs/getting_started.md](docs/getting_started.md) |
| Simulation types, batch-file knobs, submitting, continuing a run | [docs/running_simulations.md](docs/running_simulations.md) |
| Building `.data` input files | [docs/building_data_files.md](docs/building_data_files.md) |
| Analysis notebooks and the plotters that run after each job | [docs/analysis.md](docs/analysis.md) |
| Downloading output and trajectories | [docs/downloading_output.md](docs/downloading_output.md) |
| `lsync`, collaborator access, SSH keys on a cluster, merge conflicts | [docs/github_sync.md](docs/github_sync.md) |
| A new cluster account: clone, Python, copying data files | [docs/cluster_setup.md](docs/cluster_setup.md) |
| Expanse and Pod settings, partitions, allocation sizes | [docs/cluster_reference.md](docs/cluster_reference.md) |
| The run log behind the "Recent runs" table | [docs/run_log.md](docs/run_log.md) |
| LJ units, simulation parameters, barostat and piston choices | [docs/physics_reference.md](docs/physics_reference.md) |
| A specific error message | [docs/common_issues.md](docs/common_issues.md) |
| Every file in the repo, and the layout of a run folder | [docs/repository_layout.md](docs/repository_layout.md) |
| Scripts in `scripts/`, one line each | [scripts/README.md](scripts/README.md) |
| What is in `archive/` and how to read old runs | [archive/README.md](archive/README.md) |

Standalone notes at the repo root: [`lj_units_cheat_sheet.md`](lj_units_cheat_sheet.md), [`expanse_lammps_guide.md`](expanse_lammps_guide.md), [`slurm_commands_and_compiling.md`](slurm_commands_and_compiling.md), [`slab_data_file_info.md`](slab_data_file_info.md), [`simple_shear_notes.md`](simple_shear_notes.md). The plan for the random 10-bead gel is in [docs/random10_gel_plan.md](docs/random10_gel_plan.md).

The history of changes is in the git log. For questions, contact Dylan Pollard (pollard@ucsb.edu).
