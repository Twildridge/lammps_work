# scripts/ — what each file is (2026-10-11)

Rule (same as `lammps_data_files_local/`): a notebook or script lives here while a deck in `simulations/`, a doc, or a current analysis uses it. Superseded ones are in `archive/scripts/` (open them from this folder: their relative paths assume `scripts/`).

## Shell scripts (run on the clusters unless noted)

| Script | Role |
|---|---|
| `run_lammps.sh` | **The Expanse job runner.** Every `.batch` in `simulations/` calls it: finds the deck by folder name, forwards the env knobs as `-var`s, runs `mpirun lmp`, then `postprocess.sh`; with `CHAIN=1` hands over to `continue_sim.sh --chain`. |
| `postprocess.sh` | Regenerates every post-run plot of a finished run directory (calls `lib/plot_*.py`); run automatically by the runners, or by hand. |
| `continue_sim.sh` | Resume a finished run, or chain a run that hit its wall-clock budget, from the simulation folder (`docs/running_simulations.md` §5e). |
| `git_sync.sh` | What `lsync` runs: pull --rebase, stage, commit, push (`docs/github_sync.md`). |
| `install_lsync.sh` | One-time per machine: registers the `lsync` alias and nbstripout. |
| `node_git_sync.sh` | Sourced by `run_lammps.sh`: best-effort `git pull` on the compute node before a job starts. |

## Data-file builders and converters (notebooks; converter code in `lib/`)

| Notebook | Role |
|---|---|
| `slab_with_support.ipynb` | Finite gel block in a bath, with or without sheets (`includeWalls`). Current: the tall isolated gel (`tall_rho04_nowalls`). |
| `slab_with_support_periodic.ipynb` | The x/y-periodic slab — the reference 14000002 state comes from it. |
| `slab_two_pistons.ipynb` | Equilibrated periodic slab → two-piston (feed / permeate) file for `triaxial_*_two_pist` and `shear_slab_two_pist`. |
| `add_more_plates_to_gel.ipynb` | Equilibrated isolated gel → six plates + two wet pistons for `compress_slab`. |
| `pure_solvent_1.ipynb`, `pure_polymer.ipynb` | Pure-phase boxes for `solvent_phase` / `polymer_phase`. |

## Analysis notebooks (code in `lib/`)

| Notebook | Role |
|---|---|
| `triaxial_compression_single_two_pist.ipynb`, `triaxial_compression_sweep_two_pist.ipynb` | Two-piston compression: one level / whole sweep (M, G, D_c, κ). |
| `triaxial_permeation_single_two_pist.ipynb` | Two-piston permeation (Q, κ, D_c, M from the polymer displacement). |
| `shear_analysis_single.ipynb`, `shear_analysis_sweep.ipynb` | Shear G for `shear_slab` and `shear_slab_two_pist`. |
| `bulk_modulus_analysis_single.ipynb`, `bulk_modulus_analysis_sweep.ipynb` | Bulk modulus K from `compress_slab`. |
| `cross_deck_comparison.ipynb` | Figures that combine decks (D_c(φ_s), κ(φ_s), EOS). |
| `calibration_analysis.ipynb` | φ_f calibration / PMV sweep (`calibration/` holds its λ table and summary figure). |

## Folders

`lib/` all Python (analysis modules, converters, plotters — see `docs/repository_layout.md`); `tests/` deck linter, fixtures, plot and notebook smoke tests; `calibration/` the calibration outputs the notebooks read.

Moved out on 2026-10-11: the Pod and Bridges-2 runners (with every `*_pod.batch` / `*_bridges.batch` and the Pod benchmarks) → `../archive/pod_bridges/`; the LAMMPS build script and its log → `../lammps_build/`; the one-piston analysis notebooks (`triaxial_compression{,_single,_sweep}`, `triaxial_permeation`), the walled / angled / isolate / split / add_plates builders and `run_analysis_notebook.sh` → `../archive/scripts/`.
