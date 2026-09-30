# archive/ — superseded / parked tools and decks

Files here are **kept for reference, not maintained**. Nothing in `archive/` is
run by `run_lammps.sh`, `postprocess.sh` or `continue_sim.sh` (they key on
`simulations/<folder>/<folder>.lmp`), and the docs describe only what lives
under `simulations/` and `scripts/`.

## Archived simulation decks (2026-09-22): the one-piston triaxial workflow

| Folder | Came from | Why it is here |
|---|---|---|
| `simulations/triaxial_compression/` | `simulations/triaxial_compression/` | One-piston axial compression: periodic slab, single velocity-controlled piston + support, aniso-NPH Phase 0, cumulative strain sweep with auto-sized holds. Superseded by **`simulations/triaxial_compression_two_pist/`**, which keeps the strain sweep, hold sizing, observables and file conventions but holds both reservoirs at the bath pressure with NPT-pistons (drained consolidation at constant P). Contains the Expanse/Bridges/Pod batch files and `chain_after_slab.batch`. |
| `simulations/triaxial_permeation/` | `simulations/triaxial_permeation/` | One-piston permeation: force-controlled piston pushing solvent through the slab at constant volume. Superseded by **`simulations/triaxial_permeation_two_pist/`** (feed/permeate NPT-pistons, primary flux from the permeate-piston displacement). |

* The runs these decks produced are still analysed by the one-piston notebooks
  (`scripts/triaxial_compression_{single,sweep}.ipynb`, `scripts/triaxial_permeation.ipynb`,
  `scripts/triaxial_compression.ipynb`), which stay in `scripts/` — `lib/triaxial.py`
  handles one- and two-piston runs alike.
* `continue_sim.sh` still knows these folders (`final_tricomp` / `final_triperm`
  prefixes, `_c<level>` sweep tag); it works from a folder that holds the original
  SLURM output file next to the archived `.lmp` / `.batch`.
* **To run one again:** `git mv archive/simulations/<folder> simulations/` — the
  batch self-sync guard, `run_lammps.sh` and `postprocess.sh` all key on the folder
  name, so nothing else needs editing.
* `shear_slab.lmp` was re-laid-out on 2026-09-22 to follow these decks' section
  structure and conventions (see its header), so the archived compression deck is
  also the reference for reading the shear deck.

## Archived tools (2026-08): the volume-of-mixing / cavity-Widom pipeline

Retired in Phase 0 of the calibrated volume-fraction plan (2026-08), which
supersedes the old volume-of-mixing pipeline with partial-molar-volume NPT sweeps
calibrated against the notebook Voronoi estimator (see
`simulations/calibration_sweep/README.md`).

| File | Came from | Why it is here |
|---|---|---|
| `volmix_sweep.sh` | `simulations/volmix_sweep/` | Old ΔV_mix SLURM pipeline (three-volume subtraction, per-pressure slab runs). Superseded by `calibration_sweep.sh`, which **inherits its SLURM chaining, manifest, and resume-flag machinery** — read this file to understand that machinery's origin. |
| `volume_of_mixing.ipynb` | `scripts/` | Analysis notebook for the volmix_sweep data. Its ΔV_mix(P) deliverable is replaced by the calibration analysis (fit-derivative error bars instead of three-volume subtraction). |
| `volume_of_mixing.png` | `scripts/` | Output figure of the notebook above. |
| `cavity_widom.py` | `scripts/` | Cavity-biased Widom insertion for μ_ex(z). The chemical-potential route to solvent activity is parked in favor of the volume-fraction calibration. |
| `clearance_sensitivity_cell.py` | `scripts/` | Standalone notebook cell probing the bounding-box clearance artifact in the old ΔV_mix — moot once ΔV_mix comes from closed periodic NPT boxes. |

### Related gating (not deleted, just off by default)

* The cavity-Widom trajectory dumps in the triaxial decks (two-piston and the
  archived one-piston ones) are gated behind `skip_widom` (default 1 = off).
  `slab_with_support.lmp` never had them.
* `run_lammps.sh` defaults `SKIP_WIDOM=1`; `postprocess.sh` and `continue_sim.sh`
  skip cavity-Widom post-processing gracefully when the trajectory is absent.
* The solvent density profiles written to `output_files/chemical_potential/`
  are **still produced** — the production notebooks read them.

### To re-enable cavity-Widom

Run with `SKIP_WIDOM=0` (or `-var skip_widom 0` for direct LAMMPS invocation)
and copy `cavity_widom.py` back to `scripts/` — `postprocess.sh`,
`continue_sim.sh`, and `plot_lammps_log.py` still look for it there.

## Archived (2026-09-29): the fixed-volume compress_slab workflow

`compress_slab` was rebuilt on 2026-09-29 on the template of
`triaxial_compression_two_pist` (wet pistons hold the bath at P*, six
solvent-transparent plates seated onto the gel, eps_vol = 0 reference window,
one file set per level tagged `_<hold>_c<level>`, analysis in `scripts/lib/bulk.py`).
The versions it replaced are kept here as they were at commit `a98f53e`.

| File | Came from | Why it is here |
|---|---|---|
| `simulations/compress_slab/` | `simulations/compress_slab/` | The 2026-08 deck and its Expanse / Pod / Bridges batch files: six plates ON the faces of a fully periodic fixed-volume box (Langevin+NPH to P* first), `fix move` ramps to a ladder of `stage_targets`, stage stresses appended into one file over a static `poly_bulk`/`solv_bulk` group and a fixed `bulk_vol`. Reads `-var stage_targets` (`COMPRESS_STAGES`), which the run scripts no longer pass. |
| `bulk_modulus_analysis.ipynb` | `scripts/` | Its notebook: the osmotic path K(t) = -V dPi/dV over the whole run (Rg volume) and the drained K from a fit over step-gap-clustered stages. Replaced by `scripts/bulk_modulus_analysis_{single,sweep}.ipynb`. |
| `add_more_plates_to_gel.ipynb` | `scripts/` | The notebook-only converter that put the plates on the box faces and deleted the polymer within 1 sigma of them. Replaced by `scripts/add_more_plates_to_gel.py` + the notebook of the same name. |

* **To analyse a run the old deck produced**, open the archived notebook from
  `scripts/` (its paths are relative to that folder).
* `plot_lammps_log.py` still draws its normal-stress panel when a run directory
  holds `bulk_modulus_plot_data_*.dat`; the new deck does not write that file.
