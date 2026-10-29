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

* The runs these decks produced are still analysed by the one-piston notebooks, which
  moved to `archive/scripts/` on 2026-10-11 (`triaxial_compression_{single,sweep}.ipynb`,
  `triaxial_permeation.ipynb`, `triaxial_compression.ipynb`) — `lib/triaxial.py` handles
  one- and two-piston runs alike. Open them from `scripts/` (their paths assume that folder).
* `continue_sim.sh` still knows these folders (`final_tricomp` / `final_triperm`
  prefixes, `_c<level>` sweep tag); it works from a folder that holds the original
  SLURM output file next to the archived `.lmp` / `.batch`.
* **To run one again:** `git mv archive/simulations/<folder> simulations/` — the
  batch self-sync guard, `run_lammps.sh` and `postprocess.sh` all key on the folder
  name, so nothing else needs editing.
* `shear_slab.lmp` was re-laid-out on 2026-09-22 to follow these decks' section
  structure and conventions (see its header), so the archived compression deck is
  also the reference for reading the shear deck.

## Archived (2026-10-11): one-piston notebooks, old data-file builders, headless runner → `archive/scripts/`

`scripts/` was reduced to what the current pipeline uses (`scripts/README.md`). Moved here as they were:

| File | Why it is here |
|---|---|
| `scripts/triaxial_compression.ipynb` | Long-form one-piston analysis (the single/sweep notebooks were distilled from it 2026-09-02). |
| `scripts/triaxial_compression_single.ipynb`, `scripts/triaxial_compression_sweep.ipynb` | One-piston compression analysis; superseded by `scripts/triaxial_compression_{single,sweep}_two_pist.ipynb`. |
| `scripts/triaxial_permeation.ipynb` | One-piston permeation analysis; superseded by `scripts/triaxial_permeation_single_two_pist.ipynb`. |
| `scripts/add_walls_to_slab.ipynb`, `scripts/add_walls_with_angles.ipynb` | Built the `walled_*` inputs of `slab_with_flow` / the one-piston triaxial decks (Jan–Jun 2026). |
| `scripts/slab_with_support_angled.ipynb` | Angled-bond slab generator (Feb 2026 lineage; no current deck uses angle terms). |
| `scripts/isolate_gel.ipynb`, `scripts/split_gel_slab.ipynb` | Notebook front-ends of `lib/isolate_gel.py` / `lib/split_gel.py`, which stay (used by `calibration_sweep` and `slab_shear_plates_periodic.py`). |
| `scripts/add_plates_to_gel.ipynb` | Notebook front-end of `lib/add_plates_to_gel.py` (renamed `lib/add_plates_to_gel_periodic.py` on 2026-10-14, used by `slab_shear_plates_periodic.py`); it built the 2026-04 isolated-cube shear input. |
| `scripts/run_analysis_notebook.sh` | Headless notebook runner (2026-08-06); nothing calls it, `scripts/tests/run_notebook_tests.sh` covers the need. |

* All notebooks here use paths relative to `scripts/` (`lib`, `../../lammps_data_files_local/`, `../flow_data_local/`): copy one back to `scripts/` or open it from there to run it.
* The data files these builders produced are in `lammps_data_files_local/archive/` (catalog: `lammps_data_files_local/README.md`).

## Archived (2026-10-11): Pod and Bridges-2 → `archive/pod_bridges/`

Expanse is the only cluster in use. The Pod (CNSI, UCSB) and Bridges-2 (PSC) job runners,
every `*_pod.batch` / `*_bridges.batch` of the current decks (kept under
`archive/pod_bridges/simulations/<folder>/`, so they do not collide with the one-piston
decks' own Pod/Bridges batches in `archive/simulations/`), and the Pod benchmark notes
`documenting_pod_runs.md` are here as they were. `docs/cluster_reference.md` still describes
both machines. To run on one again: copy its runner back to `scripts/` and its batch back to
`simulations/<folder>/` — the batch `cd`s to `scripts/` and calls `./run_lammps_<cluster>.sh`,
and both runners still source `scripts/node_git_sync.sh`.

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

## Archived (2026-10-28): the two-batch rule, `shear_slab`, the stray data folder

Every folder under `simulations/` now keeps at most two batch files, `<folder>_single.batch` and `<folder>_sweep.batch`
(renamed from the unsuffixed files where there was one of each kind); everything else moved here, unchanged, under
`archive/simulations/<folder>/`:

| Folder | Archived | Why |
|---|---|---|
| `shear_slab/` (whole folder: deck + batch) | `shear_slab.lmp`, `shear_slab.batch` | The plate-driven shear of a periodic slab with plates on its z-faces (`slab_shear_plates_periodic.py` input). Superseded by **`simulations/shear_slab_two_pist/`** (isolated gel in the two-piston bath, grip plates). The shear notebooks still load its runs (`DECK = "shear_slab"`). |
| `shear_slab_two_pist/` | `shear_slab_two_pist_periodic.lmp` + `_periodic.batch`, `_fine.batch`, `_fine_g005.batch`, `_fine_g005_rot45.batch` | The laterally periodic two-piston shear deck (reads the compression deck's `*_two_pist.data`) and the fine-cadence displacement batches of the October D_c diagnosis (`DISP_FINE_NFREQ`; the deck no longer writes the fine or the hold-onset displacement files — D_c comes from the unload, Phase 2c). To run the periodic deck again: move both files back and keep `DECK_FILE=shear_slab_two_pist_periodic` in its batch. The rot45 run (`shear_6`) was `_fine_g005_rot45.batch`: the same run today is `shear_slab_two_pist_single.batch` with that DATANAME. |
| `slab_with_support/` | `slab_with_support.batch`, `_piston_transparent.batch`, `_isolated.batch`, `_isolated_tall.batch`, `_isolated_tall_rot45.batch` | The walled periodic slab runs (support + piston sheets), the rho10 cube and the aligned / rot45 tall gels. `slab_with_support_single.batch` (the former `_isolated_tall_random.batch`) runs any of them: set `DATANAME` (the gel from `data_file_generator.ipynb`), keep `BARO_MODE=iso`, `PIN_GEL=1` for an isolated gel. |
| `solvent_pure/` | `solvent_pure_old_no_piston_support.lmp` | The pre-sheet solvent deck. |
| `triaxial_compression_two_pist/` | `_half.batch`, `_quarter.batch` | The thinner gels (L² test of the hold relaxation, 2026-10-03/06: drainage confirmed, hold D_c thickness-independent). `_single.batch` is the former `_fine.batch` (one level, fine cadence, unload); `_sweep.batch` the former unsuffixed sweep. |
| `triaxial_permeation_two_pist/` | `_quarter.batch`, `_quarter_smalldp.batch` | The quarter gel and its small-dP ringing test (2026-10-07/08: amplitude explanation rejected). `_single.batch` is the former unsuffixed batch. |

Also removed the same day: `lammps_work/lammps_data_files_local/` (two August `..._with_six_plates_FIXED*.data` files that had been
committed by mistake; input data lives in `lammps_data_files_local/` NEXT TO the repo, never inside it) and every `.info.json`
sidecar and `_profiles.png` / `_zprofile.png` in that data folder (the converters stopped writing them; `--png` draws one on request).
