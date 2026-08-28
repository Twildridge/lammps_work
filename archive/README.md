# archive/ — superseded / parked tools

Files here are **kept for reference, not maintained**. They were retired in
Phase 0 of the calibrated volume-fraction plan (2026-08), which supersedes the
old volume-of-mixing pipeline with partial-molar-volume NPT sweeps calibrated
against the notebook Voronoi estimator (see
`simulations/calibration_sweep/README.md` once Phase 2 lands).

| File | Came from | Why it is here |
|---|---|---|
| `volmix_sweep.sh` | `simulations/volmix_sweep/` | Old ΔV_mix SLURM pipeline (three-volume subtraction, per-pressure slab runs). Superseded by `calibration_sweep.sh`, which **inherits its SLURM chaining, manifest, and resume-flag machinery** — read this file to understand that machinery's origin. |
| `volume_of_mixing.ipynb` | `scripts/` | Analysis notebook for the volmix_sweep data. Its ΔV_mix(P) deliverable is replaced by the calibration analysis (fit-derivative error bars instead of three-volume subtraction). |
| `volume_of_mixing.png` | `scripts/` | Output figure of the notebook above. |
| `cavity_widom.py` | `scripts/` | Cavity-biased Widom insertion for μ_ex(z). The chemical-potential route to solvent activity is parked in favor of the volume-fraction calibration. |
| `clearance_sensitivity_cell.py` | `scripts/` | Standalone notebook cell probing the bounding-box clearance artifact in the old ΔV_mix — moot once ΔV_mix comes from closed periodic NPT boxes. |

## Related gating (not deleted, just off by default)

* The cavity-Widom trajectory dumps in `triaxial_compression.lmp` and
  `triaxial_permeation.lmp` are gated behind `skip_widom` (default 1 = off).
  `slab_with_support.lmp` never had them.
* `run_lammps.sh` now defaults `SKIP_WIDOM=1`; `postprocess.sh` and
  `continue_sim.sh` already skip cavity-Widom post-processing gracefully when
  the trajectory is absent.
* The solvent density profiles written to `output_files/chemical_potential/`
  are **still produced** — the production notebooks read them.

## To re-enable cavity-Widom

Run with `SKIP_WIDOM=0` (or `-var skip_widom 0` for direct LAMMPS invocation)
and copy `cavity_widom.py` back to `scripts/` — `postprocess.sh`,
`continue_sim.sh`, and `plot_lammps_log.py` still look for it there.
