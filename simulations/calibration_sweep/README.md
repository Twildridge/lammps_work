# calibration_sweep — calibrated solvent volume-fraction pipeline

**Supersedes `volmix_sweep` (archived in `lammps_work/archive/`); the SLURM
chaining / manifest / resume machinery is inherited from it.** This sweep
measures partial molecular volumes and the Voronoi→thermodynamic calibration
λ(φ_p, P) from homogeneous periodic-gel NPT runs; ΔV_mix(P) over the full
pressure range is a by-product with fit-based error bars.

## The math

Thermodynamic solvent volume fraction (the one Flory–Rehner, Π, K(φ_p), D_c
are written in):

    φ_f = N_f v̄_f / V,    v̄_f = (∂V/∂N_f)_{T,P,N_p}

Euler (V extensive in N_f, N_p at fixed T,P):

    V = N_f v̄_f + N_p v̄_p   →   φ_p = 1 − φ_f exactly; no dry-reference
                                  volume; v̄_p by closure, never fit separately.

Per pressure, the analysis fits ln⟨V⟩ vs ln N_f (spline/cubic) across the
composition grid and differentiates analytically:

    φ_f^th(N_f, P) = ∂ln⟨V⟩/∂ln N_f |_P
    v̄_f = φ_f^th · V / N_f ;    v̄_p = (V − N_f v̄_f)/N_p   (Euler closure)

The same boxes get a Python-side Voronoi pass (`scripts/lib/volfrac.py`, the
identical code path as production analysis):

    φ_f^vor = Σ_{i∈solvent} v_i^voro / V

Calibration surface, joint fit over (φ_p, P) with the exact anchor λ(0,P)=1
(reservoir bins are NEVER shifted):

    λ(φ_p, P) = 1 + (a₁ + b₁P)φ_p + (a₂ + b₂P)φ_p²

ΔV_mix at scale = 1 (per-loading pure-solvent companions hold EXACTLY the
mixed box's atoms; polymer companion once per P):

    ΔV_mix(N_f, P) = ⟨V_mix⟩ − ⟨V_sol(N_f)⟩ − ⟨V_pol⟩

## The grids

All grids/counts are variables at the top of `calibration_sweep.sh` — nothing
is buried, nothing is hard-coded in any `.lmp` (fully `-var` driven via
`run_lammps.sh`).

| Variable | Default | Meaning |
|---|---|---|
| `PRESSURES` | 0.50 … 2.00 step 0.25 | pressure grid |
| `NF_GRID` | 10 values, 55800 → 20800 | exact solvent counts (geometric, ratio ≈ 0.896 from the base snapshot's equilibrium N_f = 55868; N_p = 104490 fixed) |
| `NREPS` | 2 | NPT thermal replicas per grid point |
| `CALIB_STEPS` / `_LOWP` | 300k / 500k (P ≤ 0.75) | NPT run length — verify convergence from the ⟨V⟩ trace, don't assume |
| `CALIB_FRAMES`, `CALIB_DUMP_EVERY` | 5, 2000 | all-atom Voronoi frames near each run's end (engines write frames+1) |

`calibration_sweep.sh` is the single source of truth for these. The config cell
of `scripts/calibration_analysis.ipynb` restates `NF_GRID`, `N_P`, `PRESSURES`
and the step counts and **must be updated with it** — the notebook addresses
runs by exact filename stem, so a stale grid doesn't error, it just reports the
mismatched loadings as missing. (Both this table and the notebook were left
stale by the 2026-08-29 regrid; fixed 2026-08-30.)

Base snapshot (`BASE_SNAPSHOT`): the 14M-step P=1.5
`slab_with_support_periodic` equilibration. No per-pressure slab runs — NPT
re-equilibration at each grid P erases the snapshot's origin.

**Before production:** set the low end of `NF_GRID` to span the max φ_p
reached in the c0.40 `triaxial_compression` profiles (validation needs the
calibration surface to cover the production range).

## How to run (Expanse login node)

```bash
cd ~/Documents/lammps_work/simulations/calibration_sweep

# smoke test: ONE pressure (prep + 41 runs), no auto-chaining
bash calibration_sweep.sh --only 1.5

# full sweep (7 pressures, auto-chained one pressure per batch)
bash calibration_sweep.sh

# resume from pressure index 3 (prep artifacts already on disk)
bash calibration_sweep.sh --from 3 --skip-prep
```

The one-time `[prep]` job runs `isolate_gel.py` → `adjust_solvent.py` (per
N_f) → `split_gel.py` (per N_f), all with skip-if-present logic, then the
per-pressure batches run the mixed NPT + solvent companion per (N_f, rep) and
one polymer companion per P. Chained batches auto-pass `--skip-prep`.

## Where outputs land

- Run dirs: `~/Documents/lammps_runs/calibration_sweep/<engine>_<dataname>_<interaction>_<timestamp>/`
- Manifests: `~/Documents/lammps_runs/calibration_sweep/sweep_manifest/`
  (`prep_*.path`, `p${P}_nf${NF}_rep${R}_{mixed,solvent}.workdir`,
  `p${P}_polymer.workdir`, `sweep_config_<ts>.txt`)
- Per run: `output_files/volume_data/box_dimensions_*.dat` (step lx ly lz),
  `output_files/volume_data/pressure_tensor_*.dat` (aniso-gate input),
  `traj_files/calib_*.lammpstrj` (Voronoi frames, on scratch)
- SLURM log: `calibration_sweep_<ts>.log` in this directory

## Analysis / regenerating the calibration json

`scripts/calibration_analysis.ipynb` syncs the sweep outputs to
`flow_data_local/calibration_sweep/` (paramiko + manifest, volmix style),
applies the aniso gate, fits φ_f^th and λ(φ_p,P), and writes
`scripts/calibration/calibration_lambda.json` (coefficients, covariance, full
raw table, metadata). Re-running the notebook end-to-end regenerates the json;
production notebooks load it via `lib/volfrac.py: load_calibration()`.

## Engine notes

`polymer_pure.lmp` (mixed + polymer boxes) and `solvent_pure.lmp` both run the
barostat **aniso** (2026-08) and log the full pressure tensor. The
convergence gate — Pxx≈Pyy≈Pzz within CI, stable aspect ratio, off-diagonals
flagged only — is checked in analysis, never enforced in-run. No Voronoi
computes in any `.lmp`; no silent switch to a triclinic barostat.
