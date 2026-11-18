# Repository Layout

[← back to README](../README.md)

```
lammps_work/                    ← This git repository
│
├── simulations/                ← One folder per simulation type
│   ├── slab_with_support/      ← Gel equilibration and compression (main workhorse)
│   ├── uniaxial_compression/  ← STANDARD compression workflow (two-piston, 2026-09-16; standard since 2026-09-22): NPT-pistons on both reservoirs + load piston; cumulative strain sweep
│   ├── chain_length_testing/       ← beads-per-strand sweep of the periodic random slab (2026-11-03): slab_with_support.lmp (symlink) with the Rg-plateau halt, one job per N (submit_chain_length_sweep.sh)
│   ├── uniaxial_permeation/   ← STANDARD permeation workflow (two-piston): feed/permeate NPT-pistons, one constant-dP drive, then a velocity-controlled stage and a frozen-network stage (no sweep)
│   │                              (the one-piston uniaxial_compression_one_pist/ and uniaxial_permeation_one_pist/ moved to archive/simulations/ on 2026-09-22)
│   │                              (simple_shear_plates/, the plate-driven periodic-slab shear, moved to archive/simulations/ on 2026-10-28)
│   ├── simple_shear/    ← Shear modulus on the two-piston system (the current G workflow). simple_shear.lmp (2026-10-13) = ISOLATED gel from the x_plates converter file,
│   │                              rigid x-face grips driven -/+z, core stresses, Phase 2c unload after every hold (2026-10-28: D_c read there); the periodic deck
│   │                              simple_shear_periodic.lmp + batch and the fine-cadence batches are in archive/simulations/simple_shear/ (2026-10-28)
│   ├── isotropic_compression/          ← Bulk modulus K: six plates + two wet pistons, on the template of uniaxial_compression (2026-09-29)
│   ├── pure_solvent/           ← Pure solvent, one deck: single state point (pure_solvent_single.batch) or pressure sweep / EOS (pure_solvent_sweep.batch, PHASE_SWEEP=1)
│   └── pure_polymer/           ← Pure polymer (or mixed gel), one deck: single state point (pure_polymer_single.batch) or pressure sweep / EOS (pure_polymer_sweep.batch, PHASE_SWEEP=1)
│
├── scripts/                    ← CURRENT notebooks + the shell scripts (see scripts/README.md); EVERY .py lives in scripts/lib/ (since 2026-10-01)
│   ├── README.md               ← One line per file in this folder (2026-10-11)
│   ├── run_lammps.sh           ← Job runner for Expanse (the Bridges-2 / Pod runners are in archive/pod_bridges/ since 2026-10-11)
│   ├── postprocess.sh          ← Post-run plotting (calls lib/plot_*.py; run by the job runners and continue_sim.sh)
│   ├── git_sync.sh             ← One-command GitHub sync (MacBook + clusters) — this is what `lsync` runs
│   ├── install_lsync.sh        ← One-time-per-machine setup: registers `lsync` + nbstripout (see github_sync.md)
│   ├── batch_prologue.sh       ← Sourced first by every .batch: git pull + re-exec, MPI env, modules
│   ├── node_git_sync.sh        ← Serialised fast-forward git pull on the compute node (used by batch_prologue.sh)
│   ├── lib/site.sh             ← Site settings: LAMMPS binary, paths, MPI_MAP_BY, WALL_MARGIN, modules
│   ├── lib/run_common.sh       ← Knob table + launch / verdict / chaining shared by run_lammps.sh and continue_sim.sh
│   ├── continue_sim.sh         ← Resume a finished run without redoing setup (see running_simulations.md)
│   │
│   │   # ── data-file builders / converters (notebooks; converter code in lib/) ──
│   ├── data_file_generator.ipynb       ← Finite gel block in a bath (was slab_with_support.ipynb until 2026-10-28): knobs only; NETWORK = tetrahedral | rot45 | random; code in lib/data_file_generator.py
│   ├── slab_with_support_periodic.ipynb ← xy-periodic crosslinked slab (bonds wrap x,y only; finite-z; p p p; one support+piston per z-period). The reference 14000002 slab. Knobs only since 2026-11-02 (NETWORK = tetrahedral | random); code in lib/periodic_data_file_generator.py
│   ├── slab_two_pistons_periodic.ipynb          ← CONVERTER notebook (code: lib/slab_two_pistons_periodic.py): equilibrated periodic slab → two-piston (feed/permeate) data file (input for uniaxial_*, simple_shear)
│   ├── pure_polymer.ipynb              ← Build pure polymer data file (pure_polymer sweep)
│   ├── pure_solvent_1.ipynb            ← Build pure solvent data file (pure_solvent sweep)
│   │
│   │   # ── analysis notebooks (code in lib/) ──
│   ├── uniaxial_compression_single.ipynb ← STANDARD: two-piston compression, one level: the eleven figures + wet-piston bath check + solvent expelled
│   ├── uniaxial_compression_sweep.ipynb  ← STANDARD: two-piston compression, all levels (M, G, D_c, κ vs strain) + per-level bath check
│   ├── uniaxial_permeation_single.ipynb  ← STANDARD: two-piston permeation: D_c + M from the polymer displacement, kappa from the velocity-controlled and frozen stages, the reconciliation kappa = D_c/M; then checks and the gel's state under flow (stresses, phi_s, pore sizes)
│   ├── shear_analysis_single.ipynb          ← one shear level: G (network + series estimators, increments from γ = 0), N1/N2, P_th; D_c and κ from the level's UNLOAD (2026-10-28); DECK in Config
│   ├── shear_analysis_sweep.ipynb           ← every level of a shear sweep: G vs γ, stress–strain slope, D_c (unloads), κ
│   ├── bulk_modulus_analysis_{single,sweep}.ipynb ← Bulk modulus K (code in lib/bulk.py, on top of lib/uniaxial.py)
│   ├── calibration_analysis.ipynb           ← φ_f calibration / PMV sweep analysis (lib/volfrac.py)
│   ├── cross_deck_comparison.ipynb          ← figures built from more than one deck's run (lib/compare.py)
│   ├── (one-piston analysis notebooks, the walled/angled/isolate/split/add_plates builders and run_analysis_notebook.sh → archive/scripts/, 2026-10-11;
│   │    volume_of_mixing.ipynb, cavity_widom.py, clearance_sensitivity_cell.py → archive/, 2026-08)
│   │
│   ├── lib/                    ← ALL Python code (analysis modules, converters, CLI tools). Notebooks do
│   │   │                          `sys.path.insert(0, 'lib')`; shell scripts call `$SCRIPT_DIR/lib/<tool>.py`
│   │   │   # analysis modules (imported by the notebooks)
│   │   ├── uniaxial.py                      ← all analysis code behind the uniaxial_* notebooks and the archived one-piston ones (Config, readers, Terzaghi, plateau bootstrap, G, D_c, figures)
│   │   ├── shear.py                         ← all shear analysis code for both shear decks (mirrors uniaxial.py and imports its machinery; G_STRAIN = M_STRAIN's counterpart)
│   │   ├── bulk.py                          ← all bulk-modulus analysis code (on top of uniaxial.py)
│   │   ├── volfrac.py                       ← Voronoi volume fraction + λ calibration (shared with calibration_analysis.ipynb)
│   │   ├── psd.py                           ← geometric porosity + pore-size distribution (2026-09-24)
│   │   │   # data-file converters / CLI tools (imported by the converter notebooks, or run as `python scripts/lib/<tool>.py`)
│   │   ├── slab_two_pistons_periodic.py              ← CONVERTER: equilibrated periodic slab → two-piston (feed/permeate) data file
│   │   ├── add_more_plates_to_gel.py        ← CONVERTER (2026-09-29, layouts 2026-10-14): isolated gel → plates by layout + two wet pistons (isotropic_compression / simple_shear inputs); imports slab_two_pistons_periodic
│   │   ├── data_file_generator.py           ← GENERATOR (2026-10-28, from slab_with_support.ipynb): diamond-lattice gel block in a bath; build(network='tetrahedral'|'rot45'|'random'), WWW bond switching (periodic= option 2026-11-02), isotropy check
│   │   ├── periodic_data_file_generator.py  ← GENERATOR (2026-11-02, from slab_with_support_periodic.ipynb): the x/y-periodic slab; build(network='tetrahedral'|'random'); CLI --network random --output <file>
│   │   ├── slab_shear_plates_periodic.py             ← CONVERTER (2026-09-22): periodic slab → isolate + z-normal plates on its faces = simple_shear_plates input; imports isolate_gel + add_plates_to_gel_periodic
│   │   ├── add_plates_to_gel_periodic.py             ← plate builder (module extracted from the notebook 2026-09-22)
│   │   ├── isolate_gel.py                   ← CLI: strip bath/walls and define interior control volume (calibration_sweep, archived volmix_sweep)
│   │   ├── adjust_solvent.py                ← CLI: set N_f of an isolated gel exactly (calibration_sweep)
│   │   ├── split_gel.py                     ← CLI: split isolated gel into polymer-only / solvent-only (calibration_sweep, archived volmix_sweep)
│   │   │   # post-run plotters (called by postprocess.sh / continue_sim.sh / the job runners)
│   │   ├── plot_lammps_log.py               ← Plot T, P, volume convergence from log.lammps (+ shear diagnostics)
│   │   ├── plot_compression_strain_sweep.py ← Plot stress-strain / M across a uniaxial_compression sweep
│   │   ├── plot_shear_strain_sweep.py       ← Plot stress-strain / G across a simple_shear sweep
│   │   ├── plot_stress_profiles.py          ← Plot stress and volume fraction profiles
│   │   ├── plot_piston_data.py              ← Plot piston position and velocity (multi-piston aware: one line per piston + pressure / flux panels)
│   │   └── plot_eos.py                      ← Plot P* vs ρ* for EOS sweeps (pure_solvent / pure_polymer with PHASE_SWEEP=1)
│   │
│   └── tests/                  ← lint_lmp.py (deck linter), make_fixtures.py + run_plot_tests.sh + run_notebook_tests.sh + lib_headless_test.py
│
├── lammps_build/               ← build_lammps.sh (Expanse: 22Jul2025_update3, MPI+OpenMP, CMake) + the SLURM log of the build in use + README (2026-10-11; was scripts/)
│
├── docs/                       ← Everything below this repo-layout page: setup, running,
│                                  analysis, cluster reference, troubleshooting (see README.md)
├── archive/                    ← Superseded, unmaintained (see archive/README.md):
│   ├── simulations/uniaxial_compression_one_pist/  ← one-piston compression deck + batch files (archived 2026-09-22)
│   ├── simulations/uniaxial_permeation_one_pist/   ← one-piston permeation deck + batch files (archived 2026-09-22)
│   └── volmix_sweep.sh, volume_of_mixing.ipynb, cavity_widom.py, … ← 2026-08 volmix / Widom tools
├── .claude/commands/lsync.md   ← Makes `/lsync` work in Claude Code (see github_sync.md)
│
├── README.md                   ← Short entry point — links into docs/
├── simple_shear_notes.md         ← Append-only engineering decision log (gitignored) — NOT current-state; docs/ is
├── .gitattributes              ← Runs nbstripout on every notebook commit (see github_sync.md)
├── lj_units_cheat_sheet.md     ← Unit conversions and parameter reference
├── expanse_lammps_guide.md     ← Expanse-specific setup and GitHub guide
├── slurm_commands_and_compiling.md  ← Cluster commands and LAMMPS build notes
├── slab_data_file_info.md      ← Log of every .data file that has been generated
└── (documenting_pod_runs.md → archive/pod_bridges/, 2026-10-11)
```

**What lives outside the repo** (on each machine):

```
~/Documents/
├── lammps_work/        ← This repo (git-tracked)
├── lammps_data/
│   └── input_data/     ← .data files (large; NOT in git — copy manually)
│       └── calibration_sweep/  ← calibration_sweep.sh prep files + run symlinks (~370 entries, kept out of the top level)
└── lammps_runs/        ← Timestamped output from each job (NOT in git)
    └── <folder>/                                   ← one folder per simulation type (same names as simulations/)
        └── <folder>_<dataname>_<interaction>_<YYYYMMDD_HHMMSS>/     ← one run; see "Inside a run folder" below
```

## Names changed on 2026-11-18

The simulation types were renamed for what they are. Decks, batches, notebooks, the library and
the docs use the new names; text written before that date was updated in place.

| Old | New |
|---|---|
| `triaxial_compression_two_pist` (folder, deck, batches) | `uniaxial_compression` |
| `triaxial_permeation_two_pist` | `uniaxial_permeation` |
| `shear_slab_two_pist` (also `_periodic`, `_isolated`) | `simple_shear` |
| `compress_slab` | `isotropic_compression` |
| `solvent_pure` + `solvent_phase` | `pure_solvent` (one deck; `PHASE_SWEEP=1` is the former `solvent_phase`) |
| `polymer_pure` + `polymer_phase` | `pure_polymer` (one deck; `PHASE_SWEEP=1` is the former `polymer_phase`) |
| `triaxial_compression_{single,sweep}_two_pist.ipynb`, `triaxial_permeation_single_two_pist.ipynb` | `uniaxial_compression_{single,sweep}.ipynb`, `uniaxial_permeation_single.ipynb` |
| `lib/triaxial.py`, `import triaxial as tri` | `lib/uniaxial.py`, `import uniaxial as uni` |
| `final_tricomp_*`, `final_triperm_*`, `strained_tricomp_*`, `strained_triperm_*` | `final_unicomp_*`, `final_uniperm_*`, `strained_unicomp_*`, `strained_uniperm_*` |
| archived one-piston `triaxial_compression`, `triaxial_permeation` (decks, notebooks) | `uniaxial_compression_one_pist`, `uniaxial_permeation_one_pist` |
| archived plate-driven `shear_slab` | `simple_shear_plates` |

Not renamed:

* Data files. `<stem>_two_pist.data`, `_two_pist_pm40.data` and so on keep the suffix: it is what
  tells the converted two-piston file from the periodic slab it was made from, and the files exist
  under those names on the cluster. `DATANAME` strings are unchanged everywhere.
* Local `RUN_ID` folders under `flow_data_local/` (some start with `two_pist_`).
* Run folders that already exist on the cluster. Runs made before the rename are still in
  `lammps_runs/triaxial_compression_two_pist/`, `.../triaxial_permeation_two_pist/`,
  `.../shear_slab_two_pist/` and `.../compress_slab/`; new runs go to the new folder names. The
  notebooks' sync searches both (`LEGACY_RUN_FOLDERS` in `lib/uniaxial.py`). Moving the old run
  folders into the new ones on the cluster makes that table unnecessary.
* The `two_pist` flag of `uniaxial.Config` and the `cfg.two_pist` tests in the libraries, which
  tell two-piston runs from the archived one-piston ones.
* The archived `solvent_pure` files under `archive/`.

## Inside a run folder

Layout since 2026-11-17. `run_lammps.sh` creates it (`make_run_dirs` in `scripts/lib/run_common.sh`) and writes the
same guide into every run folder as `README.txt`.

```
<folder>_<dataname>_<interaction>_<YYYYMMDD_HHMMSS>/
├── README.txt        ← this layout, so the folder explains itself
├── run_info.txt      ← what was run: SLURM job + its output file, deck, input data file, every -var given to LAMMPS, git commit
├── log.lammps        ← LAMMPS log: thermo output and the deck's own ">>>" progress lines
├── <deck>.lmp        ← the copy of the deck this run executed
├── data_files/       ← EVERY .data and .restart file of the run
│   ├── <dataname>.data                          ← the INPUT configuration (symlink into lammps_data/input_data/)
│   ├── final_<deck>_<stem>.data / .restart      ← the state at the END of the run
│   ├── final_config_<stem>_snap<N>.data         ← snapshots during an equilibration (resume points of continue_sim.sh)
│   └── strained_<deck>_<stem>_<stage>.restart   ← the strained state at the end of a level's hold, before any unload
├── output_files/     ← EVERY measured series and profile (.dat text tables; .dump for per-pair data)
│   ├── stress_data/         ← stress profiles and series, strain logs, reservoir pressures
│   ├── volume_data/         ← gel and box dimensions and volumes, polymer centre of mass
│   ├── piston_data/         ← piston / plate positions, velocities, forces, pressures
│   ├── permeation_data/     ← solvent crossing the gel (flux, expelled volume)
│   ├── displacement_data/   ← polymer displacement profiles
│   ├── pair_data/           ← per-pair and per-bond dumps (pore and bond analysis)
│   └── chemical_potential/  ← solvent density profiles
├── output_plots/     ← EVERY figure made from this run (.png)
│   ├── <stem>_<figure>.png  ← the results figures of postprocess.sh
│   └── convergence_plots/   ← run-health figures: T, P and volume convergence, diagnostics
├── traj_files/       ← symlink to scratch: trajectories (.lammpstrj). Scratch is purged; copy what must last.
└── continuation_<stamp>/    ← a continuation of this run (continue_sim.sh), same layout inside
```

**File names.** Every output file is `<quantity>_<stem>[_<stage>].<ext>`:

| Part | Meaning |
|------|---------|
| `<stem>` | `<dataname>_<epsSS>_<epsSP>_<steps>`; `<steps>` is the length of one hold (or of the production run) |
| `<deck>` (in `data_files/`) | `config` = equilibration (`slab_with_support`, `chain_length_testing`), `unicomp` = uniaxial compression, `uniperm` = uniaxial permeation, `shear`, `compress` = six-plate bulk compression; the pure-phase runs write `purepol_<stem>.data` / `puresolv_<stem>.data` |
| `_ref` | zero-load reference window, before any loading |
| `_c<strain>` / `_g<strain>` | the hold at that compression / shear level |
| `_u<strain>` | the unload after that level (re-swelling or relaxation; permeation: `_u`) |
| `_v`, `_f` | permeation: velocity-controlled stage, frozen-network stage |
| no stage | the whole run |

**The strained-state restart (2026-11-17).** The compression, shear and permeation decks write
`data_files/strained_<deck>_<stem>_<stage>.restart` at the end of every level's hold, before the unload starts
(permeation: at the end of the constant-dP drive, before the velocity-controlled, frozen and unload stages). It is the
equilibrated strained state, kept so that an unload (the trusted D_c) can be repeated without redoing the drive and the
hold. It is a binary LAMMPS restart, about 150 MB for the 770 k-atom wide gel, one per level; `SAVE_STRAINED=0` in the
batch turns it off. `read_restart` brings back atoms, velocities, bonds, box and groups, but not fixes or computes, and
no deck has a start-from-restart path yet: a repeat unload needs a deck that reads the restart and declares the unload
stage's fixes. `lmp -restart2data <file>.restart <file>.data` converts one to a text data file.

**Runs started before 2026-11-17** have their `final_*` and snapshot files at the top level of the run folder, and no
`README.txt`, `run_info.txt` or `strained_*` restart. `continue_sim.sh`, the pipeline batches and the plotters look in
`data_files/` first and then at the top level, so both layouts work.

---
