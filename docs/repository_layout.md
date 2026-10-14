# Repository Layout

[← back to README](../README.md)

```
lammps_work/                    ← This git repository
│
├── simulations/                ← One folder per simulation type
│   ├── slab_with_support/      ← Gel equilibration and compression (main workhorse)
│   ├── triaxial_compression_two_pist/  ← STANDARD compression workflow (two-piston, 2026-09-16; standard since 2026-09-22): NPT-pistons on both reservoirs + load piston; cumulative strain sweep
│   ├── triaxial_permeation_two_pist/   ← STANDARD permeation workflow (two-piston): feed/permeate NPT-pistons, one constant-dP drive (no sweep)
│   │                              (the one-piston triaxial_compression/ and triaxial_permeation/ moved to archive/simulations/ on 2026-09-22)
│   ├── shear_slab/             ← Shear modulus measurement (plate-driven xz shear; current G workflow; deck laid out like the triaxial decks since 2026-09-22)
│   ├── shear_slab_two_pist/    ← Shear modulus on the two-piston system. shear_slab_two_pist.lmp (2026-10-13) = ISOLATED gel from the six-plate converter file (plates
│   │                              deleted), rigid x-face grips driven -/+z, core stresses; shear_slab_two_pist_periodic.lmp (2026-10-03) = laterally periodic slab, z-face
│   │                              layers rigid in x/y driven +/-x (DECK_FILE in its .batch). Same notebooks, DECK = "shear_slab_two_pist_iso" / "shear_slab_two_pist"
│   ├── compress_slab/          ← Bulk modulus K: six plates + two wet pistons, on the template of triaxial_compression_two_pist (2026-09-29)
│   ├── solvent_phase/          ← Pure solvent equation of state sweep
│   ├── solvent_pure/           ← Single-state pure solvent run
│   ├── polymer_phase/          ← Pure polymer equation of state sweep
│   └── polymer_pure/           ← Single-state pure polymer run
│
├── scripts/                    ← CURRENT notebooks + the shell scripts (see scripts/README.md); EVERY .py lives in scripts/lib/ (since 2026-10-01)
│   ├── README.md               ← One line per file in this folder (2026-10-11)
│   ├── run_lammps.sh           ← Job runner for Expanse (the Bridges-2 / Pod runners are in archive/pod_bridges/ since 2026-10-11)
│   ├── postprocess.sh          ← Post-run plotting (calls lib/plot_*.py; run by the job runners and continue_sim.sh)
│   ├── git_sync.sh             ← One-command GitHub sync (MacBook + clusters) — this is what `lsync` runs
│   ├── install_lsync.sh        ← One-time-per-machine setup: registers `lsync` + nbstripout (see github_sync.md)
│   ├── node_git_sync.sh        ← Sourced by run_lammps.sh: best-effort git pull on the compute node
│   ├── continue_sim.sh         ← Resume a finished run without redoing setup (see running_simulations.md)
│   │
│   │   # ── data-file builders / converters (notebooks; converter code in lib/) ──
│   ├── slab_with_support.ipynb         ← Finite gel block in a bath, with or without sheets (includeWalls); current input: the tall isolated gel
│   ├── slab_with_support_periodic.ipynb ← xy-periodic crosslinked slab (bonds wrap x,y only; finite-z; p p p; one support+piston per z-period). The reference 14000002 slab
│   ├── slab_two_pistons.ipynb          ← CONVERTER notebook (code: lib/slab_two_pistons.py): equilibrated periodic slab → two-piston (feed/permeate) data file (input for triaxial_*_two_pist, shear_slab_two_pist)
│   ├── add_more_plates_to_gel.ipynb    ← CONVERTER notebook (code: lib/add_more_plates_to_gel.py) — six plates + two wet pistons, input for compress_slab
│   ├── pure_polymer.ipynb              ← Build pure polymer data file (polymer_phase)
│   ├── pure_solvent_1.ipynb            ← Build pure solvent data file (solvent_phase)
│   │
│   │   # ── analysis notebooks (code in lib/) ──
│   ├── triaxial_compression_single_two_pist.ipynb ← STANDARD: two-piston compression, one level: the eleven figures + wet-piston bath check + solvent expelled
│   ├── triaxial_compression_sweep_two_pist.ipynb  ← STANDARD: two-piston compression, all levels (M, G, D_c, κ vs strain) + per-level bath check
│   ├── triaxial_permeation_single_two_pist.ipynb  ← STANDARD: two-piston permeation: profile evolutions, P_feed/P_perm measured vs applied, Q_perm(t), permeability k, D_c + M from the polymer displacement (no sweep)
│   ├── shear_analysis_single.ipynb          ← one shear level: G (network + series estimators, increments from γ = 0), N1/N2, D_c, κ, P_th (shear_slab or shear_slab_two_pist output: DECK in Config)
│   ├── shear_analysis_sweep.ipynb           ← every level of a shear sweep: G vs γ, stress–strain slope, D_c, κ
│   ├── bulk_modulus_analysis_{single,sweep}.ipynb ← Bulk modulus K (code in lib/bulk.py, on top of lib/triaxial.py)
│   ├── calibration_analysis.ipynb           ← φ_f calibration / PMV sweep analysis (lib/volfrac.py)
│   ├── cross_deck_comparison.ipynb          ← figures built from more than one deck's run (lib/compare.py)
│   ├── (one-piston analysis notebooks, the walled/angled/isolate/split/add_plates builders and run_analysis_notebook.sh → archive/scripts/, 2026-10-11;
│   │    volume_of_mixing.ipynb, cavity_widom.py, clearance_sensitivity_cell.py → archive/, 2026-08)
│   │
│   ├── lib/                    ← ALL Python code (analysis modules, converters, CLI tools). Notebooks do
│   │   │                          `sys.path.insert(0, 'lib')`; shell scripts call `$SCRIPT_DIR/lib/<tool>.py`
│   │   │   # analysis modules (imported by the notebooks)
│   │   ├── triaxial.py                      ← all analysis code behind triaxial_*_two_pist / triaxial_compression_{single,sweep}.ipynb (Config, readers, Terzaghi, plateau bootstrap, G, D_c, figures)
│   │   ├── shear.py                         ← all shear analysis code for both shear decks (mirrors triaxial.py and imports its machinery; G_STRAIN = M_STRAIN's counterpart)
│   │   ├── bulk.py                          ← all bulk-modulus analysis code (on top of triaxial.py)
│   │   ├── volfrac.py                       ← Voronoi volume fraction + λ calibration (shared with calibration_analysis.ipynb)
│   │   ├── psd.py                           ← geometric porosity + pore-size distribution (2026-09-24)
│   │   │   # data-file converters / CLI tools (imported by the converter notebooks, or run as `python scripts/lib/<tool>.py`)
│   │   ├── slab_two_pistons.py              ← CONVERTER: equilibrated periodic slab → two-piston (feed/permeate) data file
│   │   ├── add_more_plates_to_gel.py        ← CONVERTER (2026-09-29): isolated gel → six plates + two wet pistons (compress_slab input); imports slab_two_pistons
│   │   ├── slab_shear_plates.py             ← CONVERTER (2026-09-22): periodic slab → isolate + z-normal plates on its faces = shear_slab input; imports isolate_gel + add_plates_to_gel
│   │   ├── add_plates_to_gel.py             ← plate builder (module extracted from the notebook 2026-09-22)
│   │   ├── isolate_gel.py                   ← CLI: strip bath/walls and define interior control volume (calibration_sweep, archived volmix_sweep)
│   │   ├── adjust_solvent.py                ← CLI: set N_f of an isolated gel exactly (calibration_sweep)
│   │   ├── split_gel.py                     ← CLI: split isolated gel into polymer-only / solvent-only (calibration_sweep, archived volmix_sweep)
│   │   │   # post-run plotters (called by postprocess.sh / continue_sim.sh / the job runners)
│   │   ├── plot_lammps_log.py               ← Plot T, P, volume convergence from log.lammps (+ shear diagnostics)
│   │   ├── plot_compression_strain_sweep.py ← Plot stress-strain / M across a triaxial_compression sweep
│   │   ├── plot_shear_strain_sweep.py       ← Plot stress-strain / G across a shear_slab sweep
│   │   ├── plot_stress_profiles.py          ← Plot stress and volume fraction profiles
│   │   ├── plot_piston_data.py              ← Plot piston position and velocity (multi-piston aware: one line per piston + pressure / flux panels)
│   │   └── plot_eos.py                      ← Plot P* vs ρ* for EOS sweeps (solvent_phase / polymer_phase)
│   │
│   └── tests/                  ← lint_lmp.py (deck linter), make_fixtures.py + run_plot_tests.sh + run_notebook_tests.sh + lib_headless_test.py
│
├── lammps_build/               ← build_lammps.sh (Expanse: 22Jul2025_update3, MPI+OpenMP, CMake) + the SLURM log of the build in use + README (2026-10-11; was scripts/)
│
├── lammps_data/                ← Reserved for small committed .data files (currently empty;
│                                  input_data/ lives outside the repo — see building_data_files.md)
│
├── docs/                       ← Everything below this repo-layout page: setup, running,
│                                  analysis, cluster reference, troubleshooting (see README.md)
├── archive/                    ← Superseded, unmaintained (see archive/README.md):
│   ├── simulations/triaxial_compression/  ← one-piston compression deck + batch files (archived 2026-09-22)
│   ├── simulations/triaxial_permeation/   ← one-piston permeation deck + batch files (archived 2026-09-22)
│   └── volmix_sweep.sh, volume_of_mixing.ipynb, cavity_widom.py, … ← 2026-08 volmix / Widom tools
├── .claude/commands/lsync.md   ← Makes `/lsync` work in Claude Code (see github_sync.md)
│
├── README.md                   ← Short entry point — links into docs/
├── shear_slab_notes.md         ← Append-only engineering decision log (gitignored) — NOT current-state; docs/ is
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
└── lammps_runs/        ← Timestamped output from each job (NOT in git)
    └── slab_*_20250101_120000/
        ├── data_files/     ← symlink to input .data file
        ├── output_files/   ← stress, volume, piston .dat files
        ├── output_plots/   ← auto-generated plots
        ├── traj_files/     ← symlink to scratch (large trajectory files)
        └── log.lammps
```

---
