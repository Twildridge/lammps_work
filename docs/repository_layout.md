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
│   ├── compress_slab/          ← Bulk modulus K measurement (plate-driven isotropic six-face compression; undergrad project, in development)
│   ├── solvent_phase/          ← Pure solvent equation of state sweep
│   ├── solvent_pure/           ← Single-state pure solvent run
│   ├── polymer_phase/          ← Pure polymer equation of state sweep
│   └── polymer_pure/           ← Single-state pure polymer run
│
├── scripts/                    ← Shell and Python scripts
│   ├── run_lammps.sh           ← Job runner for Expanse
│   ├── run_lammps_bridges.sh   ← Job runner for Bridges-2
│   ├── run_lammps_pod.sh       ← Job runner for Pod
│   ├── build_lammps.sh         ← One-shot LAMMPS build script (Expanse-style cmake)
│   ├── git_sync.sh             ← One-command GitHub sync (MacBook + clusters) — this is what `lsync` runs
│   ├── install_lsync.sh        ← One-time-per-machine setup: registers `lsync` + nbstripout (see github_sync.md)
│   ├── continue_sim.sh         ← Resume a finished run without redoing setup (see running_simulations.md)
│   ├── isolate_gel.py          ← CLI: strip bath/walls and define interior control volume (used by volmix_sweep)
│   ├── add_walls_to_slab.ipynb         ← Build slab data file (no angles)
│   ├── add_walls_with_angles.ipynb     ← Build slab data file (with angles)
│   ├── slab_with_support.ipynb         ← Build basic slab geometry
│   ├── slab_with_support_periodic.ipynb ← CURRENT: xy-periodic crosslinked slab (bonds wrap x,y only; finite-z; p p p; one support+piston per z-period; no side padding). Input for triaxial_* runs
│   ├── slab_with_support_angled.ipynb  ← Build angled-chain slab geometry
│   ├── slab_two_pistons.ipynb / .py    ← CONVERTER: equilibrated periodic slab → two-piston (feed/permeate) data file (input for triaxial_*_two_pist)
│   ├── isolate_gel.ipynb               ← Extract just the swollen gel from a run
│   ├── split_gel_slab.ipynb            ← Split a gel slab into pieces
│   ├── add_plates_to_gel.ipynb / .py   ← Attach shear plates to a gel's faces (axis 'x' = isolated cube, 'z' = periodic slab; module extracted 2026-09-22)
│   ├── slab_shear_plates.py            ← CONVERTER (2026-09-22): periodic slab → isolate + z-normal plates on its faces = shear_slab input
│   ├── add_more_plates_to_gel.ipynb    ← Variant: plates on all six faces — input for compress_slab
│   ├── pure_polymer.ipynb              ← Build pure polymer data file
│   ├── pure_solvent_1.ipynb            ← Build pure solvent data file
│   ├── triaxial_compression_single.ipynb    ← one-piston runs (deck archived 2026-09-22): one strain level — φ_s, σ/σ′ evolution (zz,xx,yy), piston P, M, G, D_c, κ
│   ├── triaxial_compression_sweep.ipynb     ← one-piston runs: the same 11 figures overlaid for every sweep level; M, G, D_c, κ vs strain
│   ├── triaxial_compression.ipynb           ← long-form original (solvent-phase stress, ss/pp virial, Widom diagnostics); the two above were distilled from it
│   ├── lib/
│   │   ├── triaxial.py                      ← all analysis code behind triaxial_compression_{single,sweep}.ipynb (Config, readers, Terzaghi, plateau bootstrap, G, D_c, figures)
│   │   └── volfrac.py                       ← Voronoi volume fraction + λ calibration (shared with calibration_analysis.ipynb)
│   ├── triaxial_permeation.ipynb            ← one-piston runs (deck archived 2026-09-22): piston/thickness/stress/density/permeate + partial-vs-ss, with Phase 1.5 reference overlays
│   ├── triaxial_compression_single_two_pist.ipynb ← STANDARD: two-piston compression, one level: the eleven figures + wet-piston bath check + solvent expelled
│   ├── triaxial_compression_sweep_two_pist.ipynb  ← STANDARD: two-piston compression, all levels (M, G, D_c, κ vs strain) + per-level bath check
│   ├── triaxial_permeation_single_two_pist.ipynb  ← STANDARD: two-piston permeation: profile evolutions, P_feed/P_perm measured vs applied, Q_perm(t), permeability k (no sweep)
│   ├── tests/                           ← lint_lmp.py (deck linter), make_fixtures.py + run_plot_tests.sh + run_notebook_tests.sh + lib_headless_test.py
│   ├── bulk_modulus_analysis.ipynb          ← Drained vs osmotic bulk modulus K
│   ├── shear_analysis.ipynb                 ← G, N1/N2, stress profiles (shear_slab output)
│   ├── (volume_of_mixing.ipynb, cavity_widom.py, clearance_sensitivity_cell.py → moved to archive/, 2026-08)
│   ├── plot_lammps_log.py      ← Plot T, P, volume convergence from log.lammps (+ shear diagnostics)
│   ├── plot_compression_strain_sweep.py ← Plot stress-strain / M across a triaxial_compression sweep
│   ├── plot_shear_strain_sweep.py ← Plot stress-strain / G across a shear_slab sweep
│   ├── split_gel.py            ← CLI: split isolated gel into polymer-only / solvent-only (used by volmix_sweep)
│   ├── plot_stress_profiles.py ← Plot stress and volume fraction profiles
│   ├── plot_piston_data.py     ← Plot piston position and velocity (multi-piston aware: one line per piston + pressure / flux panels)
│   └── plot_eos.py             ← Plot P* vs ρ* for EOS sweeps (solvent_phase / polymer_phase)
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
└── documenting_pod_runs.md     ← HPC performance benchmarks
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
