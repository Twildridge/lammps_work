# Building Input Data Files

[← back to README](../README.md)

A gel simulation starts from a `.data` file: a text file with every atom's position and type and the bond topology. The files are built on the Mac with the Jupyter notebooks in `scripts/` (a notebook is a Python document run one cell at a time; `jupyter lab` in a terminal opens them) or with the scripts in `scripts/lib/`. They are kept in `lammps_data_files_local/` next to the repo, catalogued in its `README.md`, and copied to `~/Documents/lammps_data/input_data/` on the cluster.

### Which notebook or script to use

The notebooks live in `scripts/`. The Python modules they import live in `scripts/lib/`, which the notebooks add to `sys.path` themselves. A generator builds a fresh lattice; a converter re-arranges an equilibrated final configuration for a production deck. The two converters for current decks have their own sections below.

| Goal | Tool | Output file |
|------|------|-------------|
| Isolated gel block in a solvent bath (aligned, rotated or random network) | `data_file_generator.ipynb` + `lib/data_file_generator.py` | `slab_support_5beads_tall_rho04_nowalls.data`, `..._nowalls_rot45.data`, `..._nowalls_random.data` |
| Slab periodic in x and y | `slab_with_support_periodic.ipynb` + `lib/periodic_data_file_generator.py` | `slab_support_periodic_5beads_tall_rho04_new.data` (tetrahedral), `..._rho04_random.data` |
| Periodic slab with longer strands, or the rot45 standard gel | `lib/build_beads_sweep.py` | `slab_support_periodic_<N>beads_tall_rho04_<network>.data` |
| Two-piston input for `uniaxial_compression`, `uniaxial_permeation` and the periodic `simple_shear` deck | `slab_two_pistons_periodic.ipynb` + `lib/slab_two_pistons_periodic.py` | `<input stem>_two_pist.data` (`_two_pist_pm40` with `--margin-perm 40`) |
| Isolated gel with wet pistons and load plates, for `isotropic_compression` and the isolated-gel shear deck | `lib/add_more_plates_to_gel.py` (script only) | `<input stem>_<layout>.data`, for example `_six_plates`, `_x_plates` |
| Plated periodic slab for the archived `simple_shear_plates` deck | `lib/slab_shear_plates_periodic.py` (CLI only) | `..._with_plates_film7.data` |
| Pure polymer box (no solvent), for EOS and baseline runs | `pure_polymer.ipynb` | |
| Pure solvent box, for EOS and solvent calibration | `pure_solvent.ipynb` | |

#### `data_file_generator.ipynb`: isolated gel

The notebook sets the knobs and calls `lib/data_file_generator.py`, which builds a finite gel block (diamond lattice) in a solvent bath, optionally with a frozen support sheet and a piston sheet.

- `bridge_max_f = N` (2026-11-22, default 0 = off) adds one strand to every junction with N or fewer strands, to the nearest-to-typical-length neighbour that has a free site (`bridge_dangling`). The faces of an isolated block cut strands, which leaves dangling ends (one strand) and loop links (two strands) that swell outwards as the surface fringe. On the 10-cell random 10-bead cube `bridge_max_f = 2` adds about 257 strands (1.6 % of all): junctions with one strand go from 60 to 2, with two from 426 to 116. The interior isotropy is unchanged; the skin's fourth-rank anisotropy rises from about 0.07 to 0.08–0.09, because the new strands lie along the faces. Bridged junctions fluctuate less, so more bridging trades fringe against surface mobility.
- The solvent bath is placed with the same overlap test as before (against the last 100 particles) through a ring buffer, about 20 times faster with the same beads accepted.
- `chain_placement='walk'` (2026-11-22, as the periodic generator) lays each strand as a random walk pinned to its two junctions. Use it for strands longer than 5 beads, whose beads would overlap on the straight line; `scripts/submit_cube_pipeline.sh` does so for its 10-bead cube.
- `includeWalls = False` writes an isolated gel with no sheets. Five atom types are still declared. `pistonMargin` is ignored in this mode, because a margin without sheets only wrote a vacuum frame round the box (as in the `square_rho10_nowalls` cube).
- `NETWORK = "tetrahedral" | "rot45" | "random"` selects the preset (`gen.NETWORKS`).
- The isolated gels are the starting point for `isotropic_compression` (K) and the isolated-gel shear deck. Equilibrate them with `slab_with_support_single.batch` (`BARO_MODE=iso`, `PIN_GEL=1`, `PISTON_TRANSPARENT=1`; `DATANAME` from the environment, default the random gel's stem).

**`tetrahedral`** builds the tall isolated gel `slab_support_5beads_tall_rho04_nowalls.data`:

- 9 × 9 × 21 cells (the July-8 gel), 435 k atoms.
- `preSwell = 1.10`, which gives a = 5.28 σ, the equilibrium lattice constant of this gel at P* = 1.5, 1.0_1.0 (measured on the 14000002 slab: 5.33 in-plane, 5.26 along z).
- `rhoSolv = 0.40` over the whole box. The bath at P* = 1.5 is 0.48 and the gel interior holds 0.22 solvent + 0.37 polymer, so a uniform 0.40 leaves the 15 σ pads at 15 σ.
- 15 σ of bath on every side; the bottom pad is `belowThickness + 2·supportThickness`.

**`rot45`** (`rot45_y = True` in the library) builds the same network with the diamond lattice turned 45° about y: box x = crystal [101], z = [−101], y = [010].

- The rotated lattice is periodic in the box frame with the body-centred-tetragonal cell a/√2 × a × a/√2, so `units_x` and `units_z` count those cells (3.73 σ) and `units_y` counts cubic cells.
- The preset builds `slab_support_5beads_tall_rho04_nowalls_rot45.data`: 13 × 9 × 30 cells = 48.5 × 47.5 × 112.0 σ (the aligned block is 47.5 × 47.5 × 110.9), 15257 crosslinks, 448 k atoms.
- Purpose: in the isolated-gel shear deck (x gap, z shear) an xz shear of the aligned gel is a cube-face shear of the lattice, C44, and of the rotated gel a cube-diagonal shear, (C11 − C12)/2, the constant every compression-side G measures. The two runs decide whether the shear G (0.15–0.18) and the compression G (0.04–0.06) differ by lattice anisotropy.
- A front view of the rotated gel shows an axis-aligned square mesh (strands projected along x and z). The aligned gel shows a ±45° mesh, because its strands are the cube's body diagonals.

**`random`** (`bond_switch = N`) builds `slab_support_5beads_tall_rho04_nowalls_random.data` from the rot45 geometry with `bond_switch = 3`, seed 7; see [Random networks](#random-networks-bond-switching-and-the-build-checks). The build takes about three minutes, and the output of the checks is kept in the notebook.

#### `slab_with_support_periodic.ipynb`: periodic slab

The notebook sets the knobs for `lib/periodic_data_file_generator.py` with `NETWORK = "tetrahedral" | "random"`. The library is also a CLI:

```bash
python lib/periodic_data_file_generator.py --network random --output <file>
```

- CLI defaults: `--units 9 9 23`, `--beads-per-chain 5`, `--rho-solv 0.4`, `--pre-swell 0.93`, `--chain-placement straight`. `--network` also accepts `rot45` and `random_u` (the random preset with a density-uniformity weight `w_uniform` = 0.3 in the switch acceptance).
- `random` builds `slab_support_periodic_5beads_tall_rho04_random.data`: the reference slab's 9 × 9 × 23 cells with the lattice re-wired by bond switching, periodic across x and y. Strands and repulsions across the x/y boundaries are minimum-image, and only the z faces are tethered. The recipe is that of the isolated random gel, starting from the aligned lattice.
- This file is the random-network input of the two-piston compression deck. `slab_with_support_periodic_random.batch` equilibrates it and `chain_random_periodic.batch` converts and submits; see `scripts/submit_random_periodic_pipeline.sh`.
- `rot45` is the periodic twin of the isolated rot45 gel, with the same cell counting (BCT cells along x and z, cubic along y), so the slab tiles the box. The sites are the conventional cells rotated and kept in the half-open box [0, L_x) × [0, L_y); the bonds are found by the minimum-image distance search of the aligned lattice. Checked on a 6 × 4 × 10 test slab: every interior junction four-coordinated, the (101) face junctions three-coordinated, 204 strands crossing the x/y boundaries.

#### `lib/build_beads_sweep.py`: longer strands and the rot45 standard gel

```bash
python scripts/lib/build_beads_sweep.py [--beads 10 15 20 30] [--dry-run]                 # random network
python scripts/lib/build_beads_sweep.py --network rot45 [--beads 5 15] [--dry-run]       # rot45 lattice
```

The driver calls `periodic_data_file_generator.py` and writes `slab_support_periodic_<N>beads_tall_rho04_<network>.data` to `--outdir` (default `~/Documents/lammps_data/input_data`).

- The generator's lattice constant is `1.2 (N - 1) pre_swell`, so at fixed cell counts the box grows like N and the atom count like N³. The driver holds the slab size instead, with fewer cells and a `pre_swell` that puts the lattice at the mesh the gel is expected to swell to, a_eq(N).
- Random network: ~43 σ in plane and ~111 σ thick (the 5-bead slab's size), a_eq(N) = 4.82 ((N − 1)/4)^0.59, from the equilibrated 5-bead random slab and the good-solvent exponent.
- rot45: see [The standard gel](#the-standard-gel-rot45-periodic-slab).
- `--box-xy` and `--thickness` override the size targets.
- Interior beads of longer strands go on a random walk between their two junctions (`--chain-placement walk` of the generator); on the straight line they would be 0.24–0.4 σ apart.
- The isotropy report's surface-recess line is calibrated on the 5-bead mesh (4 × 4 σ patches) and reads as voids for an 8–16 σ mesh. Ignore it there; the interior verdict is what counts.

#### Random networks: bond switching and the build checks

`NETWORK = "random"` in both generators re-wires the lattice by Wooten–Winer–Weaire bond switching on the crosslink graph (`www_bond_switch`): strand A–B is kept, A–C and B–D are replaced by A–D and B–C, until N switches per strand are accepted. This keeps every functionality, the strand count and length and the block faces, and leaves no lattice directions.

The junctions are relaxed between switch sweeps:

- Harmonic strand energy about the rest length, and a soft repulsion between non-bonded junctions.
- A weak tether of every junction to its lattice site, which spreads the stretch.
- Face junctions held in all three directions on anchors jittered by 0.25 r₀, so the lattice spacing is not re-formed. Held along the face normal only, at a quarter of the stiffness, junctions bonded to partners deep inside were pulled up to 3 σ into the block and slid sideways, leaving dents and loose loops on the surface (mean recess 1.3 σ, 17 % of 4 × 4 σ patches recessed by more than 2 σ).
- Metropolis acceptance annealed from T = 4 to 0.01 in units of the strand energy, after six randomising sweeps.
- No switch may re-create a lattice pair, and lattice strands that survive the randomising sweeps are switched away.
- A final pass re-wires every strand longer than 1.9 r₀ and frees the tethers of junctions on any that resist. The hard ceiling is 2.45 r₀ = 5.6 σ end-to-end, which puts beads 1.4 σ apart against the FENE R₀ of 1.5, after the harmonic push-off.
- `keep_size = False` (relaxing the block scale as well) is experimental and shrinks the block to half.

**Relaxation knobs.** The random presets pass `RANDOM_RELAX = dict(rest_len=1.4, rep_cut=1.3, w_rep=3.0)` as `switch_options`.

- Reason: with rest length r₀ and a repulsion only below 0.9 r₀ at weight 1, the junctions clump. A random f = 4 network wants strands ~1.5 r₀ long at the lattice's junction density, and a harmonic energy about r₀ shortens them by piling junctions together. The gel is then built with empty 5 σ voxels and 3–5 σ empty spheres (polymer-density CoV 0.38–0.41; the lattice has 0.035–0.05 and 1.8 σ), which the equilibration inherits and swells.
- With `RANDOM_RELAX` the strands end at 1.43 ± 0.21 r₀ = 2.75 σ (max 1.9 r₀) and the repulsion spaces the junctions like a soft-sphere liquid: the junction count within 2 r₀ has sd 2.0 (lattice 1.9, old relaxation 6.7).
- Periodic slab: CoV 0.09, emptiest voxel 70 % of the mean, largest empty sphere 2.1 σ, no over-long strands, interior isotropy OK, and a uniform 1.2 σ recess of the face junctions (the repulsion of the layer beneath).
- Isolated gel (13 × 9 × 30 rot45 cells): CoV 0.12, emptiest voxel 58 %, largest empty sphere 2.5 σ, isotropy OK.
- Both knobs are needed. CoV on a 9 × 9 × 8 test slab: `rest_len` 1.4 alone 0.34, `rep_cut` 1.3 alone 0.18, `rep_cut` 1.3 at weight 3 alone 0.14, `rest_len` 1.4 + `rep_cut` 1.1 at weight 3 0.19.
- Tried and rejected: `w_all` 0.5 (CoV 0.27); `rep_cut` 1.5–1.6 at weight 3–5 (0.10, but the face junctions were pushed 1.3–1.8 σ in and half the patches dented); `rest_len` 1.4 + `rep_cut` 1.2 (0.13); `w_teth` 8 (0.095, interior isotropy marginal).
- `switch_options={}` reproduces the old relaxation.
- Files built with `NETWORK = "random"` before 2026-11-03 have the voids and must be rebuilt: the random periodic slab and its uniaxial run, the isolated random gel and its shear runs, and any older beads-per-chain sweep file.
- A random f = 4 network spans less volume per strand than the open diamond lattice, so at the lattice block size it starts pre-stretched (about 28 % with the old relaxation); the equilibration run sets its real size.

**Isotropy report.** `verify_network_isotropy(<file>)` reads a written file back. It reports for all strands, for the interior (more than 7 σ, three strand lengths, inside the faces) and for the skin:

- the strand fabric tensor;
- the fourth-rank anisotropy of ⟨n n n n⟩ (the Born term of the elastic tensor);
- the Steinhardt Q₄ and Q₆ of the strand directions;
- the Born-affine Zener ratio A_xz = C₄₄ / ((C₁₁ − C₁₂)/2) in the box frame.

Each is compared with the same statistic for as many random directions. The second-rank fabric is 1/3 for every cubic lattice, so only the fourth-rank numbers discriminate.

The verdict is the interior's. The tethered, pre-stretched faces imprint an in-plane bias two to three strand lengths deep (⟨n_x²⟩ 0.27–0.31 within 4 σ of the x-faces, ⟨n_y²⟩ 0.27–0.30 at 2–8 σ from the y-faces, isotropic beyond 8 σ). It is not lattice memory: banning lattice pairs left it unchanged. The equilibration relaxes it with the pre-stretch, and the shear deck's grips and core margins sit on it.

Reference values for the isolated random preset (interior, against random directions and the lattices):

| Metric | Random gel | Random directions | Lattices |
|--------|-----------|-------------------|----------|
| Fourth-rank anisotropy | 0.017 | 0.016 ± 0.004 | 0.54 |
| Q₄ | 0.009 | 0.009 ± 0.002 | 0.51 |
| A_xz | 0.985 | 1.00 ± 0.01 | ∞ (aligned), 0 (rot45) |

All strands together fail Q₄ through the skin (0.013), and the two lattices fail everything, as they must.

The report's `surface` line measures the recess of the face junctions from the 0.2 / 99.8-percentile planes of the junction cloud and counts separately the junctions that poke more than 1 σ out of them. The relaxation pushes a few junctions out through a face (30 of 15 257 on the isolated block, 5 on the periodic slab), and measured from the outermost junction those would read as a 1.2–1.8 σ recess of every patch. Current builds: isolated gel 0.46 σ mean recess, 1.3 % of patches recessed more than 2 σ, none more than 4 σ; periodic slab 0.03 σ and no dents.

**Density report.** Every random build also prints `verify_density_uniformity(<file>, periodic=...)`:

- polymer density on 5 σ voxels more than 6 σ inside the gel: its CoV and the emptiest and fullest voxel;
- the largest empty sphere, from a 0.5 σ grid;
- the volume more than 3 σ from any bead, and how many voids it makes.

The verdict is `DENSITY OK` at CoV ≤ 0.15 and no 3 σ empty sphere. The voxel, margin and probe scale with a / 4.46, so a longer-strand mesh is judged per mesh.

#### `lib/slab_shear_plates_periodic.py`: plated periodic slab (archived deck)

This CLI converter builds the input of `simple_shear_plates.lmp`, which is archived in `archive/simulations/simple_shear_plates/` (see [`archive/README.md`](../archive/README.md)); the shear notebooks still read its runs.

```bash
python scripts/lib/slab_shear_plates_periodic.py
```

- The default input is the 14000002 slab in `lammps_data_files_local/`, the two-piston converter's input, so the shear G and the compression G come from one gel state.
- It isolates the slab (`isolate_gel.py`: x/y periodicity kept, bath solvent trimmed in z) and puts plates normal to z on its two free faces (`add_plates_to_gel_periodic.py`, `axis='z'`), which the deck drives ±x. The shear gap is the slab thickness (plate separation ~120 σ) and the sheared membrane is laterally periodic (48 × 48 σ).
- It deletes the solvent beads within `--solvent-delete` (0.87 σ, calibrated) of a plate atom, so the plates' excluded volume does not compress the gel (bulk gel pressure in the shear hold 1.50 instead of 1.53).
- Other defaults: `--clearance 4.0`, `--spacing 1.5`, `--offset 0.5`, `--cutoff 2.5`, `--box-buffer 2.0`.
- It writes `isolated_..._14000002_with_plates_film7.data`, an `_isolated.data` intermediate and a `slab_data_file_info.md` entry. The archived `simple_shear_plates.batch` names this file in `DATANAME`.

`lib/add_plates_to_gel_periodic.py` is the plate builder behind it: rigid plates (atom type 4, harmonic-bonded to surface polymer) on the faces normal to `axis` (`'x'` = the isolated-cube geometry, `'z'` = the periodic slab).

#### Archived builders

The one-piston-era notebooks are in `archive/scripts/` and described in [`archive/README.md`](../archive/README.md): `add_walls_to_slab.ipynb`, `add_walls_with_angles.ipynb`, `slab_with_support_angled.ipynb`, `isolate_gel.ipynb`, `split_gel_slab.ipynb` and `add_plates_to_gel.ipynb`.

The CLI modules `lib/isolate_gel.py` (used by `calibration_sweep` and `slab_shear_plates_periodic.py`) and `lib/split_gel.py` (used by `calibration_sweep`) stay in `scripts/lib/`. The isolated-cube shear pipeline (`isolate_gel.ipynb`, then `add_plates_to_gel.ipynb`; plates on the x-faces of a finite gel driven ±z) needs the pre-2026-09-22 `simple_shear_plates.lmp` and the pre-2026-09-23 shear notebook from git history.

### Running the isolated-gel converter from the Terminal (2026-10-14)

`lib/add_more_plates_to_gel.py` converts an equilibrated isolated-gel snapshot into a production input. It has no notebook. The terminal inside JupyterLab (**File → New → Terminal**) has the Python that runs the notebooks.

```bash
cd "<your lammps_work folder>/scripts/lib"            # quote the path, it has spaces
python add_more_plates_to_gel.py --layout x_plates \
    --input ../../../lammps_data_files_local/final_config_slab_support_5beads_tall_rho04_nowalls_1.0_1.0_8000000.data
```

`python` must be the Anaconda Python that runs JupyterLab. `which python` shows which is active and `conda activate base` switches to it; on Dylan's Mac it is `/opt/anaconda3/bin/python`. The same line works in a notebook cell with a leading `!` and the path relative to the notebook (`!python lib/add_more_plates_to_gel.py ...`).

#### What it builds

The input must be the final configuration of a `slab_with_support` run with `BARO_MODE=iso` and `PISTON_TRANSPARENT=1` (`slab_with_support_single.batch`). The gel is made whole, free beads are removed, the bath is cropped, two wet pistons close z, and the load plates of the chosen layout are parked 3 σ outside the gel faces.

| `--layout` | Plates | Deck |
|------------|--------|------|
| `six_plates` (default) | All six faces | `isotropic_compression.lmp` |
| `x_plates` | The two x-face grip plates | `simple_shear_isolated.lmp` |
| `z_plates` | Support + load plate (a uniaxial-style isolated-gel input) | No deck yet |
| `no_plates` | None | |

The 11-type table is fixed and absent plates are empty, so any deck reads any layout. A `six_plates` file also runs in the isolated-gel shear deck, which deletes the four plates it does not use.

#### Options

`--help` lists every option. The main ones:

- `--bath-margin-xy` (12 σ) and `--reservoir-z` (15 σ) set the bath kept round the gel.
- `--reservoir-feed` and `--reservoir-perm` set the two reservoirs separately (default `--reservoir-z`), for a deeper feed side or more room under the permeate piston.
- `--plate-park-gap` (3 σ) and `--plate-overhang` (6 σ) place the plates.
- `--no-log` skips the `slab_data_file_info.md` entry.
- `--png` also draws `<stem>_profiles.png`, the density profiles along z and x. Use it to check that the gel sits between the plates and pistons.

#### Output and checks

The script writes `<input stem>_<layout>.data` next to the input. It prints the geometry, the validation (`RESULT: OK`, or the problems, in which case nothing is written), the gel-unchanged check, the `slab_data_file_info.md` entry it appended and, last, `DATANAME for the batch file: <stem>`.

#### Submitting

Copy the `.data` file to `~/Documents/lammps_data/input_data/` on the cluster.

- `six_plates`: `DATANAME=<stem> sbatch isotropic_compression_single.batch`. Both `isotropic_compression` batches read `DATANAME` from the environment; the default is the square rho-1.0 cube.
- `x_plates`: put the stem in `simple_shear_isolated_single.batch` (its `DATANAME` line).

The periodic-slab converters (`slab_two_pistons_periodic.py`, `slab_shear_plates_periodic.py`) run the same way; `slab_two_pistons_periodic.ipynb` is the notebook front-end of the first.

### The two-piston converter (2026-09-16)

`slab_two_pistons_periodic.ipynb` (code in `lib/slab_two_pistons_periodic.py`) is not a lattice generator. It starts from an equilibrated periodic `slab_with_support` snapshot and re-arranges the box into the two-piston (feed / permeate) geometry: the old piston is deleted, the permeate reservoir is padded to ~10 σ, and three new sheets are added (feed piston type 5, permeate piston type 6, load piston type 7) with vacuum margins, for 7 atom types in all. It is the required input of `uniaxial_compression` and `uniaxial_permeation`.

#### Why it starts from an equilibrated slab

The two-piston decks hold each reservoir at its pressure with the NPT-piston scheme of Marioni et al. (J. Membr. Sci. 738 (2026) 124837, Eq. 3), which acts in z only. `lx` and `ly` are fixed for the whole run, so the scheme could never swell a fresh `pre_swell = 0.93` lattice laterally. Starting from the aniso-NPH slab (σ_p,xx/σ_p,zz = 1.0005), the gel arrives with zero transverse network stress and its equilibrium `lx`, `ly` and thickness by construction. The NPT-piston phase in the decks is then a z-settle of the reservoirs, not a swelling equilibration.

#### Input

The notebook default is `final_config_slab_support_periodic_…_14000002.data`, the 2026-09-07 piston-transparent rerun. If that file is not on the Mac, the notebook falls back to the 14000000 snapshot as a smoke-test input only and says so.

#### Options (Config cell or CLI)

- `permeate_thickness` (10 σ).
- `feed_thickness` (None = keep).
- `margin_perm` (None = feed + 5 σ): vacuum below the permeate piston.
- `margin_feed` (15 σ): vacuum above the feed piston. The feed piston rises by ~strain·L₀ under compression, so raise it for sweeps beyond ε ≈ 0.10; the converter prints the limit. The rot45 gels need 30 σ: at 15 σ a 0.10 strain on these taller gels halts at ε ≈ 0.08–0.09 (see `running_simulations.md`).
- `piston_clearance` (1.0 σ).
- `load_piston_frac` (0.5).
- `sheet_source`: `support` (default) copies the input sheet's pattern, 23,316 beads; `hex` rebuilds at `sheet_spacing` (0.2 σ).

#### Output and checks

- `<input stem>_two_pist.data` next to the input, and an entry appended to `slab_data_file_info.md`. `--png` also draws the z-density histogram.
- A file converted with `--margin-perm 40` is named `<stem>_two_pist_pm40` (`chain_random_periodic.batch` with `MARGIN_PERM=40` applies the name; by hand, pass `--output`).
- Gel and solvent coordinates are copied verbatim with one uniform z shift; the converter checks this.

#### Which shear deck reads which file

`simulations/simple_shear/` holds two decks:

- The laterally periodic deck (`simple_shear_periodic.lmp`, the folder default; `simple_shear.lmp` and `simple_shear_single.batch` point to it) needs no data file of its own. It reads the `<stem>_two_pist.data` or `_two_pist_pm40.data` file unchanged: the support and the load piston stay parked and the shear is applied through the polymer face layers. Anything built for `uniaxial_compression` can be sheared as is.
- The isolated-gel deck (`simple_shear_isolated.lmp`, `simple_shear_isolated_single.batch`) reads the `x_plates` file of `lib/add_more_plates_to_gel.py`.

### Workflow for a new gel

1. Open `data_file_generator.ipynb` (isolated gel; pick `NETWORK`) or `slab_with_support_periodic.ipynb` (x/y-periodic slab), or use `lib/build_beads_sweep.py` for longer strands and the rot45 slab.
2. Set the unit-cell count, beads per chain, pre-swell and solvent density.
3. Run all cells. The `.data` file lands in `../../lammps_data_files_local/`; its name encodes the parameters.
4. Copy it to `~/Documents/lammps_data/input_data/` on the cluster and run the equilibration batch.
5. Convert the equilibrated final configuration with the converter the production deck needs (table above).

#### The standard gel: rot45 periodic slab

The rot45 periodic slab is the standard example gel. The random networks were too heterogeneous to serve; the rotated lattice is the orientation on which the compression side measures (C11 − C12)/2, and longer strands on the same lattice are the route towards a less anisotropic lattice gel.

`python scripts/lib/build_beads_sweep.py --network rot45 [--beads 5 15]` builds `slab_support_periodic_<N>beads_tall_rho04_rot45.data` with:

- a_eq(N) = 5.28 ((N − 1)/4)^0.59: the 5-bead lattice's equilibrium constant at P* = 1.5 with the good-solvent exponent. The build starts at a_eq, where the 14000002 slab started at 0.93 × 4.8 and swelled 18 %.
- In-plane cells for a ~48 σ box (the equilibrated 14000002 slab is 48 × 48).
- Thickness 139.7 σ = 1.15 × the reference's 23 cubic cells at a_eq(5), rounded up to whole BCT cells. The extra 15 % gives the z-profiles more bins.
- Straight strands for 5 beads, the random walk for longer ones.

| | 5 beads | 15 beads |
|---|---------|----------|
| BCT cells | 13 × 9 × 38 | 9 × 7 × 18 |
| a_eq | 5.28 | 11.06 (`pre_swell` 0.658) |
| Box (σ) | 48.5 × 47.5 × 174.9 | 70.4 × 77.4 × 173.7 |
| Lattice height (σ) | 141.9 (+16.8 % over 23 cells) | 140.7 |
| Crosslinks | 18 018 | 4 662 |
| Strands | 35 802 | 9 198 |
| Atoms | 353 825 | 660 402 |

The 15-bead cross-section is 2.4× the 5-bead slab's (`BOX_XY_ROT45_N`: ~74 σ instead of 48) so that the polymer beads per z-bin match: 124 236 polymer beads over 140.7 σ = 883 per σ, against 884 for the 5-bead file. The 15-bead gel's polymer density at its mesh is 2.4× lower, and matching the count keeps the stress and volume-fraction profiles as quiet per bin despite the fewer, floppier chains. The price is 1.9× the atoms, so the two-piston runs need about twice the wall time or twice the nodes. Cheaper 15-bead fallbacks: 6 × 4 × 18 cells at 48 σ (252 k atoms) and 7 × 5 × 18 at 55 σ (367 k, `--box-xy 55`).

With walk placement a few hundred bonds start longer than the FENE R₀ of 1.5 σ and some polymer pairs closer than 0.5 σ. The harmonic push-off of `slab_with_support.lmp` Phase 0 removes both, as it did for the random sweep files.

The lattice constant is the c* estimate. The aniso NPH equilibration sets the real in-plane size, and a gel that swells past a_eq gives a proportionally larger and more expensive two-piston box.

Pipeline:

1. Equilibrate with `NETWORK=rot45 bash submit_chain_length_sweep.sh [5 15]` in `simulations/chain_length_testing/` (`chain_length_testing_sweep.batch`: aniso NPH at P* = 1.5, Rg-plateau halt). The step budgets 14 M / 30 M are the `final_config` tags; the 2026-11-10 runs were extended to 20 M / 60 M and their stems carry those tags.
2. Convert the final configuration once per gel with `slab_two_pistons_periodic.py --margin-perm 40 --margin-feed 30` to `<stem>_two_pist_pm40`, which both two-piston decks and the periodic shear deck read.
3. Run `uniaxial_compression_single.batch` and `uniaxial_permeation_single.batch` with `DATANAME` exported.

`chain_random_periodic.batch` does the copy, conversion and submission in one job (`running_simulations.md`).

#### Catalogue

Every generated file gets an entry in [`slab_data_file_info.md`](../slab_data_file_info.md). The folder is catalogued in `lammps_data_files_local/README.md`: one row per file with the measured solvent density, provenance, the batch that used it and its status. Active files stay at the top level; obsolete lineages are in `archive/<date>_<lineage>/`. The `rho` in a raw file's name is the generator's nominal density and was wrong for some files (`rho6` = 1.0; the cubes `new_02` = 0.6, `new_03` = 1.0). The two current ones were renamed to `rho10` (= 1.0) on 2026-10-11.

---
