# Building Input Data Files

[← back to README](../README.md)

Before running a gel simulation you need a `.data` file — a text file describing every atom's position, bond topology, and atom types. These are generated on your MacBook using Jupyter notebooks in `scripts/`.

> **What is a Jupyter notebook?** An interactive Python document where you run code cells one at a time and see results immediately. Open JupyterLab from your terminal with `jupyter lab`, then navigate to the file.

### Which notebook to use

| Notebook | What it builds | When to use |
|----------|---------------|-------------|
| `add_walls_to_slab.ipynb` | Gel slab + flat support + piston (no chain angles) | Standard compression/flow runs |
| `add_walls_with_angles.ipynb` | Same but with harmonic and cosine angle terms | When running angle-restrained chains |
| `slab_with_support.ipynb` | Basic slab geometry builder | Reference / older geometry |
| `slab_with_support_angled.ipynb` | Angled-chain slab geometry | Angled geometry variants |
| `isolate_gel.ipynb` | Extracts the swollen polymer (+solvent) from a finished slab run | Pre-step for `shear_slab`; modulus analysis |
| `add_plates_to_gel.ipynb` | Attaches rigid shear plates on the x-faces of an isolated gel (atom type 4, harmonic-bonded to surface polymer) | **Required input for `shear_slab.lmp`** |
| `add_more_plates_to_gel.ipynb` | Variant that adds plates on all six faces | **Required input for `compress_slab.lmp`** (bulk modulus K) |
| `split_gel_slab.ipynb` | Splits a slab into polymer-only and solvent-only files | Isolated component analysis |
| `pure_polymer.ipynb` | Pure polymer box (no solvent) | EOS and baseline runs |
| `pure_solvent_1.ipynb` | Pure solvent box | EOS and solvent calibration |
| `slab_two_pistons.ipynb` (+ `slab_two_pistons.py`) | **CONVERTER**: an *equilibrated* `slab_with_support` periodic snapshot → the two-piston (feed / permeate) geometry: old piston deleted, permeate reservoir padded to ~10 σ, three new sheets (feed piston type 5, permeate piston type 6, dry piston type 7), vacuum margins, 7 atom types | **Required input for `triaxial_permeation_two_pist` / `triaxial_compression_two_pist`** (2026-09-16) |

### The two-piston converter (2026-09-16)

`slab_two_pistons.ipynb` is not a lattice generator: it **starts from the equilibrated**
`final_config_slab_support_periodic_…_14000002.data` (the 2026-09-07 piston-transparent rerun) and only
re-arranges the box.  The reason is physical: the two-piston decks hold each reservoir at its pressure with
the NPT-piston scheme of Marioni et al. (J. Membr. Sci. 738 (2026) 124837, Eq. 3), which acts in *z only* —
`lx`, `ly` are fixed for the whole run — so it could never swell a fresh `pre_swell = 0.93` lattice laterally.
Starting from the aniso-NPH slab (σ_p,xx/σ_p,zz = 1.0005) the gel arrives with zero transverse network stress
and its equilibrium `lx`, `ly` and thickness *by construction*; the NPT-piston phase in the decks is then a
z-settle of the reservoirs, not a swelling equilibration.  Gel and solvent coordinates are copied verbatim
(one uniform z shift; the converter checks this).  Output: `<input stem>_two_pist.data` next to the input,
a `.info.json` sidecar with every geometry number, a z-density histogram PNG, and an entry appended to
`slab_data_file_info.md`.  Knobs (Config cell / CLI): `permeate_thickness` (10 σ), `feed_thickness` (None =
keep), `margin_perm` (feed + 5 σ), `margin_feed` (15 σ — the feed piston rises by ~strain·L₀ under
compression, so raise it for sweeps beyond ε ≈ 0.10; the converter prints the limit), `piston_clearance`,
`dry_piston_frac`, `sheet_source` (`support` copies the input sheet's pattern, 23,316 beads; `hex` rebuilds
at `sheet_spacing`).  If the 14000002 file is not on your Mac the notebook falls back to the 14000000 snapshot
as a **smoke-test input only** and says so.

Typical shear-modulus pipeline: run `slab_with_support` to equilibrate → `isolate_gel.ipynb` to strip the support/piston → `add_plates_to_gel.ipynb` to attach plates → submit `shear_slab.lmp` with the `*_with_plates.data` file.

### Typical workflow for a new slab

1. Open `add_walls_to_slab.ipynb`
2. Set the unit cell dimensions (e.g. `10×10×8`), beads per chain, solvent density
3. Run all cells → generates a `.data` file in `../lammps_data_files_local/`
4. The output filename encodes all key parameters (e.g. `slab_support_5beads_tall_3.data`)
5. Copy to `~/Documents/lammps_data/input_data/` and then to the cluster

Generated file specs are logged in [`slab_data_file_info.md`](../slab_data_file_info.md) so you can always look up what was built and when.

---
