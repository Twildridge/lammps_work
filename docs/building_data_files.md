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

Typical shear-modulus pipeline: run `slab_with_support` to equilibrate → `isolate_gel.ipynb` to strip the support/piston → `add_plates_to_gel.ipynb` to attach plates → submit `shear_slab.lmp` with the `*_with_plates.data` file.

### Typical workflow for a new slab

1. Open `add_walls_to_slab.ipynb`
2. Set the unit cell dimensions (e.g. `10×10×8`), beads per chain, solvent density
3. Run all cells → generates a `.data` file in `../lammps_data_files_local/`
4. The output filename encodes all key parameters (e.g. `slab_support_5beads_tall_3.data`)
5. Copy to `~/Documents/lammps_data/input_data/` and then to the cluster

Generated file specs are logged in [`slab_data_file_info.md`](../slab_data_file_info.md) so you can always look up what was built and when.

---
