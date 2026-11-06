#!/usr/bin/env python3
"""build_benchmark_data.py -- every fresh data file the allocation benchmark needs (2026-11-05).

Existing files are never overwritten or rebuilt, so the production inputs in input_data/ are
reused as they are.

  scale 1.0  input_data/slab_support_periodic_<N>beads_tall_rho04_random.data, N = 5..30
             (N = 5 is the generator's default random lattice, 9 x 9 x 23 cells, pre_swell 0.93;
             N >= 10 uses the knobs of scripts/lib/build_beads_sweep.py at the same physical size)
  scale 0.5  input_data/allocation_benchmark/scale0.5/<same stem>.data, N = 5, 15, 30
  scale 2.0  input_data/allocation_benchmark/scale2.0/<same stem>.data, N = 5, 15, 30

The scale is the number of lattice cells, and so roughly the number of beads, relative to the
scale 1.0 slab.  The cell counts are listed explicitly because rounding a scaled box to whole cells
would leave N = 30 at 3 x 3 x 9 for "double": only 1.3 x.  A 2 x 2 x 8 network is too narrow for
physics; these files are for timing only.  Support and piston sheets scale with the area, solvent
with the volume, so the composition drifts a little with scale; collect_results.py records the
real counts.
"""
import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'scripts' / 'lib'))
from build_beads_sweep import SWITCH_SEEDS, knobs  # noqa: E402
from periodic_data_file_generator import build  # noqa: E402

ALL_N = (5, 10, 15, 20, 25, 30)
# (units_x = units_y, units_z); scale 1.0 is 1863 / 176 / 63 cells for N = 5 / 15 / 30
SCALED_CELLS = {
    0.5: {5: (7, 19), 15: (3, 10), 30: (2, 8)},    # 931 / 90 / 32 cells
    2.0: {5: (11, 31), 15: (5, 14), 30: (4, 8)},   # 3751 / 350 / 128 cells
}


def stem(n):
    return f'slab_support_periodic_{n}beads_tall_rho04_random'


def build_one(path, n, cells=None, check=True):
    if path.exists() and path.stat().st_size > 0:
        print(f'exists, kept: {path}', flush=True)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    if n == 5:
        pre_swell, placement = 0.93, 'straight'     # the generator's default random lattice
        ux, uz = cells if cells else (9, 23)
    else:
        k = knobs(n)
        pre_swell, placement = k['pre_swell'], 'walk'
        ux, uz = cells if cells else (k['units_xy'], k['units_z'])
    print(f'\n##### N = {n}: {ux} x {ux} x {uz} cells, pre_swell {pre_swell:.4f} -> {path}', flush=True)
    build('random', str(path), beads_per_chain=n, units_x=ux, units_y=ux, units_z=uz,
          pre_swell=pre_swell, chain_placement=placement,
          bond_switch_seed=SWITCH_SEEDS.get(n, 7), check_isotropy=None if check else False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input-dir', default=str(Path.home() / 'Documents/lammps_data/input_data'))
    ap.add_argument('--no-check', action='store_true', help='skip the isotropy and density reports')
    ap.add_argument('--scales', type=float, nargs='+', default=[1.0, 0.5, 2.0])
    ap.add_argument('--beads', type=int, nargs='+', default=None, help='restrict to these N')
    args = ap.parse_args()
    root = Path(args.input_dir)
    for scale in args.scales:
        if scale == 1.0:
            for n in ALL_N:
                if args.beads is None or n in args.beads:
                    build_one(root / f'{stem(n)}.data', n, check=not args.no_check)
        else:
            for n, cells in SCALED_CELLS[scale].items():
                if args.beads is None or n in args.beads:
                    build_one(root / 'allocation_benchmark' / f'scale{scale}' / f'{stem(n)}.data', n, cells,
                              check=False)


if __name__ == '__main__':
    main()
