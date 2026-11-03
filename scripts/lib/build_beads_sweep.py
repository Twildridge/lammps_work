#!/usr/bin/env python3
"""build_beads_sweep.py -- the beads-per-chain sweep of the x/y-periodic RANDOM slab (2026-11-03).

5 beads per strand is too short: the gel is elastic rather than osmotic (mixing).  This builds
one periodic random (bond-switched) slab per strand length, each at about the SAME physical size
as the 5-bead reference (slab_support_periodic_5beads_tall_rho04_random: 43.4 x 43.1 sigma in
plane, 116 sigma thick after its 14 M-step equilibration), so that the run cost per step and the
collective relaxation length stay comparable and only the mesh changes.

Knobs per N (all derived, printed in the summary):
  a_eq(N)     the lattice constant the gel is EXPECTED to equilibrate at: the 5-bead random slab
              ended at a = 43.4 / 9 = 4.82 sigma; a strand of N - 1 bonds in an athermal (WCA)
              solvent scales as (N - 1)^0.59, so a_eq = 4.82 * ((N - 1) / 4)^0.59.  The NPH aniso
              barostat corrects whatever this estimate misses; the point is to not start 2x
              pre-stretched (the lattice a = 1.2 (N - 1) would be).
  pre_swell   a_eq / (1.2 (N - 1)): the generator's scale knob that gives that a.
  units_x,y   round(43 / a_eq), the in-plane cells (box = units * a_eq); units_z round(111 / a_eq).
  placement   'walk': interior beads on a random walk pinned to both junctions (the straight line
              would put them 0.24-0.4 sigma apart).
Everything else is the reference's (rho 0.40, 4 sigma below the support, 25 sigma above the gel,
piston at 0.85, support spacing 0.2, bond_switch 3 / seed 7).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from periodic_data_file_generator import build  # noqa: E402

A_REF, N_REF, NU = 4.82, 5, 0.59          # the 5-bead random slab's equilibrated lattice constant
# bond_switch_seed per N where the default (7) left a detectable texture: the 20-bead file with
# seed 7 failed only the interior 4th-rank anisotropy (0.068 vs random 0.039 +/- 0.007); seed 11
# passes every metric (2026-11-03).
SWITCH_SEEDS = {20: 11}
BOX_XY_REF, THICK_REF = 43.0, 111.0       # its in-plane box and gel thickness (sigma)


def knobs(n_beads, box_xy=BOX_XY_REF, thick=THICK_REF):
    a_eq = A_REF * ((n_beads - 1) / (N_REF - 1)) ** NU
    pre_swell = a_eq / (1.2 * (n_beads - 1))
    units_xy = max(2, round(box_xy / a_eq))
    units_z = max(2, round(thick / a_eq))
    return dict(a_eq=a_eq, pre_swell=pre_swell, units_xy=units_xy, units_z=units_z)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--beads', type=int, nargs='+', default=[10, 15, 20, 30])
    ap.add_argument('--outdir', default=str(Path.home() / 'Documents/lammps_data/input_data'))
    ap.add_argument('--box-xy', type=float, default=BOX_XY_REF, help='target in-plane box (sigma)')
    ap.add_argument('--thickness', type=float, default=THICK_REF, help='target gel thickness (sigma)')
    ap.add_argument('--dry-run', action='store_true', help='print the knobs only')
    ap.add_argument('--no-check', action='store_true', help='skip the isotropy report')
    args = ap.parse_args()
    out = Path(args.outdir)
    rows = []
    for n in args.beads:
        k = knobs(n, args.box_xy, args.thickness)
        stem = f"slab_support_periodic_{n}beads_tall_rho04_random"
        print(f"\n##### N = {n}: a_eq {k['a_eq']:.3f}  pre_swell {k['pre_swell']:.4f}  "
              f"cells {k['units_xy']} x {k['units_xy']} x {k['units_z']}  ->  {stem}.data", flush=True)
        if args.dry_run:
            continue
        t0 = time.time()
        build('random', str(out / f'{stem}.data'), beads_per_chain=n,
              units_x=k['units_xy'], units_y=k['units_xy'], units_z=k['units_z'],
              pre_swell=k['pre_swell'], chain_placement='walk',
              bond_switch_seed=SWITCH_SEEDS.get(n, 7),
              check_isotropy=False if args.no_check else None)
        rows.append((n, k, stem, time.time() - t0))
    for n, k, stem, dt in rows:
        print(f"built N={n} in {dt:.0f} s: {stem}.data")


if __name__ == '__main__':
    main()
