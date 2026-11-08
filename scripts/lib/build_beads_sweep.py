#!/usr/bin/env python3
"""build_beads_sweep.py -- the beads-per-chain sweep of the x/y-periodic slab (2026-11-03; random
network) -- and, since 2026-11-07, of the rot45 LATTICE slab (--network rot45), the standard gel.

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

--network rot45 (2026-11-07): the diamond lattice rotated 45 degrees about y, periodic across x
and y (periodic_data_file_generator.py 'rot45'), the standard gel of the compression + permeation
programme from this date on (the random networks were too heterogeneous; a controlled lattice is
wanted, and the compression side then measures (C11 - C12)/2).  Knobs per N:
  a_eq(N)     5.28 ((N - 1)/4)^0.59 -- the 5-bead LATTICE's equilibrium lattice constant at
              P* = 1.5 (5.33 in plane / 5.26 along z on the 14000002 slab; the isolated rot45 gel
              kept 5.28) with the same good-solvent scaling.  The build starts AT a_eq (the
              14000002 slab started at 0.93 x 4.8 = 4.46 and swelled 18 %), so the aniso NPH
              equilibration only has to correct the estimate.
  units_x     round(box_xy / (a_eq/sqrt2)) BCT cells, units_y round(box_xy / a_eq) cubic cells:
              box ~48 sigma in plane for 5 beads (the equilibrated 14000002 slab's 48 x 48), ~55 sigma
              for 15 beads (BOX_XY_ROT45_N: the same atom count as the 5-bead file, see below);
  units_z     ceil(thickness / (a_eq/sqrt2)) BCT cells with thickness = 1.15 x the reference
              slab's 23 cubic cells at a_eq(5) = 139.7 sigma (the user's 2026-11-07 request:
              ~15 % taller for finer stress / volume-fraction profiles) -> 38 cells = 141.9 sigma
              for 5 beads, 18 cells = 140.8 sigma for 15 beads.
  placement   'straight' for 5 beads (the reference's), 'walk' for longer strands.
Outputs slab_support_periodic_<N>beads_tall_rho04_rot45.data.  The 5-bead rot45 build is the
reference of this family (no separate 5-bead batch: chain_length_testing_sweep.batch with
NETWORK=rot45 NBEADS=5 equilibrates it).
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
A_REF_NET = {'random': 4.82,              # per network: the 5-bead slab's equilibrated lattice constant
             'rot45': 5.28}               # the lattice at P* = 1.5 (14000002 slab; isolated rot45 gel)
# bond_switch_seed per N where the default (7) left a detectable texture: the 20-bead file with
# seed 7 failed only the interior 4th-rank anisotropy (0.068 vs random 0.039 +/- 0.007); seed 11
# passes every metric (2026-11-03).
SWITCH_SEEDS = {20: 11}
BOX_XY_REF, THICK_REF = 43.0, 111.0       # its in-plane box and gel thickness (sigma)
# rot45 lattice targets (2026-11-07): the equilibrated 14000002 slab's 48 x 48 in plane; thickness
# 1.15 x its 23 cubic cells at a_eq(5) = 1.15 x 23 x 5.28 = 139.7 sigma (rounded UP to whole cells)
BOX_XY_ROT45, THICK_ROT45 = 48.0, 1.15 * 23 * 5.28
# per-N in-plane box of the rot45 family (2026-11-07, user): the 15-bead slab is built WIDER (7 x 5 cells, ~55 sigma)
# so that it carries about the 5-bead file's atom count (~354 k) -- more strands per z-bin against the noise of
# fewer, floppier chains; the height stays the +15 % one.  A longer-strand gel at the same box would have 1/3 of
# the polymer beads and 2/3 of the atoms.  The equilibrium lattice constant is only estimated (a_eq below, the
# c* scaling); the aniso NPH equilibration sets the real in-plane size, so the converted two-piston box may come
# out 10-20 % larger (more atoms) if the gel swells past a_eq.
# 2026-11-07, later still (user): matched on POLYMER BEADS PER z-BIN instead -- the stress / volume-fraction
# profiles are binned along z, and their noise goes with the beads per bin, i.e. the polymer per unit
# height: 5-bead 125 424 / 141.9 = 884 per sigma; the 15-bead gel's polymer density at its mesh is 2.4x
# lower, so its cross-section has to be 2.4x the 5-bead one: 9 x 7 cells, ~70 x 77 sigma (~640 k atoms,
# 1.8x the 5-bead file; the two-piston runs then need ~2x the wall time or twice the nodes).  The 55-sigma
# build (7 x 5, 367 k atoms, 490 beads per sigma) is the cheaper fallback (--box-xy 55).
BOX_XY_ROT45_N = {5: 48.0, 15: 74.0}


def knobs(n_beads, box_xy=BOX_XY_REF, thick=THICK_REF, network='random'):
    """The build knobs for N beads per strand -> dict(a_eq, pre_swell, units_x, units_y, units_z, units_xy,
    placement).  random: cubic cells, box ~box_xy, thickness ~thick (rounded to the nearest cell).
    rot45: BCT cells a_eq/sqrt2 along x and z, cubic along y; the thickness is rounded UP."""
    import math
    a_eq = A_REF_NET[network] * ((n_beads - 1) / (N_REF - 1)) ** NU
    pre_swell = a_eq / (1.2 * (n_beads - 1))
    if network == 'rot45':
        cxz = a_eq / math.sqrt(2.0)
        ux, uy, uz = max(2, round(box_xy / cxz)), max(2, round(box_xy / a_eq)), max(2, math.ceil(thick / cxz - 1e-9))
        return dict(a_eq=a_eq, pre_swell=pre_swell, units_x=ux, units_y=uy, units_z=uz, units_xy=(ux, uy),
                    cell_xz=cxz, placement='straight' if n_beads <= 5 else 'walk',
                    gel=(ux * cxz, uy * a_eq, uz * cxz))
    units_xy = max(2, round(box_xy / a_eq))
    units_z = max(2, round(thick / a_eq))
    return dict(a_eq=a_eq, pre_swell=pre_swell, units_x=units_xy, units_y=units_xy, units_z=units_z, units_xy=units_xy,
                placement='walk', gel=(units_xy * a_eq, units_xy * a_eq, units_z * a_eq))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--beads', type=int, nargs='+', default=[10, 15, 20, 30])
    ap.add_argument('--network', choices=('random', 'rot45'), default='random',
                    help="random (2026-11-03 sweep) or rot45 (the lattice, 2026-11-07; defaults: --beads 5 15, box 48, thickness 139.7)")
    ap.add_argument('--outdir', default=str(Path.home() / 'Documents/lammps_data/input_data'))
    ap.add_argument('--box-xy', type=float, default=None, help='target in-plane box (sigma); default 43 (random) / 48 (rot45)')
    ap.add_argument('--thickness', type=float, default=None, help='target gel thickness (sigma); default 111 (random) / 139.7 (rot45)')
    ap.add_argument('--dry-run', action='store_true', help='print the knobs only')
    ap.add_argument('--no-check', action='store_true', help='skip the isotropy report (rot45: the density report)')
    args = ap.parse_args()
    rot = args.network == 'rot45'
    beads = args.beads if (not rot or '--beads' in sys.argv) else [5, 15]
    box_xy = args.box_xy if args.box_xy is not None else (BOX_XY_ROT45 if rot else BOX_XY_REF)
    thick = args.thickness if args.thickness is not None else (THICK_ROT45 if rot else THICK_REF)
    out = Path(args.outdir)
    rows = []
    for n in beads:
        bxy = BOX_XY_ROT45_N.get(n, box_xy) if (rot and args.box_xy is None) else box_xy
        k = knobs(n, bxy, thick, args.network)
        stem = f"slab_support_periodic_{n}beads_tall_rho04_{args.network}"
        gx, gy, gz = k['gel']
        print(f"\n##### N = {n} ({args.network}): a_eq {k['a_eq']:.3f}  pre_swell {k['pre_swell']:.4f}  "
              f"cells {k['units_x']} x {k['units_y']} x {k['units_z']}{' (BCT along x, z)' if rot else ''}  "
              f"gel {gx:.1f} x {gy:.1f} x {gz:.1f} sigma  placement {k['placement']}  ->  {stem}.data", flush=True)
        if args.dry_run:
            continue
        t0 = time.time()
        extra = {} if rot else dict(bond_switch_seed=SWITCH_SEEDS.get(n, 7))
        build(args.network, str(out / f'{stem}.data'), beads_per_chain=n,
              units_x=k['units_x'], units_y=k['units_y'], units_z=k['units_z'],
              pre_swell=k['pre_swell'], chain_placement=k['placement'],
              check_isotropy=False if args.no_check else None, **extra)
        rows.append((n, k, stem, time.time() - t0))
    for n, k, stem, dt in rows:
        print(f"built N={n} in {dt:.0f} s: {stem}.data")


if __name__ == '__main__':
    main()
