#!/usr/bin/env python3
"""
slab_shear_plates.py -- CONVERTER: equilibrated PERIODIC slab -> shear_slab input (2026-09-22)

Turns a slab_with_support periodic snapshot (network percolating in x and y,
finite in z, support + piston + bath solvent present) into the plated geometry
that simulations/shear_slab/shear_slab.lmp expects -- plates NORMAL TO z on the
slab's two free faces, driven +/-x (z = gap, x = shear, y = neutral), the same
axes as the compression decks:

  1. ISOLATE  (isolate_gel.py, unchanged): drop the support (type 4), piston
     (type 5) and the bath solvent beyond the polymer's z-extent; x and y keep
     the ORIGINAL box and periodicity (the network wraps them); z is re-boxed to
     the full polymer span + `clearance`.  Every bond is validated against the
     FENE limit under the new box.
  2. PLATE    (add_plates_to_gel.py with axis='z'): square-lattice sheets
     (type 4) `offset` outside each 1-percentile polymer surface in z, spanning
     the full periodic x-y box, harmonic-bonded to the nearest surface polymer
     bead; the z box is extended `box_buffer` beyond each plate.  Solvent beads
     within `solvent_delete` of a plate atom are REMOVED: the plates are inserted
     into the fluid, and in a periodic box the gel cannot drain from their
     excluded volume otherwise compresses it (bulk gel P 1.53 instead of 1.50 in
     the shear hold, 2026-09-22); the radius is calibrated against that pressure
     (shear_slab.lmp runs no barostat on this input -- see its Phase 1).

The result is a laterally infinite membrane sheared face-on between two plates.
(The 2026-04-28 isolated-cube input had the plates on the x-faces driven +/-z,
because plating the ends of the tall finite gel bent it; that geometry needs
the pre-2026-09-22 deck.)

Input:  final_config_slab_support_periodic_..._14000002.data  (the 2026-09-07
        piston-transparent aniso-NPH slab, same input as the two-piston converter)
Output: isolated_slab_support_periodic_..._14000002_with_plates.data (+ .info.json,
        an intermediate _isolated.data, and a slab_data_file_info.md entry)
Note:   image flags and the Velocities section of the input are dropped (the
        parsers read x y z only); bonds across x/y are minimum-image, which the
        deck's `boundary p p p` handles, and shear_slab.lmp creates velocities.

Usage:
  python slab_shear_plates.py [--input F] [--output F] [--clearance 0.5]
        [--spacing 1.5] [--offset 0.5] [--cutoff 2.5] [--box-buffer 2.0] [--no-log]
"""
import argparse
import datetime
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import isolate_gel as ig            # noqa: E402
import add_plates_to_gel as apg     # noqa: E402

DATA_DIR = HERE.parent.parent / 'lammps_data_files_local'
DEFAULT_CLEARANCE = 0.5    # sigma per face (isolate_gel's value): the film keeps the bath solvent at bath density
DEFAULT_SOLVENT_DELETE = 0.87  # sigma; CALIBRATED 2026-09-22 on the 14000002 slab: bulk gel P in the shear hold =
                               # 1.530 - 2.2e-5 * n_deleted (0 -> 1.530, 772 @0.7 sigma -> 1.513, 1943 @1.0 sigma -> 1.487);
                               # 0.87 sigma deletes ~1350 beads -> 1.500.  Re-calibrate for a different slab or plate lattice.
DEFAULT_INPUT  = str(DATA_DIR / 'final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000002.data')
DEFAULT_OUTPUT = str(DATA_DIR / 'isolated_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000002_with_plates.data')
INFO_MD = HERE.parent / 'slab_data_file_info.md'


def isolate(input_file, clearance):
    print("=" * 60); print("1. ISOLATE (isolate_gel.py)")
    atoms, bonds, box, masses = ig.parse_lammps_data(input_file)
    print(f"  Read {len(atoms)} atoms, {len(bonds)} bonds; box "
          f"{box['xhi']-box['xlo']:.3f} x {box['yhi']-box['ylo']:.3f} x {box['zhi']-box['zlo']:.3f}")
    periodic = ig.detect_percolating_dims(atoms, bonds, box)
    if periodic != {'x', 'y'}:
        raise SystemExit(f"expected a slab percolating in x and y only, got {sorted(periodic)}")
    atoms = ig.unwrap_atoms_via_bonds(atoms, bonds, box, periodic)
    atoms = ig.rotate_mobile_atoms(atoms, periodic)          # no-op for a percolating slab
    ext   = ig.find_gel_extent(atoms, clearance, 0.1, box, periodic)
    atoms = ig.remove_non_gel_atoms(atoms, ext, periodic)
    new_box = {'xlo': ext['xmin'], 'xhi': ext['xmax'], 'ylo': ext['ymin'], 'yhi': ext['ymax'],
               'zlo': ext['zmin'], 'zhi': ext['zmax']}
    ig.validate_bonds(atoms, bonds, new_box)
    return atoms, bonds, new_box, masses


def summarise(path):
    atoms, bonds, box, _ = apg.parse_lammps_data(path)
    types = {}
    for a in atoms:
        types[a['type']] = types.get(a['type'], 0) + 1
    poly_x  = [a['z'] for a in atoms if a['type'] in (1, 2)]
    plate_x = sorted({round(a['z'], 6) for a in atoms if a['type'] == 4})
    n_pb = sum(1 for b in bonds if b['type'] == 2)
    return {'n_atoms': len(atoms), 'n_bonds': len(bonds), 'n_plate_bonds': n_pb, 'types': types,
            'box': box, 'gel_z': [min(poly_x), max(poly_x)], 'gel_gap': max(poly_x) - min(poly_x),
            'plate_z': plate_x, 'lx': box['xhi'] - box['xlo'], 'ly': box['yhi'] - box['ylo'],
            'lz': box['zhi'] - box['zlo']}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input',  default=DEFAULT_INPUT)
    ap.add_argument('--output', default=DEFAULT_OUTPUT)
    ap.add_argument('--clearance',  type=float, default=DEFAULT_CLEARANCE, help='box z clearance beyond the polymer span on each free face (sigma); the film is bath solvent at bath density, so this does not set the gel pressure')
    ap.add_argument('--solvent-delete', type=float, default=DEFAULT_SOLVENT_DELETE,
                    help='remove solvent beads closer than this (sigma) to a plate atom. CALIBRATED against the bulk gel pressure in the '
                         'shear hold: with nothing deleted the plate atoms\' excluded volume in a box the gel cannot drain from read 1.53 '
                         'instead of 1.50 (2026-09-22). Re-check the P_th(bulk) printout of shear_analysis_single.ipynb after any change.')
    ap.add_argument('--spacing',    type=float, default=apg.PLATE_SPACING)
    ap.add_argument('--offset',     type=float, default=apg.PLATE_OFFSET)
    ap.add_argument('--cutoff',     type=float, default=apg.BOND_CUTOFF)
    ap.add_argument('--box-buffer', type=float, default=2.0,
                    help='box extension beyond each plate (sigma); 2.0 here vs 0.5 for the isolated cube: '
                         'the plates are position-controlled while the Phase-1a NPT box breathes, so leave room')
    ap.add_argument('--no-log', action='store_true', help='do not append to slab_data_file_info.md')
    args = ap.parse_args()

    inp, out = Path(args.input), Path(args.output)
    if not inp.exists():
        raise SystemExit(f"input not found: {inp}")
    iso = out.with_name(out.name.replace('_with_plates', '_isolated'))

    atoms, bonds, box, masses = isolate(str(inp), args.clearance)
    ig.write_lammps_data(str(iso), atoms, bonds, box, masses)
    print(f"  Wrote isolated slab -> {iso}")

    print("=" * 60); print("2. PLATE   (add_plates_to_gel.py, axis='z')")
    apg.add_plates(str(iso), str(out), spacing=args.spacing, offset=args.offset,
                   cutoff=args.cutoff, box_buffer=args.box_buffer, axis='z',
                   solvent_delete=args.solvent_delete)

    s = summarise(str(out))
    n_plate = s['types'].get(4, 0)
    info = {'input_file': str(inp), 'output_file': str(out), 'intermediate_isolated': str(iso),
            'date': datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
            'config': vars(args), 'geometry': 'plates normal to z (slab faces), driven +/-x; x,y periodic',
            **{k: s[k] for k in ('n_atoms', 'n_bonds', 'n_plate_bonds', 'types', 'box', 'gel_z', 'gel_gap', 'plate_z', 'lx', 'ly', 'lz')}}
    out.with_suffix(out.suffix + '.info.json').write_text(json.dumps(info, indent=1, default=str))
    print("=" * 60)
    print(f"SUMMARY: {s['n_atoms']} atoms ({n_plate} plate), {s['n_bonds']} bonds ({s['n_plate_bonds']} plate-polymer)")
    print(f"  box X (shear, periodic) {s['lx']:.3f}  Y {s['ly']:.3f}  Z (gap) {s['lz']:.3f}")
    print(f"  gel Z-extent {s['gel_z'][0]:.3f}..{s['gel_z'][1]:.3f}  gel_gap {s['gel_gap']:.2f}  plates at Z = {s['plate_z']}")
    print(f"  plate-to-box-face clearance: {s['plate_z'][0]-s['box']['zlo']:.2f} / {s['box']['zhi']-s['plate_z'][-1]:.2f} sigma")
    if not args.no_log:
        t = s['types']
        entry = (f"\nShear input from the PERIODIC slab  [{info['date']}, slab_shear_plates.py]:\n"
                 f"  Input: {inp.name}\n"
                 f"  Steps: isolate (support/piston/bath solvent removed; x,y periodic kept; z = polymer span + {args.clearance}) "
                 f"-> plates NORMAL TO z on the slab faces (driven +/-x in shear_slab.lmp)\n"
                 f"  Box: X (shear, periodic) {s['lx']:.3f} x Y {s['ly']:.3f} x Z (gap) {s['lz']:.3f}\n"
                 f"  Gel Z-extent {s['gel_z'][0]:.2f}..{s['gel_z'][1]:.2f} (gel_gap {s['gel_gap']:.2f}); plates at Z = "
                 f"{', '.join(f'{p:.2f}' for p in s['plate_z'])} (spacing {args.spacing}, offset {args.offset}, box_buffer {args.box_buffer}, solvent_delete {args.solvent_delete})\n"
                 f"  Crosslinks: {t.get(1,0)}   Chain beads: {t.get(2,0)}   Solvent: {t.get(3,0)}   Plates: {n_plate}   "
                 f"plate-polymer bonds: {s['n_plate_bonds']} (cutoff {args.cutoff})\n"
                 f"  Total atoms: {s['n_atoms']}   Total bonds: {s['n_bonds']}   Output: {out.name}\n")
        with open(INFO_MD, 'a') as f:
            f.write(entry)
        print(f"  Logged to {INFO_MD.name}")


if __name__ == '__main__':
    main()
