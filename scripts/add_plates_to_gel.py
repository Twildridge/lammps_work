#!/usr/bin/env python3
"""
add_plates_to_gel.py

Adds rigid graphene-like shear plates to the top and bottom of an isolated gel
data file. The plates drive shear in shear_slab.lmp via fix move linear.

Design
------
- Plate atoms are atom type 4 (mass 1.0), arranged on a square lattice
  spanning the full periodic box in xy, at z = z_surface ± PLATE_OFFSET.
- Each plate atom is bonded (harmonic, bond type 2) to the single nearest
  backbone polymer atom (type 1 or 2) within BOND_CUTOFF. If no polymer atom
  is found within the cutoff, the plate atom still exists but has no bond
  (relies on LJ for contact).
- Existing polymer FENE bonds (type 1) are unchanged.
- Output data file has 4 atom types and 2 bond types.

Bond styles needed in shear_slab.lmp:
    bond_style hybrid fene harmonic
    bond_coeff 1 fene 30.0 1.5 1.0 1.0        # polymer-polymer (unchanged)
    bond_coeff 2 harmonic 30.0 1.0             # plate-polymer

Pair coeff for plate (type 4) — add to shear_slab.lmp:
    pair_coeff 4 4 1.0 1.0 1.122               # plate-plate  (WCA)
    pair_coeff 1 4 1.0 1.0 1.122               # polymer-plate (WCA)
    pair_coeff 2 4 1.0 1.0 1.122               # crosslinker-plate (WCA)
    pair_coeff 3 4 1.0 1.0 1.122               # solvent-plate (WCA)

Usage
-----
    python add_plates_to_gel.py \\
        --input  /path/to/isolated_gel.data \\
        --output /path/to/gel_with_plates.data \\
        [--spacing 1.5] [--offset 0.5] [--cutoff 2.5] [--surface-depth 4.0]
"""

import numpy as np
import argparse
from scipy.spatial import KDTree

# ── Defaults ─────────────────────────────────────────────────────────────────
PLATE_TYPE      = 4      # new atom type for plate beads
PLATE_BOND_TYPE = 2      # harmonic bond type (type 1 = FENE, kept)
PLATE_SPACING   = 1.5    # square lattice constant (σ)
PLATE_OFFSET    = 0.5    # distance from polymer surface (σ)
BOND_CUTOFF     = 2.5    # max 3D distance for plate→polymer bond creation (σ)
SURFACE_DEPTH   = 4.0    # z-depth of surface layer to search for polymer atoms (σ)
SURFACE_PCT     = 1.0    # percentile to robustly define surface z (ignore outliers)

POLYMER_TYPES   = {1, 2} # atom types that can be bonded to plate
BOND_POLY_TYPES = {1, 2} # subset bonded to plate (backbone + crosslinker)


# ── I/O helpers ──────────────────────────────────────────────────────────────

def parse_lammps_data(filename):
    """Return atoms, bonds, box_bounds, masses from a LAMMPS molecular data file."""
    atoms, bonds, masses = [], [], {}
    box = {}

    with open(filename) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if 'xlo xhi' in line:
            p = line.split(); box['xlo'], box['xhi'] = float(p[0]), float(p[1])
        elif 'ylo yhi' in line:
            p = line.split(); box['ylo'], box['yhi'] = float(p[0]), float(p[1])
        elif 'zlo zhi' in line:
            p = line.split(); box['zlo'], box['zhi'] = float(p[0]), float(p[1])

        elif line == 'Masses':
            i += 2
            while i < len(lines) and lines[i].strip() and \
                    not lines[i].strip()[0].isalpha():
                p = lines[i].split()
                if len(p) >= 2:
                    masses[int(p[0])] = float(p[1])
                i += 1
            continue

        elif line.startswith('Atoms'):
            i += 2
            while i < len(lines) and lines[i].strip() and \
                    not lines[i].strip()[0].isalpha():
                p = lines[i].split()
                if len(p) >= 6:
                    atoms.append({
                        'id':   int(p[0]),
                        'mol':  int(p[1]),
                        'type': int(p[2]),
                        'x':    float(p[3]),
                        'y':    float(p[4]),
                        'z':    float(p[5]),
                    })
                i += 1
            continue

        elif line.startswith('Bonds'):
            i += 2
            while i < len(lines) and lines[i].strip() and \
                    not lines[i].strip()[0].isalpha():
                p = lines[i].split()
                if len(p) >= 4:
                    bonds.append({
                        'id':    int(p[0]),
                        'type':  int(p[1]),
                        'atom1': int(p[2]),
                        'atom2': int(p[3]),
                    })
                i += 1
            continue

        i += 1

    return atoms, bonds, box, masses


def write_lammps_data(filename, atoms, bonds, box, masses):
    """Write LAMMPS molecular data file. atom IDs are renumbered 1..N."""
    old2new = {a['id']: i + 1 for i, a in enumerate(atoms)}

    valid_bonds = []
    for i, b in enumerate(bonds):
        a1 = old2new.get(b['atom1'])
        a2 = old2new.get(b['atom2'])
        if a1 and a2:
            valid_bonds.append({'id': i + 1, 'type': b['type'],
                                 'atom1': a1, 'atom2': a2})

    n_types = max(a['type'] for a in atoms)
    n_btypes = max((b['type'] for b in valid_bonds), default=1)

    with open(filename, 'w') as f:
        f.write("LAMMPS data file - isolated gel with shear plates\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(valid_bonds)} bonds\n\n")
        f.write(f"{n_types} atom types\n")
        f.write(f"{n_btypes} bond types\n\n")
        f.write(f"{box['xlo']:.6f} {box['xhi']:.6f} xlo xhi\n")
        f.write(f"{box['ylo']:.6f} {box['yhi']:.6f} ylo yhi\n")
        f.write(f"{box['zlo']:.6f} {box['zhi']:.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        for t in range(1, n_types + 1):
            f.write(f"{t} {masses.get(t, 1.0):.4f}\n")

        f.write("\nAtoms\n\n")
        for i, a in enumerate(atoms, 1):
            f.write(f"{i} {a['mol']} {a['type']} "
                    f"{a['x']:.6f} {a['y']:.6f} {a['z']:.6f}\n")

        if valid_bonds:
            f.write("\nBonds\n\n")
            for b in valid_bonds:
                f.write(f"{b['id']} {b['type']} {b['atom1']} {b['atom2']}\n")


# ── Plate generation ──────────────────────────────────────────────────────────

def surface_z(polymer_z: np.ndarray, top: bool, pct: float) -> float:
    """Robust surface z using percentile (avoids outlier atoms)."""
    return float(np.percentile(polymer_z, 100 - pct if top else pct))


def make_plate(box, z_plane, spacing, mol_id):
    """
    Generate plate atoms on a square lattice in the xy-plane at z = z_plane.
    Lattice spans [xlo, xhi) × [ylo, yhi) so atoms tile perfectly under PBC.
    Returns list of atom dicts (type = PLATE_TYPE).
    """
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']

    nx = int(np.floor(lx / spacing))
    ny = int(np.floor(ly / spacing))

    # Slight inset from edges so atoms sit on a regular grid centred in box
    dx = lx / nx
    dy = ly / ny
    x0 = box['xlo'] + dx / 2.0
    y0 = box['ylo'] + dy / 2.0

    plate = []
    for ix in range(nx):
        for iy in range(ny):
            plate.append({
                'id':   None,          # assigned later
                'mol':  mol_id,
                'type': PLATE_TYPE,
                'x':    x0 + ix * dx,
                'y':    y0 + iy * dy,
                'z':    z_plane,
            })
    return plate


def create_plate_bonds(plate_atoms, poly_atoms, poly_id_offset,
                       bond_id_start, cutoff):
    """
    For each plate atom, find the nearest polymer atom within `cutoff` and
    create a harmonic bond (type PLATE_BOND_TYPE).

    poly_id_offset: index in the final atom list where polymer starts (0-based).
                    plate atoms are appended AFTER polymer, so their final IDs
                    come later — we store references by list position and
                    remap in write_lammps_data via old2new.

    Returns list of bond dicts using the ORIGINAL atom id field.
    """
    # Build KDTree of surface polymer atom positions
    poly_xyz = np.array([[a['x'], a['y'], a['z']] for a in poly_atoms])
    tree = KDTree(poly_xyz)

    bonds = []
    bid = bond_id_start
    bonded_plate = 0

    for pa in plate_atoms:
        q = np.array([pa['x'], pa['y'], pa['z']])
        dist, idx = tree.query(q, k=1)
        if dist <= cutoff:
            bonds.append({
                'id':    bid,
                'type':  PLATE_BOND_TYPE,
                'atom1': pa['id'],        # plate atom id
                'atom2': poly_atoms[idx]['id'],  # polymer atom id
            })
            bid += 1
            bonded_plate += 1

    return bonds, bonded_plate


# ── Main ─────────────────────────────────────────────────────────────────────

def add_plates(input_file, output_file,
               spacing=PLATE_SPACING, offset=PLATE_OFFSET,
               cutoff=BOND_CUTOFF, surface_depth=SURFACE_DEPTH,
               surface_pct=SURFACE_PCT):

    print("=" * 60)
    print(f"add_plates_to_gel.py")
    print(f"  Input  : {input_file}")
    print(f"  Output : {output_file}")
    print(f"  Lattice spacing : {spacing} σ")
    print(f"  Plate offset    : {offset} σ above/below gel surface")
    print(f"  Bond cutoff     : {cutoff} σ")
    print(f"  Surface depth   : {surface_depth} σ")
    print("=" * 60)

    atoms, bonds, box, masses = parse_lammps_data(input_file)
    print(f"Read {len(atoms)} atoms, {len(bonds)} bonds")

    # ── Identify polymer atoms ────────────────────────────────────────────────
    poly_all  = [a for a in atoms if a['type'] in POLYMER_TYPES]
    poly_z    = np.array([a['z'] for a in poly_all])

    z_top_surf = surface_z(poly_z, top=True,  pct=surface_pct)
    z_bot_surf = surface_z(poly_z, top=False, pct=surface_pct)
    gel_thick  = z_top_surf - z_bot_surf

    print(f"\nGel polymer z-extent: {z_bot_surf:.3f} → {z_top_surf:.3f}  "
          f"(thickness {gel_thick:.2f} σ)")

    # ── Surface polymer atoms (for bonding) ───────────────────────────────────
    surf_top_poly = [a for a in poly_all
                     if a['z'] >= z_top_surf - surface_depth
                     and a['type'] in BOND_POLY_TYPES]
    surf_bot_poly = [a for a in poly_all
                     if a['z'] <= z_bot_surf + surface_depth
                     and a['type'] in BOND_POLY_TYPES]
    print(f"Surface polymer atoms available for bonding: "
          f"{len(surf_top_poly)} (top), {len(surf_bot_poly)} (bottom)")

    # ── Generate plate atoms ──────────────────────────────────────────────────
    z_top_plate = z_top_surf + offset
    z_bot_plate = z_bot_surf - offset

    # Give plates their own molecule IDs (above any existing mol ID)
    max_mol = max(a['mol'] for a in atoms)
    plate_top = make_plate(box, z_top_plate, spacing, mol_id=max_mol + 1)
    plate_bot = make_plate(box, z_bot_plate, spacing, mol_id=max_mol + 2)

    print(f"\nTop plate   : z = {z_top_plate:.3f}  |  {len(plate_top)} atoms")
    print(f"Bottom plate: z = {z_bot_plate:.3f}  |  {len(plate_bot)} atoms")
    print(f"Plate separation (surface-to-surface): {z_top_plate - z_bot_plate:.3f} σ")

    # ── Assign IDs ────────────────────────────────────────────────────────────
    # Original atoms keep their IDs; plate atoms get new sequential IDs.
    next_id = max(a['id'] for a in atoms) + 1
    for a in plate_top + plate_bot:
        a['id'] = next_id
        next_id += 1

    # Add mass for plate type
    masses[PLATE_TYPE] = 1.0

    # ── Create harmonic bonds ─────────────────────────────────────────────────
    next_bid = max((b['id'] for b in bonds), default=0) + 1

    top_bonds, n_top_bonded = create_plate_bonds(
        plate_top, surf_top_poly, 0, next_bid, cutoff)
    next_bid += len(top_bonds)

    bot_bonds, n_bot_bonded = create_plate_bonds(
        plate_bot, surf_bot_poly, 0, next_bid, cutoff)

    print(f"\nHarmonic bonds created:")
    print(f"  Top plate → polymer: {len(top_bonds)} bonds "
          f"({n_top_bonded}/{len(plate_top)} plate atoms bonded, "
          f"{100*n_top_bonded/len(plate_top):.1f}%)")
    print(f"  Bot plate → polymer: {len(bot_bonds)} bonds "
          f"({n_bot_bonded}/{len(plate_bot)} plate atoms bonded, "
          f"{100*n_bot_bonded/len(plate_bot):.1f}%)")
    print(f"  Total new bonds: {len(top_bonds) + len(bot_bonds)}")

    # ── Extend box z to contain plates ───────────────────────────────────────
    # Add 0.5σ buffer beyond the plate so PBC images don't overlap.
    new_box = dict(box)
    new_box['zlo'] = min(box['zlo'], z_bot_plate - 0.5)
    new_box['zhi'] = max(box['zhi'], z_top_plate + 0.5)

    # ── Assemble and write ────────────────────────────────────────────────────
    all_atoms = atoms + plate_top + plate_bot
    all_bonds = bonds + top_bonds + bot_bonds

    write_lammps_data(output_file, all_atoms, all_bonds, new_box, masses)

    n_plate = len(plate_top) + len(plate_bot)
    print(f"\nWrote {len(all_atoms)} atoms ({n_plate} plate), "
          f"{len(all_bonds)} bonds to:\n  {output_file}")
    print("=" * 60)

    # ── Summary for shear_slab.lmp ────────────────────────────────────────────
    print("\n── shear_slab.lmp hints ─────────────────────────────────────")
    print(f"  gel_thick (surface-to-surface) ≈ {gel_thick:.2f} σ")
    print(f"  plate separation               ≈ {z_top_plate - z_bot_plate:.2f} σ")
    print(f"  For target_strain=0.10, nsteps=285000, dt=0.005:")
    vshear = 0.10 * gel_thick / (2 * 285000 * 0.005)
    print(f"    vshear = {vshear:.6f} σ/τ")
    print("  bond_style hybrid fene harmonic")
    print("  bond_coeff 1 fene 30.0 1.5 1.0 1.0   # polymer-polymer")
    print("  bond_coeff 2 harmonic 30.0 1.0        # plate-polymer")
    print("  pair_coeff 4 4 1.0 1.0 1.122          # plate-plate (WCA)")
    print("  pair_coeff 1 4 1.0 1.0 1.122          # polymer-plate (WCA)")
    print("  pair_coeff 2 4 1.0 1.0 1.122          # crosslinker-plate (WCA)")
    print("  pair_coeff 3 4 1.0 1.0 1.122          # solvent-plate (WCA)")
    print("─" * 60)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input',  required=True, help='Input isolated gel .data file')
    p.add_argument('--output', required=True, help='Output .data file with plates')
    p.add_argument('--spacing',       type=float, default=PLATE_SPACING,
                   help=f'Plate lattice constant in σ (default {PLATE_SPACING})')
    p.add_argument('--offset',        type=float, default=PLATE_OFFSET,
                   help=f'Plate z-offset from gel surface in σ (default {PLATE_OFFSET})')
    p.add_argument('--cutoff',        type=float, default=BOND_CUTOFF,
                   help=f'Max 3D distance for plate→polymer bond (default {BOND_CUTOFF})')
    p.add_argument('--surface-depth', type=float, default=SURFACE_DEPTH,
                   help=f'Depth of surface layer searched for polymer (default {SURFACE_DEPTH})')
    args = p.parse_args()

    add_plates(
        input_file    = args.input,
        output_file   = args.output,
        spacing       = args.spacing,
        offset        = args.offset,
        cutoff        = args.cutoff,
        surface_depth = args.surface_depth,
    )


if __name__ == '__main__':
    main()
