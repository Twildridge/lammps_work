"""
add_plates_to_gel.py

(Extracted verbatim from add_plates_to_gel.ipynb on 2026-09-22 so that
slab_shear_plates.py -- the periodic-slab -> shear input converter -- can reuse
it; the notebook now imports this module.  Additions: the `box_buffer`
parameter (default 0.5 = the old hard-coded value) and `axis` (default 'x' =
the isolated-cube geometry; 'z' plates the faces of the periodic slab).)

Adds rigid graphene-like shear plates to the left and right sides (x-faces)
of an isolated gel data file. The plates drive shear in shear_slab.lmp via
fix move linear.

Design
------
- The gel is tall in z (long axis) and shorter in x and y. Plates are placed
  on the x-faces (left/right) so the shear gap is in x and the long z-axis
  runs along the plate faces. This avoids the bending stresses that arise
  when plates are placed on the top/bottom of a high-aspect-ratio gel.
- Plate atoms are atom type 4 (mass 1.0), arranged on a square lattice
  spanning the full periodic box in yz, at x = x_surface ± PLATE_OFFSET.
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
SURFACE_DEPTH   = 4.0    # x-depth of surface layer to search for polymer atoms (σ)
SURFACE_PCT     = 1.0    # percentile to robustly define surface x (ignore outliers)
BOX_BUFFER      = 0.5    # box extension beyond each plate (σ); the isolated-cube files were built with 0.5

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

def surface_x(polymer_x: np.ndarray, right: bool, pct: float) -> float:
    """Robust surface coordinate using percentile (avoids outlier atoms)."""
    return float(np.percentile(polymer_x, 100 - pct if right else pct))


OTHER_AXES = {'x': ('y', 'z'), 'y': ('x', 'z'), 'z': ('x', 'y')}


def make_plate(box, x_plane, spacing, mol_id, axis='x'):
    """
    Generate plate atoms on a square lattice in the plane normal to `axis` at
    coordinate x_plane.  The lattice spans the full box in the other two
    directions so atoms tile perfectly under PBC.
    Returns list of atom dicts (type = PLATE_TYPE).
    """
    a1, a2 = OTHER_AXES[axis]
    l1 = box[f'{a1}hi'] - box[f'{a1}lo']
    l2 = box[f'{a2}hi'] - box[f'{a2}lo']

    n1 = int(np.floor(l1 / spacing))
    n2 = int(np.floor(l2 / spacing))

    # Slight inset from edges so atoms sit on a regular grid centred in box
    d1 = l1 / n1
    d2 = l2 / n2
    c1_0 = box[f'{a1}lo'] + d1 / 2.0
    c2_0 = box[f'{a2}lo'] + d2 / 2.0

    plate = []
    for i1 in range(n1):
        for i2 in range(n2):
            plate.append({
                'id':   None,          # assigned later
                'mol':  mol_id,
                'type': PLATE_TYPE,
                axis:   x_plane,
                a1:     c1_0 + i1 * d1,
                a2:     c2_0 + i2 * d2,
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
               surface_pct=SURFACE_PCT, box_buffer=BOX_BUFFER, axis='x',
               solvent_delete=0.0):

    print("=" * 60)
    print(f"add_plates_to_gel.py")
    print(f"  Input  : {input_file}")
    print(f"  Output : {output_file}")
    print(f"  Lattice spacing : {spacing} σ")
    print(f"  Plate offset    : {offset} σ left/right of gel surface")
    print(f"  Bond cutoff     : {cutoff} σ")
    print(f"  Surface depth   : {surface_depth} σ")
    print(f"  Plate normal    : {axis}  (plates span the {''.join(OTHER_AXES[axis])} box faces)")
    print(f"  Solvent delete  : {solvent_delete} σ  (solvent beads closer than this to a plate atom are removed)")
    print("=" * 60)

    atoms, bonds, box, masses = parse_lammps_data(input_file)
    print(f"Read {len(atoms)} atoms, {len(bonds)} bonds")

    # ── Identify polymer atoms ────────────────────────────────────────────────
    poly_all  = [a for a in atoms if a['type'] in POLYMER_TYPES]
    poly_x    = np.array([a[axis] for a in poly_all])

    x_right_surf = surface_x(poly_x, right=True,  pct=surface_pct)
    x_left_surf  = surface_x(poly_x, right=False, pct=surface_pct)
    gel_gap      = x_right_surf - x_left_surf   # gap direction (axis)

    print(f"\nGel polymer {axis}-extent: {x_left_surf:.3f} → {x_right_surf:.3f}  "
          f"(gap {gel_gap:.2f} σ)")

    # ── Surface polymer atoms (for bonding) ───────────────────────────────────
    surf_right_poly = [a for a in poly_all
                       if a[axis] >= x_right_surf - surface_depth
                       and a['type'] in BOND_POLY_TYPES]
    surf_left_poly  = [a for a in poly_all
                       if a[axis] <= x_left_surf + surface_depth
                       and a['type'] in BOND_POLY_TYPES]
    print(f"Surface polymer atoms available for bonding: "
          f"{len(surf_right_poly)} (hi/+{axis}), {len(surf_left_poly)} (lo/-{axis})")

    # ── Generate plate atoms ──────────────────────────────────────────────────
    x_right_plate = x_right_surf + offset
    x_left_plate  = x_left_surf  - offset

    # Give plates their own molecule IDs (above any existing mol ID)
    max_mol = max(a['mol'] for a in atoms)
    plate_right = make_plate(box, x_right_plate, spacing, mol_id=max_mol + 1, axis=axis)
    plate_left  = make_plate(box, x_left_plate,  spacing, mol_id=max_mol + 2, axis=axis)

    print(f"\nHi plate (+{axis}): {axis} = {x_right_plate:.3f}  |  {len(plate_right)} atoms")
    print(f"Lo plate (-{axis}): {axis} = {x_left_plate:.3f}  |  {len(plate_left)} atoms")
    print(f"Plate separation (surface-to-surface): {x_right_plate - x_left_plate:.3f} σ")

    # ── Assign IDs ────────────────────────────────────────────────────────────
    # Original atoms keep their IDs; plate atoms get new sequential IDs.
    next_id = max(a['id'] for a in atoms) + 1
    for a in plate_right + plate_left:
        a['id'] = next_id
        next_id += 1

    # Add mass for plate type
    masses[PLATE_TYPE] = 1.0

    # ── Create harmonic bonds ─────────────────────────────────────────────────
    next_bid = max((b['id'] for b in bonds), default=0) + 1

    right_bonds, n_right_bonded = create_plate_bonds(
        plate_right, surf_right_poly, 0, next_bid, cutoff)
    next_bid += len(right_bonds)

    left_bonds, n_left_bonded = create_plate_bonds(
        plate_left, surf_left_poly, 0, next_bid, cutoff)

    print(f"\nHarmonic bonds created:")
    print(f"  Hi plate → polymer: {len(right_bonds)} bonds "
          f"({n_right_bonded}/{len(plate_right)} plate atoms bonded, "
          f"{100*n_right_bonded/len(plate_right):.1f}%)")
    print(f"  Lo plate → polymer: {len(left_bonds)} bonds "
          f"({n_left_bonded}/{len(plate_left)} plate atoms bonded, "
          f"{100*n_left_bonded/len(plate_left):.1f}%)")
    print(f"  Total new bonds: {len(right_bonds) + len(left_bonds)}")

    # ── Extend box x to contain plates ───────────────────────────────────────
    # box_buffer beyond each plate: keeps the plate's periodic image clear of
    # the other plate, and leaves room for the Phase-1a NPT box to shrink
    # without a box face crossing a (position-controlled) plate.
    new_box = dict(box)
    new_box[f'{axis}lo'] = min(box[f'{axis}lo'], x_left_plate  - box_buffer)
    new_box[f'{axis}hi'] = max(box[f'{axis}hi'], x_right_plate + box_buffer)

    # ── Delete solvent overlapping the plates ─────────────────────────────────
    # The plates are inserted INTO the fluid (at the 1-percentile surface).  In a
    # fixed periodic box the gel cannot drain from, the volume the plate atoms
    # exclude compresses the fluid (+0.03 in the bulk gel pressure on the
    # 14000002 slab with nothing deleted).  Removing the solvent beads whose
    # centres lie within `solvent_delete` of a plate atom (periodic images in the
    # in-plane directions included) gives that volume back; the radius is
    # calibrated by the bulk pressure in the shear hold (see slab_shear_plates.py).
    n_del = 0
    if solvent_delete > 0.0:
        a1, a2 = OTHER_AXES[axis]
        L1 = box[f'{a1}hi'] - box[f'{a1}lo']; L2 = box[f'{a2}hi'] - box[f'{a2}lo']
        pts = []
        for pa in plate_right + plate_left:
            for i1 in (-1, 0, 1):
                for i2 in (-1, 0, 1):
                    q = {axis: pa[axis], a1: pa[a1] + i1 * L1, a2: pa[a2] + i2 * L2}
                    pts.append([q['x'], q['y'], q['z']])
        ptree = KDTree(np.array(pts))
        keep = []
        for a in atoms:
            if a['type'] == 3 and ptree.query([a['x'], a['y'], a['z']], k=1)[0] < solvent_delete:
                n_del += 1
                continue
            keep.append(a)
        atoms = keep
        print(f"\nDeleted {n_del} solvent beads within {solvent_delete} σ of a plate atom "
              f"({n_del / (len(plate_right) + len(plate_left)):.2f} per plate atom)")

    # ── Assemble and write ────────────────────────────────────────────────────
    all_atoms = atoms + plate_right + plate_left
    all_bonds = bonds + right_bonds + left_bonds

    write_lammps_data(output_file, all_atoms, all_bonds, new_box, masses)

    n_plate = len(plate_right) + len(plate_left)
    print(f"\nWrote {len(all_atoms)} atoms ({n_plate} plate), "
          f"{len(all_bonds)} bonds to:\n  {output_file}")
    print("=" * 60)

    # ── Summary for shear_slab.lmp ────────────────────────────────────────────
    print("\n── shear_slab.lmp hints ─────────────────────────────────────")
    print(f"  gel_gap  ({axis} surface-to-surface) ≈ {gel_gap:.2f} σ")
    print(f"  plate separation ({axis})            ≈ {x_right_plate - x_left_plate:.2f} σ")
    print(f"  For target_strain=0.10, nsteps=285000, dt=0.005:")
    vshear = 0.10 * gel_gap / (2 * 285000 * 0.005)
    print(f"    vshear = {vshear:.6f} σ/τ  (plates move in ±{'z' if axis == 'x' else 'x'})")
    print("  bond_style hybrid fene harmonic")
    print("  bond_coeff 1 fene 30.0 1.5 1.0 1.0   # polymer-polymer")
    print("  bond_coeff 2 harmonic 30.0 1.0        # plate-polymer")
    print("  pair_coeff 4 4 1.0 1.0 1.122          # plate-plate (WCA)")
    print("  pair_coeff 1 4 1.0 1.0 1.122          # polymer-plate (WCA)")
    print("  pair_coeff 2 4 1.0 1.0 1.122          # crosslinker-plate (WCA)")
    print("  pair_coeff 3 4 1.0 1.0 1.122          # solvent-plate (WCA)")
    print("─" * 60)


# ── CLI ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Add shear plates to the x-faces of an isolated gel data file.")
    parser.add_argument('--input',  required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--spacing', type=float, default=PLATE_SPACING)
    parser.add_argument('--offset',  type=float, default=PLATE_OFFSET)
    parser.add_argument('--cutoff',  type=float, default=BOND_CUTOFF)
    parser.add_argument('--surface-depth', type=float, default=SURFACE_DEPTH)
    parser.add_argument('--box-buffer', type=float, default=BOX_BUFFER)
    parser.add_argument('--axis', default='x', choices=('x', 'y', 'z'), help="plate normal ('x' isolated cube; 'z' periodic slab)")
    parser.add_argument('--solvent-delete', type=float, default=0.0, help='remove solvent beads closer than this (sigma) to a plate atom')
    a = parser.parse_args()
    add_plates(a.input, a.output, a.spacing, a.offset, a.cutoff, a.surface_depth, box_buffer=a.box_buffer, axis=a.axis, solvent_delete=a.solvent_delete)
