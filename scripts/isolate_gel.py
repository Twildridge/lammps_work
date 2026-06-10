#!/usr/bin/env python3
"""
isolate_gel.py
==============
CLI wrapper for isolating an equilibrated gel slab from a slab_with_support
LAMMPS data file.

Removes support (type 4) and piston (type 5) atoms, trims any solvent/polymer
outside a tight bounding box around the polymer network, and writes a new data
file containing only gel-internal atoms (types 1=crosslink, 2=chain, 3=solvent).

Usage
-----
  python3 isolate_gel.py --input <final_config.data> --output <isolated.data>
                         [--clearance 0.2] [--percentile 0.1]

Output atom types
-----------------
  1  Crosslink   (unchanged)
  2  Chain bead  (unchanged)
  3  Solvent     (gel-internal only; bath solvent trimmed)
"""

import argparse
import numpy as np

try:
    from scipy.spatial import ConvexHull
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False
    print("WARNING: scipy not available — skipping slab rotation (bounding box will be axis-aligned)")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
POLYMER_TYPES = {1, 2}
SOLVENT_TYPES = {3}
SUPPORT_TYPES = {4}
PISTON_TYPES  = {5}
ROTATE_TYPES  = POLYMER_TYPES | SOLVENT_TYPES
BOX_CLEARANCE = 0.2


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def find_min_bounding_rect_angle(hull_points):
    n = len(hull_points)
    min_area  = float('inf')
    best_angle = 0.0
    for i in range(n):
        edge  = hull_points[(i + 1) % n] - hull_points[i]
        angle = np.arctan2(edge[1], edge[0])
        c, s  = np.cos(-angle), np.sin(-angle)
        rot   = np.column_stack([
            c * hull_points[:, 0] - s * hull_points[:, 1],
            s * hull_points[:, 0] + c * hull_points[:, 1],
        ])
        area = (rot[:, 0].max() - rot[:, 0].min()) * \
               (rot[:, 1].max() - rot[:, 1].min())
        if area < min_area:
            min_area   = area
            best_angle = angle
    return best_angle


def rotate_mobile_atoms(atoms):
    if not _SCIPY_AVAILABLE:
        print("  Skipping rotation (scipy unavailable) — using axis-aligned bounding box")
        return atoms
    poly = [a for a in atoms if a['type'] in POLYMER_TYPES]
    if not poly:
        raise ValueError("No polymer atoms found in input file.")
    xy  = np.array([[a['x'], a['y']] for a in poly])
    com = xy.mean(axis=0)
    hull  = ConvexHull(xy)
    angle = find_min_bounding_rect_angle(xy[hull.vertices])
    print(f"  Slab rotation angle (MABR): {np.degrees(angle):.2f}°")
    c, s = np.cos(-angle), np.sin(-angle)
    for a in atoms:
        if a['type'] in ROTATE_TYPES:
            x = a['x'] - com[0];  y = a['y'] - com[1]
            a['x'] = c*x - s*y + com[0]
            a['y'] = s*x + c*y + com[1]
    return atoms


def find_gel_extent(atoms, clearance, percentile=0.1):
    poly = [a for a in atoms if a['type'] in POLYMER_TYPES]
    xs = np.array([a['x'] for a in poly])
    ys = np.array([a['y'] for a in poly])
    zs = np.array([a['z'] for a in poly])
    xmin, xmax = np.percentile(xs, percentile),     np.percentile(xs, 100-percentile)
    ymin, ymax = np.percentile(ys, percentile),     np.percentile(ys, 100-percentile)
    zmin, zmax = np.percentile(zs, percentile),     np.percentile(zs, 100-percentile)
    print(f"  Polymer extent  x: {xmin:.2f}–{xmax:.2f}  "
          f"y: {ymin:.2f}–{ymax:.2f}  z: {zmin:.2f}–{zmax:.2f}")
    return dict(xmin=xmin-clearance, xmax=xmax+clearance,
                ymin=ymin-clearance, ymax=ymax+clearance,
                zmin=zmin-clearance, zmax=zmax+clearance)


def remove_non_gel_atoms(atoms, ext):
    kept = []
    n_walls = n_outside = 0
    for a in atoms:
        if a['type'] in SUPPORT_TYPES | PISTON_TYPES:
            n_walls += 1;  continue
        if (a['x'] < ext['xmin'] or a['x'] > ext['xmax'] or
            a['y'] < ext['ymin'] or a['y'] > ext['ymax'] or
            a['z'] < ext['zmin'] or a['z'] > ext['zmax']):
            n_outside += 1;  continue
        kept.append(a)
    print(f"  Removed {n_walls} support/piston atoms, "
          f"{n_outside} atoms outside gel bounds")
    return kept


# ---------------------------------------------------------------------------
# LAMMPS I/O
# ---------------------------------------------------------------------------
def parse_lammps_data(path):
    atoms = [];  bonds = [];  box = {};  masses = {}
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if 'xlo xhi' in line:
            p = line.split();  box['xlo'] = float(p[0]);  box['xhi'] = float(p[1])
        elif 'ylo yhi' in line:
            p = line.split();  box['ylo'] = float(p[0]);  box['yhi'] = float(p[1])
        elif 'zlo zhi' in line:
            p = line.split();  box['zlo'] = float(p[0]);  box['zhi'] = float(p[1])
        elif line == 'Masses':
            i += 2
            while i < len(lines) and lines[i].strip() \
                    and not lines[i].strip()[0].isalpha():
                p = lines[i].split()
                if len(p) >= 2:
                    try: masses[int(p[0])] = float(p[1])
                    except ValueError: pass
                i += 1
            continue
        elif line.startswith('Atoms'):
            i += 2
            while i < len(lines) and lines[i].strip() \
                    and not lines[i].strip()[0].isalpha():
                p = lines[i].split()
                if len(p) >= 6:
                    try:
                        atoms.append({'id': int(p[0]), 'mol': int(p[1]),
                                      'type': int(p[2]),
                                      'x': float(p[3]), 'y': float(p[4]),
                                      'z': float(p[5])})
                    except ValueError: pass
                i += 1
            continue
        elif line.startswith('Bonds'):
            i += 2
            while i < len(lines) and lines[i].strip() \
                    and not lines[i].strip()[0].isalpha():
                p = lines[i].split()
                if len(p) >= 4:
                    try:
                        bonds.append({'id': int(p[0]), 'type': int(p[1]),
                                      'atom1': int(p[2]), 'atom2': int(p[3])})
                    except ValueError: pass
                i += 1
            continue
        i += 1
    return atoms, bonds, box, masses


def write_lammps_data(path, atoms, bonds, box, masses):
    old2new = {a['id']: i+1 for i, a in enumerate(atoms)}
    valid_bonds = [
        {'id': i+1, 'type': b['type'],
         'atom1': old2new[b['atom1']], 'atom2': old2new[b['atom2']]}
        for i, b in enumerate(bonds)
        if b['atom1'] in old2new and b['atom2'] in old2new
    ]
    types_present = sorted({a['type'] for a in atoms})
    n_atom_types  = max(types_present) if types_present else 3

    with open(path, 'w') as f:
        f.write("LAMMPS data file — isolated gel (support/piston/bath solvent removed)\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(valid_bonds)} bonds\n\n")
        f.write(f"{n_atom_types} atom types\n")
        f.write("1 bond types\n\n")
        f.write(f"{box['xlo']:.6f} {box['xhi']:.6f} xlo xhi\n")
        f.write(f"{box['ylo']:.6f} {box['yhi']:.6f} ylo yhi\n")
        f.write(f"{box['zlo']:.6f} {box['zhi']:.6f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for t in range(1, n_atom_types + 1):
            f.write(f"{t} {masses.get(t, 1.0):.1f}\n")
        f.write("\nAtoms\n\n")
        for i, a in enumerate(atoms, 1):
            f.write(f"{i} {a['mol']} {a['type']} "
                    f"{a['x']:.6f} {a['y']:.6f} {a['z']:.6f}\n")
        if valid_bonds:
            f.write("\nBonds\n\n")
            for b in valid_bonds:
                f.write(f"{b['id']} {b['type']} {b['atom1']} {b['atom2']}\n")


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def isolate_gel(input_file, output_file, clearance=BOX_CLEARANCE, percentile=0.1):
    print(f"Input:     {input_file}")
    print(f"Output:    {output_file}")
    print(f"Clearance: {clearance}  |  Percentile: {percentile}")
    atoms, bonds, box, masses = parse_lammps_data(input_file)
    print(f"Read {len(atoms)} atoms, {len(bonds)} bonds")
    atoms = rotate_mobile_atoms(atoms)
    ext   = find_gel_extent(atoms, clearance, percentile)
    atoms = remove_non_gel_atoms(atoms, ext)
    new_box = {k: ext[k] for k in ('xmin','xmax','ymin','ymax','zmin','zmax')}
    new_box = {'xlo': ext['xmin'], 'xhi': ext['xmax'],
               'ylo': ext['ymin'], 'yhi': ext['ymax'],
               'zlo': ext['zmin'], 'zhi': ext['zmax']}
    lx = new_box['xhi']-new_box['xlo']
    ly = new_box['yhi']-new_box['ylo']
    lz = new_box['zhi']-new_box['zlo']
    print(f"  New box: Lx={lx:.2f}  Ly={ly:.2f}  Lz={lz:.2f}  "
          f"V={lx*ly*lz:.1f} σ³")
    write_lammps_data(output_file, atoms, bonds, new_box, masses)
    print(f"Wrote {len(atoms)} atoms → {output_file}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Isolate gel from slab_with_support data file "
                    "(removes bath solvent, support, piston).")
    parser.add_argument('--input',      required=True,
                        help="Input LAMMPS data file (final_config_*.data)")
    parser.add_argument('--output',     required=True,
                        help="Output LAMMPS data file (isolated_*.data)")
    parser.add_argument('--clearance',  type=float, default=BOX_CLEARANCE,
                        help=f"Box clearance around gel (default {BOX_CLEARANCE})")
    parser.add_argument('--percentile', type=float, default=0.1,
                        help="Percentile for outlier exclusion (default 0.1)")
    args = parser.parse_args()
    isolate_gel(args.input, args.output, args.clearance, args.percentile)
