#!/usr/bin/env python3
"""
isolate_gel.py
==============
CLI wrapper for isolating an equilibrated gel slab from a slab_with_support
LAMMPS data file.

Removes support (type 4) and piston (type 5) atoms, defines a control volume
strictly inside the polymer network using a percentile inset, and writes a new
data file containing only the atoms within that interior region.

The control volume is the inner (cv_percentile, 100-cv_percentile) percentile
of the polymer atom distribution along each axis.  No outward clearance is
added, so the box is guaranteed to lie inside the gel regardless of pressure.
This removes the variability in bath-solvent retention that arose when the
old clearance-based planes coincided with the gel face.

Usage
-----
  python3 isolate_gel.py --input <final_config.data> --output <isolated.data>
                         [--cv-percentile 5.0]

Output atom types
-----------------
  1  Crosslink   (unchanged)
  2  Chain bead  (unchanged)
  3  Solvent     (gel-interior only; bath solvent and surface atoms trimmed)
"""

import argparse

try:
    import numpy as np
    _NUMPY_AVAILABLE = True
except ImportError:
    _NUMPY_AVAILABLE = False
    print("WARNING: numpy not available — using pure-Python fallback (no rotation)")

try:
    from scipy.spatial import ConvexHull
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False

from collections import deque


# ---------------------------------------------------------------------------
# Pure-Python percentile (matches numpy's linear interpolation default)
# ---------------------------------------------------------------------------
def _percentile(data, pct):
    s = sorted(data)
    n = len(s)
    if n == 0:
        raise ValueError("Empty sequence for percentile")
    idx = (pct / 100.0) * (n - 1)
    lo = int(idx)
    hi = lo + 1
    if hi >= n:
        return s[-1]
    return s[lo] * (1 - (idx - lo)) + s[hi] * (idx - lo)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
POLYMER_TYPES = {1, 2}
SOLVENT_TYPES = {3}
SUPPORT_TYPES = {4}
PISTON_TYPES  = {5}
ROTATE_TYPES  = POLYMER_TYPES | SOLVENT_TYPES
CV_PERCENTILE = 3.0   # inner 94 % of polymer distribution per axis


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def unwrap_atoms_via_bonds(atoms, bonds, box):
    """
    Unwrap polymer atoms using BFS over the bond graph so that the
    entire gel network lands in one continuous piece (no periodic splits).
    Then fold solvent atoms to the nearest image of the polymer COM.

    Without this step, a gel that straddles a periodic boundary (common
    at lower pressures / larger boxes) produces a convex hull that spans
    the whole simulation box, corrupting the MABR angle and the gel-extent
    bounding box — leaving the gel rotated and retaining bath solvent.
    """
    if not _NUMPY_AVAILABLE:
        print("  Skipping PBC unwrap (numpy unavailable)")
        return atoms

    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']

    atom_by_id   = {a['id']: a for a in atoms}
    polymer_ids  = {a['id'] for a in atoms if a['type'] in POLYMER_TYPES}

    # Build adjacency list for polymer only
    adj = {pid: [] for pid in polymer_ids}
    for b in bonds:
        a1, a2 = b['atom1'], b['atom2']
        if a1 in polymer_ids and a2 in polymer_ids:
            adj[a1].append(a2)
            adj[a2].append(a1)

    # BFS: place each bonded neighbour in the same image as its parent
    visited = set()
    seed    = next(iter(polymer_ids))
    queue   = deque([seed])
    visited.add(seed)

    while queue:
        curr_id = queue.popleft()
        curr    = atom_by_id[curr_id]
        for nb_id in adj[curr_id]:
            if nb_id in visited:
                continue
            nb = atom_by_id[nb_id]
            nb['x'] -= lx * round((nb['x'] - curr['x']) / lx)
            nb['y'] -= ly * round((nb['y'] - curr['y']) / ly)
            nb['z'] -= lz * round((nb['z'] - curr['z']) / lz)
            visited.add(nb_id)
            queue.append(nb_id)

    if len(visited) < len(polymer_ids):
        print(f"  WARNING: bond graph has {len(polymer_ids) - len(visited)} "
              f"disconnected polymer atoms — they keep their wrapped positions")

    # Fold solvent to nearest image of polymer COM
    poly_atoms = [atom_by_id[pid] for pid in polymer_ids]
    cx = float(np.mean([a['x'] for a in poly_atoms]))
    cy = float(np.mean([a['y'] for a in poly_atoms]))
    cz = float(np.mean([a['z'] for a in poly_atoms]))

    for a in atoms:
        if a['type'] in SOLVENT_TYPES:
            a['x'] -= lx * round((a['x'] - cx) / lx)
            a['y'] -= ly * round((a['y'] - cy) / ly)
            a['z'] -= lz * round((a['z'] - cz) / lz)

    print(f"  Unwrapped {len(visited)} polymer atoms; "
          f"folded solvent near COM ({cx:.1f}, {cy:.1f}, {cz:.1f})")
    return atoms


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
    if not (_NUMPY_AVAILABLE and _SCIPY_AVAILABLE):
        print("  Skipping rotation (numpy/scipy unavailable) — using axis-aligned bounding box")
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


def find_gel_extent(atoms, cv_percentile=CV_PERCENTILE):
    """
    Return a control volume defined by the inner (cv_percentile, 100-cv_percentile)
    percentile of polymer atom positions along each axis.  The bounds are strictly
    inside the gel — no outward clearance is added — so bath solvent can never
    leak in regardless of how the gel face sits relative to the box edge.
    """
    poly = [a for a in atoms if a['type'] in POLYMER_TYPES]
    xs = [a['x'] for a in poly]
    ys = [a['y'] for a in poly]
    zs = [a['z'] for a in poly]
    if _NUMPY_AVAILABLE:
        pct = lambda arr, p: float(np.percentile(arr, p))
    else:
        pct = _percentile
    xmin, xmax = pct(xs, cv_percentile), pct(xs, 100 - cv_percentile)
    ymin, ymax = pct(ys, cv_percentile), pct(ys, 100 - cv_percentile)
    zmin, zmax = pct(zs, cv_percentile), pct(zs, 100 - cv_percentile)
    lx = xmax - xmin;  ly = ymax - ymin;  lz = zmax - zmin
    print(f"  Control volume  x: {xmin:.2f}–{xmax:.2f}  "
          f"y: {ymin:.2f}–{ymax:.2f}  z: {zmin:.2f}–{zmax:.2f}")
    print(f"  CV dimensions   Lx={lx:.2f}  Ly={ly:.2f}  Lz={lz:.2f}  "
          f"V={lx*ly*lz:.1f} σ³  (cv_percentile={cv_percentile})")
    return dict(xmin=xmin, xmax=xmax,
                ymin=ymin, ymax=ymax,
                zmin=zmin, zmax=zmax)


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
          f"{n_outside} atoms outside control volume")
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
def isolate_gel(input_file, output_file, cv_percentile=CV_PERCENTILE):
    print(f"Input:        {input_file}")
    print(f"Output:       {output_file}")
    print(f"cv_percentile: {cv_percentile}  (inner {100 - 2*cv_percentile:.0f}% of polymer distribution)")
    atoms, bonds, box, masses = parse_lammps_data(input_file)
    print(f"Read {len(atoms)} atoms, {len(bonds)} bonds")
    atoms = unwrap_atoms_via_bonds(atoms, bonds, box)
    atoms = rotate_mobile_atoms(atoms)
    ext   = find_gel_extent(atoms, cv_percentile)
    atoms = remove_non_gel_atoms(atoms, ext)
    new_box = {'xlo': ext['xmin'], 'xhi': ext['xmax'],
               'ylo': ext['ymin'], 'yhi': ext['ymax'],
               'zlo': ext['zmin'], 'zhi': ext['zmax']}
    write_lammps_data(output_file, atoms, bonds, new_box, masses)
    print(f"Wrote {len(atoms)} atoms → {output_file}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Isolate gel interior from slab_with_support data file "
                    "(removes bath solvent, support, piston, and gel surface layer).")
    parser.add_argument('--input',         required=True,
                        help="Input LAMMPS data file (final_config_*.data)")
    parser.add_argument('--output',        required=True,
                        help="Output LAMMPS data file (isolated_*.data)")
    parser.add_argument('--cv-percentile', type=float, default=CV_PERCENTILE,
                        help=f"Percentile inset for control volume on each axis "
                             f"(default {CV_PERCENTILE}). "
                             f"E.g. 3.0 keeps the inner 94%% of the polymer distribution.")
    args = parser.parse_args()
    isolate_gel(args.input, args.output, args.cv_percentile)
