#!/usr/bin/env python3
"""
isolate_gel.py
==============
CLI wrapper for isolating an equilibrated gel slab from a slab_with_support
LAMMPS data file.

Removes support (type 4) and piston (type 5) atoms, trims any solvent/polymer
outside a tight bounding box around the polymer network, and writes a new data
file containing only gel-internal atoms (types 1=crosslink, 2=chain, 3=solvent).

This isolated file is now used as a *fully periodic mixed-gel* input: it is
re-equilibrated under NPT (via polymer_pure.lmp) so that V_mix is a true
thermodynamic box volume, measured on the same footing as the pure-species
references. The box written here is only an initial condition.

Why the redesign (volume-of-mixing pipeline, 2026-06):
  * V_mix used to be THIS geometric bounding box, while V_pol/V_sol were NPT
    box volumes — subtracting two different kinds of volume is not a valid
    ΔV_mix. V_mix is now an NPT volume too.
  * Clearance was shown (sensitivity sweep) to be a near-constant additive
    offset (~0.018 in ΔV_mix/(V_sol+V_pol) at every P*), not the source of the
    pressure trend. So clearance here no longer feeds the result; it is purely
    a PBC-safe initial gap and is raised 0.2 -> 0.5σ so periodic images of edge
    atoms start ~1σ apart (no self-overlap at t=0).

Usage
-----
  python3 isolate_gel.py --input <final_config.data> --output <isolated.data>
                         [--clearance 0.5] [--percentile 0.1]

Output atom types
-----------------
  1  Crosslink   (unchanged)
  2  Chain bead  (unchanged)
  3  Solvent     (gel-internal only; bath solvent trimmed)
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
# PBC-safe initial gap (per face). NOT a result-affecting parameter anymore:
# V_mix is set by the downstream NPT run, so this only places periodic images
# of edge atoms ~1σ apart at t=0 to avoid self-overlap. (Was 0.2 when this box
# WAS the reported V_mix.)
BOX_CLEARANCE = 0.5
# Maximum extent of the FENE bond used by the engines (bond_coeff 1 30.0 1.5 1.0 1.0).
# Any bond longer than this in the written file is a broken network, not a stretched one.
FENE_R0_MAX = 1.5


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def detect_percolating_dims(atoms, bonds, box):
    """
    Return the set of dimensions ('x','y','z') in which the polymer network
    percolates, i.e. in which at least one bond crosses the periodic boundary.

    This distinguishes the two geometries this script can be handed:

      * A FINITE gel blob floating in a bath (no bond crosses any boundary).
        It can legitimately be unwrapped, rotated onto its principal axes and
        re-boxed to a tight bounding box — the original operations here.

      * A gel that is PERIODIC (percolating) in one or more directions, e.g.
        the slab_support_periodic snapshots, which percolate in x and y and
        are finite only in z.

    A percolating direction must be left completely alone. Unwrapping a
    network that closes on itself through the boundary is ill-defined: BFS
    over the bond graph assigns images consistently until it reaches a loop
    that wraps the box, and that loop-closure bond is then left stretched by
    a full box length. Rotating or re-boxing such a direction breaks every
    boundary-crossing bond in the same way, because the periodic images no
    longer tile the new box.
    """
    lengths = {'x': box['xhi'] - box['xlo'],
               'y': box['yhi'] - box['ylo'],
               'z': box['zhi'] - box['zlo']}
    pos = {a['id']: a for a in atoms}
    crossings = {'x': 0, 'y': 0, 'z': 0}
    for b in bonds:
        a1 = pos.get(b['atom1']);  a2 = pos.get(b['atom2'])
        if a1 is None or a2 is None:
            continue
        for d in ('x', 'y', 'z'):
            if abs(a2[d] - a1[d]) > lengths[d] / 2.0:
                crossings[d] += 1
    periodic = {d for d in ('x', 'y', 'z') if crossings[d] > 0}
    desc = ", ".join(f"{d}:{crossings[d]}" for d in ('x', 'y', 'z'))
    print(f"  Boundary-crossing bonds ({desc})")
    if periodic:
        print(f"  Network PERCOLATES in {sorted(periodic)} — those directions "
              f"keep the original box, coordinates and periodicity")
    else:
        print("  Network is finite in all directions — free to unwrap/rotate/re-box")
    return periodic


def unwrap_atoms_via_bonds(atoms, bonds, box, periodic_dims=frozenset()):
    """
    Unwrap polymer atoms using BFS over the bond graph so that the
    entire gel network lands in one continuous piece (no periodic splits).
    Then fold solvent atoms to the nearest image of the polymer COM.

    Without this step, a gel that straddles a periodic boundary (common
    at lower pressures / larger boxes) produces a convex hull that spans
    the whole simulation box, corrupting the MABR angle and the gel-extent
    bounding box — leaving the gel rotated and retaining bath solvent.

    Dimensions listed in `periodic_dims` are skipped: the network closes on
    itself through those boundaries, so there is no consistent unwrapping and
    the wrapped coordinates are already the correct ones.
    """
    if not _NUMPY_AVAILABLE:
        print("  Skipping PBC unwrap (numpy unavailable)")
        return atoms

    free = [d for d in ('x', 'y', 'z') if d not in periodic_dims]
    if not free:
        print("  Skipping PBC unwrap — network percolates in all three directions")
        return atoms
    if periodic_dims:
        print(f"  Unwrapping only in {free} (periodic: {sorted(periodic_dims)})")

    L = {'x': box['xhi'] - box['xlo'],
         'y': box['yhi'] - box['ylo'],
         'z': box['zhi'] - box['zlo']}

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
            for d in free:
                nb[d] -= L[d] * round((nb[d] - curr[d]) / L[d])
            visited.add(nb_id)
            queue.append(nb_id)

    if len(visited) < len(polymer_ids):
        print(f"  WARNING: bond graph has {len(polymer_ids) - len(visited)} "
              f"disconnected polymer atoms — they keep their wrapped positions")

    # Fold solvent to nearest image of polymer COM (free dimensions only)
    poly_atoms = [atom_by_id[pid] for pid in polymer_ids]
    com = {d: float(np.mean([a[d] for a in poly_atoms])) for d in ('x', 'y', 'z')}

    for a in atoms:
        if a['type'] in SOLVENT_TYPES:
            for d in free:
                a[d] -= L[d] * round((a[d] - com[d]) / L[d])

    print(f"  Unwrapped {len(visited)} polymer atoms; folded solvent near COM "
          f"({com['x']:.1f}, {com['y']:.1f}, {com['z']:.1f})")
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


def rotate_mobile_atoms(atoms, periodic_dims=frozenset()):
    if not (_NUMPY_AVAILABLE and _SCIPY_AVAILABLE):
        print("  Skipping rotation (numpy/scipy unavailable) — using axis-aligned bounding box")
        return atoms
    # The MABR rotation is about the z axis, so it moves x and y. A rotation is
    # only a symmetry of the system when both are non-periodic: rotating a
    # percolating direction leaves the periodic images no longer tiling the
    # box, which stretches every boundary-crossing bond past the FENE limit.
    if periodic_dims & {'x', 'y'}:
        print(f"  Skipping rotation — network percolates in "
              f"{sorted(periodic_dims & {'x', 'y'})} (rotation is not a symmetry there)")
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


def find_gel_extent(atoms, clearance, percentile=0.1,
                    box=None, periodic_dims=frozenset()):
    """
    Bounding box around the polymer, per dimension.

    A percolating dimension keeps the ORIGINAL box bounds exactly. Shrinking
    it would change that direction's box length, so bonds that cross the
    boundary — which are correct only under the original length — would all
    break. Only the finite directions are re-boxed.

    A finite dimension spans EVERY polymer atom, not a percentile of them.
    The box has to contain the whole network for two reasons: any polymer atom
    left outside would be wrapped to the opposite face on read_data, and
    deleting it instead severs its bonds and leaves dangling chain ends. So the
    only thing trimmed off a finite face is bath solvent, which is exactly what
    "isolating the gel" means here.

    `percentile` no longer cuts the network; it is now purely diagnostic,
    reporting how far the outermost beads reach beyond the bulk of the surface.
    The extra box length that buys is harmless — this box is only the starting
    configuration, and V_mix comes from the downstream NPT run.
    """
    poly = [a for a in atoms if a['type'] in POLYMER_TYPES]
    if _NUMPY_AVAILABLE:
        pct = lambda arr, p: float(np.percentile(arr, p))
    else:
        pct = _percentile

    ext = {}
    for d in ('x', 'y', 'z'):
        if d in periodic_dims:
            if box is None:
                raise ValueError("box is required to preserve periodic dimensions")
            lo, hi = box[f'{d}lo'], box[f'{d}hi']
            ext[f'{d}min'], ext[f'{d}max'] = lo, hi
            print(f"  {d}: {lo:.2f}–{hi:.2f} (periodic — original box kept)")
        else:
            vals = [a[d] for a in poly]
            lo, hi = min(vals), max(vals)
            ext[f'{d}min'], ext[f'{d}max'] = lo - clearance, hi + clearance
            plo, phi = pct(vals, percentile), pct(vals, 100 - percentile)
            n_out = sum(1 for v in vals if v < plo or v > phi)
            print(f"  {d}: {lo:.2f}–{hi:.2f} (finite — full polymer span, "
                  f"+{clearance} clearance)")
            print(f"     surface roughness: {n_out} polymer atoms beyond "
                  f"p{percentile}/p{100-percentile} ({plo:.2f}–{phi:.2f}) — kept, not trimmed")
    return ext


def remove_non_gel_atoms(atoms, ext, periodic_dims=frozenset()):
    # Whitelist gel types (1,2,3): drop support/piston AND any other stray type
    # so the isolated mixed file is always exactly polymer+solvent. This keeps
    # the mixed-gel NPT atom counts identical to the split pure-species files.
    #
    # Percolating dimensions are not tested: their bounds are the original box,
    # so every atom is inside by construction, and a floating-point comparison
    # at the boundary could otherwise delete an atom and orphan its bonds.
    # The polymer network is NEVER trimmed. find_gel_extent sizes the box to
    # contain every polymer atom, so the only thing a face cuts away is bath
    # solvent. Deleting a bead instead would orphan its bonds — write_lammps_data
    # silently drops any bond whose partner is gone, so the file would still load
    # and run, just with severed chain ends nobody asked for.
    gel_types = POLYMER_TYPES | SOLVENT_TYPES
    test_dims = [d for d in ('x', 'y', 'z') if d not in periodic_dims]
    kept = []
    n_nongel = n_outside = 0
    for a in atoms:
        if a['type'] not in gel_types:
            n_nongel += 1;  continue
        if (a['type'] not in POLYMER_TYPES and
                any(a[d] < ext[f'{d}min'] or a[d] > ext[f'{d}max'] for d in test_dims)):
            n_outside += 1;  continue
        kept.append(a)

    n_poly_in  = sum(1 for a in atoms if a['type'] in POLYMER_TYPES)
    n_poly_out = sum(1 for a in kept  if a['type'] in POLYMER_TYPES)
    if n_poly_in != n_poly_out:
        raise ValueError(f"{n_poly_in - n_poly_out} polymer atoms were dropped; "
                         f"the network must be preserved intact")
    print(f"  Removed {n_nongel} non-gel atoms (support/piston/other), "
          f"{n_outside} bath solvent atoms outside gel bounds")
    print(f"  Kept all {n_poly_out} polymer atoms (network intact)")
    return kept


def validate_bonds(atoms, bonds, box, max_bond=FENE_R0_MAX):
    """
    Post-condition check: every surviving bond must be shorter than the FENE
    maximum extent under the minimum-image convention of the NEW box.

    A violation means the geometry pipeline broke the network (see
    detect_percolating_dims). Downstream this surfaces only as an opaque
    'Bond atom missing in image check' at LAMMPS runtime, hours later and
    after the job has been recorded as submitted, so fail here instead.
    """
    L = {'x': box['xhi'] - box['xlo'],
         'y': box['yhi'] - box['ylo'],
         'z': box['zhi'] - box['zlo']}
    pos = {a['id']: a for a in atoms}
    worst = 0.0
    n_bad = 0
    n_orphan = 0
    for b in bonds:
        a1 = pos.get(b['atom1']);  a2 = pos.get(b['atom2'])
        if a1 is None or a2 is None:
            n_orphan += 1
            continue
        r2 = 0.0
        for d in ('x', 'y', 'z'):
            delta = a2[d] - a1[d]
            delta -= L[d] * round(delta / L[d])
            r2 += delta * delta
        r = r2 ** 0.5
        worst = max(worst, r)
        if r > max_bond:
            n_bad += 1
    print(f"  Bond check: longest min-image bond = {worst:.3f} σ "
          f"(FENE limit {max_bond}); {len(bonds) - n_orphan} bonds kept, "
          f"{n_orphan} severed")
    if n_orphan:
        raise ValueError(
            f"{n_orphan} bond(s) lost an endpoint to trimming. The network must "
            f"come through whole — a severed bond is silently dropped by "
            f"write_lammps_data and leaves a dangling chain end in the run.")
    if n_bad:
        raise ValueError(
            f"{n_bad} bond(s) exceed the FENE maximum extent of {max_bond} σ "
            f"(longest {worst:.3f} σ). The isolated configuration is not a "
            f"valid FENE network and would fail in LAMMPS. This usually means "
            f"a percolating direction was unwrapped, rotated or re-boxed.")


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
    periodic = detect_percolating_dims(atoms, bonds, box)
    atoms = unwrap_atoms_via_bonds(atoms, bonds, box, periodic)
    atoms = rotate_mobile_atoms(atoms, periodic)
    ext   = find_gel_extent(atoms, clearance, percentile, box, periodic)
    atoms = remove_non_gel_atoms(atoms, ext, periodic)
    new_box = {'xlo': ext['xmin'], 'xhi': ext['xmax'],
               'ylo': ext['ymin'], 'yhi': ext['ymax'],
               'zlo': ext['zmin'], 'zhi': ext['zmax']}
    lx = new_box['xhi']-new_box['xlo']
    ly = new_box['yhi']-new_box['ylo']
    lz = new_box['zhi']-new_box['zlo']
    print(f"  New box: Lx={lx:.2f}  Ly={ly:.2f}  Lz={lz:.2f}  "
          f"V={lx*ly*lz:.1f} σ³")
    # Validate BEFORE writing: a broken network must never reach the disk, or
    # the split/adjust steps downstream will happily propagate it into every
    # per-loading input file.
    validate_bonds(atoms, bonds, new_box)
    write_lammps_data(output_file, atoms, bonds, new_box, masses)
    print(f"Wrote {len(atoms)} atoms → {output_file}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Isolate gel interior from slab_with_support data file "
                    "(removes bath solvent, support, piston, and gel surface layer).")
    parser.add_argument('--input',      required=True,
                        help="Input LAMMPS data file (final_config_*.data)")
    parser.add_argument('--output',     required=True,
                        help="Output LAMMPS data file (isolated_*.data)")
    parser.add_argument('--clearance',  type=float, default=BOX_CLEARANCE,
                        help=f"PBC-safe initial gap per face (default {BOX_CLEARANCE}); "
                             f"only an initial condition — V_mix comes from the NPT run")
    parser.add_argument('--percentile', type=float, default=0.1,
                        help="Percentile for outlier exclusion (default 0.1)")
    args = parser.parse_args()
    isolate_gel(args.input, args.output, args.clearance, args.percentile)
