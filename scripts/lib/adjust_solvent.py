#!/usr/bin/env python3
"""
adjust_solvent.py
=================
CLI tool for setting the solvent count N_f of an isolated mixed-gel LAMMPS
data file (the output of isolate_gel.py) to an exact target value.

Used by the calibration sweep (2026-08 plan) to realize a fixed grid of
composition targets from ONE network topology: the φ_p grid is set by DIRECTLY
adding/deleting solvent, decoupled from the swelling curve, so every pressure
sees identical composition coverage.

  * target N_f < current: randomly chosen solvent beads (type 3) are deleted.
  * target N_f > current: new solvent beads are inserted at random positions
    at least --min-dist from every existing atom (PBC-aware). If a placement
    stalls, the clearance is relaxed in 0.05σ steps down to a floor of 0.7σ —
    residual soft overlaps at that scale are absorbed by the downstream
    engines' gentle-start stages (polymer_pure.lmp: harmonic pre-relax +
    nve/limit; solvent_pure.lmp: minimize + nve/limit).

The polymer (types 1, 2) and all bonds are untouched, so N_p is constant
across the whole composition grid — the Euler-closure requirement.

Output file: <input stem>_nf<TARGET>.data (or --output). Run split_gel.py on
it afterwards for the per-loading _solvent_only ΔV_mix companion.

Usage
-----
  python3 adjust_solvent.py --input <isolated.data> --target-nf <N>
                            [--output <file.data>] [--seed 12345]
                            [--min-dist 0.9]

Atom types (unchanged from isolate_gel.py output)
-------------------------------------------------
  1  Crosslink   (untouched)
  2  Chain bead  (untouched)
  3  Solvent     (deleted or inserted to reach the target)
"""

import argparse
import os
import random

try:
    import numpy as np
    _NUMPY_AVAILABLE = True
except ImportError:
    _NUMPY_AVAILABLE = False

try:
    from scipy.spatial import cKDTree
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
POLYMER_TYPES = {1, 2}
SOLVENT_TYPE  = 3
MIN_DIST      = 0.9    # default insertion clearance (σ) from any existing atom
MIN_DIST_FLOOR = 0.7   # hard floor for the progressive relaxation — deeper
                       # overlaps risk force blow-ups even under nve/limit
BATCH         = 4096   # candidate positions drawn per rejection-sampling round
MAX_ROUNDS    = 200    # rounds per clearance level before relaxing it


# ---------------------------------------------------------------------------
# LAMMPS I/O (same format/conventions as isolate_gel.py)
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


def write_lammps_data(path, title, atoms, bonds, box, masses):
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
        f.write(f"{title}\n\n")
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
# Deletion
# ---------------------------------------------------------------------------
def delete_solvent(atoms, n_delete, rng):
    solvent_ids = [a['id'] for a in atoms if a['type'] == SOLVENT_TYPE]
    if n_delete > len(solvent_ids):
        raise ValueError(f"Cannot delete {n_delete} solvent beads — "
                         f"only {len(solvent_ids)} present")
    doomed = set(rng.sample(solvent_ids, n_delete))
    kept = [a for a in atoms if a['id'] not in doomed]
    print(f"  Deleted {n_delete} randomly chosen solvent beads")
    return kept


# ---------------------------------------------------------------------------
# Insertion
# ---------------------------------------------------------------------------
def insert_solvent(atoms, n_insert, box, rng, min_dist=MIN_DIST):
    """Insert n_insert solvent beads at random PBC-aware positions ≥ min_dist
    from every existing atom (and every already-accepted insertion). The
    clearance relaxes progressively (0.05σ steps, floor MIN_DIST_FLOOR) if the
    box is too crowded — downstream gentle-start stages absorb the resulting
    soft overlaps."""
    if not (_NUMPY_AVAILABLE and _SCIPY_AVAILABLE):
        raise RuntimeError("Insertion requires numpy + scipy (cKDTree). "
                           "Deletion-only use works without them.")

    lo = np.array([box['xlo'], box['ylo'], box['zlo']])
    L  = np.array([box['xhi'] - box['xlo'],
                   box['yhi'] - box['ylo'],
                   box['zhi'] - box['zlo']])
    # Wrap existing coordinates into [0, L) — the isolated box is periodic, and
    # cKDTree's boxsize metric needs canonical images.
    existing = np.array([[a['x'], a['y'], a['z']] for a in atoms])
    existing = np.mod(existing - lo, L)
    tree = cKDTree(existing, boxsize=L)

    np_rng   = np.random.default_rng(rng.getrandbits(64))
    accepted = np.empty((0, 3))
    rounds_at_level = 0

    while len(accepted) < n_insert:
        cand = np_rng.random((BATCH, 3)) * L
        # clearance vs existing atoms
        d_exist = tree.query(cand, k=1)[0]
        cand = cand[d_exist >= min_dist]
        # clearance vs already-accepted insertions from earlier batches
        if len(accepted) and len(cand):
            delta = np.abs(cand[:, None, :] - accepted[None, :, :])
            delta = np.minimum(delta, L - delta)
            d_acc = np.sqrt((delta**2).sum(axis=2)).min(axis=1)
            cand = cand[d_acc >= min_dist]
        # clearance within the batch itself: greedily drop the higher-indexed
        # member of every too-close pair (survivors are mutually separated)
        if len(cand) > 1:
            pairs = cKDTree(cand, boxsize=L).query_pairs(min_dist)
            drop = np.zeros(len(cand), dtype=bool)
            for i, j in sorted(pairs):
                if not drop[i] and not drop[j]:
                    drop[j] = True
            cand = cand[~drop]
        n_take = min(len(cand), n_insert - len(accepted))
        if n_take:
            accepted = np.vstack([accepted, cand[:n_take]])
            rounds_at_level = 0
        else:
            rounds_at_level += 1
            if rounds_at_level >= MAX_ROUNDS:
                if min_dist - 0.05 < MIN_DIST_FLOOR - 1e-9:
                    raise RuntimeError(
                        f"Could not place bead {len(accepted)+1}/{n_insert} even at "
                        f"the clearance floor {MIN_DIST_FLOOR}σ — box too crowded. "
                        f"Lower the target N_f or pre-expand the box.")
                min_dist = round(min_dist - 0.05, 2)
                rounds_at_level = 0
                print(f"  WARNING: placement stalling — clearance relaxed to {min_dist}σ "
                      f"({len(accepted)}/{n_insert} placed)")

    next_id  = max(a['id']  for a in atoms) + 1
    next_mol = max(a['mol'] for a in atoms) + 1
    out = accepted + lo   # back to the file's coordinate origin
    new_atoms = [{'id': next_id + i, 'mol': next_mol + i, 'type': SOLVENT_TYPE,
                  'x': float(p[0]), 'y': float(p[1]), 'z': float(p[2])}
                 for i, p in enumerate(out)]
    print(f"  Inserted {n_insert} solvent beads (final clearance {min_dist}σ)")
    return atoms + new_atoms


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def adjust_solvent(input_file, target_nf, output_file=None, seed=12345,
                   min_dist=MIN_DIST):
    if output_file is None:
        stem, ext = os.path.splitext(input_file)
        output_file = f"{stem}_nf{target_nf}{ext or '.data'}"

    print(f"Input:      {input_file}")
    print(f"Output:     {output_file}")
    print(f"Target N_f: {target_nf}  |  Seed: {seed}  |  Min dist: {min_dist}σ")

    atoms, bonds, box, masses = parse_lammps_data(input_file)
    n_f = sum(1 for a in atoms if a['type'] == SOLVENT_TYPE)
    n_p = sum(1 for a in atoms if a['type'] in POLYMER_TYPES)
    n_other = len(atoms) - n_f - n_p
    if n_other:
        raise ValueError(f"Input has {n_other} atoms of unexpected type — "
                         f"expected an isolate_gel.py output (types 1,2,3 only)")
    print(f"Read {len(atoms)} atoms ({n_p} polymer, {n_f} solvent), "
          f"{len(bonds)} bonds")

    rng = random.Random(seed)
    delta = target_nf - n_f
    if delta < 0:
        atoms = delete_solvent(atoms, -delta, rng)
    elif delta > 0:
        atoms = insert_solvent(atoms, delta, box, rng, min_dist)
    else:
        print("  Target equals current N_f — atoms unchanged")

    n_f_out = sum(1 for a in atoms if a['type'] == SOLVENT_TYPE)
    assert n_f_out == target_nf, f"N_f bookkeeping error: {n_f_out} != {target_nf}"

    title = (f"LAMMPS data file — solvent-adjusted mixed gel "
             f"(N_f {n_f} -> {target_nf}, N_p {n_p}, seed {seed})")
    write_lammps_data(output_file, title, atoms, bonds, box, masses)
    print(f"Wrote {len(atoms)} atoms ({n_p} polymer, {n_f_out} solvent) "
          f"→ {output_file}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Set the solvent count N_f of an isolated mixed-gel data "
                    "file to an exact target (random delete / non-overlapping "
                    "random insert). Polymer and bonds are untouched.")
    parser.add_argument('--input',     required=True,
                        help="Input LAMMPS data file (isolated_*.data from isolate_gel.py)")
    parser.add_argument('--target-nf', required=True, type=int,
                        help="Exact solvent bead count for the output file")
    parser.add_argument('--output',    default=None,
                        help="Output data file (default: <input stem>_nf<TARGET>.data)")
    parser.add_argument('--seed',      type=int, default=12345,
                        help="RNG seed for delete choice / insert positions (default 12345)")
    parser.add_argument('--min-dist',  type=float, default=MIN_DIST,
                        help=f"Insertion clearance from any existing atom in σ "
                             f"(default {MIN_DIST}; auto-relaxes to {MIN_DIST_FLOOR} "
                             f"floor if the box is crowded)")
    args = parser.parse_args()
    adjust_solvent(args.input, args.target_nf, args.output, args.seed,
                   args.min_dist)
