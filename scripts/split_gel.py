"""
split_gel.py
============
Standalone CLI wrapper around the split_gel function from split_gel_slab.ipynb.

Reads an equilibrated LAMMPS data file (atom_style molecular, 5 atom types)
and writes two new data files:

  <stem>_polymer_only.data  — types 1 (crosslink) and 2 (chain bead) only
                               positions unchanged; bonds preserved
  <stem>_solvent_only.data  — types 3 (solvent), 4 (support), 5 (piston) only
                               positions unchanged; no bonds
                               type remapping: 3→1, 4→2, 5→3

Usage
-----
  python split_gel.py <input.data> [--output-dir DIR] [--output-stem STEM]

  If --output-dir is given, output files are written there (keeping the
  input basename as stem).  If --output-stem is given directly, that overrides
  both the directory and input filename.  Default: same directory as input.

Atom types in the original file
--------------------------------
  1  Crosslink  (polymer)
  2  Chain bead (polymer)
  3  Solvent
  4  Bottom support (frozen)
  5  Top piston  (mobile)
"""

import sys
import os
import re
import argparse

try:
    import numpy as np
    _NUMPY = True
except ImportError:
    _NUMPY = False


# ---------------------------------------------------------------------------
# Percentile helper (numpy or pure-Python fallback)
# ---------------------------------------------------------------------------
def _pct(data, p):
    if _NUMPY:
        return float(np.percentile(data, p))
    s = sorted(data)
    n = len(s)
    if n == 0:
        raise ValueError("Empty sequence")
    idx = (p / 100.0) * (n - 1)
    lo = int(idx)
    hi = min(lo + 1, n - 1)
    return s[lo] * (1 - (idx - lo)) + s[hi] * (idx - lo)


def trim_to_cv(atoms, cv_percentile, label=""):
    """
    Trim *atoms* to the inner (cv_percentile, 100-cv_percentile) percentile
    of their own x/y/z distribution.  Returns (trimmed_atoms, cv_box).
    cv_box is a dict with xlo/xhi/ylo/yhi/zlo/zhi keys.
    """
    if not atoms or cv_percentile <= 0:
        xs = [a['x'] for a in atoms]
        ys = [a['y'] for a in atoms]
        zs = [a['z'] for a in atoms]
        box = {'xlo': min(xs), 'xhi': max(xs),
               'ylo': min(ys), 'yhi': max(ys),
               'zlo': min(zs), 'zhi': max(zs)}
        return atoms, box

    xs = [a['x'] for a in atoms]
    ys = [a['y'] for a in atoms]
    zs = [a['z'] for a in atoms]
    xlo = _pct(xs, cv_percentile);  xhi = _pct(xs, 100 - cv_percentile)
    ylo = _pct(ys, cv_percentile);  yhi = _pct(ys, 100 - cv_percentile)
    zlo = _pct(zs, cv_percentile);  zhi = _pct(zs, 100 - cv_percentile)
    kept = [a for a in atoms
            if xlo <= a['x'] <= xhi and ylo <= a['y'] <= yhi and zlo <= a['z'] <= zhi]
    cv_box = {'xlo': xlo, 'xhi': xhi, 'ylo': ylo, 'yhi': yhi, 'zlo': zlo, 'zhi': zhi}
    lx = xhi - xlo;  ly = yhi - ylo;  lz = zhi - zlo
    tag = f" ({label})" if label else ""
    print(f"  CV trim {cv_percentile}%{tag}: {len(kept)}/{len(atoms)} atoms  "
          f"Lx={lx:.2f} Ly={ly:.2f} Lz={lz:.2f} V={lx*ly*lz:.1f} σ³")
    return kept, cv_box


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_data_file(path):
    """
    Parse a LAMMPS data file (atom_style molecular).
    Returns a dict with keys:
        header    : dict of scalar counts/bounds
        masses    : list of (type_id, mass, comment)
        atoms     : list of dicts {id, mol, type, x, y, z}
        bonds     : list of dicts {id, type, atom1, atom2}
    """
    with open(path) as f:
        lines = f.readlines()

    header = {
        'n_atoms': 0, 'n_bonds': 0,
        'n_angles': 0, 'n_dihedrals': 0, 'n_impropers': 0,
        'n_atom_types': 0, 'n_bond_types': 0,
        'xlo': 0.0, 'xhi': 0.0,
        'ylo': 0.0, 'yhi': 0.0,
        'zlo': 0.0, 'zhi': 0.0,
        'title': ''
    }
    masses = []
    atoms  = []
    bonds  = []

    section = None
    for i, raw in enumerate(lines):
        line = raw.strip()

        if i == 0:
            header['title'] = line
            continue

        if not line or line.startswith('#'):
            continue

        line = line.split('#')[0].strip()
        if not line:
            continue

        if re.match(r'^[A-Za-z][A-Za-z /]+$', line):
            section = line
            continue

        if section is None:
            for kw, key in [
                ('atoms',        'n_atoms'),
                ('bonds',        'n_bonds'),
                ('angles',       'n_angles'),
                ('dihedrals',    'n_dihedrals'),
                ('impropers',    'n_impropers'),
                ('atom types',   'n_atom_types'),
                ('bond types',   'n_bond_types'),
            ]:
                if line.endswith(kw):
                    header[key] = int(line.split()[0])
                    break

            m = re.match(r'^([-\d.eE+]+)\s+([-\d.eE+]+)\s+xlo\s+xhi', line)
            if m:
                header['xlo'], header['xhi'] = float(m.group(1)), float(m.group(2))
            m = re.match(r'^([-\d.eE+]+)\s+([-\d.eE+]+)\s+ylo\s+yhi', line)
            if m:
                header['ylo'], header['yhi'] = float(m.group(1)), float(m.group(2))
            m = re.match(r'^([-\d.eE+]+)\s+([-\d.eE+]+)\s+zlo\s+zhi', line)
            if m:
                header['zlo'], header['zhi'] = float(m.group(1)), float(m.group(2))
            continue

        if section == 'Masses':
            parts = line.split('#', 1)
            nums  = parts[0].split()
            comment = parts[1].strip() if len(parts) > 1 else ''
            masses.append((int(nums[0]), float(nums[1]), comment))

        elif section == 'Atoms':
            parts = line.split()
            atoms.append({
                'id':  int(parts[0]),
                'mol': int(parts[1]),
                'type': int(parts[2]),
                'x': float(parts[3]),
                'y': float(parts[4]),
                'z': float(parts[5]),
            })

        elif section == 'Bonds':
            parts = line.split()
            bonds.append({
                'id':    int(parts[0]),
                'type':  int(parts[1]),
                'atom1': int(parts[2]),
                'atom2': int(parts[3]),
            })

    return {'header': header, 'masses': masses, 'atoms': atoms, 'bonds': bonds}


def write_data_file(path, title, box, atom_types_info, atoms, bonds=None):
    """Write a LAMMPS data file (atom_style molecular)."""
    n_bond_types = len(set(b['type'] for b in bonds)) if bonds else 0

    with open(path, 'w') as f:
        f.write(f"{title}\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(bonds) if bonds else 0} bonds\n")
        f.write("0 angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n\n")
        f.write(f"{len(atom_types_info)} atom types\n")
        if n_bond_types:
            f.write(f"{n_bond_types} bond types\n")
        f.write("\n")
        f.write(f"{box['xlo']:.6f} {box['xhi']:.6f} xlo xhi\n")
        f.write(f"{box['ylo']:.6f} {box['yhi']:.6f} ylo yhi\n")
        f.write(f"{box['zlo']:.6f} {box['zhi']:.6f} zlo zhi\n")
        f.write("\nMasses\n\n")
        for tid, mass, label in atom_types_info:
            comment = f"  # {label}" if label else ''
            f.write(f"{tid} {mass:.1f}{comment}\n")

        f.write("\nAtoms\n\n")
        for a in atoms:
            f.write(f"{a['id']} {a['mol']} {a['type']} "
                    f"{a['x']:.6f} {a['y']:.6f} {a['z']:.6f}\n")

        if bonds:
            f.write("\nBonds\n\n")
            for b in bonds:
                f.write(f"{b['id']} {b['type']} {b['atom1']} {b['atom2']}\n")


# ---------------------------------------------------------------------------
# Main splitting logic
# ---------------------------------------------------------------------------

CV_PERCENTILE = 8.0   # inner 84 % of each species' distribution per axis


def split_gel(input_path, output_stem=None, polymer_stem=None, solvent_stem=None,
              cv_percentile=CV_PERCENTILE):
    """
    Split a slab-with-support data file into polymer-only and solvent-only files.

    An 8 % percentile trim is applied independently to each species so that the
    output files sample the homogeneous interior — avoiding the solvent-rich shell
    around the finite gel and the polymer-depleted periphery.

    Parameters
    ----------
    input_path    : path to the combined .data file
    output_stem   : base stem for both outputs (overridden per-species below)
    polymer_stem  : override stem for the polymer_only file (e.g. different directory)
    solvent_stem  : override stem for the solvent_only file (e.g. different directory)
    cv_percentile : percentile inset applied to each species (default 8.0)

    If polymer_stem / solvent_stem are None, output_stem is used for both.
    If output_stem is also None, defaults to input path without extension.
    """
    print(f"Reading: {input_path}")
    print(f"CV trim: {cv_percentile}% per face per species")
    data = parse_data_file(input_path)

    header = data['header']
    atoms  = data['atoms']
    bonds  = data['bonds']

    if output_stem is None:
        output_stem = os.path.splitext(input_path)[0]

    # Per-species stems (fall back to shared output_stem if not specified)
    if polymer_stem is None:
        polymer_stem = output_stem
    if solvent_stem is None:
        solvent_stem = output_stem

    orig_masses = {t: (m, c) for t, m, c in data['masses']}

    # -----------------------------------------------------------------------
    # 1.  POLYMER-ONLY  (types 1 and 2; no remapping)
    # -----------------------------------------------------------------------
    POLYMER_TYPES = {1, 2}

    poly_atoms_raw  = [a for a in atoms if a['type'] in POLYMER_TYPES]

    # Trim to inner CV before renumbering so bonds to trimmed atoms are dropped
    poly_atoms_raw, poly_box = trim_to_cv(poly_atoms_raw, cv_percentile, "polymer")

    poly_ids_old    = {a['id'] for a in poly_atoms_raw}
    poly_bonds_raw  = [b for b in bonds
                       if b['atom1'] in poly_ids_old and b['atom2'] in poly_ids_old]

    poly_id_map = {a['id']: new_id for new_id, a in enumerate(poly_atoms_raw, start=1)}

    poly_atoms = [{
        'id':   poly_id_map[a['id']],
        'mol':  a['mol'],
        'type': a['type'],
        'x': a['x'], 'y': a['y'], 'z': a['z'],
    } for a in poly_atoms_raw]

    poly_bonds = [{
        'id':    new_bid,
        'type':  b['type'],
        'atom1': poly_id_map[b['atom1']],
        'atom2': poly_id_map[b['atom2']],
    } for new_bid, b in enumerate(poly_bonds_raw, start=1)]

    poly_type_info = [(t, *orig_masses.get(t, (1.0, ''))) for t in sorted(POLYMER_TYPES)]

    poly_out = f"{polymer_stem}_polymer_only.data"
    write_data_file(
        poly_out,
        title="LAMMPS data file — polymer only (crosslinks + chain beads)",
        box=poly_box,
        atom_types_info=poly_type_info,
        atoms=poly_atoms,
        bonds=poly_bonds,
    )
    print(f"Wrote: {poly_out}  ({len(poly_atoms)} atoms, {len(poly_bonds)} bonds)")

    # -----------------------------------------------------------------------
    # 2.  SOLVENT-ONLY  (type 3 only; support/piston already stripped by isolate_gel)
    # -----------------------------------------------------------------------
    SOLVENT_TYPES   = {3, 4, 5}
    TYPE_REMAP      = {3: 1, 4: 2, 5: 3}
    TYPE_LABELS_NEW = {
        1: 'Solvent (was type 3)',
        2: 'Bottom support — frozen (was type 4)',
        3: 'Top piston — mobile (was type 5)',
    }

    solv_atoms_raw = [a for a in atoms if a['type'] in SOLVENT_TYPES]

    # Trim solvent to its own inner CV for homogeneous density
    solv_atoms_raw, solv_box = trim_to_cv(solv_atoms_raw, cv_percentile, "solvent")

    solv_id_map = {a['id']: new_id for new_id, a in enumerate(solv_atoms_raw, start=1)}

    present_orig_types = sorted({a['type'] for a in solv_atoms_raw})
    present_new_types  = sorted({TYPE_REMAP[t] for t in present_orig_types})

    solv_atoms = [{
        'id':   solv_id_map[a['id']],
        'mol':  a['mol'],
        'type': TYPE_REMAP[a['type']],
        'x': a['x'], 'y': a['y'], 'z': a['z'],
    } for a in solv_atoms_raw]

    OLD_FROM_NEW = {v: k for k, v in TYPE_REMAP.items()}
    solv_type_info = [
        (new_t, orig_masses.get(OLD_FROM_NEW[new_t], (1.0, ''))[0], TYPE_LABELS_NEW[new_t])
        for new_t in present_new_types
    ]

    solv_out = f"{solvent_stem}_solvent_only.data"
    write_data_file(
        solv_out,
        title="LAMMPS data file — solvent only",
        box=solv_box,
        atom_types_info=solv_type_info,
        atoms=solv_atoms,
        bonds=None,
    )
    print(f"Wrote: {solv_out}  ({len(solv_atoms)} atoms)")

    return poly_out, solv_out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split a LAMMPS slab-with-support data file into polymer-only "
                    "and solvent-only files."
    )
    parser.add_argument("input", help="Input .data file path")
    parser.add_argument(
        "--output-dir", default=None,
        help="Directory for both output files (default: same directory as input)"
    )
    parser.add_argument(
        "--output-stem", default=None,
        help="Override output filename stem for both files (takes precedence over --output-dir)"
    )
    parser.add_argument(
        "--polymer-dir", default=None,
        help="Directory for the polymer_only file (overrides --output-dir for polymer)"
    )
    parser.add_argument(
        "--solvent-dir", default=None,
        help="Directory for the solvent_only file (overrides --output-dir for solvent)"
    )
    parser.add_argument(
        "--cv-percentile", type=float, default=CV_PERCENTILE,
        help=f"Percentile trim applied independently to each species (default {CV_PERCENTILE}). "
             f"E.g. 8.0 keeps the inner 84%% of each species' distribution."
    )
    args = parser.parse_args()

    basename = os.path.splitext(os.path.basename(args.input))[0]

    # Shared stem (fallback for both species)
    if args.output_stem is not None:
        shared_stem = args.output_stem
    else:
        out_dir = args.output_dir if args.output_dir else os.path.dirname(args.input)
        shared_stem = os.path.join(out_dir, basename)

    # Per-species stems override shared stem when a specific dir is given
    polymer_stem = os.path.join(args.polymer_dir, basename) if args.polymer_dir else None
    solvent_stem = os.path.join(args.solvent_dir, basename) if args.solvent_dir else None

    split_gel(args.input, shared_stem, polymer_stem=polymer_stem, solvent_stem=solvent_stem,
              cv_percentile=args.cv_percentile)
