#!/usr/bin/env python3
"""
add_more_plates_to_gel.py
=========================
CONVERTER: an equilibrated slab_with_support snapshot of an ISOLATED gel (a finite
network floating in its solvent bath, not touching its periodic image)  ->  the
six-plate + two-wet-piston data file for compress_slab.lmp   (2026-09-29).

The bulk-modulus counterpart of slab_two_pistons.py, laid out like it.  It replaces
the notebook-only add_more_plates_to_gel code of 2026-07/08, which put the six
plates ON the faces of the periodic box and deleted every polymer bead within
1 sigma of them (the network was cut to fit the box).  Here nothing is cut:

  * the gel is made whole (bond-graph unwrap) and centred;
  * the bath is CROPPED to the margins asked for, so the production run does not
    carry the large equilibration bath (x, y stay periodic; the solvent pairs the
    new periodic seam brings closer than overlap_min are removed);
  * six solvent-TRANSPARENT load plates are PARKED plate_park_gap outside the six
    bounding-box faces of the gel -- compress_slab.lmp seats them onto the gel
    (Phase 1.25) exactly as triaxial_compression_two_pist.lmp seats its load
    piston and support;
  * z is closed by two wet NPT-pistons (feed above, permeate below) with vacuum
    margins, so the bath is held at P_target for the whole run and the solvent the
    gel expels is measured by the piston displacement.

Geometry produced (z from bottom to top; x and y are periodic, bath all round):

    zlo = 0
    | vacuum margin (margin_perm)
    | permeate piston (type 6)                        full-width sheet
    | bath (reservoir_z)
    |   plate z-lo (type 4)   <- the "support" of the two-piston decks
    |   gel, with plates x-lo/x-hi (8/9) and y-lo/y-hi (10/11) on its sides
    |   plate z-hi (type 7)   <- the "load piston" of the two-piston decks
    | bath (reservoir_z)
    | feed piston (type 5)                            full-width sheet
    | vacuum margin (margin_feed)
    zhi

Atom types written (11): 1 crosslink, 2 chain bead, 3 solvent, 4 plate z-lo,
5 feed piston, 6 permeate piston, 7 plate z-hi, 8 plate x-lo, 9 plate x-hi,
10 plate y-lo, 11 plate y-hi.  Types 4/5/6/7 keep the meaning they have in the
two-piston decks (lower loading plate / feed / permeate / upper loading plate), so
lib/triaxial.py reads the z planes of a compress_slab trajectory unchanged.
One type per plate: the deck groups them by type (the old mol-ID lookup returned
EMPTY groups whenever leftover wall atoms carried higher mol IDs than the gel).

Usage (CLI)
-----------
  python3 add_more_plates_to_gel.py --input <final_config_...data> [--output <...>]
        [--bath-margin-xy 12] [--reservoir-z 15] [--plate-park-gap 3]
        [--plate-overhang 6] [--plate-spacing 0.5] [--piston-spacing 0.5]
        [--margin-feed 10] [--margin-perm 10] [--no-log] [--no-png]
  --self-check-only <file> validates an existing output file without rewriting it.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys
from collections import deque
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import slab_two_pistons as s2p  # noqa: E402  (hex sheet, json helper, info log)

try:
    from scipy.spatial import cKDTree
    _SCIPY = True
except ImportError:  # pragma: no cover
    _SCIPY = False

# ---------------------------------------------------------------------------
# Constants (mirror compress_slab.lmp)
# ---------------------------------------------------------------------------
T_CROSSLINK, T_CHAIN, T_SOLVENT = 1, 2, 3
T_FEED, T_PERM = 5, 6
PLATE_TYPE = {'zlo': 4, 'zhi': 7, 'xlo': 8, 'xhi': 9, 'ylo': 10, 'yhi': 11}
PLATE_ORDER = ('zlo', 'zhi', 'xlo', 'xhi', 'ylo', 'yhi')
POLYMER_TYPES = (T_CROSSLINK, T_CHAIN)
WALL_TYPES = (4, 5, 6, 7, 8, 9, 10, 11)
N_TYPES = 11
FENE_R0_MAX = 1.5
AX = {'x': 0, 'y': 1, 'z': 2}
TYPE_NAMES = {1: 'crosslink', 2: 'chain bead', 3: 'solvent', 4: 'plate z-lo', 5: 'feed piston',
              6: 'permeate piston', 7: 'plate z-hi', 8: 'plate x-lo', 9: 'plate x-hi',
              10: 'plate y-lo', 11: 'plate y-hi'}


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
@dataclass
class Config:
    """Every knob of the converter (the notebook's Config cell builds one)."""
    input_file: str
    output_file: str = None            # None -> <input stem>_six_plates.data next to the input
    bath_margin_xy: float = 12.0       # sigma of bath between a PARKED side plate and the periodic box face
                                       # (the gel sees its image across 2*(park gap + margin))
    reservoir_z: float = 15.0          # sigma of bath between a parked z plate and its wet piston
    plate_park_gap: float = 3.0        # plates are parked this far outside the gel bounding-box faces
    plate_overhang: float = 6.0        # each plate extends this far beyond the gel BB in its two in-plane directions
    plate_spacing: float = 0.5         # square-lattice spacing of the load plates (snapped to the plate size)
    piston_spacing: float = 0.5        # hex spacing of the wet pistons (snapped to tile lx, ly)
    piston_clearance: float = 1.0      # each wet piston plane this far outside the outermost solvent bead
    margin_feed: float = 10.0          # vacuum above the feed piston (the expelled solvent raises it)
    margin_perm: float = 10.0          # vacuum below the permeate piston
    overlap_min: float = 0.8           # solvent closer than this to another bead across the new seam is removed
    drop_types: tuple = (4, 5)         # wall types of the equilibration run, deleted (old support / piston)
    drop_fragments: bool = True        # delete polymer not bonded to the main network (free corner beads of the lattice)
    density_window: float = 4.0        # bath-density window starts this far outside the gel BB
    log_info: bool = True              # append the entry to slab_data_file_info.md
    make_png: bool = True

    def __post_init__(self):
        if self.output_file is None:
            p = Path(self.input_file)
            self.output_file = str(p.with_name(p.stem + '_six_plates.data'))
        self.drop_types = tuple(int(t) for t in self.drop_types)
        assert self.plate_overhang < self.plate_park_gap + min(self.bath_margin_xy, self.reservoir_z), \
            'plate_overhang must be smaller than plate_park_gap + the bath margin (plates must stay inside the box)'


# ---------------------------------------------------------------------------
# Fast array parser (the snapshot of a 12x12x12 gel holds ~1.5 M atoms)
# ---------------------------------------------------------------------------
def read_data(path):
    """LAMMPS molecular data file -> dict(ids, mol, typ, xyz, vel, bonds, box, n_types).
    Image flags, when present, are ignored: the gel is unwrapped over its bond graph."""
    with open(path) as f:
        lines = f.read().split('\n')
    box, n_atoms, n_bonds, n_types = {}, 0, 0, 0
    sec = {}
    for i, raw in enumerate(lines):
        t = raw.strip()
        if not t:
            continue
        if t.endswith(' atoms'):
            n_atoms = int(t.split()[0])
        elif t.endswith(' bonds'):
            n_bonds = int(t.split()[0])
        elif t.endswith('atom types'):
            n_types = int(t.split()[0])
        elif 'xlo xhi' in t or 'ylo yhi' in t or 'zlo zhi' in t:
            p = t.split()
            box[p[2]], box[p[3]] = float(p[0]), float(p[1])
        elif t.split()[0] in ('Atoms', 'Velocities', 'Bonds', 'Masses'):
            sec[t.split()[0]] = i

    def block(name, n, ncol, dtype):
        if name not in sec or n == 0:
            return None
        rows = [l.split('#')[0].split()[:ncol] for l in lines[sec[name] + 2: sec[name] + 2 + n]]
        return np.array(rows, dtype=dtype)

    A = block('Atoms', n_atoms, 6, float)
    if A is None or len(A) != n_atoms:
        raise ValueError(f'{path}: could not read the Atoms section ({n_atoms} atoms declared)')
    ids = A[:, 0].astype(np.int64)
    out = dict(ids=ids, mol=A[:, 1].astype(np.int64), typ=A[:, 2].astype(np.int64), xyz=A[:, 3:6].copy(),
               box=box, n_types=n_types)
    V = block('Velocities', n_atoms, 4, float)
    vel = np.zeros((n_atoms, 3))
    if V is not None and len(V) == n_atoms:
        order = np.argsort(ids)
        pos = order[np.searchsorted(ids[order], V[:, 0].astype(np.int64))]
        vel[pos] = V[:, 1:4]
    out['vel'] = vel
    out['has_vel'] = V is not None
    B = block('Bonds', n_bonds, 4, np.int64)
    out['bonds'] = B if B is not None else np.zeros((0, 4), np.int64)
    return out


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def unwrap_gel(xyz, is_poly, bonds_idx, L):
    """Make the network whole: breadth-first walk over the bond graph, each bonded
    neighbour placed in the image nearest its parent.  Returns (xyz, n_crossing, comp):
    n_crossing = bonds that crossed a periodic boundary in the input, comp = bond-graph
    component of every polymer bead (-1 for the rest).  Raises if
    the network closes on itself through a boundary (a gel that reaches its own
    periodic image is not isolated and cannot be boxed by six plates)."""
    xyz = xyz.copy()
    n = len(xyz)
    adj = [[] for _ in range(n)]
    for a, b in bonds_idx:
        adj[a].append(b)
        adj[b].append(a)
    d0 = xyz[bonds_idx[:, 1]] - xyz[bonds_idx[:, 0]]
    n_cross = int(np.any(np.abs(d0) > 0.5 * L[None, :], axis=1).sum())
    seen = np.zeros(n, bool)
    comp = -np.ones(n, np.int64)
    poly = np.flatnonzero(is_poly)
    n_comp = 0
    for seed in poly:
        if seen[seed]:
            continue
        seen[seed] = True
        comp[seed] = n_comp
        q = deque([seed])
        while q:
            c = q.popleft()
            for nb in adj[c]:
                if seen[nb]:
                    continue
                xyz[nb] -= L * np.round((xyz[nb] - xyz[c]) / L)
                seen[nb] = True
                comp[nb] = n_comp
                q.append(nb)
        n_comp += 1
    d = np.linalg.norm(xyz[bonds_idx[:, 1]] - xyz[bonds_idx[:, 0]], axis=1)
    if len(d) and d.max() > FENE_R0_MAX:
        raise ValueError(f'{int((d > FENE_R0_MAX).sum())} bond(s) are longer than {FENE_R0_MAX} sigma after the unwrap '
                         f'(longest {d.max():.2f}): the network closes on itself through a periodic boundary, i.e. the '
                         'gel has reached its own image.  This converter needs an ISOLATED gel.')
    return xyz, n_cross, comp


def plate_lattice(lo_a, hi_a, lo_b, hi_b, spacing):
    """Square lattice over [lo_a, hi_a] x [lo_b, hi_b], spacing snapped so the rows
    fill the plate edge to edge.  Returns (n, 2) in-plane coordinates."""
    na = max(2, int(round((hi_a - lo_a) / spacing)) + 1)
    nb = max(2, int(round((hi_b - lo_b) / spacing)) + 1)
    A, B = np.meshgrid(np.linspace(lo_a, hi_a, na), np.linspace(lo_b, hi_b, nb), indexing='ij')
    return np.column_stack([A.ravel(), B.ravel()])


def make_plate(face, bb_lo, bb_hi, cfg):
    """Beads of one parked load plate.  face = 'zhi' etc.; bb_lo/bb_hi = gel bounding box."""
    n = AX[face[0]]
    a, b = [k for k in range(3) if k != n]
    plane = (bb_hi[n] + cfg.plate_park_gap) if face.endswith('hi') else (bb_lo[n] - cfg.plate_park_gap)
    ab = plate_lattice(bb_lo[a] - cfg.plate_overhang, bb_hi[a] + cfg.plate_overhang,
                       bb_lo[b] - cfg.plate_overhang, bb_hi[b] + cfg.plate_overhang, cfg.plate_spacing)
    p = np.empty((len(ab), 3))
    p[:, n] = plane
    p[:, a], p[:, b] = ab[:, 0], ab[:, 1]
    return p, plane


def seam_overlaps(xyz, typ, lx, ly, r):
    """Indices of the solvent beads to delete so that no pair involving a solvent
    bead is closer than r THROUGH the new x,y periodicity (z is finite).  Only pairs
    the new seam creates count -- pairs already that close inside the box are the
    input's own (e.g. the bonded neighbours of an unrelaxed lattice) and are left
    alone.  Greedy: of a solvent-solvent pair the second bead goes; of a
    solvent-polymer pair the solvent goes.  Returns (indices, n_seam_pairs, n_pp),
    n_pp = polymer-polymer pairs across the seam (the margin is then too small)."""
    if not _SCIPY:
        raise RuntimeError('scipy is required for the seam overlap removal')
    p = xyz.copy()
    p[:, 0] = np.mod(p[:, 0], lx)
    p[:, 1] = np.mod(p[:, 1], ly)
    p[:, 2] = p[:, 2] - p[:, 2].min() + 1.0e3
    pairs = cKDTree(p, boxsize=[lx, ly, 1.0e6]).query_pairs(r, output_type='ndarray')
    if len(pairs):
        direct = np.linalg.norm(p[pairs[:, 0]] - p[pairs[:, 1]], axis=1)
        pairs = pairs[direct >= r]
    drop = np.zeros(len(typ), bool)
    n_pp = 0
    for i, j in pairs:
        si, sj = typ[i] == T_SOLVENT, typ[j] == T_SOLVENT
        if not (si or sj):
            n_pp += 1
            continue
        if drop[i] or drop[j]:
            continue
        drop[j if sj else i] = True
    return np.flatnonzero(drop), len(pairs), n_pp


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------
def write_data(path, header_line, ids, mols, typ, xyz, vel, bonds, box):
    """LAMMPS data file: Atoms (molecular, image flags 0 0 0), Velocities, Bonds."""
    header_line = ' '.join(str(header_line).split())
    with open(path, 'w') as f:
        f.write(header_line + '\n\n')
        f.write(f'{len(ids)} atoms\n{len(bonds)} bonds\n0 angles\n0 dihedrals\n0 impropers\n\n')
        f.write(f'{N_TYPES} atom types\n1 bond types\n\n')
        for d in 'xyz':
            f.write(f"{box[d + 'lo']:.10f} {box[d + 'hi']:.10f} {d}lo {d}hi\n")
        f.write('\nMasses\n\n')
        for t in range(1, N_TYPES + 1):
            f.write(f'{t} 1.0  # {TYPE_NAMES[t]}' + ('  (deck overrides with piston_mass)' if t in (T_FEED, T_PERM) else '') + '\n')
        f.write('\nAtoms # molecular\n\n')
        np.savetxt(f, np.column_stack([ids, mols, typ, xyz, np.zeros((len(ids), 3))]),
                   fmt='%d %d %d %.10f %.10f %.10f %d %d %d')
        f.write('\nVelocities\n\n')
        np.savetxt(f, np.column_stack([ids, vel]), fmt='%d %.10g %.10g %.10g')
        if len(bonds):
            f.write('\nBonds\n\n')
            np.savetxt(f, np.column_stack([np.arange(1, len(bonds) + 1), bonds]), fmt='%d %d %d %d')


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------
def convert(cfg: Config):
    """Run the whole conversion.  Returns a dict `info` (geometry, counts, paths)."""
    t_start = _dt.datetime.now()
    inp = Path(cfg.input_file)
    if not inp.exists():
        raise FileNotFoundError(f'input data file not found: {inp}')
    print('=' * 78)
    print(f'add_more_plates_to_gel converter  ({t_start:%Y-%m-%d %H:%M})')
    print(f'  input : {inp}')
    print(f'  output: {cfg.output_file}')
    print('=' * 78)

    # ---- (a) parse, drop the equilibration walls ----------------------------
    D = read_data(str(inp))
    box = D['box']
    L = np.array([box['xhi'] - box['xlo'], box['yhi'] - box['ylo'], box['zhi'] - box['zlo']])
    typ_in = D['typ']
    counts_in = {int(t): int((typ_in == t).sum()) for t in np.unique(typ_in)}
    print(f'  read {len(typ_in)} atoms, {len(D["bonds"])} bonds; box {L[0]:.3f} x {L[1]:.3f} x {L[2]:.3f}')
    print('  input type counts: ' + ', '.join(f'{t}={n}' for t, n in counts_in.items()))
    keep = np.isin(typ_in, (1, 2, 3))
    n_dropped = int((~keep).sum())
    stray = sorted(set(np.unique(typ_in[~keep]).tolist()) - set(cfg.drop_types))
    if stray:
        raise ValueError(f'input holds atom types {stray} that are neither gel (1,2,3) nor in drop_types {cfg.drop_types}')
    ids, mol, typ, xyz, vel = D['ids'][keep], D['mol'][keep], typ_in[keep], D['xyz'][keep], D['vel'][keep]
    print(f'  deleted {n_dropped} wall beads of the equilibration run (types {list(cfg.drop_types)})')
    index_of = {int(a): i for i, a in enumerate(ids)}
    try:
        bonds_idx = np.array([[index_of[int(a)], index_of[int(b)]] for a, b in D['bonds'][:, 2:4]], dtype=np.int64)
    except KeyError as e:
        raise ValueError(f'a bond references deleted atom {e} -- polymer must not be in drop_types') from e
    btype = D['bonds'][:, 1]

    # ---- (b) make the gel whole, fold the solvent around it --------------------
    is_poly = np.isin(typ, POLYMER_TYPES)
    is_solv = typ == T_SOLVENT
    xyz, n_cross, comp = unwrap_gel(xyz, is_poly, bonds_idx, L)
    # ---- free fragments (sol fraction) -------------------------------------------
    # slab_with_support.ipynb placed (until 2026-09-29) a crosslink bead on every lattice
    # site inside the cube, including the corner sites no chain reaches: those beads are
    # free monomers that wander through the bath.  They are polymer by TYPE, so they would
    # set the bounding box the plates are parked and seated on; only the main network is kept.
    sizes = np.bincount(comp[is_poly])
    main = int(np.argmax(sizes))
    frag = is_poly & (comp != main)
    n_frag_beads, n_frag = int(frag.sum()), int(len(sizes) - 1)
    if n_frag:
        if not cfg.drop_fragments:
            raise ValueError(f'{n_frag} free polymer fragment(s) ({n_frag_beads} beads) are not bonded to the main network; '
                             'the plates would seat on them.  Set drop_fragments=True.')
        keepm = ~frag
        sel = np.flatnonzero(keepm)
        remap = -np.ones(len(typ), np.int64)
        remap[sel] = np.arange(len(sel))
        ids, mol, typ, xyz, vel = ids[sel], mol[sel], typ[sel], xyz[sel], vel[sel]
        kb = keepm[bonds_idx[:, 0]] & keepm[bonds_idx[:, 1]]
        bonds_idx, btype = remap[bonds_idx[kb]], btype[kb]
        is_poly, is_solv = np.isin(typ, POLYMER_TYPES), typ == T_SOLVENT
        print(f'  removed {n_frag} free fragment(s), {n_frag_beads} polymer bead(s) not bonded to the main network '
              f'(fragment sizes: {sorted(np.delete(sizes, main).tolist(), reverse=True)[:8]})')
    com = xyz[is_poly].mean(axis=0)
    xyz[is_solv] -= L * np.round((xyz[is_solv] - com) / L)
    bb_lo, bb_hi = xyz[is_poly].min(axis=0), xyz[is_poly].max(axis=0)
    L_bb = bb_hi - bb_lo
    L_rg = 2.0 * np.sqrt(3.0 * xyz[is_poly].var(axis=0))
    ctr = 0.5 * (bb_lo + bb_hi)
    print(f'  gel made whole ({n_cross} bonds crossed a boundary in the input); COM ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})')
    print(f'  gel bounding box  L_bb = {L_bb[0]:.2f} x {L_bb[1]:.2f} x {L_bb[2]:.2f}   '
          f'Rg-based L_rg = {L_rg[0]:.2f} x {L_rg[1]:.2f} x {L_rg[2]:.2f}')
    gap_img = L - L_bb
    print(f'  distance of the gel to its own image in the input: x {gap_img[0]:.2f}  y {gap_img[1]:.2f}  z {gap_img[2]:.2f} sigma')
    if gap_img.min() < 2.0 * FENE_R0_MAX:
        print('  WARNING: the gel is within 3 sigma of its periodic image -- it was NOT swelling freely in the input')

    # ---- (c) crop the bath ------------------------------------------------------
    half_want = 0.5 * L_bb + cfg.plate_park_gap + np.array([cfg.bath_margin_xy, cfg.bath_margin_xy, cfg.reservoir_z])
    half = np.minimum(half_want, 0.5 * L)
    cropped = half_want < 0.5 * L - 1e-9
    for k, d in enumerate('xyz'):
        if not cropped[k]:
            print(f'  NOTE: the input box is too small in {d} for the margin asked for '
                  f'({2 * half_want[k]:.2f} > {L[k]:.2f}): keeping the whole input length '
                  f'(margin {0.5 * (L[k] - L_bb[k]) - cfg.plate_park_gap:.2f} sigma instead)')
    # solvent is folded around the COM, the crop is centred on the BB centre: refold around it
    xyz[is_solv] -= L * np.round((xyz[is_solv] - ctr) / L)
    inside = np.all(np.abs(xyz - ctr) < half, axis=1) | is_poly
    n_crop = int((~inside).sum())
    rho_in = float(is_solv.sum()) / float(np.prod(L))
    sel = np.flatnonzero(inside)
    remap = -np.ones(len(typ), np.int64)
    remap[sel] = np.arange(len(sel))
    ids, mol, typ, xyz, vel = ids[sel], mol[sel], typ[sel], xyz[sel], vel[sel]
    bonds_idx = remap[bonds_idx]
    assert bonds_idx.min() >= 0, 'a polymer bead was cropped -- this cannot happen'
    is_poly, is_solv = np.isin(typ, POLYMER_TYPES), typ == T_SOLVENT
    print(f'  bath cropped to {2 * half[0]:.2f} x {2 * half[1]:.2f} x {2 * half[2]:.2f}: removed {n_crop} solvent beads '
          f'({int(is_solv.sum())} kept)')
    xyz = xyz - (ctr - half)                       # box origin at 0: gel centred, x,y in [0, 2 half)
    bb_lo, bb_hi = bb_lo - (ctr - half), bb_hi - (ctr - half)
    lx, ly = 2.0 * half[0], 2.0 * half[1]

    # ---- (d) the new periodic seam --------------------------------------------
    n_seam = 0
    if cropped[0] or cropped[1]:
        kill, n_pairs, n_pp = seam_overlaps(xyz, typ, lx, ly, cfg.overlap_min)
        if n_pp:
            raise RuntimeError(f'{n_pp} polymer-polymer pair(s) closer than {cfg.overlap_min} sigma across the new seam: '
                               'bath_margin_xy is too small')
        n_seam = len(kill)
        m = np.ones(len(typ), bool)
        m[kill] = False
        sel = np.flatnonzero(m)
        remap = -np.ones(len(typ), np.int64)
        remap[sel] = np.arange(len(sel))
        ids, mol, typ, xyz, vel = ids[sel], mol[sel], typ[sel], xyz[sel], vel[sel]
        bonds_idx = remap[bonds_idx]
        is_poly, is_solv = np.isin(typ, POLYMER_TYPES), typ == T_SOLVENT
        print(f'  new x,y seam: {n_pairs} pair(s) closer than {cfg.overlap_min} sigma -> removed {n_seam} solvent beads')

    # ---- bath density away from the gel --------------------------------------------
    out_gel = is_solv & np.any((xyz < bb_lo - cfg.density_window) | (xyz > bb_hi + cfg.density_window), axis=1)
    v_tot = float(np.prod(2.0 * half))
    v_in = float(np.prod(np.minimum(L_bb + 2.0 * cfg.density_window, 2.0 * half)))
    rho_bulk = float(out_gel.sum()) / max(v_tot - v_in, 1e-9)
    print(f'  bath number density ({cfg.density_window:g} sigma clear of the gel BB): rho = {rho_bulk:.4f} /sigma^3 '
          f'(input box average incl. the gel: {rho_in:.4f})')

    # ---- (e) wet pistons + six plates ----------------------------------------------
    sz = xyz[is_solv, 2]
    z_perm, z_feed = float(sz.min()) - cfg.piston_clearance, float(sz.max()) + cfg.piston_clearance
    sheet, (nx, ny, sdx, sdy) = s2p.make_sheet_hex(lx, ly, cfg.piston_spacing)
    new_pos, new_typ, new_mol = [], [], []
    mol_next = int(mol.max()) + 1
    for zp, t in ((z_feed, T_FEED), (z_perm, T_PERM)):
        new_pos.append(np.column_stack([sheet, np.full(len(sheet), zp)]))
        new_typ.append(np.full(len(sheet), t))
        new_mol.append(np.full(len(sheet), mol_next))
        mol_next += 1
    print(f'  wet pistons: hex sheet {nx} x {ny} = {len(sheet)} beads each (dx={sdx:.4f}, row pitch={sdy:.4f}); '
          f'permeate z={z_perm:.3f}, feed z={z_feed:.3f}')
    planes, n_plate = {}, {}
    for face in PLATE_ORDER:
        p, plane = make_plate(face, bb_lo, bb_hi, cfg)
        planes[face], n_plate[face] = float(plane), len(p)
        new_pos.append(p)
        new_typ.append(np.full(len(p), PLATE_TYPE[face]))
        new_mol.append(np.full(len(p), mol_next))
        mol_next += 1
    print('  load plates (parked): ' + '  '.join(f'{f} @ {planes[f]:.2f} ({n_plate[f]})' for f in PLATE_ORDER))

    n_old = len(typ)
    all_xyz = np.vstack([xyz] + new_pos)
    all_typ = np.concatenate([typ] + new_typ).astype(np.int64)
    all_mol = np.concatenate([mol] + new_mol).astype(np.int64)
    all_vel = np.vstack([vel, np.zeros((len(all_typ) - n_old, 3))])
    zlo_new = z_perm - cfg.margin_perm
    all_xyz[:, 2] -= zlo_new
    new_box = dict(xlo=0.0, xhi=lx, ylo=0.0, yhi=ly, zlo=0.0, zhi=z_feed + cfg.margin_feed - zlo_new)
    zs = lambda z: z - zlo_new      # noqa: E731
    bb_lo_o = np.array([bb_lo[0], bb_lo[1], zs(bb_lo[2])])
    bb_hi_o = np.array([bb_hi[0], bb_hi[1], zs(bb_hi[2])])
    planes_o = {f: (zs(v) if f[0] == 'z' else v) for f, v in planes.items()}
    wet = np.isin(all_typ, (T_FEED, T_PERM))
    all_xyz[wet, 0] = np.mod(all_xyz[wet, 0], lx)
    all_xyz[wet, 1] = np.mod(all_xyz[wet, 1], ly)
    mob = all_typ <= 3
    all_xyz[mob, 0] = np.mod(all_xyz[mob, 0], lx)       # only solvent can sit on the far face after the crop
    all_xyz[mob, 1] = np.mod(all_xyz[mob, 1], ly)
    bonds = np.column_stack([btype, bonds_idx + 1])
    new_ids = np.arange(1, len(all_typ) + 1)

    A = lx * ly
    V_bb = float(np.prod(L_bb))
    info = dict(
        input_file=str(inp), output_file=str(cfg.output_file), date=t_start.strftime('%Y-%m-%d %H:%M'),
        config=asdict(cfg), box_in=dict(box), box_out=new_box, lx=lx, ly=ly, area=A,
        gel_bb_lo=bb_lo_o.tolist(), gel_bb_hi=bb_hi_o.tolist(), L_bb=L_bb.tolist(), L_rg=L_rg.tolist(),
        V_bb=V_bb, V_rg=float(np.prod(L_rg)), image_gap_in=gap_img.tolist(),
        plates=planes_o, n_plate=n_plate, z_feed_piston=zs(z_feed), z_perm_piston=zs(z_perm),
        n_sheet=len(sheet), rho_bulk=rho_bulk, n_fragments_removed=n_frag, n_fragment_beads_removed=n_frag_beads, n_cropped=n_crop, n_seam_removed=n_seam, n_walls_dropped=n_dropped,
        bath_margin_xy_actual=float(half[0] - 0.5 * L_bb[0] - cfg.plate_park_gap),
        reservoir_z_actual=float(half[2] - 0.5 * L_bb[2] - cfg.plate_park_gap),
        n_atoms=int(len(all_typ)), n_bonds=int(len(bonds)),
        counts={int(t): int((all_typ == t).sum()) for t in range(1, N_TYPES + 1)},
        # the wet pistons back out by (expelled solvent)/A: the largest VOLUMETRIC strain the two
        # vacuum margins can absorb before a piston reaches a box face (2 sigma guard each)
        max_strain_vol_margins=float((cfg.margin_feed + cfg.margin_perm - 4.0) * A / V_bb),
    )

    # ---- (f) validate BEFORE writing ------------------------------------------------
    val = validate(all_typ, all_xyz, bonds, new_box, n_old, cfg.overlap_min)
    info['validation'] = val
    if not val['ok']:
        raise RuntimeError('validation FAILED -- file not written:\n  ' + '\n  '.join(val['problems']))

    cfgs = ' '.join(f'{k}={v}' for k, v in asdict(cfg).items() if k not in ('input_file', 'output_file'))
    header = (f'six-plate + two-wet-piston data file from {inp.name} via add_more_plates_to_gel.py on {info["date"]}; '
              'types 1 xl 2 chain 3 solv 4 plate_zlo 5 feed_piston 6 perm_piston 7 plate_zhi 8 plate_xlo 9 plate_xhi '
              f'10 plate_ylo 11 plate_yhi; {cfgs}; gel_bb_lo={bb_lo_o.round(3).tolist()} gel_bb_hi={bb_hi_o.round(3).tolist()} '
              f'rho_bulk={rho_bulk:.4f}')
    write_data(cfg.output_file, header, new_ids, all_mol, all_typ, all_xyz, all_vel, bonds, new_box)
    with open(str(cfg.output_file) + '.info.json', 'w') as f:
        json.dump(s2p._jsonable(info), f, indent=1)
    print(f'  wrote {len(all_typ)} atoms, {len(bonds)} bonds -> {cfg.output_file}')
    print(f'  sidecar: {cfg.output_file}.info.json')
    print_layout(info)
    if cfg.make_png:
        info['png'] = plot_profiles(all_typ, all_xyz, new_box, info, str(cfg.output_file).replace('.data', '_profiles.png'))
    entry = info_log_entry(info)
    print('\n--- slab_data_file_info.md entry ---\n' + entry)
    if cfg.log_info:
        s2p.append_info_log(entry)
    print(f'\nconversion done in {(_dt.datetime.now() - t_start).total_seconds():.1f} s')
    return info


# ---------------------------------------------------------------------------
# Validation (also the CLI self-check)
# ---------------------------------------------------------------------------
def validate(typ, xyz, bonds, box, n_old=0, overlap_min=0.8, verbose=True):
    """Checks: (1) no INTERACTING pair closer than overlap_min that the CONVERSION made --
    a pair with a new wall bead (index >= n_old) or a pair of mobile beads that is close
    only through the new x,y periodicity; solvent-plate and wall-wall pairs are off in the
    deck and exempt, and pairs already that close in the input are reported, not failed
    (the deck's Phase 0 minimises); (2) every bond shorter than the
    FENE limit with NO periodic wrap (the gel is whole); (3) all eleven types present;
    (4) solvent and polymer between the wet pistons; (5) every plate inside the box
    and outside the gel; (6) no atom outside the box.  Returns dict(ok, problems, ...)."""
    say = print if verbose else (lambda *a, **k: None)
    problems = []
    say('\n--- validation ---')
    lx, ly = box['xhi'] - box['xlo'], box['yhi'] - box['ylo']
    n_close, worst = 0, np.inf
    if _SCIPY:
        p = xyz.copy()
        p[:, 0] = np.mod(p[:, 0] - box['xlo'], lx)
        p[:, 1] = np.mod(p[:, 1] - box['ylo'], ly)
        p[:, 2] = p[:, 2] - box['zlo'] + 1.0e3
        # plate beads that overhang the box wrap onto the far side; they interact with polymer only
        pairs = cKDTree(p, boxsize=[lx, ly, 1.0e6]).query_pairs(overlap_min, output_type='ndarray')
        n_pre = 0
        if len(pairs):
            ta, tb = typ[pairs[:, 0]], typ[pairs[:, 1]]
            wall = np.isin(ta, WALL_TYPES) & np.isin(tb, WALL_TYPES)
            plate = (4, 7, 8, 9, 10, 11)
            solv_plate = ((ta == T_SOLVENT) & np.isin(tb, plate)) | ((tb == T_SOLVENT) & np.isin(ta, plate))
            is_new = (pairs[:, 0] >= n_old) | (pairs[:, 1] >= n_old)
            seam = np.linalg.norm(p[pairs[:, 0]] - p[pairs[:, 1]], axis=1) >= overlap_min
            rel = ~wall & ~solv_plate & (is_new | seam)
            n_pre = int((~wall & ~solv_plate & ~is_new & ~seam).sum())
            n_close = int(rel.sum())
            if n_close:
                d = np.linalg.norm(p[pairs[rel, 0]] - p[pairs[rel, 1]], axis=1)
                d = np.minimum(d, 1e9)
                worst = float(d.min())
                problems.append(f'{n_close} interacting pair(s) closer than {overlap_min} sigma (closest {worst:.3f}); e.g. types '
                                + ', '.join(f'{typ[i]}-{typ[j]}' for i, j in pairs[rel][:5]))
        say(f'  interacting pairs closer than {overlap_min} made by the conversion: {n_close}'
            + (f'   ({n_pre} already in the input, left alone)' if n_pre else ''))
    else:
        say('  (scipy missing -> overlap check skipped)')
    if len(bonds):
        d = np.linalg.norm(xyz[bonds[:, 2] - 1] - xyz[bonds[:, 1] - 1], axis=1)
        say(f'  bonds: longest {d.max():.3f} sigma (FENE limit {FENE_R0_MAX}); none may wrap')
        if d.max() > FENE_R0_MAX:
            problems.append(f'{int((d > FENE_R0_MAX).sum())} bond(s) longer than {FENE_R0_MAX} sigma without wrapping')
    counts = {int(t): int((typ == t).sum()) for t in range(1, N_TYPES + 1)}
    say('  type counts: ' + ', '.join(f'{t} {TYPE_NAMES[t]}={n}' for t, n in counts.items()))
    for t in range(1, N_TYPES + 1):
        if counts[t] == 0:
            problems.append(f'no atoms of type {t} ({TYPE_NAMES[t]})')
    if counts[T_FEED] and counts[T_PERM]:
        z_feed, z_perm = float(xyz[typ == T_FEED, 2].mean()), float(xyz[typ == T_PERM, 2].mean())
        mz = xyz[typ <= 3, 2]
        n_out = int(((mz < z_perm) | (mz > z_feed)).sum())
        say(f'  mobile z in [{mz.min():.3f}, {mz.max():.3f}] vs wet pistons [{z_perm:.3f}, {z_feed:.3f}]: {n_out} outside')
        if n_out:
            problems.append(f'{n_out} mobile beads outside [permeate piston, feed piston]')
        pl = np.isin(typ, (4, 7, 8, 9, 10, 11))
        if pl.any():
            pz = xyz[pl, 2]
            if pz.min() <= z_perm or pz.max() >= z_feed:
                problems.append('a load plate reaches past a wet piston: reservoir_z too small for plate_overhang')
            px, py = xyz[pl, 0], xyz[pl, 1]
            if px.min() < box['xlo'] or px.max() > box['xhi'] or py.min() < box['ylo'] or py.max() > box['yhi']:
                problems.append('a load plate reaches past the periodic box: bath_margin_xy too small for plate_overhang')
        poly = np.isin(typ, POLYMER_TYPES)
        if poly.any() and pl.any():
            lo, hi = xyz[poly].min(axis=0), xyz[poly].max(axis=0)
            for face, t in PLATE_TYPE.items():
                k = AX[face[0]]
                c = float(xyz[typ == t, k].mean())
                gap = (c - hi[k]) if face.endswith('hi') else (lo[k] - c)
                if gap < overlap_min:
                    problems.append(f'plate {face} is {gap:.2f} sigma from the gel face (parked plates must be outside the gel)')
    out = ((xyz[:, 2] < box['zlo']) | (xyz[:, 2] > box['zhi']))
    n_box = int(out.sum())
    say(f'  atoms outside the box in z: {n_box}')
    if n_box:
        problems.append(f'{n_box} atoms outside the box in z')
    ok = not problems
    say('  RESULT: ' + ('OK' if ok else 'PROBLEMS:\n    ' + '\n    '.join(problems)))
    return dict(ok=ok, problems=problems, n_close=n_close, closest=worst, counts=counts, n_outside_box=n_box)


def print_layout(info):
    print('\n--- layout (output coordinates, box origin at 0) ---')
    b = info['box_out']
    print(f"  box  {b['xhi']:.3f} x {b['yhi']:.3f} x {b['zhi']:.3f}   (x, y periodic bath; z closed by the wet pistons)")
    lo, hi = info['gel_bb_lo'], info['gel_bb_hi']
    for k, d in enumerate('xyz'):
        print(f"  {d}: plate {d}lo {info['plates'][d + 'lo']:8.3f} | gel BB [{lo[k]:8.3f}, {hi[k]:8.3f}] "
              f"(L_bb {info['L_bb'][k]:.2f}, L_rg {info['L_rg'][k]:.2f}) | plate {d}hi {info['plates'][d + 'hi']:8.3f}")
    print(f"  z: permeate piston {info['z_perm_piston']:.3f} | feed piston {info['z_feed_piston']:.3f} | "
          f"bath margin x,y {info['bath_margin_xy_actual']:.2f}, reservoir z {info['reservoir_z_actual']:.2f} sigma")
    print(f"  the vacuum margins absorb the solvent expelled up to a VOLUMETRIC strain of ~{info['max_strain_vol_margins']:.2f}; "
          'raise margin_feed / margin_perm for deeper sweeps')


def plot_profiles(typ, xyz, box, info, out_png):
    """Number density along z and along x per species, plate / piston planes marked."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    L = [box['xhi'], box['yhi'], box['zhi']]
    for ax, d in zip(axes, 'zx'):
        k = AX[d]
        a, b = [j for j in range(3) if j != k]
        edges = np.arange(0.0, L[k] + 0.5, 0.5)
        zc = 0.5 * (edges[1:] + edges[:-1])
        # density inside the column through the gel centre (the deck's core column), not the box average
        lo, hi = info['gel_bb_lo'], info['gel_bb_hi']
        col = (xyz[:, a] > lo[a]) & (xyz[:, a] < hi[a]) & (xyz[:, b] > lo[b]) & (xyz[:, b] < hi[b])
        area = (hi[a] - lo[a]) * (hi[b] - lo[b])
        for t, c, lab in (((1, 2), 'tab:orange', 'polymer'), ((3,), 'tab:blue', 'solvent')):
            h, _ = np.histogram(xyz[col & np.isin(typ, t), k], bins=edges)
            ax.plot(zc, h / (area * 0.5), color=c, lw=1.3, label=lab)
        for f in (d + 'lo', d + 'hi'):
            ax.axvline(info['plates'][f], color='k', ls='--', lw=1.5)
        if d == 'z':
            ax.axvline(info['z_feed_piston'], color='tab:green', ls='-', lw=2, label='feed piston')
            ax.axvline(info['z_perm_piston'], color='tab:purple', ls='-', lw=2, label='permeate piston')
        ax.set_xlabel(f'{d} (sigma)')
        ax.set_ylabel('number density in the gel column (1/sigma^3)')
        ax.set_title(f'along {d}  (dashed: parked plates)', fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle(Path(info['output_file']).name, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    print(f'  density profiles -> {out_png}')
    return out_png


def info_log_entry(info):
    c, cfg, b = info['counts'], info['config'], info['box_out']
    n_pl = sum(info['n_plate'].values())
    return '\n'.join([
        f'Converted six-plate + two-wet-piston data file  [{info["date"]}, add_more_plates_to_gel.py]:',
        f'  Input: {Path(info["input_file"]).name}'
        + (f'   ({info["n_fragment_beads_removed"]} free polymer bead(s) not bonded to the network removed)'
           if info['n_fragment_beads_removed'] else ''),
        f"  Box: {b['xhi']:.4f} x {b['yhi']:.4f} x {b['zhi']:.3f}   (bath cropped: x,y margin {info['bath_margin_xy_actual']:.1f}, "
        f"z reservoir {info['reservoir_z_actual']:.1f} sigma; {info['n_cropped']} solvent cropped, {info['n_seam_removed']} removed at the seam)",
        f"  Gel: L_bb {info['L_bb'][0]:.2f} x {info['L_bb'][1]:.2f} x {info['L_bb'][2]:.2f}, "
        f"L_rg {info['L_rg'][0]:.2f} x {info['L_rg'][1]:.2f} x {info['L_rg'][2]:.2f}; bath rho = {info['rho_bulk']:.4f}",
        f"  Plates: parked {cfg['plate_park_gap']} sigma outside the BB faces, overhang {cfg['plate_overhang']}, spacing "
        f"{cfg['plate_spacing']} ({n_pl} beads in six plates); wet pistons {info['n_sheet']} beads each, "
        f"z = {info['z_perm_piston']:.2f} / {info['z_feed_piston']:.2f}, vacuum margins {cfg['margin_perm']} / {cfg['margin_feed']}",
        f'  Crosslinks: {c[1]}   Chain beads: {c[2]}   Solvent: {c[3]}   Plates (4,7-11): {n_pl}   '
        f'Feed piston: {c[5]}   Permeate piston: {c[6]}',
        f'  Total atoms: {info["n_atoms"]}   Total bonds: {info["n_bonds"]}   Output: {Path(info["output_file"]).name}',
        '',
    ])


def self_check(path, overlap_min=0.8):
    D = read_data(path)
    order = np.argsort(D['ids'])
    assert np.array_equal(D['ids'][order], np.arange(1, len(order) + 1)), 'atom ids are not 1..N'
    typ = D['typ'][order]
    n_mob = int((typ <= 3).sum())
    assert np.all(typ[:n_mob] <= 3), 'mobile beads must come first in a six-plate data file'
    return validate(typ, D['xyz'][order], D['bonds'][:, 1:4], D['box'], n_old=n_mob, overlap_min=overlap_min)


def check_gel_unchanged(input_file, output_file, drop_types=(4, 5)):
    """The network must come through as a RIGID copy: every bond vector identical to the
    input's minimum-image bond vector.  Returns (ok, message)."""
    a, b = read_data(input_file), read_data(output_file)
    La = np.array([a['box'][d + 'hi'] - a['box'][d + 'lo'] for d in 'xyz'])
    ia = {int(i): k for k, i in enumerate(a['ids'])}
    pa = a['xyz'][[ia[int(i)] for i in a['bonds'][:, 3]]] - a['xyz'][[ia[int(i)] for i in a['bonds'][:, 2]]]
    pa -= La * np.round(pa / La)
    ib = {int(i): k for k, i in enumerate(b['ids'])}
    pb = b['xyz'][[ib[int(i)] for i in b['bonds'][:, 3]]] - b['xyz'][[ib[int(i)] for i in b['bonds'][:, 2]]]
    n_poly = (int(np.isin(a['typ'], POLYMER_TYPES).sum()), int(np.isin(b['typ'], POLYMER_TYPES).sum()))
    if len(pa) != len(pb):
        return False, (f'bond count changed {len(pa)} -> {len(pb)}: a BONDED fragment was removed '
                       f'(polymer beads {n_poly[0]} -> {n_poly[1]})')
    dev = float(np.abs(pa - pb).max()) if len(pa) else 0.0
    ok = dev < 1e-5
    return ok, (f'{len(pa)} bond vectors identical to the input within {dev:.1e} sigma; polymer beads {n_poly[0]} -> {n_poly[1]}'
                + (f' ({n_poly[0] - n_poly[1]} unbonded free bead(s) removed)' if n_poly[0] != n_poly[1] else ''))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Convert an equilibrated isolated-gel snapshot into the six-plate + '
                                             'two-wet-piston data file for compress_slab.lmp.')
    ap.add_argument('--input')
    ap.add_argument('--output', default=None)
    ap.add_argument('--bath-margin-xy', type=float, default=12.0)
    ap.add_argument('--reservoir-z', type=float, default=15.0)
    ap.add_argument('--plate-park-gap', type=float, default=3.0)
    ap.add_argument('--plate-overhang', type=float, default=6.0)
    ap.add_argument('--plate-spacing', type=float, default=0.5)
    ap.add_argument('--piston-spacing', type=float, default=0.5)
    ap.add_argument('--piston-clearance', type=float, default=1.0)
    ap.add_argument('--margin-feed', type=float, default=10.0)
    ap.add_argument('--margin-perm', type=float, default=10.0)
    ap.add_argument('--keep-fragments', action='store_true',
                    help='refuse (instead of deleting) polymer not bonded to the main network')
    ap.add_argument('--no-log', action='store_true', help='do not append to slab_data_file_info.md')
    ap.add_argument('--no-png', action='store_true')
    ap.add_argument('--self-check-only', metavar='DATAFILE', help='validate an existing six-plate data file and exit')
    args = ap.parse_args()
    if args.self_check_only:
        sys.exit(0 if self_check(args.self_check_only)['ok'] else 1)
    if not args.input:
        ap.error('--input is required')
    cfg = Config(input_file=args.input, output_file=args.output, bath_margin_xy=args.bath_margin_xy,
                 reservoir_z=args.reservoir_z, plate_park_gap=args.plate_park_gap, plate_overhang=args.plate_overhang,
                 plate_spacing=args.plate_spacing, piston_spacing=args.piston_spacing,
                 piston_clearance=args.piston_clearance, margin_feed=args.margin_feed, margin_perm=args.margin_perm,
                 drop_fragments=not args.keep_fragments, log_info=not args.no_log, make_png=not args.no_png)
    convert(cfg)
    ok, msg = check_gel_unchanged(cfg.input_file, cfg.output_file)
    print(('  gel-unchanged check: OK -- ' if ok else '  gel-unchanged check: FAILED -- ') + msg)
    sys.exit(0 if ok else 1)
