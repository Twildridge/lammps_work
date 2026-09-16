#!/usr/bin/env python3
"""
slab_two_pistons.py
===================
CONVERTER: an equilibrated slab_with_support PERIODIC snapshot  ->  a two-piston
(feed / permeate) data file for triaxial_permeation_two_pist.lmp and
triaxial_compression_two_pist.lmp   (2026-09-16).

Why a converter and not a fresh lattice
---------------------------------------
The two-piston decks control pressure with the NPT-piston scheme of Marioni et
al., J. Membr. Sci. 738 (2026) 124837 (Eq. 3): each reservoir is closed by a
damped, force-loaded piston sheet that holds it at a prescribed pressure.  That
scheme acts in z ONLY (lx, ly are fixed for the whole run), so it can never swell
a fresh pre_swell=0.93 lattice laterally.  Starting from the slab_with_support
final configuration (aniso NPH at P*, piston transparent to solvent, sig_p,xx/zz
= sig_p,yy/zz = 1.0005) the gel arrives with zero transverse network stress and
its equilibrium lx, ly, gel dimensions by construction; the NPT-piston phase in
the new decks is then only a z-settle of the two reservoirs, not a swelling
equilibration.

Geometry produced (z from bottom to top; support-BELOW-gel convention kept so
every existing z-binning / analysis convention carries over -- it is the paper's
figure rotated by 90 degrees):

    zlo = 0
    | vacuum margin (permeate side, margin_perm)
    | permeate piston   (type 6, P_perm)                 <- NEW sheet
    | permeate reservoir  = existing below-support solvent, PADDED to permeate_thickness
    | support sheet     (type 4, frozen, WCA to polymer only -- unchanged)
    | gap | gel slab (types 1,2 + internal solvent) | feed reservoir (existing top solvent)
    |     [dry piston, type 7, parked at dry_piston_frac of the feed reservoir height]
    | feed piston       (type 5, P_feed)                 <- NEW sheet
    | vacuum margin (feed side, margin_feed)
    zhi

The old type-5 piston of the input is DELETED.  lx, ly are UNCHANGED (that is the
whole point of the converter).  Gel and solvent coordinates are identical to the
input up to one uniform z shift (and the z "roll", see unwrap_z, which is a no-op
for the standard input).

Atom types written (7): 1 crosslink, 2 chain bead, 3 solvent, 4 support,
5 feed piston, 6 permeate piston, 7 dry piston.  Masses are 1.0 for all seven;
the decks override 5/6/7 with `mass 5 ${piston_mass}` etc. after read_data so one
data file serves any piston mass.

Reuses isolate_gel.py's parser (parse_lammps_data) and bond validator
(validate_bonds); the writer is extended here because the two-piston file needs a
Velocities section, 7 atom types and a fully documented first line.

Usage (CLI)
-----------
  python3 slab_two_pistons.py --input <final_config_...data> [--output <...>]
        [--permeate-thickness 10] [--feed-thickness None] [--margin-perm auto]
        [--margin-feed 15] [--piston-clearance 1.0] [--dry-piston-frac 0.5]
        [--sheet-source support|hex] [--sheet-spacing 0.2] [--seed 42] [--no-log]
  The validation (self-check) always runs; --self-check-only validates an
  existing output file without rewriting it.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys
from dataclasses import dataclass, asdict, field
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import isolate_gel  # noqa: E402  (parser + FENE bond validator)

try:
    from scipy.spatial import cKDTree
    _SCIPY = True
except ImportError:  # pragma: no cover
    _SCIPY = False

# ---------------------------------------------------------------------------
# Constants (mirror the decks)
# ---------------------------------------------------------------------------
T_CROSSLINK, T_CHAIN, T_SOLVENT, T_SUPPORT = 1, 2, 3, 4
T_FEED, T_PERM, T_DRY = 5, 6, 7
POLYMER_TYPES = (T_CROSSLINK, T_CHAIN)
OVERLAP_MIN = 0.8            # sigma; same rejection radius as slab_with_support_periodic.ipynb
FENE_R0_MAX = isolate_gel.FENE_R0_MAX
TYPE_NAMES = {1: 'crosslink', 2: 'chain bead', 3: 'solvent', 4: 'support',
              5: 'feed piston', 6: 'permeate piston', 7: 'dry piston'}
DEFAULT_INPUT = ('../../lammps_data_files_local/'
                 'final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000002.data')
SMOKE_INPUT = ('../../lammps_data_files_local/'
               'final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000000.data')


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
@dataclass
class Config:
    """Every knob of the converter (the notebook's Config cell builds one)."""
    input_file: str = DEFAULT_INPUT
    output_file: str = None            # None -> <input stem>_two_pist.data next to the input
    permeate_thickness: float = 10.0   # sigma of solvent between support and permeate piston (padded up)
    feed_thickness: object = None      # None = keep the existing top-reservoir solvent; a number trims/pads to it
    margin_perm: object = None         # vacuum below the permeate piston; None -> feed_thickness + 5
    margin_feed: float = 15.0          # vacuum above the feed piston
    piston_clearance: float = 1.0      # wet piston plane this far outside the outermost solvent bead
    dry_piston_frac: float = 0.5       # dry piston z as a fraction of the feed reservoir height (from gel top)
    sheet_source: str = 'support'      # 'support': copy the input support's (x,y) pattern (identical count);
                                       # 'hex': fresh hex sheet at sheet_spacing (make_sheet logic)
    sheet_spacing: float = 0.2         # hex sheet spacing (only for sheet_source='hex'); snapped to tile lx, ly
    seed: int = 42
    density_window: float = 4.0        # bulk-density window starts this far above the gel top (sigma)
    log_info: bool = True              # append the entry to slab_data_file_info.md
    make_png: bool = True

    def __post_init__(self):
        if self.output_file is None:
            p = Path(self.input_file)
            self.output_file = str(p.with_name(p.stem + '_two_pist.data'))
        if self.feed_thickness is not None:
            self.feed_thickness = float(self.feed_thickness)
        if self.margin_perm is not None:
            self.margin_perm = float(self.margin_perm)
        assert self.sheet_source in ('support', 'hex'), "sheet_source must be 'support' or 'hex'"


# ---------------------------------------------------------------------------
# Parsing helpers (velocities are not read by isolate_gel's parser)
# ---------------------------------------------------------------------------
def parse_velocities(path):
    """Velocities section -> dict id -> (vx, vy, vz); empty if the file has none."""
    vel = {}
    with open(path) as f:
        lines = f.readlines()
    i = 0
    n = len(lines)
    while i < n:
        if lines[i].strip().startswith('Velocities'):
            i += 2
            while i < n and lines[i].strip() and not lines[i].strip()[0].isalpha():
                p = lines[i].split()
                if len(p) >= 4:
                    try:
                        vel[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
                    except ValueError:
                        pass
                i += 1
            break
        i += 1
    return vel


def parse_image_flags(path):
    """Atoms section -> dict id -> (ix, iy, iz) when the file carries image flags
    (write_data does); empty otherwise.  The network bonds across the x,y PBC, so
    LAMMPS needs the x,y flags of bonded atoms to stay consistent ("Inconsistent
    image flags" warning otherwise, and compute com/gyration/displace use them)."""
    flags = {}
    with open(path) as f:
        lines = f.readlines()
    i, n = 0, len(lines)
    while i < n:
        if lines[i].strip().startswith('Atoms'):
            i += 2
            while i < n and lines[i].strip() and not lines[i].strip()[0].isalpha():
                p = lines[i].split()
                if len(p) >= 9:
                    try:
                        flags[int(p[0])] = (int(p[6]), int(p[7]), int(p[8]))
                    except ValueError:
                        pass
                i += 1
            break
        i += 1
    return flags


def _arrays(atoms):
    ids = np.array([a['id'] for a in atoms], dtype=np.int64)
    typ = np.array([a['type'] for a in atoms], dtype=np.int64)
    xyz = np.array([[a['x'], a['y'], a['z']] for a in atoms], dtype=float)
    return ids, typ, xyz


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def gel_extent(zp):
    """Bounding-box, percentile and Rg-based z-extent of the polymer (decks' conventions)."""
    com = float(zp.mean())
    rg2 = float(np.mean((zp - com) ** 2))          # Rg_zz^2 (unit masses)
    L_rg = 2.0 * np.sqrt(3.0 * rg2)                # L = 2*sqrt(3)*sqrt(Rg_zz^2)  (uniform prism)
    return dict(bb_lo=float(zp.min()), bb_hi=float(zp.max()),
                p001=float(np.percentile(zp, 0.1)), p999=float(np.percentile(zp, 99.9)),
                com=com, L_rg=L_rg, rg_lo=com - 0.5 * L_rg, rg_hi=com + 0.5 * L_rg,
                L_bb=float(zp.max() - zp.min()))


def unwrap_z(xyz, typ, box, support_z, permeate_thickness):
    """Roll the periodic z coordinate so that the box reads
    below-support solvent | support | gel | feed solvent  with NO solvent crossing z.

    The periodic slab_with_support box has solvent on both sides of the z boundary
    (the top solvent and the below-support solvent are ONE reservoir through the
    PBC).  We cut that reservoir at a plane at most permeate_thickness below the
    support: if the input already has <= permeate_thickness of solvent below the
    support (the standard case) nothing moves; otherwise the excess below-support
    solvent is rolled to the top of the feed.  Returns (xyz, shift_applied)."""
    zlo, zhi = box['zlo'], box['zhi']
    Lz = zhi - zlo
    below = support_z - zlo
    if below <= permeate_thickness + 1e-9:
        print(f'  unwrap_z: {below:.2f} sigma of solvent below the support (<= permeate_thickness) -- no roll needed')
        return xyz, 0.0
    cut = support_z - permeate_thickness
    z = xyz[:, 2]
    znew = zlo + np.mod(z - cut, Lz)
    xyz = xyz.copy()
    xyz[:, 2] = znew
    print(f'  unwrap_z: {below:.2f} sigma below the support > permeate_thickness -> rolled the cut plane '
          f'z={cut:.2f} to zlo ({below - permeate_thickness:.2f} sigma moved to the top of the feed)')
    return xyz, cut - zlo


def make_sheet_hex(lx, ly, spacing):
    """Hex sheet tiling lx, ly exactly (same logic as make_sheet() in
    slab_with_support_periodic.ipynb: nx = round(lx/spacing), even row count)."""
    nx = max(1, round(lx / spacing))
    snap_dx = lx / nx
    row_pitch = spacing * np.sqrt(3.0)
    ny = max(2, round(ly / row_pitch))
    if ny % 2:
        ny += 1
    snap_dy = ly / ny
    pts = []
    for j in range(ny):
        y = j * snap_dy
        xshift = snap_dx / 2 if (j % 2) else 0.0
        for i in range(nx):
            pts.append(((i * snap_dx + xshift) % lx, y))
    return np.array(pts), (nx, ny, snap_dx, snap_dy)


def sheet_xy(cfg, box, support_xy):
    """(x, y) pattern for the three new sheets."""
    lx, ly = box['xhi'] - box['xlo'], box['yhi'] - box['ylo']
    if cfg.sheet_source == 'support':
        # The input support was dilated affinely with the box by the aniso barostat, so
        # its hex pattern still tiles lx, ly exactly and its bead count (23,316 for the
        # 9x9 box) is what the decks' N_pist assumptions refer to.  Copying it gives
        # pistons IDENTICAL in count and pattern to the existing sheets.
        xy = support_xy.copy()
        # in-row spacing: sort the beads of the lowest row (same y) by x
        y0 = np.round(xy[:, 1], 6)
        row = np.sort(xy[y0 == y0.min(), 0])
        spacing = float(np.median(np.diff(row))) if len(row) > 1 else np.nan
        print(f'  sheets: copying the support pattern ({len(xy)} beads/sheet; x-spacing now {spacing:.4f} sigma '
              f'after the box swelled; 0.2 sigma in the generator)')
        return xy, dict(source='support', n=len(xy), spacing=spacing)
    xy, (nx, ny, sdx, sdy) = make_sheet_hex(lx, ly, cfg.sheet_spacing)
    print(f'  sheets: fresh hex sheet {nx} x {ny} = {len(xy)} beads, snapped dx={sdx:.4f} dy(row)={sdy:.4f}')
    return xy, dict(source='hex', n=len(xy), spacing=sdx, nx=nx, ny=ny)


# ---------------------------------------------------------------------------
# Overlap-rejected random insertion (0.8 sigma, x,y periodic, z finite)
# ---------------------------------------------------------------------------
def _tree(xyz, box):
    lx, ly = box['xhi'] - box['xlo'], box['yhi'] - box['ylo']
    p = xyz.copy()
    p[:, 0] = np.mod(p[:, 0] - box['xlo'], lx)
    p[:, 1] = np.mod(p[:, 1] - box['ylo'], ly)
    p[:, 2] = p[:, 2] - box['zlo'] + 1.0e3         # finite z: shift positive, huge periodic length
    return cKDTree(p, boxsize=[lx, ly, 1.0e6]), (lx, ly)


def _wrap_like(xyz, box):
    lx, ly = box['xhi'] - box['xlo'], box['yhi'] - box['ylo']
    p = xyz.copy()
    p[:, 0] = np.mod(p[:, 0] - box['xlo'], lx)
    p[:, 1] = np.mod(p[:, 1] - box['ylo'], ly)
    p[:, 2] = p[:, 2] - box['zlo'] + 1.0e3
    return p


def insert_solvent(xyz_all, box, zlo, zhi, n_want, rng, label=''):
    """Place n_want solvent beads uniformly in [zlo, zhi) x full x,y, rejecting any
    candidate closer than OVERLAP_MIN to ANY existing atom (x,y minimum image) or
    to an already accepted candidate.  Returns the (m, 3) accepted positions."""
    if n_want <= 0:
        return np.zeros((0, 3))
    if not _SCIPY:
        raise RuntimeError('scipy is required for the overlap-rejected insertion')
    lx, ly = box['xhi'] - box['xlo'], box['yhi'] - box['ylo']
    accepted = []
    n_acc = 0
    attempts = 0
    max_attempts = 200 * n_want + 10000
    tree, _ = _tree(xyz_all, box)
    while n_acc < n_want and attempts < max_attempts:
        m = min(max(2000, 4 * (n_want - n_acc)), 50000)
        cand = np.column_stack([box['xlo'] + rng.random(m) * lx,
                                box['ylo'] + rng.random(m) * ly,
                                zlo + rng.random(m) * (zhi - zlo)])
        attempts += m
        q = _wrap_like(cand, box)
        # (1) no overlap with anything already present (existing atoms + accepted beads)
        nn = tree.query_ball_point(q, r=OVERLAP_MIN, return_length=True)
        ok = np.asarray(nn) == 0
        cand, q = cand[ok], q[ok]
        if not len(cand):
            continue
        # (2) no overlap within the batch (greedy: drop the second of every close pair)
        ct = cKDTree(q, boxsize=[lx, ly, 1.0e6])
        pairs = ct.query_pairs(OVERLAP_MIN, output_type='ndarray')
        drop = np.zeros(len(cand), bool)
        for a, b in pairs:
            if not drop[a] and not drop[b]:
                drop[b] = True
        cand = cand[~drop]
        take = cand[:n_want - n_acc]
        accepted.append(take)
        n_acc += len(take)
        # rebuild the tree with the new beads so the next batch respects them
        xyz_all = np.vstack([xyz_all, take])
        tree, _ = _tree(xyz_all, box)
    out = np.vstack(accepted) if accepted else np.zeros((0, 3))
    print(f'  insert_solvent{label}: placed {len(out)}/{n_want} beads in z=[{zlo:.2f}, {zhi:.2f}) '
          f'({attempts} candidates tried, {OVERLAP_MIN} sigma rejection)')
    if len(out) < n_want:
        print(f'  WARNING: only {len(out)} of {n_want} beads could be placed')
    return out



def fill_to_density(xyz, is_solv, box, z_lo, z_hi, rho_bulk, rng, label='', bin_w=0.5, empty_frac=0.5):
    """Bring the solvent in [z_lo, z_hi) up to rho_bulk by inserting ONLY into the
    z-intervals that are (nearly) empty.  A plain 0.8-sigma rejection over the whole
    range would also drop beads INSIDE the existing liquid (~40 % of random points in
    a rho~0.43 WCA fluid are > 0.8 sigma from every atom) and over-densify it while
    under-filling the pad (seen on the first converter run, 2026-09-16).  Existing
    solvent is therefore never touched: bins of width bin_w whose density is below
    empty_frac*rho_bulk are merged into intervals and each interval is filled to
    rho_bulk*A*height (its own existing beads subtracted).  Returns (m, 3)."""
    A = (box['xhi'] - box['xlo']) * (box['yhi'] - box['ylo'])
    edges = np.arange(z_lo, z_hi + bin_w * 0.5, bin_w)
    if edges[-1] < z_hi - 1e-9:
        edges = np.append(edges, z_hi)
    sz = xyz[is_solv, 2]
    cnt, _ = np.histogram(sz, bins=edges)
    dens = cnt / (A * np.diff(edges))
    empty = dens < empty_frac * rho_bulk
    intervals = []
    i = 0
    while i < len(empty):
        if empty[i]:
            j = i
            while j + 1 < len(empty) and empty[j + 1]:
                j += 1
            intervals.append((edges[i], edges[j + 1]))
            i = j + 1
        else:
            i += 1
    print(f'  fill_to_density{label}: {len(intervals)} empty interval(s) in z=[{z_lo:.2f}, {z_hi:.2f}): '
          + ', '.join(f'[{a:.2f},{b:.2f}]' for a, b in intervals))
    out = []
    cur = xyz
    for a, b in intervals:
        n_exist = int(((sz >= a) & (sz < b)).sum())
        n_want = int(round(rho_bulk * A * (b - a))) - n_exist
        if n_want <= 0:
            continue
        p = insert_solvent(cur, box, a, b, n_want, rng, label)
        if len(p):
            out.append(p)
            cur = np.vstack([cur, p])
    return np.vstack(out) if out else np.zeros((0, 3))

# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------
def write_two_pist_data(path, header_line, ids, mols, typ, xyz, vel, bonds, box, n_types=7, img=None):
    """LAMMPS data file: Atoms (molecular: id mol type x y z), Velocities, Bonds.
    `header_line` must be ONE line (the data-file format allows a single comment)."""
    header_line = ' '.join(str(header_line).split())
    with open(path, 'w') as f:
        f.write(header_line + '\n\n')
        f.write(f'{len(ids)} atoms\n{len(bonds)} bonds\n0 angles\n0 dihedrals\n0 impropers\n\n')
        f.write(f'{n_types} atom types\n1 bond types\n\n')
        f.write(f"{box['xlo']:.10f} {box['xhi']:.10f} xlo xhi\n")
        f.write(f"{box['ylo']:.10f} {box['yhi']:.10f} ylo yhi\n")
        f.write(f"{box['zlo']:.10f} {box['zhi']:.10f} zlo zhi\n\n")
        f.write('Masses\n\n')
        for t in range(1, n_types + 1):
            f.write(f'{t} 1.0  # {TYPE_NAMES.get(t, "")}' + ('  (deck overrides with piston_mass)' if t >= 5 else '') + '\n')
        f.write('\nAtoms # molecular\n\n')
        for i in range(len(ids)):
            f.write(f'{ids[i]} {mols[i]} {typ[i]} {xyz[i, 0]:.10f} {xyz[i, 1]:.10f} {xyz[i, 2]:.10f}'
                    + (f' {img[i, 0]} {img[i, 1]} 0\n' if img is not None else '\n'))
        f.write('\nVelocities\n\n')
        for i in range(len(ids)):
            f.write(f'{ids[i]} {vel[i, 0]:.10g} {vel[i, 1]:.10g} {vel[i, 2]:.10g}\n')
        if len(bonds):
            f.write('\nBonds\n\n')
            for k, (bt, a1, a2) in enumerate(bonds, 1):
                f.write(f'{k} {bt} {a1} {a2}\n')


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------
def convert(cfg: Config):
    """Run the whole conversion.  Returns a dict `info` (geometry, counts, paths)."""
    t_start = _dt.datetime.now()
    inp = Path(cfg.input_file)
    if not inp.exists():
        raise FileNotFoundError(f'input data file not found: {inp}')
    is_smoke = inp.name.endswith('_14000000.data')
    print('=' * 78)
    print(f'slab_two_pistons converter  ({t_start:%Y-%m-%d %H:%M})')
    print(f'  input : {inp}')
    print(f'  output: {cfg.output_file}')
    if is_smoke:
        print('  NOTE: this is the 14000000 snapshot (piston VISIBLE to solvent, Nose-Hoover T~1.24).\n'
              '        It is a SMOKE-TEST input only; the production input is the 2026-09-07\n'
              '        PISTON_TRANSPARENT=1 rerun ..._14000002.data (pull it from Expanse).')
    print('=' * 78)
    rng = np.random.default_rng(cfg.seed)

    # ---- (a) parse ------------------------------------------------------
    atoms, bonds_in, box, masses = isolate_gel.parse_lammps_data(str(inp))
    vel_in = parse_velocities(str(inp))
    img_in = parse_image_flags(str(inp))
    ids, typ, xyz = _arrays(atoms)
    mols = np.array([a['mol'] for a in atoms], dtype=np.int64)
    print(f'  read {len(ids)} atoms, {len(bonds_in)} bonds, {len(vel_in)} velocities; '
          f'box x[{box["xlo"]:.3f},{box["xhi"]:.3f}] y[{box["ylo"]:.3f},{box["yhi"]:.3f}] z[{box["zlo"]:.3f},{box["zhi"]:.3f}]')
    counts_in = {int(t): int((typ == t).sum()) for t in np.unique(typ)}
    print('  input type counts: ' + ', '.join(f'{t}={n}' for t, n in counts_in.items()))
    lx, ly = box['xhi'] - box['xlo'], box['yhi'] - box['ylo']
    A = lx * ly

    # ---- (b) identify groups, delete the old piston -----------------------
    is_pist_old = typ == T_FEED
    old_piston_z = float(xyz[is_pist_old, 2].mean()) if is_pist_old.any() else np.nan
    n_old_piston = int(is_pist_old.sum())
    keep = ~is_pist_old
    ids, typ, xyz, mols = ids[keep], typ[keep], xyz[keep], mols[keep]
    print(f'  deleted the old piston: {n_old_piston} type-5 beads at z={old_piston_z:.3f}')
    is_sup = typ == T_SUPPORT
    support_z = float(xyz[is_sup, 2].mean())
    sup_spread = float(xyz[is_sup, 2].std())
    assert sup_spread < 1e-3, f'support sheet is not flat (z std {sup_spread})'
    support_xy = xyz[is_sup, :2].copy()

    # ---- z unwrap (roll) so the box reads permeate | support | gel | feed ---
    xyz, roll = unwrap_z(xyz, typ, box, support_z, cfg.permeate_thickness)
    support_z = float(xyz[is_sup, 2].mean())
    is_poly = np.isin(typ, POLYMER_TYPES)
    is_solv = typ == T_SOLVENT
    G = gel_extent(xyz[is_poly, 2])
    below_layer = (is_solv) & (xyz[:, 2] < support_z)
    print(f'  support plane z = {support_z:.3f}; below-support solvent: {int(below_layer.sum())} beads in '
          f'[{xyz[below_layer, 2].min():.2f}, {support_z:.2f}] = {support_z - xyz[below_layer, 2].min():.2f} sigma')
    print(f'  gel z-extent: BB [{G["bb_lo"]:.2f}, {G["bb_hi"]:.2f}] (L_bb={G["L_bb"]:.2f}); '
          f'p0.1/p99.9 [{G["p001"]:.2f}, {G["p999"]:.2f}]; Rg-based [{G["rg_lo"]:.2f}, {G["rg_hi"]:.2f}] (L_rg={G["L_rg"]:.2f})')
    gel_top = G['bb_hi']                     # outermost polymer bead: safe for piston placement
    top_layer = is_solv & (xyz[:, 2] > gel_top)
    z_top_solv = float(xyz[is_solv, 2].max())
    print(f'  top (feed) solvent above the gel BB top: {int(top_layer.sum())} beads, up to z={z_top_solv:.2f} '
          f'({z_top_solv - gel_top:.2f} sigma above the gel top)')

    # ---- bulk solvent density from the feed reservoir, away from gel & old piston ----
    w_lo, w_hi = gel_top + cfg.density_window, box['zhi'] - 1.0
    if np.isfinite(old_piston_z) and (w_lo - 1.5 < old_piston_z < w_hi + 1.5):
        # old piston (visible to solvent in the 14000000 file) sits in the window: shrink it
        if old_piston_z > 0.5 * (w_lo + w_hi):
            w_hi = min(w_hi, old_piston_z - 1.5)
        else:
            w_lo = max(w_lo, old_piston_z + 1.5)
    n_win = int((is_solv & (xyz[:, 2] >= w_lo) & (xyz[:, 2] < w_hi)).sum())
    rho_bulk = n_win / (A * (w_hi - w_lo))
    print(f'  bulk solvent number density (feed window z=[{w_lo:.2f}, {w_hi:.2f}], {n_win} beads): rho = {rho_bulk:.4f} /sigma^3')

    # ---- (c) permeate padding ------------------------------------------
    z_perm_lo = support_z - cfg.permeate_thickness
    in_perm = is_solv & (xyz[:, 2] >= z_perm_lo) & (xyz[:, 2] < support_z)
    n_target = int(round(rho_bulk * A * cfg.permeate_thickness))
    print(f'  permeate reservoir target: {cfg.permeate_thickness} sigma x A={A:.2f} x rho={rho_bulk:.4f} = {n_target} beads; '
          f'existing {int(in_perm.sum())} in the region (empty sub-intervals are filled to rho_bulk)')
    # Pad ONLY below the existing layer.  Filling interior voids (e.g. the 1-2 sigma slot a
    # solvent-VISIBLE old piston leaves in the 14000000 smoke-test file) by random insertion
    # packs the beads toward the slot centre (both sides are liquid-bounded) and leaves a
    # local over-density; such a void is reported here instead and closes during the
    # NPT-piston settle (the pistons regulate the reservoir pressure).  The 14000002 input
    # (transparent piston) has no void.
    z_below_min = float(xyz[below_layer, 2].min())
    A_bin = 0.5
    e = np.arange(z_below_min, support_z + 1e-9, A_bin)
    if len(e) > 1:
        cnt, _ = np.histogram(xyz[below_layer, 2], bins=e)
        dens = cnt / (A * np.diff(e))
        n_void = int((dens < 0.5 * rho_bulk).sum())
        if n_void:
            print(f'  NOTE: {n_void} x {A_bin}-sigma bin(s) inside the existing below-support layer are nearly empty '
                  f'(old piston slot); left for the settle to close')
    pad_perm = fill_to_density(xyz, is_solv, box, z_perm_lo, z_below_min, rho_bulk, rng, ' [permeate pad]')
    n_perm_final = int(in_perm.sum()) + len(pad_perm)
    print(f'  permeate reservoir now holds {n_perm_final} beads = {n_perm_final / (A * cfg.permeate_thickness):.4f} /sigma^3 '
          f'(target {rho_bulk:.4f}; the NPT-piston settle absorbs the remainder)')
    xyz_solv_new = [pad_perm]

    # ---- feed trim / pad (only when feed_thickness is given) -----------------
    n_trim_feed = 0
    pad_feed = np.zeros((0, 3))
    if cfg.feed_thickness is not None:
        z_feed_hi = gel_top + cfg.feed_thickness
        above = is_solv & (xyz[:, 2] > z_feed_hi)
        if above.any():
            n_trim_feed = int(above.sum())
            keep = ~above
            ids, typ, xyz, mols = ids[keep], typ[keep], xyz[keep], mols[keep]
            is_poly, is_solv, is_sup = np.isin(typ, POLYMER_TYPES), typ == T_SOLVENT, typ == T_SUPPORT
            print(f'  feed trim: removed {n_trim_feed} solvent beads above z={z_feed_hi:.2f}')
        z_top_now = float(xyz[is_solv, 2].max())
        if z_feed_hi > z_top_now + 0.5:
            pad_feed = fill_to_density(np.vstack([xyz, pad_perm]), np.concatenate([is_solv, np.ones(len(pad_perm), bool)]),
                                       box, gel_top, z_feed_hi, rho_bulk, rng, ' [feed pad]')
            xyz_solv_new.append(pad_feed)

    # ---- (d) three sheets -------------------------------------------------
    solv_z_all = np.concatenate([xyz[is_solv, 2]] + [p[:, 2] for p in xyz_solv_new if len(p)])
    z_min_solv, z_max_solv = float(solv_z_all.min()), float(solv_z_all.max())
    z_perm_piston = z_min_solv - cfg.piston_clearance
    z_feed_piston = z_max_solv + cfg.piston_clearance
    z_dry = gel_top + cfg.dry_piston_frac * (z_feed_piston - gel_top)
    # the dry piston must not overlap POLYMER (solvent passes through it): nudge up if needed
    zp = xyz[is_poly, 2]
    nudge = 0
    while np.any(np.abs(zp - z_dry) < OVERLAP_MIN) and nudge < 40:
        z_dry += 0.25
        nudge += 1
    if nudge:
        print(f'  dry piston plane nudged up {0.25 * nudge:.2f} sigma to clear stray polymer beads')
    sheet, sheet_info = sheet_xy(cfg, box, support_xy)
    n_sheet = len(sheet)
    feed_thick_actual = z_feed_piston - gel_top
    margin_perm = cfg.margin_perm if cfg.margin_perm is not None else feed_thick_actual + 5.0
    print(f'  pistons: permeate z={z_perm_piston:.3f}  feed z={z_feed_piston:.3f}  dry z={z_dry:.3f}   '
          f'(feed reservoir {feed_thick_actual:.2f} sigma from gel BB top to the feed piston)')
    print(f'  margins: permeate side {margin_perm:.2f} sigma (auto = feed + 5 unless given), feed side {cfg.margin_feed:.2f} sigma')

    # ---- (e) assemble, shift so zlo = 0, wrap x,y, renumber ------------------
    new_pos, new_typ = [], []
    for p in xyz_solv_new:
        if len(p):
            new_pos.append(p)
            new_typ.append(np.full(len(p), T_SOLVENT))
    for zpl, t in ((z_feed_piston, T_FEED), (z_perm_piston, T_PERM), (z_dry, T_DRY)):
        new_pos.append(np.column_stack([sheet[:, 0] + box['xlo'] * 0.0, sheet[:, 1], np.full(n_sheet, zpl)]))
        new_typ.append(np.full(n_sheet, t))
    new_pos = np.vstack(new_pos)
    new_typ = np.concatenate(new_typ).astype(np.int64)
    # the copied support pattern is already in box coordinates; a hex sheet starts at 0 -> add xlo/ylo
    if cfg.sheet_source == 'hex':
        m = new_typ >= T_FEED
        new_pos[m, 0] += box['xlo']
        new_pos[m, 1] += box['ylo']

    n_old = len(ids)
    all_xyz = np.vstack([xyz, new_pos])
    all_typ = np.concatenate([typ, new_typ])
    mol_next = int(mols.max()) + 1
    all_mol = np.concatenate([mols, np.arange(mol_next, mol_next + len(new_typ))])
    # velocities: keep types 1-3 (and support, which is 0); new atoms 0
    vel = np.zeros((len(all_typ), 3))
    for i in range(n_old):
        v = vel_in.get(int(ids[i]))
        if v is not None and all_typ[i] in (1, 2, 3):
            vel[i] = v
    # image flags: x,y carried over for the original atoms (bonded network crosses the
    # x,y PBC), z reset to 0 everywhere (the new box is finite in z and nothing crosses it)
    img = np.zeros((len(all_typ), 3), dtype=np.int64)
    for i in range(n_old):
        fl = img_in.get(int(ids[i]))
        if fl is not None:
            img[i, 0], img[i, 1] = fl[0], fl[1]
    # id remap (contiguous), bonds
    old2new = {int(o): i + 1 for i, o in enumerate(ids)}
    new_ids = np.arange(1, len(all_typ) + 1)
    bonds = []
    n_orphan = 0
    for b in bonds_in:
        a1, a2 = old2new.get(b['atom1']), old2new.get(b['atom2'])
        if a1 is None or a2 is None:
            n_orphan += 1
            continue
        bonds.append((b['type'], a1, a2))
    assert n_orphan == 0, f'{n_orphan} bonds lost an endpoint -- polymer must never be trimmed'
    # z shift so zlo = 0; x,y wrapped into the (unchanged) lateral box
    zlo_new = z_perm_piston - margin_perm
    zhi_new = z_feed_piston + cfg.margin_feed
    all_xyz[:, 2] -= zlo_new
    all_xyz[:, 0] = box['xlo'] + np.mod(all_xyz[:, 0] - box['xlo'], lx)
    all_xyz[:, 1] = box['ylo'] + np.mod(all_xyz[:, 1] - box['ylo'], ly)
    new_box = dict(xlo=box['xlo'], xhi=box['xhi'], ylo=box['ylo'], yhi=box['yhi'],
                   zlo=0.0, zhi=zhi_new - zlo_new)
    zs = lambda z: z - zlo_new    # noqa: E731  (input z -> output z)

    # ---- geometry record --------------------------------------------------
    layout = [
        ('zhi (box top)', new_box['zhi']),
        ('feed piston (type 5)', zs(z_feed_piston)),
        ('top of feed solvent', zs(z_max_solv)),
        ('dry piston (type 7)', zs(z_dry)),
        ('gel top (BB / Rg)', f'{zs(G["bb_hi"]):.3f} / {zs(G["rg_hi"]):.3f}'),
        ('gel bottom (BB / Rg)', f'{zs(G["bb_lo"]):.3f} / {zs(G["rg_lo"]):.3f}'),
        ('support (type 4)', zs(support_z)),
        ('bottom of permeate solvent', zs(z_min_solv)),
        ('permeate piston (type 6)', zs(z_perm_piston)),
        ('zlo (box bottom)', 0.0),
    ]
    info = dict(
        input_file=str(inp), output_file=str(cfg.output_file), date=t_start.strftime('%Y-%m-%d %H:%M'),
        smoke_test_input=is_smoke, config=asdict(cfg), lx=lx, ly=ly, area=A,
        box_in=dict(box), box_out=new_box, z_shift=-zlo_new, z_roll=roll,
        old_piston_z_in=old_piston_z, n_old_piston=n_old_piston,
        support_z=zs(support_z), gel_bb=(zs(G['bb_lo']), zs(G['bb_hi'])), gel_rg=(zs(G['rg_lo']), zs(G['rg_hi'])),
        L_bb=G['L_bb'], L_rg=G['L_rg'],
        z_feed_piston=zs(z_feed_piston), z_perm_piston=zs(z_perm_piston), z_dry_piston=zs(z_dry),
        feed_thickness_actual=feed_thick_actual, permeate_thickness=cfg.permeate_thickness,
        margin_perm=margin_perm, margin_feed=cfg.margin_feed,
        rho_bulk=rho_bulk, density_window=(zs(w_lo), zs(w_hi)),
        n_added_permeate=int(len(pad_perm)), n_added_feed=int(len(pad_feed)), n_trimmed_feed=n_trim_feed,
        sheet=sheet_info, n_sheet=n_sheet, n_atoms=int(len(all_typ)), n_bonds=len(bonds),
        counts={int(t): int((all_typ == t).sum()) for t in range(1, 8)},
        layout=layout,
        # the feed piston rises by (expelled solvent)/A under compression: this is the largest
        # gel strain the feed-side margin can absorb before the piston reaches the box top
        max_strain_feed_margin=(cfg.margin_feed - 2.0) / G['L_bb'],
    )

    # ---- (f) validate BEFORE writing ------------------------------------------
    val = validate(all_typ, all_xyz, bonds, new_box, n_old, info)
    info['validation'] = val
    if not val['ok']:
        raise RuntimeError('validation FAILED -- file not written:\n  ' + '\n  '.join(val['problems']))

    # ---- write -----------------------------------------------------------------
    cfgs = ' '.join(f'{k}={v}' for k, v in asdict(cfg).items())
    header = (f'two-piston (feed/permeate) data file from {inp.name} via slab_two_pistons.py on {info["date"]}; '
              f'types 1 xl 2 chain 3 solv 4 support 5 feed_piston 6 perm_piston 7 dry_piston; '
              f'{cfgs}; z: support={info["support_z"]:.3f} gel_bb=[{info["gel_bb"][0]:.3f},{info["gel_bb"][1]:.3f}] '
              f'gel_rg=[{info["gel_rg"][0]:.3f},{info["gel_rg"][1]:.3f}] perm_piston={info["z_perm_piston"]:.3f} '
              f'dry_piston={info["z_dry_piston"]:.3f} feed_piston={info["z_feed_piston"]:.3f}; '
              f'rho_bulk={rho_bulk:.4f} n_added_perm={len(pad_perm)}' + (' SMOKE-TEST-INPUT(14000000)' if is_smoke else ''))
    write_two_pist_data(cfg.output_file, header, new_ids, all_mol, all_typ, all_xyz, vel, bonds, new_box,
                        img=(img if img_in else None))
    print(f'  image flags: x,y carried over from the input for {n_old} atoms ({len(img_in)} read), z reset to 0'
          if img_in else '  image flags: none in the input -> none written')
    with open(str(cfg.output_file) + '.info.json', 'w') as f:
        json.dump(_jsonable(info), f, indent=1)
    print(f'  wrote {len(all_typ)} atoms, {len(bonds)} bonds -> {cfg.output_file}')
    print(f'  sidecar: {cfg.output_file}.info.json')

    print_layout(info)
    if cfg.make_png:
        info['png'] = plot_z_histogram(all_typ, all_xyz, new_box, info, str(cfg.output_file).replace('.data', '_zprofile.png'))
    entry = info_log_entry(info)
    print('\n--- slab_data_file_info.md entry ---\n' + entry)
    if cfg.log_info:
        append_info_log(entry)
    dt = (_dt.datetime.now() - t_start).total_seconds()
    print(f'\nconversion done in {dt:.1f} s')
    return info


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return o


# ---------------------------------------------------------------------------
# Validation (also the CLI self-check)
# ---------------------------------------------------------------------------
def validate(typ, xyz, bonds, box, n_old, info=None, verbose=True):
    """Checks: (1) no pair closer than 0.8 sigma between NEW beads and anything they
    interact with (solvent-dry piston and wall-wall pairs are OFF in the decks and are
    exempt); (2) every FENE bond < 1.5 under x,y minimum image (isolate_gel.validate_bonds);
    (3) atom-type counts; (4) no solvent outside [permeate piston, feed piston];
    (5) no atom outside the box.  Returns dict(ok, problems, ...)."""
    problems = []
    say = print if verbose else (lambda *a, **k: None)
    say('\n--- validation ---')
    lx, ly = box['xhi'] - box['xlo'], box['yhi'] - box['ylo']
    # (1) overlaps
    n_close = 0
    worst = np.inf
    if _SCIPY:
        p = _wrap_like(xyz, box)
        tree = cKDTree(p, boxsize=[lx, ly, 1.0e6])
        pairs = tree.query_pairs(OVERLAP_MIN, output_type='ndarray')
        if len(pairs):
            a, b = pairs[:, 0], pairs[:, 1]
            ta, tb = typ[a], typ[b]
            is_new = (a >= n_old) | (b >= n_old)
            wall = np.isin(ta, (4, 5, 6, 7)) & np.isin(tb, (4, 5, 6, 7))          # wall-wall: off
            dry_solv = ((ta == T_DRY) & (tb == T_SOLVENT)) | ((tb == T_DRY) & (ta == T_SOLVENT))  # dry-solvent: off
            sup_solv = ((ta == T_SUPPORT) & (tb == T_SOLVENT)) | ((tb == T_SUPPORT) & (ta == T_SOLVENT))  # off
            rel = is_new & ~wall & ~dry_solv & ~sup_solv
            n_close = int(rel.sum())
            if n_close:
                d = np.linalg.norm(p[a[rel]] - p[b[rel]], axis=1)
                worst = float(d.min())
                bad = pairs[rel][:5]
                problems.append(f'{n_close} interacting pair(s) closer than {OVERLAP_MIN} sigma involving new beads '
                                f'(closest {worst:.3f}); e.g. types ' + ', '.join(f'{typ[i]}-{typ[j]}' for i, j in bad))
            n_pre = int((~is_new & ~wall & ~dry_solv & ~sup_solv).sum())
            say(f'  overlaps < {OVERLAP_MIN}: {n_close} involving new beads (interacting pairs)'
                + (f'; {n_pre} pre-existing close pairs in the input (untouched)' if n_pre else ''))
    else:
        say('  (scipy missing -> overlap check skipped)')
    # (2) bonds
    atoms_l = [dict(id=i + 1, x=xyz[i, 0], y=xyz[i, 1], z=xyz[i, 2]) for i in range(len(typ)) if typ[i] in POLYMER_TYPES]
    bonds_l = [dict(atom1=a1, atom2=a2) for (_, a1, a2) in bonds]
    try:
        isolate_gel.validate_bonds(atoms_l, bonds_l, box, FENE_R0_MAX)
    except ValueError as e:
        problems.append(f'bond check failed: {e}')
    # no bond may cross the (now finite) z direction
    zc = 0
    Lz = box['zhi'] - box['zlo']
    pos = {a['id']: a for a in atoms_l}
    for b in bonds_l:
        if abs(pos[b['atom2']]['z'] - pos[b['atom1']]['z']) > 0.5 * Lz:
            zc += 1
    if zc:
        problems.append(f'{zc} bond(s) cross the z boundary')
    say(f'  bonds crossing z: {zc}')
    # (3) counts
    counts = {int(t): int((typ == t).sum()) for t in range(1, 8)}
    say('  type counts: ' + ', '.join(f'{t} {TYPE_NAMES[t]}={n}' for t, n in counts.items()))
    for t in (5, 6, 7):
        if counts[t] == 0:
            problems.append(f'no atoms of type {t} ({TYPE_NAMES[t]})')
    if not (counts[5] == counts[6] == counts[7] == counts[4]):
        say(f'  NOTE: sheet counts differ (support {counts[4]}, pistons {counts[5]}/{counts[6]}/{counts[7]})')
    # (4) solvent between the wet pistons
    z_feed = float(xyz[typ == T_FEED, 2].mean()) if counts[5] else np.inf
    z_perm = float(xyz[typ == T_PERM, 2].mean()) if counts[6] else -np.inf
    sz = xyz[typ == T_SOLVENT, 2]
    n_out = int(((sz < z_perm) | (sz > z_feed)).sum())
    say(f'  solvent z in [{sz.min():.3f}, {sz.max():.3f}] vs pistons [{z_perm:.3f}, {z_feed:.3f}]: {n_out} outside')
    if n_out:
        problems.append(f'{n_out} solvent beads outside [permeate piston, feed piston]')
    # polymer must also lie between the wet pistons and above the support
    pz = xyz[np.isin(typ, POLYMER_TYPES), 2]
    n_pout = int(((pz < z_perm) | (pz > z_feed)).sum())
    if n_pout:
        problems.append(f'{n_pout} polymer beads outside [permeate piston, feed piston]')
    # (5) inside the box
    out = ((xyz[:, 0] < box['xlo']) | (xyz[:, 0] >= box['xhi']) | (xyz[:, 1] < box['ylo']) | (xyz[:, 1] >= box['yhi'])
           | (xyz[:, 2] < box['zlo']) | (xyz[:, 2] > box['zhi']))
    n_box = int(out.sum())
    say(f'  atoms outside the box: {n_box}')
    if n_box:
        problems.append(f'{n_box} atoms outside the box')
    # gel coordinates identical to the input up to a uniform z shift: checked by the caller
    ok = not problems
    say('  RESULT: ' + ('OK' if ok else 'PROBLEMS:\n    ' + '\n    '.join(problems)))
    return dict(ok=ok, problems=problems, n_close=n_close, closest=worst, counts=counts,
                n_solvent_outside=n_out, n_outside_box=n_box, n_bonds_cross_z=zc)


def print_layout(info):
    print('\n--- z layout (output coordinates, zlo = 0) ---')
    print(f'  {"plane":<30s} {"z (sigma)":>12s}')
    for name, z in info['layout']:
        zs = f'{z:.3f}' if isinstance(z, (int, float, np.floating)) else str(z)
        print(f'  {name:<30s} {zs:>12s}')
    print(f'  lx = {info["lx"]:.4f}  ly = {info["ly"]:.4f}  (unchanged)   lz = {info["box_out"]["zhi"]:.3f}')
    print(f'  gel L_bb = {info["L_bb"]:.3f}  L_rg = {info["L_rg"]:.3f};  feed reservoir {info["feed_thickness_actual"]:.2f} sigma; '
          f'permeate reservoir {info["permeate_thickness"]:.2f} sigma')
    print(f'  feed-side margin {info["margin_feed"]:.1f} sigma absorbs the solvent expelled up to a compression strain of '
          f'~{info["max_strain_feed_margin"]:.2f} (piston must stay 2 sigma below zhi); raise margin_feed for deeper sweeps')


def plot_z_histogram(typ, xyz, box, info, out_png):
    """z number-density histogram per atom type (PNG next to the data file)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    A = info['area']
    edges = np.arange(box['zlo'], box['zhi'] + 0.5, 0.5)
    fig, ax = plt.subplots(figsize=(11, 5))
    cols = {1: 'tab:orange', 2: 'tab:red', 3: 'tab:blue', 4: 'k', 5: 'tab:green', 6: 'tab:purple', 7: 'tab:brown'}
    for t in range(1, 8):
        z = xyz[typ == t, 2]
        if not len(z):
            continue
        h, _ = np.histogram(z, bins=edges)
        dens = h / (A * 0.5)
        if t in (4, 5, 6, 7):
            ax.axvline(float(z.mean()), color=cols[t], lw=2, ls='--', label=f'{t} {TYPE_NAMES[t]} (z={z.mean():.1f})')
        else:
            ax.plot(0.5 * (edges[1:] + edges[:-1]), dens, color=cols[t], lw=1.3, label=f'{t} {TYPE_NAMES[t]}')
    ax.set_xlabel('z (sigma)')
    ax.set_ylabel('number density (1/sigma^3)')
    ax.set_title(Path(info['output_file']).name, fontsize=9)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    print(f'  z-density histogram -> {out_png}')
    return out_png


def info_log_entry(info):
    c = info['counts']
    cfg = info['config']
    lines = [
        f'Converted two-piston (feed/permeate) data file  [{info["date"]}, slab_two_pistons.py]:',
        f'  Input: {Path(info["input_file"]).name}' + ('   (SMOKE-TEST input: 14000000 snapshot, not the production 14000002)' if info['smoke_test_input'] else ''),
        f'  Box: {info["lx"]:.4f} x {info["ly"]:.4f} x {info["box_out"]["zhi"]:.3f}   (lx, ly unchanged from the input; zlo = 0)',
        f'  z: perm piston {info["z_perm_piston"]:.2f} | support {info["support_z"]:.2f} | gel BB [{info["gel_bb"][0]:.2f}, {info["gel_bb"][1]:.2f}] '
        f'(L_bb {info["L_bb"]:.2f}, L_rg {info["L_rg"]:.2f}) | dry piston {info["z_dry_piston"]:.2f} | feed piston {info["z_feed_piston"]:.2f}',
        f'  Reservoirs: permeate {info["permeate_thickness"]:.1f} sigma (padded +{info["n_added_permeate"]} beads at rho={info["rho_bulk"]:.4f}), '
        f'feed {info["feed_thickness_actual"]:.2f} sigma; margins perm {info["margin_perm"]:.1f} / feed {info["margin_feed"]:.1f} sigma (vacuum)',
        f'  Sheets: {info["n_sheet"]} beads each ({info["sheet"]["source"]} pattern, spacing {info["sheet"]["spacing"]:.4f}); '
        f'piston_clearance {cfg["piston_clearance"]}, dry_piston_frac {cfg["dry_piston_frac"]}, seed {cfg["seed"]}',
        f'  Crosslinks: {c[1]}   Chain beads: {c[2]}   Solvent: {c[3]}   Support: {c[4]}   Feed piston: {c[5]}   Permeate piston: {c[6]}   Dry piston: {c[7]}',
        f'  Total atoms: {info["n_atoms"]}   Total bonds: {info["n_bonds"]}   Output: {Path(info["output_file"]).name}',
        '',
    ]
    return '\n'.join(lines)


def append_info_log(entry, path=None):
    path = Path(path) if path else _HERE.parent / 'slab_data_file_info.md'
    with open(path, 'a') as f:
        f.write('\n' + entry)
    print(f'  appended to {path}')


# ---------------------------------------------------------------------------
# Standalone self-check of an existing output file
# ---------------------------------------------------------------------------
def self_check(path):
    atoms, bonds_in, box, _ = isolate_gel.parse_lammps_data(path)
    ids, typ, xyz = _arrays(atoms)
    old2new = {int(o): i + 1 for i, o in enumerate(ids)}
    bonds = [(b['type'], old2new[b['atom1']], old2new[b['atom2']]) for b in bonds_in]
    # treat every atom as "new" so every interacting pair is checked
    return validate(typ, xyz, bonds, box, n_old=0)


def check_gel_unchanged(input_file, output_file, tol=1e-6):
    """Confirm lx, ly are unchanged and gel (types 1,2) coordinates are identical to
    the input up to a uniform z shift (and any z roll).  Returns (ok, message)."""
    a_in, _, b_in, _ = isolate_gel.parse_lammps_data(input_file)
    a_out, _, b_out, _ = isolate_gel.parse_lammps_data(output_file)
    pin = np.array([[a['x'], a['y'], a['z']] for a in a_in if a['type'] in POLYMER_TYPES])
    pout = np.array([[a['x'], a['y'], a['z']] for a in a_out if a['type'] in POLYMER_TYPES])
    msgs = []
    ok = True
    for d in ('x', 'y'):
        if abs((b_in[d + 'hi'] - b_in[d + 'lo']) - (b_out[d + 'hi'] - b_out[d + 'lo'])) > tol:
            ok = False
            msgs.append(f'l{d} changed')
    if len(pin) != len(pout):
        return False, f'polymer count changed {len(pin)} -> {len(pout)}'
    dz = pout[:, 2] - pin[:, 2]
    shift = float(np.median(dz))
    resid = np.abs(dz - shift)
    dxy = np.abs(pout[:, :2] - pin[:, :2])
    lx, ly = b_in['xhi'] - b_in['xlo'], b_in['yhi'] - b_in['ylo']
    dxy[:, 0] = np.minimum(dxy[:, 0], lx - dxy[:, 0])
    dxy[:, 1] = np.minimum(dxy[:, 1], ly - dxy[:, 1])
    if resid.max() > 1e-5 or dxy.max() > 1e-5:
        ok = False
        msgs.append(f'gel coordinates differ beyond a uniform z shift (max |dz-shift| {resid.max():.2e}, max dxy {dxy.max():.2e})')
    msgs.append(f'lx={lx:.4f} ly={ly:.4f} unchanged; gel z shift = {shift:+.4f} (max deviation {resid.max():.1e}, max xy deviation {dxy.max():.1e})')
    return ok, '; '.join(msgs)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _none_or_float(s):
    return None if str(s).lower() in ('none', 'auto', '') else float(s)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Convert an equilibrated slab_with_support periodic snapshot into a '
                                             'two-piston (feed/permeate) data file.')
    ap.add_argument('--input', default=DEFAULT_INPUT)
    ap.add_argument('--output', default=None)
    ap.add_argument('--permeate-thickness', type=float, default=10.0)
    ap.add_argument('--feed-thickness', type=_none_or_float, default=None)
    ap.add_argument('--margin-perm', type=_none_or_float, default=None)
    ap.add_argument('--margin-feed', type=float, default=15.0)
    ap.add_argument('--piston-clearance', type=float, default=1.0)
    ap.add_argument('--dry-piston-frac', type=float, default=0.5)
    ap.add_argument('--sheet-source', choices=('support', 'hex'), default='support')
    ap.add_argument('--sheet-spacing', type=float, default=0.2)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--no-log', action='store_true', help='do not append to slab_data_file_info.md')
    ap.add_argument('--no-png', action='store_true')
    ap.add_argument('--self-check-only', metavar='DATAFILE', help='validate an existing two-piston data file and exit')
    args = ap.parse_args()
    if args.self_check_only:
        r = self_check(args.self_check_only)
        sys.exit(0 if r['ok'] else 1)
    cfg = Config(input_file=args.input, output_file=args.output, permeate_thickness=args.permeate_thickness,
                 feed_thickness=args.feed_thickness, margin_perm=args.margin_perm, margin_feed=args.margin_feed,
                 piston_clearance=args.piston_clearance, dry_piston_frac=args.dry_piston_frac,
                 sheet_source=args.sheet_source, sheet_spacing=args.sheet_spacing, seed=args.seed,
                 log_info=not args.no_log, make_png=not args.no_png)
    info = convert(cfg)
    ok, msg = check_gel_unchanged(cfg.input_file, cfg.output_file)
    print(('  gel-unchanged check: OK -- ' if ok else '  gel-unchanged check: FAILED -- ') + msg)
    sys.exit(0 if ok else 1)
