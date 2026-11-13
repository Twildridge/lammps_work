"""Crosslinker fluctuations of a gel network from a LAMMPS trajectory (2026-11-12).

How far the crosslinkers wander about their mean positions, measured against the phantom-network
prediction for Gaussian strands -- the question being whether the network behaves as a lattice of
stretched (non-Gaussian) strands or as a thermally fluctuating entropic network.

Everything comes from the trajectory alone (`id type mol x y z`, wrapped): the bond topology is rebuilt
from the frames (a bonded pair never exceeds the FENE limit, in any frame), so no data file or bond dump
has to be local.

Per strand (crosslinker i -> crosslinker j, n bonds), with R(t) the end-to-end vector:
    Rbar   = <R>_t                      the mean (elastically active) vector
    dR2    = <|R - Rbar|^2>_t           the fluctuation
    Gamma  = <dR2> / <|Rbar|^2>         averaged over strands
Phantom network of Gaussian strands, functionality f, in its reference state:
    <dR2> = (2/f) N b^2,  <Rbar^2> = (1 - 2/f) N b^2   ->   Gamma_ph = 2/(f - 2)        (1 for f = 4)
Gamma << Gamma_ph: mean vectors stretched far beyond the fluctuation (lattice-like, non-Gaussian strands);
Gamma ~ Gamma_ph: Gaussian reference state.  No N b^2 enters the ratio.

Per crosslinker, with u(t) its displacement from its own time-mean position (crosslinker-COM drift
removed):
    delta  = sqrt(<u^2>) / <|Rbar|>     Lindemann ratio on the crosslinker sub-lattice
    phantom: <u^2> = (f-1)/(f(f-2)) N b^2   ->   delta_ph = sqrt(f-1)/(f-2)             (0.866 for f = 4)

Chain extension: x = sqrt(<R^2>)/(n b), b the measured mean bond length.  Gaussian strands give an
isotropic modulus on any cubic lattice (sum of R_i R_j is a second-rank tensor); the anisotropy enters
with the finite-extensibility correction, (3/5) x^2 to leading order for a freely jointed chain.
"""
from __future__ import annotations

import io
import itertools
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

XLINK_TYPE, BEAD_TYPE = 1, 2
R_BOND_MAX = 1.3          # a FENE bond (R0 = 1.5, mean 0.97) is inside this in EVERY frame; chance contacts are not


def read_polymer_frames(traj, types=(XLINK_TYPE, BEAD_TYPE), max_frames=None):
    """-> dict(ts, ids, typ, pos (n_frames, n, 3), box (n_frames, 3, 2)); atoms sorted by id, identical in
    every frame (frames whose polymer ids differ from the first are dropped)."""
    import pandas as pd
    ts, P, B, ids, typ = [], [], [], None, None
    with open(traj) as f:
        line = f.readline()
        while line:
            if line.startswith('ITEM: TIMESTEP'):
                t = int(f.readline()); f.readline(); n = int(f.readline()); f.readline()
                box = np.array([f.readline().split()[:2] for _ in range(3)], float)
                cols = f.readline().split()[2:]
                # the atom block goes through a buffer: read_csv on the open handle reads ahead past the
                # block and the next frame's header is lost
                df = pd.read_csv(io.StringIO(''.join(itertools.islice(f, n))), sep=r'\s+', header=None, names=cols, engine='c')
                df = df[df['type'].isin(types)].sort_values('id')
                i = df['id'].to_numpy()
                if ids is None:
                    ids, typ = i, df['type'].to_numpy()
                if np.array_equal(i, ids):
                    ts.append(t); B.append(box); P.append(df[['x', 'y', 'z']].to_numpy(float))
                if max_frames and len(ts) >= max_frames:
                    break
            line = f.readline()
    if not ts:
        return None
    return dict(ts=np.array(ts), ids=ids, typ=typ, pos=np.array(P), box=np.array(B))


def _min_image(d, L):
    return d - L * np.round(d / L)


def bonds_from_frames(F):
    """Bonded pairs (index arrays i, j) = pairs inside R_BOND_MAX in every frame."""
    L0 = F['box'][0, :, 1] - F['box'][0, :, 0]
    p0 = np.mod(F['pos'][0] - F['box'][0, :, 0], L0)
    p0[p0 >= L0] = 0.0
    pairs = cKDTree(p0, boxsize=L0).query_pairs(R_BOND_MAX, output_type='ndarray')
    i, j = pairs[:, 0], pairs[:, 1]
    keep = np.ones(len(i), bool)
    for k in range(len(F['ts'])):
        L = F['box'][k, :, 1] - F['box'][k, :, 0]
        d = _min_image(F['pos'][k][i] - F['pos'][k][j], L)
        keep &= (d ** 2).sum(axis=1) < R_BOND_MAX ** 2
    return i[keep], j[keep]


def strands_from_bonds(typ, bi, bj):
    """Walk from every crosslinker along chain beads to the next crosslinker.
    -> (a, b, nbonds) index arrays, each strand once (a < b); loops (a == b) and dangling ends dropped;
    the bond-degree array; and the atom-index path of every strand (a, beads..., b)."""
    n = len(typ)
    deg = np.bincount(np.r_[bi, bj], minlength=n)
    order = np.argsort(np.r_[bi, bj], kind='stable')
    nbr = np.r_[bj, bi][order]
    start = np.r_[0, np.cumsum(deg)]
    isx = typ == XLINK_TYPE
    A, Bq, N, paths = [], [], [], []
    for x in np.flatnonzero(isx):
        for first in nbr[start[x]:start[x + 1]]:
            prev, cur, nb, path = x, first, 1, [x, first]
            while not isx[cur]:
                nn = nbr[start[cur]:start[cur + 1]]
                nxt = nn[nn != prev]
                if len(nn) != 2 or len(nxt) != 1:      # dangling end, or a spurious extra contact
                    cur = -1
                    break
                prev, cur, nb = cur, nxt[0], nb + 1
                path.append(cur)
            if cur > x:
                A.append(x); Bq.append(cur); N.append(nb); paths.append(path)
    return np.array(A, int), np.array(Bq, int), np.array(N, int), deg, paths


def analyse(traj, f_target=None, max_frames=None):
    """Crosslinker-fluctuation statistics of one trajectory (see the module docstring) -> dict, or None
    when the file is absent / holds fewer than 3 usable frames."""
    traj = Path(traj)
    if not traj.exists():
        return None
    F = read_polymer_frames(traj, max_frames=max_frames)
    if F is None or len(F['ts']) < 3:
        return None
    nt = len(F['ts'])
    typ, pos = F['typ'], F['pos']
    Lbox = F['box'][:, :, 1] - F['box'][:, :, 0]
    bi, bj = bonds_from_frames(F)
    a, b, nb, deg, _ = strands_from_bonds(typ, bi, bj)
    isx = typ == XLINK_TYPE
    f = int(f_target or np.bincount(deg[isx]).argmax())
    # interior strands: both crosslinkers carry the full functionality
    full = (deg[a] == f) & (deg[b] == f)
    a, b, nb = a[full], b[full], nb[full]

    # strand vectors
    Rv = np.array([_min_image(pos[k][b] - pos[k][a], Lbox[k]) for k in range(nt)])     # (nt, ns, 3)
    Rbar = Rv.mean(axis=0)
    c = nt / (nt - 1.0)                                    # the mean is estimated from the same nt frames
    dR2_ax = ((Rv - Rbar) ** 2).mean(axis=0) * c           # (ns, 3)
    Rbar2 = (Rbar ** 2).sum(axis=1) - dR2_ax.sum(axis=1) / nt   # unbiased |<R>|^2
    R2 = (Rv ** 2).sum(axis=2).mean(axis=0)
    blen = np.mean([np.linalg.norm(_min_image(pos[k][bi] - pos[k][bj], Lbox[k]), axis=1).mean() for k in range(nt)])

    # crosslinker displacements about their own mean, crosslinker-COM drift removed
    xi = np.flatnonzero(isx & (deg == f))
    U = np.zeros((nt, len(xi), 3))
    for k in range(1, nt):
        U[k] = U[k - 1] + _min_image(pos[k][xi] - pos[k - 1][xi], Lbox[k])
    U -= U.mean(axis=1, keepdims=True)
    u2_ax = ((U - U.mean(axis=0)) ** 2).mean(axis=0) * c   # (nx, 3)

    def _boot(num, den, nboot=200, seed=0):
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, len(num), (nboot, len(num)))
        return np.percentile(num[idx].mean(axis=1) / den[idx].mean(axis=1), [2.5, 97.5])

    dR2 = dR2_ax.sum(axis=1)
    Gamma = dR2.mean() / Rbar2.mean()
    d_nn = np.sqrt(Rbar2.mean())
    x_ext = np.sqrt(R2.mean()) / (np.mean(nb) * blen)
    return dict(
        traj=traj.name, n_frames=nt, ts=F['ts'], f=f,
        n_xlink=int(isx.sum()), n_xlink_full=len(xi), n_strands=len(a),
        deg_hist=np.bincount(deg[isx]), nbonds=np.bincount(nb), nb_mean=float(np.mean(nb)), bond_len=float(blen),
        Gamma=float(Gamma), Gamma_ci=_boot(dR2, Rbar2), Gamma_ph=2.0 / (f - 2),
        Gamma_ax=dR2_ax.mean(axis=0) / (Rbar2.mean() / 3.0),          # per-axis fluctuation, same normalisation
        dR2=float(dR2.mean()), Rbar=float(d_nn), R_rms=float(np.sqrt(R2.mean())),
        lindemann=float(np.sqrt(u2_ax.sum(axis=1).mean()) / d_nn), lindemann_ph=np.sqrt(f - 1.0) / (f - 2),
        u2_ax=u2_ax.mean(axis=0),
        x_ext=float(x_ext), x_mean=float(d_nn / (np.mean(nb) * blen)), fjc_corr=float(0.6 * x_ext ** 2),
        strand_Gamma=dR2 / np.maximum(Rbar2, 1e-12), strand_x=np.sqrt(R2) / (nb * blen),
    )


def report(S, label=''):
    if S is None:
        print(f'crosslinker fluctuations{" " + label if label else ""}: no trajectory')
        return
    print(f"crosslinker fluctuations {label}: {S['n_frames']} frames (steps {S['ts'][0]}..{S['ts'][-1]}), "
          f"f = {S['f']}, {S['n_strands']} interior strands of {S['nb_mean']:.3g} bonds, b = {S['bond_len']:.3f}\n"
          f"  crosslinker degree histogram {S['deg_hist'].tolist()}\n"
          f"  strand:  |Rbar| = {S['Rbar']:.3f}   sqrt<dR^2> = {np.sqrt(S['dR2']):.3f}   "
          f"Gamma = <dR^2>/<Rbar^2> = {S['Gamma']:.3f}  [{S['Gamma_ci'][0]:.3f}, {S['Gamma_ci'][1]:.3f}]   "
          f"phantom Gaussian {S['Gamma_ph']:.3g}\n"
          f"  crosslinker Lindemann ratio sqrt<u^2>/|Rbar| = {S['lindemann']:.3f}   phantom Gaussian {S['lindemann_ph']:.3f}\n"
          f"  extension sqrt<R^2>/(n b) = {S['x_ext']:.3f}   leading finite-extensibility correction (3/5)x^2 = {S['fjc_corr']:.3f}")


# ---------------------------------------------------------------------------------------------------
#  Network modes (2026-11-12): crosslinker motion, chord stretch / swing, jump-rope, and the frozen
#  fraction of the polymer density f(q)
# ---------------------------------------------------------------------------------------------------
def _unwrap(F):
    """Unwrapped polymer positions (nt, n, 3), polymer-COM drift removed (frame 0 keeps its coordinates)."""
    pos, Lbox = F['pos'], F['box'][:, :, 1] - F['box'][:, :, 0]
    U = np.zeros_like(pos)
    for k in range(1, len(pos)):
        U[k] = U[k - 1] + _min_image(pos[k] - pos[k - 1], Lbox[k])
    U -= U.mean(axis=1, keepdims=True)
    return pos[0][None] + U


def frozen_fraction(X, q=None, groups=None):
    """Single-bead frozen fraction f_s(q) (the non-ergodicity / Lamb-Mossbauer factor of the self
    intermediate scattering function): f_s(q) = < |<exp(i q.u_j)>_t|^2 >_j, u_j the bead's displacement
    about its mean position, averaged over the three box axes; the finite-frame bias of the time mean is
    removed.  f_s -> 1: beads are localised on the scale 2 pi/q (a fixed frame at that scale); f_s -> 0:
    they explore more than that scale (the structure at that scale is thermally renewed).  A harmonic
    solid follows exp(-q^2 <u^2>_axis).  Works for ordered and random networks alike (no Bragg sampling).
    groups: {name: boolean mask over atoms} -> one curve each."""
    nt = X.shape[0]
    q = np.logspace(np.log10(0.25), np.log10(12.0), 28) if q is None else np.asarray(q, float)
    U = X - X.mean(axis=0)
    out = {k: np.zeros(len(q)) for k in (groups or {'all': None})}
    for iq, qq in enumerate(q):
        m2 = (np.abs(np.exp(1j * qq * U).mean(axis=0)) ** 2).mean(axis=1)       # (natoms,), mean over the 3 axes
        fs = m2 - (1.0 - m2) / (nt - 1.0)
        for k, msk in (groups or {'all': None}).items():
            out[k][iq] = fs.mean() if msk is None else fs[msk].mean()
    return dict(q=q, **out)


def network_modes(traj, max_frames=None, map_half=11.0, map_thick=3.0):
    """Decomposition of the network's thermal motion in the trajectory (strands of the modal length,
    both crosslinkers fully connected).  Per strand, about its mean chord Rbar (unit vector e):
        stretch   <(dR.e)^2>                    chord length fluctuation
        swing     <|dR|^2 - (dR.e)^2>/2         chord direction fluctuation, per transverse axis
        rope      <|d - dbar|^2>/2              jump-rope: transverse offset d of the middle bead from the
                                                instantaneous chord, fluctuating part, per axis
        xl        <u^2>/3 of its two crosslinkers
    Instantaneous (every strand, every frame): chord extension x = |R|/(n b) and the middle bead's
    offset h = |d|, for the joint histogram against the inextensible-rope bound h <= (n b/2) sqrt(1 - x^2).
    Also per-bead <u^2>, a slice of mean positions for the map, and frozen_fraction()."""
    traj = Path(traj)
    if not traj.exists():
        return None
    F = read_polymer_frames(traj, max_frames=max_frames)
    if F is None or len(F['ts']) < 3:
        return None
    nt, typ = len(F['ts']), F['typ']
    Lbox = F['box'][:, :, 1] - F['box'][:, :, 0]
    bi, bj = bonds_from_frames(F)
    a, b, nb, deg, paths = strands_from_bonds(typ, bi, bj)
    isx = typ == XLINK_TYPE
    f = int(np.bincount(deg[isx]).argmax())
    n = int(np.bincount(nb).argmax())
    keep = (deg[a] == f) & (deg[b] == f) & (nb == n)
    P = np.array([p for p, k in zip(paths, keep) if k])            # (ns, n + 1) atom indices
    a, b = P[:, 0], P[:, -1]
    X = _unwrap(F)                                                 # (nt, natoms, 3)
    c = nt / (nt - 1.0)
    blen = float(np.mean(np.linalg.norm(_min_image(X[:, bi] - X[:, bj], Lbox[:, None, :]), axis=2)))

    Rv = X[:, b] - X[:, a]                                         # unwrapped, so no image correction
    Rv = _min_image(Rv, Lbox[:, None, :])
    Rbar = Rv.mean(axis=0)
    e = Rbar / np.linalg.norm(Rbar, axis=1, keepdims=True)
    dR = Rv - Rbar
    par = (dR * e).sum(axis=2)
    stretch = (par ** 2).mean(axis=0) * c
    swing = ((dR ** 2).sum(axis=2) - par ** 2).mean(axis=0) * c / 2.0

    mid = P[:, n // 2]
    dm = _min_image(X[:, mid] - X[:, a], Lbox[:, None, :])
    Rn = np.linalg.norm(Rv, axis=2)
    d = dm - (dm * Rv).sum(axis=2, keepdims=True) * Rv / (Rn ** 2)[..., None]     # offset from the instantaneous chord
    h = np.linalg.norm(d, axis=2)
    rope = ((d - d.mean(axis=0)) ** 2).sum(axis=2).mean(axis=0) * c / 2.0
    bow = np.linalg.norm(d.mean(axis=0), axis=1)                    # static bend of the mean shape

    u2 = ((X - X.mean(axis=0)) ** 2).sum(axis=2).mean(axis=0) * c   # per atom, 3D
    xl = 0.5 * (u2[a] + u2[b]) / 3.0

    amp = dict(xl=xl, stretch=stretch, swing=swing, rope=rope)
    keys = list(amp)
    corr = np.corrcoef(np.array([amp[k] for k in keys]))

    # the map: mean positions in a slab, thin along y, centred on the gel
    M = X.mean(axis=0)
    cen = np.median(M[isx], axis=0)
    xi = float(np.linalg.norm(Rbar, axis=1).mean())
    if map_half is None:                       # the same number of meshes in view whatever the strand length
        map_half = 3.5 * xi
    sel = lambda idx: ((np.abs(M[idx, 1] - cen[1]) < map_thick / 2) & (np.abs(M[idx, 0] - cen[0]) < map_half)
                       & (np.abs(M[idx, 2] - cen[2]) < map_half))
    inpl = np.sqrt((((X - M) ** 2)[:, :, [0, 2]]).sum(axis=2).mean(axis=0) * c)     # rms in-plane (x, z) displacement
    sP = P[sel(P[:, 0]) | sel(P[:, -1])]
    atoms = np.flatnonzero(sel(np.arange(len(typ))))
    fq = frozen_fraction(X, groups=dict(xl=isx & (deg == f), bead=~isx))
    # instantaneous coupling inside a strand: does the middle bead come in when the chord lengthens?
    dl, dh = Rn - Rn.mean(axis=0), h - h.mean(axis=0)
    rope_stretch = float((dl * dh).sum() / np.sqrt((dl ** 2).sum() * (dh ** 2).sum()))
    return dict(
        traj=traj.name, n_frames=nt, ts=F['ts'], f=f, n=n, n_strands=len(P), bond_len=blen, Rbar=xi,
        amp={k: float(v.mean()) for k, v in amp.items()}, keys=keys, corr=corr, bow=float(np.sqrt((bow ** 2).mean())),
        stiffening=float(swing.mean() / stretch.mean()), rope_stretch=rope_stretch,
        x_inst=(Rn / (n * blen)).ravel(), h_inst=(h / (n * blen)).ravel(),
        u2_xl=float(u2[isx & (deg == f)].mean()), u2_bead=float(u2[~isx].mean()), u2_axis=u2 / 3.0,
        fq=fq, map=dict(M=M, r=inpl, isx=isx, strands=sP, atoms=atoms, cen=cen, half=map_half, thick=map_thick),
    )


def report_modes(S, label=''):
    if S is None:
        print(f'network modes{" " + label if label else ""}: no trajectory')
        return
    A = S['amp']
    print(f"network modes {label}: {S['n_frames']} frames, {S['n_strands']} strands of {S['n']} bonds, mean chord {S['Rbar']:.3f}\n"
          f"  mean-square amplitude per axis (sigma^2):  crosslinker {A['xl']:.3f}   chord stretch {A['stretch']:.3f}   "
          f"chord swing {A['swing']:.3f}   jump-rope (middle bead) {A['rope']:.3f}   static bow {S['bow']:.3f} sigma\n"
          f"  swing/stretch = {S['stiffening']:.2f}  (1 for Gaussian strands; = d ln f/d ln R of the strand, > 1 strain-stiffened)\n"
          f"  correlation across strands [{', '.join(S['keys'])}]:\n" + np.array2string(S['corr'], precision=2, prefix='  '))
    print(f"  instantaneous correlation of chord length and middle-bead offset within a strand: {S['rope_stretch']:+.2f}  "
          f"(-1 = an inextensible rope trading slack for length)")
    for k, nm in (('xl', 'crosslinkers'), ('bead', 'chain beads')):
        print(f"  {nm}: rms displacement {np.sqrt(S['u2_' + k]):.2f} sigma = {np.sqrt(S['u2_' + k]) / S['Rbar']:.2f} of the mean chord;  "
              f"frozen fraction at the mesh scale f_s(2 pi/chord) = {mesh_frozen(S, k):.2f}")


def mesh_frozen(S, k='bead'):
    """f_s at q = 2 pi / (mean chord): the share of beads still localised on the mesh scale."""
    return float(np.interp(2 * np.pi / S['Rbar'], S['fq']['q'], S['fq'][k]))
