"""psd.py -- geometric porosity and pore-size distribution of the polymer network
(2026-09-24), in-house, no external binary.

Used by lib/triaxial.py (add_perm_psd / fig_perm_psd) on the SAME trajectory frames
the Voronoi pass streams (the polymer coordinates are kept by _volume_fractions, so
the 240 MB dump is never re-read for this).

Method
------
Polymer beads only (types 1, 2; solvent 3 and the walls 4-7 are dropped), wrapped into
the box.  A regular grid of spacing ~h covers the membrane region (full x, y; a z
window handed in by the caller), periodic in x and y only.  For every grid point the
distance to the nearest bead SURFACE

    d = r_nn - sigma_bead/2                       (cKDTree, boxsize = (Lx, Ly, 2 Lz))

is the radius of the largest empty sphere CENTRED there.

* geometric porosity per z-bin  eps_g(r_probe) = fraction of grid points with
  d >= r_probe (default r_probe = 0.5 = the solvent bead radius): the volume a
  solvent-sized probe can reach.  It is a geometric quantity, not the thermodynamic
  solvent fraction phi_s (which the lambda-calibrated Voronoi gives).
* pore size per void point by the Gelb-Gubbins covering step: a point's pore
  diameter is that of the LARGEST sphere that contains it (not the sphere centred on
  it).  Implemented on the grid by descending radius levels r_k: the points within r_k
  of any centre with d >= r_k are covered by a sphere of radius >= r_k (an exact
  Euclidean distance transform per level, periodic in x, y), so D(p) = 2 max r_k.  The
  resolution in D is one histogram bin (2 dr).
* per z-bin: mean and median pore diameter over the void points, and a histogram of
  D (counts) the caller sums into regions (membrane interior, feed / permeate halves).
  The histogram is volume-weighted (every grid point is one volume element).
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage

POLYMER_TYPES = (1, 2)


def _axis(L, h):
    """n points spaced L/n ~ h over [0, L) (exact periodic wrap)."""
    n = max(int(round(L / h)), 1)
    return np.arange(n) * (L / n), L / n


def distance_field(box, types, xyz, z_lo, z_hi, h=0.25, sigma_bead=1.0, polymer_types=POLYMER_TYPES, workers=-1):
    """Distance from every grid point to the nearest polymer bead surface.
    -> dict(x, y, z (absolute grid coordinates), d [nx, ny, nz] float32, hx, hy, hz)."""
    types = np.asarray(types)
    xyz = np.asarray(xyz, float)
    keep = np.isin(types, polymer_types)
    if keep.sum() == 0:
        raise ValueError('no polymer beads in this frame')
    origin = np.array([box['x'][0], box['y'][0], box['z'][0]])
    L = np.array([box[k][1] - box[k][0] for k in 'xyz'])
    p = xyz[keep] - origin
    p[:, :2] %= L[:2]                                   # periodic in x, y
    p[:, 2] = np.clip(p[:, 2], 0.0, np.nextafter(L[2], 0.0))
    z0, z1 = max(z_lo - origin[2], 0.0), min(z_hi - origin[2], L[2])
    if z1 <= z0:
        raise ValueError(f'empty z window [{z_lo}, {z_hi}]')
    xs, hx = _axis(L[0], h)
    ys, hy = _axis(L[1], h)
    nz = max(int(round((z1 - z0) / h)), 1)
    zs = z0 + (np.arange(nz) + 0.5) * ((z1 - z0) / nz)
    hz = (z1 - z0) / nz
    # z is NOT periodic: a pseudo-period of 2 Lz never wraps a real pair (all z in [0, Lz))
    tree = cKDTree(p, boxsize=[L[0], L[1], 2.0 * L[2]])
    G = np.stack(np.meshgrid(xs, ys, zs, indexing='ij'), axis=-1).reshape(-1, 3)
    r, _ = tree.query(G, k=1, workers=workers)
    d = (r - 0.5 * sigma_bead).astype(np.float32).reshape(len(xs), len(ys), len(zs))
    return dict(x=xs + origin[0], y=ys + origin[1], z=zs + origin[2], d=d, hx=hx, hy=hy, hz=hz)


def covering_radius(d, h, r_probe=0.5, dr=0.125, r_cap=None):
    """Gelb-Gubbins covering on the grid.  R[p] = radius of the largest sphere
    (centred anywhere on the grid, radius d(c)) that contains p, for void points
    (d >= r_probe); NaN elsewhere.  Levels every `dr` from r_probe up to max d (or
    `r_cap`: spheres larger than that -- the reservoirs -- are counted at r_cap, which
    keeps the level count bounded); each level is one Euclidean distance transform
    (periodic in x, y).  h may be a 3-tuple (hx, hy, hz)."""
    d = np.asarray(d, np.float32)
    hv = np.broadcast_to(np.asarray(h, float), (3,))
    void = d >= r_probe
    R = np.where(void, d, np.nan).astype(np.float32)   # every void point holds its own sphere
    if not void.any():
        return R
    r_max = float(d.max()) if r_cap is None else min(float(d.max()), float(r_cap))
    R = np.minimum(R, np.float32(r_max))
    levels = np.arange(r_probe + dr, r_max + 1e-9, dr)
    for r in levels[::-1]:                              # largest spheres first
        S = d >= r
        if not S.any():
            continue
        pad = int(np.ceil(r / hv[:2].min())) + 1
        Sp = np.pad(~S, ((pad, pad), (pad, pad), (0, 0)), mode='wrap')
        dist = ndimage.distance_transform_edt(Sp, sampling=hv)[pad:-pad, pad:-pad, :]
        cov = void & (dist <= r) & ~(R >= r)
        R[cov] = r
    return R


def psd_frame(box, types, xyz, z_lo, z_hi, bin_edges, h=0.25, r_probe=0.5, d_edges=None,
              sigma_bead=1.0, polymer_types=POLYMER_TYPES):
    """Porosity and pore-size statistics of one frame on the caller's z-bins.
    bin_edges: the notebook's z-bin edges (absolute z).  d_edges: pore-DIAMETER
    histogram edges (default 2 r_probe .. 8 sigma in steps of 0.25).
    -> dict(zc, por, d_mean, d_med, n_grid, n_void, hist [nz, nd], d_edges, ts=None)
    with NaN in bins the grid does not cover."""
    F = distance_field(box, types, xyz, z_lo, z_hi, h, sigma_bead, polymer_types)
    d, zs = F['d'], F['z']
    if d_edges is None:
        d_edges = np.arange(2.0 * r_probe, 8.0 + 1e-9, 0.25)
    d_edges = np.asarray(d_edges, float)
    R = covering_radius(d, (F['hx'], F['hy'], F['hz']), r_probe, dr=0.5 * (d_edges[1] - d_edges[0]),
                        r_cap=0.5 * d_edges[-1])          # pores wider than the histogram range land in its last bin
    D = 2.0 * R                                          # pore diameter per void point (NaN = solid)
    edges = np.asarray(bin_edges, float)
    nz = len(edges) - 1
    zc = 0.5 * (edges[:-1] + edges[1:])
    bi = np.searchsorted(edges, zs, side='right') - 1     # z-bin of every grid plane
    por = np.full(nz, np.nan); dm = np.full(nz, np.nan); dmed = np.full(nz, np.nan)
    n_grid = np.zeros(nz, int); n_void = np.zeros(nz, int)
    hist = np.zeros((nz, len(d_edges) - 1), float)
    for b in np.unique(bi):
        if b < 0 or b >= nz:
            continue
        sl = bi == b
        Db = D[:, :, sl].ravel()
        v = Db[np.isfinite(Db)]
        n_grid[b], n_void[b] = Db.size, v.size
        por[b] = v.size / Db.size
        if v.size:
            dm[b], dmed[b] = float(v.mean()), float(np.median(v))
            hist[b] = np.histogram(np.clip(v, d_edges[0], d_edges[-1] - 1e-6), bins=d_edges)[0]
    return dict(zc=zc, por=por, d_mean=dm, d_med=dmed, n_grid=n_grid, n_void=n_void, hist=hist, d_edges=d_edges,
                z_grid=(float(zs[0]), float(zs[-1])), h=(F['hx'], F['hy'], F['hz']))


def region_density(hist, n_void, mask, d_edges):
    """Volume-weighted probability density of the pore diameter over the z-bins in
    `mask` (counts summed, normalised to unit area) -> (D centres, density); NaN
    density if the region holds no void point."""
    c = hist[mask].sum(axis=0)
    tot = c.sum()
    w = np.diff(d_edges)
    dens = c / (tot * w) if tot > 0 else np.full(len(c), np.nan)
    return 0.5 * (d_edges[:-1] + d_edges[1:]), dens
