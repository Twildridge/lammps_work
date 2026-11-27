#!/usr/bin/env python3
"""volfrac.py — single source of truth for solvent volume-fraction estimators
and the Voronoi→thermodynamic calibration (2026-08 plan).

Imported by calibration_analysis.ipynb AND both production notebooks
(uniaxial_compression_single.ipynb, uniaxial_permeation_single.ipynb) the same way
cavity_widom.py used to be — a local import from scripts/:

    import sys; sys.path.insert(0, str(Path('lib').resolve()))   # or scripts/lib
    import volfrac

The math (verify against simulations/calibration_sweep/README.md)
-----------------------------------------------------------------
Thermodynamic solvent volume fraction — the one Flory–Rehner, Π, K(φ_p), D_c
are written in:

    φ_f = N_f v̄_f / V,    v̄_f = (∂V/∂N_f)_{T,P,N_p}

Euler (V extensive in N_f, N_p at fixed T,P):

    V = N_f v̄_f + N_p v̄_p   →   φ_p = 1 − φ_f exactly; no dry-reference
                                  volume; v̄_p by closure, never fit separately.

Per pressure the calibration fits ln⟨V⟩ vs ln N_f over the composition grid
and differentiates analytically:

    φ_f^th(N_f, P) = ∂ln⟨V⟩/∂ln N_f |_P
    v̄_f = φ_f^th · V / N_f ;    v̄_p = (V − N_f v̄_f)/N_p     (Euler closure)

The SPATIAL estimator is the periodic Voronoi pass in this module (Python
voro++ only — NO tessellation in any .lmp):

    φ_f^vor = Σ_{i∈solvent} v_i^voro / V        (per box, or per z-bin)

Calibration surface, per-pressure quadratics with the exact anchor λ(0,P)=1,
so reservoir bins are NEVER shifted; a(P), b(P) linearly interpolated between
sweep pressures (a single polynomial-in-P surface was too stiff to track the
isotherms):

    λ(φ_p, P) = 1 + a(P)·φ_p + b(P)·φ_p²

Application per z-bin, per frame (solvent primary; polymer by complement):

    φ_p^vor(z) = 1 − φ_f^vor(z)
    φ_f^cal(z) = min(1, λ(φ_p^vor(z), P_local) · φ_f^vor(z))
    φ_p^cal(z) = 1 − φ_f^cal(z)

Coefficients live in scripts/calibration/calibration_lambda.json, written by
calibration_analysis.ipynb (coefficients, covariance, full raw table,
metadata). load_calibration() reads it; lambda_of()/phi_calibrated() apply it.

Conventions
-----------
Monodisperse (plain, not radical) tessellation — all beads are σ=1. Mobile
types are (1, 2, 3) = crosslink, chain, solvent; walls (4, 5) are excluded by
default, matching the convention formerly referenced as
longitudinal_modulus_analysis.ipynb (that notebook is now
isotropic_compression.ipynb). Norm 'bin' divides by the geometric bin volume
(absolute, the W/V_solv-honest choice); 'mobile' divides by the summed cell
volume in the bin (exact saturation φ_f+φ_p=1 — use for the CALIBRATED φ).
"""

from pathlib import Path
import json

import numpy as np

SOLVENT_TYPE = 3
MOBILE_TYPES = (1, 2, 3)

# Default location of the calibration artifact, relative to this file
# (scripts/lib/volfrac.py → scripts/calibration/calibration_lambda.json).
CALIBRATION_JSON = Path(__file__).resolve().parent.parent / 'calibration' / 'calibration_lambda.json'


# ---------------------------------------------------------------------------
# Trajectory streaming
# ---------------------------------------------------------------------------
def stream_traj_frames(traj_file, want_ts=None):
    """{ts: (box, types[N], xyz[N,3])} — ALL atom types, streamed so a 100 MB
    dump is never held in memory. box = {'x': (lo,hi), 'y': ..., 'z': ...}.

    want_ts=None reads every frame (calibration dumps are small); otherwise
    only the requested timesteps are materialized."""
    want = None if want_ts is None else {int(t) for t in np.atleast_1d(want_ts)}
    out = {}
    with open(traj_file) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith('ITEM: TIMESTEP'):
                continue
            ts = int(f.readline()); f.readline(); n = int(f.readline()); f.readline()
            box = {k: tuple(map(float, f.readline().split()[:2])) for k in 'xyz'}
            cols = f.readline().split()[2:]
            if want is not None and ts not in want:
                for _ in range(n):
                    f.readline()
                continue
            ti, xi, yi, zi = (cols.index(c) for c in ('type', 'x', 'y', 'z'))
            a = np.empty((n, 4))
            for j in range(n):
                p = f.readline().split()
                a[j] = (float(p[ti]), float(p[xi]), float(p[yi]), float(p[zi]))
            out[ts] = (box, a[:, 0].astype(int), a[:, 1:])
    return out


# ---------------------------------------------------------------------------
# Binning and tessellation
# ---------------------------------------------------------------------------
def _frame_bins(box, bin_width):
    """(origin, L, edges, centers, widths, V_bin_k) on the frame's own z grid —
    origin = box zlo, width = bin_width, i.e. identical to the LAMMPS chunks."""
    origin = np.array([box['x'][0], box['y'][0], box['z'][0]])
    L      = np.array([box[k][1] - box[k][0] for k in 'xyz'])
    nbz    = int(np.ceil(L[2] / bin_width))
    edges  = np.minimum(origin[2] + np.arange(nbz + 1) * bin_width, origin[2] + L[2])
    widths = np.diff(edges)
    return origin, L, edges, 0.5 * (edges[:-1] + edges[1:]), widths, L[0] * L[1] * widths


# tess keeps every cell's full geometry (~50 kB per cell), so one container over a
# million-atom frame needs tens of GB.  Above TESS_BLOCK_MIN_N points the frame is
# tessellated in z-pencils instead: each pencil is its own container holding the pencil's
# points plus a halo of periodic neighbours in x and y (z stays periodic over the whole box),
# and only the volumes of the pencil's own points are kept.  The cells are the same as the
# single container's as long as no cell reaches further than TESS_HALO sideways.
TESS_BLOCK_MIN_N = 60000      # at or below this: one periodic container (the original path)
TESS_BLOCK_W     = 16.0       # target pencil width in x and y [sigma]
TESS_HALO        = 4.0        # halo around each pencil [sigma]
TESS_WORKERS     = 4          # pencils tessellated in parallel (each holds ~1.5 GB); 1 = serial


def _pencil_volumes(args):
    """volumes of the first n_own points of q (pencil-local coordinates; z periodic)."""
    q, lim_hi, zlo, n_own = args
    import tess
    pad = 1e-6                                       # voro++ rejects a point lying on a wall to round-off
    cont = tess.Container(q, limits=((-pad, -pad, zlo), (lim_hi[0] + pad, lim_hi[1] + pad, lim_hi[2])),
                          periodic=(False, False, True))
    return np.array([cont[i].volume() for i in range(n_own)])


def _cell_volumes_blocked(p, origin, L, w=None, halo=None, workers=None):
    """Voronoi cell volumes of the points p (wrapped into the periodic box origin, L), from
    z-pencils with a halo -- bounded memory, see the note above."""
    w = TESS_BLOCK_W if w is None else w
    halo = TESS_HALO if halo is None else halo
    workers = TESS_WORKERS if workers is None else workers
    nb = [max(1, int(L[k] // w)) for k in (0, 1)]
    bw = [L[k] / nb[k] for k in (0, 1)]
    if min(nb) < 3 or min(bw) <= 2 * halo:           # box too narrow for pencils + halo
        import tess
        cont = tess.Container(p, limits=(tuple(origin), tuple(origin + L)), periodic=True)
        return np.array([c.volume() for c in cont])
    rel = [(p[:, k] - origin[k]) % L[k] for k in (0, 1)]
    idx = [np.minimum((rel[k] / bw[k]).astype(int), nb[k] - 1) for k in (0, 1)]
    zlo, zhi = float(origin[2]), float(origin[2] + L[2])
    tasks, owners = [], []
    for i in range(nb[0]):
        dx = (rel[0] - i * bw[0] + halo) % L[0]          # pencil-local x, halo included
        inx = np.flatnonzero(dx < bw[0] + 2 * halo)
        for j in range(nb[1]):
            dy = (rel[1][inx] - j * bw[1] + halo) % L[1]
            s = inx[dy < bw[1] + 2 * halo]
            own = (idx[0][s] == i) & (idx[1][s] == j)
            if not own.any():
                continue
            s = np.concatenate([s[own], s[~own]])        # the pencil's own points first
            q = np.column_stack([dx[s], (rel[1][s] - j * bw[1] + halo) % L[1], p[s, 2]])
            tasks.append((q, (bw[0] + 2 * halo, bw[1] + 2 * halo, zhi), zlo, int(own.sum())))
            owners.append(s[:int(own.sum())])
    res = None
    if workers > 1 and len(tasks) > 1:
        try:
            import multiprocessing as mp
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn')) as ex:
                res = list(ex.map(_pencil_volumes, tasks))
        except Exception as e:
            print(f'    (parallel tessellation unavailable: {type(e).__name__}: {e} -> serial)')
    if res is None:
        res = [_pencil_volumes(t) for t in tasks]
    vol = np.full(len(p), np.nan)
    for s, v in zip(owners, res):
        vol[s] = v
    return vol


def _tessellate(box, types, xyz, mobile_only=True, mobile_types=MOBILE_TYPES):
    """Periodic voro++ tessellation of one frame → (points, types_kept, volumes,
    origin, L). Shared by the per-bin and whole-box estimators so the two can
    never drift apart. Drops exact duplicate positions (voro++ fails on them)
    and warns if the cell volumes do not close the box."""
    import tess
    origin = np.array([box['x'][0], box['y'][0], box['z'][0]])
    L      = np.array([box[k][1] - box[k][0] for k in 'xyz'])
    keep = np.isin(types, mobile_types) if mobile_only else np.ones(len(types), bool)
    t_k  = np.asarray(types)[keep]
    p    = origin + (np.asarray(xyz)[keep] - origin) % L        # wrap into primary box
    _, ui = np.unique(np.round(p, 6), axis=0, return_index=True)
    if len(ui) < len(p):
        print(f'    warning: dropping {len(p)-len(ui)} duplicate positions')
        ui = np.sort(ui); p, t_k = p[ui], t_k[ui]
    if len(p) <= TESS_BLOCK_MIN_N:
        cont = tess.Container(p, limits=(tuple(origin), tuple(origin + L)), periodic=True)
        vol  = np.array([c.volume() for c in cont])
        del cont
    else:
        vol = _cell_volumes_blocked(p, origin, L)
    _err = abs(vol.sum() - L.prod()) / L.prod()
    if _err > 0.01:
        print(f'    warning: sum(V_cell) off box volume by {100*_err:.1f}%')
    return p, t_k, vol, origin, L


def phi_voronoi_frame(box, types, xyz, bin_width, mobile_only=True, norm='bin',
                      solvent_type=SOLVENT_TYPE, mobile_types=MOBILE_TYPES):
    """Solvent Voronoi volume fraction per z-bin for one frame → (centers, φ_s).

    Full 3-D PERIODIC tessellation via voro++. Makes no incompressibility
    assumption — it measures the space the polymer actually leaves the solvent.
    norm='bin'    → φ_s = Σ V_cell(solvent) / V_bin              [absolute]
    norm='mobile' → φ_s = Σ V_cell(solvent) / Σ V_cell(in bin)   [exact
                     saturation; the choice for the CALIBRATED φ]
    norm='both'   → (centers, φ_bin, φ_mobile) from ONE tessellation."""
    p, t_k, vol, origin, L = _tessellate(box, types, xyz, mobile_only, mobile_types)
    _, _, edges, centers, widths, V_bin_k = _frame_bins(box, bin_width)
    nbz = len(centers)
    bi  = np.clip(((p[:, 2] - origin[2]) / bin_width).astype(int), 0, nbz - 1)
    s   = (t_k == solvent_type)
    V_s = np.bincount(bi[s], weights=vol[s], minlength=nbz)[:nbz]
    V_t = np.bincount(bi,    weights=vol,    minlength=nbz)[:nbz]
    den_mob = np.where(V_t > 0, V_t, np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        if norm == 'both':
            return centers, V_s / V_bin_k, V_s / den_mob
        return centers, V_s / (V_bin_k if norm == 'bin' else den_mob)


def phi_voronoi_box(box, types, xyz, mobile_only=True,
                    solvent_type=SOLVENT_TYPE, mobile_types=MOBILE_TYPES):
    """Whole-box solvent Voronoi fraction for one frame → scalar φ_f^vor.

    The calibration-sweep estimator: one number per (homogeneous, fully
    periodic) box, φ_f^vor = Σ_{i∈solvent} v_i / V_box. With mobile_only and a
    box containing only types 1,2,3 the tessellation closes the box exactly,
    so this is identical to the 'mobile' norm."""
    p, t_k, vol, origin, L = _tessellate(box, types, xyz, mobile_only, mobile_types)
    return float(vol[t_k == solvent_type].sum() / L.prod())


def phi_voronoi_traj(traj_file, want_ts=None, verbose=True, **kw):
    """Stream a calibration dump → (ts array, φ_f^vor per frame array).

    Convenience wrapper for the calibration notebook: mean ± SE over the
    calib_frames+1 near-end frames is the per-run Voronoi measurement."""
    frames = stream_traj_frames(traj_file, want_ts)
    ts = sorted(frames)
    phis = []
    for t in ts:
        box, typ, xyz = frames[t]
        phis.append(phi_voronoi_box(box, typ, xyz, **kw))
        if verbose:
            print(f'    ts {t}: phi_f^vor = {phis[-1]:.4f}')
    return np.array(ts, float), np.array(phis)


# ---------------------------------------------------------------------------
# Calibration: load and apply λ(φ_p, P)
# ---------------------------------------------------------------------------
def load_calibration(path=None):
    """Read calibration_lambda.json → dict with (at least) 'coeffs':
    {'pressures','a','b'} (per-sweep-pressure quadratic coefficients),
    'per_pressure', 'raw_table', 'metadata'. Raises with a pointed message if
    the artifact has not been generated yet or uses the retired global-surface
    schema."""
    p = Path(path) if path is not None else CALIBRATION_JSON
    if not p.exists():
        raise FileNotFoundError(
            f"Calibration artifact not found: {p}\n"
            f"Run scripts/calibration_analysis.ipynb end-to-end to generate it "
            f"(see simulations/calibration_sweep/README.md).")
    with open(p) as f:
        calib = json.load(f)
    coeffs = calib.get('coeffs', {})
    if 'a1' in coeffs:
        raise KeyError(
            f"calibration json uses the retired global-surface schema (a1..b2); "
            f"re-run scripts/calibration_analysis.ipynb to regenerate ({p})")
    for k in ('pressures', 'a', 'b'):
        if k not in coeffs:
            raise KeyError(f"calibration json missing coeffs['{k}'] ({p})")
    return calib


def lambda_of(phi_p, P, calib):
    """λ(φ_p, P) = 1 + a(P)·φ_p + b(P)·φ_p².

    a, b are fit per sweep pressure (the global polynomial-in-P surface was
    too stiff — reduced χ² ≈ 35 and the P=0.5 isotherm missed entirely) and
    linearly interpolated in P here; P outside the sweep range clamps to the
    nearest sweep pressure. P may be a scalar or an array broadcastable
    against phi_p (per-bin P_local).

    The exact anchor λ(0, P) = 1 is structural: a polymer-free (reservoir) bin
    is never shifted, and the correction scales with local polymer content.
    In φ_p the quadratic extrapolates beyond the sweep's data window
    (φ_p ≈ 0.55–0.80), pinned to the anchor at φ_p = 0."""
    c = calib['coeffs']
    Ps = np.asarray(c['pressures'], float)
    a = np.interp(P, Ps, np.asarray(c['a'], float))
    b = np.interp(P, Ps, np.asarray(c['b'], float))
    phi_p = np.asarray(phi_p, float)
    return 1.0 + a * phi_p + b * phi_p**2


def phi_calibrated(phi_f_vor, P, calib):
    """Apply the calibration to a Voronoi solvent fraction (scalar or array):

        φ_p^vor = 1 − φ_f^vor
        φ_f^cal = min(1, λ(φ_p^vor, P) · φ_f^vor)

    P may be a scalar or an array broadcastable against phi_f_vor (per-bin
    P_local). Returns φ_f^cal; the polymer fraction is 1 − φ_f^cal, always by
    complement, never fit independently. NaNs pass through."""
    phi_f_vor = np.asarray(phi_f_vor, float)
    lam = lambda_of(1.0 - phi_f_vor, P, calib)
    return np.minimum(1.0, lam * phi_f_vor)


# ---------------------------------------------------------------------------
# Calibration-run tabular readers (box_dimensions / pressure_tensor)
# ---------------------------------------------------------------------------
def read_box_dimensions(path):
    """box_dimensions_*.dat → (step, lx, ly, lz) arrays. Columns: step lx ly lz."""
    d = np.loadtxt(path, comments='#', ndmin=2)
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]


def read_pressure_tensor(path):
    """pressure_tensor_*.dat → (step, P[ n,6 ]) with columns pxx pyy pzz pxy pxz pyz
    (window-averaged every volume_freq steps; the aniso-gate input)."""
    d = np.loadtxt(path, comments='#', ndmin=2)
    return d[:, 0], d[:, 1:7]


def aniso_gate(steps, ptens, lx, ly, lz, tail_frac=0.5,
               normal_tol=0.05, aspect_tol=0.10, offdiag_tol=0.05):
    """The convergence gate (checked in ANALYSIS, never enforced in-run).

    Over the last tail_frac of the run:
      * normal-stress isotropy: max pairwise |⟨Pii⟩−⟨Pjj⟩| / |⟨P⟩| ≤ normal_tol
      * aspect-ratio drift of lx/lz and ly/lz > aspect_tol — FLAGGED only:
        under `fix npt aniso` the box shape is a soft/diffusive mode (a pure
        fluid has no shear restoring force at all, so lx/lz random-walks tens
        of percent while ⟨V⟩ stays converged); shape drift alone therefore
        cannot invalidate a volume average
      * off-diagonals: |⟨Pij⟩| / |⟨P⟩| ≤ offdiag_tol — FLAGGED only (no uni)

    Returns dict(passed, flags[list of str], stats). Production averages are
    only valid when 'passed'; the hard criterion is normal-stress isotropy —
    aspect-ratio and off-diagonal flags do not fail the gate."""
    n  = len(steps)
    k0 = int(n * (1 - tail_frac))
    pt = ptens[k0:]
    diag = pt[:, :3].mean(axis=0)
    off  = pt[:, 3:].mean(axis=0)
    pbar = diag.mean()
    flags, passed = [], True

    span = np.ptp(diag) / abs(pbar) if pbar else np.inf
    if span > normal_tol:
        passed = False
        flags.append(f'normal-stress anisotropy {100*span:.1f}% > {100*normal_tol:.0f}%')

    m  = min(len(lx), n)
    ax = (lx[:m] / lz[:m])[int(m * (1 - tail_frac)):]
    ay = (ly[:m] / lz[:m])[int(m * (1 - tail_frac)):]
    for name, a in (('lx/lz', ax), ('ly/lz', ay)):
        drift = np.ptp(a) / a.mean()
        if drift > aspect_tol:
            flags.append(f'aspect ratio {name} drift {100*drift:.1f}% > {100*aspect_tol:.0f}% — flagged only')

    for name, v in zip(('Pxy', 'Pxz', 'Pyz'), off):
        if abs(v) / abs(pbar) > offdiag_tol:
            flags.append(f'off-diagonal {name} = {v:.4f} ({100*abs(v)/abs(pbar):.1f}% of P) — flagged only')

    return {'passed': passed, 'flags': flags,
            'stats': {'diag_mean': diag.tolist(), 'offdiag_mean': off.tolist(),
                      'pbar': float(pbar)}}


def block_average_volume(steps, lx, ly, lz, tail_frac=0.5, n_blocks=5):
    """⟨V⟩ ± SE from the box_dimensions trace: last tail_frac of the run split
    into n_blocks non-overlapping blocks (block means → mean and standard
    error). Verify the tail is equilibrated from the trace itself — the gate
    and a look at V(t) — don't assume."""
    V  = lx * ly * lz
    k0 = int(len(V) * (1 - tail_frac))
    tail = V[k0:]
    blocks = np.array_split(tail, n_blocks)
    bm = np.array([b.mean() for b in blocks if len(b)])
    return float(bm.mean()), float(bm.std(ddof=1) / np.sqrt(len(bm)))
