#!/usr/bin/env python3
"""volfrac.py — single source of truth for solvent volume-fraction estimators
and the Voronoi→thermodynamic calibration (2026-08 plan).

Imported by calibration_analysis.ipynb AND both production notebooks
(triaxial_compression.ipynb, triaxial_permeation.ipynb) the same way
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

Calibration surface, joint fit over (φ_p, P) with the exact anchor λ(0,P)=1,
so reservoir bins are NEVER shifted:

    λ(φ_p, P) = 1 + (a₁ + b₁P)·φ_p + (a₂ + b₂P)·φ_p²

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
bulk_modulus_analysis.ipynb). Norm 'bin' divides by the geometric bin volume
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
    cont = tess.Container(p, limits=(tuple(origin), tuple(origin + L)), periodic=True)
    vol  = np.array([c.volume() for c in cont])
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
    {'a1','b1','a2','b2'}, 'covariance', 'raw_table', 'metadata'. Raises with
    a pointed message if the artifact has not been generated yet."""
    p = Path(path) if path is not None else CALIBRATION_JSON
    if not p.exists():
        raise FileNotFoundError(
            f"Calibration artifact not found: {p}\n"
            f"Run scripts/calibration_analysis.ipynb end-to-end to generate it "
            f"(see simulations/calibration_sweep/README.md).")
    with open(p) as f:
        calib = json.load(f)
    for k in ('a1', 'b1', 'a2', 'b2'):
        if k not in calib.get('coeffs', {}):
            raise KeyError(f"calibration json missing coeffs['{k}'] ({p})")
    return calib


def lambda_of(phi_p, P, calib):
    """λ(φ_p, P) = 1 + (a₁ + b₁P)·φ_p + (a₂ + b₂P)·φ_p².

    The exact anchor λ(0, P) = 1 is structural: a polymer-free (reservoir) bin
    is never shifted, and the correction scales with local polymer content."""
    c = calib['coeffs']
    phi_p = np.asarray(phi_p, float)
    return 1.0 + (c['a1'] + c['b1'] * P) * phi_p + (c['a2'] + c['b2'] * P) * phi_p**2


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
               normal_tol=0.05, aspect_tol=0.02, offdiag_tol=0.05):
    """The convergence gate (checked in ANALYSIS, never enforced in-run).

    Over the last tail_frac of the run:
      * normal-stress isotropy: max pairwise |⟨Pii⟩−⟨Pjj⟩| / |⟨P⟩| ≤ normal_tol
      * aspect-ratio stability: drift of lx/lz and ly/lz ≤ aspect_tol (relative)
      * off-diagonals: |⟨Pij⟩| / |⟨P⟩| ≤ offdiag_tol — FLAGGED only (no tri)

    Returns dict(passed, flags[list of str], stats). Production averages are
    only valid when 'passed'; off-diagonal flags do not fail the gate."""
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
            passed = False
            flags.append(f'aspect ratio {name} drift {100*drift:.1f}% > {100*aspect_tol:.0f}%')

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
