"""Cross-deck comparison -- figures built from MORE THAN ONE deck's run (2026-10-03).

Notebook: scripts/cross_deck_comparison.ipynb.  Each deck is loaded with its own library
(lib/triaxial.py for the two-piston compression and permeation decks, lib/shear.py for
shear_slab) through `load_deck`, which returns a plain dict
    dict(kind='compression' | 'permeation' | 'shear', cfg, R, levels | P)
so a comparison function takes the decks it needs and nothing here re-implements a loader.

First comparison: D_c and kappa as functions of the solvent volume fraction phi_s.

  Compression sweep (uniform states).  Level i is a uniformly compressed gel at phi_s,i with the
  network stress sigma'_i and the displacement-strain secant modulus M_i = dsigma'_i/eps_disp.
  The levels give the equilibrium curve sigma'(phi_s) (monotone PCHIP, Monte-Carlo CI), which the
  closures and fig_states (b) use.  Their hold-relaxation D_c (tri.fit_Dc) is a LOAD measurement
  and is NOT used here (2026-10-21): the hold starts from the ramp's face skins plus an undrained
  interior strain (the 2026-10-20 amplitude-spectrum analysis), and its D_c moves with the mode
  count, trim and window (0.05-0.085 on run 7).  It stays in the compression notebooks, labelled LOAD.

  Unload (free re-swelling after a hold; deck Phase 2c, files _u<lvl>).  The ONLY compression-side
  D_c in this comparison: the re-swelling profile fit of every level that has an unload
  (unload_states), with its Rg thickness-trace D_c as a check, and kappa = D_c / M_i with that
  level's secant M.  Closure (ii) uses this D_c (constant at one level, interpolated over several).

  Permeation (one run, a range of phi_s at once).  In the steady state the network is at rest and
  the solvent flux q = Q/A is the same at every z.  Darcy, the uniform total stress
  (dp/dz = dsigma'/dz) and a network stress that depends on z through phi_s alone give
  q = kappa |dsigma'/dphi_s| |dphi_s/dz|; with phi_p L = const (polymer conserved, lateral box
  fixed) a strain increment is dphi_p/phi_p, so kappa phi_p |dsigma'/dphi_s| = kappa M = D_c and

      D_c(phi_s) = q phi_p / |dphi_s/dz| ,        kappa(phi_s) = D_c(phi_s) / M_perm .

  D_c needs the flux and the slope of the steady phi_s(z) profile only -- no modulus, no fit of a
  transient.  M_perm is the permeation run's own displacement-derived modulus dP L/(2|u_F|) (a
  membrane average), so the local kappa carries the local D_c and one M.  The membrane-average
  numbers of the permeation notebook (Darcy kappa = Q L/(A dP), the consolidation-fit D_c) are
  drawn as bands over the phi_s range the membrane spans.

  No modulus here comes from the volume fraction.  The tangent modulus of the fitted sigma'(phi_s)
  curve, M_cur = phi_p dsigma'/dphi_p (2026-10-03), was dropped on 2026-10-15: it rests on the
  volume-fraction profiles, which are not trusted to that level yet.  Every M is a displacement-
  strain secant of its own deck.

  Closures (print_summary, fig_closure): the network stress the compression curve assigns to the
  phi_s at the support against the applied dP (it must carry all of it), and the flux the
  compression D_c(phi_s) predicts for the measured profile against the measured one.
"""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.interpolate import PchipInterpolator

import triaxial as tri
from triaxial import WONG, sig, smart_legend

PHI_LABEL = {'cal': r'$\lambda$-calibrated Voronoi', 'vor': 'Voronoi', 'mf': 'mass fraction'}


@dataclass
class Config:
    """Knobs of the comparison (the decks keep their own Config objects)."""
    NAME: str = 'comparison'              # plots go to <base_dir>/plots/comparison/<NAME>
    base_dir: str = '../../flow_data_local'
    PHI_KIND: str = 'cal'                 # which solvent volume fraction: 'cal' (lambda-calibrated Voronoi) | 'vor' | 'mf'
    PERM_TRIM: float = 6.0                # sigma dropped at each end of the membrane before the phi_s(z) slope is taken
                                          # (the contact layer on the support, the diffuse feed face)
    PERM_N_SEG: int = 3                   # the profile is cut into this many z-segments, one local D_c / kappa each (the
                                          # Voronoi-based profiles come from the few steady trajectory frames: more
                                          # segments than ~3 leave the slope of a segment unresolved; 'mf' carries 4-5)
    PERM_POLY_DEG: int = 2                # degree of the polynomial phi_s(z) behind the continuous D_c(phi_s) curve
    N_MC: int = 400                       # Monte-Carlo draws for the CI of the compression curve sigma'(phi_s)
    # ---- the solvent-measures figure of the permeation membrane (fig_solvent_measures, 2026-10-08) ----
    CONC_MARGIN: float = 5.0              # sigma of reservoir kept on either side of the membrane (the reservoir values of
                                          # the mass fraction -- > 1 on the feed side -- stay visible, the rest is trimmed)
    CONC_FRAMES: bool = True              # draw the individual steady Voronoi frames of c_s faintly (the frame-to-frame noise
                                          # under the 95 % band; the reference keeps its band only, its frames are not stored)
    CONC_DIFF_PANEL: bool = True          # a second panel with the two other measures minus the calibrated fraction
    ci_level: float = 0.95
    save: bool = True

    @property
    def PLOT_DIR(self):
        return Path(self.base_dir) / 'plots' / 'comparison' / self.NAME


def _save(fig, ccfg, stem):
    if ccfg.save:
        ccfg.PLOT_DIR.mkdir(parents=True, exist_ok=True)
        out = ccfg.PLOT_DIR / f'{stem}.png'
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print('saved', out)
    plt.show()
    return fig


# ===========================================================================
#  1. LOADING
# ===========================================================================
def sync_deck(cfg, force=False):
    """Pull one deck's files from Expanse with its own library's sync."""
    if _is_shear(cfg):
        import shear
        return shear.sync_from_expanse(cfg, force=force)
    return tri.sync_from_expanse(cfg, force=force)


def _is_shear(cfg):
    return type(cfg).__module__ == 'shear'


def load_deck(cfg, verbose=True):
    """Load one run with the library its Config belongs to.
    tri.Config(mode='compression') -> dict(kind='compression', cfg, R, levels)   (M, G, D_c, phi per level)
    tri.Config(mode='permeation')  -> dict(kind='permeation',  cfg, R, P)        (flux, phi profile, D_c, M)
    shear.Config                   -> dict(kind='shear',       cfg, R, levels)   (G, N1/N2, D_c per level)"""
    if _is_shear(cfg):
        import shear
        R = shear.load_reference(cfg)
        levels = [L for L in (shear.load_level(cfg, R, lvl) for lvl in cfg.STRAINS) if L is not None]
        levels.sort(key=lambda L: L['gamma'])
        return dict(kind='shear', cfg=cfg, R=R, levels=levels)
    R = tri.load_reference(cfg)
    if cfg.mode == 'permeation':
        P = tri.load_permeation(cfg, R, verbose=verbose)
        tri.add_perm_volume_fractions(cfg, R, P)
        if cfg.PSD_ENABLE and cfg.VOR_ENABLE:
            tri.add_perm_psd(cfg, R, P)                       # probe porosity (fig_solvent_measures); cfg.PSD_R_PROBE
        tri.add_perm_displacement(cfg, R, P, verbose=verbose)
        return dict(kind='permeation', cfg=cfg, R=R, P=P)
    levels = [L for L in (tri.load_level(cfg, R, lvl, verbose=verbose) for lvl in cfg.COMP_LEVELS) if L is not None]
    levels.sort(key=lambda L: L['eps'])
    tri.add_volume_fractions(cfg, R, levels)
    return dict(kind='compression', cfg=cfg, R=R, levels=levels)


# ===========================================================================
#  2. COMPRESSION: uniform states -> sigma'(phi_s) (equation of state); the hold D_c is NOT used (2026-10-21)
# ===========================================================================
def _interior_mean(prof, mask, ci):
    """(mean, CI half-width) of a (mean, lo, hi) bin profile over the interior bins: the bins'
    scatter about the interior mean (t-interval) and the mean per-bin frame CI, in quadrature."""
    m = np.asarray(prof[0], float)[mask]
    m = m[np.isfinite(m)]
    h_bins = stats.t.ppf(0.5 + ci / 2, max(len(m) - 1, 1)) * np.std(m, ddof=1) / np.sqrt(len(m)) if len(m) > 1 else 0.0
    h_frames = float(np.nanmean((np.asarray(prof[2], float) - np.asarray(prof[1], float))[mask]) / 2) / np.sqrt(max(len(m), 1))
    return float(np.mean(m)), float(np.hypot(h_bins, h_frames))


def compression_states(ccfg, comp):
    """The sweep as a table of uniform states: the reference (sigma' = 0) and one row per level.
    -> dict of arrays (phi, phi_h, sig, sig_h, eps, M_sec, M_sec_lo, M_sec_hi, label) ordered by increasing
    compression, row 0 = the reference.  M_sec is the level's displacement-strain secant dsigma'/eps_disp
    (L['M_net']).  No D_c here (2026-10-21): the levels' hold-relaxation D_c is a LOAD measurement and stays
    in the sweep notebook; the comparison's compression-side D_c comes from the unload levels only
    (unload_states)."""
    key = 'phi_' + ccfg.PHI_KIND
    R = comp['R']
    if R.get(key) is None:
        raise ValueError(f"the compression reference has no {key} profile (VOR_ENABLE / calibration artifact?)")
    p0, h0 = _interior_mean(R[key], R['interior'], ccfg.ci_level)
    rows = [dict(label='ref', phi=p0, phi_h=h0, sig=0.0, sig_h=0.0, eps=0.0, M_sec=np.nan, M_sec_lo=np.nan, M_sec_hi=np.nan)]
    for L in comp['levels']:
        if L.get(key) is None:
            print(f"  NOTE: level {L['lvl']} has no {key} profile -> left out of the comparison")
            continue
        p, h = _interior_mean(L[key], L['interior'], ccfg.ci_level)
        rows.append(dict(label=str(L['lvl']), phi=p, phi_h=h, sig=float(L['dsig_net']),
                         sig_h=float((L['M_net_hi'] - L['M_net_lo']) / 2 * L['eps_M']), eps=float(L['eps_M']),
                         M_sec=float(L['M_net']), M_sec_lo=float(L['M_net_lo']), M_sec_hi=float(L['M_net_hi'])))
    S = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    if not np.all(np.diff(S['phi']) < 0):
        raise ValueError('phi_s does not fall monotonically along the sweep: ' + ', '.join(f'{p:.3f}' for p in S['phi']))
    return S


def _curve(phi_s, sg):
    """Monotone (PCHIP) sigma'(phi_p) through the states -> the sigma'(phi_s) callable."""
    f = PchipInterpolator(1.0 - np.asarray(phi_s, float), np.asarray(sg, float), extrapolate=True)
    return lambda ps: f(1.0 - np.asarray(ps, float))


def compression_curve(ccfg, S, seed=0):
    """sigma'(phi_s) from the states S, with a Monte-Carlo CI (the levels' sigma' and phi_s drawn from
    their CIs, the curve rebuilt each time).  Used by the closures (permeation_local,
    fig_closure) and drawn in fig_states (b).  -> dict(sig(phi), sig_ci(phi) -> (lo, hi), n_draws)."""
    f_sig = _curve(S['phi'], S['sig'])
    rng = np.random.default_rng(seed)
    z = stats.norm.ppf(0.5 + ccfg.ci_level / 2)
    draws = []
    for _ in range(int(ccfg.N_MC)):
        ph = S['phi'] + rng.normal(0.0, 1.0, len(S['phi'])) * S['phi_h'] / z
        sg = S['sig'] + rng.normal(0.0, 1.0, len(S['sig'])) * S['sig_h'] / z
        if np.all(np.diff(ph) < 0):
            draws.append(_curve(ph, np.maximum.accumulate(sg)))          # sigma' is non-decreasing with compression
    q = [100 * (0.5 - ccfg.ci_level / 2), 100 * (0.5 + ccfg.ci_level / 2)]

    def sig_ci(ps):
        if not draws:
            v = f_sig(ps)
            return v, v
        A = np.array([d(ps) for d in draws])
        return tuple(np.percentile(A, q, axis=0))

    return dict(sig=f_sig, sig_ci=sig_ci, n_draws=len(draws))


# ===========================================================================
#  3. PERMEATION: the steady phi_s(z) profile + the flux -> D_c(phi_s), kappa(phi_s)
# ===========================================================================
def permeation_local(ccfg, perm, C):
    """Local D_c and kappa along the steady membrane (see the module docstring).
    -> dict: z, phi, phi_lo, phi_hi (the steady profile inside the fit range), q, q_h, dP,
       poly (phi_s(z) polynomial), zz / phi_fit / Dc_fit / kappa_fit (the continuous curves),
       seg = dict(z, phi, slope, slope_h, Dc, Dc_lo, Dc_hi, kappa, kappa_lo, kappa_hi) (one row per segment; the
       bounds come from the slope's regression CI and the flux CI, and the upper ones are inf when the slope CI reaches 0),
       phi_bot / phi_top (the fitted phi_s at the two ends of the fit range), the membrane-average
       numbers of the permeation notebook (kappa_darcy, Dc_prof, Dc_trace, M_perm) and the closures.
       Every local kappa is the local D_c over M_perm, the run's own displacement-derived modulus."""
    P, R, cfg = perm['P'], perm['R'], perm['cfg']
    key = 'phi_' + ccfg.PHI_KIND
    if P.get(key) is None:
        raise ValueError(f"the permeation run has no steady {key} profile (VOR_ENABLE / calibration artifact / trajectories?)")
    fl = P.get('flux') or {}
    if 'Q_N' not in fl:
        raise ValueError('the permeation run has no N-slope flux (permeate_count / permeation file)')
    NS = fl['Q_N']
    A = P['area']
    q, q_h = NS['mean'] / A, 0.5 * (NS['hi'] - NS['lo']) / A
    dP = fl.get('dP_meas', fl.get('dP_app'))
    z = np.asarray(P['z'], float)
    m, lo, hi = (np.asarray(a, float) for a in P[key])
    sel = P['in_mem'] & (z >= P['z_mem_lo'] + ccfg.PERM_TRIM) & (z <= P['z_mem_hi'] - ccfg.PERM_TRIM) & np.isfinite(m)
    zs, ps = z[sel], m[sel]
    if len(zs) < 4 * ccfg.PERM_N_SEG:
        raise ValueError(f'only {len(zs)} membrane bins inside the fit range -- lower PERM_N_SEG or PERM_TRIM')
    out = dict(z=zs, phi=ps, phi_lo=lo[sel], phi_hi=hi[sel], q=float(q), q_h=float(q_h), dP=float(dP), A=float(A),
               z_lo=float(zs[0]), z_hi=float(zs[-1]), kind=ccfg.PHI_KIND,
               phi_ref=_interior_mean(R[key], R['interior'], ccfg.ci_level)[0] if R.get(key) is not None else np.nan)
    # ---- the permeation notebook's membrane averages; its M divides every local kappa below ----
    kd = (fl.get('k') or {}).get('N_measured') or (fl.get('k') or {}).get('N_applied')
    out['kappa_darcy'] = (kd['k'], kd['lo'], kd['hi']) if kd else None
    F = P.get('Dc')
    out['Dc_prof'] = float(F['Dc']) if F else np.nan
    out['Dc_trace'] = float(F['trace']['Dc']) if (F and F.get('trace')) else np.nan
    Mp = P.get('M_perm')
    # the permeation M is displacement-derived: dP L/(2|u_F|) from the steady free-face displacement (perm_M_estimates),
    # a membrane average -- the only modulus this run measures (no volume-fraction-derived modulus since 2026-10-15)
    has_M = bool(Mp and Mp.get('primary'))
    M_perm = float(Mp[Mp['primary']]['M']) if has_M else np.nan
    M_lo, M_hi = (float(Mp[Mp['primary']]['lo']), float(Mp[Mp['primary']]['hi'])) if has_M else (np.nan, np.nan)
    out['M_perm'], out['M_perm_ci'], out['M_perm_how'] = M_perm, (M_lo, M_hi), (Mp['primary'] if has_M else None)
    if not has_M:
        print("  NOTE: the permeation run has no displacement-derived M (P['M_perm']) -> the local kappa is left nan")
    # ---- continuous: polynomial phi_s(z) ----
    poly = np.polynomial.Polynomial.fit(zs, ps, ccfg.PERM_POLY_DEG)
    zz = np.linspace(zs[0], zs[-1], 200)
    pf, dpf = poly(zz), poly.deriv()(zz)
    with np.errstate(divide='ignore', invalid='ignore'):
        Dc_fit = np.where(dpf > 0, q * (1.0 - pf) / dpf, np.nan)
    out.update(poly=poly, zz=zz, phi_fit=pf, Dc_fit=Dc_fit, kappa_fit=Dc_fit / M_perm,
               phi_bot=float(pf[0]), phi_top=float(pf[-1]))
    # ---- segments: a straight line per z-segment, slope +/- its regression CI ----
    seg = {k: [] for k in ('z', 'phi', 'slope', 'slope_h', 'Dc', 'Dc_lo', 'Dc_hi', 'kappa', 'kappa_lo', 'kappa_hi')}
    for idx in np.array_split(np.arange(len(zs)), ccfg.PERM_N_SEG):
        r = stats.linregress(zs[idx], ps[idx])
        t = stats.t.ppf(0.5 + ccfg.ci_level / 2, max(len(idx) - 2, 1))
        pm = float(np.mean(ps[idx]))
        sh = t * r.stderr
        D = q * (1.0 - pm) / r.slope if r.slope > 0 else np.nan
        D_lo = (q - q_h) * (1.0 - pm) / (r.slope + sh) if r.slope > 0 else np.nan
        D_hi = (q + q_h) * (1.0 - pm) / (r.slope - sh) if r.slope > sh else np.inf
        for k, v in (('z', float(np.mean(zs[idx]))), ('phi', pm), ('slope', r.slope), ('slope_h', sh), ('Dc', D), ('Dc_lo', D_lo), ('Dc_hi', D_hi),
                     ('kappa', D / M_perm), ('kappa_lo', D_lo / M_hi), ('kappa_hi', D_hi / M_lo)):
            seg[k].append(float(v))
    out['seg'] = {k: np.array(v) for k, v in seg.items()}
    # ---- closures ----
    #  (i) the network stress across the fit range, read off the compression curve, against the share of dP
    #      that falls inside the fit range if the pore pressure were linear in z
    out['sig_z'] = C['sig'](pf)                                        # sigma'(z) the compression curve assigns to phi_s(z)
    out['dsig'] = float(out['sig_z'][0] - out['sig_z'][-1])
    L_mem = float(P['z_mem_hi'] - P['z_mem_lo'])
    out['L_mem'], out['L_fit'] = L_mem, float(zs[-1] - zs[0])
    out['kappa_local_mean'] = float(q * out['L_fit'] / out['dsig']) if out['dsig'] > 0 else np.nan   # harmonic mean of kappa(z)
    return out


def flux_from_unload(RW, PL):
    """The flux the UNLOAD D_c predicts for the measured steady profile (closure ii): q L = integral of
    D_c/(1 - phi_s) dphi_s over the phi_s range of the fit.  One unload level: D_c constant at its value;
    several: linear in phi_s between them, held at the end values.  -> (q_pred, number of unload levels used);
    (nan, 0) without an unload."""
    if RW is None or not np.isfinite(RW['Dc']).any():
        return np.nan, 0
    ok = np.isfinite(RW['Dc'])
    order = np.argsort(RW['phi'][ok])
    ph, Dc = RW['phi'][ok][order], RW['Dc'][ok][order]
    x = np.linspace(PL['phi_bot'], PL['phi_top'], 400)
    D = np.interp(x, ph, Dc) if len(ph) > 1 else np.full_like(x, Dc[0])
    integ = float(np.sum(0.5 * (D[1:] / (1 - x[1:]) + D[:-1] / (1 - x[:-1])) * np.diff(x)))
    return integ / PL['L_fit'], int(ok.sum())


def unload_states(ccfg, decks):
    """Every compression level WITH an unload (deck Phase 2c, files _u<lvl>) across the given compression decks: the
    free re-swelling profile-fit D_c -- the trusted D_c and the only compression-side D_c of this comparison
    (2026-10-21) -- its Rg thickness-trace D_c as a check, the level's phi_s, its displacement-strain secant M and
    kappa = D_c/M (CI from the M CI).  -> dict of arrays (run_id, lvl, eps, phi, phi_h, Dc, Dc_R2, tau1,
    Dc_trace_rg, M_sec, M_sec_lo, M_sec_hi, M_pist, kappa, kappa_lo, kappa_hi) ordered by decreasing phi_s, or None
    when no level has an unload.  The hold D_c of these levels is NOT carried: it is a load measurement (the
    compression notebooks keep it)."""
    key = 'phi_' + ccfg.PHI_KIND
    rows, seen = [], set()
    for d in decks:
        if d is None or d.get('kind') != 'compression' or d['cfg'].RUN_ID in seen:
            continue
        seen.add(d['cfg'].RUN_ID)
        for L in d['levels']:
            Pf = (L.get('unload') or {}).get('prof')
            if Pf is None or L.get(key) is None:
                continue
            p, h = _interior_mean(L[key], L['interior'], ccfg.ci_level)
            Dc = float(Pf['Dc'])
            M, Mlo, Mhi = float(L['M_net']), float(L['M_net_lo']), float(L['M_net_hi'])
            Tr = (L['unload'] or {}).get('trace_rg')
            rows.append(dict(run_id=d['cfg'].RUN_ID, lvl=str(L['lvl']), eps=float(L['eps_M']), phi=float(p), phi_h=float(h),
                             Dc=Dc, Dc_R2=float(Pf['R2']), tau1=float(Pf['tau1']), Dc_trace_rg=float(Tr['Dc']) if Tr else np.nan,
                             M_sec=M, M_sec_lo=Mlo, M_sec_hi=Mhi, M_pist=float(L.get('M_pist', np.nan)),
                             kappa=Dc / M, kappa_lo=Dc / Mhi, kappa_hi=Dc / Mlo))
    if not rows:
        return None
    rows.sort(key=lambda r: -r['phi'])
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def phi_dependence(ccfg, comp, perm, unload=None):
    """Everything for the D_c(phi_s) / kappa(phi_s) comparison: -> dict(S, C, PL, RW, q_pred, n_unload).
    comp: the compression sweep (the states and the equation of state sigma'(phi_s); its hold D_c is not used);
    perm: the permeation run;  unload: a compression deck, or a list of them, whose levels have an unload -- the
    ONLY source of a compression-side D_c here (2026-10-21).  Levels of `comp` with an unload count too."""
    decks = [comp] + ([unload] if isinstance(unload, dict) else list(unload or []))
    S = compression_states(ccfg, comp)
    C = compression_curve(ccfg, S)
    PL = permeation_local(ccfg, perm, C)
    RW = unload_states(ccfg, decks)
    if RW is None:
        print("  NOTE: no compression level with an unload -> no compression-side D_c in this comparison "
              "(the holds' D_c is a load measurement and is not used; pass a run with _u<lvl> files as `unload`)")
    q_pred, n_u = flux_from_unload(RW, PL)
    return dict(S=S, C=C, PL=PL, RW=RW, q_pred=q_pred, n_unload=n_u, comp=comp, perm=perm, unload=decks[1:])


# ===========================================================================
#  4. SUMMARY
# ===========================================================================
def print_summary(ccfg, X):
    S, C, PL, RW = X['S'], X['C'], X['PL'], X.get('RW')
    print(f"D_c and kappa vs the solvent volume fraction  (phi_s = {PHI_LABEL[ccfg.PHI_KIND]})")
    print(f"\nCOMPRESSION STATES  {X['comp']['cfg'].RUN_ID}:  M = the level's displacement-strain secant dsigma'/eps_disp (the sweep "
          f"notebook's);  sigma'(phi_s) = PCHIP through the states ({C['n_draws']} MC draws), for the closures.")
    print("  The levels' hold D_c (a LOAD measurement) is not used anywhere in this comparison; the sweep notebook keeps it.")
    print("  level   phi_s            eps     sigma'    M [CI]")
    for i in range(len(S['phi'])):
        if i == 0:
            print(f"  {S['label'][i]:<6s}  {S['phi'][i]:.4f} ± {S['phi_h'][i]:.4f}  0.0000  0.0000    --")
            continue
        print(f"  {S['label'][i]:<6s}  {S['phi'][i]:.4f} ± {S['phi_h'][i]:.4f}  {S['eps'][i]:.4f}  {S['sig'][i]:.4f}    "
              f"{S['M_sec'][i]:.3f} [{S['M_sec_lo'][i]:.3f}, {S['M_sec_hi'][i]:.3f}]")
    sg = PL['seg']
    print(f"\nPERMEATION  {X['perm']['cfg'].RUN_ID}:  steady profile over z in [{PL['z_lo']:.0f}, {PL['z_hi']:.0f}] ({len(PL['z'])} bins, "
          f"{ccfg.PERM_TRIM:g} sigma trimmed at each end);  q = Q/A = {PL['q']:.3e} ± {PL['q_h']:.1e};  dP = {PL['dP']:.4f}")
    print(f"  phi_s: {PL['phi_bot']:.3f} (support side) -> {PL['phi_top']:.3f} (feed side);  zero-flux reference {PL['phi_ref']:.3f}")
    how = {'lag': 'u_F from the per-atom (Lagrangian) steady displacement of the free face, paired with the zero-flux reference',
           'prof': 'u_F from the steady-profile parabola', 'bb': 'u_F from the bounding-box thickness change',
           'trace': 'u_F from the thickness-trace asymptote'}.get(PL['M_perm_how'], str(PL['M_perm_how']))
    print(f"  M_perm = {PL['M_perm']:.4f} [{PL['M_perm_ci'][0]:.4f}, {PL['M_perm_ci'][1]:.4f}] = dP L/(2|u_F|) of the permeation run itself "
          f"({how}; PERM_M_PRIMARY = '{PL['M_perm_how']}') divides every local kappa")
    print("  segment   z      phi_s    dphi_s/dz [CI half]      D_c = q phi_p/(dphi_s/dz) [CI]     kappa = D_c/M_perm [CI]")
    for i in range(len(sg['z'])):
        print(f"  {i + 1:<7d}  {sg['z'][i]:6.1f}  {sg['phi'][i]:.4f}   {sg['slope'][i]:.3e} ± {sg['slope_h'][i]:.1e}   {sg['Dc'][i]:.4f} [{sg['Dc_lo'][i]:.4f}, {sg['Dc_hi'][i]:.4f}]"
              f"          {sg['kappa'][i]:.4f} [{sg['kappa_lo'][i]:.4f}, {sg['kappa_hi'][i]:.4f}]"
              + ('   <-- slope not resolved (CI reaches 0): lower bound only' if not np.isfinite(sg['Dc_hi'][i]) else ''))
    kd = PL['kappa_darcy']
    print("  membrane averages (permeation notebook):  "
          + (f"Darcy kappa = Q L/(A dP) = {kd[0]:.4f} [{kd[1]:.4f}, {kd[2]:.4f}]"
             + (f", kappa_Darcy x M_perm = {kd[0] * PL['M_perm']:.4f};  " if np.isfinite(PL['M_perm']) else ';  ') if kd else '')
          + f"consolidation-fit D_c = {PL['Dc_prof']:.4f} (profile) / {PL['Dc_trace']:.4f} (trace)")
    if RW is not None:
        print("\nUNLOAD (free re-swelling after a hold) -- the compression-side D_c of this comparison (the trusted estimate):")
        print("  run                                         level  eps_disp  phi_s            D_c profile  R^2    D_c Rg trace  M_sec [CI]             kappa = D_c/M_sec [CI]")
        for i in range(len(RW['phi'])):
            print(f"  {RW['run_id'][i]:<43s} {RW['lvl'][i]:<6s} {RW['eps'][i]:.4f}    {RW['phi'][i]:.4f} ± {RW['phi_h'][i]:.4f}  {RW['Dc'][i]:.4f}       {RW['Dc_R2'][i]:.3f}  "
                  f"{RW['Dc_trace_rg'][i]:.4f}        {RW['M_sec'][i]:.3f} [{RW['M_sec_lo'][i]:.3f}, {RW['M_sec_hi'][i]:.3f}]   "
                  f"{RW['kappa'][i]:.4f} [{RW['kappa_lo'][i]:.4f}, {RW['kappa_hi'][i]:.4f}]")
        if kd:
            print(f"  against the permeation Darcy kappa {kd[0]:.4f} [{kd[1]:.4f}, {kd[2]:.4f}]")
    else:
        print("\nUNLOAD: no loaded compression level has an unload -> no compression-side D_c / kappa in this comparison")
    print("\nCLOSURES")
    share = PL['dP'] * PL['L_fit'] / PL['L_mem']
    print(f"  (i)  network stress across the fit range, from the compression curve at the measured phi_s: "
          f"sigma'({PL['phi_bot']:.3f}) - sigma'({PL['phi_top']:.3f}) = {PL['dsig']:.4f}\n"
          f"       against dP = {PL['dP']:.4f} over the whole membrane ({PL['L_mem']:.0f} sigma; {share:.4f} if the pressure fell "
          f"linearly over the {PL['L_fit']:.0f} sigma of the fit range): ratio {PL['dsig'] / share:.2f}")
    print(f"       -> harmonic-mean local kappa over the fit range = q L_fit/dsigma' = {PL['kappa_local_mean']:.4f}"
          + (f"  (Darcy, whole membrane: {kd[0]:.4f})" if kd else ''))
    if X.get('n_unload', 0) > 0:
        how_u = (f"constant at the unload value {RW['Dc'][0]:.4f} (one level, phi_s {RW['phi'][0]:.3f})" if X['n_unload'] == 1
                 else f"interpolated over {X['n_unload']} unload levels")
        print(f"  (ii) flux the UNLOAD D_c predicts for this profile ({how_u}): q = {X['q_pred']:.3e}  vs measured {PL['q']:.3e}  "
              f"(ratio {X['q_pred'] / PL['q']:.2f})")
    else:
        print("  (ii) flux from the compression-side D_c: skipped (no unload level)")


def print_phi_kinds(ccfg, comp, perm, unload=None, kinds=('cal', 'vor', 'mf')):
    """The same comparison under each phi_s definition, one line per permeation segment: how much
    of the result hangs on the choice (and on the noise of the Voronoi profiles).  The compression-side
    D_c it compares with is the UNLOAD one (constant at one level, interpolated over several)."""
    from dataclasses import replace
    print('phi_s definition   segment phi_s    D_c local [CI]              kappa local    |  unload D_c (compression side)   ratio')
    for kind in kinds:
        try:
            X = phi_dependence(replace(ccfg, PHI_KIND=kind), comp, perm, unload)
        except ValueError as e:
            print(f'  {kind:<4s} skipped: {e}')
            continue
        sg, RW = X['PL']['seg'], X.get('RW')
        for i in range(len(sg['z'])):
            if RW is not None and np.isfinite(RW['Dc']).any():
                ok = np.isfinite(RW['Dc'])
                order = np.argsort(RW['phi'][ok])
                Dc_u = float(np.interp(sg['phi'][i], RW['phi'][ok][order], RW['Dc'][ok][order]))
                tail = f"{Dc_u:.4f}{' (one level, constant)' if ok.sum() == 1 else '':<24s} {sg['Dc'][i] / Dc_u:.2f}"
            else:
                tail = 'n/a'
            print(f"  {kind:<16s} {sg['phi'][i]:.3f}            {sg['Dc'][i]:.4f} [{sg['Dc_lo'][i]:.4f}, {sg['Dc_hi'][i]:.4f}]   {sg['kappa'][i]:.4f}         |  {tail}")
        print(f"  {kind:<16s} closures: sigma' drop {X['PL']['dsig']:.4f} (linear share {X['PL']['dP'] * X['PL']['L_fit'] / X['PL']['L_mem']:.4f});  "
              + (f"flux predicted (unload D_c)/measured {X['q_pred'] / X['PL']['q']:.2f}" if X.get('n_unload', 0) else 'flux closure skipped (no unload)'))


# ===========================================================================
#  5. FIGURES -- the permeation membrane on its own (2026-10-08)
# ===========================================================================
def fig_thermo_pressure_norm(ccfg, perm):
    """The permeation notebook's normalised thermodynamic-pressure figure (tri.fig_thermo_pressure_norm:
    solvent, polymer and total traces with the back pressure removed, / dP_ext) with its two dash-dot
    guides dP_th = (4/3)(G/M) dP_ext, one per shear-modulus estimate (cfg.G_REF = the shear decks' G,
    cfg.G_COMP_REF = the lateral network stress of the compression holds; M = this run's).  The deck's
    own library draws and saves it under the deck's plot folder; a copy goes to the comparison folder."""
    cfg, R, P = perm['cfg'], perm['R'], perm['P']
    if cfg.G_REF is None and cfg.G_COMP_REF is None:
        print('NOTE: set G_REF / G_COMP_REF in the permeation Config for the two dP_th guides')
    fig = tri.fig_thermo_pressure_norm(cfg, R, P)
    if fig is not None and ccfg.save:
        ccfg.PLOT_DIR.mkdir(parents=True, exist_ok=True)
        out = ccfg.PLOT_DIR / 'perm_thermo_pressure_norm.png'
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print('copied to', out)
    return fig


def _band(ax, x, prof, color, ls, lw, label, alpha=0.22, scale=1.0, zorder=3):
    """A (mean, lo, hi) bin profile times `scale`: 95 % band + line.  Returns the line handle."""
    m, lo, hi = (np.asarray(a, float) * scale for a in prof)
    ax.fill_between(x, lo, hi, color=color, alpha=alpha, lw=0, zorder=zorder)
    h, = ax.plot(x, m, ls=ls, color=color, lw=lw, zorder=zorder + 1, label=label)
    return h


def _in_window(z, lo, hi):
    """1 inside [lo, hi], nan outside: multiplies a profile to trim it to the window."""
    return np.where((np.asarray(z, float) >= lo) & (np.asarray(z, float) <= hi), 1.0, np.nan)


def fig_solvent_measures(ccfg, perm):
    """The analogue of the permeation notebook's figure 7a for the membrane community: three measures
    of "how much solvent is where", reference (zero flux, dashed) and steady permeation (solid), each
    with its 95 % band over frames / snapshots:

      * solvent CONCENTRATION  c_s = phi_s^cal rho_s,0 / M_s  (left axis, beads per sigma^3; M_s = 1 in LJ
        units, rho_s,0 = the bulk reservoir density at zero flux) -- phi_s^cal is the lambda-calibrated Voronoi
        fraction, the thermodynamic solvent volume fraction;
      * the solvent MASS fraction  rho_s(z)/rho_s,0  (= the NUMBER fraction, all beads of mass 1) -- the local
        density relative to the bulk; it is NOT 1 in the feed reservoir, whose solvent is denser under the feed
        pressure, and it counts solvent packed against the polymer that the calibration assigns to the network;
      * the PROBE POROSITY  eps_g  (lib/psd.py; a solvent-sized probe of radius cfg.PSD_R_PROBE, the fraction of
        space its centre can occupy, as PoreBlazer in Marioni et al.) -- geometric, below the thermodynamic
        fraction because the excluded shell around every bead is counted as solid.

    The fractions sit on the right axis, scaled by rho_s,0 so the three curves are directly comparable (the
    calibrated curve reads phi_s^cal there).  The individual steady Voronoi frames of c_s are drawn faintly
    (ccfg.CONC_FRAMES) to show the frame-to-frame noise.  The z range is the membrane plus ccfg.CONC_MARGIN sigma
    of reservoir on either side.  (b) (ccfg.CONC_DIFF_PANEL): the mass fraction and eps_g minus phi_s^cal."""
    cfg, R, P = perm['cfg'], perm['R'], perm['P']
    rho0 = float(R.get('rho_s0', np.nan))
    if not np.isfinite(rho0):
        print('solvent-measures figure skipped (no reference solvent density -> no rho_s,0)')
        return None
    if P.get('phi_cal') is None:
        print('solvent-measures figure skipped (no calibrated phi_s: VOR_ENABLE + calibration artifact needed)')
        return None
    S, S0 = P.get('psd'), R.get('psd')
    if S is None:
        print('  NOTE: no steady probe porosity (PSD_ENABLE / tessellated frames) -> eps_g left out')
    z = np.asarray(R['z'], float)
    zx = tri.zn(R, z)
    # ---- z window: membrane (reference and steady) + the margin ----
    bw = cfg.binWidth
    zf = P['shade_hi'] if P.get('shade_hi') is not None else R['z_feed']          # the final feed-piston plane
    lo_z = min(R['z_gel_lo'], P['z_mem_lo']) - ccfg.CONC_MARGIN
    hi_z = min(zf + bw, max(R['z_gel_hi'], P['z_mem_hi']) + ccfg.CONC_MARGIN)       # feed side: up to the piston (no vacuum)
    win = _in_window(z, lo_z, hi_z)
    ncol = 2 if ccfg.CONC_DIFF_PANEL else 1
    fig, axes = plt.subplots(1, ncol, figsize=(10.5 * ncol + 2, 7), constrained_layout=True,
                             gridspec_kw=(dict(width_ratios=[1.45, 1]) if ncol == 2 else None), squeeze=False)
    ax = axes[0, 0]
    C_CAL, C_MF, C_POR = WONG['green'], WONG['orange'], WONG['blue']
    rp = cfg.PSD_R_PROBE
    # ---- reference (dashed) ----
    if R.get('phi_cal') is not None:
        _band(ax, zx, [np.asarray(a) * win for a in R['phi_cal']], C_CAL, '--', 1.8, r'$c_s$ reference (zero flux)', scale=rho0, alpha=0.18)
    if R.get('phi_mf') is not None:
        _band(ax, zx, [np.asarray(a) * win for a in R['phi_mf']], C_MF, '--', 1.8, r'$\rho_s/\rho_{s,0}$ reference', scale=rho0, alpha=0.18)
    if S0 is not None:
        _band(ax, zx, [np.asarray(a) * win for a in S0['por']], C_POR, '--', 1.8, r'$\epsilon_g$ reference', scale=rho0, alpha=0.18)
    # ---- steady (solid) ----
    vf = P.get('vf')
    if ccfg.CONC_FRAMES and vf is not None and vf.get('phi_cal') is not None and P.get('vor_steady_ts') is not None:
        sm = np.isin(vf['ts'], P['vor_steady_ts'])
        for row in np.asarray(vf['phi_cal'], float)[sm]:
            ax.plot(zx, row * win * rho0, '-', color=C_CAL, lw=0.8, alpha=0.45, zorder=4)
        n_fr = int(sm.sum())
    else:
        n_fr = len(P.get('vor_steady_ts', []))
    h_cal = _band(ax, zx, [np.asarray(a) * win for a in P['phi_cal']], C_CAL, '-', 2.8,
                  rf'$c_s=\phi_s^{{\rm cal}}\,\rho_{{s,0}}/M_s$  steady' + (f' (faint: each of the {n_fr} frames)' if ccfg.CONC_FRAMES else f' ({n_fr} frames)'), scale=rho0)
    h_mf = h_por = None
    if P.get('phi_mf') is not None:
        n_sn = int((np.asarray(P['mf_ts']) >= P['halt_ts']).sum()) if 'mf_ts' in P else 0
        h_mf = _band(ax, zx, [np.asarray(a) * win for a in P['phi_mf']], C_MF, '-', 2.4,
                     r'mass (= number) fraction $\rho_s(z)/\rho_{s,0}$  steady' + (f' ({n_sn} snapshots)' if n_sn else ''), scale=rho0)
    if S is not None:
        h_por = _band(ax, zx, [np.asarray(a) * win for a in S['por']], C_POR, '-', 2.4,
                      rf'probe porosity $\epsilon_g$  steady ($r_{{\rm probe}}={rp:g}\,\sigma$)', scale=rho0)
    ax.axhline(rho0, color='k', ls=':', lw=1.2, alpha=0.6, zorder=1)
    tri.shade_gel(ax, R, P)
    tri.mark_walls(ax, R, P)
    tri.finish_axes(ax, r'concentration  $c_s=\phi_s\,\rho_{s,0}/M_s$   (beads/$\sigma^3$)',
                    r'(a) concentration $c_s$ vs mass fraction $\rho_s/\rho_{s,0}$ vs probe porosity $\epsilon_g$', R=R)
    ax.set_xlim(sorted([tri.zn(R, lo_z), tri.zn(R, hi_z)]))
    ax.set_ylim(0, rho0 * 1.42)                                           # headroom for the legend and the numbers box
    ax2 = ax.secondary_yaxis('right', functions=(lambda c: c / rho0, lambda f: f * rho0))
    ax2.set_ylabel(r'fraction  $c_s/\rho_{s,0}$  (--)')
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
    # ---- the numbers: interior means and the reservoir values of the mass fraction ----
    def _im(D, k, mask):
        return np.nanmean(np.asarray(D[k][0], float)[mask]) if D.get(k) is not None else np.nan
    im = P['interior']
    txt = ['steady, gel interior:',
           rf'  $c_s$ = {_im(P, "phi_cal", im) * rho0:.3f} $\sigma^{{-3}}$  ($\phi_s^{{\rm cal}}$ = {_im(P, "phi_cal", im):.3f})',
           rf'  $\rho_s/\rho_{{s,0}}$ = {_im(P, "phi_mf", im):.3f}' + (rf'    $\epsilon_g$ = {np.nanmean(S["por"][0][im]):.3f}' if S is not None else '')]
    if P.get('phi_mf') is not None:
        feed = (z >= P['z_mem_hi'] + bw) & (z <= zf - bw)
        perm_res = (z <= P['z_mem_lo'] - bw) & (z >= R['z_perm'] + bw) if np.isfinite(R.get('z_perm', np.nan)) else np.zeros_like(z, bool)
        txt.append(r'reservoirs, $\rho_s/\rho_{s,0}$:  feed ' + (f'{_im(P, "phi_mf", feed):.3f}' if feed.any() else '--')
                   + (rf',  permeate {_im(P, "phi_mf", perm_res):.3f}' if perm_res.any() else ''))
        txt.append(rf'  ($\rho_{{s,0}}$ = {rho0:.3f}: bulk reservoir at zero flux)')
    tri.annotate_box(ax, '\n'.join(txt), loc='upper left', fontsize=10)
    # ---- (b) differences from the calibrated fraction ----
    if ccfg.CONC_DIFF_PANEL:
        bx = axes[0, 1]
        cal = np.asarray(P['phi_cal'][0], float)
        cal_h = 0.5 * (np.asarray(P['phi_cal'][2], float) - np.asarray(P['phi_cal'][1], float))
        mem = P['in_mem'] & np.isfinite(cal) & (cal > cfg.PHI_FLOOR)
        wm = np.where(mem, 1.0, np.nan)
        bx.axhline(0.0, color=C_CAL, lw=2.0, zorder=2, label=r'$\phi_s^{\rm cal}$ (zero line)')
        bx.fill_between(zx, -cal_h * wm, cal_h * wm, color=C_CAL, alpha=0.18, lw=0, zorder=1, label=r'95 % band of $\phi_s^{\rm cal}$')
        for D, k, col, lab in ((P, 'phi_mf', C_MF, r'$\rho_s/\rho_{s,0}-\phi_s^{\rm cal}$  (solvent the calibration assigns to the network)'),
                               (S, 'por', C_POR, r'$\epsilon_g-\phi_s^{\rm cal}$  (excluded shell around the beads, $r_{\rm probe}$)')):
            if D is None or D.get(k) is None:
                continue
            m, lo, hi = (np.asarray(a, float) for a in D[k])
            bx.fill_between(zx, (lo - cal) * wm, (hi - cal) * wm, color=col, alpha=0.2, lw=0, zorder=2)
            bx.plot(zx, (m - cal) * wm, '-', color=col, lw=2.4, zorder=3, label=lab + f'  (interior mean {np.nanmean((m - cal)[P["interior"]]):+.3f})')
        tri.shade_gel(bx, R, P)
        tri.mark_walls(bx, R, P)
        tri.finish_axes(bx, r'difference from $\phi_s^{\rm cal}$  (--)', r'(b) steady membrane: $\rho_s/\rho_{s,0}$ and $\epsilon_g$ minus $\phi_s^{\rm cal}$', R=R)
        bx.set_xlim(sorted([tri.zn(R, P['z_mem_lo'] - 0.5 * ccfg.CONC_MARGIN), tri.zn(R, P['z_mem_hi'] + 0.5 * ccfg.CONC_MARGIN)]))
        smart_legend(bx, fontsize=10)
    fig.suptitle(rf'Three measures of the solvent across the membrane, zero flux (dashed) and steady permeation (solid):  '
                 rf'$c_s=\phi_s^{{\rm cal}}\rho_{{s,0}}/M_s$ ($\lambda$-calibrated Voronoi),  $\rho_s/\rho_{{s,0}}$ (mass = number fraction),  '
                 rf'$\epsilon_g$ (probe porosity)' + '\n'
                 rf'$\rho_{{s,0}}$ = {rho0:.3f} $\sigma^{{-3}}$ (bulk reservoir at zero flux), $M_s$ = 1, $r_{{\rm probe}}$ = {rp:g} $\sigma$  |  {cfg.sim_name}',
                 fontsize=12, fontweight='bold')
    return _save(fig, ccfg, 'perm_solvent_measures')


# ===========================================================================
#  6. FIGURES -- D_c and kappa against phi_s
# ===========================================================================
def _phi_axis(ax, ccfg):
    ax.set_xlabel(rf"solvent volume fraction  $\phi_s$  ({PHI_LABEL[ccfg.PHI_KIND]})")
    ax.invert_xaxis()                                    # compression to the right
    ax.grid(alpha=0.3)


def fig_states(ccfg, X):
    """(a) the steady permeation phi_s(z) with its fits and the compression levels' phi_s;
    (b) the compression curve sigma'(phi_s); (c) the modulus M of both decks, every one from a displacement
    strain (sweep secants, permeation dP L/(2|u_F|), the re-swelling run's secant); the range the membrane spans."""
    S, C, PL = X['S'], X['C'], X['PL']
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(24, 6.5), constrained_layout=True)
    a.fill_between(PL['z'], PL['phi_lo'], PL['phi_hi'], color='0.6', alpha=0.3, lw=0)
    a.plot(PL['z'], PL['phi'], 'o', color='k', ms=4, label='permeation, steady profile')
    a.plot(PL['zz'], PL['phi_fit'], '-', color=WONG['vermillion'], lw=2.5, label=f"polynomial (degree {ccfg.PERM_POLY_DEG})")
    for idx in np.array_split(np.arange(len(PL['z'])), ccfg.PERM_N_SEG)[1:]:
        a.axvline(0.5 * (PL['z'][idx[0]] + PL['z'][idx[0] - 1]), color='0.6', ls=':', lw=1.2)
    a.axhline(PL['phi_ref'], color=WONG['blue'], ls='--', lw=1.8, label=f"zero-flux reference  {sig(PL['phi_ref'])}")
    for i in range(1, len(S['phi'])):
        if PL['phi_bot'] - 0.03 < S['phi'][i] < PL['phi_top'] + 0.03:
            a.axhline(S['phi'][i], color=tri.level_color(i - 1), lw=1.2, alpha=0.8)
            a.annotate(rf"$\varepsilon$ = {S['label'][i]}", (PL['z'][-1], S['phi'][i]), ha='right', va='bottom', fontsize=10, color=tri.level_color(i - 1))
    a.set(xlabel=r'$z$  ($\sigma$)   [support $\leftarrow$ $\rightarrow$ feed]', ylabel=r'$\phi_s$')
    a.set_title(r'(a) steady membrane: $\phi_s(z)$, segments, compression levels', fontsize=14)
    a.grid(alpha=0.3)
    smart_legend(a, fontsize=10)
    pp = np.linspace(S['phi'][-1], S['phi'][0], 300)
    lo, hi = C['sig_ci'](pp)
    b.fill_between(pp, lo, hi, color=WONG['blue'], alpha=0.2, lw=0)
    b.plot(pp, C['sig'](pp), '-', color=WONG['blue'], lw=2.2, label='PCHIP through the levels')
    b.plot(S['phi'], S['sig'], 'o', color=WONG['blue'], ms=9, label='compression levels')
    b.errorbar(S['phi'], S['sig'], xerr=S['phi_h'], yerr=S['sig_h'], fmt='none', color=WONG['blue'], capsize=3)
    b.plot([PL['phi_bot']], [PL['dP']], '*', color=WONG['vermillion'], ms=18, zorder=5,
           label=rf"permeation: $\Delta P$ = {sig(PL['dP'])} at the support-side $\phi_s$")
    b.set_ylim(-0.02, max(4 * PL['dP'], 0.3))
    b.set_ylabel(r"network stress  $\sigma'_{zz}$")
    b.set_title(r"(b) compression: $\sigma'(\phi_s)$", fontsize=14)
    # (c) every modulus from a displacement strain; the volume-fraction tangent M_cur was dropped on 2026-10-15
    c.errorbar(S['phi'][1:], S['M_sec'][1:], yerr=[S['M_sec'][1:] - S['M_sec_lo'][1:], S['M_sec_hi'][1:] - S['M_sec'][1:]], fmt='s',
               color=WONG['blue'], ms=8, capsize=4, label=r"compression levels: secant $\Delta\sigma'/\varepsilon_{\rm disp}$")
    if np.isfinite(PL['M_perm']):                      # the permeation M is displacement-derived too: dP L/(2|u_F|)
        c.fill_between([PL['phi_bot'], PL['phi_top']], PL['M_perm_ci'][0], PL['M_perm_ci'][1], color=WONG['vermillion'], alpha=0.2, lw=0)
        _bands(c, PL, PL['M_perm'], WONG['vermillion'], rf"permeation: $\Delta P\,L/(2|u_F|)$ = {sig(PL['M_perm'])}")
    RW = X.get('RW')
    if RW is not None:
        for i in range(len(RW['phi'])):
            c.errorbar([RW['phi'][i]], [RW['M_sec'][i]], yerr=[[RW['M_sec'][i] - RW['M_sec_lo'][i]], [RW['M_sec_hi'][i] - RW['M_sec'][i]]], fmt='*',
                       color=WONG['green'], ms=16, capsize=4, zorder=6,
                       label=f"unload level {RW['lvl'][i]} ({RW['run_id'][i]}): secant $\\Delta\\sigma'/\\varepsilon_{{\\rm disp}}$ = {sig(RW['M_sec'][i])}")
    c.set_ylabel(r"modulus  $M$")
    c.set_ylim(bottom=0)
    c.set_title('(c) modulus $M$ of both decks, each from its displacement strain', fontsize=14)
    for ax in (b, c):
        ax.axvspan(PL['phi_bot'], PL['phi_top'], color=WONG['orange'], alpha=0.15, label=r'$\phi_s$ range of the membrane')
        _phi_axis(ax, ccfg)
    smart_legend(b, fontsize=10)
    c.legend(loc='upper left', fontsize=9)
    fig.suptitle(f"States  |  compression {X['comp']['cfg'].RUN_ID}  |  permeation {X['perm']['cfg'].RUN_ID}", fontsize=12, fontweight='bold')
    return _save(fig, ccfg, 'states_phi')


def _seg_points(ax, ccfg, sg, key):
    """The z-segment values with their asymmetric CI; a segment whose slope CI reaches 0 keeps
    its lower bar and is drawn hollow (its upper bound is unbounded)."""
    y, lo, hi = sg[key], sg[key + '_lo'], sg[key + '_hi']
    res = np.isfinite(hi)
    lab = f"permeation, local: {ccfg.PERM_N_SEG} z-segments (line fits)"
    if res.any():
        ax.errorbar(sg['phi'][res], y[res], yerr=[y[res] - lo[res], hi[res] - y[res]], fmt='D', color=WONG['vermillion'], ms=9, capsize=5, lw=2, label=lab)
    if (~res).any():
        ax.errorbar(sg['phi'][~res], y[~res], yerr=[y[~res] - lo[~res], 0 * y[~res]], fmt='D', mfc='none', color=WONG['vermillion'], ms=9, capsize=5, lw=2,
                    label='... slope not resolved: lower bound only')


def _bands(ax, PL, y, col, lab):
    ax.plot([PL['phi_bot'], PL['phi_top']], [y, y], '-', color=col, lw=3, alpha=0.8, solid_capstyle='butt', label=lab)


def fig_Dc(ccfg, X):
    """D_c against phi_s: the UNLOAD levels (free re-swelling profile fit, stars; Rg thickness trace hollow) -- the only
    compression-side D_c here (2026-10-21) -- the permeation steady state (local: q phi_p / (dphi_s/dz), segments + the
    polynomial curve) and the permeation notebook's membrane averages as bars over the membrane's phi_s range."""
    PL, RW = X['PL'], X.get('RW')
    sg = PL['seg']
    fig, ax = plt.subplots(figsize=(11, 7), constrained_layout=True)
    ax.plot(PL['phi_fit'], PL['Dc_fit'], '-', color=WONG['vermillion'], lw=2.2, alpha=0.8,
            label=rf"permeation, local: $q\,\phi_p/(\mathrm{{d}}\phi_s/\mathrm{{d}}z)$, degree-{ccfg.PERM_POLY_DEG} profile")
    _seg_points(ax, ccfg, sg, 'Dc')
    if np.isfinite(PL['Dc_prof']):
        _bands(ax, PL, PL['Dc_prof'], WONG['green'], f"permeation, transient fit (profile): {sig(PL['Dc_prof'])}")
    if np.isfinite(PL['Dc_trace']):
        _bands(ax, PL, PL['Dc_trace'], WONG['reddishpurple'], f"permeation, transient fit (trace): {sig(PL['Dc_trace'])}")
    if PL['kappa_darcy'] and np.isfinite(PL['M_perm']):
        _bands(ax, PL, PL['kappa_darcy'][0] * PL['M_perm'], '0.35', rf"permeation, Darcy $\kappa \times M$ = {sig(PL['kappa_darcy'][0] * PL['M_perm'])}")
    if RW is not None:
        for i in range(len(RW['phi'])):
            ax.errorbar([RW['phi'][i]], [RW['Dc'][i]], xerr=[RW['phi_h'][i]], fmt='*', color=WONG['blue'], ms=18, capsize=4, zorder=6,
                        label=f"UNLOAD, free re-swelling profile fit: {sig(RW['Dc'][i])}  [{RW['run_id'][i]}, level {RW['lvl'][i]}]")
            if np.isfinite(RW['Dc_trace_rg'][i]):
                ax.plot([RW['phi'][i]], [RW['Dc_trace_rg'][i]], '*', color=WONG['blue'], mfc='none', ms=15, mew=1.5, zorder=6,
                        label=f"... its Rg thickness trace: {sig(RW['Dc_trace_rg'][i])}")
    else:
        ax.text(0.02, 0.03, 'no compression level with an unload: no compression-side $D_c$ drawn', transform=ax.transAxes, fontsize=11, color='0.4')
    ax.set_ylabel(r'$D_c$  ($\sigma^2/\tau$)')
    ax.set_yscale('log')
    ax.set_title(r'Cooperative diffusivity $D_c(\phi_s)$: unload (free re-swelling) vs steady permeation  --  no hold (load) $D_c$', fontsize=15)
    _phi_axis(ax, ccfg)
    smart_legend(ax, fontsize=10)
    return _save(fig, ccfg, 'Dc_vs_phi')


def fig_kappa(ccfg, X):
    """kappa against phi_s: the UNLOAD levels (free re-swelling D_c over the level's displacement-strain secant M,
    stars) -- the only compression-side kappa here (2026-10-21) -- the permeation local kappa = D_c(phi_s)/M_perm and
    the Darcy membrane average."""
    PL, RW = X['PL'], X.get('RW')
    sg = PL['seg']
    fig, ax = plt.subplots(figsize=(11, 7), constrained_layout=True)
    ax.plot(PL['phi_fit'], PL['kappa_fit'], '-', color=WONG['vermillion'], lw=2.2, alpha=0.8,
            label=rf"permeation, local: $q\,\phi_p/(\mathrm{{d}}\phi_s/\mathrm{{d}}z)\;/\;M_{{\rm perm}}$, degree-{ccfg.PERM_POLY_DEG} profile;  "
                  rf"$M_{{\rm perm}} = \Delta P\,L/(2|u_F|)$ = {sig(PL['M_perm'])}")
    _seg_points(ax, ccfg, sg, 'kappa')
    kd = PL['kappa_darcy']
    if kd:
        ax.fill_between([PL['phi_bot'], PL['phi_top']], kd[1], kd[2], color='0.35', alpha=0.2, lw=0)
        _bands(ax, PL, kd[0], '0.35', rf"permeation, Darcy $QL/(A\,\Delta P)$ = {sig(kd[0])} (membrane average)")
    if RW is not None:
        for i in range(len(RW['phi'])):
            ax.errorbar([RW['phi'][i]], [RW['kappa'][i]], yerr=[[RW['kappa'][i] - RW['kappa_lo'][i]], [RW['kappa_hi'][i] - RW['kappa'][i]]], fmt='*',
                        color=WONG['blue'], ms=18, capsize=4, zorder=6,
                        label=rf"UNLOAD $D_c$ / the level's secant $M$ = {sig(RW['kappa'][i])}  [{RW['run_id'][i]}, level {RW['lvl'][i]}]")
    else:
        ax.text(0.02, 0.03, 'no compression level with an unload: no compression-side $\kappa$ drawn', transform=ax.transAxes, fontsize=11, color='0.4')
    ax.set_ylabel(r'$\kappa$  ($\sigma^5/(\epsilon\,\tau)$)')
    ax.set_yscale('log')
    ax.set_title(r'Hydraulic permeability $\kappa(\phi_s)$: unload (free re-swelling) vs steady permeation  --  no hold (load) $D_c$', fontsize=15)
    _phi_axis(ax, ccfg)
    smart_legend(ax, fontsize=10)
    return _save(fig, ccfg, 'kappa_vs_phi')


def fig_closure(ccfg, X):
    """(a) the network stress the compression curve assigns to the measured phi_s(z), against the
    linear pressure drop; (b) measured flux against the one the compression D_c(phi_s) predicts."""
    S, C, PL = X['S'], X['C'], X['PL']
    P = X['perm']['P']
    fig, (a, b) = plt.subplots(1, 2, figsize=(17, 6.5), constrained_layout=True, gridspec_kw=dict(width_ratios=[2, 1]))
    lo, hi = C['sig_ci'](PL['phi'])
    a.fill_between(PL['z'], lo, hi, color=WONG['blue'], alpha=0.2, lw=0)
    # NOT a measurement of the membrane's network stress (2026-10-15): the pore pressure under flow is not measured yet
    # (it needs the solvent chemical-potential profile), so sigma'(z) here is the value the compression EQUILIBRIUM curve
    # assigns to the measured phi_s(z), and the dashed line is the pore-pressure drop a uniform kappa WOULD give -- an
    # assumption drawn for scale, to be replaced by the measured profile when the permeation notebook has it.
    a.plot(PL['z'], C['sig'](PL['phi']), 'o', color=WONG['blue'], ms=4, label=r"$\sigma'_\mathrm{compression}(\phi_s(z))$: the equilibrium curve's value at the measured $\phi_s$, bin by bin (not measured)")
    a.plot(PL['zz'], PL['sig_z'], '-', color=WONG['blue'], lw=2.2, label='the same on the fitted profile')
    zl = np.array([P['z_mem_lo'], P['z_mem_hi']])
    a.plot(zl, [PL['dP'], 0.0], 'k--', lw=1.8, label=rf"ASSUMED linear pore-pressure drop (uniform $\kappa$): $\Delta P$ = {sig(PL['dP'])} at the support, 0 at the feed face")
    a.set(xlabel=r'$z$  ($\sigma$)   [support $\leftarrow$ $\rightarrow$ feed]', ylabel=r"network stress  $\sigma'_{zz}$  (inferred)")
    a.set_title(r"(a) network stress IMPLIED by the compression curve at the measured $\phi_s(z)$ vs an assumed linear $p_{\rm pore}$ drop", fontsize=13)
    a.grid(alpha=0.3)
    smart_legend(a, fontsize=10)
    b.bar([0], [PL['q']], yerr=[PL['q_h']], color=WONG['vermillion'], capsize=8, width=0.6)
    if np.isfinite(X['q_pred']):
        b.bar([1], [X['q_pred']], color=WONG['blue'], width=0.6)
        b.set_title(f"(b) flux: predicted (unload $D_c$, {X['n_unload']} level{'s' if X['n_unload'] > 1 else ''})/measured = {X['q_pred'] / PL['q']:.2f}", fontsize=14)
    else:
        b.text(1, 0.5 * PL['q'], 'no unload level:\nno prediction', ha='center', fontsize=11, color='0.4')
        b.set_title('(b) flux: no compression-side $D_c$ (no unload level)', fontsize=14)
    b.set_xticks([0, 1])
    b.set_xticklabels(['measured\n$Q/A$', 'from the UNLOAD $D_c$\n' + r'$\frac{1}{L}\int D_c\,\mathrm{d}\phi_s/\phi_p$'], fontsize=12)
    b.set_ylabel(r'flux  $q$  ($\sigma/\tau$)')
    b.grid(axis='y', alpha=0.3)
    return _save(fig, ccfg, 'closure')
