"""Cross-deck comparison -- figures built from MORE THAN ONE deck's run (2026-10-03).

Notebook: scripts/cross_deck_comparison.ipynb.  Each deck is loaded with its own library
(lib/triaxial.py for the two-piston compression and permeation decks, lib/shear.py for
shear_slab) through `load_deck`, which returns a plain dict
    dict(kind='compression' | 'permeation' | 'shear', cfg, R, levels | P)
so a comparison function takes the decks it needs and nothing here re-implements a loader.

First comparison: D_c and kappa as functions of the solvent volume fraction phi_s.

  Compression sweep (uniform states).  Level i is a uniformly compressed gel at phi_s,i with the
  network stress sigma'_i and the hold-relaxation D_c,i (tri.fit_Dc).  The levels give the
  equilibrium curve sigma'(phi_s), and with it the modulus that belongs to a small disturbance
  of the state phi_s IN ITS CURRENT CONFIGURATION (the frame fit_Dc works in: it uses the
  compressed thickness).  With the polymer conserved and the lateral box fixed, phi_p L = const,
  so an increment of compressive strain is  de = -dL/L = dphi_p/phi_p  and

      M_cur(phi_s) = phi_p dsigma'/dphi_p ,        kappa(phi_s) = D_c / M_cur .

  (The sweep notebook's kappa = D_c/M divides by the SECANT modulus from eps = 0; the two agree
  at small strain and part as the gel stiffens.)

  Permeation (one run, a range of phi_s at once).  In the steady state the network is at rest and
  the solvent flux q = Q/A is the same at every z.  Darcy, the uniform total stress
  (dp/dz = dsigma'/dz) and the same sigma'(phi_s) give  q = kappa |dsigma'/dphi_s| |dphi_s/dz|, i.e.

      D_c(phi_s) = kappa M_cur = q phi_p / |dphi_s/dz| ,   kappa(phi_s) = D_c(phi_s) / M_cur(phi_s) .

  D_c needs the flux and the slope of the steady phi_s(z) profile only -- no modulus, no fit of a
  transient; kappa borrows M_cur from the compression curve.  The membrane-average numbers of the
  permeation notebook (Darcy kappa = Q L/(A dP), the consolidation-fit D_c) are drawn as bands
  over the phi_s range the membrane spans.

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
    N_MC: int = 400                       # Monte-Carlo draws for the CI of the compression curve sigma'(phi_s) and M_cur
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
        tri.add_perm_displacement(cfg, R, P, verbose=verbose)
        return dict(kind='permeation', cfg=cfg, R=R, P=P)
    levels = [L for L in (tri.load_level(cfg, R, lvl, verbose=verbose) for lvl in cfg.COMP_LEVELS) if L is not None]
    levels.sort(key=lambda L: L['eps'])
    tri.add_volume_fractions(cfg, R, levels)
    return dict(kind='compression', cfg=cfg, R=R, levels=levels)


# ===========================================================================
#  2. COMPRESSION: uniform states -> sigma'(phi_s), M_cur(phi_s), D_c(phi_s), kappa(phi_s)
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
    """The sweep as a table of uniform states: the reference (sigma' = 0, no hold -> no D_c) and
    one row per level.  -> dict of arrays (phi, phi_h, sig, sig_h, eps, Dc, M_sec, R2, label)
    ordered by increasing compression, row 0 = the reference."""
    key = 'phi_' + ccfg.PHI_KIND
    R = comp['R']
    if R.get(key) is None:
        raise ValueError(f"the compression reference has no {key} profile (VOR_ENABLE / calibration artifact?)")
    p0, h0 = _interior_mean(R[key], R['interior'], ccfg.ci_level)
    rows = [dict(label='ref', phi=p0, phi_h=h0, sig=0.0, sig_h=0.0, eps=0.0, Dc=np.nan, M_sec=np.nan, R2=np.nan)]
    for L in comp['levels']:
        if L.get(key) is None:
            print(f"  NOTE: level {L['lvl']} has no {key} profile -> left out of the comparison")
            continue
        p, h = _interior_mean(L[key], L['interior'], ccfg.ci_level)
        F = L.get('Dc')
        rows.append(dict(label=str(L['lvl']), phi=p, phi_h=h, sig=float(L['dsig_net']),
                         sig_h=float((L['M_net_hi'] - L['M_net_lo']) / 2 * L['eps_M']), eps=float(L['eps_M']),
                         Dc=float(F['Dc']) if F else np.nan, M_sec=float(L['M_net']), R2=float(F['R2']) if F else np.nan))
    S = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    if not np.all(np.diff(S['phi']) < 0):
        raise ValueError('phi_s does not fall monotonically along the sweep: ' + ', '.join(f'{p:.3f}' for p in S['phi']))
    return S


def _curve(phi_s, sg):
    """Monotone (PCHIP) sigma'(phi_p) through the states -> (sigma'(phi_s), M_cur(phi_s)) callables."""
    f = PchipInterpolator(1.0 - np.asarray(phi_s, float), np.asarray(sg, float), extrapolate=True)
    df = f.derivative()
    return (lambda ps: f(1.0 - np.asarray(ps, float))), (lambda ps: (1.0 - np.asarray(ps, float)) * df(1.0 - np.asarray(ps, float)))


def compression_curve(ccfg, S, seed=0):
    """sigma'(phi_s) and M_cur(phi_s) = phi_p dsigma'/dphi_p from the states S, with a Monte-Carlo
    CI (the levels' sigma' and phi_s drawn from their CIs, the curve rebuilt each time).
    -> dict(sig(phi), M(phi), sig_ci(phi) -> (lo, hi), M_ci(phi) -> (lo, hi), M_lvl, M_lvl_lo, M_lvl_hi,
            M_strain (the same modulus from the strain: (1 - eps) dsigma'/deps), kappa, kappa_lo, kappa_hi)."""
    f_sig, f_M = _curve(S['phi'], S['sig'])
    rng = np.random.default_rng(seed)
    z = stats.norm.ppf(0.5 + ccfg.ci_level / 2)
    draws = []
    for _ in range(int(ccfg.N_MC)):
        ph = S['phi'] + rng.normal(0.0, 1.0, len(S['phi'])) * S['phi_h'] / z
        sg = S['sig'] + rng.normal(0.0, 1.0, len(S['sig'])) * S['sig_h'] / z
        if np.all(np.diff(ph) < 0):
            draws.append(_curve(ph, np.maximum.accumulate(sg)))          # sigma' is non-decreasing with compression
    q = [100 * (0.5 - ccfg.ci_level / 2), 100 * (0.5 + ccfg.ci_level / 2)]

    def _ci(which):
        def g(ps):
            if not draws:
                v = (f_sig, f_M)[which](ps)
                return v, v
            A = np.array([d[which](ps) for d in draws])
            return tuple(np.percentile(A, q, axis=0))
        return g

    C = dict(sig=f_sig, M=f_M, sig_ci=_ci(0), M_ci=_ci(1), n_draws=len(draws))
    C['M_lvl'] = f_M(S['phi'])
    C['M_lvl_lo'], C['M_lvl_hi'] = C['M_ci'](S['phi'])
    C['M_strain'] = (1.0 - S['eps']) * np.gradient(S['sig'], S['eps'])
    with np.errstate(invalid='ignore', divide='ignore'):
        C['kappa'] = S['Dc'] / C['M_lvl']
        C['kappa_lo'], C['kappa_hi'] = S['Dc'] / C['M_lvl_hi'], S['Dc'] / C['M_lvl_lo']
        C['kappa_sec'] = S['Dc'] / S['M_sec']
    return C


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
       numbers of the permeation notebook (kappa_darcy, Dc_prof, Dc_trace, M_perm) and the closures."""
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
    # ---- continuous: polynomial phi_s(z) ----
    poly = np.polynomial.Polynomial.fit(zs, ps, ccfg.PERM_POLY_DEG)
    zz = np.linspace(zs[0], zs[-1], 200)
    pf, dpf = poly(zz), poly.deriv()(zz)
    with np.errstate(divide='ignore', invalid='ignore'):
        Dc_fit = np.where(dpf > 0, q * (1.0 - pf) / dpf, np.nan)
    out.update(poly=poly, zz=zz, phi_fit=pf, Dc_fit=Dc_fit, kappa_fit=Dc_fit / C['M'](pf),
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
        Ml, Mh = C['M_ci'](pm)
        for k, v in (('z', float(np.mean(zs[idx]))), ('phi', pm), ('slope', r.slope), ('slope_h', sh), ('Dc', D), ('Dc_lo', D_lo), ('Dc_hi', D_hi),
                     ('kappa', D / float(C['M'](pm))), ('kappa_lo', D_lo / float(Mh)), ('kappa_hi', D_hi / float(Ml))):
            seg[k].append(float(v))
    out['seg'] = {k: np.array(v) for k, v in seg.items()}
    # ---- the permeation notebook's membrane averages ----
    kd = (fl.get('k') or {}).get('N_measured') or (fl.get('k') or {}).get('N_applied')
    out['kappa_darcy'] = (kd['k'], kd['lo'], kd['hi']) if kd else None
    F = P.get('Dc')
    out['Dc_prof'] = float(F['Dc']) if F else np.nan
    out['Dc_trace'] = float(F['trace']['Dc']) if (F and F.get('trace')) else np.nan
    Mp = P.get('M_perm')
    out['M_perm'] = float(Mp[Mp['primary']]['M']) if (Mp and Mp.get('primary')) else np.nan
    # ---- closures ----
    #  (i) the network stress across the fit range, read off the compression curve, against the share of dP
    #      that falls inside the fit range if the pore pressure were linear in z
    out['sig_z'] = C['sig'](pf)                                        # sigma'(z) the compression curve assigns to phi_s(z)
    out['dsig'] = float(out['sig_z'][0] - out['sig_z'][-1])
    L_mem = float(P['z_mem_hi'] - P['z_mem_lo'])
    out['L_mem'], out['L_fit'] = L_mem, float(zs[-1] - zs[0])
    out['kappa_local_mean'] = float(q * out['L_fit'] / out['dsig']) if out['dsig'] > 0 else np.nan   # harmonic mean of kappa(z)
    return out


def flux_from_compression(S, C, PL):
    """The flux the COMPRESSION D_c(phi_s) predicts for the measured steady profile:
    q L = integral of D_c/(1 - phi_s) dphi_s over the phi_s range of the fit.  D_c is interpolated
    linearly between the levels and held at the first level's value above it (no hold at eps = 0).
    -> (q_pred, fraction of the phi_s range that lies above the first level)."""
    ok = np.isfinite(S['Dc'])
    ph, Dc = S['phi'][ok][::-1], S['Dc'][ok][::-1]                     # increasing phi_s
    x = np.linspace(PL['phi_bot'], PL['phi_top'], 400)
    D = np.interp(x, ph, Dc)
    integ = float(np.sum(0.5 * (D[1:] / (1 - x[1:]) + D[:-1] / (1 - x[:-1])) * np.diff(x)))
    above = float(np.clip((PL['phi_top'] - ph[-1]) / (PL['phi_top'] - PL['phi_bot']), 0.0, 1.0))
    return integ / PL['L_fit'], above


def phi_dependence(ccfg, comp, perm):
    """Everything for the D_c(phi_s) / kappa(phi_s) comparison: -> dict(S, C, PL, q_pred, q_above)."""
    S = compression_states(ccfg, comp)
    C = compression_curve(ccfg, S)
    PL = permeation_local(ccfg, perm, C)
    q_pred, above = flux_from_compression(S, C, PL)
    return dict(S=S, C=C, PL=PL, q_pred=q_pred, q_above=above, comp=comp, perm=perm)


# ===========================================================================
#  4. SUMMARY
# ===========================================================================
def print_summary(ccfg, X):
    S, C, PL = X['S'], X['C'], X['PL']
    print(f"D_c and kappa vs the solvent volume fraction  (phi_s = {PHI_LABEL[ccfg.PHI_KIND]})")
    print(f"\nCOMPRESSION  {X['comp']['cfg'].RUN_ID}:  uniform states; M_cur = phi_p dsigma'/dphi_p (PCHIP through the states, "
          f"{C['n_draws']} MC draws); kappa = D_c/M_cur")
    print("  level   phi_s            eps     sigma'    M_secant   M_cur [CI]              M_cur(strain)   D_c      R^2    kappa [CI]                 kappa(secant)")
    for i in range(len(S['phi'])):
        if i == 0:
            print(f"  {S['label'][i]:<6s}  {S['phi'][i]:.4f} ± {S['phi_h'][i]:.4f}  0.0000  0.0000    --         {C['M_lvl'][i]:.3f} [{C['M_lvl_lo'][i]:.3f}, {C['M_lvl_hi'][i]:.3f}]"
                  f"   {C['M_strain'][i]:.3f}           --")
            continue
        print(f"  {S['label'][i]:<6s}  {S['phi'][i]:.4f} ± {S['phi_h'][i]:.4f}  {S['eps'][i]:.4f}  {S['sig'][i]:.4f}    {S['M_sec'][i]:.3f}      "
              f"{C['M_lvl'][i]:.3f} [{C['M_lvl_lo'][i]:.3f}, {C['M_lvl_hi'][i]:.3f}]   {C['M_strain'][i]:.3f}           {S['Dc'][i]:.4f}   {S['R2'][i]:.3f}  "
              f"{C['kappa'][i]:.4f} [{C['kappa_lo'][i]:.4f}, {C['kappa_hi'][i]:.4f}]   {C['kappa_sec'][i]:.4f}")
    sg = PL['seg']
    print(f"\nPERMEATION  {X['perm']['cfg'].RUN_ID}:  steady profile over z in [{PL['z_lo']:.0f}, {PL['z_hi']:.0f}] ({len(PL['z'])} bins, "
          f"{ccfg.PERM_TRIM:g} sigma trimmed at each end);  q = Q/A = {PL['q']:.3e} ± {PL['q_h']:.1e};  dP = {PL['dP']:.4f}")
    print(f"  phi_s: {PL['phi_bot']:.3f} (support side) -> {PL['phi_top']:.3f} (feed side);  zero-flux reference {PL['phi_ref']:.3f}")
    print("  segment   z      phi_s    dphi_s/dz [CI half]      D_c = q phi_p/(dphi_s/dz) [CI]     kappa = D_c/M_cur [CI]")
    for i in range(len(sg['z'])):
        print(f"  {i + 1:<7d}  {sg['z'][i]:6.1f}  {sg['phi'][i]:.4f}   {sg['slope'][i]:.3e} ± {sg['slope_h'][i]:.1e}   {sg['Dc'][i]:.4f} [{sg['Dc_lo'][i]:.4f}, {sg['Dc_hi'][i]:.4f}]"
              f"          {sg['kappa'][i]:.4f} [{sg['kappa_lo'][i]:.4f}, {sg['kappa_hi'][i]:.4f}]"
              + ('   <-- slope not resolved (CI reaches 0): lower bound only' if not np.isfinite(sg['Dc_hi'][i]) else ''))
    kd = PL['kappa_darcy']
    print("  membrane averages (permeation notebook):  "
          + (f"Darcy kappa = Q L/(A dP) = {kd[0]:.4f} [{kd[1]:.4f}, {kd[2]:.4f}];  " if kd else '')
          + f"consolidation-fit D_c = {PL['Dc_prof']:.4f} (profile) / {PL['Dc_trace']:.4f} (trace);  M = {PL['M_perm']:.4f}")
    print("\nCLOSURES")
    share = PL['dP'] * PL['L_fit'] / PL['L_mem']
    print(f"  (i)  network stress across the fit range, from the compression curve at the measured phi_s: "
          f"sigma'({PL['phi_bot']:.3f}) - sigma'({PL['phi_top']:.3f}) = {PL['dsig']:.4f}\n"
          f"       against dP = {PL['dP']:.4f} over the whole membrane ({PL['L_mem']:.0f} sigma; {share:.4f} if the pressure fell "
          f"linearly over the {PL['L_fit']:.0f} sigma of the fit range): ratio {PL['dsig'] / share:.2f}")
    print(f"       -> harmonic-mean local kappa over the fit range = q L_fit/dsigma' = {PL['kappa_local_mean']:.4f}"
          + (f"  (Darcy, whole membrane: {kd[0]:.4f})" if kd else ''))
    print(f"  (ii) flux the compression D_c(phi_s) predicts for this profile: q = {X['q_pred']:.3e}  vs measured {PL['q']:.3e}  "
          f"(ratio {X['q_pred'] / PL['q']:.2f}; {X['q_above']:.0%} of the phi_s range lies above the first level, where D_c is held at its value)")


def print_phi_kinds(ccfg, comp, perm, kinds=('cal', 'vor', 'mf')):
    """The same comparison under each phi_s definition, one line per permeation segment: how much
    of the result hangs on the choice (and on the noise of the Voronoi profiles)."""
    from dataclasses import replace
    print('phi_s definition   segment phi_s    D_c local [CI]              kappa local    |  compression D_c at that phi_s   ratio')
    for kind in kinds:
        try:
            X = phi_dependence(replace(ccfg, PHI_KIND=kind), comp, perm)
        except ValueError as e:
            print(f'  {kind:<4s} skipped: {e}')
            continue
        S, sg = X['S'], X['PL']['seg']
        ok = np.isfinite(S['Dc'])
        for i in range(len(sg['z'])):
            Dc_c = float(np.interp(sg['phi'][i], S['phi'][ok][::-1], S['Dc'][ok][::-1]))
            print(f"  {kind:<16s} {sg['phi'][i]:.3f}            {sg['Dc'][i]:.4f} [{sg['Dc_lo'][i]:.4f}, {sg['Dc_hi'][i]:.4f}]   {sg['kappa'][i]:.4f}         |  "
                  f"{Dc_c:.4f}{' (first level, held)' if sg['phi'][i] > S['phi'][ok].max() else '':<22s} {sg['Dc'][i] / Dc_c:.2f}")
        print(f"  {kind:<16s} closures: sigma' drop {X['PL']['dsig']:.4f} (linear share {X['PL']['dP'] * X['PL']['L_fit'] / X['PL']['L_mem']:.4f});  "
              f"flux predicted/measured {X['q_pred'] / X['PL']['q']:.2f}")


# ===========================================================================
#  5. FIGURES
# ===========================================================================
def _phi_axis(ax, ccfg):
    ax.set_xlabel(rf"solvent volume fraction  $\phi_s$  ({PHI_LABEL[ccfg.PHI_KIND]})")
    ax.invert_xaxis()                                    # compression to the right
    ax.grid(alpha=0.3)


def fig_states(ccfg, X):
    """(a) the steady permeation phi_s(z) with its fits and the compression levels' phi_s;
    (b) the compression curve sigma'(phi_s) and (c) M_cur(phi_s), with the range the membrane spans."""
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
    for ax, f, ci, y, yl, tt in ((b, C['sig'], C['sig_ci'], S['sig'], r"network stress  $\sigma'_{zz}$", r"(b) compression: $\sigma'(\phi_s)$"),
                                 (c, C['M'], C['M_ci'], C['M_lvl'], r"$M_\mathrm{cur}=\phi_p\,\mathrm{d}\sigma'/\mathrm{d}\phi_p$", r'(c) current-frame tangent modulus')):
        lo, hi = ci(pp)
        ax.fill_between(pp, lo, hi, color=WONG['blue'], alpha=0.2, lw=0)
        ax.plot(pp, f(pp), '-', color=WONG['blue'], lw=2.2, label='PCHIP through the levels')
        ax.plot(S['phi'], y, 'o', color=WONG['blue'], ms=9, label='compression levels')
        ax.axvspan(PL['phi_bot'], PL['phi_top'], color=WONG['orange'], alpha=0.15, label=r'$\phi_s$ range of the membrane')
        ax.set_ylabel(yl)
        ax.set_yscale('log' if ax is c else 'linear')
        ax.set_title(tt, fontsize=14)
        _phi_axis(ax, ccfg)
    b.errorbar(S['phi'], S['sig'], xerr=S['phi_h'], yerr=S['sig_h'], fmt='none', color=WONG['blue'], capsize=3)
    b.plot([PL['phi_bot']], [PL['dP']], '*', color=WONG['vermillion'], ms=18, zorder=5,
           label=rf"permeation: $\Delta P$ = {sig(PL['dP'])} at the support-side $\phi_s$")
    b.set_ylim(-0.02, max(4 * PL['dP'], 0.3))
    c.plot(S['phi'][1:], S['M_sec'][1:], 's', color='0.45', ms=7, label=r"secant $M=\sigma'/\varepsilon$ (sweep notebook)")
    for ax in (b, c):
        smart_legend(ax, fontsize=10)
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
    """D_c against phi_s: the compression holds (one point per level), the permeation steady state
    (local: q phi_p / (dphi_s/dz), segments + the polynomial curve) and the permeation notebook's
    membrane averages as bars over the membrane's phi_s range."""
    S, C, PL = X['S'], X['C'], X['PL']
    sg = PL['seg']
    fig, ax = plt.subplots(figsize=(11, 7), constrained_layout=True)
    ok = np.isfinite(S['Dc'])
    ax.plot(S['phi'][ok], S['Dc'][ok], 'o-', color=WONG['blue'], ms=10, lw=2, label='compression: hold relaxation fit, one level each')
    ax.plot(PL['phi_fit'], PL['Dc_fit'], '-', color=WONG['vermillion'], lw=2.2, alpha=0.8,
            label=rf"permeation, local: $q\,\phi_p/(\mathrm{{d}}\phi_s/\mathrm{{d}}z)$, degree-{ccfg.PERM_POLY_DEG} profile")
    _seg_points(ax, ccfg, sg, 'Dc')
    if np.isfinite(PL['Dc_prof']):
        _bands(ax, PL, PL['Dc_prof'], WONG['green'], f"permeation, transient fit (profile): {sig(PL['Dc_prof'])}")
    if np.isfinite(PL['Dc_trace']):
        _bands(ax, PL, PL['Dc_trace'], WONG['reddishpurple'], f"permeation, transient fit (trace): {sig(PL['Dc_trace'])}")
    if PL['kappa_darcy'] and np.isfinite(PL['M_perm']):
        _bands(ax, PL, PL['kappa_darcy'][0] * PL['M_perm'], '0.35', rf"permeation, Darcy $\kappa \times M$ = {sig(PL['kappa_darcy'][0] * PL['M_perm'])}")
    ax.set_ylabel(r'$D_c$  ($\sigma^2/\tau$)')
    ax.set_yscale('log')
    ax.set_title(r'Cooperative diffusivity $D_c(\phi_s)$: compression holds vs steady permeation', fontsize=15)
    _phi_axis(ax, ccfg)
    smart_legend(ax, fontsize=10)
    return _save(fig, ccfg, 'Dc_vs_phi')


def fig_kappa(ccfg, X):
    """kappa against phi_s: compression D_c/M_cur per level (and the sweep notebook's D_c/M_secant),
    the permeation local kappa = D_c(phi_s)/M_cur(phi_s), and the Darcy membrane average."""
    S, C, PL = X['S'], X['C'], X['PL']
    sg = PL['seg']
    fig, ax = plt.subplots(figsize=(11, 7), constrained_layout=True)
    ok = np.isfinite(C['kappa'])
    ax.errorbar(S['phi'][ok], C['kappa'][ok], yerr=[C['kappa'][ok] - C['kappa_lo'][ok], C['kappa_hi'][ok] - C['kappa'][ok]],
                fmt='o-', color=WONG['blue'], ms=10, lw=2, capsize=5, label=r'compression: $D_c/M_\mathrm{cur}$')
    ax.plot(S['phi'][ok], C['kappa_sec'][ok], 's--', color='0.45', ms=7, lw=1.2, label=r'compression: $D_c/M_\mathrm{secant}$ (sweep notebook)')
    ax.plot(PL['phi_fit'], PL['kappa_fit'], '-', color=WONG['vermillion'], lw=2.2, alpha=0.8,
            label=rf"permeation, local: $q/(|\mathrm{{d}}\sigma'/\mathrm{{d}}\phi_s|\,\mathrm{{d}}\phi_s/\mathrm{{d}}z)$, degree-{ccfg.PERM_POLY_DEG} profile")
    _seg_points(ax, ccfg, sg, 'kappa')
    kd = PL['kappa_darcy']
    if kd:
        ax.fill_between([PL['phi_bot'], PL['phi_top']], kd[1], kd[2], color='0.35', alpha=0.2, lw=0)
        _bands(ax, PL, kd[0], '0.35', rf"permeation, Darcy $QL/(A\,\Delta P)$ = {sig(kd[0])} (membrane average)")
    ax.set_ylabel(r'$\kappa$  ($\sigma^5/(\epsilon\,\tau)$)')
    ax.set_yscale('log')
    ax.set_title(r'Hydraulic permeability $\kappa(\phi_s)$: compression vs steady permeation', fontsize=15)
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
    a.plot(PL['z'], C['sig'](PL['phi']), 'o', color=WONG['blue'], ms=4, label=r"$\sigma'_\mathrm{compression}(\phi_s(z))$, bin by bin")
    a.plot(PL['zz'], PL['sig_z'], '-', color=WONG['blue'], lw=2.2, label='the same on the fitted profile')
    zl = np.array([P['z_mem_lo'], P['z_mem_hi']])
    a.plot(zl, [PL['dP'], 0.0], 'k--', lw=1.8, label=rf"linear pore-pressure drop: $\Delta P$ = {sig(PL['dP'])} at the support, 0 at the feed face")
    a.set(xlabel=r'$z$  ($\sigma$)   [support $\leftarrow$ $\rightarrow$ feed]', ylabel=r"network stress  $\sigma'_{zz}$")
    a.set_title(r"(a) does the membrane's $\phi_s(z)$ carry $\Delta P$ on the compression curve?", fontsize=14)
    a.grid(alpha=0.3)
    smart_legend(a, fontsize=10)
    b.bar([0], [PL['q']], yerr=[PL['q_h']], color=WONG['vermillion'], capsize=8, width=0.6)
    b.bar([1], [X['q_pred']], color=WONG['blue'], width=0.6)
    b.set_xticks([0, 1])
    b.set_xticklabels(['measured\n$Q/A$', 'from compression\n' + r'$\frac{1}{L}\int D_c\,\mathrm{d}\phi_s/\phi_p$'], fontsize=12)
    b.set_ylabel(r'flux  $q$  ($\sigma/\tau$)')
    b.set_title(f"(b) flux: predicted/measured = {X['q_pred'] / PL['q']:.2f}", fontsize=14)
    b.grid(axis='y', alpha=0.3)
    return _save(fig, ccfg, 'closure')
