"""Compression HOLD (load) D_c fits, archived 2026-11-17.

Everything lib/triaxial.py computed or drew from the consolidation transient of a compression hold, moved here
unchanged when the compression D_c became the UNLOAD (free re-swelling) fit only:

    fit_Dc                  held-slab modal fit of u_z(zeta, t)/L over the hold (tau_1 = L^2/(4 pi^2 D_c))
    fit_Dc_stress           D_c from the load-piston stress trace
    hold_adequacy,          was the hold long enough (fitted D_c vs the deck's sizing constant)
    print_hold_check
    fit_Dc_timevarying,     D_c free in time (hold and unload) / D_c(phi_s) from the strain PDE, with their figures
    fit_Dc_phi
    fig_Dc                  the hold consolidation figure (profiles, interior strain, stress trace, unload traces)
    fig_Dc_sweep_load       the sweep figure with the hold fits (the former fig_Dc_sweep)
    add_load_Dc             the block load_level ran: L['Dc'], L['Dc_stress'], L['kappa_hold'], L['kappa_stress']

Why it was dropped: the hold starts from the ramp's face skins plus an undrained interior strain, so its D_c moves
with what is assumed about the onset state (0.05-0.085 on run 7 against 0.12 from the re-swelling of the same gel).

To use it again:
    import sys; sys.path.insert(0, '../archive/scripts/lib')       # from scripts/
    import triaxial_load_Dc as ldc
    ldc.add_load_Dc(cfg, R, L); ldc.fig_Dc(cfg, R, L); ldc.print_hold_check(cfg, [L])
It runs against the current lib/triaxial.py (every name of that module is imported below); it is not covered by
the test scripts apart from tests/dc_convention_test.py, which uses fit_Dc.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'scripts' / 'lib'))
import triaxial as tri  # noqa: E402

globals().update({k: v for k, v in vars(tri).items() if not k.startswith('__')})


def _w_modes(zh, kk):
    """Displacement modes of a slab HELD between two drained plates at the same bath pressure
    (2026-09-28).  With both faces pinned (u''(0) = u''(1) = -q/D_c) and the pore pressure equal at
    both faces (eps(0) = eps(1), uniform total stress) the strain modes are periodic, cos/sin(2 m pi
    zeta), so u relaxes in sin(pi kk zeta) (antisymmetric) and cos(pi kk zeta) - 1 (symmetric) with
    kk = 2 m, both zero at both plates, decaying as exp(-(pi kk)^2 D_c t/L^2).  Returns
    (len(zh), 2 len(kk)): the sin block, then the cos block.  The previous shapes
    (2 zeta - 1) + cos(k pi zeta), odd k, at k^2 pi^2 D_c/L^2 -- Terzaghi's constant-load series
    bent to vanish at both plates -- are not solutions of the held-slab equation and returned ~5x
    the true D_c on the exact relaxation (synthetic test, 2026-09-28)."""
    zh = np.asarray(zh, float)[:, None]
    return np.hstack([np.sin(np.pi * zh * kk[None, :]), np.cos(np.pi * zh * kk[None, :]) - 1.0])


def window_decay(lam, t, W):
    """exp(-lam t) averaged over the ave/chunk window [t - W, t] (clipped at t = 0, the hold
    onset) -- what a block-averaged snapshot tagged with its window END holds; W = 0 (or t = 0):
    the instantaneous value.  lam: array of rates, t: scalar."""
    lam = np.asarray(lam, float)
    w = min(W, t)
    if not w > 0:
        return np.exp(-lam * t)
    return (np.exp(-lam * (t - w)) - np.exp(-lam * t)) / (lam * w)


def fit_Dc(cfg, R, disp, free_offset=False):
    """Consolidation fit of D_c to u_z/L on the polymer domain during a hold between two
    drained plates: modes sin(2 m pi zeta), cos(2 m pi zeta) - 1 at 4 m^2 pi^2 D_c/L^2, free
    amplitudes (see _w_modes and the D_c notes in the notebooks).  Returns a dict or None.

    The domain is bounded by the HELD plate planes.  DL is the cumulative closure of
    the plate gap since the seated reference; f_sup is the share of it the support
    took (0 for the top-only drive of runs before 2026-09-18, 1/2 for the symmetric
    drive since).  Only the affine end state / plotted IC depend on f_sup: the fit
    itself models the hold-referenced u_dat with both faces pinned, so the odd-mode
    (centre-symmetric) expansion is the same in both cases -- the symmetric drive
    just makes the hold-onset state actually symmetric about the gel centre.

    The snapshots are ave/chunk block averages over the whole output interval, so each
    mode's decay is averaged over that window (DC_WINDOW_AVG, 2026-10-03); F['W'] is the
    window in LJ time and the model curves (F['T'], F['u_model']) are averaged the same way.

    Fit window (2026-10-09): F['Dc'] is fitted on the first DC_FIT_TAU1 x tau_1 of the hold
    (tau_1 of the fit itself, iterated), because a real hold is not one poroelastic transient:
    a slower process follows it and, fitted whole with these modes, pulls D_c down (quarter gel:
    0.109 over 6 tau_1, 0.054 whole).  F['Dc_all'] is the whole-hold value, F['stab'] a scan
    over windows and mode counts, F['stab_flag'] marks a > DC_STAB_FLAG disagreement.  The
    stress-trace D_c (fit_Dc_stress) is the independent check load_level prints next to it.

    free_offset (2026-10-09, fig_v0_check): project the per-snapshot rigid translation out of the data
    and the modes (centre both over the fit bins) -- algebraically the fit of the STRAIN PDE
    d eps/dt = D_c d2 eps/dz2, in which the barycentric velocity v0(t) = phi_s v_s + phi_p v_p
    (uniform in z) has dropped out, against the displacement fit above, whose modes fix v0(t)
    through the pinned-face boundary conditions."""
    if disp is None or not np.isfinite(R['z_support']):
        return None
    ts, z, Nc = disp['ts'], disp['z'], disp['Nc']
    z_sup = float(disp.get('z_supp_held', R['z_support']))
    if not np.isfinite(z_sup):
        z_sup = float(R['z_support'])
    span = disp['z_pist_held'] - z_sup
    L = disp.get('L_bb', span - 2.0)
    gap = 0.5 * (span - L)
    z_perm, z_feed = z_sup + gap, disp['z_pist_held'] - gap
    DL_pist = R['z_piston'] - disp['z_pist_held']      # piston travel DOWN since the reference
    DL_sup = z_sup - R['z_support']                    # support travel UP since the reference
    DL = DL_pist + DL_sup                              # total plate-gap closure
    if not (L > 0 and DL > 0):
        return None
    f_sup = float(DL_sup / DL)
    zeta = (z - z_perm) / L
    uhat = disp['uz'] / L
    idx = np.where((Nc.min(axis=0) > cfg.Ncount_min) & (zeta > 0) & (zeta < 1))[0]
    # DC_TRIM_BINS counts COARSE bins (binWidth): the same sigma of gel is dropped at each face on the fine file
    n_trim = int(round(cfg.DC_TRIM_BINS * cfg.binWidth / abs(float(z[1] - z[0])))) if len(z) > 1 else cfg.DC_TRIM_BINS
    if len(idx) < 4 + 2 * n_trim:
        return None
    if n_trim:
        idx = idx[n_trim:-n_trim]
    zf = zeta[idx]
    t_lj = (ts - disp['t_hold']) * cfg.dt_lj
    early_all = np.where(t_lj <= cfg.DC_FRAC_EARLY * t_lj[-1])[0]
    if len(early_all) < 2:
        return None
    if np.nanmax(np.abs(disp['uz'][:, idx])) > 0.5 * DL:
        print('  WARNING: |u_dat| ~ DL -- the disp file does not look hold-referenced.')
    if not cfg.DC_FREE_AMPS:
        print('  NOTE: DC_FREE_AMPS=False (Terzaghi load-control amplitudes) has no meaning for the held-slab '
              'modes (2026-09-28); fitting free amplitudes')
    kk = 2.0 * np.arange(1, cfg.DC_N_MODES + 1)            # even wavenumbers: m = 1..N, 2N amplitudes (sin + cos blocks)
    W = float(np.min(np.diff(ts))) * cfg.dt_lj if (cfg.DC_WINDOW_AVG and len(ts) > 1) else 0.0
    W = float(disp.get('W', W))                            # a subsampled snapshot set keeps the deck's averaging window (fig_Dc_skin)
    Y = uhat[:, idx]                                       # (n_snap, n_bins) data

    def decay(Dc, t, k=kk):                                # per mode, window-averaged like the data; t scalar or array
        return _decay_mat((np.pi * k) ** 2 * Dc / L ** 2, t, W)

    def _fit(sel, k):
        """(Dc, A, R2) of the modal fit (wavenumbers k) to the snapshots `sel`: amplitudes by linear least
        squares at each trial D_c, D_c by bounded Brent on the residual."""
        Wm = _w_modes(zf, k)                               # (n_bins, 2K)
        if free_offset:                                    # strain-PDE fit: rigid translation projected out
            Wm = Wm - Wm.mean(axis=0, keepdims=True)
            y = (Y[sel] - Y[sel].mean(axis=1, keepdims=True)).ravel()
        else:
            y = Y[sel].ravel()

        def X(Dc):
            d = np.tile(decay(Dc, t_lj[sel], k), (1, 2)) - 1.0      # (n_sel, 2K)
            return (Wm[None, :, :] * d[:, None, :]).reshape(-1, Wm.shape[1])

        def resid(Dc):
            Xd = X(Dc)
            return float(np.sum((Xd @ _lsq_amps(Xd, y) - y) ** 2))

        # coarse log grid first (short windows have secondary minima), then bounded Brent in the best cell
        grid = np.geomspace(cfg.DC_BOUNDS[0], cfg.DC_BOUNDS[1], 61)
        j = int(np.argmin([resid(g) for g in grid]))
        Dc = float(minimize_scalar(resid, bounds=(grid[max(j - 1, 0)], grid[min(j + 1, len(grid) - 1)]), method='bounded').x)
        Xd = X(Dc)
        A = _lsq_amps(Xd, y)
        ss_t = np.sum((y - np.mean(y)) ** 2)
        R2 = float(1.0 - np.sum((Xd @ A - y) ** 2) / ss_t) if ss_t > 1e-30 else np.nan
        return Dc, A, R2

    def _window(Dc, n_tau):
        """snapshots within the first n_tau tau_1 of the hold (at least 3, never beyond early_all)"""
        sel = early_all[t_lj[early_all] <= n_tau * L ** 2 / (4.0 * np.pi ** 2 * Dc)]
        return early_all[:3] if len(sel) < 3 else sel

    # ---- the whole-hold fit (the pre-2026-10-09 result) and the self-consistent early window ----
    Dc_all, A_all, R2_all = _fit(early_all, kk)
    n_tau = float(cfg.DC_FIT_TAU1)
    early, Dc, A, R2 = early_all, Dc_all, A_all, R2_all
    if n_tau > 0:
        seen = {len(early_all)}
        for _ in range(12):
            sel = _window(Dc, n_tau)
            if len(sel) == len(early) and np.array_equal(sel, early):
                break
            if len(sel) in seen:                           # oscillating between two window lengths: keep the current fit
                break
            seen.add(len(sel))
            early, (Dc, A, R2) = sel, _fit(sel, kk)
    tau1 = L ** 2 / (4.0 * np.pi ** 2 * Dc)
    stab = []
    if cfg.DC_STABILITY:
        for nt in (1.0, 2.0, 4.0, 10.0):
            sel = _window(Dc, nt)
            d, _, r = _fit(sel, kk)
            stab.append(dict(what=f'{nt:g} tau_1', n=int(len(sel)), N=int(cfg.DC_N_MODES), Dc=d, R2=r))
        stab.append(dict(what='whole hold', n=int(len(early_all)), N=int(cfg.DC_N_MODES), Dc=Dc_all, R2=R2_all))
        for N in (3, 8):
            if N != cfg.DC_N_MODES:
                d, _, r = _fit(early, 2.0 * np.arange(1, N + 1))
                stab.append(dict(what=f'window, {N} modes', n=int(len(early)), N=N, Dc=d, R2=r))
    stab_flag = bool(n_tau > 0 and len(early) < len(early_all) and abs(Dc_all / Dc - 1.0) > cfg.DC_STAB_FLAG)

    def T(zh, t):                                          # window-averaged like the snapshots; instantaneous at t = 0 (the IC)
        return _w_modes(zh, kk) @ (A * np.tile(decay(Dc, t)[0], 2))

    def T_all(zh, t):
        return _w_modes(zh, kk) @ (A_all * np.tile(decay(Dc_all, t)[0], 2))

    # absolute u_z/L (referenced to the pre-drive state): affine end state
    # (DL/L)(f_sup - zeta) -- the support face moved UP by DL_sup, the piston face
    # DOWN by DL_pist -- plus the fitted transient
    u_model = lambda zh, t: (DL / L) * (f_sup - zh) + T(zh, t)
    u_IC = lambda zh: u_model(zh, 0.0)
    u_tr = lambda zh, t: T(zh, t) - T(zh, 0.0)             # u_z/L since the hold onset: what the file holds
    u_tr_all = lambda zh, t: T_all(zh, t) - T_all(zh, 0.0)
    beta = np.nan
    hold_T = float(t_lj[-1])
    return dict(Dc=Dc, A=A, beta=beta, R2=R2, L=L, DL=DL, DL_pist=DL_pist, DL_sup=DL_sup, f_sup=f_sup,
                gap=gap, z_perm=z_perm, z_feed=z_feed, z_sup=z_sup, z_pist=disp['z_pist_held'],
                zeta=zeta, idx=idx, zf=zf, uhat=uhat, early=early, t_lj=t_lj, ts=ts, t_hold=float(disp['t_hold']),
                T=T, u_model=u_model, u_IC=u_IC, u_tr=u_tr, u_tr_all=u_tr_all, kk=kk, hold_T=hold_T, W=W,
                fine=bool(disp.get('fine', False)),
                shown=early[::max(1, int(np.ceil(len(early) / max(cfg.DC_PLOT_MAX, 1))))],
                # 2026-10-09: the fit window, the whole-hold fit and the stability scan
                fit_tau1=n_tau, tau1=tau1, window_T=float(t_lj[early[-1]]), n_tau1_fit=float(t_lj[early[-1]] / tau1),
                early_all=early_all, Dc_all=Dc_all, A_all=A_all, R2_all=R2_all, stab=stab, stab_flag=stab_flag,
                free_offset=bool(free_offset), u_mean=Y.mean(axis=1),           # <u>/L over the fit bins, every snapshot
                hold_check=hold_adequacy(cfg, L, hold_T, Dc))


def interior_strain(F, U):
    """-d(u_z/L)/d zeta over the fit domain by least squares, per snapshot: the interior strain change since
    the hold onset (compression positive).  U: (n_snap, n_bins) of u_z/L on F['zf']."""
    m1 = F['zf'] - F['zf'].mean()
    return -(np.atleast_2d(U) * m1[None, :]).sum(axis=1) / np.sum(m1 ** 2)


def fit_Dc_stress(cfg, steps, P, L, t_hold):
    """D_c from the load-piston stress trace of a hold (2026-10-09), independent of the displacement profiles.
    The load piston is solvent-transparent, so P = F_load/A is the network stress at the held face, M eps(L, t);
    in the held slab eps(L, t) - eps_bar relaxes in the cos(2 m pi zeta) modes, every one at 4 m^2 pi^2 D_c/L^2:
        P(t) = P_inf + a sum_m g_m exp(-4 m^2 pi^2 D_c t/L^2),   g_m = exp(-(2 m pi delta/L)^2 / 2),
    the hold-onset state being a Gaussian skin of width delta at each face (the ramp compresses the faces only).
    (a) the series (300 modes) is fitted for (P_inf, a, D_c, delta) on t >= DC_STRESS_TMIN; (b) a single
    exponential P_inf + B exp(-t/tau) on the tail t >= DC_STRESS_TAIL_TAU1 x tau_1, where only m = 1 survives
    whatever the onset state.  Standard errors are scaled by the residual autocorrelation time.  Returns a dict;
    ok = False with `why` when the trace cannot resolve tau_1 (too few blocks, or tau_1 shorter than 2 blocks or
    2 DC_STRESS_TMIN -- the quarter gel)."""
    steps = np.asarray(steps, float)
    P = np.asarray(P, float)
    t = (steps - t_hold) * cfg.dt_lj
    keep = (t > 0) & np.isfinite(P)
    t, P = t[keep], P[keep]
    out = dict(ok=False, why='', tail_ok=False, t=t, P=P, L=L, t_min=float(cfg.DC_STRESS_TMIN))
    if len(t) < 40:
        out['why'] = 'fewer than 40 stress blocks in the hold'
        return out
    block = float(np.median(np.diff(t)))
    out['block'] = block
    w = t >= cfg.DC_STRESS_TMIN
    if w.sum() < 40:
        out['why'] = f'fewer than 40 blocks after t_min = {cfg.DC_STRESS_TMIN:.0f} tau'
        return out
    m = np.arange(1, 301, dtype=float)

    def series(tt, P_inf, a, Dc, delta):
        g = np.exp(-0.5 * (2.0 * m * np.pi * delta / L) ** 2)
        return P_inf + a * (np.exp(-4.0 * (m[None, :] ** 2) * np.pi ** 2 * Dc * np.atleast_1d(tt)[:, None] / L ** 2) * g[None, :]).sum(axis=1)

    P_end = float(np.mean(P[t >= 0.75 * t[-1]]))
    p0 = (P_end, max(float(P[w][0] - P_end), 1e-4), cfg.DC_SLOW_REF / 4.0, 2.0)
    try:
        popt, pcov = curve_fit(series, t[w], P[w], p0=p0, maxfev=20000,
                               bounds=((-np.inf, 0.0, 1e-4, 0.3), (np.inf, np.inf, 10.0, 0.5 * L)))
    except (RuntimeError, ValueError) as e:
        out['why'] = f'series fit did not converge ({e})'
        return out
    res = P[w] - series(t[w], *popt)
    tau_ac = autocorr_time(res)
    se = np.sqrt(np.abs(np.diag(pcov))) * np.sqrt(max(tau_ac, 1.0))
    Dc, delta = float(popt[2]), float(popt[3])
    tau1 = L ** 2 / (4.0 * np.pi ** 2 * Dc)
    out.update(P_inf=float(popt[0]), a=float(popt[1]), Dc=Dc, Dc_se=float(se[2]), delta=delta, delta_se=float(se[3]),
               tau1=tau1, rms=float(np.sqrt(np.mean(res ** 2))), n=int(w.sum()), tau_ac=float(tau_ac),
               model=lambda tt, p=popt: series(tt, *p))
    # resolved only if the fitted mode is slower than the blocks AND most of the relaxation is still to come at
    # t_min: on the quarter gel (tau_1 ~ 160 tau) the trace has relaxed before t_min and the series would fit the
    # slower secondary process instead
    P_0 = float(np.mean(P[:2]))
    frac_after = float((series(np.array([cfg.DC_STRESS_TMIN]), *popt)[0] - popt[0]) / max(P_0 - popt[0], 1e-12))
    out['frac_after_tmin'] = frac_after
    if tau1 < 2.0 * block or tau1 < 2.0 * cfg.DC_STRESS_TMIN:
        out['why'] = f'tau_1 = {tau1:.0f} tau is not resolved (block {block:.0f} tau, t_min {cfg.DC_STRESS_TMIN:.0f} tau)'
        return out
    if frac_after < 0.3:
        out['why'] = (f'only {frac_after:.0%} of the stress relaxation is left at t_min = {cfg.DC_STRESS_TMIN:.0f} tau '
                      f'(the transient precedes the fitted window; the series then fits the slower process: tau_1 = {tau1:.0f})')
        return out
    out['ok'] = True
    # ---- the tail: m = 1 only ----
    t_tail = max(cfg.DC_STRESS_TAIL_TAU1 * tau1, cfg.DC_STRESS_TMIN)
    wt = t >= t_tail
    out.update(t_tail=float(t_tail), n_tail=int(wt.sum()))
    if wt.sum() >= 20 and (t[-1] - t_tail) > 1.5 * tau1:
        f = lambda tt, P_inf, B, tau: P_inf + B * np.exp(-tt / tau)
        B0 = float(series(np.array([t_tail]), *popt)[0] - popt[0])
        try:
            pt, ct = curve_fit(f, t[wt], P[wt], p0=(popt[0], max(B0, 1e-4), tau1), maxfev=20000,
                               bounds=((-np.inf, 0.0, block), (np.inf, np.inf, 50.0 * (t[-1] - t_tail))))
            rt = P[wt] - f(t[wt], *pt)
            st = np.sqrt(np.abs(np.diag(ct))) * np.sqrt(max(autocorr_time(rt), 1.0))
            tau_t = float(pt[2])
            out.update(tau_tail=tau_t, tau_tail_se=float(st[2]), Dc_tail=L ** 2 / (4.0 * np.pi ** 2 * tau_t),
                       Dc_tail_se=L ** 2 / (4.0 * np.pi ** 2 * tau_t ** 2) * float(st[2]), B_tail=float(pt[1]),
                       P_inf_tail=float(pt[0]), model_tail=lambda tt, p=pt: f(tt, *p),
                       tail_ok=bool(np.isfinite(tau_t) and tau_t < (t[-1] - t_tail) and pt[1] > 0))
        except (RuntimeError, ValueError):
            pass
    return out


def hold_residual(T, tau1, f):
    """(end residual, mean residual over the last fraction f of a hold of length T)."""
    end = float(np.exp(-T / tau1))
    avg = float((tau1 / (f * T)) * (np.exp(-(1.0 - f) * T / tau1) - np.exp(-T / tau1)))
    return end, avg


def n_tau_needed(tau1, target, f):
    n = np.arange(1.0, 30.001, 0.05)
    ok = [x for x in n if hold_residual(x * tau1, tau1, f)[1] <= target]
    return float(ok[0]) if ok else float('nan')


def hold_adequacy(cfg, L, hold_T, Dc_fit):
    """Was the hold long enough?  'fit': tau_1 = L^2/(4 pi^2 D_c) for the fitted D_c (held
    slab, first mode sin(2 pi zeta); 2026-09-28).  'slow': the DECK's own sizing formula,
    L^2/(pi^2 Dc_est) with Dc_est = DC_SLOW_REF -- the hold it actually prescribed (its
    Dc_est = 0.17 in that formula is the hold of D_c = 0.0425 in the held slab: adequate up to
    eps ~ 0.2, short at higher strain where the fitted D_c drops faster than the deck's (1-eps)^2 sizing).
    residual = mean excess stress over the last plateau_frac."""
    out = {}
    for tag, Dx, fac in (('fit', Dc_fit, 4.0), ('slow', cfg.DC_SLOW_REF, 1.0)):
        tau1 = L * L / (fac * np.pi ** 2 * Dx)
        end, avg = hold_residual(hold_T, tau1, cfg.plateau_frac)
        out[tag] = dict(Dc=Dx, tau1=tau1, end=end, avg=avg,
                        need=n_tau_needed(tau1, cfg.DC_TARGET_RESID, cfg.plateau_frac),
                        ok=avg <= cfg.DC_TARGET_RESID)
    return out


def _dc_time_panels(cfg, axc, axd, F, S, title_lvl=''):
    """The two time-scale panels of fig_Dc (2026-10-09).  (c) the interior strain change since the hold
    onset -- every snapshot, the windowed fit and the whole-hold fit; (d) the load-piston stress trace with
    the held-slab series and the tail exponential (fit_Dc_stress), or why it is absent."""
    t, zf = F['t_lj'], F['zf']
    s_dat = interior_strain(F, F['uhat'][:, F['idx']])
    td = np.geomspace(max(t[0] / 5.0, 1.0), t[-1], 300)       # starts before the first snapshot: the model's early shape
    s_win = np.array([interior_strain(F, F['u_tr'](zf, ti))[0] for ti in td])
    axc.semilogx(t, s_dat, 'o', color='k', ms=4 if len(t) < 60 else 2.5, alpha=0.7, label='snapshots (all)')
    axc.semilogx(td, s_win, '-', color=WONG['vermillion'], lw=2.4,
                 label=rf"fit to the first {F['n_tau1_fit']:.1f} $\tau_1$ ({len(F['early'])} snaps): $D_c={sig(F['Dc'])}$, $\tau_1={F['tau1']:.0f}\,\tau$")
    if F['fit_tau1'] > 0 and len(F['early']) < len(F['early_all']):
        s_all = np.array([interior_strain(F, F['u_tr_all'](zf, ti))[0] for ti in td])
        axc.semilogx(td, s_all, '--', color=WONG['blue'], lw=2.2,
                     label=rf"fit to the whole hold ({len(F['early_all'])} snaps): $D_c={sig(F['Dc_all'])}$")
        axc.axvspan(t[0], F['window_T'], color=WONG['vermillion'], alpha=0.07)
    axc.set(xlabel=r'hold time  ($\tau$)', ylabel=r'interior strain change  $-\partial(u_z/L)/\partial\zeta$')
    axc.xaxis.set_minor_formatter(NullFormatter())
    axc.set_title(r'(c) the transient in time' + (' -- two time scales' if F['stab_flag'] else ''), fontsize=14,
                  color=WONG['vermillion'] if F['stab_flag'] else 'k')
    axc.grid(alpha=0.3)
    smart_legend(axc, fontsize=10)
    if S is None:
        axd.text(0.5, 0.5, 'no load-piston trace (piston_force_avg) or DC_STRESS off', ha='center', va='center', transform=axd.transAxes)
        axd.set_title('(d) load-piston stress trace', fontsize=14)
        return
    tt, P = S['t'], S['P']
    axd.loglog(tt, P, '.', color='0.6', ms=3, label=f"$F_\\mathrm{{load}}/A$, blocks of {S.get('block', np.nan):.0f} $\\tau$")
    if S['ok']:
        axd.loglog(td[td >= S['t_min']], S['model'](td[td >= S['t_min']]), '-', color=WONG['reddishpurple'], lw=2.4,
                   label=rf"held-slab series from {S['t_min']:.0f} $\tau$: $D_c={sig(S['Dc'])}\pm{sig(S['Dc_se'], 1)}$, $\tau_1={S['tau1']:.0f}\,\tau$")
        if S['tail_ok']:
            tt2 = td[td >= S['t_tail']]
            axd.loglog(tt2, S['model_tail'](tt2), '-', color=WONG['skyblue'], lw=2.4,
                       label=rf"tail from {S['t_tail']:.0f} $\tau$: $\tau={S['tau_tail']:.0f}\pm{S['tau_tail_se']:.0f}$ $\to$ $D_c={sig(S['Dc_tail'])}$")
        # the profile fit's tau_1 on the same tail, amplitude fitted: is the stress slower or faster than the profiles?
        wt = tt >= S.get('t_tail', S['t_min'])
        if wt.sum() >= 10:
            Xp = np.column_stack([np.ones(wt.sum()), np.exp(-tt[wt] / F['tau1'])])
            Ap = np.linalg.lstsq(Xp, P[wt], rcond=None)[0]
            tt3 = td[td >= tt[wt][0]]
            axd.loglog(tt3, Ap[0] + Ap[1] * np.exp(-tt3 / F['tau1']), ':', color=WONG['vermillion'], lw=2.2,
                       label=rf"decaying at the profile fit's $\tau_1={F['tau1']:.0f}\,\tau$ (amplitude fitted)")
        lo = np.nanpercentile(P, 0.5)
        axd.set_ylim(max(lo * 0.8, 1e-4), np.nanpercentile(P[tt < max(10 * S['t_min'], tt[0])], 99.5) * 1.3)
        axd.set_title(rf"(d) stress trace: $D_c={sig(S['Dc'])}$ vs profiles {sig(F['Dc'])}  (ratio {F['Dc'] / S['Dc']:.2f})", fontsize=14)
    else:
        axd.set_title('(d) stress trace: $\\tau_1$ not resolved by this trace', fontsize=14)
        annotate_box(axd, S['why'].replace('; ', ';\n').replace(' (', '\n('), loc='lower left', fontsize=10, color='0.3')
    axd.set(xlabel=r'hold time  ($\tau$)', ylabel=r'$F_\mathrm{load}/A$  (LJ)')
    axd.xaxis.set_minor_formatter(NullFormatter())
    axd.yaxis.set_minor_formatter(NullFormatter())
    ylo, yhi = axd.get_ylim()
    axd.set_yticks([v for v in (0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0) if ylo <= v <= yhi])
    axd.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    axd.grid(alpha=0.3, which='both')
    axd.legend(loc='upper right', fontsize=10)


def _unload_trace_panel(cfg, axe, F, Uf, S):
    """(e) of fig_Dc (2026-10-09; Rg trace added 2026-10-15): the thickness traces of the free re-swelling -- bounding box
    (the deck's L_bb, the outermost beads) and Rg (bulk-weighted) -- each with its one-parameter series fit, the free-asymptote
    fit of the bb trace, and the series drawn at the hold's and the re-swelling profile D_c values.  The profile fit is the
    primary D_c; the bb trace reads high because the outermost beads spring back within ~100 tau of the retraction."""
    T, Tr, Pf, L = Uf['trace'], Uf.get('trace_rg'), Uf['prof'], Uf['L']
    if T is None and Tr is None:
        axe.text(0.5, 0.5, 'no thickness trace (gel_dimensions_{bb,rg}_*_u<lvl>)', ha='center', va='center', transform=axe.transAxes)
        return
    T0 = T if T is not None else Tr
    tt = T0['t']
    td = np.geomspace(max(tt[tt > 0][0], 1.0), tt[-1], 400)
    if T is not None:
        axe.semilogx(tt[tt > 0], T['dL'][tt > 0], '.', color='0.6', ms=2, label=r'$L_\mathrm{bb}(t) - L_\mathrm{hold}$ (bounding box, every print)')
        axe.semilogx(td, T['model'](td), '-', color=WONG['vermillion'], lw=2.4,
                     label=rf"bb series, $\Delta L_\infty$ fixed = {T['dL_inf']:.2f} $\sigma$: $D_c={sig(T['Dc'])}\pm{sig(T['Dc_se'], 1)}$ (reads high: surface beads)")
        axe.semilogx(td, T['model_free'](td), '--', color=WONG['vermillion'], lw=1.4, alpha=0.8,
                     label=rf"bb series, $\Delta L_\infty$ free = {T['dL_inf_free']:.2f}: $D_c={sig(T['Dc_free'])}$")
    if Tr is not None:
        tr = Tr['t']
        axe.semilogx(tr[tr > 0], Tr['dL'][tr > 0], '.', color=WONG['skyblue'], ms=2, alpha=0.6, label=r'$L_\mathrm{Rg}(t) - L_\mathrm{hold,Rg}$ (Rg thickness)')
        axe.semilogx(td, Tr['model'](td), '-', color=WONG['blue'], lw=2.2,
                     label=rf"Rg series, $\Delta L_\infty$ fixed = {Tr['dL_inf']:.2f} $\sigma$: $D_c={sig(Tr['Dc'])}\pm{sig(Tr['Dc_se'], 1)}$")
    dLi = T0['dL_inf']
    if Pf is not None:
        axe.semilogx(td, dLi * _unload_series(td, Pf['Dc'], L), '-', color='k', lw=1.8, alpha=0.8,
                     label=rf"bb series at the PRIMARY (profile) $D_c={sig(Pf['Dc'])}$")
    for D_, lab, ls in ((F['Dc'], "hold, windowed profile", ':'), (F['Dc_all'], "hold, whole-hold profile", '-.')):
        axe.semilogx(td, dLi * _unload_series(td, D_, L), ls, color='0.3', lw=1.4, label=rf"bb series at the {lab} $D_c={sig(D_)}$")
    if S is not None and S.get('ok'):
        axe.semilogx(td, dLi * _unload_series(td, S['Dc'], L), ls='--', color=WONG['green'], lw=1.4,
                     label=rf"bb series at the hold's stress-trace $D_c={sig(S['Dc'])}$")
    axe.axhline(dLi, color='k', lw=0.8, alpha=0.5)
    axe.axvline(Uf['retract_T'], color='k', lw=0.8, ls=':', alpha=0.6)
    axe.text(Uf['retract_T'] * 1.1, 0.93 * dLi, 'plates stopped', fontsize=9, color='0.3')
    axe.set(xlabel=r'time since retraction start  ($\tau$)', ylabel=r'$\Delta L(t)$  ($\sigma$)')
    axe.xaxis.set_minor_formatter(NullFormatter())
    axe.set_title(r'(e) re-swelling thickness traces ($\tau_1=L^2/\pi^2 D_c$): checks on the profile fit -- '
                  + (rf"bb {sig(T['Dc'])}, " if T is not None else '') + (rf"Rg {sig(Tr['Dc'])}, " if Tr is not None else '')
                  + (rf"profiles {sig(Pf['Dc'])}" if Pf is not None else ''), fontsize=13)
    axe.grid(alpha=0.3)
    axe.legend(loc='lower right', fontsize=9)


# ---------------------------------------------------------------------------
#  D_c(t): the same modal problems with the diffusivity free in time (2026-10-18).
#  The strain equation is linear with a time-only coefficient, so D_c(t) enters only through the
#  diffusion time theta(t) = int_0^t D_c dt': every mode decays as exp(-lam_k theta), lam_k = mu_k^2/L^2.
#  Fit: D_c piecewise constant over DC_TV_BINS equal time bins of the identifiable window (the first
#  DC_TV_TAU1 tau_1), theta(t) its running integral, amplitudes A_k shared by all snapshots (linear least
#  squares inside the residual), log D_k by least squares with standard errors from the Jacobian.  A free
#  theta_i per snapshot (at the fitted amplitudes) is drawn behind it as the unregularised reading; a
#  monotone-pooled version of those theta_i was tried first and its flat runs made d theta/dt useless.
#  Window averaging of the snapshots is not modelled here: each snapshot is placed at its window centre.
# ---------------------------------------------------------------------------
def _theta_free(Wm, Y, lam, A, theta0, th_hi):
    """per-snapshot diffusion time at FIXED amplitudes (bounded Brent) with its standard error from the model's
    sensitivity to theta_i -- the unregularised reading drawn as dots behind the piecewise fit"""
    theta = np.asarray(theta0, float).copy(); se = np.empty(len(Y))
    for i in range(len(Y)):
        def r(th, i=i):
            return float(np.sum((Wm[i] @ (A * (np.exp(-lam * th) - 1.0)) - Y[i]) ** 2))
        theta[i] = minimize_scalar(r, bounds=(0.0, th_hi), method='bounded').x
        g = Wm[i] @ (A * lam * np.exp(-lam * theta[i]))
        S = float(np.sum(g ** 2)); n = len(Y[i])
        se[i] = np.sqrt(r(theta[i]) / max(n - 1, 1) / S) if S > 0 else np.inf
    return theta, se


def _dc_bins_fit(Wm, Y, lam, t, D0, K):
    """D_c piecewise constant over K equal time bins on [0, t_max]: theta(t) = sum_k D_k x (overlap of [0, t] with bin k).
    u_i(zeta) = sum_k A_k w_k(zeta) [exp(-lam_k theta_i) - 1], amplitudes by linear least squares inside the residual
    (variable projection), log D_k by least_squares.  Returns edges, D (K), se of log D, A, R2, theta at the snapshots."""
    from scipy.optimize import least_squares
    t = np.asarray(t, float)
    edges = np.linspace(0.0, t.max() * (1 + 1e-9), K + 1)
    over = np.clip(t[:, None] - edges[None, :-1], 0.0, None) - np.clip(t[:, None] - edges[None, 1:], 0.0, None)   # (n, K)
    y = np.concatenate(Y)

    def X(theta):
        return np.vstack([Wm[i] * (np.exp(-lam * theta[i]) - 1.0)[None, :] for i in range(len(Y))])

    def resid(logD):
        Xd = X(over @ np.exp(logD))
        return Xd @ _lsq_amps(Xd, y) - y

    r = least_squares(resid, np.full(K, np.log(D0)), bounds=(np.log(1e-6), np.log(50.0)), x_scale='jac', max_nfev=300)
    D = np.exp(r.x)
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac) * (np.sum(r.fun ** 2) / max(len(r.fun) - K, 1))
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        se = np.full(K, np.nan)
    theta = over @ D
    Xd = X(theta); A = _lsq_amps(Xd, y)
    ss = np.sum((y - y.mean()) ** 2)
    R2 = float(1.0 - np.sum((Xd @ A - y) ** 2) / ss) if ss > 1e-30 else np.nan
    return edges, D, se, A, R2, theta


def fit_Dc_timevarying(cfg, L, n_max=300, K=None):
    """D_c(t) of the hold (modes of fit_Dc) and of the free re-swelling (modes of fit_Dc_unload) with the
    diffusivity piecewise constant over K equal time bins of the identifiable window (see the section comment).
    Stores and returns L['Dc_tv'] = {'hold': {...}, 'unload': {...}} with the bin edges, D (per bin) and its
    standard error, theta at the snapshots, the per-snapshot free theta_i (dots), the constant-D_c reference."""
    K = int(K or cfg.DC_TV_BINS)
    n_id = float(cfg.DC_TV_TAU1)                        # snapshots within n_id tau_1 of the constant fit: beyond, the modes
    out = {}                                            # have decayed and theta is not identifiable
    F = L.get('Dc')
    if F is not None:
        sel = F['early_all'][F['t_lj'][F['early_all']] <= n_id * F['tau1']]
        sel = sel[::max(1, int(np.ceil(len(sel) / n_max)))]
        t = np.maximum(F['t_lj'][sel] - 0.5 * F['W'], 0.0)
        lam = np.tile((np.pi * F['kk'] / F['L']) ** 2, 2)
        Wm = [_w_modes(F['zf'], F['kk'])] * len(sel)
        Y = [F['uhat'][i][F['idx']] for i in sel]
        edges, D, se, A, R2, theta = _dc_bins_fit(Wm, Y, lam, t, F['Dc'], K)
        th_i, se_i = _theta_free(Wm, Y, lam, A, theta, 3.0 * theta.max() + 8.0 / lam.min())
        out['hold'] = dict(t=t, edges=edges, D=D, se_logD=se, theta=theta, theta_i=th_i, se_i=se_i, Dc_const=F['Dc'], tau1=F['tau1'],
                           R2=R2, R2_const=F['R2'], n=len(sel), K=K, L=F['L'], n_id=n_id, label=r'hold, $\tau_1 = L^2/(4\pi^2 D_c)$')
    Uf = L.get('unload')
    Pf = (Uf or {}).get('prof')
    if Pf is not None:
        sel = np.where(Pf['t_lj'] <= n_id * Pf['tau1'])[0]
        sel = sel[::max(1, int(np.ceil(len(sel) / n_max)))]
        t = np.maximum(Pf['t_lj'][sel] - 0.5 * Pf['W'], 0.0)
        lam = (np.pi * Pf['k'] / Uf['L']) ** 2
        Wm = [np.cos(np.pi * np.outer(Pf['zl'][i], Pf['k'])) for i in sel]
        Y = [Pf['yl'][i] for i in sel]
        edges, D, se, A, R2, theta = _dc_bins_fit(Wm, Y, lam, t, Pf['Dc'], K)
        th_i, se_i = _theta_free(Wm, Y, lam, A, theta, 3.0 * theta.max() + 8.0 / lam.min())
        out['unload'] = dict(t=t, edges=edges, D=D, se_logD=se, theta=theta, theta_i=th_i, se_i=se_i, Dc_const=Pf['Dc'], tau1=Pf['tau1'],
                             R2=R2, R2_const=Pf['R2'], n=len(sel), K=K, L=Uf['L'], n_id=n_id, label=r're-swelling, $\tau_1 = L^2/(\pi^2 D_c)$')
    L['Dc_tv'] = out
    for k, d in out.items():
        print(f"  D_c(t) {k:9s}: {d['n']} snapshots within {d['n_id']:g} tau_1, {d['K']} bins of {d['edges'][1]:.0f} tau, R^2 {d['R2']:.3f} "
              f"(constant D_c {d['Dc_const']:.4f}: {d['R2_const']:.3f});  D_c per bin: "
              + ', '.join(f"{v:.3f}" for v in d['D']) + f"  (x/ {np.exp(np.nanmedian(d['se_logD'])):.2f} typical)")
    return out


def fig_Dc_timevarying(cfg, R, L):
    """Third D_c figure of the single-level notebook (2026-10-18): (a) the diffusion time theta(t) of the hold and of the
    re-swelling (the per-snapshot theta_i as dots, the binned-D_c fit as a line) against the constant-D_c lines, each
    transient's tau_1 marked; (b) the CUMULATIVE ratio theta_i/(D_c t) against t/tau_1 (dots; 1 at every t = constant
    D_c; a secant, so it smooths), and the LOCAL ratio D_k/D_c of the binned fit (faint steps = the slope of the line
    in (a) over D_c).  The steps carry a caveat: the bin heights trade off against the free amplitudes (they change
    with the number of bins at the same R^2) and the late bins are weakly determined once the modes have decayed
    (the printed x/ errors are conditional on the amplitudes), so theta is what the data determine, D_k is a reading.
    2026-10-18: the "mean to 1.5 tau_1" lines were dropped (a mean of the secant ratio over a window says nothing
    the dots do not)."""
    TV = L.get('Dc_tv') or fit_Dc_timevarying(cfg, L)
    if not TV:
        print('D_c(t) figure skipped (no profile fits)')
        return None
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(18, 6.5), constrained_layout=True)
    cols = {'hold': WONG['blue'], 'unload': WONG['vermillion']}
    for k, d in TV.items():
        c = cols[k]
        ok = np.isfinite(d['se_i']) & (d['se_i'] < 0.5 * np.maximum(d['theta_i'], 1e-9) + 0.05 * d['theta'].max())
        axa.errorbar(d['t'][ok], d['theta_i'][ok], yerr=d['se_i'][ok], fmt='o', color=c, ms=2.5, alpha=0.3, elinewidth=0.6, capsize=0,
                     label=f"{d['label']}: $\\theta_i$ per snapshot at the fitted amplitudes ($\\pm$ s.e.; {np.sum(~ok)} unresolved left out)")
        tt = np.linspace(0, d['t'].max(), 300)
        over = np.clip(tt[:, None] - d['edges'][None, :-1], 0.0, None) - np.clip(tt[:, None] - d['edges'][None, 1:], 0.0, None)
        axa.plot(tt, over @ d['D'], '-', color=c, lw=2.2, label=f"{d['K']}-bin $D_c(t)$ fit ($R^2$ {d['R2']:.3f}; first {d['n_id']:g} $\\tau_1$)")
        axa.plot(tt, d['Dc_const'] * tt, '--', color=c, lw=1.5, label=f"constant $D_c$ = {sig(d['Dc_const'])} ($R^2$ {d['R2_const']:.3f})")
        axa.axvline(d['tau1'], color=c, ls=':', lw=1.0, label=f"$\\tau_1$ = {d['tau1']:.0f} $\\tau$ ($L$ = {d['L']:.1f}, {d['n']} snapshots)")
        tn = d['t'] / d['tau1']
        ratio = d['theta_i'] / (d['Dc_const'] * np.maximum(d['t'], 1e-9))
        rse = d['se_i'] / (d['Dc_const'] * np.maximum(d['t'], 1e-9))
        okb = ok & (tn > 0.05)
        axb.errorbar(tn[okb], ratio[okb], yerr=rse[okb], fmt='o', color=c, ms=3, alpha=0.6, elinewidth=0.7, capsize=0,
                     label=f"{d['label']}: cumulative $\\theta_i/(D_c t)$, constant-fit $D_c$ = {sig(d['Dc_const'])} ($t > 0.05\\,\\tau_1$)")
        se_typ = np.exp(np.nanmedian(d['se_logD'][np.isfinite(d['se_logD'])])) if np.isfinite(d['se_logD']).any() else np.nan
        _lab = f"{k.replace('unload', 're-swelling')}: {d['K']}-bin local $D_c(t)/D_c$ (slope of the (a) line; s.e. x/ {se_typ:.2f} at fixed amplitudes)"
        _vals, _edges = d['D'] / d['Dc_const'], d['edges'] / d['tau1']
        if hasattr(axb, 'stairs'):
            axb.stairs(_vals, _edges, baseline=None, color=c, lw=1.8, alpha=0.4, label=_lab)
        else:   # matplotlib < 3.4 (the Expanse anaconda3/2021.05 module): step plot of the same bins
            axb.step(_edges, np.append(_vals, _vals[-1]), where='post', color=c, lw=1.8, alpha=0.4, label=_lab)
    axb.axhline(1.0, color='k', ls='--', lw=1.2, label=r'$\theta = D_c t$ (constant $D_c$)')
    axb.axhspan(0.8, 1.2, color='0.5', alpha=0.12, label=r'$\pm 20\%$')
    axb.axvline(1.0, color='0.5', ls=':', lw=1)
    axa.set_xlabel(r'time since the onset  ($\tau$)'); axa.set_ylabel(r'$\theta(t) = \int_0^t D_c\,\mathrm{d}t^\prime$  ($\sigma^2$)')
    axa.set_title(r'(a) diffusion time: every mode decays as $e^{-\mu_k^2\theta/L^2}$, amplitudes shared by all snapshots', fontsize=12)
    axb.set_xlabel(r"$t/\tau_1$  (each transient's own $\tau_1$)"); axb.set_ylabel(r'$\theta_i\,/\,(D_c\,t)$')
    axb.set_title(r'(b) against the constant-$D_c$ line: cumulative $\theta_i/(D_c t)$ (dots) and local $\mathrm{d}\theta/\mathrm{d}t\,/D_c$ (steps); 1 everywhere = constant $D_c$', fontsize=12)
    axb.set_ylim(0, 1.6); axb.set_xlim(0, float(cfg.DC_TV_TAU1))
    axa.grid(alpha=0.3); axa.legend(fontsize=8.5, loc='upper left')
    hb, lb = axb.get_legend_handles_labels()
    order = [i for i, l in enumerate(lb) if 'cumulative' in l] + [i for i, l in enumerate(lb) if 'local' in l] \
        + [i for i, l in enumerate(lb) if 'cumulative' not in l and 'local' not in l]
    axb.grid(alpha=0.3); axb.legend([hb[i] for i in order], [lb[i] for i in order], fontsize=8.5, loc='upper left')
    fig.suptitle(f"Cooperative diffusivity with $D_c$ free in time (level _c{L['lvl']})  |  {cfg.sim_name}\n"
                 "hold: $\\sin, \\cos(2m\\pi\\zeta)$ modes between the held plates;  re-swelling: $\\cos(k\\pi\\zeta)$ modes, faces free",
                 fontsize=12, fontweight='bold')
    return _save(fig, cfg, 'Dc_timevarying', L['lvl'])


# ---------------------------------------------------------------------------
#  D_c(phi_s): both compression transients with ONE state-dependent diffusivity (2026-10-18).
#  v_p = q(t) + kappa dsigma'/dz with sigma' = sigma'(eps) gives du/dt = q(t) + D_c(eps) d2u/dz2,
#  D_c = kappa M_t (tangent modulus), and for the strain eps = -du/dz the conservative form
#      d eps/dt = d/dz [ D_c(eps) d eps/dz ],      phi_s = phi_s0 - (1 - phi_s0) eps      (S129).
#  D_c(phi_s) = D_ref exp[a (phi_s - phi_s0)]: D_ref at the swollen reference, a the sensitivity.
#  Hold: eps periodic (pinned faces, both drained into one bath), mean eps = DL/L fixed.  Re-swelling:
#  eps = 0 at both free faces.  Each transient's initial strain profile is free (periodic modes about
#  DL/L; a uniform value + sine modes) and is fitted together with (D_ref, a); the observable is the
#  displacement since the onset, u(zeta,t) - u(zeta,0) = -int_0^zeta [eps(t) - eps(0)] dzeta' (plus a
#  constant removed with the COM for the re-swelling).  Explicit finite differences.  (D_ref, a) are
#  shared in the joint fit and free per transient in the separate fits; the three are compared.
# ---------------------------------------------------------------------------
def _eps_pde(eps0, Dfun, L, t_out, periodic):
    """d eps/dt = d/dzeta [D(eps) d eps/dzeta] / L^2 on zeta in [0, 1] (len(eps0) points), explicit Euler at the
    stability limit, sampled at t_out (LJ, increasing).  periodic: zeta_j = j/n, j = 0..n-1;  else Dirichlet
    eps = 0 at both ends.  Returns (len(t_out), n)."""
    eps = np.asarray(eps0, float).copy()
    n = len(eps)
    dz = 1.0 / n if periodic else 1.0 / (n - 1)
    out = np.empty((len(t_out), n))
    t = 0.0
    if not periodic:
        eps[0] = eps[-1] = 0.0
    for m, tm in enumerate(t_out):
        while t < tm:
            if periodic:
                e1 = np.roll(eps, -1)
                Dm = Dfun(0.5 * (eps + e1))
                flux = Dm * (e1 - eps) / dz                     # between j and j+1
                deriv = (flux - np.roll(flux, 1)) / (dz * L ** 2)
            else:
                Dm = Dfun(0.5 * (eps[1:] + eps[:-1]))
                flux = Dm * (eps[1:] - eps[:-1]) / dz
                deriv = np.zeros(n)
                deriv[1:-1] = (flux[1:] - flux[:-1]) / (dz * L ** 2)
            dt = min(0.4 * dz ** 2 * L ** 2 / max(float(Dm.max()), 1e-12), tm - t)
            eps = eps + dt * deriv
            t += dt
        out[m] = eps
    return out


def _u_from_eps(eps_t, eps0, periodic):
    """u(zeta,t) - u(zeta,0) = -int_0^zeta [eps(t) - eps(0)] on the solver grid (trapezoid)."""
    d = eps_t - eps0[None, :]
    n = d.shape[1]
    if periodic:                                                # grid j/n, close the integral at zeta = 1
        d = np.hstack([d, d[:, :1]])
    dz = 1.0 / (d.shape[1] - 1)
    u = -np.hstack([np.zeros((d.shape[0], 1)), np.cumsum(0.5 * (d[:, 1:] + d[:, :-1]) * dz, axis=1)])
    return np.linspace(0.0, 1.0, d.shape[1]), u


def fit_Dc_phi(cfg, L, n_grid=81, n_max=80, n_modes=6):
    """One D_c(phi_s) for the hold and the free re-swelling (see the section comment).  Stores and returns
    L['Dc_phi'] with the joint fit ('joint'), the separate fits ('hold', 'unload'), each with D_ref, a, their
    standard errors, D_c at the swollen reference and at the hold state, R2 per transient, and the model
    profiles at the shown snapshots."""
    from scipy.optimize import least_squares
    F = L.get('Dc')
    Pf = (L.get('unload') or {}).get('prof')
    if F is None or Pf is None:
        print('D_c(phi_s) fit skipped (needs the hold profile fit and the re-swelling profile fit)')
        return None
    Uf = L['unload']
    phi0 = float(L['phi_cal'][0]) if L.get('phi_cal') is not None and False else float(cfg.PHI_S0)
    phip0 = 1.0 - phi0
    n_id = float(cfg.DC_TV_TAU1)
    # ---- data: hold (Eulerian bins zf, hold-referenced u/L) and re-swelling (material zeta per snapshot, COM removed) ----
    selh = F['early_all'][F['t_lj'][F['early_all']] <= n_id * F['tau1']]
    selh = selh[::max(1, int(np.ceil(len(selh) / n_max)))]
    th = np.maximum(F['t_lj'][selh] - 0.5 * F['W'], 0.0)
    Yh = [F['uhat'][i][F['idx']] for i in selh]
    zh = F['zf']
    selu = np.where(Pf['t_lj'] <= n_id * Pf['tau1'])[0]
    selu = selu[::max(1, int(np.ceil(len(selu) / n_max)))]
    tu = np.maximum(Pf['t_lj'][selu] - 0.5 * Pf['W'], 0.0)
    Yu = [Pf['yl'][i] for i in selu]
    Zu = [Pf['zl'][i] for i in selu]
    Lh, Lu = float(F['L']), float(Uf['L'])
    eps_h = float(F['DL'] / F['L'])                               # mean strain of the held gel
    eps_u0 = float(Uf['dL_inf_ref'] / Uf['L_hold'])             # uniform strain released by the re-swelling
    zg_p = np.arange(n_grid) / n_grid                            # periodic grid (hold)
    zg_d = np.linspace(0.0, 1.0, n_grid)                         # Dirichlet grid (re-swelling)
    mm = np.arange(1, n_modes + 1)
    Mh = np.hstack([np.cos(2 * np.pi * np.outer(zg_p, mm)), np.sin(2 * np.pi * np.outer(zg_p, mm))])   # (n, 2N)
    Mu = np.sin(np.pi * np.outer(zg_d, mm))                                                          # (n, N)
    wh, wu = 1.0 / np.sqrt(sum(len(y) for y in Yh)), 1.0 / np.sqrt(sum(len(y) for y in Yu))

    def Dfun_of(logD, a):
        return lambda e: np.exp(logD) * np.exp(-a * phip0 * e)   # phi_s - phi_s0 = -phi_p0 eps

    def model_hold(logD, a, amps):
        eps0 = eps_h + Mh @ amps
        E = _eps_pde(eps0, Dfun_of(logD, a), Lh, th, True)
        zz, U = _u_from_eps(E, eps0, True)
        return eps0, E, [np.interp(zh, zz, U[i]) for i in range(len(th))]

    def model_unload(logD, a, amps):
        eps0 = amps[0] + Mu @ amps[1:]
        E = _eps_pde(eps0, Dfun_of(logD, a), Lu, tu, False)
        zz, U = _u_from_eps(E, eps0, False)
        return eps0, E, [np.interp(Zu[i], zz, U[i]) - np.mean(np.interp(Zu[i], zz, U[i])) for i in range(len(tu))]

    def res_hold(logD, a, amps):
        return np.concatenate([(m - y - np.mean(m - y) * 0.0) for m, y in zip(model_hold(logD, a, amps)[2], Yh)]) * wh

    def res_unload(logD, a, amps):
        return np.concatenate([(m - (y - np.mean(y))) for m, y in zip(model_unload(logD, a, amps)[2], Yu)]) * wu

    nh, nu = 2 * n_modes, n_modes + 1
    x0h = np.zeros(nh); x0u = np.concatenate([[eps_u0], np.zeros(n_modes)])

    def _fit(which):
        if which == 'joint':
            f = lambda x: np.concatenate([res_hold(x[0], x[1], x[2:2 + nh]), res_unload(x[0], x[1], x[2 + nh:])])
            x0 = np.concatenate([[np.log(0.5 * (F['Dc'] + Pf['Dc'])), 0.0], x0h, x0u])
        elif which == 'hold':
            f = lambda x: res_hold(x[0], x[1], x[2:])
            x0 = np.concatenate([[np.log(F['Dc']), 0.0], x0h])
        else:
            f = lambda x: res_unload(x[0], x[1], x[2:])
            x0 = np.concatenate([[np.log(Pf['Dc']), 0.0], x0u])
        lo = np.full(len(x0), -np.inf); hi = np.full(len(x0), np.inf)
        lo[0], hi[0] = np.log(cfg.DC_BOUNDS[0]), np.log(cfg.DC_BOUNDS[1]); lo[1], hi[1] = -60.0, 60.0
        r = least_squares(f, x0, bounds=(lo, hi), x_scale='jac', max_nfev=400)
        J = r.jac
        try:
            cov = np.linalg.inv(J.T @ J) * (np.sum(r.fun ** 2) / max(len(r.fun) - len(x0), 1))
            se = np.sqrt(np.clip(np.diag(cov), 0, None))
        except np.linalg.LinAlgError:
            cov = None; se = np.full(len(x0), np.nan)
        return r, cov, se

    def _r2(res, Y, w):
        y = np.concatenate([yy - (np.mean(yy) if w is wu else 0.0) for yy in Y])
        ss = np.sum((y - y.mean()) ** 2)
        return float(1.0 - np.sum((res / w) ** 2) / ss) if ss > 1e-30 else np.nan

    out = dict(phi0=phi0, eps_h=eps_h, eps_u0=eps_u0, n_id=n_id, th=th, tu=tu, zh=zh, Yh=Yh, Yu=Yu, Zu=Zu, selh=selh, selu=selu,
               Dc_const=dict(hold=F['Dc'], unload=Pf['Dc']), R2_const=dict(hold=F['R2'], unload=Pf['R2']))
    for which in ('joint', 'hold', 'unload'):
        r, cov, se = _fit(which)
        x = r.x
        d = dict(logD=float(x[0]), a=float(x[1]), se_logD=float(se[0]), se_a=float(se[1]), cov=cov[:2, :2] if cov is not None else None,
                 D_ref=float(np.exp(x[0])), nfev=int(r.nfev), status=int(r.status))
        d['D_hold'] = d['D_ref'] * np.exp(-d['a'] * phip0 * eps_h)   # at the held gel's mean strain
        if which in ('joint', 'hold'):
            amps = x[2:2 + nh] if which == 'joint' else x[2:]
            eps0, E, U = model_hold(x[0], x[1], amps)
            d['hold'] = dict(eps0=eps0, E=E, U=U, amps=amps, R2=_r2(res_hold(x[0], x[1], amps), Yh, wh))
        if which in ('joint', 'unload'):
            amps = x[2 + nh:] if which == 'joint' else x[2:]
            eps0, E, U = model_unload(x[0], x[1], amps)
            d['unload'] = dict(eps0=eps0, E=E, U=U, amps=amps, R2=_r2(res_unload(x[0], x[1], amps), Yu, wu))
        out[which] = d
    L['Dc_phi'] = out
    print(f"  D_c(phi_s) = D_ref exp[a (phi_s - {phi0:.3f})], phi_s = phi_s0 - (1 - phi_s0) eps;  held gel eps = {eps_h:.4f} (phi_s {phi0 - phip0 * eps_h:.3f}), "
          f"re-swelling from eps = {eps_u0:.4f}")
    for which in ('joint', 'hold', 'unload'):
        d = out[which]
        print(f"    {which:7s}: D_ref = {d['D_ref']:.4f} (x/ {np.exp(d['se_logD']):.3f}),  a = {d['a']:+.1f} +/- {d['se_a']:.1f}  ->  D_c(hold state) = {d['D_hold']:.4f}"
              + (f";  R^2 hold {d['hold']['R2']:.3f}" if 'hold' in d else '') + (f", re-swelling {d['unload']['R2']:.3f}" if 'unload' in d else '')
              + f"   [constant-D_c fits: hold {F['Dc']:.4f} ({F['R2']:.3f}), re-swelling {Pf['Dc']:.4f} ({Pf['R2']:.3f})]")
    return out


def fig_Dc_phi(cfg, R, L):
    """Fourth D_c figure of the single-level notebook (2026-10-18): (a) D_c(phi_s) of the joint and the separate fits
    with the constant-D_c values at each transient's phi_s, (b) hold and (c) re-swelling profiles against the joint model."""
    X = L.get('Dc_phi') or fit_Dc_phi(cfg, L)
    if X is None:
        return None
    phi0, phip0 = X['phi0'], 1.0 - X['phi0']
    fig, (axa, axb, axc) = plt.subplots(1, 3, figsize=(21, 6.5), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.35, 1, 1]))
    phi = np.linspace(phi0 - phip0 * max(X['eps_h'], X['eps_u0']) * 1.3, phi0 + 0.02, 200)
    cols = dict(joint='k', hold=WONG['blue'], unload=WONG['vermillion'])
    for which, ls in (('joint', '-'), ('hold', '--'), ('unload', '--')):
        d = X[which]
        D = d['D_ref'] * np.exp(d['a'] * (phi - phi0))
        axa.plot(phi, D, ls, color=cols[which], lw=2.2 if which == 'joint' else 1.6,
                 label=f"{which} fit: $D_{{ref}}$ = {sig(d['D_ref'])}, a = {d['a']:+.1f} $\\pm$ {d['se_a']:.1f}"
                       + (f"  ($R^2$ hold {d['hold']['R2']:.3f}, re-swelling {d['unload']['R2']:.3f})" if which == 'joint' else f"  ($R^2$ {d[which]['R2']:.3f})"))
        if d['cov'] is not None and which == 'joint':
            g = np.vstack([np.ones_like(phi), phi - phi0])          # d logD / d(logD_ref, a)
            sd = np.sqrt(np.einsum('ij,jk,ik->i', g.T, d['cov'], g.T))
            axa.fill_between(phi, D * np.exp(-sd), D * np.exp(sd), color='k', alpha=0.08)
    ph_h = phi0 - phip0 * X['eps_h']
    ph_u = phi0 - 0.5 * phip0 * X['eps_u0']
    e0 = X['joint']['hold']['eps0']                                   # the fitted initial strain profile of the hold: skins to core
    axa.errorbar([ph_h], [X['Dc_const']['hold']], xerr=[[max(phip0 * (e0.max() - X['eps_h']), 0)], [max(phip0 * (X['eps_h'] - e0.min()), 0)]], fmt='s', color=WONG['blue'], ms=9, capsize=4, lw=1.2,
                 label=f"constant-$D_c$ fit of the hold, {sig(X['Dc_const']['hold'])}, at the held gel's mean $\\phi_s$ (bar: skin to core at the onset)")
    axa.errorbar([ph_u], [X['Dc_const']['unload']], xerr=[[ph_u - (phi0 - phip0 * X['eps_u0'])], [phi0 - ph_u]], fmt='D', color=WONG['vermillion'], ms=9, capsize=4, lw=1.2,
                 label=f"constant-$D_c$ fit of the re-swelling, {sig(X['Dc_const']['unload'])}, at the midpoint of the $\\phi_s$ it sweeps (bar: that range)")
    axa.axvline(phi0, color='0.5', ls=':', lw=1)
    axa.text(phi0, axa.get_ylim()[0] if False else X['joint']['D_ref'] * 1.25, r'$\phi_s^0$', ha='center', fontsize=10, color='0.4')
    axa.set_xlabel(r'$\phi_s = \phi_s^0 - (1-\phi_s^0)\,\varepsilon$'); axa.set_ylabel(r'$D_c$  ($\sigma^2/\tau$)'); axa.set_yscale('log')
    axa.set_title(r'(a) $D_c(\phi_s) = D_{ref}\,e^{a(\phi_s-\phi_s^0)}$ from $\partial_t\varepsilon = \partial_z[D_c(\varepsilon)\,\partial_z\varepsilon]$', fontsize=12)
    axa.yaxis.set_major_formatter(FormatStrFormatter('%.2f')); axa.yaxis.set_minor_formatter(FormatStrFormatter('%.2f'))
    axa.set_xlim(phi.min(), phi.max())                            # the hold's skin bar may run further left: clipped
    axa.grid(alpha=0.3, which='both'); axa.legend(fontsize=7.5, loc='upper left')
    J = X['joint']
    for ax, key, t_, Y, Z, zz_grid, ttl in ((axb, 'hold', X['th'], X['Yh'], None, None, '(b) hold: data vs the joint $D_c(\\phi_s)$ model'),
                                            (axc, 'unload', X['tu'], X['Yu'], X['Zu'], None, '(c) re-swelling: data (COM removed) vs the joint model')):
        shown = np.unique(np.linspace(0, len(t_) - 1, 8).astype(int))
        norm = Normalize(vmin=t_[0], vmax=t_[-1]); cmap = plt.cm.viridis
        for i in shown:
            c = cmap(norm(t_[i]))
            if key == 'hold':
                ax.plot(X['zh'], Y[i], 'o', color=c, ms=3, alpha=0.5)
                ax.plot(X['zh'], J['hold']['U'][i], '-', color=c, lw=1.8)
            else:
                ax.plot(Z[i], Y[i] - np.mean(Y[i]), 'o', color=c, ms=3, alpha=0.5)
                ax.plot(Z[i], J['unload']['U'][i], '-', color=c, lw=1.8)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
        fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02).set_label(r'time since the onset ($\tau$)')
        ax.set_xlabel(r'$\zeta$'); ax.set_ylabel(r'$[u_z(\zeta,t) - u_z(\zeta,0)]/L$'); ax.set_title(ttl, fontsize=12); ax.grid(alpha=0.3)
    fig.suptitle(f"Cooperative diffusivity as a function of the solvent volume fraction, hold and re-swelling fitted together (level _c{L['lvl']})  |  {cfg.sim_name}\n"
                 f"joint: $D_c$ = {sig(J['D_ref'])} at $\\phi_s^0$ = {phi0:.3f}, {sig(J['D_hold'])} at the held state;  separate fits: hold {sig(X['hold']['D_hold'])}, re-swelling {sig(X['unload']['D_ref'])} at $\\phi_s^0$",
                 fontsize=12, fontweight='bold')
    return _save(fig, cfg, 'Dc_phi_fit', L['lvl'])


def fig_Dc(cfg, R, L):
    """Consolidation fit of u_z(zeta,t)/L: (a) data, (b) data + model; (c) the interior strain change vs hold
    time with the windowed and the whole-hold fits; (d) the load-piston stress trace and its own D_c."""
    F = L.get('Dc')
    if F is None:
        print('D_c figure skipped (no fit)')
        return None
    zff = np.linspace(0.0, 1.0, 400)
    dlL = F['DL'] / F['L']
    u0f = F['u_IC'](zff)
    u0_b = F['u_IC'](F['zf'])
    Uf = L.get('unload')
    if Uf is not None:
        fig_Dc_reswelling(cfg, R, L)                       # the PRIMARY D_c figure first (2026-10-15)
        fig = plt.figure(figsize=(18, 19), constrained_layout=True)
        gs = fig.add_gridspec(3, 2)
        axl, axr, axc, axd = (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]))
        axe = fig.add_subplot(gs[2, :])
        _unload_trace_panel(cfg, axe, F, Uf, L.get('Dc_stress'))
    else:
        fig, ((axl, axr), (axc, axd)) = plt.subplots(2, 2, figsize=(18, 13), constrained_layout=True)
    _dc_time_panels(cfg, axc, axd, F, L.get('Dc_stress'))
    norm = Normalize(vmin=F['ts'][F['early'][0]], vmax=F['ts'][F['early'][-1]])
    cmap = plt.cm.viridis
    for i in F['shown']:
        c = cmap(norm(F['ts'][i]))
        axl.plot(F['zf'], u0_b + F['uhat'][i][F['idx']], 'o-', color=c, ms=3, alpha=0.6)
        axr.plot(F['zf'], u0_b + F['uhat'][i][F['idx']], 'o', color=c, ms=3, alpha=0.35)
        axr.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=2.0)
    fs = F['f_sup']                       # support's share of the closure (0 top-only, 1/2 symmetric)
    for ax in (axl, axr):
        ax._tri_has_colorbar = True
        ax.plot([0, 1], [fs * dlL, (fs - 1.0) * dlL], 'k:', lw=1.8,
                label=r'affine ($t\to\infty$):  $(\Delta L/L)(f_\mathrm{sup}-\zeta)$')
        ax.plot(zff, u0f, '-', color='0.45', lw=1.6, label=r'IC: fitted hold-onset state')
        ax.plot(0, fs * dlL, 's', color=WONG['blue'], ms=9, zorder=5,
                label=r'BC: $u_z(0,t)=+\Delta L_\mathrm{sup}/L$ (held support)')
        ax.plot(1, (fs - 1.0) * dlL, 'D', color=WONG['vermillion'], ms=9, zorder=5,
                label=r'BC: $u_z(1,t)=-\Delta L_\mathrm{pist}/L$ (held piston)')
        ax.set(xlabel=flip_hint(R, r'$\zeta=(z-z_\mathrm{perm})/L$', nl=True), ylabel=r'$u_z/L$', xlim=_zlim(R))
        ax.grid(alpha=0.3)
        smart_legend(ax, fontsize=11)
    axl.set_title(r'(a) $u_z(\zeta,t)/L$ -- fitted hold snapshots (data + fitted IC offset)', fontsize=15)
    lbl = '' if cfg.DC_FREE_AMPS else rf"$\beta=p_0/M={sig(F['beta'])}$, "
    axr.set_title(rf"(b) consolidation fit, first {F['n_tau1_fit']:.1f} $\tau_1$: $D_c={sig(F['Dc'])}\ \sigma^2/\tau$, " + lbl + rf"$R^2={sig(F['R2'])}$", fontsize=15)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=[axl, axr], fraction=0.015, pad=0.04).set_label('timestep')
    fig.suptitle(f'Cooperative diffusivity, LOAD measurement: hold consolidation fit (level _c{L["lvl"]}'
                 + ('; the UNLOAD fit above is the trusted estimate' if Uf is not None else '; no unload in this run') + f')  |  {cfg.sim_name}  |  '
                 f'plate-gap closure $\\Delta L={sig(F["DL"])}\\,\\sigma$: support {fs:.0%} / piston {1 - fs:.0%}'
                 + ('  (symmetric drive)' if abs(fs - 0.5) < 0.05 else '  (top-only drive)' if fs < 0.05 else ''),
                 fontsize=12, fontweight='bold')
    return _save(fig, cfg, 'Dc_consolidation_fit', L['lvl'])


def print_hold_check(cfg, levels):
    """7b of the original notebook: was each level held long enough?"""
    hd = [L for L in levels if L.get('Dc') is not None]
    if not hd:
        return
    print(f'\nHOLD-ADEQUACY CHECK  (fit: tau_1 = L^2/(4 pi^2 D_c), held slab; slow: the deck\'s sizing formula '
          f'L^2/(pi^2 Dc_est); residual = mean excess stress over the last '
          f'{cfg.plateau_frac:.0%} of the hold, i.e. the window M is read from)')
    for L in hd:
        F = L['Dc']
        print(f"  level _c{L['lvl']}:  L = {F['L']:.1f} sigma   hold T = {F['hold_T']:.0f} tau "
              f"= {F['hold_T'] / cfg.dt_lj / 1e6:.2f}M steps   (D_c fitted on the first {F['n_tau1_fit']:.1f} tau_1 = "
              f"{len(F['early'])} of {len(F['early_all'])} snapshots; whole hold {F['Dc_all']:.3e}"
              + (f"; stress trace {L['Dc_stress']['Dc']:.3e}" if (L.get('Dc_stress') or {}).get('ok') else '')
              + (f"; unload trace {L['unload']['trace']['Dc']:.3e}" if (L.get('unload') or {}).get('trace') else '') + ')')
        for tag, h in F['hold_check'].items():
            flag = '' if h['ok'] else '   <-- TOO SHORT'
            print(f"     {tag:<4s} D_c={h['Dc']:.3f}:  tau_1 = {h['tau1']:.0f} tau = {h['tau1'] / cfg.dt_lj / 1e6:.2f}M steps"
                  f"  |  held {F['hold_T'] / h['tau1']:.2f} tau_1  ->  residual {h['avg'] * 100:5.2f}% "
                  f"(end {h['end'] * 100:5.2f}%)  |  {cfg.DC_TARGET_RESID:.0%} needs {h['need']:.1f} tau_1 "
                  f"= {h['need'] * h['tau1'] / cfg.dt_lj / 1e6:.1f}M steps{flag}")
        if 'M_pist' in L:
            mixed = (L.get('M_pist_ref') == 'absent')
            mp, mn = (L['M_pist_abs'], L['M_net_abs']) if mixed else (L['M_pist'], L['M_net'])
            print(f"     observed M_piston/M_network - 1 = {mp / mn - 1:+.1%}   vs   "
                  f"predicted unrelaxed excess (slow D_c) = {F['hold_check']['slow']['avg']:+.1%}"
                  + ('   [both ABSOLUTE here: no piston_force_avg_ref, so the increment ratio is not like-for-like]' if mixed else ''))
    print(f'  (triaxial_compression.lmp sizes each hold as n_tau_hold * L^2/(pi^2 Dc_est) from the live compressed '
          f'BB thickness and Dc_est = {cfg.DC_SLOW_REF:.2f}; the held slab relaxes as L^2/(4 pi^2 D_c), so that is the hold of '
          f'D_c = {cfg.DC_SLOW_REF / 4:.4f} -- adequate while the fitted D_c stays above it, but the fitted D_c falls with strain '
          f'faster than the (1-eps)^2 sizing assumes, so the higher levels are held too short: compare the observed '
          f'piston/network excess above.)')


def fig_Dc_sweep_load(cfg, R, levels):
    """(a) D_c vs applied strain; then one consolidation-fit panel per level."""
    hd = [L for L in levels if L.get('Dc') is not None]
    if not hd:
        print('D_c figure skipped (no level produced a fit)')
        return None
    n = len(hd)
    fig, axes = plt.subplots(1, n + 1, figsize=(7 + 6.5 * n, 6), constrained_layout=True)
    hu = [L for L in hd if (L.get('unload') or {}).get('prof')]
    fig.suptitle(f'Cooperative diffusivity across the sweep  |  {cfg.sim_name}  |  LOAD measurements (hold consolidation, stress trace)'
                 + ('; stars = UNLOAD (free re-swelling), the trusted estimate' if hu else '; this run has no unload, so no trusted (unload) D_c'),
                 fontsize=13, fontweight='bold')
    ax = axes[0]
    e = np.array([L['eps'] for L in hd])
    v = np.array([L['Dc']['Dc'] for L in hd])
    ax.plot(e, v, 'o-', color=WONG['reddishpurple'], lw=2, ms=8,
            label=rf"LOAD: hold profile fit, first {cfg.DC_FIT_TAU1:g} $\tau_1$ of the hold" if cfg.DC_FIT_TAU1 > 0 else 'LOAD: hold profile fit')
    if hu:
        ax.plot([L['eps'] for L in hu], [L['unload']['prof']['Dc'] for L in hu], '*', color=WONG['blue'], ms=17, zorder=6,
                label='UNLOAD: free re-swelling profile fit (the trusted $D_c$)')
    for L in hd:   # R^2 label under each point, clipped inside the axes
        ax.annotate(f"$R^2$={sig(L['Dc']['R2'])}", (L['eps'], L['Dc']['Dc']), textcoords='offset points',
                    xytext=(0, -12), ha='center', va='top', fontsize=10, color='0.35',
                    annotation_clip=True)
    if cfg.DC_FIT_TAU1 > 0:            # 2026-10-09: the whole-hold fit and the stress-trace D_c alongside
        ax.plot(e, [L['Dc']['Dc_all'] for L in hd], 'o:', mfc='none', color=WONG['reddishpurple'], lw=1.5, ms=8,
                label='LOAD: hold profile fit, whole hold')
    hs = [L for L in hd if (L.get('Dc_stress') or {}).get('ok')]
    if hs:
        ax.errorbar([L['eps'] for L in hs], [L['Dc_stress']['Dc'] for L in hs], yerr=[L['Dc_stress']['Dc_se'] for L in hs],
                    fmt='^--', color=WONG['green'], lw=1.8, ms=9, capsize=5, label='LOAD: load-piston stress trace (held-slab series)')
    ax.margins(x=0.12, y=0.15)
    ax.axhline(cfg.DC_SLOW_REF / 4.0, color='0.4', ls=':', lw=1.5,
               label=f"deck hold-sizing $D_c$ = {sig(cfg.DC_SLOW_REF / 4)}  ($Dc_{{est}}$ = {sig(cfg.DC_SLOW_REF)} in the deck's $L^2/\\pi^2 D_c$ formula)")
    ax.set_xlabel(r'applied strain  $\varepsilon$')
    ax.set_ylabel(r'$D_c$  ($\sigma^2/\tau$)')
    ax.set_title(r'(a) $D_c$ vs applied strain: LOAD (hold) fits' + (', UNLOAD stars' if hu else ''), fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=11)
    zff = np.linspace(0, 1, 300)
    for k, L in enumerate(hd):
        ax = axes[k + 1]
        F = L['Dc']
        i_lvl = levels.index(L)
        cmap = plt.cm.viridis
        norm = Normalize(vmin=F['ts'][F['early'][0]], vmax=F['ts'][F['early'][-1]])
        u0_b = F['u_IC'](F['zf'])
        for i in F['shown']:
            c = cmap(norm(F['ts'][i]))
            ax.plot(F['zf'], u0_b + F['uhat'][i][F['idx']], 'o', color=c, ms=2.5, alpha=0.35)
            ax.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=1.5)
        dlL, fs = F['DL'] / F['L'], F['f_sup']
        ax.plot([0, 1], [fs * dlL, (fs - 1.0) * dlL], 'k:', lw=1.6)
        ax.set(xlabel=flip_hint(R, r'$\zeta$'), ylabel=r'$u_z/L$', xlim=_zlim(R))
        ax.set_title(fr"({'bcdefgh'[k]}) $\varepsilon={L['lvl']}$, LOAD (hold):  $D_c={sig(F['Dc'])}$, $R^2={sig(F['R2'])}$",
                     fontsize=14, color=level_color(i_lvl))
        ax.grid(alpha=0.3)
    return _save(fig, cfg, 'sweep_Dc')


def add_load_Dc(cfg, R, L, verbose=True):
    """The hold fits load_level made until 2026-11-17, on a level L it has already loaded: L['Dc'] (fit_Dc),
    L['kappa_hold'], L['Dc_stress'] (fit_Dc_stress), L['kappa_stress'].  Returns L['Dc'] (None when there is no fit)."""
    say = print if verbose else (lambda *a, **k: None)
    L['Dc'] = F = fit_Dc(cfg, R, load_disp(cfg, R, L['lvl']))
    L['Dc_stress'] = None
    if F is None:
        say('  NOTE: no hold D_c (missing disp_z_polymer / piston_position, no support plane, or too few bins)')
        return None
    say(f"  LOAD D_c (hold profile fit) = {F['Dc']:.4e} sigma^2/tau  (R^2 = {F['R2']:.3f};  L = {F['L']:.2f}, hold = {F['hold_T'] / F['tau1']:.1f} tau_1; "
        f"first {F['n_tau1_fit']:.1f} tau_1 = {len(F['early'])} of {len(F['early_all'])} snapshots);  whole hold {F['Dc_all']:.4e}")
    L['kappa_hold'] = _kappa_from(L, F['Dc'])
    if cfg.DC_STRESS and 'pfa_P' in L:
        S = L['Dc_stress'] = fit_Dc_stress(cfg, L['pfa_step'], L['pfa_P'], F['L'], F['t_hold'])
        if S['ok']:
            say(f"  LOAD D_c (stress trace) = {S['Dc']:.4e} +/- {S['Dc_se']:.1e}")
            L['kappa_stress'] = {key: S['Dc'] / L[Mk] for key, Mk in (('net', 'M_net'), ('pist', 'M_pist')) if Mk in L}
        else:
            say(f"  LOAD D_c (stress trace): not resolved -- {S['why']}")
    return F


# ---------------------------------------------------------------------------
#  For reference: the D_c block of load_level as it stood before 2026-11-17 (the unload part is still in
#  lib/triaxial.py; add_load_Dc above is the hold part of it).
# ---------------------------------------------------------------------------
    #     # ---- D_c and kappa ---------------------------------------------------
    #     L['Dc'] = fit_Dc(cfg, R, load_disp(cfg, R, lvl))
    #     if L['Dc'] is None:
    #         say('  NOTE: no D_c (missing disp_z_polymer / piston_position, no support plane, or too few bins)')
    #     else:
    #         F = L['Dc']
    #         say(f"  D_c = {F['Dc']:.4e} sigma^2/tau  (R^2 = {F['R2']:.3f};  L = {F['L']:.2f}, "
    #             f"DL/L = {F['DL'] / F['L']:.4f} [support share {F['f_sup']:.2f}], hold = {F['hold_T']:.0f} tau = {F['hold_T'] / F['tau1']:.1f} tau_1; "
    #             f"fitted: first {F['n_tau1_fit']:.1f} tau_1 = {len(F['early'])} of {len(F['early_all'])} {'FINE' if F['fine'] else 'coarse'} snapshots, "
    #             f"first at {F['t_lj'][F['early'][0]]:.0f} tau, block {F['W']:.0f} tau)")
    #         if F['fit_tau1'] > 0 and len(F['early']) < len(F['early_all']):
    #             say(f"      whole hold: D_c = {F['Dc_all']:.4e} (R^2 = {F['R2_all']:.3f})"
    #                 + (f"   <-- {F['Dc'] / F['Dc_all']:.2f}x the windowed value: a slower process after the poroelastic transient "
    #                    f"(fig_Dc panel c); kappa below uses the windowed D_c" if F['stab_flag'] else '   (consistent)'))
    #         if F.get('stab'):
    #             say('      stability: ' + '   '.join(f"{s['what']} ({s['n']} snaps): {s['Dc']:.3e}" for s in F['stab']))
    #         L['kappa_hold'] = _kappa_from(L, F['Dc'])
    #         say('  kappa = D_c(hold)/M:  ' + '   '.join(f"{k}: {v['k']:.4e} [{v['lo']:.4e}, {v['hi']:.4e}]"
    #                                                    for k, v in L['kappa_hold'].items()))
    #         # ---- the stress trace: an independent D_c (2026-10-09) ----
    #         L['Dc_stress'] = None
    #         if cfg.DC_STRESS and 'pfa_P' in L:
    #             S = fit_Dc_stress(cfg, L['pfa_step'], L['pfa_P'], F['L'], F['t_hold'])
    #             L['Dc_stress'] = S
    #             if S['ok']:
    #                 say(f"  D_c (stress trace) = {S['Dc']:.4e} +/- {S['Dc_se']:.1e}  (held-slab series from {S['t_min']:.0f} tau, "
    #                     f"{S['n']} blocks of {S['block']:.0f} tau, skin {S['delta']:.1f} sigma, rms {S['rms']:.4f}; tau_1 = {S['tau1']:.0f} tau)"
    #                     + (f";  tail exponential from {S['t_tail']:.0f} tau: tau = {S['tau_tail']:.0f} +/- {S['tau_tail_se']:.0f} "
    #                        f"-> D_c = {S['Dc_tail']:.4e}" if S['tail_ok'] else ';  tail exponential not credible')
    #                     + f"   [profile/stress = {F['Dc'] / S['Dc']:.2f}]")
    #                 L['kappa_stress'] = {key: S['Dc'] / L[Mk] for key, Mk in (('net', 'M_net'), ('pist', 'M_pist')) if Mk in L}
    #                 say('      kappa from the stress-trace D_c:  ' + '   '.join(f"{k}: {v:.4e}" for k, v in L['kappa_stress'].items()))
    #             else:
    #                 say(f"  D_c (stress trace): not resolved -- {S['why']}")
    #         # ---- the unload / free re-swelling (deck Phase 2c), the skin-free control (2026-10-09) ----
    #         L['unload'] = None
    #         if cfg.DC_UNLOAD:
    #             Uf = fit_Dc_unload(cfg, R, load_unload(cfg, R, lvl, F), F)
    #             L['unload'] = Uf
    #             if Uf is not None:
    #                 T, Pf = Uf['trace'], Uf['prof']
    #                 say(f"  UNLOAD (free re-swelling after the hold; retraction {Uf['retract_T']:.0f} tau; L = {Uf['L']:.2f} = {cfg.UNLOAD_L!r} "
    #                     f"[hold {Uf['L_hold']:.2f}, reference {Uf['L_ref']:.2f}: {Uf['L_ref_src']}]; tau_1 = L^2/(pi^2 D_c), faces free)")
    #                 for T_, nm in ((T, 'bounding-box thickness trace'), (Uf.get('trace_rg'), 'Rg thickness trace')):
    #                     if T_ is not None:
    #                         say(f"    D_c ({nm}, dL_inf fixed = {T_['dL_inf']:.2f} sigma) = {T_['Dc']:.4e} +/- {T_['Dc_se']:.1e}  "
    #                             f"(R^2 {T_['R2']:.3f}, rms {T_['rms']:.3f} sigma; tau_1 = {T_['tau1']:.0f} tau, {T_['t'][-1] / T_['tau1']:.1f} tau_1 recorded, "
    #                             f"{100 * T_['relaxed']:.1f}% relaxed);  free asymptote: D_c = {T_['Dc_free']:.4e} +/- {T_['Dc_free_se']:.1e}, "
    #                             f"dL_inf = {T_['dL_inf_free']:.2f} +/- {T_['dL_inf_free_se']:.2f}")
    #                 if T is not None and Pf is not None and T['Dc'] > 1.3 * Pf['Dc']:
    #                     say("    (the bounding-box trace reads high: L_bb follows the outermost beads, which spring back within ~100 tau of the "
    #                         "retraction -- ~3 sigma of the 12 before the bulk has moved -- so it is not a bulk re-swelling measure; the profile fit is)")
    #                 if Pf is not None:
    #                     say(f"    D_c (profiles, first {Pf['n_tau1_fit']:.1f} tau_1 = {len(Pf['early'])} of {Pf['n_all']} {'FINE' if Uf['fine'] else 'coarse'} snapshots) = "
    #                         f"{Pf['Dc']:.4e} (R^2 {Pf['R2']:.3f});  whole unload {Pf['Dc_all']:.4e}"
    #                         + ('   <-- differ' if Pf['stab_flag'] else '') + ('   stability: ' + '  '.join(f"{x['what']} ({x['n']}): {x['Dc']:.3e}" for x in Pf['stab']) if Pf['stab'] else ''))
    #                 ref = Pf['Dc'] if Pf is not None else T['Dc']
    #                 say(f"    vs the hold:  windowed profile {F['Dc'] / ref:.2f}x,  whole hold {F['Dc_all'] / ref:.2f}x"
    #                     + (f",  stress trace {L['Dc_stress']['Dc'] / ref:.2f}x" if (L.get('Dc_stress') or {}).get('ok') else '')
    #                     + "  the re-swelling profile D_c  (1.0 = the hold's transient is the same D_c; < 1 = the hold reads slower)")
    #         # ---- the PRIMARY D_c (DC_PRIMARY, 2026-10-15) and kappa = D_c/M from it ----
    #         L['Dc_primary'] = primary_Dc(cfg, L)
    #         Dp = L['Dc_primary']
    #         L['kappa'] = _kappa_from(L, Dp['Dc'])
    #         say(f"  PRIMARY D_c = {Dp['Dc']:.4e} ({Dp['src']});  kappa = D_c/M:  "
    #             + '   '.join(f"{k}: {v['k']:.4e} [{v['lo']:.4e}, {v['hi']:.4e}]" for k, v in L['kappa'].items()))
