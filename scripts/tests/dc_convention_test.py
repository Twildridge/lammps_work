"""One consolidation equation, one D_c, both D_c fitters (2026-10-03).

The compression fit (tri.fit_Dc: slab HELD between two drained plates) and the permeation fit
(tri.fit_perm_Dc: one face pinned on the support, one free, pressure drop prescribed) were each
checked on synthetic data built from their OWN mode series.  Here neither series is used: a
finite-difference solver integrates

    du/dt = q(t) + D_c d2u/dz2        (u: network displacement, q: total flux, uniform in z)

with the PHYSICAL boundary conditions of each experiment, and both fitters read the result.

    held slab    u(0) = a(t), u(L) = b(t) (plates: ramp, then held);  u'(0) = u'(L)
                 (same bath pressure at both faces -> same network stress, total stress uniform)
    permeation   u(0) = 0 (support);  u'(L) = 0 (free feed face, sigma' = 0);
                 u'(0) = -dP(t)/M (the support carries the whole pressure drop from t = 0+)
    unload       u'(0) = u'(L) = s (both faces free after the plates retract: the network stress
                 there is zero, i.e. the strain returns to the reference value, s = dL_inf/L above
                 the uniform compressed state);  mean u = 0 (the gel's COM stays)   [2026-10-09]

q(t) is the extra unknown that the third condition fixes.  Eigenmodes in ONE convention
(L = the FULL thickness between the faces in both problems, never a half-thickness):

    held slab    strain modes cos/sin(2 m pi z/L)   rate 4 m^2 pi^2 D_c/L^2   tau_1 = L^2/(4 pi^2 D_c)
    permeation   strain modes sin(k pi z/L)         rate   k^2 pi^2 D_c/L^2   tau_1 = L^2/(  pi^2 D_c)
    unload       strain modes sin(k pi z/L)         rate   k^2 pi^2 D_c/L^2   tau_1 = L^2/(  pi^2 D_c)
                 (the SAME class as the permeation transient: strain pinned at the faces, solvent crossing
                 them; the hold is the other class -- no flux through the plates, strain free at the faces)

The factor 4 between the two tau_1 is physical (drainage path L/2 against L), not a convention:
both fits return the D_c of the equation above, and D_c = kappa M in both.  The test asserts that
each fitter recovers the solver's D_c and that the solver's steady permeation flux is Darcy's
(D_c/M) dP/L.  The snapshots are built the way the decks write them (block averages over the whole
output interval, binned by CURRENT z), so the test also covers the fits' handling of both.

Usage: python dc_convention_test.py [D_c]
"""
import sys, io, contextlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
import triaxial as tri

DT_LJ = 0.005            # tau per MD step (the decks' timestep_prod)
SAMPLE = 5000            # MD steps between stored solver states (the thickness-trace cadence)
BIN = 2.0                # sigma, the decks' z-bin
N_BIN = 300              # polymer beads per populated bin (> Config.Ncount_min)


def solve(D, L, n, steps, mode, drive, dt=0.5):
    """Backward-Euler integration of du/dt = q + D u'' on n + 1 nodes over `steps` MD steps.
    mode 'held': drive(t) -> (u(0), u(L));  mode 'perm': drive(t) -> dP(t)/M.
    Returns (Z, sample steps, u[sample, node], q[sample]); the state is stored every SAMPLE steps."""
    h = L / n
    A = np.zeros((n + 2, n + 2))                     # unknowns u_0 .. u_n, q
    for i in range(1, n):
        A[i, i - 1:i + 2] = np.array([-1.0, 2.0, -1.0]) * D / h ** 2
        A[i, i] += 1.0 / dt
        A[i, n + 1] = -1.0
    lo = np.array([-3.0, 4.0, -1.0]) / (2 * h)       # one-sided u'(0)
    hi = np.array([1.0, -4.0, 3.0]) / (2 * h)        # one-sided u'(L)
    if mode == 'free':                               # unload: u'(0) = u'(L) = s(t), mean u = 0 fixes q
        A[0, 0:3] = lo
        A[n, n - 2:n + 1] = hi
        A[n + 1, 0:n + 1] = 1.0 / (n + 1)
    else:
        A[0, 0] = 1.0                                # u(0) prescribed in both problems
    if mode == 'held':
        A[n, n] = 1.0                                # u(L) prescribed
        A[n + 1, 0:3] = lo                           # u'(0) - u'(L) = 0
        A[n + 1, n - 2:n + 1] -= hi
    elif mode == 'perm':
        A[n, n - 2:n + 1] = hi                       # u'(L) = 0
        A[n + 1, 0:3] = lo                           # u'(0) = -dP(t)/M
    Ai = np.linalg.inv(A)
    per = int(round(SAMPLE * DT_LJ / dt))
    assert abs(per * dt - SAMPLE * DT_LJ) < 1e-9
    x = np.zeros(n + 2)
    U, Q, S = [x[:n + 1].copy()], [0.0], [0]
    for k in range(1, int(steps // SAMPLE) * per + 1):
        t = k * dt
        rhs = np.zeros(n + 2)
        rhs[1:n] = x[1:n] / dt
        if mode == 'held':
            rhs[0], rhs[n] = drive(t)
        elif mode == 'free':
            rhs[0] = rhs[n] = drive(t)
        else:
            rhs[n + 1] = -drive(t)
        x = Ai @ rhs
        if k % per == 0:
            U.append(x[:n + 1].copy()); Q.append(x[n + 1]); S.append(k // per * SAMPLE)
    return np.linspace(0.0, L, n + 1), np.array(S), np.array(U), np.array(Q)


def eulerian_bins(zc, z_face, Z, u, shift=0.0):
    """What fix ave/chunk writes for one configuration: the displacement `u + shift` of the
    material points Z, read at the bin centres zc of their CURRENT position z_face + Z + u + shift.
    Returns (uz[bin], populated[bin])."""
    z_now = z_face + Z + u + shift
    inside = (zc > z_now[0]) & (zc < z_now[-1])
    return np.where(inside, np.interp(zc, z_now, u + shift), 0.0), inside


def window_snapshots(zc, z_face, Z, S, U, ends, stride, shift=None, ref=None):
    """Block averages over the windows (end - stride, end] of the stored states, as the decks'
    ave/chunk (nevery x nrepeat = nfreq) writes them; `ref` is subtracted (the displace/atom reset)."""
    uz, Nc = [], []
    for e in ends:
        sel = np.where((S > e - stride) & (S <= e))[0]
        prof = [eulerian_bins(zc, z_face, Z, U[i] - (0.0 if ref is None else ref), 0.0 if shift is None else shift[i]) for i in sel]
        pop = np.all([p[1] for p in prof], axis=0)
        uz.append(np.where(pop, np.mean([p[0] for p in prof], axis=0), 0.0))
        Nc.append(np.where(pop, N_BIN, 0))
    return np.array(uz), np.array(Nc)


def solver_check(D, L=20.0, dP=0.1, M=0.42):
    """The solver against the two closed-form steady states of the permeation problem (a thin slab
    run to 12 tau_1): u_F = -dP L/(2M) and the Darcy flux |q| = (D/M) dP/L.  Returns |q|/Darcy."""
    steps = int(np.ceil(12 * L ** 2 / (np.pi ** 2 * D) / DT_LJ / SAMPLE)) * SAMPLE
    Z, S, U, Q = solve(D, L, 100, steps, 'perm', lambda t: dP / M)
    assert abs(U[-1][-1] / (-dP * L / (2 * M)) - 1.0) < 0.01, 'steady u_F is not -dP L/(2M)'
    r = float(-Q[-1] / ((D / M) * dP / L))
    assert abs(r - 1.0) < 0.01, 'steady flux is not Darcy with kappa = D_c/M'
    return r


def _cfg(**kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return tri.Config(DATANAME='dc_convention', INTERACTION='1.0_1.0', RUN_ID='dc_convention', two_pist=True, **kw)


def held_slab(D, split=0.5, L=113.6, closure=12.8, v=0.04, hold=8_400_000, n_snap=20, noise=0.01, seed=1):
    """The compression deck's level: plates close `closure` sigma at v (the load piston's share
    `split` coming DOWN, the rest the support going UP), then a hold of `hold` steps read as
    n_snap block-averaged, hold-referenced profiles.  Returns (fit, tau_1 of the equation)."""
    rng = np.random.default_rng(seed)
    t_ramp = closure / v
    ramp_steps = int(np.ceil(t_ramp / DT_LJ / SAMPLE)) * SAMPLE
    drive = lambda t: ((1.0 - split) * closure * min(t / t_ramp, 1.0), -split * closure * min(t / t_ramp, 1.0))
    Z, S, U, _ = solve(D, L, 400, ramp_steps + hold, 'held', drive)
    assert np.allclose(U[-1], np.interp(Z, [0, L], drive(1e9)), atol=2e-3 * closure), 'hold not relaxed to uniform strain'
    i_hold = int(np.argmax(S >= ramp_steps))
    stride = hold // n_snap
    ends = ramp_steps + stride * np.arange(1, n_snap + 1)
    z_sup0, gap = 30.0, 1.0
    z_sup = z_sup0 + (1.0 - split) * closure                  # held support plane
    z_face = z_sup + gap                                      # pinned polymer face
    zc = np.arange(0.5 * BIN, z_face + L + 20.0, BIN)
    # the fit's frame is the held configuration, so place the material points there (u(t_hold) removed)
    uz, Nc = window_snapshots(zc, z_face, Z, S, U, ends, stride, ref=U[i_hold])
    uz = uz + np.where(Nc > 0, rng.normal(0.0, noise, uz.shape), 0.0)
    cfg = _cfg(mode='compression', COMP_LEVELS=['0.10'], DC_N_MODES=5, DC_TRIM_BINS=2)
    R = dict(z_support=z_sup0, z_piston=z_sup0 + 2 * gap + L + closure)
    disp = dict(ts=ends.astype(float), z=zc, Nc=Nc, uz=uz, z_pist_held=z_face + L + gap, z_supp_held=z_sup,
                t_hold=float(ramp_steps), L_bb=L)
    return tri.fit_Dc(cfg, R, disp), L ** 2 / (4 * np.pi ** 2 * D)


def permeation(D, rigid_drop=3.0, L0=126.4, dP=0.1, M=0.42, onset=1_000_000, n_ramp=100, ramp_steps=10_000,
               prod=9_800_000, stride=50_000, noise=0.01, seed=2):
    """The permeation deck's run: dP ramped in n_ramp stairs of ramp_steps from `onset`, then held;
    disp_z_polymer every `stride` steps (block-averaged, reset at the ramp start) and the thickness
    trace every SAMPLE steps.  rigid_drop > 0: the network starts that far above the contact plane and
    falls onto it during the ramp (the u_0 that PERM_RIGID = 'contact' removes).
    Returns (fit, tau_1 of the equation)."""
    rng = np.random.default_rng(seed)
    t_on = onset * DT_LJ
    stair = ramp_steps * DT_LJ
    s_full = dP / M
    drive = lambda t: s_full * min(max(np.floor((t - t_on) / stair) + 1.0, 0.0), n_ramp) / n_ramp
    Z, S, U, _ = solve(D, L0, 400, onset + n_ramp * ramp_steps + prod, 'perm', drive)
    z_sup, gap = 51.0, 1.0
    z_P = z_sup + gap                                         # contact plane (R['z_support'] + PERM_GAP)
    u0 = -rigid_drop * np.clip((S - onset) / (n_ramp * ramp_steps), 0.0, 1.0)
    zc = np.arange(0.5 * BIN, z_P + rigid_drop + L0 + 20.0, BIN)
    ends = onset + stride * np.arange(1, (S[-1] - onset) // stride + 1)
    uz, Nc = window_snapshots(zc, z_P + rigid_drop, Z, S, U, ends, stride, shift=u0)
    uz = uz + np.where(Nc > 0, rng.normal(0.0, noise, uz.shape), 0.0)
    tj = (np.arange(n_ramp) * ramp_steps) * DT_LJ
    fz = dict(t0=float(onset), t_end=float(onset + (n_ramp - 1) * ramp_steps), full=dP, n=n_ramp, tj=tj,
              wj=np.full(n_ramp, 1.0 / n_ramp), source='synthetic')
    disp = dict(ts=ends.astype(float), z=zc, Nc=Nc, uz=uz, stride=float(stride), forcing=fz, reset_mode='ramp_start',
                t_reset=float(onset), bb_step=S.astype(float), bb_L=L0 + U[:, -1] + rng.normal(0.0, noise, len(S)))
    cfg = _cfg(mode='permeation', PERM_GAP=gap, PERM_RIGID='contact' if rigid_drop > 0 else 'none',
               PERM_DC_N_MODES=8, PERM_TRACE_SKIP=1_500_000)
    P = dict(halt_ts=int(S[-1] - 0.25 * (S[-1] - onset)))
    F = tri.fit_perm_Dc(cfg, dict(z_support=z_sup), P, disp)
    return F, L0 ** 2 / (np.pi ** 2 * D)


def unload(D, L_hold=113.6, s=0.11, retract=60_000, hold=25_000_000, stride=25_000, noise=0.01, seed=3):
    """The unload of the fine batch: plates retract over `retract` steps (the faces are free from the start,
    they detach within a few tau), then a re-swelling hold of `hold` steps; the thickness trace every SAMPLE
    steps from the retraction start, the displacement profile every `stride` steps reset when the plates
    stop.  Returns (fit, tau_1 of the equation)."""
    rng = np.random.default_rng(seed)
    Z, S, U, _ = solve(D, L_hold, 400, retract + hold, 'free', lambda t: s)
    i_free = int(np.argmax(S >= retract))
    assert abs((U[-1][-1] - U[-1][0]) / (s * L_hold) - 1.0) < 0.01, 'unload not relaxed to the reference thickness'
    z_bot = 31.0                                              # bottom face of the compressed gel at the end of the hold
    zc = np.arange(0.5 * BIN, z_bot + L_hold * (1.0 + s) + 20.0, BIN)
    ends = retract + stride * np.arange(1, hold // stride + 1)
    uz, Nc = window_snapshots(zc, z_bot, Z, S, U, ends, stride, ref=U[i_free])
    uz = uz + np.where(Nc > 0, rng.normal(0.0, noise, uz.shape), 0.0)
    Ud = dict(fine=True, ts=ends.astype(float), z=zc, Nc=Nc, uz=uz, t_retract=0.0, t_free=float(retract), z_pist_free=np.nan,
              bb_step=S.astype(float), bb_L=L_hold + U[:, -1] - U[:, 0] + rng.normal(0.0, 3 * noise, len(S)),
              L_ref=L_hold * (1.0 + s), L_ref_src='synthetic', L_hold=L_hold, Z_bot=z_bot, Z_top=z_bot + L_hold)
    cfg = _cfg(mode='compression', COMP_LEVELS=['0.10'], DC_N_MODES=5, DC_TRIM_BINS=2)
    F_hold = dict(L=L_hold, z_perm=z_bot, z_feed=z_bot + L_hold, Dc=D, Dc_all=D)
    return tri.fit_Dc_unload(cfg, dict(z_support=z_bot - 1.0), Ud, F_hold), L_hold ** 2 / (np.pi ** 2 * D)


def main(D=0.10, tol=0.10):
    print(f'consolidation equation du/dt = q(t) + D_c u\'\' solved by finite differences with D_c = {D} sigma^2/tau')
    print(f'solver check (permeation steady state): u_F = -dP L/(2M) within 1%, |q|/[(D_c/M) dP/L] = {solver_check(D):.4f}\n')
    out = {}
    for split in (0.5, 1.0):
        F, tau1 = held_slab(D, split)
        n = len(F['kk'])
        fam = 'sin (strain symmetric)' if np.abs(F['A'][:n]).max() > np.abs(F['A'][n:]).max() else 'cos - 1 (strain antisymmetric)'
        out[f'held slab, drive_split {split}'] = F['Dc']
        print(f"held slab (fit_Dc), drive_split {split}:  D_c = {F['Dc']:.4f}  ({F['Dc'] / D:.3f} x true, R^2 {F['R2']:.4f});  "
              f"L = {F['L']:.1f} (full thickness), tau_1 = L^2/(4 pi^2 D_c) = {tau1:.0f} tau;  largest amplitude in the {fam} block")
    for drop in (0.0, 3.0):
        F, tau1 = permeation(D, drop)
        T = F['trace']
        out[f'permeation profile, rigid drop {drop}'] = F['Dc']
        out[f'permeation trace, rigid drop {drop}'] = T['Dc']
        print(f"permeation (fit_perm_Dc), rigid drop {drop} sigma (PERM_RIGID = {F['rigid']}):  profile D_c = {F['Dc']:.4f}  "
              f"({F['Dc'] / D:.3f} x true, R^2 {F['R2']:.4f});  trace D_c = {T['Dc']:.4f}  ({T['Dc'] / D:.3f} x true, R^2 {T['R2']:.4f});  "
              f"L_0 = {F['L0']:.1f} (full thickness), tau_1 = L_0^2/(pi^2 D_c) = {tau1:.0f} tau")
    Fu, tau1 = unload(D)
    T, Pf = Fu['trace'], Fu['prof']
    out['unload trace (dL_inf fixed)'] = T['Dc']
    out['unload trace (dL_inf free)'] = T['Dc_free']
    out['unload profile'] = Pf['Dc']
    print(f"unload (fit_Dc_unload):  thickness trace D_c = {T['Dc']:.4f} ({T['Dc'] / D:.3f} x true, R^2 {T['R2']:.4f}; dL_inf free: "
          f"{T['Dc_free']:.4f}, dL_inf {T['dL_inf_free']:.2f} vs {Fu['dL_inf_ref']:.2f});  profile D_c = {Pf['Dc']:.4f} "
          f"({Pf['Dc'] / D:.3f} x true, R^2 {Pf['R2']:.4f}, first {Pf['n_tau1_fit']:.1f} tau_1; whole {Pf['Dc_all']:.4f});  "
          f"L = {Fu['L']:.1f} (compressed thickness), tau_1 = L^2/(pi^2 D_c) = {tau1:.0f} tau")
    bad = {k: v for k, v in out.items() if abs(v / D - 1.0) > tol}
    print(f"\nlargest deviation from the true D_c: {max(abs(v / D - 1.0) for v in out.values()):.1%}  (tolerance {tol:.0%})")
    assert not bad, f'fitters off by more than {tol:.0%}: {bad}'
    print('PASS: all three fitters return the D_c of the same equation -- no factor between their conventions')


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.10)
