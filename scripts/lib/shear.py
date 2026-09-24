"""shear.py -- shared analysis library for the shear_slab notebooks (2026-09-23).

Used by  scripts/shear_analysis_single.ipynb  (one shear-strain level)  and
         scripts/shear_analysis_sweep.ipynb   (every level of a strain sweep).

The shear counterpart of lib/triaxial.py, laid out section for section like it
so the two notebooks read like the compression ones and -- more important -- so
G is estimated the way the compression notebooks estimate M and G:

  * the eps = 0 reference of the compression deck is the gamma = 0 REFERENCE
    window of shear_slab.lmp (Phase 1.5, _ref files): every G is an INCREMENT
    from that reading (G_SUBTRACT_REF), i.e. a slope, like M and G there;
  * M_network <-> G_network:  plateau-averaged polymer sigma_xz(z) profile over
    the last plateau_frac of the hold, mean over the wall-trimmed interior bins,
    CI = t-interval over the bins (+) the reference's own interval in quadrature;
  * M_piston  <-> G_series:   the fine block-averaged bulk sigma_p,xz(t) series
    (stress_series_*, the analogue of piston_force_avg) read over the longest
    drift-free trailing window with the same circular block bootstrap;
  * the strain is the PRESCRIBED plate-based gamma the deck records (plate
    x-displacement / plate separation), averaged over the plateau window;
  * D_c from the transverse (even-sine) relaxation fit of u_x(z,t) during the
    hold, kappa = D_c/G  (= k/eta, the same kappa as D_c/M in compression);
  * P_th = -(1/3) tr(sigma^t) in the bulk is the P* = 1.5 check.

Shared machinery (style, readers, block bootstrap, plateau window, formatting,
smart_legend, annotate_box, the Expanse puller) is imported from triaxial.py;
nothing statistical is re-implemented here.

Layout of this file
    0. style                         (triaxial's)
    1. Config                        knobs, path builders, NSTEPS auto-detect, _g<level> tags
    2. file readers                  tensor / series / profile files of shear_slab.lmp
    3. reference state               load_reference   (gamma = 0, shared by every level)
    4. one strain level              load_level       (stresses, plateau, G, N1/N2, P_th, D_c, kappa)
    5. Expanse sync                  sync_from_expanse (via triaxial.sync_pull)
    6. plotting primitives           evolution profiles in z/L_z with the plate planes
    7. figures, single level         fig_strain ... fig_thermo_pressure
    8. figures, sweep                fig_*_sweep
    9. summaries                     print_summary, print_hold_check

Geometry (shear_slab.lmp since 2026-09-22): z = gap (plate normal), x = shear
direction, y neutral; the network is periodic in x and y; profiles are z-binned
over the BULK region (plate_excl = 3 sigma inside each plate plane) in reduced
coordinates z/L_z.
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import triaxial as tri  # noqa: E402

# ===========================================================================
#  0. STYLE  (triaxial's palette and rcParams)
# ===========================================================================
WONG, EVO_CMAP, GEL_SHADE = tri.WONG, tri.EVO_CMAP, tri.GEL_SHADE
level_color, setup_style = tri.level_color, tri.setup_style
sig, fmt_step, fmt_val_unc, fmt_ci, fmt_mu = tri.sig, tri.fmt_step, tri.fmt_val_unc, tri.fmt_ci, tri.fmt_mu
mean_ci, block_bootstrap_ci, plateau_window = tri.mean_ci, tri.block_bootstrap_ci, tri.plateau_window
smart_legend, annotate_box, robust_ylim, rolling_mean, subsample = (tri.smart_legend, tri.annotate_box,
                                                                     tri.robust_ylim, tri.rolling_mean, tri.subsample)
COMP = ('xx', 'yy', 'zz', 'xy', 'xz', 'yz')          # column order of every tensor/profile file
CI = {c: i for i, c in enumerate(COMP)}
SERIES_COLS = ('sp_xx', 'sp_yy', 'sp_zz', 'sp_xz', 'ss_xx', 'ss_yy', 'ss_zz', 'ss_xz')


# ===========================================================================
#  1. CONFIG
# ===========================================================================
@dataclass
class Config:
    """Every knob of the shear analysis.  The notebook Config cell builds one of
    these; nothing else in the notebook needs editing to switch runs."""
    # ---- which run -------------------------------------------------------
    DATANAME: str
    INTERACTION: str                      # "epsSS_epsSP"
    RUN_ID: str                           # local folder under flow_data_local/{shear,plots/shear}
    STRAINS: list = field(default_factory=list)   # shear-strain targets as STRINGS ("0.1"), = STRAINS_LIST
                                          # in shear_slab.batch (files are tagged _g<level>)
    NSTEPS: object = None                 # the <steps> tag of the file names (= NSTEPS of the batch, one tag for every
                                          # level); None -> resolved from the files on disk, an int pins it
    base_dir: str = '../../flow_data_local'
    # ---- profile / window knobs -----------------------------------------
    binWidth: float = 2.0                 # z-bin (sigma); must match shear_slab.lmp
    plate_excl: float = 3.0               # bulk region starts this far inside each plate plane (shear_slab.lmp plate_excl)
    wall_margin: float = 4.0              # "interior" = bins >= this far inside each plate plane (as triaxial's wall_margin)
    n_curves: int = 10                    # evolution curves drawn per profile
    Ncount_min: int = 200                 # min polymer atoms per z-bin to count as gel (profiles, D_c domain)
    dt_lj: float = 0.005
    ci_level: float = 0.95
    plateau_frac: float = 0.25            # FIXED trailing fraction of the hold for the profile plateau means
    plateau_frac_auto: float = 0.45       # LONGEST candidate trailing window for plateau_window (series -> G_series)
    roll_win: int = 11                    # rolling-mean window (samples) for the series plots
    P_BARO: float = 1.5                   # bath pressure of the compression runs = the P_th check line
    # ---- G as an increment from the gamma = 0 reference --------------------
    G_SUBTRACT_REF: bool = True           # G = (stress - its own gamma = 0 reading) / gamma for BOTH estimators
                                          # (needs the _ref files of Phase 1.5; older runs fall back to absolute, flagged)
    # ---- D_c transverse relaxation fit ----------------------------------
    DC_N_MODES: int = 5
    DC_TRIM_BINS: int = 1
    DC_FRAC_EARLY: float = 1.0
    DC_BOUNDS: tuple = (1e-6, 1.0)
    DC_TARGET_RESID: float = 0.01
    DC_SLOW_REF: float = 0.03             # sigma^2/tau; D_c(shear) ~ (G/M) D_c(compression) ~ 0.16 x 0.17, for the hold check
    # ---- Expanse ---------------------------------------------------------
    EXPANSE_HOST: str = 'login.expanse.sdsc.edu'
    EXPANSE_USER: str = 'dpollard'
    RUNS_ROOT: str = '/home/dpollard/Documents/lammps_runs/shear_slab'
    TRAJ_ROOT: str = '/expanse/lustre/scratch/dpollard/temp_project/lammps_trajectories'
    # ---- derived ---------------------------------------------------------
    DATA_DIR: Path = field(init=False)
    PLOT_DIR: Path = field(init=False)
    TRAJ_DIR: Path = field(init=False)
    sim_name: str = field(init=False)
    mode: str = field(init=False, default='shear')

    def __post_init__(self):
        assert self.STRAINS, 'STRAINS is empty -- list at least one level, e.g. ["0.1"]'
        self.STRAINS = [str(l) for l in self.STRAINS]
        base = Path(self.base_dir)
        self.DATA_DIR = base / 'shear' / self.RUN_ID
        self.PLOT_DIR = base / 'plots' / 'shear' / self.RUN_ID
        self.TRAJ_DIR = base / 'traj_files.nosync'
        for d in (self.DATA_DIR, self.PLOT_DIR, self.TRAJ_DIR):
            d.mkdir(parents=True, exist_ok=True)
        self._tags = {}
        if self.NSTEPS is not None:
            self.NSTEPS = int(self.NSTEPS)
        t = self.tag_for()
        self.sim_name = f'{self.DATANAME}_{self.INTERACTION}' + (f'_{t}' if t else '')
        print(f'Config (shear): {self.sim_name}  |  levels {self.STRAINS}\n'
              f'  data  {self.DATA_DIR}\n  plots {self.PLOT_DIR}\n'
              f'  file tag (NSTEPS): {t if t else "not local yet"}')

    # ---- <steps> tag: one tag for every file of a shear run --------------
    def _glob_tag(self):
        pat = re.compile(rf'^stress_tensor_polymer(?:_ref)?_{re.escape(self.DATANAME)}_'
                         rf'{re.escape(self.INTERACTION)}_(\d+)(?:_g[0-9.]+)?\.dat$')
        hits = sorted({int(m.group(1)) for f in self.DATA_DIR.glob('stress_tensor_polymer*.dat')
                       if (m := pat.match(f.name)) and f.stat().st_size > 0})
        return str(hits[-1]) if hits else None

    def tag_for(self, lvl=None):
        if isinstance(self.NSTEPS, int):
            return str(self.NSTEPS)
        if 'tag' not in self._tags:
            t = self._glob_tag()
            if t is None:
                return None
            self._tags['tag'] = t
        return self._tags['tag']

    def gsuf(self, lvl):
        return '' if lvl is None else f'_g{lvl}'

    def stem(self, lvl=None):
        t = self.tag_for()
        return f'{self.DATANAME}_{self.INTERACTION}_{t if t is not None else "*"}{self.gsuf(lvl)}'

    def path(self, name, lvl=None, ext='dat'):
        """<name>_<stem>[_g<lvl>].dat; reference files are named with their own '_ref' (lvl=None)."""
        return self.DATA_DIR / f'{name}_{self.stem(lvl)}.{ext}'

    def plot(self, stem, lvl=None):
        t = self.tag_for() or 'untagged'
        return self.PLOT_DIR / f'{stem}_{self.DATANAME}_{self.INTERACTION}_{t}{self.gsuf(lvl)}.png'


# ===========================================================================
#  2. FILE READERS
# ===========================================================================
def read_tensor(path):
    """fix ave/time scalar file (step + 6 components) -> (steps, array[n, 6]); None if absent."""
    a = tri.load2c(path, 7)
    if a is None:
        return None, None
    return a[:, 0].astype(int), a[:, 1:7]


def read_series(path):
    """stress_series file (step + sp_xx sp_yy sp_zz sp_xz ss_xx ss_yy ss_zz ss_xz) -> dict of arrays."""
    a = tri.load2c(path, 9)
    if a is None:
        return None
    d = {'step': a[:, 0].astype(int)}
    for j, c in enumerate(SERIES_COLS):
        d[c] = a[:, 1 + j]
    return d


def read_profiles(path):
    """stress_profile_z file (ave/chunk: chunk, z/Lz, Ncount, xx yy zz xy xz yz) ->
    dict(ts, zf, N[n, nb], stack{comp: [n, nb]}); None if absent/empty."""
    if not Path(path).exists():
        return None
    snaps = tri.read_ave_chunk_file(path)
    snaps = [s for s in snaps if s[1].shape[1] >= 9]
    if not snaps:
        return None
    nb = min(s[1].shape[0] for s in snaps)
    ts = np.array([s[0] for s in snaps])
    zf = snaps[-1][1][:nb, 1]
    N = np.array([s[1][:nb, 2] for s in snaps])
    stack = {c: np.array([s[1][:nb, 3 + j] for s in snaps]) for j, c in enumerate(COMP)}
    return dict(ts=ts, zf=zf, N=N, stack=stack)


def read_disp(path):
    """disp_x_polymer file (ave/chunk: chunk, z/Lz, Ncount, u_x) -> dict(ts, zf, N, ux)."""
    if not Path(path).exists():
        return None
    snaps = [s for s in tri.read_ave_chunk_file(path) if s[1].shape[1] >= 4]
    if len(snaps) < 2:
        return None
    nb = min(s[1].shape[0] for s in snaps)
    return dict(ts=np.array([s[0] for s in snaps]), zf=snaps[-1][1][:nb, 1],
                N=np.array([s[1][:nb, 2] for s in snaps]), ux=np.array([s[1][:nb, 3] for s in snaps]))


def _first_col(path, ncol, col):
    a = tri.load2c(path, ncol)
    return None if a is None else a[:, col]


# ===========================================================================
#  3. REFERENCE STATE  (gamma = 0; shared by every level)
# ===========================================================================
def _geometry(cfg, R, lvl_for_fallback=None):
    """Box, plate planes (reduced z), bulk and interior masks on the profile grid."""
    B = tri.load2c(cfg.path('box_dimensions'), 4)
    if B is None and lvl_for_fallback is not None:
        B = tri.load2c(cfg.path('box_dimensions', lvl_for_fallback), 4)
    if B is not None:
        R['LX'], R['LY'], R['LZ'] = float(np.mean(B[-5:, 1])), float(np.mean(B[-5:, 2])), float(np.mean(B[-5:, 3]))
    else:
        R['LX'] = R['LY'] = R['LZ'] = np.nan
    R['AREA'] = R['LX'] * R['LY']
    # plate separation and planes: shear_strain col 3 (plate_sep); plate_pressure (absolute planes)
    ps = _first_col(cfg.path('shear_strain'), 4, 2)
    if ps is None and lvl_for_fallback is not None:
        ps = _first_col(cfg.path('shear_strain', lvl_for_fallback), 4, 2)
    R['plate_sep'] = float(ps[0]) if ps is not None else np.nan
    PP = tri.load2c(cfg.path('plate_pressure'), 5)
    R['z_top_abs'], R['z_bot_abs'] = (float(PP[0, 1]), float(PP[0, 2])) if PP is not None else (np.nan, np.nan)
    if np.isfinite(R['z_top_abs']) and not np.isfinite(R['plate_sep']):
        R['plate_sep'] = R['z_top_abs'] - R['z_bot_abs']
    # masks on the profile grid: bulk = populated bins; the plate planes sit plate_excl outside
    # the bulk edges, so interior = bins >= (wall_margin - plate_excl) inside the bulk edges
    zf, N = R['zf'], R['N_ref']
    bulk = N.min(axis=0) > cfg.Ncount_min
    R['in_bulk'] = bulk
    if bulk.any():
        i_lo, i_hi = int(np.argmax(bulk)), int(len(bulk) - 1 - np.argmax(bulk[::-1]))
        dz = cfg.binWidth / R['LZ'] if np.isfinite(R['LZ']) else (zf[1] - zf[0])
        R['zf_bot'] = float(zf[i_lo] - 0.5 * dz - cfg.plate_excl / R['LZ']) if np.isfinite(R['LZ']) else np.nan
        R['zf_top'] = float(zf[i_hi] + 0.5 * dz + cfg.plate_excl / R['LZ']) if np.isfinite(R['LZ']) else np.nan
        extra = max(0, int(np.ceil((cfg.wall_margin - cfg.plate_excl) / cfg.binWidth - 1e-9)))
        interior = bulk.copy()
        interior[:i_lo + extra] = False
        interior[len(bulk) - extra:] = False if extra else interior[len(bulk) - extra:]
        interior[i_hi + 1 - extra:] = False
        R['interior'] = interior
    else:
        R['zf_bot'] = R['zf_top'] = np.nan
        R['interior'] = bulk
    return R


def load_reference(cfg, verbose=True):
    """The gamma = 0 reference state: geometry + the _ref stress tensors, profiles and
    fine series of Phase 1.5.  Runs older than 2026-09-23 have no _ref files: the
    geometry is then taken from the first level and every G is ABSOLUTE (flagged)."""
    say = print if verbose else (lambda *a, **k: None)
    R = dict(have_ref=False)
    first = cfg.STRAINS[0]
    prof = read_profiles(cfg.path('stress_profile_z_polymer_ref'))
    src = 'reference'
    if prof is None:
        prof = read_profiles(cfg.path('stress_profile_z_polymer', first))
        src = f'level _g{first} (no _ref files)'
        if prof is None:
            raise FileNotFoundError(f'no stress_profile_z_polymer file for the reference or level {first} in {cfg.DATA_DIR}'
                                    ' -- run the sync cell')
    R['zf'], R['N_ref'] = prof['zf'], prof['N']
    _geometry(cfg, R, first)
    say(f"reference geometry from the {src} profile: {len(R['zf'])} bins, bulk {int(R['in_bulk'].sum())} bins, "
        f"interior {int(R['interior'].sum())} bins;  Lz = {R['LZ']:.2f}, plate_sep = {R['plate_sep']:.2f} sigma, "
        f"plates at z/Lz = {R['zf_bot']:.3f} / {R['zf_top']:.3f}")
    # ---- reference stresses (only when the _ref files exist) ----------------
    tp_ts, tp = read_tensor(cfg.path('stress_tensor_polymer_ref'))
    ts_ts, tsv = read_tensor(cfg.path('stress_tensor_solvent_ref'))
    if tp is None or src != 'reference':
        say('  NOTE: no gamma = 0 reference files (pre-2026-09-23 run, or not synced) -> G will be ABSOLUTE')
        R['ref'] = None
        return R
    R['have_ref'] = True
    Rf = dict(ts=tp_ts, tensor_p=tp, tensor_s=tsv)
    pp = prof
    sp = read_profiles(cfg.path('stress_profile_z_solvent_ref'))
    Rf['prof_p'], Rf['prof_s'] = pp, sp
    im = R['interior']
    # polymer sigma_xz: interior mean per snapshot -> mean + t-interval over the snapshots
    per_snap = np.array([np.nanmean(pp['stack']['xz'][i][im]) for i in range(len(pp['ts']))])
    m, lo, hi = mean_ci(per_snap, cfg.ci_level)
    Rf['sp_xz_int'], Rf['sp_xz_int_half'] = float(m), float(0.5 * (hi - lo))
    if sp is not None:
        tot = pp['stack']['xz'] + sp['stack']['xz']
        per_t = np.array([np.nanmean(tot[i][im]) for i in range(len(pp['ts']))])
        m, lo, hi = mean_ci(per_t, cfg.ci_level)
        Rf['st_xz_int'], Rf['st_xz_int_half'] = float(m), float(0.5 * (hi - lo))
    # profile means with bands (mean over the reference snapshots, per bin)
    Rf['prof_mean'] = {c: mean_ci(pp['stack'][c], cfg.ci_level) for c in COMP}
    if sp is not None:
        Rf['prof_tot_mean'] = {c: mean_ci(pp['stack'][c] + sp['stack'][c], cfg.ci_level) for c in COMP}
    # bulk tensor means (P_th, N1, N2 of the reference)
    T = tp + (tsv if tsv is not None else 0.0)
    Rf['Pth'] = mean_ci(-(T[:, CI['xx']] + T[:, CI['yy']] + T[:, CI['zz']]) / 3.0, cfg.ci_level)
    Rf['N1_p'] = mean_ci(tp[:, CI['xx']] - tp[:, CI['yy']], cfg.ci_level)
    Rf['N2_p'] = mean_ci(tp[:, CI['yy']] - tp[:, CI['zz']], cfg.ci_level)
    # fine series: block bootstrap of the bulk polymer sigma_xz (the analogue of P_ref)
    S = read_series(cfg.path('stress_series_ref'))
    Rf['series'] = S
    if S is not None and len(S['step']) >= 2:
        m, lo, hi, blk, tau = block_bootstrap_ci(S['sp_xz'], cfg.ci_level)
        Rf['S_ref'], Rf['S_ref_lo'], Rf['S_ref_hi'] = m, lo, hi
    else:
        Rf['S_ref'] = Rf['sp_xz_int']; Rf['S_ref_lo'] = Rf['sp_xz_int'] - Rf['sp_xz_int_half']; Rf['S_ref_hi'] = Rf['sp_xz_int'] + Rf['sp_xz_int_half']
    R['ref'] = Rf
    say(f"  reference (gamma = 0, {len(pp['ts'])} snapshots): <sigma_p,xz>_int = {Rf['sp_xz_int']:+.5f} ± {Rf['sp_xz_int_half']:.5f}"
        f"   series <sigma_p,xz> = {Rf['S_ref']:+.5f} [{Rf['S_ref_lo']:+.5f}, {Rf['S_ref_hi']:+.5f}]"
        f"   P_th(bulk) = {Rf['Pth'][0]:.4f} [{Rf['Pth'][1]:.4f}, {Rf['Pth'][2]:.4f}]  (P_bath {cfg.P_BARO})")
    return R


# ===========================================================================
#  4. ONE STRAIN LEVEL
# ===========================================================================
def _plateau_mean(x):
    return float(np.nanmean(x)) if np.size(x) else np.nan


def fit_Dc(cfg, R, disp, t_hold):
    """Transverse relaxation fit of u_x(z,t)/L on the bulk bins during the hold.
    The plates are frozen, so the hold-referenced displacement vanishes at both
    plates and, by antisymmetry of the shear, at the gap centre: even sine modes
    sin(2 pi k zhat), zhat = (z - z_bot)/plate_sep, each relaxing as
    exp(-4 pi^2 k^2 D_c t / L^2) with free amplitudes (the hold-onset state is
    fitted).  The slowest mode has tau_1 = L^2/(4 pi^2 D_c)."""
    if disp is None or not np.isfinite(R['plate_sep']):
        return None
    L = float(R['plate_sep'])
    ok = disp['N'].min(axis=0) > cfg.Ncount_min
    idx = np.where(ok)[0]
    if len(idx) < 4 + 2 * cfg.DC_TRIM_BINS:
        return None
    i_lo, i_hi = idx[0], idx[-1]
    n_gel = i_hi - i_lo + 1
    # bin centres mapped onto the plate-to-plate gap: the populated bins span
    # [plate_excl, L - plate_excl] in absolute distance from the bottom plate
    scale = (L - 2.0 * cfg.plate_excl) / (n_gel * cfg.binWidth)
    zhat_all = (cfg.plate_excl + (np.arange(len(ok)) - i_lo + 0.5) * cfg.binWidth * scale) / L
    if cfg.DC_TRIM_BINS:
        idx = idx[cfg.DC_TRIM_BINS:-cfg.DC_TRIM_BINS]
    zf = zhat_all[idx]
    uhat = disp['ux'] / L
    t_lj = (disp['ts'] - t_hold) * cfg.dt_lj
    early = np.where((t_lj > 0) & (t_lj <= cfg.DC_FRAC_EARLY * t_lj[-1]))[0]
    if len(early) < 2:
        return None
    kk = np.arange(1, cfg.DC_N_MODES + 1, dtype=float)

    def basis(zh, Dc, t):
        return (1.0 - np.exp(-4.0 * np.pi ** 2 * kk ** 2 * Dc * t / L ** 2))[None, :] * np.sin(2.0 * np.pi * np.outer(zh, kk))

    def amps(Dc):
        X = np.vstack([basis(zf, Dc, t_lj[i]) for i in early])
        y = np.concatenate([uhat[i][idx] for i in early])
        return np.linalg.lstsq(X, y, rcond=None)[0]

    def resid(Dc):
        A = amps(Dc)
        return float(sum(np.sum((basis(zf, Dc, t_lj[i]) @ A - uhat[i][idx]) ** 2) for i in early))

    Dc = float(minimize_scalar(resid, bounds=cfg.DC_BOUNDS, method='bounded').x)
    A = amps(Dc)
    y_all = np.concatenate([uhat[i][idx] for i in early])
    p_all = np.concatenate([basis(zf, Dc, t_lj[i]) @ A for i in early])
    ss_t = np.sum((y_all - y_all.mean()) ** 2)
    R2 = float(1.0 - np.sum((y_all - p_all) ** 2) / ss_t) if ss_t > 1e-30 else np.nan
    u_model = lambda zh, t: basis(np.asarray(zh, float), Dc, t) @ A
    u_inf = lambda zh: np.sin(2.0 * np.pi * np.outer(np.asarray(zh, float), kk)) @ A
    hold_T = float(t_lj[-1])
    check = {}
    for tag, Dx in (('fit', Dc), ('slow', cfg.DC_SLOW_REF)):
        tau1 = L * L / (4.0 * np.pi ** 2 * Dx)
        end, avg = tri.hold_residual(hold_T, tau1, cfg.plateau_frac)
        check[tag] = dict(Dc=Dx, tau1=tau1, end=end, avg=avg, need=tri.n_tau_needed(tau1, cfg.DC_TARGET_RESID, cfg.plateau_frac),
                          ok=avg <= cfg.DC_TARGET_RESID)
    return dict(Dc=Dc, A=A, R2=R2, L=L, zf=zf, idx=idx, uhat=uhat, early=early, t_lj=t_lj, ts=disp['ts'],
                u_model=u_model, u_inf=u_inf, kk=kk, hold_T=hold_T, hold_check=check, at_bound=(Dc >= 0.999 * cfg.DC_BOUNDS[1]))


def load_level(cfg, R, lvl, verbose=True):
    """Everything for ONE shear-strain level `lvl` (string, e.g. "0.1"): strain history,
    bulk stress tensors, z-profiles, the fine sigma_p,xz series + plateau, G (profile and
    series estimators, increments from gamma = 0), N1/N2, P_th, the D_c fit and kappa.
    Returns a dict L, or None if the core files are missing."""
    say = print if verbose else (lambda *a, **k: None)
    L = dict(lvl=lvl, gamma_target=float(lvl))
    say(f'\n=== level _g{lvl}  (target shear strain {float(lvl):.4f}) ===')
    ts, tp = read_tensor(cfg.path('stress_tensor_polymer', lvl))
    if tp is None:
        say(f'  level {lvl}: stress_tensor_polymer file missing -- skipping level')
        return None
    _, tsv = read_tensor(cfg.path('stress_tensor_solvent', lvl))
    L['ts'], L['tensor_p'], L['tensor_s'] = ts, tp, tsv
    t0, t1 = float(ts[0]), float(ts[-1])
    nfreq = float(ts[1] - ts[0]) if len(ts) > 1 else 0.0
    L['t_hold'] = t0 - nfreq                     # the hold started one averaging epoch before the first row
    L['halt_ts'] = int(t1 - cfg.plateau_frac * (t1 - t0))
    L['plat'] = ts >= L['halt_ts']
    if L['plat'].sum() < 1:
        L['plat'][-1] = True
    say(f'  hold: {len(ts)} tensor snapshots, steps {int(t0)} -> {int(t1)} (hold start ~{int(L["t_hold"])});  '
        f'plateau window: steps >= {L["halt_ts"]} (last {cfg.plateau_frac:.0%}, {int(L["plat"].sum())} snapshots)')

    # ---- strain: plate-based gamma (prescribed), surface-COM gamma (slip diagnostic) ----
    SS = tri.load2c(cfg.path('shear_strain', lvl), 4)
    if SS is not None:
        L['strain_ts'], L['gamma_ts'] = SS[:, 0], SS[:, 3]
        L['plate_sep'] = float(SS[0, 2])
        pl = SS[:, 0] >= L['halt_ts']
        L['gamma'] = float(np.mean(SS[pl, 3])) if pl.any() else float(SS[-1, 3])
    else:
        L['gamma'] = L['gamma_target']
        L['plate_sep'] = R['plate_sep']
        say('  NOTE: no shear_strain file -> gamma = the target value')
    SU = tri.load2c(cfg.path('shear_strain_surface', lvl), 3)
    if SU is not None:
        L['surf_ts'], L['gamma_surf_ts'] = SU[:, 0], SU[:, 2]
        pl = SU[:, 0] >= L['halt_ts']
        L['gamma_surf'] = float(np.mean(SU[pl, 2])) if pl.any() else float(SU[-1, 2])
    g = L['gamma']
    say(f"  strain: held plate-based gamma = {g:.5f} (G denominator; target {L['gamma_target']}), "
        f"surface-COM gamma = {L.get('gamma_surf', np.nan):.5f}, plate_sep = {L['plate_sep']:.2f}")

    # ---- profiles (bulk z-bins) ---------------------------------------------
    pp = read_profiles(cfg.path('stress_profile_z_polymer', lvl))
    sp = read_profiles(cfg.path('stress_profile_z_solvent', lvl))
    L['prof_p'], L['prof_s'] = pp, sp
    im = R['interior']
    if pp is not None:
        pl_p = pp['ts'] >= L['halt_ts']
        if pl_p.sum() < 1:
            pl_p[-1] = True
        L['prof_plat'] = pl_p
        L['prof_p_plat'] = {c: np.nanmean(pp['stack'][c][pl_p], axis=0) for c in COMP}
        if sp is not None:
            L['prof_t'] = {c: pp['stack'][c] + sp['stack'][c] for c in COMP}
            L['prof_t_plat'] = {c: np.nanmean(L['prof_t'][c][pl_p], axis=0) for c in COMP}
            L['Pth_prof'] = -(L['prof_t']['xx'] + L['prof_t']['yy'] + L['prof_t']['zz']) / 3.0

    # ---- fine series + plateau window (the analogue of the piston force) ----
    S = read_series(cfg.path('stress_series', lvl))
    L['series'] = S
    if S is not None:
        hold = S['step'] >= L['t_hold']
        if hold.sum() >= 4:
            L['PF'] = plateau_window(S['step'][hold], S['sp_xz'][hold], cfg.plateau_frac_auto, cfg.ci_level)
            p = L['PF']
            say(f"  series plateau <sigma_p,xz> = {p['mean']:+.5f} [{p['lo']:+.5f}, {p['hi']:+.5f}]  "
                f"(auto window last {p['frac']:.0%} of the hold, n={p['n']}, block={p['block']}, tau~{p['tau']:.1f})"
                + ('  DRIFT WARNING: no drift-free window -- extend the hold' if p.get('warn') else ''))
    else:
        say('  NOTE: no stress_series file -> G_series skipped')

    # ---- G: network (profile) and series estimators, increments from gamma = 0 ----
    Rf = R.get('ref')
    sub = bool(cfg.G_SUBTRACT_REF) and Rf is not None
    L['G_ref'] = 'measured' if sub else ('absent' if cfg.G_SUBTRACT_REF else 'off')
    L['G'] = {}
    if pp is not None:
        ref_v = Rf['sp_xz_int'] if sub else 0.0
        ref_h = Rf['sp_xz_int_half'] if sub else 0.0
        bins = L['prof_p_plat']['xz'][im]
        bins = bins[np.isfinite(bins)]
        Ga, Ga_lo, Ga_hi = mean_ci(bins / g, cfg.ci_level)
        m, lo, hi = mean_ci((bins - ref_v) / g, cfg.ci_level)
        half = np.sqrt((0.5 * (hi - lo)) ** 2 + (ref_h / g) ** 2)
        L['G']['net'] = dict(G=float(m), lo=float(m - half), hi=float(m + half), abs=float(Ga), abs_lo=float(Ga_lo),
                             abs_hi=float(Ga_hi), ref=ref_v, nbins=len(bins), sigma=float(np.mean(bins)))
        if sp is not None:
            tb = L['prof_t_plat']['xz'][im]
            tb = tb[np.isfinite(tb)]
            rv = Rf['st_xz_int'] if (sub and 'st_xz_int' in Rf) else 0.0
            rh = Rf['st_xz_int_half'] if (sub and 'st_xz_int_half' in Rf) else 0.0
            m, lo, hi = mean_ci((tb - rv) / g, cfg.ci_level)
            half = np.sqrt((0.5 * (hi - lo)) ** 2 + (rh / g) ** 2)
            Ta, Ta_lo, Ta_hi = mean_ci(tb / g, cfg.ci_level)
            L['G']['tot'] = dict(G=float(m), lo=float(m - half), hi=float(m + half), ref=rv, sigma=float(np.mean(tb)),
                                 abs=float(Ta), abs_lo=float(Ta_lo), abs_hi=float(Ta_hi),
                                 solvent_share=float(np.mean(L['prof_t_plat']['xz'][im] - L['prof_p_plat']['xz'][im]) / max(abs(np.mean(tb)), 1e-30)))
    if 'PF' in L:
        p = L['PF']
        rv = Rf['S_ref'] if sub else 0.0
        rh = 0.5 * (Rf['S_ref_hi'] - Rf['S_ref_lo']) if sub else 0.0
        half = np.sqrt((0.5 * (p['hi'] - p['lo'])) ** 2 + rh ** 2) / g
        L['G']['ser'] = dict(G=(p['mean'] - rv) / g, lo=(p['mean'] - rv) / g - half, hi=(p['mean'] - rv) / g + half,
                             abs=p['mean'] / g, abs_lo=p['lo'] / g, abs_hi=p['hi'] / g, ref=rv, sigma=p['mean'])
    how = 'increment from gamma = 0' if sub else ('ABSOLUTE (no _ref files)' if cfg.G_SUBTRACT_REF else 'absolute')
    for key, name in (('net', 'G_network (profile, interior bins)'), ('ser', 'G_series   (plateau window)'),
                      ('tot', 'G_total    (sigma^t_xz, poroelastic check)')):
        Gd = L['G'].get(key)
        if Gd:
            say(f"  {name} = {Gd['G']:.4f} [{Gd['lo']:.4f}, {Gd['hi']:.4f}]  ({how}"
                + (f"; ref {Gd['ref']:+.5f} subtracted; absolute {Gd.get('abs', np.nan):.4f}" if sub else '')
                + (f"; {Gd['nbins']} bins" if 'nbins' in Gd else '')
                + (f"; solvent share of sigma_xz {Gd['solvent_share']:+.1%}" if 'solvent_share' in Gd else '') + ')')

    # ---- normal stress differences and P_th from the bulk tensors ----------
    tt = tp + (tsv if tsv is not None else 0.0)
    pl = L['plat']
    L['N1_p'] = mean_ci((tp[:, CI['xx']] - tp[:, CI['yy']])[pl], cfg.ci_level)
    L['N2_p'] = mean_ci((tp[:, CI['yy']] - tp[:, CI['zz']])[pl], cfg.ci_level)
    L['N1_t'] = mean_ci((tt[:, CI['xx']] - tt[:, CI['yy']])[pl], cfg.ci_level)
    L['N2_t'] = mean_ci((tt[:, CI['yy']] - tt[:, CI['zz']])[pl], cfg.ci_level)
    L['Pth_ts'] = -(tt[:, CI['xx']] + tt[:, CI['yy']] + tt[:, CI['zz']]) / 3.0
    L['Pth'] = mean_ci(L['Pth_ts'][pl], cfg.ci_level)
    L['sxz_p_ts'] = tp[:, CI['xz']]
    L['sxz_t_ts'] = tt[:, CI['xz']]
    say(f"  bulk P_th (plateau) = {L['Pth'][0]:.4f} [{L['Pth'][1]:.4f}, {L['Pth'][2]:.4f}]  (P_bath {cfg.P_BARO});  "
        f"N1_p/sigma_p,xz = {L['N1_p'][0] / max(abs(np.nanmean(L['sxz_p_ts'][pl])), 1e-30):+.3f}, "
        f"N2_p/sigma_p,xz = {L['N2_p'][0] / max(abs(np.nanmean(L['sxz_p_ts'][pl])), 1e-30):+.3f}")

    # ---- D_c and kappa --------------------------------------------------------
    L['Dc'] = fit_Dc(cfg, R, read_disp(cfg.path('disp_x_polymer', lvl)), L['t_hold'])
    if L['Dc'] is None:
        say('  NOTE: no D_c (missing disp_x_polymer, no plate_sep, or too few bins)')
    else:
        F = L['Dc']
        say(f"  D_c (shear) = {F['Dc']:.4e} sigma^2/tau  (R^2 = {F['R2']:.3f}; L = {F['L']:.2f}, hold = {F['hold_T']:.0f} tau)"
            + ('   AT THE FIT BOUND -- not converged' if F['at_bound'] else ''))
        L['kappa'] = {}
        for key in ('net', 'ser'):
            Gd = L['G'].get(key)
            if Gd and Gd['G'] > 0 and Gd['lo'] > 0:
                L['kappa'][key] = dict(k=F['Dc'] / Gd['G'], lo=F['Dc'] / Gd['hi'], hi=F['Dc'] / Gd['lo'])
        if L['kappa']:
            say('  kappa = D_c/G:  ' + '   '.join(f"{k}: {v['k']:.4e} [{v['lo']:.4e}, {v['hi']:.4e}]" for k, v in L['kappa'].items()))
    return L


# ===========================================================================
#  5. EXPANSE SYNC
# ===========================================================================
_REF_DAT = ('stress_tensor_polymer_ref', 'stress_tensor_solvent_ref', 'stress_profile_z_polymer_ref',
            'stress_profile_z_solvent_ref', 'stress_series_ref')
_RUN_DAT = ('shear_strain', 'plate_pressure', 'box_dimensions', 'gel_dimensions_rg', 'gel_volume_rg', 'gel_volume_bb')
_PROD_DAT = ('stress_tensor_polymer', 'stress_tensor_solvent', 'stress_profile_z_polymer', 'stress_profile_z_solvent',
             'stress_series', 'shear_strain', 'shear_strain_surface', 'disp_x_polymer', 'gel_dimensions_rg',
             'gel_dimensions_bb', 'gel_volume_rg', 'gel_volume_bb', 'box_dimensions', 'polymer_com')
_REQUIRED = ('stress_tensor_polymer', 'stress_profile_z_polymer', 'shear_strain')


def sync_files(cfg, levels=None):
    levels = cfg.STRAINS if levels is None else [str(l) for l in levels]
    data = [cfg.path(n) for n in _REF_DAT + _RUN_DAT]
    req = []
    for l in levels:
        data += [cfg.path(n, l) for n in _PROD_DAT]
        req += [cfg.path(n, l) for n in _REQUIRED]
    return data, [], req


def sync_from_expanse(cfg, levels=None, force=False):
    """Pull every file the shear notebooks read from Expanse in ONE login (the
    triaxial puller); the _ref files are optional for runs older than 2026-09-23."""
    data, traj, req = sync_files(cfg, levels)
    optional = {cfg.path(n).name for n in _REF_DAT}
    tri.sync_pull(cfg, data, traj, req, optional, force, refresh=lambda: sync_files(cfg, levels))


# ===========================================================================
#  6. PLOTTING PRIMITIVES
# ===========================================================================
def mark_plates(ax, R):
    for zf in (R.get('zf_bot'), R.get('zf_top')):
        if zf is not None and np.isfinite(zf):
            ax.axvline(zf, color='k', ls='-.', lw=1.5, alpha=0.85, zorder=4)


def shade_bulk(ax, R):
    b = R['in_bulk']
    if b.any():
        ax.axvspan(float(R['zf'][b].min()), float(R['zf'][b].max()), **GEL_SHADE)


def finish_axes(ax, ylabel, title):
    ax.set_xlabel(r'$z/L$  (gap direction)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    lo = R_XLIM[0] if R_XLIM else 0.0
    ax.set_xlim(lo, R_XLIM[1] if R_XLIM else 1.0)


R_XLIM = None


def _xlim_from(R):
    global R_XLIM
    if np.isfinite(R.get('zf_bot', np.nan)) and np.isfinite(R.get('zf_top', np.nan)):
        w = R['zf_top'] - R['zf_bot']
        R_XLIM = (max(0.0, R['zf_bot'] - 0.05 * w), min(1.0, R['zf_top'] + 0.05 * w))
    else:
        R_XLIM = (0.0, 1.0)


def _save(fig, cfg, stem, lvl=None):
    return tri._save(fig, cfg, stem, lvl)


def post_hold(cfg, ts, stack):
    return subsample(np.asarray(ts), np.asarray(stack), cfg.n_curves)


def plot_evolution(ax, cfg, R, ts, stack, ylabel, title, ref=None, band=None, colorbar=True, legend=True):
    """Time-coloured profiles (cividis), final curve bold black, bulk shaded, plate
    planes dash-dot; ref=(mean, lo, hi) draws the gamma = 0 reference dashed."""
    ts = np.asarray(ts)
    stack = np.asarray(stack)
    norm = Normalize(vmin=ts.min(), vmax=ts.max())
    cmap = plt.get_cmap(EVO_CMAP)
    zx = R['zf']
    b = R['in_bulk']
    if ref is not None:
        rm, rlo, rhi = ref
        ax.fill_between(zx[b], np.asarray(rlo)[b], np.asarray(rhi)[b], color=WONG['skyblue'], alpha=0.25, lw=0, zorder=1)
        ax.plot(zx[b], np.asarray(rm)[b], '--', color=WONG['blue'], lw=2.2, alpha=0.9, zorder=2, label=r'reference ($\gamma=0$)')
    for i in range(len(ts)):
        last = (i == len(ts) - 1)
        c = 'k' if last else cmap(norm(ts[i]))
        if band is not None:
            ax.fill_between(zx[b], (stack[i] - band[i])[b], (stack[i] + band[i])[b], color=c, alpha=(0.20 if last else 0.06), lw=0)
        ax.plot(zx[b], stack[i][b], '-', color=c, lw=(3.5 if last else 1.6), alpha=(1.0 if last else 0.75),
                zorder=(5 if last else 3), label=('final (plateau)' if last else None))
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    shade_bulk(ax, R)
    mark_plates(ax, R)
    finish_axes(ax, ylabel, title)
    if colorbar:
        sm = plt.cm.ScalarMappable(cmap=EVO_CMAP, norm=norm)
        sm.set_array([])
        ax.figure.colorbar(sm, ax=ax, fraction=0.046, pad=0.02).set_label('timestep')
        ax._tri_has_colorbar = True
    if ref is not None and legend:
        smart_legend(ax, fontsize=12)


def _final_profile(ax, R, zx_mask, m, lo, hi, color, label):
    b = R['in_bulk']
    ax.fill_between(R['zf'][b], np.asarray(lo)[b], np.asarray(hi)[b], color=color, alpha=0.25, lw=0, zorder=3)
    ax.plot(R['zf'][b], np.asarray(m)[b], '-', color=color, lw=2.8, zorder=4, label=label)


def _ref_profile(ax, R, key, comp, tot=False):
    Rf = R.get('ref')
    if Rf is None:
        return None
    src = Rf.get('prof_tot_mean' if tot else 'prof_mean')
    if src is None:
        return None
    m, lo, hi = src[comp]
    b = R['in_bulk']
    ax.fill_between(R['zf'][b], lo[b], hi[b], color='0.5', alpha=0.25, lw=0, zorder=1)
    ax.plot(R['zf'][b], m[b], '--', color='k', lw=2.0, alpha=0.9, zorder=2, label=r'reference ($\gamma=0$)')
    return m


def level_handles(levels, ref=False):
    h = [Line2D([0], [0], color=level_color(i), lw=3, label=fr'$\gamma={L["gamma"]:.3f}$') for i, L in enumerate(levels)]
    if ref:
        h.insert(0, Line2D([0], [0], color='k', ls='--', lw=2, label=r'reference ($\gamma=0$)'))
    return h


# ===========================================================================
#  7. FIGURES -- SINGLE LEVEL
# ===========================================================================
def fig_strain(cfg, R, levels, stem='strain_diagnostic'):
    """gamma vs step: solid plate-based gamma (prescribed), dashed surface-COM gamma
    (slip diagnostic), dotted target; shaded = plateau window."""
    fig, ax = plt.subplots(figsize=(10, 6.5), constrained_layout=True)
    any_ = False
    for i, L in enumerate(levels):
        col = level_color(i) if len(levels) > 1 else WONG['blue']
        if 'gamma_ts' not in L:
            continue
        any_ = True
        ax.plot(L['strain_ts'], L['gamma_ts'], '-', color=col, lw=2.2, label=fr'$\gamma$ plate-based (target {L["lvl"]})')
        if 'gamma_surf_ts' in L:
            ax.plot(L['surf_ts'], L['gamma_surf_ts'], '--', color=col, lw=1.5, alpha=0.7, label=fr'$\gamma_{{\rm surf}}$ (polymer surface COM)')
        ax.axvspan(L['halt_ts'], float(L['strain_ts'][-1]), color=col, alpha=0.06)
        ax.axhline(L['gamma_target'], color=col, ls=':', lw=1.0, alpha=0.6)
        ax.annotate(f"held: {sig(L['gamma'], 4)}  surf {sig(L.get('gamma_surf', np.nan), 4)}",
                    (L['strain_ts'][-1], L['gamma_ts'][-1]), textcoords='offset points', xytext=(-6, 9), ha='right',
                    va='bottom', fontsize=10, color=col, bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.8))
    ax.set_xlabel('time step')
    ax.set_ylabel(r'shear strain  $\gamma = \Delta x_{\rm plates}/L_{\rm plates}$')
    ax.set_title('Strain diagnostic: solid = plate-based $\\gamma$ (prescribed), dashed = surface COM,\n'
                 'dotted = target, shaded = plateau window', fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=10)
    if not any_:
        plt.close(fig)
        print('strain diagnostic skipped (no shear_strain files)')
        return None
    return _save(fig, cfg, stem, levels[0]['lvl'] if len(levels) == 1 else None)


def _stress_evo_panels(cfg, R, L, kind, stem, suptitle):
    """1 x 4 evolution panels (xz, xx, yy, zz) of kind 't' (total) or 'p' (polymer)."""
    if L.get('prof_p') is None or (kind == 't' and L.get('prof_t') is None):
        print(f'{stem} skipped (profile files missing)')
        return None
    comps = ('xz', 'xx', 'yy', 'zz')
    fig, axes = plt.subplots(1, 4, figsize=(30, 6.5), constrained_layout=True)
    fig.suptitle(suptitle, fontsize=13, fontweight='bold')
    stacks = L['prof_t'] if kind == 't' else L['prof_p']['stack']
    for k, (ax, comp) in enumerate(zip(axes, comps)):
        ts, ev = post_hold(cfg, L['prof_p']['ts'], stacks[comp])
        Rf = R.get('ref')
        ref = None
        if Rf is not None:
            src = Rf.get('prof_tot_mean' if kind == 't' else 'prof_mean')
            if src is not None:
                ref = src[comp]
        lab = (r'$\sigma^t_{%s}$' % comp) if kind == 't' else (r'$\sigma_{p,%s}$' % comp)
        title = f'({"abcd"[k]}) ' + ('total' if kind == 't' else 'polymer (network)') + ' ' + lab
        if kind == 't' and comp in ('xx', 'yy', 'zz'):
            Pb = cfg.P_BARO
            ev = -ev / Pb
            ref = None if ref is None else tuple(-np.asarray(r) / Pb for r in ref)
            ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6, zorder=1)
            lab = r'$-$' + lab + r'$/P_{\rm bath}$'
        plot_evolution(ax, cfg, R, ts, ev, lab + '$(z,t)$', title, ref=ref, legend=False)
        fin = ev[-1][R['interior']]
        robust_ylim(ax, list(ev) + ([ref[0]] if ref is not None else []), zmask=R['in_bulk'], pad=0.3,
                    include_zero=(comp == 'xz'))
        h, lb = ax.get_legend_handles_labels()
        note = 'plateau mean, interior = ' + fmt_mu(fin)
        h.append(Patch(alpha=0, label=note)); lb.append(note)
        smart_legend(ax, handles=h, labels=lb, fontsize=12)
    return _save(fig, cfg, stem, L['lvl'])


def fig_total_stress(cfg, R, L):
    return _stress_evo_panels(cfg, R, L, 't', 'total_stress_evolution',
                              f'Total stress: shear $\\sigma^t_{{xz}}$ and normal components / $P_{{\\rm bath}}={sig(cfg.P_BARO)}$, '
                              f'reference -> sheared (hold from step {fmt_step(L["t_hold"])}; plateau from {fmt_step(L["halt_ts"])})  |  {cfg.sim_name}')


def fig_partial_stress(cfg, R, L):
    """solvent | polymer | total sigma_xz(z) evolutions (the poroelastic split of the shear stress)."""
    if L.get('prof_p') is None or L.get('prof_s') is None:
        print('partial-stress figure skipped (profile files missing)')
        return None
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    fig.suptitle(f'Shear stress $\\sigma_{{xz}}(z)$: solvent, polymer and total  |  $\\gamma = {sig(L["gamma"], 4)}$  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    for k, (ax, st, lab, title) in enumerate(zip(axes,
            (L['prof_s']['stack']['xz'], L['prof_p']['stack']['xz'], L['prof_t']['xz']),
            (r'$\sigma_{s,xz}$', r'$\sigma_{p,xz}$', r'$\sigma^t_{xz}$'),
            ('(a) solvent partial (-> 0 at rest)', '(b) polymer partial = network shear stress', '(c) total'))):
        ts, ev = post_hold(cfg, L['prof_p']['ts'], st)
        plot_evolution(ax, cfg, R, ts, ev, lab + '$(z,t)$', title, legend=False)
        robust_ylim(ax, list(ev), zmask=R['in_bulk'], pad=0.3)
        note = 'plateau mean, interior = ' + fmt_mu(ev[-1][R['interior']])
        h, lb = ax.get_legend_handles_labels()
        h.append(Patch(alpha=0, label=note)); lb.append(note)
        smart_legend(ax, handles=h, labels=lb, fontsize=12)
    return _save(fig, cfg, 'partial_stress_evolution', L['lvl'])


def fig_network_stress(cfg, R, L):
    """Network shear stress sigma_p,xz(z): the gamma = 0 REFERENCE (dashed, grey band, ~0)
    and the FINAL plateau-averaged state (solid, 95 % band) -- the shear analogue of the
    compression notebooks' network-stress figure; the interior mean over gamma is G."""
    if L.get('prof_p') is None:
        print('network-stress figure skipped (no polymer profile)')
        return None
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(17, 6.5), constrained_layout=True)
    fig.suptitle(f'Network shear stress: reference and final equilibrated state  |  $\\gamma = {sig(L["gamma"], 4)}$  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    for ax, comp, title in ((axA, 'xz', r'(a) polymer $\sigma_{p,xz}(z)$  (G = interior mean / $\gamma$)'),
                            (axB, 'xx', r'(b) polymer $\sigma_{p,xx}(z)$  (normal, for $N_1$)')):
        rm = _ref_profile(ax, R, 'prof_mean', comp)
        m, lo, hi = mean_ci(L['prof_p']['stack'][comp][L['prof_plat']], cfg.ci_level)
        _final_profile(ax, R, None, m, lo, hi, WONG['blue'], f'final (plateau, {int(L["prof_plat"].sum())} snapshots)')
        ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
        shade_bulk(ax, R); mark_plates(ax, R)
        finish_axes(ax, r'$\sigma_{p,%s}(z)$' % comp, title)
        robust_ylim(ax, [m] + ([rm] if rm is not None else []), zmask=R['in_bulk'], pad=0.25)
        h, lb = ax.get_legend_handles_labels()
        note = 'final mean, interior = ' + fmt_mu(m[R['interior']]) + (f"  ->  $G_{{\\rm net}}$ = {sig(L['G']['net']['G'])}" if comp == 'xz' and 'net' in L['G'] else '')
        h.append(Patch(alpha=0, label=note)); lb.append(note)
        smart_legend(ax, handles=h, labels=lb, fontsize=12)
    return _save(fig, cfg, 'network_stress_final', L['lvl'])


def fig_series(cfg, R, L):
    """Bulk polymer shear stress sigma_p,xz(t) from the fine block-averaged series
    (drive + hold), linear and log |sigma|, with the auto plateau window -- the
    analogue of the compression notebooks' piston-pressure history."""
    S = L.get('series')
    if S is None:
        print('series figure skipped (no stress_series file)')
        return None
    st, x = S['step'], S['sp_xz']
    xr = rolling_mean(x, cfg.roll_win)
    fig, (axL, axG) = plt.subplots(1, 2, figsize=(17, 6), constrained_layout=True)
    fig.suptitle(f'Bulk polymer shear stress history (block-averaged every {int(st[1] - st[0]) if len(st) > 1 else "?"} steps)  |  '
                 f'$\\gamma = {sig(L["gamma"], 4)}$  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    for ax in (axL, axG):
        ax.axvline(L['t_hold'], color=WONG['blue'], ls='--', lw=1.6, alpha=0.7, label=f'hold start (step {fmt_step(L["t_hold"])})')
        if 'PF' in L:
            ax.axvspan(L['PF']['step0'], float(st[-1]), color=WONG['green'], alpha=0.10, label='plateau window (auto)')
        ax.set_xlabel('time step'); ax.grid(alpha=0.3)
    axL.plot(st, x, '-', color=WONG['vermillion'], lw=1.0, alpha=0.35, label=r'$\sigma_{p,xz}$ (block-avg)')
    axL.plot(st, xr, '-', color=WONG['vermillion'], lw=2.6, alpha=0.95, label=f'rolling mean ({cfg.roll_win})')
    if 'ss_xz' in S:
        axL.plot(st, S['ss_xz'], '-', color=WONG['skyblue'], lw=1.2, alpha=0.8, label=r'$\sigma_{s,xz}$ (solvent)')
    Rf = R.get('ref')
    if Rf is not None:
        axL.axhline(Rf['S_ref'], color='k', ls=':', lw=1.2, alpha=0.7, label=f"reference $\\gamma=0$: {sig(Rf['S_ref'], 2)}")
    axL.axhline(0, color='k', ls='--', lw=0.8, alpha=0.4)
    axL.set_ylabel(r'$\sigma_{xz}$  (bulk, LJ)'); axL.set_title('(a) shear stress')
    smart_legend(axL, fontsize=12)
    pos = np.abs(x) > 0
    axG.plot(st[pos], np.log(np.abs(x[pos])), '-', color=WONG['vermillion'], lw=1.0, alpha=0.35, label=r'$\ln|\sigma_{p,xz}|$')
    posr = np.abs(xr) > 0
    axG.plot(st[posr], np.log(np.abs(xr[posr])), '-', color=WONG['vermillion'], lw=2.6, alpha=0.95, label=f'rolling mean ({cfg.roll_win})')
    axG.set_ylabel(r'$\ln|\sigma_{p,xz}|$'); axG.set_title('(b) log shear stress (relaxation view)')
    smart_legend(axG, fontsize=12)
    if 'PF' in L:
        p = L['PF']
        annotate_box(axL, f"plateau $\\langle\\sigma_{{p,xz}}\\rangle$ = {fmt_val_unc(p['mean'], 0.5 * (p['hi'] - p['lo']))}\n"
                          f"last {p['frac']:.0%} of the hold, n={p['n']}, block={p['block']}", loc='upper right', fontsize=12)
    return _save(fig, cfg, 'shear_stress_history', L['lvl'])


def _pt(ax, x, d, marker, color, ms, label, hollow=False, key='G'):
    lo, hi = (d['lo'], d['hi']) if key == 'G' else (d['abs_lo'], d['abs_hi'])
    v = d[key] if key == 'G' else d['abs']
    ax.errorbar([x], [v], yerr=[[v - lo], [hi - v]], fmt=marker, ms=ms, color=color, capsize=(5 if hollow else 8),
                lw=(1.5 if hollow else 2.5), alpha=(0.7 if hollow else 1.0), label=label, **({'mfc': 'none'} if hollow else {}))


def fig_G(cfg, R, L):
    """Shear modulus, two panes (the layout of the compression notebooks' M figure):
    (a) the increment estimates next to the absolute stress/gamma values they replace
    (hollow) and the total-stress estimate (poroelastic check); (b) network and series
    estimates alone."""
    if not L['G']:
        print('G figure skipped (no profile / series files)')
        return None
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    ci = int(cfg.ci_level * 100)
    sub = L['G_ref'] == 'measured'
    fig.suptitle('Shear modulus' + (' (increment from $\\gamma=0$)' if sub else ' (absolute: no $\\gamma=0$ reference files)')
                 + f'   |   {cfg.RUN_ID}   |   $\\gamma = {sig(L["gamma"], 4)}$', fontsize=14, fontweight='bold')
    spec = (('net', 'o', WONG['blue'], 'network'), ('ser', 's', WONG['vermillion'], 'series'))
    for ax, pane in ((axA, 'a'), (axB, 'b')):
        for k, (key, mk, col, name) in enumerate(spec):
            d = L['G'].get(key)
            if d is None:
                continue
            _pt(ax, k, d, mk, col, 13, f"{name}  $G = {sig(d['G'])}$\n{ci}% CI [{sig(d['lo'])}, {sig(d['hi'])}]")
            ax.axhline(d['G'], color=col, ls='--', lw=1.2, alpha=0.5)
            if pane == 'a' and sub and 'abs' in d:
                _pt(ax, k + 0.15, d, mk, col, 10, f"absolute $\\sigma/\\gamma = {sig(d['abs'])}$\n(ref {d['ref']:+.5f} subtracted)", hollow=True, key='abs')
        if pane == 'a' and 'tot' in L['G']:
            d = L['G']['tot']
            ax.errorbar([2], [d['G']], yerr=[[d['G'] - d['lo']], [d['hi'] - d['G']]], fmt='D', ms=10, color=WONG['green'],
                        capsize=6, lw=1.8, mfc='none', label=f"total $\\sigma^t_{{xz}}/\\gamma = {sig(d['G'])}$\n(solvent share {d['solvent_share']:+.1%})")
        ax.set_xticks([0, 1, 2] if pane == 'a' else [0, 1])
        ax.set_xticklabels(([r'network' + '\n' + r'$\Delta\langle\sigma_{p,xz}\rangle_{\rm int}/\gamma$',
                             r'series' + '\n' + r'$\Delta\langle\sigma_{p,xz}\rangle_{\rm plateau}/\gamma$']
                            + (['total\n' + r'$\Delta\sigma^t_{xz}/\gamma$'] if pane == 'a' else [])), fontsize=13)
        ax.set_ylabel(r'$G$  (LJ units)')
        ax.set_xlim(-0.5, 2.6 if pane == 'a' else 1.5)
        ax.grid(axis='y', alpha=0.3)
        ax.set_title('(a) with the absolute values (hollow) and the total-stress check' if pane == 'a'
                     else '(b) shear modulus, two estimates', fontsize=13)
        smart_legend(ax, fontsize=11)
    return _save(fig, cfg, 'G_comparison', L['lvl'])


def fig_normal(cfg, R, L):
    """N1 = sigma_xx - sigma_yy and N2 = sigma_yy - sigma_zz vs step (polymer and total),
    plateau shaded -- the linearity / isotropy check (both should vanish)."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    tp, tt = L['tensor_p'], L['tensor_p'] + (L['tensor_s'] if L['tensor_s'] is not None else 0.0)
    fig.suptitle(f'Normal stress differences over the hold  |  $\\gamma = {sig(L["gamma"], 4)}$  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    for ax, (a, b, name, Np, Nt) in zip(axes, (('xx', 'yy', r'$N_1=\sigma_{xx}-\sigma_{yy}$', L['N1_p'], L['N1_t']),
                                               ('yy', 'zz', r'$N_2=\sigma_{yy}-\sigma_{zz}$', L['N2_p'], L['N2_t']))):
        ax.plot(L['ts'], tp[:, CI[a]] - tp[:, CI[b]], '-o', ms=4, color=WONG['blue'], lw=2, label=f'polymer: plateau {fmt_val_unc(Np[0], 0.5 * (Np[2] - Np[1]))}')
        ax.plot(L['ts'], tt[:, CI[a]] - tt[:, CI[b]], ':s', ms=4, color=WONG['reddishpurple'], lw=2, label=f'total: plateau {fmt_val_unc(Nt[0], 0.5 * (Nt[2] - Nt[1]))}')
        ax.axvspan(L['halt_ts'], float(L['ts'][-1]), color=WONG['green'], alpha=0.10, label='plateau window')
        ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
        ax.set_xlabel('time step'); ax.set_ylabel(name + '  (bulk, LJ)'); ax.set_title(name); ax.grid(alpha=0.3)
        smart_legend(ax, fontsize=12)
    return _save(fig, cfg, 'normal_stress_differences', L['lvl'])


def fig_Dc(cfg, R, L):
    """Transverse relaxation fit of u_x(zhat, t)/L: data (left) and data + model (right)."""
    F = L.get('Dc')
    if F is None:
        print('D_c figure skipped (no fit)')
        return None
    zff = np.linspace(0, 1, 400)
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(18, 7), constrained_layout=True)
    norm = Normalize(vmin=F['ts'][F['early'][0]], vmax=F['ts'][F['early'][-1]])
    cmap = plt.cm.viridis
    for i in F['early']:
        c = cmap(norm(F['ts'][i]))
        axl.plot(F['zf'], F['uhat'][i][F['idx']], 'o-', color=c, ms=3, alpha=0.6)
        axr.plot(F['zf'], F['uhat'][i][F['idx']], 'o', color=c, ms=3, alpha=0.35)
        axr.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=2.0)
    for ax in (axl, axr):
        ax._tri_has_colorbar = True
        ax.plot(zff, F['u_inf'](zff), 'k:', lw=1.8, label=r'fitted $t\to\infty$ state')
        ax.axhline(0, color='0.5', lw=1, alpha=0.6)
        ax.set(xlabel=r'$\hat z=(z-z_{\rm bot})/L_{\rm plates}$', ylabel=r'$u_x/L_{\rm plates}$', xlim=(0, 1))
        ax.grid(alpha=0.3)
        smart_legend(ax, fontsize=11)
    axl.set_title(r'$u_x(\hat z,t)/L$ -- hold snapshots (hold-referenced)', fontsize=15)
    axr.set_title(rf"Even-sine relaxation fit: $D_c={sig(F['Dc'])}\ \sigma^2/\tau$, $R^2={sig(F['R2'])}$"
                  + ('  (AT BOUND)' if F['at_bound'] else ''), fontsize=15)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    fig.colorbar(sm, ax=[axl, axr], fraction=0.015, pad=0.04).set_label('timestep')
    fig.suptitle(f'Cooperative diffusivity (shear) fit, level _g{L["lvl"]}  |  {cfg.sim_name}  |  plates frozen, '
                 f'$L_{{\\rm plates}}={sig(F["L"])}\\,\\sigma$', fontsize=12, fontweight='bold')
    return _save(fig, cfg, 'Dc_shear_fit', L['lvl'])


def fig_kappa(cfg, R, L):
    """kappa = D_c/G for the network and the series G (CI propagated from G)."""
    K = L.get('kappa')
    if not K:
        print('kappa figure skipped (no D_c or no G)')
        return None
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    for k, (key, lab, col, mk) in enumerate((('net', r'$D_c/G_\mathrm{network}$', WONG['blue'], 'o'),
                                             ('ser', r'$D_c/G_\mathrm{series}$', WONG['vermillion'], 's'))):
        v = K.get(key)
        if v is None:
            continue
        ax.errorbar([k], [v['k']], yerr=[[v['k'] - v['lo']], [v['hi'] - v['k']]], fmt=mk, ms=13, color=col, capsize=8, lw=2.5,
                    label=f"{lab} = {sig(v['k'])}")
    ax.set_xticks([0, 1]); ax.set_xticklabels(['network $G$', 'series $G$'], fontsize=16)
    ax.set_ylabel(r'$\kappa = D_c/G$  (LJ: $\sigma^5/(\epsilon\,\tau)$)')
    ax.set_title(f'Hydraulic permeability $\\kappa = D_c/G$\n$D_c = {sig(L["Dc"]["Dc"])}$  |  $\\gamma = {sig(L["gamma"], 4)}$', fontsize=15)
    ax.set_xlim(-0.5, 1.5); ax.grid(axis='y', alpha=0.3)
    smart_legend(ax, fontsize=13)
    return _save(fig, cfg, 'kappa', L['lvl'])


def fig_thermo_pressure(cfg, R, L):
    """P_th = -(1/3) tr(sigma^t): (a) bulk value vs step with the reference and P_bath,
    (b) profile evolution over the hold -- the P* = 1.5 check for the shear box (the
    thermo `press` of the deck is a box average over vacuum and is NOT this)."""
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(19, 6.5), constrained_layout=True)
    fig.suptitle(f'Thermodynamic pressure $P_{{\\rm th}}=-\\frac{{1}}{{3}}\\mathrm{{tr}}(\\sigma^t)$ in the bulk  |  $\\gamma = {sig(L["gamma"], 4)}$  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    axA.plot(L['ts'], L['Pth_ts'], '-o', ms=4, lw=2, color=WONG['blue'], label=f"bulk $P_{{\\rm th}}$: plateau {fmt_val_unc(L['Pth'][0], 0.5 * (L['Pth'][2] - L['Pth'][1]))}")
    axA.axvspan(L['halt_ts'], float(L['ts'][-1]), color=WONG['green'], alpha=0.10, label='plateau window')
    axA.axhline(cfg.P_BARO, color='k', ls=':', lw=1.4, label=f'$P_{{\\rm bath}} = {sig(cfg.P_BARO)}$ (compression runs)')
    Rf = R.get('ref')
    if Rf is not None:
        axA.axhline(Rf['Pth'][0], color='0.4', ls='--', lw=1.4, label=f"reference $\\gamma=0$: {sig(Rf['Pth'][0])}")
    axA.set_xlabel('time step'); axA.set_ylabel(r'$P_{\rm th}$  (LJ)'); axA.set_title('(a) bulk value over the hold'); axA.grid(alpha=0.3)
    smart_legend(axA, fontsize=12)
    if L.get('Pth_prof') is not None:
        ts, ev = post_hold(cfg, L['prof_p']['ts'], L['Pth_prof'])
        ref = None
        if Rf is not None and Rf.get('prof_tot_mean') is not None:
            pm = Rf['prof_tot_mean']
            ref = tuple(-(pm['xx'][j] + pm['yy'][j] + pm['zz'][j]) / 3.0 for j in range(3))
            ref = (ref[0], np.minimum(ref[1], ref[2]), np.maximum(ref[1], ref[2]))
        plot_evolution(axB, cfg, R, ts, ev, r'$P_{\rm th}(z,t)$  (LJ)', '(b) profile: reference -> hold -> plateau', ref=ref, legend=False)
        axB.axhline(cfg.P_BARO, color='k', ls=':', lw=1.2, alpha=0.6)
        robust_ylim(axB, list(ev) + ([ref[0]] if ref is not None else []), zmask=R['in_bulk'], pad=0.45, include_zero=False)
        h, lb = axB.get_legend_handles_labels()
        note = 'plateau mean, interior = ' + fmt_mu(np.nanmean(L['Pth_prof'][L['prof_plat']], axis=0)[R['interior']])
        h.append(Patch(alpha=0, label=note)); lb.append(note)
        smart_legend(axB, handles=h, labels=lb, fontsize=12)
    else:
        axB.text(0.5, 0.5, 'solvent profile\nnot found', ha='center', va='center', transform=axB.transAxes)
    return _save(fig, cfg, 'thermo_pressure', L['lvl'])


# ===========================================================================
#  8. FIGURES -- SWEEP
# ===========================================================================
def _final_overlay(ax, cfg, R, levels, get_final, ref, ylabel, title, include_zero=True):
    b = R['in_bulk']
    finals = []
    if ref is not None:
        m, lo, hi = ref
        ax.fill_between(R['zf'][b], lo[b], hi[b], color='0.5', alpha=0.2, lw=0, zorder=1)
        ax.plot(R['zf'][b], m[b], '--', color='k', lw=2.0, alpha=0.9, zorder=2)
        finals.append(m)
    for i, L in enumerate(levels):
        m = get_final(L)
        if m is None:
            continue
        ax.plot(R['zf'][b], m[b], '-', color=level_color(i), lw=2.8, zorder=4)
        finals.append(m)
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.4)
    shade_bulk(ax, R); mark_plates(ax, R)
    finish_axes(ax, ylabel, title)
    robust_ylim(ax, finals, zmask=b, pad=0.25, include_zero=include_zero)


def fig_network_stress_sweep(cfg, R, levels):
    """Final plateau sigma_p,xz(z) and sigma_p,xx(z) of every level + the reference."""
    hp = [L for L in levels if L.get('prof_p') is not None]
    if not hp:
        print('network-stress sweep figure skipped (no profiles)')
        return None
    fig, axes = plt.subplots(1, 2, figsize=(17, 6.5), constrained_layout=True)
    fig.suptitle(f'Network shear stress, final state of every level  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    Rf = R.get('ref')
    for ax, comp, title in ((axes[0], 'xz', r'(a) polymer $\sigma_{p,xz}(z)$'), (axes[1], 'xx', r'(b) polymer $\sigma_{p,xx}(z)$')):
        ref = Rf['prof_mean'][comp] if Rf is not None else None
        _final_overlay(ax, cfg, R, hp, lambda L, c=comp: L['prof_p_plat'][c], ref, r'$\sigma_{p,%s}(z)$' % comp, title)
        smart_legend(ax, handles=level_handles(hp, ref=Rf is not None), fontsize=11)
    return _save(fig, cfg, 'sweep_network_stress')


def fig_total_stress_sweep(cfg, R, levels):
    """Final plateau total stress of every level: sigma^t_xz and -sigma^t_ii / P_bath."""
    hp = [L for L in levels if L.get('prof_t') is not None]
    if not hp:
        print('total-stress sweep figure skipped (no solvent profiles)')
        return None
    comps = ('xz', 'xx', 'yy', 'zz')
    fig, axes = plt.subplots(1, 4, figsize=(30, 6.5), constrained_layout=True)
    fig.suptitle(f'Total stress, final state of every level (normal components / $P_{{\\rm bath}}={sig(cfg.P_BARO)}$)  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    Rf = R.get('ref')
    for k, (ax, comp) in enumerate(zip(axes, comps)):
        s = 1.0 if comp == 'xz' else -1.0 / cfg.P_BARO
        ref = None
        if Rf is not None and Rf.get('prof_tot_mean') is not None:
            m, lo, hi = Rf['prof_tot_mean'][comp]
            ref = (s * m, np.minimum(s * lo, s * hi), np.maximum(s * lo, s * hi))
        lab = (r'$\sigma^t_{xz}$' if comp == 'xz' else r'$-\sigma^t_{%s}/P_{\rm bath}$' % comp)
        _final_overlay(ax, cfg, R, hp, lambda L, c=comp, s=s: s * L['prof_t_plat'][c], ref, lab + '$(z)$',
                       f'({"abcd"[k]}) total {lab}', include_zero=(comp == 'xz'))
        if comp != 'xz':
            ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6)
        smart_legend(ax, handles=level_handles(hp, ref=ref is not None), fontsize=11)
    return _save(fig, cfg, 'sweep_total_stress')


def fig_series_sweep(cfg, R, levels):
    """sigma_p,xz(t) of every level from the fine series (drive + hold), plateau dotted."""
    hs = [L for L in levels if L.get('series') is not None]
    if not hs:
        print('series sweep figure skipped')
        return None
    fig, ax = plt.subplots(figsize=(13, 6.5), constrained_layout=True)
    for i, L in enumerate(hs):
        S = L['series']
        ax.plot(S['step'], rolling_mean(S['sp_xz'], cfg.roll_win), '-', color=level_color(i), lw=2.0,
                label=fr"$\gamma={sig(L['gamma'], 4)}$" + (f":  plateau {sig(L['PF']['mean'], 3)}" if 'PF' in L else ''))
        if 'PF' in L:
            ax.hlines(L['PF']['mean'], L['PF']['step0'], S['step'][-1], colors=level_color(i), linestyles=':', lw=1.5)
    Rf = R.get('ref')
    if Rf is not None:
        ax.axhline(Rf['S_ref'], color='k', ls='--', lw=1.2, alpha=0.7, label=f"reference $\\gamma=0$: {sig(Rf['S_ref'], 2)}")
    ax.axhline(0, color='k', ls='--', lw=0.8, alpha=0.4)
    ax.set_xlabel('time step'); ax.set_ylabel(r'$\langle\sigma_{p,xz}\rangle_{\rm bulk}$  (LJ)')
    ax.set_title(f'Bulk polymer shear stress histories (rolling mean over {cfg.roll_win} blocks; dotted = plateau)  |  {cfg.sim_name}', fontsize=13)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=11)
    return _save(fig, cfg, 'sweep_shear_stress_history')


def fig_G_sweep(cfg, R, levels):
    """G vs held strain, two panes as fig_G: (a) increments with the absolute values hollow
    and the total-stress check, (b) the network and series estimates alone."""
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    sub = all(L['G_ref'] == 'measured' for L in levels)
    fig.suptitle('Shear modulus across the sweep' + (' (increment from $\\gamma=0$)' if sub else ' (absolute)') + f'   |   {cfg.sim_name}',
                 fontsize=13, fontweight='bold')

    def _ser(ax, key, marker, color, label, hollow=False, dx=0.0, val='G'):
        Ls = [L for L in levels if key in L['G']]
        if not Ls:
            return
        g = np.array([L['gamma'] for L in Ls])
        v = np.array([L['G'][key][val if val == 'G' else 'abs'] for L in Ls])
        lo = np.array([L['G'][key]['lo' if val == 'G' else 'abs_lo'] for L in Ls])
        hi = np.array([L['G'][key]['hi' if val == 'G' else 'abs_hi'] for L in Ls])
        kw = dict(mfc='none', alpha=0.7, lw=1.5, ms=8, capsize=4, ls=':') if hollow else dict(lw=2, ms=10, capsize=6, ls='-')
        ax.errorbar(g + dx, v, yerr=[v - lo, hi - v], fmt=marker, color=color, label=label, **kw)

    for ax, pane in ((axA, 'a'), (axB, 'b')):
        _ser(ax, 'net', 'o', WONG['blue'], r'network  $\Delta\langle\sigma_{p,xz}\rangle_{\rm int}/\gamma$')
        _ser(ax, 'ser', 's', WONG['vermillion'], r'series  $\Delta\langle\sigma_{p,xz}\rangle_{\rm plateau}/\gamma$')
        if pane == 'a':
            if sub:
                _ser(ax, 'net', 'o', WONG['blue'], 'network, absolute', hollow=True, dx=0.002, val='abs')
                _ser(ax, 'ser', 's', WONG['vermillion'], 'series, absolute', hollow=True, dx=-0.002, val='abs')
            Ls = [L for L in levels if 'tot' in L['G']]
            if Ls:
                g = np.array([L['gamma'] for L in Ls]); v = np.array([L['G']['tot']['G'] for L in Ls])
                ax.errorbar(g, v, yerr=[v - [L['G']['tot']['lo'] for L in Ls], [L['G']['tot']['hi'] for L in Ls] - v],
                            fmt='D', mfc='none', color=WONG['green'], ms=8, lw=1.5, capsize=4, ls=':', label=r'total $\sigma^t_{xz}/\gamma$ (check)')
        ax.set_xlabel(r'held shear strain  $\gamma$'); ax.set_ylabel(r'$G$  (LJ units)'); ax.grid(alpha=0.3)
        ax.set_title('(a) with the absolute values (hollow) and the total-stress check' if pane == 'a' else '(b) shear modulus per level, two estimates', fontsize=13)
        smart_legend(ax, fontsize=11)
    return _save(fig, cfg, 'sweep_G')


def fig_stress_strain_sweep(cfg, R, levels):
    """The stress-strain curve itself: plateau sigma_p,xz (profile interior and series) vs
    held gamma, INCLUDING the gamma = 0 reading (hollow), with a least-squares line through
    each series (slope = offset-free G, as the compression notebooks' figure 7b) and the
    through-origin fit for comparison."""
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    Rf = R.get('ref')
    fits = []

    def _series(key, color, marker, name, ref_v=None, ref_h=None):
        Ls = [L for L in levels if key in L['G']]
        if not Ls:
            return
        g = np.array([L['gamma'] for L in Ls])
        s = np.array([L['G'][key]['sigma'] for L in Ls])
        err = np.array([0.5 * (L['G'][key].get('abs_hi', L['G'][key]['hi']) - L['G'][key].get('abs_lo', L['G'][key]['lo'])) * L['gamma'] for L in Ls])
        ax.errorbar(g, s, yerr=err, fmt=marker + '-', lw=2, ms=9, color=color, capsize=6, label=name)
        gg, ss = list(g), list(s)
        if ref_v is not None and np.isfinite(ref_v):
            ax.errorbar([0.0], [ref_v], yerr=ref_h, fmt=marker, ms=9, mfc='none', color=color, capsize=6)
            gg, ss = [0.0] + gg, [ref_v] + ss
        if len(gg) >= 2:
            slope, icpt = np.polyfit(gg, ss, 1)
            xs = np.linspace(0, max(gg) * 1.05, 20)
            ax.plot(xs, slope * xs + icpt, ls=':', lw=1.5, color=color, alpha=0.8)
            G0 = float(np.sum(np.array(gg) * np.array(ss)) / np.sum(np.array(gg) ** 2))
            fits.append((name.split()[0], slope, G0, len(gg)))

    _series('net', WONG['blue'], 'o', 'network $\\langle\\sigma_{p,xz}\\rangle_{\\rm int}$ (plateau)',
            Rf['sp_xz_int'] if Rf else None, Rf['sp_xz_int_half'] if Rf else None)
    _series('ser', WONG['vermillion'], 's', 'series $\\langle\\sigma_{p,xz}\\rangle$ (plateau)',
            Rf['S_ref'] if Rf else None, 0.5 * (Rf['S_ref_hi'] - Rf['S_ref_lo']) if Rf else None)
    ax.plot([], [], 'o', mfc='none', color='0.4', label=r'hollow = $\gamma=0$ reading')
    ax.set_title('Stress vs strain: least-squares slopes incl. $\\gamma=0$ (dotted) and through-origin:  '
                 + ',  '.join(f"$G_{{\\rm {n[:4]}}}\\approx{sig(s)}$ / {sig(g0)} ({k} pts)" for n, s, g0, k in fits) + '\n' + cfg.sim_name, fontsize=11)
    ax.set_xlabel(r'held shear strain  $\gamma$'); ax.set_ylabel(r'plateau $\sigma_{p,xz}$  (LJ)')
    ax.grid(alpha=0.3); ax.set_xlim(left=-0.005)
    smart_legend(ax, fontsize=12)
    L_fits = {n: (s, g0) for n, s, g0, k in fits}
    return _save(fig, cfg, 'sweep_stress_strain'), L_fits


def fig_normal_sweep(cfg, R, levels):
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    g = np.array([L['gamma'] for L in levels])
    for key, lab, col, mk in (('N1_p', r'$N_1$ polymer', WONG['blue'], 'o'), ('N2_p', r'$N_2$ polymer', WONG['vermillion'], 's'),
                              ('N1_t', r'$N_1$ total', WONG['blue'], 'o'), ('N2_t', r'$N_2$ total', WONG['vermillion'], 's')):
        v = np.array([L[key][0] for L in levels]); lo = np.array([L[key][1] for L in levels]); hi = np.array([L[key][2] for L in levels])
        hollow = key.endswith('_t')
        ax.errorbar(g, v, yerr=[v - lo, hi - v], fmt=mk + ('--' if hollow else '-'), ms=8, lw=1.5, color=col, capsize=5,
                    label=lab, **({'mfc': 'none'} if hollow else {}))
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    ax.set_xlabel(r'held shear strain  $\gamma$'); ax.set_ylabel('normal stress difference  (LJ)')
    ax.set_title('Normal stress differences vs strain (plateau means; both -> 0 for a linear isotropic network)', fontsize=12)
    ax.grid(alpha=0.3); smart_legend(ax, fontsize=11)
    return _save(fig, cfg, 'sweep_normal_stress')


def fig_Dc_sweep(cfg, R, levels):
    hd = [L for L in levels if L.get('Dc') is not None]
    if not hd:
        print('D_c sweep figure skipped (no fits)')
        return None
    n = len(hd)
    fig, axes = plt.subplots(1, n + 1, figsize=(7 + 6.5 * n, 6), constrained_layout=True)
    fig.suptitle(f'Cooperative diffusivity (shear) across the sweep  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    ax = axes[0]
    g = np.array([L['gamma'] for L in hd]); v = np.array([L['Dc']['Dc'] for L in hd])
    ax.plot(g, v, 'o-', color=WONG['reddishpurple'], lw=2, ms=8)
    for L in hd:
        ax.annotate(f"$R^2$={sig(L['Dc']['R2'])}", (L['gamma'], L['Dc']['Dc']), textcoords='offset points', xytext=(0, -12),
                    ha='center', va='top', fontsize=10, color='0.35', annotation_clip=True)
    ax.margins(x=0.12, y=0.15)
    ax.axhline(cfg.DC_SLOW_REF, color='0.4', ls=':', lw=1.5, label=f'slow reference $D_c$ = {sig(cfg.DC_SLOW_REF)}')
    ax.set_xlabel(r'held shear strain  $\gamma$'); ax.set_ylabel(r'$D_c$  ($\sigma^2/\tau$)'); ax.set_title(r'(a) $D_c$ vs strain', fontsize=15)
    ax.grid(alpha=0.3); smart_legend(ax, fontsize=11)
    zff = np.linspace(0, 1, 300)
    for k, L in enumerate(hd):
        ax = axes[k + 1]; F = L['Dc']
        cmap = plt.cm.viridis; norm = Normalize(vmin=F['ts'][F['early'][0]], vmax=F['ts'][F['early'][-1]])
        for i in F['early']:
            c = cmap(norm(F['ts'][i]))
            ax.plot(F['zf'], F['uhat'][i][F['idx']], 'o', color=c, ms=2.5, alpha=0.35)
            ax.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=1.5)
        ax.set(xlabel=r'$\hat z$', ylabel=r'$u_x/L$', xlim=(0, 1))
        ax.set_title(fr"({'bcdefgh'[k]}) $\gamma={sig(L['gamma'], 3)}$:  $D_c={sig(F['Dc'])}$, $R^2={sig(F['R2'])}$", fontsize=14,
                     color=level_color(levels.index(L)))
        ax.grid(alpha=0.3)
    return _save(fig, cfg, 'sweep_Dc')


def fig_kappa_sweep(cfg, R, levels):
    hk = [L for L in levels if L.get('kappa')]
    if not hk:
        print('kappa sweep figure skipped')
        return None
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    for key, lab, col, mk in (('net', r'$D_c/G_\mathrm{network}$', WONG['blue'], 'o'), ('ser', r'$D_c/G_\mathrm{series}$', WONG['vermillion'], 's')):
        Ls = [L for L in hk if key in L['kappa']]
        if not Ls:
            continue
        g = np.array([L['gamma'] for L in Ls]); k = np.array([L['kappa'][key]['k'] for L in Ls])
        ax.errorbar(g, k, yerr=[k - [L['kappa'][key]['lo'] for L in Ls], [L['kappa'][key]['hi'] for L in Ls] - k],
                    fmt=mk + '-', ms=10, lw=2, color=col, capsize=6, label=lab)
    ax.set_xlabel(r'held shear strain  $\gamma$'); ax.set_ylabel(r'$\kappa = D_c/G$  (LJ)')
    ax.set_title(r'Hydraulic permeability $\kappa = D_c/G$ vs strain', fontsize=15); ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=12)
    return _save(fig, cfg, 'sweep_kappa')


def fig_thermo_pressure_sweep(cfg, R, levels):
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    for i, L in enumerate(levels):
        ax.plot(L['ts'], L['Pth_ts'], '-o', ms=3, lw=1.8, color=level_color(i),
                label=fr"$\gamma={sig(L['gamma'], 4)}$: plateau {fmt_val_unc(L['Pth'][0], 0.5 * (L['Pth'][2] - L['Pth'][1]))}")
    ax.axhline(cfg.P_BARO, color='k', ls=':', lw=1.4, label=f'$P_{{\\rm bath}} = {sig(cfg.P_BARO)}$')
    Rf = R.get('ref')
    if Rf is not None:
        ax.axhline(Rf['Pth'][0], color='0.4', ls='--', lw=1.4, label=f"reference $\\gamma=0$: {sig(Rf['Pth'][0])}")
    ax.set_xlabel('time step'); ax.set_ylabel(r'$P_{\rm th}$  (bulk, LJ)')
    ax.set_title(f'Thermodynamic pressure $P_{{\\rm th}}=-\\frac{{1}}{{3}}\\mathrm{{tr}}(\\sigma^t)$ over every hold  |  {cfg.sim_name}', fontsize=13)
    ax.grid(alpha=0.3); smart_legend(ax, fontsize=11)
    return _save(fig, cfg, 'sweep_thermo_pressure')


# ===========================================================================
#  9. SUMMARIES
# ===========================================================================
def print_hold_check(cfg, levels):
    hd = [L for L in levels if L.get('Dc') is not None]
    if not hd:
        return
    print(f'\nHOLD-ADEQUACY CHECK  (tau_1 = L^2/(4 pi^2 D_c) for the slowest antisymmetric mode; residual = mean excess '
          f'stress over the last {cfg.plateau_frac:.0%} of the hold, i.e. the window G is read from)')
    for L in hd:
        F = L['Dc']
        print(f"  level _g{L['lvl']}:  L = {F['L']:.1f} sigma   hold T = {F['hold_T']:.0f} tau = {F['hold_T'] / cfg.dt_lj / 1e6:.2f}M steps")
        for tag, h in F['hold_check'].items():
            flag = '' if h['ok'] else '   <-- TOO SHORT'
            print(f"     {tag:<4s} D_c={h['Dc']:.3f}:  tau_1 = {h['tau1']:.0f} tau = {h['tau1'] / cfg.dt_lj / 1e6:.2f}M steps"
                  f"  |  held {F['hold_T'] / h['tau1']:.2f} tau_1  ->  residual {h['avg'] * 100:5.2f}% (end {h['end'] * 100:5.2f}%)"
                  f"  |  {cfg.DC_TARGET_RESID:.0%} needs {h['need']:.1f} tau_1 = {h['need'] * h['tau1'] / cfg.dt_lj / 1e6:.1f}M steps{flag}")
    print('  (shear_slab.lmp holds a FLAT nsteps per level; size NSTEPS in the batch from this check.)')


def print_summary(cfg, levels):
    ci = int(cfg.ci_level * 100)
    how = 'G = increment from the gamma = 0 reference' if cfg.G_SUBTRACT_REF else 'G = absolute stress / gamma'
    if cfg.G_SUBTRACT_REF and any(L['G_ref'] == 'absent' for L in levels):
        how += '; ABSOLUTE for levels without _ref files'
    print(f'\nSUMMARY  ({cfg.sim_name}; {ci}% CIs; {how})')
    print(f"{'target':>7s} {'gamma':>8s} {'G_net':>21s} {'G_ser':>21s} {'G_tot':>21s} {'P_th':>7s} {'N1/sxz':>7s} {'D_c':>10s} {'kappa_net':>10s}")
    for L in levels:
        def ci_(d):
            return f"{d['G']:.4f} [{d['lo']:.4f},{d['hi']:.4f}]" if d else 'n/a'
        sxz = abs(np.nanmean(L['sxz_p_ts'][L['plat']]))
        dc = f"{L['Dc']['Dc']:.3e}" if L.get('Dc') else 'n/a'
        kn = f"{L['kappa']['net']['k']:.3e}" if L.get('kappa', {}).get('net') else 'n/a'
        print(f"{L['gamma_target']:7.4f} {L['gamma']:8.5f} {ci_(L['G'].get('net')):>21s} {ci_(L['G'].get('ser')):>21s} "
              f"{ci_(L['G'].get('tot')):>21s} {L['Pth'][0]:7.4f} {L['N1_p'][0] / max(sxz, 1e-30):+7.3f} {dc:>10s} {kn:>10s}")
    print_hold_check(cfg, levels)
