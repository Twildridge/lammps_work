"""triaxial.py -- shared analysis library for the triaxial-compression notebooks.

Used by  scripts/triaxial_compression_single.ipynb  (one strain level)  and
         scripts/triaxial_compression_sweep.ipynb   (every level of a sweep).

Everything the two notebooks need lives here so the definitions cannot drift
between them: file readers, the Terzaghi network/pore split, the plateau
(block-bootstrap) estimators behind M, the G estimate from the lateral network
stress, the consolidation fit behind D_c, the lambda-calibrated Voronoi volume
fraction (via lib/volfrac.py) and every figure.  The notebooks only set a
`Config`, call `load_reference` / `load_level` / `add_volume_fractions`, and
then call one `fig_*` function per figure.

Layout of this file
    0. style + palette              setup_style, WONG, level_color
    1. Config                       all knobs, path builders, NSTEPS auto-detect
    2. file readers                 fix ave/time, ave/chunk, print, dumps
    3. statistics                   mean_ci, block bootstrap, plateau_window
    4. reference state              load_reference   (eps = 0, shared by all levels)
    5. one strain level             load_level       (stress, Terzaghi, piston, M, G, D_c, kappa)
    6. volume fractions             add_volume_fractions (mass-fraction / Voronoi / lambda-calibrated)
    7. Expanse sync                 sync_from_expanse
    8. plotting primitives          zn, shade_gel, mark_walls, plot_evolution, ...
    9. figures, single level        fig_strain ... fig_kappa
   10. figures, sweep               fig_*_sweep
   11. summaries                    print_summary, print_hold_check

Physics conventions (see the Notes section at the end of either notebook):
  * total stress sigma^t = sigma_p + sigma_s (group stress/atom, kinetic term included)
  * Terzaghi: sigma' = sigma^t - p_pore, p_pore read per curve from the flat
    reservoir at z/Lz ~ baseline_zf (one scalar per stress component)
  * uniaxial strain (fixed lateral box):  sigma'_zz = M eps,  sigma'_xx = lambda eps,
    lambda = M - 2G  ->  sigma'_zz/sigma'_xx = M/(M - 2G)  ->  G = (sigma'_zz - sigma'_xx)/(2 eps)
  * M_network = <sigma'_zz>_membrane / eps_applied ;  M_piston = <F_z/A>_plateau / eps_applied
  * D_c from the two-sided consolidation fit of the polymer displacement u_z(z,t)
  * kappa = D_c / M   (Darcy permeability over viscosity, k/eta, in LJ units)
"""
from __future__ import annotations

import re
import sys
import stat
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import volfrac  # noqa: E402  (scripts/lib/volfrac.py -- Voronoi + lambda calibration)


# ===========================================================================
#  0. STYLE
# ===========================================================================
WONG = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'vermillion': '#D55E00',
        'skyblue': '#56B4E9', 'yellow': '#F0E442', 'reddishpurple': '#CC79A7', 'black': '#000000'}
EVO_CMAP = 'cividis'                                   # CVD-safe sequential map for time
GEL_SHADE = dict(color='0.6', alpha=0.15, zorder=0)
LEVEL_COLORS = [WONG[k] for k in ('blue', 'vermillion', 'green', 'reddishpurple',
                                  'orange', 'skyblue', 'yellow', 'black')]
COMP_LABEL = {'zz': r'zz', 'xx': r'xx', 'yy': r'yy'}
COMPONENTS = ('zz', 'xx', 'yy')
_Z95 = 1.959963985


def level_color(i):
    return LEVEL_COLORS[i % len(LEVEL_COLORS)]


def setup_style():
    """rcParams shared by every notebook figure (CMU Serif, large labels)."""
    plt.rcParams.update({
        'font.family': 'CMU Serif', 'mathtext.fontset': 'cm', 'mathtext.rm': 'CMU Serif',
        'font.size': 20, 'axes.titlesize': 22, 'axes.labelsize': 25,
        'xtick.labelsize': 23, 'ytick.labelsize': 23, 'legend.fontsize': 23,
        'figure.titlesize': 22, 'axes.unicode_minus': False,
    })


# ===========================================================================
#  1. CONFIG
# ===========================================================================
@dataclass
class Config:
    """Every knob of the analysis.  The notebook Config cell builds one of these;
    nothing else in the notebook needs editing to switch runs."""
    # ---- which run -------------------------------------------------------
    DATANAME: str
    INTERACTION: str                      # "epsSS_epsSP"
    RUN_ID: str                           # local folder under flow_data_local/{compression,plots}
    COMP_LEVELS: list                     # applied-strain targets, as STRINGS ("0.10"), = STRAIN_TARGETS in the .batch
    NSTEPS: int = None                    # None -> auto-detect the largest NSTEPS present in DATA_DIR
    base_dir: str = '../../flow_data_local'
    # ---- profile / window knobs -----------------------------------------
    binWidth: float = 2.0                 # coarse z-bin (sigma); must match triaxial_compression.lmp
    n_curves: int = 10                    # evolution curves drawn per profile
    Ncount_min: int = 200                 # min polymer atoms per z-bin for the D_c fit domain
    dt_lj: float = 0.005                  # LJ timestep (D_c time axis)
    ci_level: float = 0.95
    plateau_frac: float = 0.25            # FIXED trailing fraction of the hold for profile plateau means
    plateau_frac_auto: float = 0.45       # LONGEST candidate trailing window for plateau_window (piston -> M_piston)
    gel_thresh: float = 0.05              # gel/membrane = bins where |sigma_p,zz| > gel_thresh * max
    flat_tol: float = 0.15                # "flat inside gel" if relative linear trend < this
    wall_margin: float = 4.0              # sigma trimmed off both gel ends for interior means
    baseline_zf: float = 0.95             # reservoir baseline window centre (z/Lz) for p_pore
    baseline_zf_half: float = 0.04        # its half-width
    roll_win: int = 21                    # rolling-mean window (samples) for the piston plots
    # ---- volume fractions (lib/volfrac.py) -------------------------------
    VOR_ENABLE: bool = True               # False -> mass-fraction only (no trajectory pass)
    VOR_MOBILE_ONLY: bool = True          # tessellate types 1,2,3 only
    VOR_NORM: str = 'bin'                 # 'bin' (absolute) | 'mobile' (saturating) for the raw phi^vor
    P_CAL: float = 1.5                    # P_local for lambda(phi_p, P); reference-state P_target
    PHI_FLOOR: float = 0.02
    REF_VOR_FRAMES: int = 3               # reference frames tessellated (~20 s each)
    VOR_MAX_FRAMES: int = 4               # plateau frames tessellated per level (~20 s each)
    # ---- G from the lateral network stress -------------------------------
    G_SUBTRACT_REF: bool = True           # use increments relative to the eps = 0 reference state
    # ---- D_c consolidation fit ------------------------------------------
    DC_N_MODES: int = 5
    DC_KMAX_IC: int = 199
    DC_FRAC_EARLY: float = 1.0
    DC_TRIM_BINS: int = 2
    DC_FREE_AMPS: bool = True
    DC_BOUNDS: tuple = (1e-6, 1.0)
    DC_SLOW_REF: float = 0.17             # sigma^2/tau, slow collective D_c (free-swelling equilibration)
    DC_TARGET_RESID: float = 0.01
    # ---- Expanse ---------------------------------------------------------
    EXPANSE_HOST: str = 'login.expanse.sdsc.edu'
    EXPANSE_USER: str = 'dpollard'
    RUNS_ROOT: str = '/home/dpollard/Documents/lammps_runs/triaxial_compression'
    TRAJ_ROOT: str = '/expanse/lustre/scratch/dpollard/temp_project/lammps_trajectories'
    # ---- derived (filled by __post_init__) -------------------------------
    DATA_DIR: Path = field(init=False)
    PLOT_DIR: Path = field(init=False)
    TRAJ_DIR: Path = field(init=False)
    sim_name: str = field(init=False)

    def __post_init__(self):
        assert self.COMP_LEVELS, 'COMP_LEVELS is empty -- list at least one level, e.g. ["0.10"]'
        self.COMP_LEVELS = [str(l) for l in self.COMP_LEVELS]
        base = Path(self.base_dir)
        self.DATA_DIR = base / 'compression' / self.RUN_ID
        self.PLOT_DIR = base / 'plots' / 'compression' / self.RUN_ID
        self.TRAJ_DIR = base / 'traj_files.nosync'
        for d in (self.DATA_DIR, self.PLOT_DIR, self.TRAJ_DIR):
            d.mkdir(parents=True, exist_ok=True)
        if self.NSTEPS is None:
            pat = re.compile(rf'^sigmazz_polymer(?:_ref)?_{re.escape(self.DATANAME)}_'
                             rf'{re.escape(self.INTERACTION)}_(\d+)(?:_c[\d.]+)?\.dat$')
            hits = sorted({int(m.group(1)) for f in self.DATA_DIR.glob('sigmazz_polymer*.dat')
                           if (m := pat.match(f.name))})
            if not hits:
                raise ValueError('NSTEPS is None and no sigmazz_polymer_*.dat in DATA_DIR yet '
                                 '-- set NSTEPS explicitly, or run sync_from_expanse first.')
            self.NSTEPS = hits[-1]
            print(f'auto-detected NSTEPS = {self.NSTEPS} from {len(hits)} candidate(s): {hits}')
        self.sim_name = f'{self.DATANAME}_{self.INTERACTION}_{self.NSTEPS}'
        print(f'Config: {self.sim_name}  |  levels {self.COMP_LEVELS}\n'
              f'  data  {self.DATA_DIR}\n  plots {self.PLOT_DIR}\n  traj  {self.TRAJ_DIR}')

    # ---- path builders: *_ref files carry no level tag, production files do ----
    def tag(self, lvl):
        return '' if lvl is None else f'_c{lvl}'

    def path(self, name, lvl=None, ext='dat'):
        return self.DATA_DIR / f'{name}_{self.sim_name}{self.tag(lvl)}.{ext}'

    def traj(self, name, lvl=None):
        return self.TRAJ_DIR / f'{name}_{self.sim_name}{self.tag(lvl)}.lammpstrj'

    def plot(self, stem, lvl=None):
        return self.PLOT_DIR / f'{stem}_{self.sim_name}{self.tag(lvl)}.png'


# ===========================================================================
#  2. FILE READERS
# ===========================================================================
def read_print_file(filepath, col_names=None):
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            rows.append([float(v) for v in line.split()])
    arr = np.array(rows)
    if col_names is None:
        col_names = [f'col_{i}' for i in range(arr.shape[1])]
    return {name: arr[:, i] for i, name in enumerate(col_names)}


def read_ave_time_file(filepath):
    """fix ave/time mode vector -> list of (timestep, bin_idx, values)."""
    out = []
    with open(filepath) as f:
        lines = [l for l in f if not l.startswith('#') and l.strip()]
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) == 2:
            ts, nrows = int(parts[0]), int(parts[1])
            vals = []
            for j in range(1, nrows + 1):
                if i + j < len(lines):
                    vp = lines[i + j].split()
                    if len(vp) == 2:
                        vals.append(float(vp[1]))
            if vals:
                out.append((ts, np.arange(1, len(vals) + 1), np.array(vals)))
            i += nrows + 1
        else:
            i += 1
    return out


def read_ave_chunk_file(filepath):
    """fix ave/chunk -> list of (timestep, array[rows, cols]); cols [chunk, Coord1, Ncount, val...]."""
    snaps = []
    with open(filepath) as f:
        lines = [l for l in f if l.strip() and not l.startswith('#')]
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) in (2, 3):
            try:
                ts, nch = int(parts[0]), int(parts[1])
            except ValueError:
                i += 1
                continue
            rows = []
            for j in range(1, nch + 1):
                if i + j < len(lines):
                    rows.append([float(v) for v in lines[i + j].split()])
            if rows:
                snaps.append((ts, np.array(rows)))
            i += nch + 1
        else:
            i += 1
    return snaps


def read_strain_file(filepath):
    """strain_zz: step L_initial L_current -> (steps, eps = (L0 - L)/L0)."""
    arr = np.atleast_2d(np.loadtxt(filepath, comments='#'))
    ts = arr[:, 0].astype(int)
    L0, L = arr[:, 1], arr[:, 2]
    return ts, (L0 - L) / L0


def load2c(path, mincols=2):
    """np.loadtxt as a 2-D array, or None if the file is missing/too narrow."""
    try:
        a = np.atleast_2d(np.loadtxt(path, comments='#'))
    except Exception:
        return None
    return a if (a.ndim == 2 and a.shape[1] >= mincols and a.size) else None


def read_box(dumpfile):
    """Box bounds + edge lengths from the first frame of any LAMMPS dump."""
    with open(dumpfile) as f:
        head = [next(f) for _ in range(9)]
    b = {k: tuple(map(float, head[5 + i].split()[:2])) for i, k in enumerate('xyz')}
    b.update(lx=b['x'][1] - b['x'][0], ly=b['y'][1] - b['y'][0], lz=b['z'][1] - b['z'][0])
    return b


def wall_z_first_frame(traj, types=(4, 5)):
    """Mean z of each atom type in the FIRST frame of a dump trajectory."""
    zc = {t: [] for t in types}
    if not Path(traj).exists():
        return {t: np.nan for t in types}
    with open(traj) as f:
        cols = None
        for line in f:
            if line.startswith('ITEM: ATOMS'):
                cols = line.split()[2:]
                break
        ti, zi = cols.index('type'), cols.index('z')
        for line in f:
            if line.startswith('ITEM:'):
                break
            p = line.split()
            t = int(float(p[ti]))
            if t in zc:
                zc[t].append(float(p[zi]))
    return {t: (np.mean(v) if v else np.nan) for t, v in zc.items()}


# ===========================================================================
#  3. STATISTICS
# ===========================================================================
def mean_ci(stack, ci=0.95):
    """(mean, lo, hi) per column via a t-interval; all-NaN columns stay NaN quietly."""
    stack = np.asarray(stack, float)
    n = stack.shape[0]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        m = np.nanmean(stack, axis=0)
        if n < 2:
            return m, m, m
        se = stats.sem(stack, axis=0, nan_policy='omit')
    half = np.asarray(se) * stats.t.ppf(0.5 + ci / 2, df=n - 1)
    return m, m - half, m + half


def rolling_mean(y, win):
    """Centered moving average with a shrinking edge window (plotting only)."""
    y = np.asarray(y, float)
    n = len(y)
    if win <= 1 or n == 0:
        return y.copy()
    half = win // 2
    csum = np.concatenate(([0.0], np.cumsum(y)))
    out = np.empty(n)
    for k in range(n):
        lo, hi = max(0, k - half), min(n, k + half + 1)
        out[k] = (csum[hi] - csum[lo]) / (hi - lo)
    return out


def subsample(ts, stack, k):
    """Evenly pick up to k snapshots (keep order; always include the last)."""
    ts = np.asarray(ts)
    stack = np.asarray(stack)
    n = len(ts)
    idx = np.arange(n) if n <= k else np.unique(np.linspace(0, n - 1, k).round().astype(int))
    return ts[idx], stack[idx]


def autocorr_time(x):
    """Integrated autocorrelation time (samples), Sokal automatic windowing."""
    x = np.asarray(x, float)
    n = len(x)
    if n < 4:
        return 1.0
    x = x - x.mean()
    var = np.dot(x, x) / n
    if var <= 0:
        return 1.0
    fx = np.fft.rfft(x, n=2 * n)
    acf = np.fft.irfft(fx * np.conj(fx))[:n].real / (var * n)
    tau = 1.0
    for W in range(1, n):
        tau = 1.0 + 2.0 * np.sum(acf[1:W + 1])
        if W >= 5.0 * tau:
            break
    return float(max(tau, 1.0))


def block_bootstrap_ci(x, ci=0.95, n_boot=2000, block=None, seed=12345):
    """Circular block-bootstrap CI for the MEAN of an autocorrelated series.
    Returns (mean, lo, hi, block_len, tau); blocks ~2*tau long preserve the
    within-block correlation so the CI reflects ~N/block effective samples."""
    x = np.asarray(x, float)
    n = len(x)
    if n < 2:
        m = float(x.mean()) if n else np.nan
        return m, np.nan, np.nan, 1, 1.0
    tau = autocorr_time(x)
    if block is None:
        block = int(np.clip(np.ceil(2.0 * tau), 1, max(1, n // 2)))
    rng = np.random.default_rng(seed)
    nblk = int(np.ceil(n / block))
    starts = rng.integers(0, n, size=(n_boot, nblk))
    idx = (starts[:, :, None] + np.arange(block)[None, None, :]) % n
    idx = idx.reshape(n_boot, -1)[:, :n]
    means = x[idx].mean(axis=1)
    lo, hi = np.percentile(means, [100 * (1 - ci) / 2, 100 * (1 + ci) / 2])
    return float(x.mean()), float(lo), float(hi), int(block), float(tau)


def plateau_window(steps, x, plateau_frac_auto=0.45, ci=0.95, fracs=None):
    """Longest DRIFT-FREE trailing window of (steps, x), block-bootstrapped.
    Candidates are plateau_frac_auto * (9,8,...,2)/9 of the series; a window
    passes when its two halves agree within their block-bootstrap CIs.  Returns
    dict(mean, lo, hi, step0, n, block, tau, frac[, warn]) -- the ONE plateau
    P every M_piston in the notebooks is read from."""
    steps = np.asarray(steps, float)
    x = np.asarray(x, float)
    if fracs is None:
        fracs = plateau_frac_auto * np.arange(9, 1, -1) / 9.0
    t0, t1 = float(steps[0]), float(steps[-1])
    for f in fracs:
        sel = steps >= t1 - f * (t1 - t0)
        xw, sw = x[sel], steps[sel]
        if len(xw) < 20:
            continue
        h = len(xw) // 2
        m1, lo1, hi1, *_ = block_bootstrap_ci(xw[:h], ci=ci)
        m2, lo2, hi2, *_ = block_bootstrap_ci(xw[h:], ci=ci)
        if (lo1 <= m2 <= hi1) or (lo2 <= m1 <= hi2):
            m, lo, hi, blk, tau = block_bootstrap_ci(xw, ci=ci)
            return dict(mean=m, lo=lo, hi=hi, step0=float(sw[0]), n=len(xw),
                        block=blk, tau=tau, frac=float(f))
    sel = steps >= t1 - fracs[-1] * (t1 - t0)
    m, lo, hi, blk, tau = block_bootstrap_ci(x[sel], ci=ci)
    return dict(mean=m, lo=lo, hi=hi, step0=float(steps[sel][0]), n=int(sel.sum()),
                block=blk, tau=tau, frac=float(fracs[-1]), warn=True)


def fmt_val_unc(v, u):
    """value +/- uncertainty, uncertainty rounded to 2 sig figs."""
    v = float(v)
    if np.isfinite(u) and u > 0:
        dec = int(np.clip(1 - np.floor(np.log10(u)), 0, 6))
        return f'{v:.{dec}f} ± {u:.{dec}f}'
    return f'{v:.3g}'


def fmt_mu(vals):
    """mean +/- std over finite values."""
    a = np.asarray(vals, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return 'n/a'
    return fmt_val_unc(np.mean(a), np.std(a))


# ===========================================================================
#  4. REFERENCE STATE  (eps = 0; shared by every level)
# ===========================================================================
def _component_stacks(cfg, comp, lvl=None, ref=False):
    """sigma<comp>_{polymer,solvent}[_ref] files -> dict(ts, bins, p, s, t) or None."""
    suffix = '_ref' if ref else ''
    fp = cfg.path(f'sigma{comp}_polymer{suffix}', None if ref else lvl)
    fs = cfg.path(f'sigma{comp}_solvent{suffix}', None if ref else lvl)
    if not (fp.exists() and fs.exists()):
        return None
    sp, ss = read_ave_time_file(fp), read_ave_time_file(fs)
    if not sp or not ss:
        return None
    n = min(len(sp), len(ss))
    P = np.array([x[2] for x in sp[:n]])
    S = np.array([x[2] for x in ss[:n]])
    return dict(ts=np.array([x[0] for x in sp[:n]], float), bins=sp[0][1], p=P, s=S, t=P + S)


def baseline_mask(z, cfg):
    """Flat far-reservoir bins used for the pore-pressure baseline (drops the
    half-empty extreme-edge bin)."""
    zf = (z - z.min()) / (z.max() - z.min())
    return ((zf >= cfg.baseline_zf - cfg.baseline_zf_half)
            & (zf <= cfg.baseline_zf + cfg.baseline_zf_half) & (zf < 0.995))


def terzaghi_split(tot_stack, bw, sd_bin, ci=0.95):
    """sigma' = sigma^t - p_pore per snapshot; p_pore = mean over the baseline bins.
    Returns (net_stack, pore_val, pore_half, net_half)."""
    tot = np.asarray(tot_stack, float)
    net = np.zeros_like(tot)
    pore = np.zeros(len(tot))
    pore_h = np.zeros(len(tot))
    for i in range(len(tot)):
        v = tot[i][bw]
        v = v[np.isfinite(v)]
        p0 = float(np.nanmean(v)) if len(v) else 0.0
        se = float(stats.sem(v)) if len(v) > 1 else 0.0
        tcr = stats.t.ppf(0.5 + ci / 2, df=max(len(v) - 1, 1))
        pore[i] = p0
        pore_h[i] = tcr * se
        net[i] = tot[i] - p0
    net_half = np.sqrt((_Z95 * np.asarray(sd_bin)[None, :]) ** 2 + pore_h[:, None] ** 2)
    return net, pore, pore_h, net_half


def load_density(path):
    """solvent_density_z: cols chunk, Coord1, Ncount, n_dens, m_dens -> (ts, z, n, m)."""
    snaps = read_ave_chunk_file(path)
    ts = np.array([s[0] for s in snaps], float)
    z = snaps[0][1][:, 1]
    return ts, z, np.array([s[1][:, 3] for s in snaps]), np.array([s[1][:, 4] for s in snaps])


def load_reference(cfg):
    """Everything that belongs to the uncompressed (eps = 0) state:
    box geometry, wall planes, reference stress profiles (zz, and xx/yy when the
    run wrote them), gel bounds, reference network stress per component and the
    reference solvent density.  Returns a dict R."""
    R = {}
    # ---- geometry from a dump header (box fixed; only the piston moves) ----
    src = cfg.traj('traj_ref')
    if not src.exists():
        alts = [cfg.traj('traj_stress', l) for l in cfg.COMP_LEVELS if cfg.traj('traj_stress', l).exists()]
        if not alts:
            raise FileNotFoundError(f'no trajectory found for the box header: {src.name} '
                                    '(run the sync cell -- traj_ref is needed once).')
        src = alts[0]
    box = read_box(src)
    R.update(BOX=box, Z_LO=box['z'][0], Z_HI=box['z'][1], LX=box['lx'], LY=box['ly'], LZ=box['lz'])
    R['AREA'] = box['lx'] * box['ly']
    R['V_BIN'] = R['AREA'] * cfg.binWidth
    wz = wall_z_first_frame(cfg.traj('traj_ref'), (4, 5))
    R['z_support'], R['z_piston'] = wz.get(4, np.nan), wz.get(5, np.nan)
    print(f'box (fixed): Lx={box["lx"]:.2f} Ly={box["ly"]:.2f} Lz={box["lz"]:.2f}  |  '
          f'A={R["AREA"]:.2f}  V_bin={R["V_BIN"]:.2f}')
    print(f'support (type 4) z = {R["z_support"]:.2f}  |  reference piston (type 5) z = {R["z_piston"]:.2f}')

    # ---- reference stress profiles, per component ----------------------
    R['stress'] = {}
    for comp in COMPONENTS:
        S = _component_stacks(cfg, comp, ref=True)
        if S is None:
            if comp == 'zz':
                raise FileNotFoundError('reference sigmazz_{polymer,solvent}_ref files missing -- sync first.')
            print(f'NOTE: no reference sigma{comp} files -> {comp} panels / G skipped')
            continue
        R['stress'][comp] = S
    zz = R['stress']['zz']
    R['z'] = R['Z_LO'] + zz['bins'] * cfg.binWidth - cfg.binWidth / 2.0
    z = R['z']
    R['n_ref'] = len(zz['ts'])
    print(f'reference stress: {R["n_ref"]} snapshots on {len(z)} z-bins '
          f'[{z.min():.1f}, {z.max():.1f}]; components {sorted(R["stress"])}')

    # ---- gel bounds from the reference polymer zz stress ---------------
    pm = np.abs(np.nanmean(zz['p'], axis=0))
    pmax = float(pm.max())
    gel = (pm > cfg.gel_thresh * pmax) if pmax > 0 else np.zeros(len(z), bool)
    R['z_gel_lo'] = float(z[gel].min()) if gel.any() else z[0]
    R['z_gel_hi'] = float(z[gel].max()) if gel.any() else z[-1]
    R['in_gel'] = (z >= R['z_gel_lo']) & (z <= R['z_gel_hi'])
    R['interior'] = (R['in_gel'] & (z >= R['z_gel_lo'] + cfg.wall_margin)
                     & (z <= (R['z_piston'] if np.isfinite(R['z_piston']) else R['z_gel_hi']) - cfg.wall_margin))
    print(f'gel interior (reference): z in [{R["z_gel_lo"]:.1f}, {R["z_gel_hi"]:.1f}]  ({R["in_gel"].sum()} bins)')

    # ---- reference totals, Terzaghi split per component ----------------
    R['bw'] = baseline_mask(z, cfg)
    for comp, S in R['stress'].items():
        S['t_m'], S['t_lo'], S['t_hi'] = mean_ci(S['t'], cfg.ci_level)
        S['p_m'] = np.nanmean(S['p'], axis=0)
        S['s_m'] = np.nanmean(S['s'], axis=0)
        S['sd_bin'] = np.nanstd(S['t'], axis=0)
        net, pore, pore_h, _ = terzaghi_split(S['t'], R['bw'], S['sd_bin'], cfg.ci_level)
        S['net'] = net
        S['net_m'], S['net_lo'], S['net_hi'] = mean_ci(net, cfg.ci_level)
        S['pore'] = float(np.mean(pore))
        S['net_interior'] = float(np.nanmean(S['net_m'][R['interior']])) if R['interior'].any() else np.nan
    print('reference pore baseline p_pore  ' + '  '.join(
        f'{c}: {S["pore"]:.4f}' for c, S in R['stress'].items()))
    print('reference network stress in gel interior  ' + '  '.join(
        f"sigma'_{c}: {S['net_interior']:+.4f}" for c, S in R['stress'].items()))

    # ---- reference solvent density + mass-fraction volume fraction -----
    fd = cfg.path('solvent_density_z_ref')
    if fd.exists():
        rts, dz, _, rm = load_density(fd)
        rho_m, rho_lo, rho_hi = mean_ci(rm, cfg.ci_level)
        rmax = float(np.nanmax(rho_m))
        res = rho_m >= 0.85 * rmax
        R['rho_s0'] = float(np.nanmean(rho_m[res]))
        R.update(dens_ts=rts, dens_z=dz, dens_m=rm)
        R['phi_mf'] = mean_ci(np.array([np.interp(z, dz, row) for row in rm]) / R['rho_s0'], cfg.ci_level)
        print(f'reference density: {len(rts)} snapshots; rho_s,0 (bulk reservoir) = {R["rho_s0"]:.4f}')
    else:
        R['rho_s0'] = np.nan
        R['phi_mf'] = None
        print(f'NOTE: {fd.name} missing -> reference mass-fraction phi skipped')
    R['phi_vor'] = R['phi_cal'] = None       # filled by add_volume_fractions
    return R


# ===========================================================================
#  5. ONE STRAIN LEVEL
# ===========================================================================
def _piston_z_func(cfg, R, lvl):
    f = cfg.path('piston_position', lvl)
    if f.exists():
        pp = np.atleast_2d(np.loadtxt(f, comments='#'))
        pt, pz = pp[:, 0], pp[:, 1]
        return (lambda t: float(np.interp(t, pt, pz))), (pt, pz)
    return (lambda t: R['z_piston']), None


def load_disp(cfg, R, lvl, halt_ts=None):
    """disp_z_polymer + piston_position + gel BB for one level -> dict/None."""
    f = cfg.path('disp_z_polymer', lvl)
    if not f.exists():
        return None
    snaps = read_ave_chunk_file(f)
    if not snaps:
        return None
    d = dict(ts=np.array([s[0] for s in snaps], float), z=snaps[0][1][:, 1],
             Nc=np.array([s[1][:, 2] for s in snaps]), uz=np.array([s[1][:, 3] for s in snaps]))
    fp = cfg.path('piston_position', lvl)
    if fp.exists():
        pp = np.atleast_2d(np.loadtxt(fp, comments='#'))
        d['z_pist_held'] = float(pp[-1, 1])
        d['t_hold'] = float(pp[np.argmax(np.isclose(pp[:, 1], pp[-1, 1])), 0])
    else:
        d['z_pist_held'] = R['z_piston']
        d['t_hold'] = float(d['ts'][0])
    fb = cfg.path('gel_dimensions_bb', lvl)
    if fb.exists():
        bb = np.atleast_2d(np.loadtxt(fb, comments='#'))
        h = bb[:, 0] >= d['t_hold'] + 0.1 * (bb[-1, 0] - d['t_hold'])
        rows = bb[h] if h.any() else bb[-1:]
        d['L_bb'] = float(np.median(rows[:, 3]))
    return d


def _w_modes(zh, kk):
    return (2.0 * zh - 1.0)[:, None] + np.cos(np.pi * zh[:, None] * kk[None, :])


def fit_Dc(cfg, R, disp):
    """Two-sided consolidation fit of D_c to u_z/L on the polymer domain
    (see the D_c notes in the notebooks).  Returns a dict or None."""
    if disp is None or not np.isfinite(R['z_support']):
        return None
    ts, z, Nc = disp['ts'], disp['z'], disp['Nc']
    span = disp['z_pist_held'] - R['z_support']
    L = disp.get('L_bb', span - 2.0)
    gap = 0.5 * (span - L)
    z_perm, z_feed = R['z_support'] + gap, disp['z_pist_held'] - gap
    DL = R['z_piston'] - disp['z_pist_held']
    if not (L > 0 and DL > 0):
        return None
    zeta = (z - z_perm) / L
    uhat = disp['uz'] / L
    idx = np.where((Nc.min(axis=0) > cfg.Ncount_min) & (zeta > 0) & (zeta < 1))[0]
    if len(idx) < 4 + 2 * cfg.DC_TRIM_BINS:
        return None
    if cfg.DC_TRIM_BINS:
        idx = idx[cfg.DC_TRIM_BINS:-cfg.DC_TRIM_BINS]
    zf = zeta[idx]
    t_lj = (ts - disp['t_hold']) * cfg.dt_lj
    early = np.where(t_lj <= cfg.DC_FRAC_EARLY * t_lj[-1])[0]
    if len(early) < 2:
        return None
    if np.nanmax(np.abs(disp['uz'][:, idx])) > 0.5 * DL:
        print('  WARNING: |u_dat| ~ DL -- the disp file does not look hold-referenced.')
    kk = (2.0 * np.arange(1, cfg.DC_N_MODES + 1) - 1.0 if cfg.DC_FREE_AMPS
          else np.arange(1, cfg.DC_KMAX_IC + 1, 2.0))

    def design(Dc, t):
        dec = np.exp(-(np.pi * kk) ** 2 * Dc * t / L ** 2)
        if cfg.DC_FREE_AMPS:
            return _w_modes(zf, kk) * (dec - 1.0)[None, :]
        return (_w_modes(zf, kk) @ (8.0 / (np.pi * kk) ** 2 * dec))[:, None]

    def amps(Dc):
        X = np.vstack([design(Dc, t_lj[i]) for i in early])
        y = np.concatenate([uhat[i][idx] for i in early])
        return np.linalg.lstsq(X, y, rcond=None)[0]

    def resid(Dc):
        A = amps(Dc)
        return float(sum(np.sum((design(Dc, t_lj[i]) @ A - uhat[i][idx]) ** 2) for i in early))

    Dc = float(minimize_scalar(resid, bounds=cfg.DC_BOUNDS, method='bounded').x)
    A = amps(Dc)
    y_all = np.concatenate([uhat[i][idx] for i in early])
    p_all = np.concatenate([design(Dc, t_lj[i]) @ A for i in early])
    ss_t = np.sum((y_all - np.mean(y_all)) ** 2)
    R2 = float(1.0 - np.sum((y_all - p_all) ** 2) / ss_t) if ss_t > 1e-30 else np.nan

    def T(zh, t):
        dec = np.exp(-(np.pi * kk) ** 2 * Dc * t / L ** 2)
        if cfg.DC_FREE_AMPS:
            return _w_modes(zh, kk) @ (A * dec)
        return A[0] * (_w_modes(zh, kk) @ (8.0 / (np.pi * kk) ** 2 * dec))

    u_model = lambda zh, t: -(DL / L) * zh + T(zh, t)
    u_IC = lambda zh: u_model(zh, 0.0)
    beta = np.nan if cfg.DC_FREE_AMPS else float(A[0])
    hold_T = float(t_lj[-1])
    return dict(Dc=Dc, A=A, beta=beta, R2=R2, L=L, DL=DL, gap=gap, z_perm=z_perm, z_feed=z_feed,
                zeta=zeta, idx=idx, zf=zf, uhat=uhat, early=early, t_lj=t_lj, ts=ts,
                T=T, u_model=u_model, u_IC=u_IC, kk=kk, hold_T=hold_T,
                hold_check=hold_adequacy(cfg, L, hold_T, Dc))


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
    """Was the hold long enough?  tau_1 = L^2/(pi^2 D_c) for the fitted and the
    slow reference D_c; residual = mean excess stress over the last plateau_frac."""
    out = {}
    for tag, Dx in (('fit', Dc_fit), ('slow', cfg.DC_SLOW_REF)):
        tau1 = L * L / (np.pi ** 2 * Dx)
        end, avg = hold_residual(hold_T, tau1, cfg.plateau_frac)
        out[tag] = dict(Dc=Dx, tau1=tau1, end=end, avg=avg,
                        need=n_tau_needed(tau1, cfg.DC_TARGET_RESID, cfg.plateau_frac),
                        ok=avg <= cfg.DC_TARGET_RESID)
    return out


def load_level(cfg, R, lvl, verbose=True):
    """Everything for ONE applied-strain level `lvl` (string, e.g. "0.10"):
    production stress stacks (zz + xx/yy when present), Terzaghi split, piston
    history + plateau, M (network & piston), G (from xx and yy), the D_c fit
    and kappa = D_c/M.  Returns a dict L, or None if the core files are missing."""
    say = print if verbose else (lambda *a, **k: None)
    z = R['z']
    L = dict(lvl=lvl, eps=float(lvl), z=z)
    say(f'\n=== level _c{lvl}  (applied strain {float(lvl):.3f}) ===')

    # ---- production stress stacks per component -------------------------
    L['stress'] = {}
    for comp in COMPONENTS:
        S = _component_stacks(cfg, comp, lvl=lvl)
        if S is None:
            if comp == 'zz':
                say(f'  level {lvl}: sigmazz files missing -- skipping level')
                return None
            say(f'  NOTE: no sigma{comp} files for level {lvl} -> {comp} panels / G_{comp[0]} skipped')
            continue
        L['stress'][comp] = S
    zz = L['stress']['zz']
    ts = zz['ts']
    L['ts'] = ts
    t0, t1 = float(ts[0]), float(ts[-1])
    L['halt_ts'] = int(t1 - cfg.plateau_frac * (t1 - t0))     # start of the plateau window
    L['evol_from'] = int(t0)                                  # start of the hold (evolution plots)
    say(f'  production: {len(ts)} snapshots, steps {int(t0)} -> {int(t1)};  '
        f'plateau window: steps >= {L["halt_ts"]} (last {cfg.plateau_frac:.0%})')

    # ---- piston position, membrane bounds ------------------------------
    L['piston_z_at'], L['piston_pos'] = _piston_z_func(cfg, R, lvl)
    L['z_pist'] = L['piston_z_at'](t1)
    spf = np.abs(zz['p'][-1])
    thr = cfg.gel_thresh * float(np.nanmax(spf))
    mem = spf > thr
    L['z_mem_lo'] = float(z[mem].min()) if mem.any() else R['z_gel_lo']
    L['z_mem_hi'] = float(z[mem].max()) if mem.any() else R['z_gel_hi']
    L['in_mem'] = (z >= L['z_mem_lo']) & (z <= L['z_mem_hi'])
    L['interior'] = (L['in_mem'] & (z >= L['z_mem_lo'] + cfg.wall_margin)
                     & (z <= L['z_pist'] - cfg.wall_margin))
    say(f'  membrane (from final polymer stress): z in [{L["z_mem_lo"]:.1f}, {L["z_mem_hi"]:.1f}] '
        f'({L["in_mem"].sum()} bins);  held piston z = {L["z_pist"]:.2f}')

    # ---- Terzaghi split per component ----------------------------------
    bw = R['bw']
    for comp, S in L['stress'].items():
        Rs = R['stress'].get(comp)
        sd_bin = Rs['sd_bin'] if Rs is not None else None
        if sd_bin is None or Rs['t'].shape[0] < 2 or not np.any(sd_bin > 0):
            sd_bin = np.full(len(z), float(np.nanstd(S['t'][-1][bw])))
        S['net'], S['pore'], S['pore_half'], S['net_half'] = terzaghi_split(S['t'], bw, sd_bin, cfg.ci_level)
        S['ref_net'] = (Rs['net_interior'] if (cfg.G_SUBTRACT_REF and Rs is not None) else 0.0)
        S['net_mem'] = np.array([np.nanmean(S['net'][i][L['in_mem']]) for i in range(len(ts))])
    say(f"  pore pressure (zz baseline): {zz['pore'][0]:.4f} -> {zz['pore'][-1]:.4f};  "
        f"network sigma'_zz in membrane: {zz['net_mem'][0]:.4f} -> {zz['net_mem'][-1]:.4f}")

    # ---- strains ---------------------------------------------------------
    fs = cfg.path('strain_zz', lvl)
    if fs.exists():
        s_ts, s_eps = read_strain_file(fs)
        plat = s_ts >= L['halt_ts']
        L['eps_rg'] = float(np.mean(s_eps[plat])) if plat.any() else float(s_eps[-1])
        L['strain_ts'], L['strain_eps'] = s_ts, s_eps
    bb = load2c(cfg.path('gel_dimensions_bb', lvl), 4)
    if bb is not None and bb[0, 3] != 0:
        pl = bb[:, 0] >= L['halt_ts']
        L['eps_bb'] = float(np.mean((bb[0, 3] - bb[pl, 3]) / bb[0, 3])) if pl.any() else float((bb[0, 3] - bb[-1, 3]) / bb[0, 3])
    ps = load2c(cfg.path('strain_piston', lvl), 4)
    if ps is not None:
        L['eps_boundary'] = float(np.median(ps[ps[:, 0] >= L['halt_ts'], 3]))
    say(f"  strain: applied {L['eps']:.4f} (M denominator) | measured Rg {L.get('eps_rg', np.nan):.4f}  "
        f"BB {L.get('eps_bb', np.nan):.4f}  boundary(diag) {L.get('eps_boundary', np.nan):.4f}")

    # ---- solvent density / mass fraction ---------------------------------
    fd = cfg.path('solvent_density_z', lvl)
    if fd.exists() and np.isfinite(R.get('rho_s0', np.nan)):
        d_ts, d_z, _, d_m = load_density(fd)
        L.update(dens_ts=d_ts, dens_z=d_z, dens_m=d_m)
        mf = np.array([np.interp(z, d_z, row) for row in d_m]) / R['rho_s0']
        L['mf_stack'] = mf
        pl = d_ts >= L['halt_ts']
        L['phi_mf'] = mean_ci(mf[pl] if pl.any() else mf[-1:], cfg.ci_level)
    L['phi_vor'] = L['phi_cal'] = None

    # ---- piston force / pressure + plateau -------------------------------
    area = R['AREA']
    B = load2c(cfg.path('box_dimensions', lvl), 3)
    if B is not None:
        area = float(np.mean(B[-5:, 1] * B[-5:, 2]))
    L['area'] = area
    ff, ffa = cfg.path('piston_force', lvl), cfg.path('piston_force_avg', lvl)
    if ff.exists():
        pf = read_print_file(ff, ['step', 'Fz'])
        L['pf_step'], L['pf_F'] = pf['step'].astype(int), pf['Fz']
        L['pf_P'] = pf['Fz'] / area
    if ffa.exists():
        pfa = read_print_file(ffa, ['step', 'Fz'])
        L['pfa_step'], L['pfa_P'] = pfa['step'].astype(int), pfa['Fz'] / area
    if 'pfa_P' in L:
        L['PF'] = plateau_window(L['pfa_step'], L['pfa_P'], cfg.plateau_frac_auto, cfg.ci_level)
        L['PF']['src'] = 'LMP block-avg'
    elif 'pf_P' in L:
        L['PF'] = plateau_window(L['pf_step'], L['pf_P'], cfg.plateau_frac_auto, cfg.ci_level)
        L['PF']['src'] = 'raw'
    if 'PF' in L:
        p = L['PF']
        say(f"  piston plateau <P> = {p['mean']:.4f} [{p['lo']:.4f}, {p['hi']:.4f}]  "
            f"({p['src']}; auto window last {p['frac']:.0%}, n={p['n']}, block={p['block']}, tau~{p['tau']:.1f})"
            + ('  DRIFT WARNING: no drift-free window -- extend the hold' if p.get('warn') else ''))
    else:
        say('  NOTE: no piston_force file -> M_piston skipped')

    # ---- M: network and piston ------------------------------------------
    eps = L['eps']
    mn = zz['net'][-1][L['in_mem']] / eps
    mn = mn[np.isfinite(mn)]
    L['M_net'], L['M_net_lo'], L['M_net_hi'] = mean_ci(mn, cfg.ci_level)
    L['M_net_nbins'] = len(mn)
    if 'PF' in L:
        p = L['PF']
        L['M_pist'], L['M_pist_lo'], L['M_pist_hi'] = p['mean'] / eps, p['lo'] / eps, p['hi'] / eps
        L['P_final'] = p['mean']
    say(f"  M_network = {L['M_net']:.4f} [{L['M_net_lo']:.4f}, {L['M_net_hi']:.4f}]"
        + (f"   M_piston = {L['M_pist']:.4f} [{L['M_pist_lo']:.4f}, {L['M_pist_hi']:.4f}]"
           f"   ratio {L['M_pist'] / L['M_net']:.4f}" if 'M_pist' in L else ''))

    # ---- G from the lateral network stress -------------------------------
    # uniaxial strain, fixed lateral box:  sigma'_zz = M eps,  sigma'_xx = (M - 2G) eps
    #   -> G = (sigma'_zz - sigma'_xx) / (2 eps),  sigma'_zz/sigma'_xx = M/(M - 2G)
    # increments relative to the eps = 0 reference when G_SUBTRACT_REF (the
    # lateral network stress need not vanish in the reference state).
    L['G'] = {}
    dzz = zz['net'] - zz['ref_net']
    for comp in ('xx', 'yy'):
        S = L['stress'].get(comp)
        if S is None:
            continue
        dxx = S['net'] - S['ref_net']
        g_bins = (dzz[-1] - dxx[-1])[L['in_mem']] / (2.0 * eps)
        g_bins = g_bins[np.isfinite(g_bins)]
        Gm, Glo, Ghi = mean_ci(g_bins, cfg.ci_level)
        lam = float(np.nanmean(dxx[-1][L['in_mem']])) / eps
        with np.errstate(invalid='ignore', divide='ignore'):
            ratio = np.array([np.nanmean(dzz[i][L['in_mem']]) / np.nanmean(dxx[i][L['in_mem']])
                              for i in range(len(ts))])
        L['G'][comp] = dict(G=Gm, lo=Glo, hi=Ghi, lam=lam, ratio=ratio, nbins=len(g_bins),
                            G_pist=(0.5 * (L['M_pist'] - lam) if 'M_pist' in L else np.nan),
                            ratio_final=float(ratio[-1]))
        say(f"  G from {comp}: {Gm:.4f} [{Glo:.4f}, {Ghi:.4f}]   lambda_{comp} = {lam:.4f}   "
            f"sigma'_zz/sigma'_{comp} (final) = {ratio[-1]:.3f}"
            + (f"   [G via M_piston = {L['G'][comp]['G_pist']:.4f}]" if 'M_pist' in L else '')
            + (f"   (ref sigma'_{comp} subtracted: {S['ref_net']:+.4f})" if cfg.G_SUBTRACT_REF else ''))

    # ---- D_c and kappa ---------------------------------------------------
    L['Dc'] = fit_Dc(cfg, R, load_disp(cfg, R, lvl))
    if L['Dc'] is None:
        say('  NOTE: no D_c (missing disp_z_polymer / piston_position, no support plane, or too few bins)')
    else:
        F = L['Dc']
        say(f"  D_c = {F['Dc']:.4e} sigma^2/tau  (R^2 = {F['R2']:.3f};  L = {F['L']:.2f}, "
            f"DL/L = {F['DL'] / F['L']:.4f}, hold = {F['hold_T']:.0f} tau)")
        L['kappa'] = {}
        for key, Mk in (('net', 'M_net'), ('pist', 'M_pist')):
            if Mk in L:
                L['kappa'][key] = dict(k=F['Dc'] / L[Mk], lo=F['Dc'] / L[Mk + '_hi'], hi=F['Dc'] / L[Mk + '_lo'])
        say('  kappa = D_c/M:  ' + '   '.join(f"{k}: {v['k']:.4e} [{v['lo']:.4e}, {v['hi']:.4e}]"
                                             for k, v in L['kappa'].items()))
    return L


# ===========================================================================
#  6. VOLUME FRACTIONS  (mass fraction / Voronoi / lambda-calibrated)
# ===========================================================================
_VF_CACHE = {}


def _volume_fractions(cfg, R, traj, ts_want, label=''):
    """One streamed pass over `traj` -> dict(ts, phi_vor, phi_vor_mobile, phi_cal)
    on the R['z'] grid, cached on (file, frames, knobs)."""
    ts_want = [int(t) for t in np.atleast_1d(ts_want)]
    key = (str(traj), tuple(ts_want), cfg.VOR_MOBILE_ONLY, cfg.VOR_NORM, cfg.binWidth, cfg.P_CAL)
    if key in _VF_CACHE:
        print(f'    {label}cached: {len(_VF_CACHE[key]["ts"])} frame(s), no re-read')
        return _VF_CACHE[key]
    print(f'    {label}streaming {Path(traj).name} for {len(ts_want)} frame(s): Voronoi (~20 s/frame) ...')
    fr = volfrac.stream_traj_frames(traj, ts_want)
    missing = [t for t in ts_want if t not in fr]
    if missing:
        print(f'    NOTE: no traj frame at {missing} -- dropped')
    ts_ok = [t for t in ts_want if t in fr]
    out = {'ts': np.array(ts_ok, float), 'phi_vor': None, 'phi_vor_mobile': None, 'phi_cal': None}
    if not ts_ok:
        return out
    pv, pm = [], []
    for t in ts_ok:
        box, typ, xyz = fr[t]
        zc, ph_b, ph_m = volfrac.phi_voronoi_frame(box, typ, xyz, cfg.binWidth,
                                                   mobile_only=cfg.VOR_MOBILE_ONLY, norm='both')
        ph = ph_b if cfg.VOR_NORM == 'bin' else ph_m
        pv.append(np.interp(R['z'], zc, ph, left=np.nan, right=np.nan))
        pm.append(np.interp(R['z'], zc, ph_m, left=np.nan, right=np.nan))
        print(f'      ts {t}: phi^vor max = {np.nanmax(pv[-1]):.3f}')
    del fr
    out['phi_vor'] = np.array(pv)
    out['phi_vor_mobile'] = np.array(pm)
    if R.get('CALIB') is not None:
        out['phi_cal'] = volfrac.phi_calibrated(out['phi_vor_mobile'], cfg.P_CAL, R['CALIB'])
    _VF_CACHE[key] = out
    return out


def add_volume_fractions(cfg, R, levels=()):
    """Voronoi + lambda-calibrated solvent volume fractions for the reference
    state (R) and each level in `levels`.  Reference: REF_VOR_FRAMES frames
    evenly over the reference window (mean + CI).  Level: VOR_MAX_FRAMES frames
    inside the plateau window (plateau mean + CI).  Mass-fraction profiles were
    already set by load_reference / load_level."""
    if 'CALIB' not in R:
        try:
            R['CALIB'] = volfrac.load_calibration()
            print(f'calibration loaded ({volfrac.CALIBRATION_JSON.name}): lambda coeffs {R["CALIB"]["coeffs"]}')
        except FileNotFoundError:
            R['CALIB'] = None
            print('NOTE: no calibration artifact -> phi_cal skipped '
                  '(generate scripts/calibration/calibration_lambda.json with calibration_analysis.ipynb)')
    if not cfg.VOR_ENABLE:
        print('VOR_ENABLE=False -> Voronoi / calibrated phi skipped')
        return
    # ---- reference ------------------------------------------------------
    tr = cfg.traj('traj_ref')
    if R.get('phi_vor') is None:
        if tr.exists() and 'dens_ts' in R:
            want, _ = subsample(R['dens_ts'], R['dens_ts'][:, None], cfg.REF_VOR_FRAMES)
            vf = _volume_fractions(cfg, R, tr, want, 'reference: ')
            if vf['phi_vor'] is not None:
                R['phi_vor'] = mean_ci(vf['phi_vor'], cfg.ci_level)
                R['phi_vor_frames'] = vf['ts']
                if vf['phi_cal'] is not None:
                    R['phi_cal'] = mean_ci(vf['phi_cal'], cfg.ci_level)
        else:
            print(f'NOTE: {tr.name} or reference density missing -> reference Voronoi skipped')
    # ---- levels ---------------------------------------------------------
    for L in levels:
        if L is None or L.get('phi_vor') is not None:
            continue
        tp = cfg.traj('traj_stress', L['lvl'])
        if not tp.exists():
            print(f'NOTE: {tp.name} missing -> Voronoi skipped for level {L["lvl"]}')
            continue
        plat = L['ts'][L['ts'] >= L['halt_ts']]
        if not len(plat):
            plat = L['ts'][-1:]
        want, _ = subsample(plat, plat[:, None], cfg.VOR_MAX_FRAMES)
        vf = _volume_fractions(cfg, R, tp, want, f'level {L["lvl"]}: ')
        if vf['phi_vor'] is not None:
            L['phi_vor'] = mean_ci(vf['phi_vor'], cfg.ci_level)
            L['phi_vor_frames'] = vf['ts']
            if vf['phi_cal'] is not None:
                L['phi_cal'] = mean_ci(vf['phi_cal'], cfg.ci_level)
    # ---- printed in-gel means ------------------------------------------
    def _line(tag, D, mask):
        parts = []
        for k, lab in (('phi_mf', 'mass-frac'), ('phi_vor', 'Voronoi'), ('phi_cal', 'calibrated')):
            if D.get(k) is not None:
                parts.append(f'{lab} {fmt_mu(D[k][0][mask])}')
        print(f'  {tag:<22s} ' + '   '.join(parts))
    print('in-gel solvent volume fractions (interior, wall_margin trimmed):')
    _line('reference (eps = 0)', R, R['interior'])
    for L in levels:
        if L is not None:
            _line(f'compressed eps={L["eps"]:.2f}', L, L['interior'])


# ===========================================================================
#  7. EXPANSE SYNC
# ===========================================================================
_PROD_DAT = ('sigmazz_polymer', 'sigmazz_solvent', 'sigmaxx_polymer', 'sigmaxx_solvent',
             'sigmayy_polymer', 'sigmayy_solvent', 'solvent_density_z', 'disp_z_polymer',
             'strain_zz', 'strain_piston', 'piston_position', 'piston_force', 'piston_force_avg',
             'box_dimensions', 'gel_dimensions_bb', 'gel_dimensions_rg', 'polymer_com', 'gel_edges')
_REF_DAT = ('sigmazz_polymer_ref', 'sigmazz_solvent_ref', 'sigmaxx_polymer_ref', 'sigmaxx_solvent_ref',
            'sigmayy_polymer_ref', 'sigmayy_solvent_ref', 'solvent_density_z_ref')
_REQUIRED_PROD = ('sigmazz_polymer', 'sigmazz_solvent', 'solvent_density_z', 'strain_zz',
                  'piston_force', 'box_dimensions', 'gel_dimensions_bb', 'disp_z_polymer')


def sync_files(cfg, levels=None):
    """(data_files, traj_files, required) the notebooks read for `levels`
    (default: every level in cfg.COMP_LEVELS)."""
    levels = cfg.COMP_LEVELS if levels is None else [str(l) for l in levels]
    data = [cfg.path(n) for n in _REF_DAT]
    traj = [cfg.traj('traj_ref')]
    req = [cfg.path('sigmazz_polymer_ref'), cfg.path('sigmazz_solvent_ref'), cfg.path('solvent_density_z_ref')]
    for l in levels:
        data += [cfg.path(n, l) for n in _PROD_DAT]
        traj += [cfg.traj('traj_stress', l)]
        req += [cfg.path(n, l) for n in _REQUIRED_PROD]
    return data, traj, req


def sync_from_expanse(cfg, levels=None, force=False):
    """Pull every file the notebooks read from Expanse in ONE login (password +
    TOTP prompts).  Logs in when ANY data (.dat) file is missing locally, or when
    force=True (which also refreshes the large trajectories).  Files that a
    previous sync looked for and did NOT find on the cluster are remembered in
    DATA_DIR/.sync_absent.json so an old run that never wrote them does not
    prompt for a login every time; force=True clears that memory."""
    import paramiko
    import getpass
    import json
    data_files, traj_files, required = sync_files(cfg, levels)
    absent_f = cfg.DATA_DIR / '.sync_absent.json'
    absent = set()
    if absent_f.exists() and not force:
        try:
            absent = set(json.loads(absent_f.read_text()))
        except Exception:
            absent = set()
    missing_dat = [f for f in data_files if not f.exists()]
    missing_new = [f for f in missing_dat if f.name not in absent]
    missing_all = missing_dat + [f for f in traj_files if not f.exists()]
    missing_req = [f for f in required if not f.exists()]
    if not force and not missing_new:
        known = len(missing_dat) - len(missing_new)
        print(f'All data files present locally' + (f' except {known} known absent on the cluster' if known else '')
              + f' ({len(missing_all)} target file(s) missing in total) -- skipping Expanse login.  '
              f'FORCE_SYNC=True re-checks everything.')
        if missing_req:
            print('  WARNING: required files missing: ' + ', '.join(f.name for f in missing_req))
        return
    why = 'force=True' if (force and not missing_new) else f'{len(missing_new)} data file(s) missing'
    print(f'Syncing from Expanse ({why}); {len(missing_all)} of {len(data_files) + len(traj_files)} target files missing locally.')
    if missing_new:
        print('  missing: ' + ', '.join(f.name for f in missing_new[:6]) + (' ...' if len(missing_new) > 6 else ''))
    bn = lambda paths: ' '.join(f'"{p.name}"' for p in paths)
    stage = f'{cfg.RUNS_ROOT}/triaxial_stage'
    script = r"""
set -u
RUNS="__RUNS__"; TRAJ="__TRAJ__"; STAGE="__STAGE__"
rm -rf "$STAGE"; mkdir -p "$STAGE/data" "$STAGE/traj"
declare -A NEWEST_D
while IFS= read -r line; do p=${line#* }; b=${p##*/}
  [ -z "${NEWEST_D[$b]:-}" ] && NEWEST_D[$b]="$p"
done < <(find "$RUNS" -name '*.dat' -not -path '*/triaxial_stage/*' -printf '%T@ %p\n' 2>/dev/null | sort -rn)
for B in __DATA_BN__; do S="${NEWEST_D[$B]:-}"; [ -n "$S" ] && cp -p "$S" "$STAGE/data/" 2>/dev/null || true; done
declare -A NEWEST_T
while IFS= read -r line; do p=${line#* }; b=${p##*/}
  [ -z "${NEWEST_T[$b]:-}" ] && NEWEST_T[$b]="$p"
done < <(find "$TRAJ" "$RUNS" -name '*.lammpstrj' -not -path '*/triaxial_stage/*' -printf '%T@ %p\n' 2>/dev/null | sort -rn)
for B in __TRAJ_BN__; do S="${NEWEST_T[$B]:-}"; [ -n "$S" ] && cp -p "$S" "$STAGE/traj/" 2>/dev/null || true; done
echo "  staged: $(ls "$STAGE/data" 2>/dev/null | wc -l) data, $(ls "$STAGE/traj" 2>/dev/null | wc -l) traj"
"""
    script = (script.replace('__RUNS__', cfg.RUNS_ROOT).replace('__TRAJ__', cfg.TRAJ_ROOT)
              .replace('__STAGE__', stage).replace('__DATA_BN__', bn(data_files))
              .replace('__TRAJ_BN__', bn(traj_files)))
    password = getpass.getpass(f'Expanse password for {cfg.EXPANSE_USER}: ')
    totp = getpass.getpass('TOTP / verification code: ')

    def auth_handler(title, instructions, prompt_list):
        return [password if 'password' in p.strip().lower() else totp for p, _ in prompt_list]

    import socket
    print(f'Connecting to {cfg.EXPANSE_HOST} ...')
    if 'expanse' in socket.gethostname().lower() or 'sdsc' in socket.gethostname().lower():
        print('  NOTE: this notebook seems to be running ON Expanse; the login node may not accept an SSH\n'
              '  connection from here.  Set SYNC=False and point cfg.base_dir at the run directory instead.')
    try:
        sock = socket.create_connection((cfg.EXPANSE_HOST, 22), timeout=30)
    except (socket.timeout, OSError) as e:
        raise ConnectionError(f'cannot reach {cfg.EXPANSE_HOST}:22 ({type(e).__name__}: {e}).  '
                              'Check network / VPN / that the host name resolves, then re-run this cell.') from e
    transport = paramiko.Transport(sock)
    transport.banner_timeout = 60
    transport.set_keepalive(30)
    try:
        transport.connect()
        transport.auth_interactive(cfg.EXPANSE_USER, auth_handler)
    except paramiko.AuthenticationException as e:
        transport.close()
        raise ConnectionError('Expanse rejected the login: wrong password, or the TOTP code expired '
                              '(codes last 30 s -- have the code ready before running the cell).') from e
    except (paramiko.SSHException, OSError) as e:
        transport.close()
        raise ConnectionError(f'SSH handshake with {cfg.EXPANSE_HOST} failed ({type(e).__name__}: {e}).  '
                              'Repeated failed logins get the connection dropped for a while; wait a minute and retry.') from e
    if not transport.is_authenticated():
        transport.close()
        raise ConnectionError('Expanse login did not complete (transport not authenticated).')
    ssh = paramiko.SSHClient()
    ssh._transport = transport
    print('Step 1 -- staging files on Expanse (find over the runs tree; can take a minute)...')
    try:
        _, stdout, stderr = ssh.exec_command('bash -s', get_pty=False, timeout=600)
        stdout.channel.sendall(script.encode())
        stdout.channel.shutdown_write()
        out, err = stdout.read().decode(), stderr.read().decode()
    except (socket.timeout, OSError, paramiko.SSHException) as e:
        transport.close()
        raise ConnectionError(f'the staging command on Expanse failed ({type(e).__name__}: {e}).') from e
    print(out)
    if err.strip():
        print('  (stderr) ' + err.strip().replace('\n', '\n  (stderr) '))
    print('Step 2 -- downloading via SFTP (skips files already present)...')
    sftp = ssh.open_sftp()

    def pull(remote_dir, local_dir):
        local_dir = Path(local_dir)
        local_dir.mkdir(parents=True, exist_ok=True)
        try:
            entries = sftp.listdir_attr(remote_dir)
        except FileNotFoundError:
            return
        for e in entries:
            rp, lp = f'{remote_dir}/{e.filename}', local_dir / e.filename
            if stat.S_ISDIR(e.st_mode):
                pull(rp, lp)
                continue
            if lp.exists() and lp.stat().st_mtime >= e.st_mtime:
                continue
            sftp.get(rp, str(lp))
    pull(f'{stage}/data', cfg.DATA_DIR)
    pull(f'{stage}/traj', cfg.TRAJ_DIR)
    sftp.close()
    ssh.close()
    still = [f.name for f in data_files if not f.exists()]
    absent_f.write_text(json.dumps(sorted(set(still)), indent=1))
    if still:
        print(f'Sync complete; {len(still)} data file(s) not found on the cluster (remembered in {absent_f.name}): '
              + ', '.join(still[:6]) + (' ...' if len(still) > 6 else ''))
    else:
        print('Sync complete; every data file present.')


# ===========================================================================
#  8. PLOTTING PRIMITIVES
# ===========================================================================
def zn(R, z):
    """Fractional box height z/Lz."""
    return (np.asarray(z, float) - R['Z_LO']) / R['LZ']


def shade_gel(ax, R, L=None):
    lo, hi = (L['z_mem_lo'], L['z_mem_hi']) if L is not None else (R['z_gel_lo'], R['z_gel_hi'])
    ax.axvspan(zn(R, lo), zn(R, hi), **GEL_SHADE)


def mark_walls(ax, R, L=None, ts=None):
    """Support (solid) + piston (dash-dot): reference position when L is None,
    else per evolution timestep (final dark, earlier faded)."""
    if np.isfinite(R['z_support']):
        ax.axvline(zn(R, R['z_support']), color='k', ls='-', lw=1.5, alpha=0.85, zorder=4)
    if L is None or ts is None:
        pistons = [(R['z_piston'] if L is None else L['z_pist'], True)]
    else:
        ts = np.asarray(ts)
        pistons = [(L['piston_z_at'](t), i == len(ts) - 1) for i, t in enumerate(ts)]
    for pz, last in pistons:
        if np.isfinite(pz):
            ax.axvline(zn(R, pz), color=('0.15' if last else '0.7'), ls='-.',
                       lw=(1.6 if last else 1.0), alpha=(0.9 if last else 0.35), zorder=(4 if last else 3))


def finish_axes(ax, ylabel, title):
    ax.set_xlabel(r'$z/L_z$')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.grid(alpha=0.3)


def annotate_box(ax, text, loc='lower left', fontsize=14, color='k'):
    x, ha = (0.02, 'left') if 'left' in loc else (0.98, 'right')
    y, va = (0.03, 'bottom') if 'lower' in loc else (0.97, 'top')
    ax.text(x, y, text, transform=ax.transAxes, va=va, ha=ha, fontsize=fontsize, color=color,
            bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.85))


def flat_inside(curve, z, mask, tol):
    """Flat = small linear TREND across the masked bins, relative to the mean."""
    v = np.asarray(curve)[mask]
    zz = np.asarray(z)[mask]
    m = np.isfinite(v)
    v, zz = v[m], zz[m]
    if len(v) < 3:
        return False
    slope = np.polyfit(zz, v, 1)[0]
    rise = abs(slope) * (zz.max() - zz.min())
    return (rise / max(abs(np.mean(v)), 1e-9)) < tol


def robust_ylim(ax, curves, zmask=None, pad=0.12, qlo=2, qhi=98, include_zero=True):
    vals = []
    for c in curves:
        c = np.asarray(c, float)
        if zmask is not None:
            c = c[zmask]
        c = c[np.isfinite(c)]
        if c.size:
            vals.append(c)
    if not vals:
        return
    v = np.concatenate(vals)
    lo, hi = np.percentile(v, [qlo, qhi])
    if include_zero:
        lo, hi = min(lo, 0.0), max(hi, 0.0)
    if hi <= lo:
        hi = lo + 1.0
    d = (hi - lo) * pad
    ax.set_ylim(lo - d, hi + d)


def post_halt(cfg, L, ts, stack):
    """Evolution-plot subsample over the whole hold (transient -> plateau)."""
    ts = np.asarray(ts)
    stack = np.asarray(stack)
    m = ts >= L['evol_from']
    if m.sum() < 2:
        m = np.zeros(len(ts), bool)
        m[-min(len(ts), cfg.n_curves):] = True
    return subsample(ts[m], stack[m], cfg.n_curves)


def plot_reference(ax, R, z, m, lo, hi, color, ylabel, title, annotate=True):
    zx = zn(R, z)
    ax.fill_between(zx, lo, hi, color=color, alpha=0.25, lw=0, zorder=2)
    ax.plot(zx, m, '-', color=color, lw=2.5, zorder=3)
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    shade_gel(ax, R)
    mark_walls(ax, R)
    finish_axes(ax, ylabel, title)
    if annotate:
        ig = R['in_gel']
        annotate_box(ax, f'mean in gel = {fmt_mu(np.asarray(m)[ig])}\nmean out gel = {fmt_mu(np.asarray(m)[~ig])}',
                     fontsize=15)


def plot_evolution(ax, cfg, R, L, z, ts, stack, ylabel, title, ref=None, band=None, colorbar=True,
                   annotate=True):
    """Time-coloured profiles (cividis); final curve bold black; membrane shaded.
    ref=(mean, lo, hi) draws the eps = 0 reference profile (dashed, band) as the
    starting point of the evolution.  band = (n, nz) half-widths for 95 % bands."""
    ts = np.asarray(ts)
    stack = np.asarray(stack)
    norm = Normalize(vmin=ts.min(), vmax=ts.max())
    cmap = plt.get_cmap(EVO_CMAP)
    zx = zn(R, z)
    if ref is not None:
        rm, rlo, rhi = ref
        ax.fill_between(zx, rlo, rhi, color=WONG['skyblue'], alpha=0.25, lw=0, zorder=1)
        ax.plot(zx, rm, '--', color=WONG['blue'], lw=2.2, alpha=0.9, zorder=2,
                label=r'reference ($\varepsilon=0$)')
    for i in range(len(ts)):
        last = (i == len(ts) - 1)
        c = 'k' if last else cmap(norm(ts[i]))
        if band is not None:
            ax.fill_between(zx, stack[i] - band[i], stack[i] + band[i], color=c,
                            alpha=(0.20 if last else 0.06), lw=0, zorder=(4 if last else 2))
        ax.plot(zx, stack[i], '-', color=c, lw=(3.5 if last else 1.6),
                alpha=(1.0 if last else 0.75), zorder=(5 if last else 3),
                label=('final (plateau)' if last else None))
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    shade_gel(ax, R, L)
    mark_walls(ax, R, L, ts)
    finish_axes(ax, ylabel, title)
    if colorbar:
        sm = plt.cm.ScalarMappable(cmap=EVO_CMAP, norm=norm)
        sm.set_array([])
        cb = ax.figure.colorbar(sm, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label('timestep')
    if annotate:
        interior = L['interior']
        if flat_inside(stack[-1], z, interior, cfg.flat_tol):
            annotate_box(ax, 'final (plateau)\nmean in gel = ' + fmt_mu(stack[-1][interior]))
        else:
            ax.text(0.02, 0.03, 'final not flat inside\n(means omitted)', transform=ax.transAxes,
                    va='bottom', ha='left', fontsize=13, color='0.35')
    if ref is not None:
        ax.legend(fontsize=12, loc='upper right')


def overlay_levels(ax, R, levels, get_z, get_ts, get_stack, cfg, ref=None, autoscale_mask=None,
                   ylabel='', title=''):
    """Sweep overlay: colour = level; within a level the evolution ramps faint ->
    bold with the final curve opaque.  ref=(mean, lo, hi) is drawn dashed black."""
    finals = []
    zx_ref = zn(R, R['z'])
    if ref is not None:
        rm, rlo, rhi = ref
        ax.fill_between(zx_ref, rlo, rhi, color='0.5', alpha=0.2, lw=0, zorder=1)
        ax.plot(zx_ref, rm, '--', color='k', lw=2.0, alpha=0.9, zorder=2)
        finals.append(np.asarray(rm))
    for i, L in enumerate(levels):
        st = get_stack(L)
        if st is None:
            continue
        ts, stack = post_halt(cfg, L, get_ts(L), st)
        zx = zn(R, get_z(L))
        col = level_color(i)
        nc = len(ts)
        for j in range(nc):
            last = (j == nc - 1)
            ramp = 0.20 + 0.80 * (j / max(nc - 1, 1))
            ax.plot(zx, stack[j], '-', color=col, lw=(3.0 if last else 1.1),
                    alpha=(1.0 if last else 0.55 * ramp), zorder=(5 if last else 3))
            if last:
                finals.append(stack[j])
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.4)
    if np.isfinite(R['z_support']):
        ax.axvline(zn(R, R['z_support']), color='k', ls='-', lw=1.5, alpha=0.85, zorder=4)
    for i, L in enumerate(levels):
        ax.axvline(zn(R, L['z_pist']), color=level_color(i), ls='-.', lw=1.2, alpha=0.6, zorder=4)
    finish_axes(ax, ylabel, title)
    if finals:
        robust_ylim(ax, finals, zmask=autoscale_mask, pad=0.15)


def level_handles(levels, ref=False):
    h = [Line2D([0], [0], color=level_color(i), lw=3, label=fr'$\varepsilon={L["eps"]:.2f}$')
         for i, L in enumerate(levels)]
    if ref:
        h.insert(0, Line2D([0], [0], color='k', ls='--', lw=2, label=r'reference ($\varepsilon=0$)'))
    return h


def _save(fig, cfg, stem, lvl=None):
    out = cfg.plot(stem, lvl)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print('saved', out)
    plt.show()
    return fig


# ===========================================================================
#  9. FIGURES -- SINGLE LEVEL
# ===========================================================================
def fig_strain(cfg, R, levels, stem='strain_diagnostic'):
    """Strain vs step: solid eps_Rg, dashed eps_BB, faint eps_piston (diagnostic);
    shaded plateau window.  Works for one level or a sweep (colour = level)."""
    fig, ax = plt.subplots(figsize=(10, 6.5), constrained_layout=True)
    any_ = False
    for i, L in enumerate(levels):
        lvl = L['lvl']
        col = level_color(i) if len(levels) > 1 else WONG['blue']
        S = load2c(cfg.path('strain_zz', lvl), 3)
        if S is None:
            continue
        any_ = True
        st, L0, Lrg = S[:, 0], S[:, 1], S[:, 2]
        eps_rg = (L0 - Lrg) / L0
        ax.plot(st, eps_rg, '-', color=col, lw=2.2, label=fr'$\varepsilon_{{Rg}}$ (target {lvl})')
        BB = load2c(cfg.path('gel_dimensions_bb', lvl), 4)
        if BB is not None and BB[0, 3] != 0:
            ax.plot(BB[:, 0], (BB[0, 3] - BB[:, 3]) / BB[0, 3], '--', color=col, lw=1.5, alpha=0.7,
                    label=fr'$\varepsilon_{{BB}}$ (target {lvl})')
        PS = load2c(cfg.path('strain_piston', lvl), 4)
        if PS is not None:
            ax.plot(PS[:, 0], PS[:, 3], '-', color=col, lw=3.2, alpha=0.30,
                    label=r'$\varepsilon_\mathrm{piston}$ (diag)')
        ax.axvspan(L['halt_ts'], st[-1], color=col, alpha=0.06)
        ax.axhline(L['eps'], color=col, ls=':', lw=1.0, alpha=0.6)
        ax.annotate(f'plateau: Rg {L.get("eps_rg", np.nan):.3f}  BB {L.get("eps_bb", np.nan):.3f}',
                    (st[-1], eps_rg[-1]), textcoords='offset points', xytext=(-6, 9), ha='right',
                    va='bottom', fontsize=10, color=col,
                    bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.8))
    ax.set_xlabel('step')
    ax.set_ylabel(r'compression strain  $\varepsilon=(L_0-L)/L_0$')
    ax.set_title('Strain diagnostic: solid $\\varepsilon_{Rg}$, dashed $\\varepsilon_{BB}$,\n'
                 'dotted = applied target, shaded = plateau window', fontsize=15)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=10, loc='best')
    if not any_:
        plt.close(fig)
        print('strain diagnostic skipped (no strain_zz files)')
        return None
    # one level -> tag the file with it (keeps clear of the long-form notebook's untagged file)
    return _save(fig, cfg, stem, levels[0]['lvl'] if len(levels) == 1 else None)


def _phi_panel(ax, cfg, R, D, L=None, title='', bands=True):
    zx = zn(R, R['z'])
    drawn = []
    for key, lab, col in (('phi_mf', r'$\phi_s^{\rm mf}=\rho_s/\rho_{s,0}$', WONG['blue']),
                          ('phi_vor', r'$\phi_s^{\rm vor}$ (Voronoi)', WONG['vermillion']),
                          ('phi_cal', r'$\phi_s^{\rm cal}=\lambda\,\phi_s^{\rm vor}$', WONG['green'])):
        if D.get(key) is None:
            continue
        m, lo, hi = D[key]
        if bands:
            ax.fill_between(zx, lo, hi, color=col, alpha=0.22, lw=0, zorder=2)
        ax.plot(zx, m, '-', lw=2.4, color=col, label=lab, zorder=3)
        drawn.append((lab, m))
    ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6)
    ax.axhline(cfg.PHI_FLOOR, color='r', ls=':', lw=1.0, alpha=0.5)
    shade_gel(ax, R, L)
    mark_walls(ax, R, L)
    finish_axes(ax, r'$\phi_s$', title)
    ax.set_ylim(0, 1.15)
    ax.legend(fontsize=12, loc='lower left')
    mask = L['interior'] if L is not None else R['interior']
    txt = '\n'.join(f'{"mf" if "mf" in lab else ("vor" if "vor" in lab else "cal")}: {fmt_mu(m[mask])}'
                    for lab, m in drawn)
    if txt:
        annotate_box(ax, 'in-gel means\n' + txt, loc='lower right', fontsize=12)


def fig_volfrac(cfg, R, L):
    """Solvent volume fraction (mass fraction, Voronoi, lambda-calibrated):
    (a) reference eps = 0 (mean + 95 % CI across frames), (b) compressed plateau."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 6.5), constrained_layout=True)
    fig.suptitle(f'Solvent volume fraction: reference vs compressed  ($P_{{\\rm CAL}}={cfg.P_CAL}$)  '
                 f'|  {cfg.sim_name}', fontsize=13, fontweight='bold')
    _phi_panel(axes[0], cfg, R, R, None, r'(a) reference, $\varepsilon=0$')
    _phi_panel(axes[1], cfg, R, L, L, fr'(b) compressed, $\varepsilon={L["eps"]:.2f}$ (plateau mean)')
    return _save(fig, cfg, 'volfrac_profiles', L['lvl'])


def _stress_evo_panels(cfg, R, L, kind, stem, suptitle):
    """1 x 3 evolution panels (zz, xx, yy) for kind='t' (total) or 'net' (network)."""
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    fig.suptitle(suptitle, fontsize=13, fontweight='bold')
    for ax, comp in zip(axes, COMPONENTS):
        S = L['stress'].get(comp)
        Rs = R['stress'].get(comp)
        lab = (r'$\sigma^{t}_{%s}$' % comp) if kind == 't' else (r"$\sigma'_{%s}$" % comp)
        title = f'({"abc"[COMPONENTS.index(comp)]}) ' + (
            r'total $\sigma^{t}_{%s}$' % comp if kind == 't' else r"network $\sigma'_{%s}=\sigma^t_{%s}-p_{\rm pore}$" % (comp, comp))
        if S is None:
            ax.text(0.5, 0.5, f'sigma{comp} files\nnot found', ha='center', va='center', transform=ax.transAxes)
            finish_axes(ax, lab, title)
            continue
        ts, ev = post_halt(cfg, L, L['ts'], S[kind])
        band = None
        if kind == 'net':
            _, band = post_halt(cfg, L, L['ts'], S['net_half'])
        ref = None
        if Rs is not None:
            ref = (Rs['t_m'], Rs['t_lo'], Rs['t_hi']) if kind == 't' else (Rs['net_m'], Rs['net_lo'], Rs['net_hi'])
        plot_evolution(ax, cfg, R, L, R['z'], ts, ev, lab + '$(z,t)$', title, ref=ref, band=band)
        if kind == 'net':
            mask = (R['z'] >= L['z_mem_lo'] + cfg.wall_margin)
            robust_ylim(ax, list(ev), zmask=mask, pad=0.15)
            annotate_box(ax, 'final (plateau)\nmean in membrane = ' + fmt_mu(ev[-1][L['in_mem']]))
    return _save(fig, cfg, stem, L['lvl'])


def fig_total_stress(cfg, R, L):
    return _stress_evo_panels(cfg, R, L, 't', 'total_stress_evolution',
                              f'Total stress evolution, reference -> compressed (steps >= {L["evol_from"]}; '
                              f'plateau >= {L["halt_ts"]})  |  {cfg.sim_name}')


def fig_network_stress(cfg, R, L):
    return _stress_evo_panels(cfg, R, L, 'net', 'network_stress_evolution',
                              f"Network stress evolution (Terzaghi), reference -> compressed  |  "
                              f"p_pore(final) = {fmt_val_unc(L['stress']['zz']['pore'][-1], L['stress']['zz']['pore_half'][-1])}"
                              f"  |  {cfg.sim_name}")


def fig_partial_stress(cfg, R, L):
    """One axis: solvent partial (blues), polymer partial (oranges) and total
    (greys) sigma_zz evolutions, reference dashed, final curves bold."""
    zz, Rz = L['stress']['zz'], R['stress']['zz']
    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    zx = zn(R, R['z'])
    fam = (('s', r'solvent $\sigma_{s,zz}$', 'Blues'), ('p', r'polymer $\sigma_{p,zz}$', 'Oranges'),
           ('t', r'total $\sigma^{t}_{zz}$', 'Greys'))
    handles = []
    for key, lab, cmap_name in fam:
        cmap = plt.get_cmap(cmap_name)
        ts, ev = post_halt(cfg, L, L['ts'], zz[key])
        ref = Rz[key + '_m']
        base = cmap(0.85)
        ax.plot(zx, ref, '--', color=base, lw=1.8, alpha=0.9, zorder=2)
        n = len(ts)
        for i in range(n):
            last = (i == n - 1)
            ax.plot(zx, ev[i], '-', color=(base if last else cmap(0.3 + 0.5 * i / max(n - 1, 1))),
                    lw=(3.2 if last else 1.2), alpha=(1.0 if last else 0.6), zorder=(5 if last else 3))
        handles.append(Line2D([0], [0], color=base, lw=3, label=lab))
    handles.append(Line2D([0], [0], color='0.4', ls='--', lw=2, label=r'reference ($\varepsilon=0$)'))
    handles.append(Line2D([0], [0], color='0.4', lw=1.2, alpha=0.6, label='hold (faint = early)'))
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    shade_gel(ax, R, L)
    mark_walls(ax, R, L)
    finish_axes(ax, r'$\sigma_{zz}(z,t)$  (LJ)',
                f'Partial and total $\\sigma_{{zz}}$: reference -> compressed ($\\varepsilon={L["eps"]:.2f}$)')
    ax.legend(handles=handles, fontsize=12, loc='best')
    return _save(fig, cfg, 'partial_stress_evolution', L['lvl'])


def fig_piston(cfg, R, L):
    """Piston pressure P = F_z/A vs step, linear + log, with the auto plateau window."""
    if 'pf_P' not in L:
        print('piston figure skipped (no piston_force file)')
        return None
    steps, P = L['pf_step'], L['pf_P']
    P_roll = rolling_mean(P, cfg.roll_win)
    peak = int(steps[int(P.argmax())])
    fig, (axL, axG) = plt.subplots(1, 2, figsize=(17, 6), constrained_layout=True)
    fig.suptitle(f'Piston pressure history ($A = {L["area"]:.1f}\\,\\sigma^2$, $F_z = P\\,A$)  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    for ax in (axL, axG):
        ax.axvline(peak, color=WONG['blue'], ls='--', lw=1.6, alpha=0.7, label=f'peak $P$ (t={peak})')
        if 'PF' in L:
            ax.axvspan(L['PF']['step0'], float(steps[-1]), color=WONG['green'], alpha=0.10,
                       label='plateau window (auto)')
        ax.set_xlabel('step')
        ax.grid(alpha=0.3)
    axL.plot(steps, P, '-', color=WONG['vermillion'], lw=1.0, alpha=0.30, label=r'$P=F_z/A$ (raw)')
    axL.plot(steps, P_roll, '-', color=WONG['vermillion'], lw=2.6, alpha=0.95, label=f'rolling mean ({cfg.roll_win})')
    if 'pfa_P' in L:
        axL.plot(L['pfa_step'], L['pfa_P'], '-', color='k', lw=1.8, alpha=0.8, label='LMP block-avg')
    axL.axhline(0, color='k', ls='--', lw=0.8, alpha=0.4)
    axL.set_ylabel(r'$P = F_z/A$  (LJ / $\sigma^2$)')
    axL.set_title('(a) piston pressure')
    axL.legend(fontsize=12, loc='upper center')     # upper right is the plateau annotation box
    pos = P > 0
    axG.plot(steps[pos], np.log(P[pos]), '-', color=WONG['vermillion'], lw=1.0, alpha=0.30, label=r'$\ln P$ (raw)')
    posr = P_roll > 0
    axG.plot(steps[posr], np.log(P_roll[posr]), '-', color=WONG['vermillion'], lw=2.6, alpha=0.95,
             label=f'rolling mean ({cfg.roll_win})')
    if 'pfa_P' in L:
        posa = L['pfa_P'] > 0
        axG.plot(L['pfa_step'][posa], np.log(L['pfa_P'][posa]), '-', color='k', lw=1.8, alpha=0.8, label='LMP block-avg')
    axG.set_ylabel(r'$\ln P$')
    axG.set_title('(b) log piston pressure (relaxation view)')
    axG.legend(fontsize=12, loc='best')
    if 'PF' in L:
        p = L['PF']
        annotate_box(axL, f"plateau $\\langle P\\rangle$ = {fmt_val_unc(p['mean'], 0.5 * (p['hi'] - p['lo']))}\n"
                          f"last {p['frac']:.0%}, n={p['n']}, block={p['block']}", loc='upper right', fontsize=12)
    return _save(fig, cfg, 'piston_pressure_history', L['lvl'])


def fig_ratio(cfg, R, L):
    """<sigma'_zz>/<sigma'_xx> and <sigma'_zz>/<sigma'_yy> in the membrane vs step
    (increments relative to the reference when G_SUBTRACT_REF)."""
    if not L['G']:
        print('ratio figure skipped (no sigmaxx / sigmayy files)')
        return None
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    for comp, col, ls in (('xx', WONG['blue'], '-'), ('yy', WONG['vermillion'], '--')):
        G = L['G'].get(comp)
        if G is None:
            continue
        ax.plot(L['ts'], G['ratio'], ls, color=col, lw=2.4, marker='o', ms=5,
                label=fr"$\sigma'_{{zz}}/\sigma'_{{{comp}}}$   final = {G['ratio_final']:.3f}")
    ax.axvspan(L['halt_ts'], float(L['ts'][-1]), color=WONG['green'], alpha=0.10, label='plateau window')
    ax.axhline(1.0, color='k', ls=':', lw=1.0, alpha=0.6)
    ax.set_xlabel('step')
    ax.set_ylabel(r"$\langle\sigma'_{zz}\rangle_{\rm mem}/\langle\sigma'_{ii}\rangle_{\rm mem}$")
    ax.set_title(r"Network stress anisotropy $\sigma'_{zz}/\sigma'_{ii} = M/(M-2G)$"
                 + ('  (relative to $\\varepsilon=0$)' if cfg.G_SUBTRACT_REF else ''), fontsize=15)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=13, loc='best')
    return _save(fig, cfg, 'network_stress_ratio', L['lvl'])


def fig_M(cfg, R, L):
    """Longitudinal modulus: network vs piston estimate with 95 % CIs."""
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    ci = int(cfg.ci_level * 100)
    ax.errorbar([0], [L['M_net']], yerr=[[L['M_net'] - L['M_net_lo']], [L['M_net_hi'] - L['M_net']]],
                fmt='o', ms=13, color=WONG['blue'], capsize=8, lw=2.5,
                label=f"network  $M = {L['M_net']:.3f}$\n{ci}% CI [{L['M_net_lo']:.3f}, {L['M_net_hi']:.3f}]")
    ax.axhline(L['M_net'], color=WONG['blue'], ls='--', lw=1.2, alpha=0.5)
    if 'M_pist' in L:
        ax.errorbar([1], [L['M_pist']], yerr=[[L['M_pist'] - L['M_pist_lo']], [L['M_pist_hi'] - L['M_pist']]],
                    fmt='s', ms=13, color=WONG['vermillion'], capsize=8, lw=2.5,
                    label=f"piston  $M = {L['M_pist']:.3f}$\n{ci}% CI [{L['M_pist_lo']:.3f}, {L['M_pist_hi']:.3f}]")
        ax.axhline(L['M_pist'], color=WONG['vermillion'], ls='--', lw=1.2, alpha=0.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([r"network ($\sigma'_{zz}/\varepsilon$)", r'piston ($P/\varepsilon$)'], fontsize=16)
    ax.set_ylabel(r'$M$  (LJ units)')
    ax.set_title(f'Longitudinal modulus\n{cfg.RUN_ID}  |  $\\varepsilon = {L["eps"]:.3f}$', fontsize=15)
    ax.set_xlim(-0.5, 1.5)
    ax.grid(axis='y', alpha=0.3)
    ax.legend(fontsize=13, loc='best')
    return _save(fig, cfg, 'M_comparison', L['lvl'])


def fig_G(cfg, R, L):
    """Shear modulus from the lateral network stress, G = (sigma'_zz - sigma'_ii)/(2 eps),
    one estimate from xx and one from yy (should agree by symmetry)."""
    if not L['G']:
        print('G figure skipped (no sigmaxx / sigmayy files)')
        return None
    fig, ax = plt.subplots(figsize=(7.5, 6), constrained_layout=True)
    ci = int(cfg.ci_level * 100)
    xs, vals = [], []
    for k, (comp, col, mk) in enumerate((('xx', WONG['blue'], 'o'), ('yy', WONG['vermillion'], 's'))):
        G = L['G'].get(comp)
        if G is None:
            continue
        ax.errorbar([k], [G['G']], yerr=[[G['G'] - G['lo']], [G['hi'] - G['G']]], fmt=mk, ms=13, color=col,
                    capsize=8, lw=2.5, label=f"from {comp}:  $G = {G['G']:.3f}$\n{ci}% CI [{G['lo']:.3f}, {G['hi']:.3f}]")
        if np.isfinite(G['G_pist']):
            ax.plot([k], [G['G_pist']], marker='x', ms=12, mew=2.5, color=col, ls='none', alpha=0.8,
                    label=f'from {comp} with $M_\\mathrm{{piston}}$: {G["G_pist"]:.3f}')
        xs.append(k)
        vals.append(G['G'])
    if vals:
        ax.axhline(np.mean(vals), color='0.3', ls='--', lw=1.2, alpha=0.6,
                   label=f'mean of network estimates = {np.mean(vals):.3f}')
    ax.set_xticks([0, 1])
    ax.set_xticklabels([r"$G_x=(\sigma'_{zz}-\sigma'_{xx})/2\varepsilon$", r"$G_y=(\sigma'_{zz}-\sigma'_{yy})/2\varepsilon$"],
                       fontsize=14)
    ax.set_ylabel(r'$G$  (LJ units)')
    ax.set_title(f'Shear modulus from network-stress anisotropy\n{cfg.RUN_ID}  |  $\\varepsilon = {L["eps"]:.3f}$',
                 fontsize=15)
    ax.set_xlim(-0.5, 1.5)
    ax.grid(axis='y', alpha=0.3)
    ax.legend(fontsize=11, loc='best')
    return _save(fig, cfg, 'G_estimate', L['lvl'])


def fig_Dc(cfg, R, L):
    """Consolidation fit of u_z(zeta,t)/L: data (left) and data + model (right)."""
    F = L.get('Dc')
    if F is None:
        print('D_c figure skipped (no fit)')
        return None
    zff = np.linspace(0.0, 1.0, 400)
    dlL = F['DL'] / F['L']
    u0f = F['u_IC'](zff)
    u0_b = F['u_IC'](F['zf'])
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(18, 7), constrained_layout=True)
    norm = Normalize(vmin=F['ts'][F['early'][0]], vmax=F['ts'][F['early'][-1]])
    cmap = plt.cm.viridis
    for i in F['early']:
        c = cmap(norm(F['ts'][i]))
        axl.plot(F['zf'], u0_b + F['uhat'][i][F['idx']], 'o-', color=c, ms=3, alpha=0.6)
        axr.plot(F['zf'], u0_b + F['uhat'][i][F['idx']], 'o', color=c, ms=3, alpha=0.35)
        axr.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=2.0)
    for ax in (axl, axr):
        ax.plot([0, 1], [0, -dlL], 'k:', lw=1.8, label=r'affine ($t\to\infty$):  $-(\Delta L/L)\,\zeta$')
        ax.plot(zff, u0f, '-', color='0.45', lw=1.6, label=r'IC: fitted hold-onset state')
        ax.plot(0, 0, 's', color=WONG['blue'], ms=9, zorder=5, label=r'BC: $u_z(0,t)=0$ (support)')
        ax.plot(1, -dlL, 'D', color=WONG['vermillion'], ms=9, zorder=5, label=r'BC: $u_z(1,t)=-\Delta L/L$ (piston)')
        ax.set(xlabel=r'$\zeta=(z-z_\mathrm{perm})/L$', ylabel=r'$u_z/L$', xlim=(0, 1))
        ax.grid(alpha=0.3)
        ax.legend(fontsize=11)
    axl.set_title(r'$u_z(\zeta,t)/L$ -- hold snapshots (data + fitted IC offset)', fontsize=15)
    lbl = '' if cfg.DC_FREE_AMPS else rf"$\beta=p_0/M={F['beta']:.3f}$, "
    axr.set_title(rf"Consolidation fit: $D_c={F['Dc']:.3e}\ \sigma^2/\tau$, " + lbl + rf"$R^2={F['R2']:.3f}$", fontsize=15)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=[axl, axr], fraction=0.015, pad=0.04).set_label('timestep')
    fig.suptitle(f'Cooperative diffusivity fit (level _c{L["lvl"]})  |  {cfg.sim_name}', fontsize=12, fontweight='bold')
    return _save(fig, cfg, 'Dc_consolidation_fit', L['lvl'])


def fig_kappa(cfg, R, L):
    """kappa = D_c / M for the network and the piston M (CI propagated from M)."""
    K = L.get('kappa')
    if not K:
        print('kappa figure skipped (no D_c or no M)')
        return None
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    for k, (key, lab, col, mk) in enumerate((('net', r'$D_c/M_\mathrm{network}$', WONG['blue'], 'o'),
                                             ('pist', r'$D_c/M_\mathrm{piston}$', WONG['vermillion'], 's'))):
        v = K.get(key)
        if v is None:
            continue
        ax.errorbar([k], [v['k']], yerr=[[v['k'] - v['lo']], [v['hi'] - v['k']]], fmt=mk, ms=13, color=col,
                    capsize=8, lw=2.5, label=f"{lab} = {v['k']:.3e}")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['network $M$', 'piston $M$'], fontsize=16)
    ax.set_ylabel(r'$\kappa = D_c/M$  (LJ: $\sigma^5/(\epsilon\,\tau)$)')
    ax.set_title(f'Hydraulic permeability $\\kappa = D_c/M$\n$D_c = {L["Dc"]["Dc"]:.3e}$  |  '
                 f'$\\varepsilon = {L["eps"]:.3f}$', fontsize=15)
    ax.set_xlim(-0.5, 1.5)
    ax.grid(axis='y', alpha=0.3)
    ax.legend(fontsize=13, loc='best')
    return _save(fig, cfg, 'kappa', L['lvl'])


# ===========================================================================
#  10. FIGURES -- SWEEP
# ===========================================================================
def fig_volfrac_sweep(cfg, R, levels):
    """(a) mass fraction, (b) Voronoi, (c) lambda-calibrated: reference dashed
    black + one plateau-mean curve per level."""
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    fig.suptitle(f'Solvent volume fraction across the sweep ($P_{{\\rm CAL}}={cfg.P_CAL}$)  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    zx = zn(R, R['z'])
    for ax, (key, title) in zip(axes, (('phi_mf', r'(a) mass fraction $\phi_s^{\rm mf}=\rho_s/\rho_{s,0}$'),
                                       ('phi_vor', r'(b) Voronoi $\phi_s^{\rm vor}$'),
                                       ('phi_cal', r'(c) $\lambda$-calibrated $\phi_s^{\rm cal}$'))):
        txt = []
        if R.get(key) is not None:
            m, lo, hi = R[key]
            ax.fill_between(zx, lo, hi, color='0.5', alpha=0.2, lw=0)
            ax.plot(zx, m, '--', color='k', lw=2.0)
            txt.append(f'ref: {fmt_mu(m[R["interior"]])}')
        for i, L in enumerate(levels):
            if L.get(key) is None:
                continue
            m, lo, hi = L[key]
            ax.fill_between(zx, lo, hi, color=level_color(i), alpha=0.15, lw=0)
            ax.plot(zx, m, '-', color=level_color(i), lw=2.4)
            ax.axvline(zn(R, L['z_pist']), color=level_color(i), ls='-.', lw=1.2, alpha=0.6)
            txt.append(f'{L["eps"]:.2f}: {fmt_mu(m[L["interior"]])}')
        ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6)
        if np.isfinite(R['z_support']):
            ax.axvline(zn(R, R['z_support']), color='k', lw=1.5, alpha=0.85)
        finish_axes(ax, r'$\phi_s$', title)
        ax.set_ylim(0, 1.15)
        ax.legend(handles=level_handles(levels, ref=True), fontsize=11, loc='lower left')
        if txt:
            annotate_box(ax, 'in-gel means\n' + '\n'.join(txt), loc='lower right', fontsize=11)
        else:
            ax.text(0.5, 0.5, 'unavailable', ha='center', va='center', transform=ax.transAxes)
    return _save(fig, cfg, 'sweep_volfrac_profiles')


def _sweep_stress_panels(cfg, R, levels, kind, stem, suptitle):
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    fig.suptitle(suptitle, fontsize=13, fontweight='bold')
    mem_lo = min(L['z_mem_lo'] for L in levels)
    mask = R['z'] >= mem_lo + cfg.wall_margin
    for ax, comp in zip(axes, COMPONENTS):
        Rs = R['stress'].get(comp)
        lab = (r'$\sigma^{t}_{%s}(z,t)$' % comp) if kind == 't' else (r"$\sigma'_{%s}(z,t)$" % comp)
        title = f'({"abc"[COMPONENTS.index(comp)]}) ' + (
            r'total $\sigma^{t}_{%s}$' % comp if kind == 't' else r"network $\sigma'_{%s}$" % comp)
        if not any(comp in L['stress'] for L in levels):
            ax.text(0.5, 0.5, f'sigma{comp} files\nnot found', ha='center', va='center', transform=ax.transAxes)
            finish_axes(ax, lab, title)
            continue
        ref = None
        if Rs is not None:
            ref = (Rs['t_m'], Rs['t_lo'], Rs['t_hi']) if kind == 't' else (Rs['net_m'], Rs['net_lo'], Rs['net_hi'])
        overlay_levels(ax, R, levels, lambda L: L['z'], lambda L: L['ts'],
                       lambda L, c=comp: (L['stress'][c][kind] if c in L['stress'] else None), cfg,
                       ref=ref, autoscale_mask=(mask if kind == 'net' else None), ylabel=lab, title=title)
        ax.legend(handles=level_handles(levels, ref=Rs is not None), fontsize=11, loc='best')
    return _save(fig, cfg, stem)


def fig_total_stress_sweep(cfg, R, levels):
    return _sweep_stress_panels(cfg, R, levels, 't', 'sweep_total_stress_evolution',
                                f'Total stress evolution, all levels (faint = early hold, bold = plateau)  |  {cfg.sim_name}')


def fig_network_stress_sweep(cfg, R, levels):
    return _sweep_stress_panels(cfg, R, levels, 'net', 'sweep_network_stress_evolution',
                                f'Network stress evolution (Terzaghi), all levels  |  {cfg.sim_name}')


def fig_partial_stress_sweep(cfg, R, levels):
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    fig.suptitle(f'Partial and total $\\sigma_{{zz}}$ evolution, all levels  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    Rz = R['stress']['zz']
    for ax, (key, title) in zip(axes, (('s', r'(a) solvent partial $\sigma_{s,zz}$'),
                                       ('p', r'(b) polymer partial $\sigma_{p,zz}$'),
                                       ('t', r'(c) total $\sigma^t_{zz}$'))):
        m, lo, hi = mean_ci(Rz[key], cfg.ci_level)
        overlay_levels(ax, R, levels, lambda L: L['z'], lambda L: L['ts'],
                       lambda L, k=key: L['stress']['zz'][k], cfg, ref=(m, lo, hi),
                       ylabel=r'$\sigma_{zz}(z,t)$ (LJ)', title=title)
        ax.legend(handles=level_handles(levels, ref=True), fontsize=11, loc='best')
    return _save(fig, cfg, 'sweep_partial_stress_evolution')


def fig_piston_sweep(cfg, R, levels):
    fig, (axP, axG) = plt.subplots(1, 2, figsize=(18, 6), constrained_layout=True)
    fig.suptitle(f'Piston pressure histories, all levels  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    for i, L in enumerate(levels):
        if 'pf_P' not in L:
            continue
        col = level_color(i)
        st, P = L['pf_step'], L['pf_P']
        Pr = rolling_mean(P, cfg.roll_win)
        axP.plot(st, P, '-', color=col, lw=0.8, alpha=0.20)
        axP.plot(st, Pr, '-', color=col, lw=2.4, alpha=0.95)
        pos = Pr > 0
        axG.plot(st[pos], np.log(Pr[pos]), '-', color=col, lw=2.4, alpha=0.95)
        if 'PF' in L:
            axP.axvspan(L['PF']['step0'], float(st[-1]), color=col, alpha=0.06)
            axP.axhline(L['PF']['mean'], color=col, ls=':', lw=1.2, alpha=0.7)
    axP.axhline(0, color='k', ls='--', lw=0.8, alpha=0.4)
    axP.set_xlabel('step')
    axP.set_ylabel(r'$P = F_z/A$  (LJ / $\sigma^2$)')
    axP.set_title('(a) piston pressure (rolling mean; dotted = plateau)')
    axG.set_xlabel('step')
    axG.set_ylabel(r'$\ln P$')
    axG.set_title('(b) log piston pressure (relaxation view)')
    for ax in (axP, axG):
        ax.grid(alpha=0.3)
        ax.legend(handles=level_handles(levels), fontsize=11, loc='best')
    return _save(fig, cfg, 'sweep_piston_pressure_history')


def fig_ratio_sweep(cfg, R, levels):
    if not any(L['G'] for L in levels):
        print('ratio figure skipped (no sigmaxx / sigmayy files)')
        return None
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    for i, L in enumerate(levels):
        for comp, ls in (('xx', '-'), ('yy', '--')):
            G = L['G'].get(comp)
            if G is None:
                continue
            ax.plot(L['ts'], G['ratio'], ls, color=level_color(i), lw=2.2, marker='o', ms=4)
    ax.axhline(1.0, color='k', ls=':', lw=1.0, alpha=0.6)
    h = level_handles(levels) + [Line2D([0], [0], color='0.3', ls='-', lw=2, label=r"$\sigma'_{zz}/\sigma'_{xx}$"),
                                 Line2D([0], [0], color='0.3', ls='--', lw=2, label=r"$\sigma'_{zz}/\sigma'_{yy}$")]
    ax.legend(handles=h, fontsize=11, loc='best')
    ax.set_xlabel('step')
    ax.set_ylabel(r"$\langle\sigma'_{zz}\rangle/\langle\sigma'_{ii}\rangle$")
    ax.set_title(r"Network stress anisotropy $\sigma'_{zz}/\sigma'_{ii}=M/(M-2G)$, all levels", fontsize=15)
    ax.grid(alpha=0.3)
    return _save(fig, cfg, 'sweep_network_stress_ratio')


def fig_M_sweep(cfg, R, levels):
    """(a) M_network and M_piston vs applied strain; (b) piston stress vs strain."""
    fig, (axM, axPP) = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    fig.suptitle(f'Longitudinal modulus across the sweep  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    eps = np.array([L['eps'] for L in levels])
    Mn = np.array([L['M_net'] for L in levels])
    axM.errorbar(eps, Mn, yerr=[Mn - [L['M_net_lo'] for L in levels], [L['M_net_hi'] for L in levels] - Mn],
                 fmt='o-', ms=10, lw=2, color=WONG['blue'], capsize=6,
                 label=r"network  $\langle\sigma'_{zz}\rangle_\mathrm{mem}/\varepsilon$")
    hp = [L for L in levels if 'M_pist' in L]
    if hp:
        ep = np.array([L['eps'] for L in hp])
        Mp = np.array([L['M_pist'] for L in hp])
        axM.errorbar(ep, Mp, yerr=[Mp - [L['M_pist_lo'] for L in hp], [L['M_pist_hi'] for L in hp] - Mp],
                     fmt='s-', ms=10, lw=2, color=WONG['vermillion'], capsize=6,
                     label=r'piston  $(F_z/A)/\varepsilon$  (block-bootstrap CI)')
        Pp = np.array([L['P_final'] for L in hp])
        axPP.errorbar(ep, Pp, yerr=[Pp - [L['PF']['lo'] for L in hp], [L['PF']['hi'] for L in hp] - Pp],
                      fmt='o-', lw=2, ms=9, color=WONG['green'], capsize=6, label='plateau piston stress')
        if len(ep) >= 2:
            M_init = (Pp[1] - Pp[0]) / (ep[1] - ep[0])
            M_sec = (Pp[-1] - Pp[0]) / (ep[-1] - ep[0])
            axPP.set_title(fr'(b) $P$ vs $\varepsilon$:  $M_\mathrm{{init}}\approx{M_init:.3g}$, '
                           fr'$M_\mathrm{{secant}}\approx{M_sec:.3g}$', fontsize=15)
        else:
            axPP.set_title(r'(b) piston stress vs strain', fontsize=15)
    axM.set_xlabel(r'applied strain  $\varepsilon$')
    axM.set_ylabel(r'$M$  (LJ units)')
    axM.set_title('(a) longitudinal modulus, two estimates', fontsize=15)
    axM.grid(alpha=0.3)
    axM.legend(fontsize=12, loc='best')
    axPP.set_xlabel(r'applied strain  $\varepsilon$')
    axPP.set_ylabel(r'$P_\mathrm{piston}=\langle F_z\rangle/A$  (LJ)')
    axPP.grid(alpha=0.3)
    axPP.legend(fontsize=12, loc='best')
    return _save(fig, cfg, 'sweep_modulus')


def fig_G_sweep(cfg, R, levels):
    hg = [L for L in levels if L['G']]
    if not hg:
        print('G figure skipped (no sigmaxx / sigmayy files)')
        return None
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    for comp, col, mk in (('xx', WONG['blue'], 'o'), ('yy', WONG['vermillion'], 's')):
        Ls = [L for L in hg if comp in L['G']]
        if not Ls:
            continue
        e = np.array([L['eps'] for L in Ls])
        g = np.array([L['G'][comp]['G'] for L in Ls])
        ax.errorbar(e, g, yerr=[g - [L['G'][comp]['lo'] for L in Ls], [L['G'][comp]['hi'] for L in Ls] - g],
                    fmt=mk + '-', ms=10, lw=2, color=col, capsize=6,
                    label=fr"$G$ from {comp}:  $(\sigma'_{{zz}}-\sigma'_{{{comp}}})/2\varepsilon$")
        gp = np.array([L['G'][comp]['G_pist'] for L in Ls])
        if np.isfinite(gp).any():
            ax.plot(e, gp, 'x', ms=11, mew=2.5, color=col, alpha=0.8, label=f'from {comp} with $M_\\mathrm{{piston}}$')
    ax.set_xlabel(r'applied strain  $\varepsilon$')
    ax.set_ylabel(r'$G$  (LJ units)')
    ax.set_title('Shear modulus from network-stress anisotropy vs strain', fontsize=15)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=11, loc='best')
    return _save(fig, cfg, 'sweep_G_estimate')


def fig_Dc_sweep(cfg, R, levels):
    """(a) D_c vs applied strain; then one consolidation-fit panel per level."""
    hd = [L for L in levels if L.get('Dc') is not None]
    if not hd:
        print('D_c figure skipped (no level produced a fit)')
        return None
    n = len(hd)
    fig, axes = plt.subplots(1, n + 1, figsize=(7 + 6.5 * n, 6), constrained_layout=True)
    fig.suptitle(f'Cooperative diffusivity across the sweep  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    ax = axes[0]
    e = np.array([L['eps'] for L in hd])
    v = np.array([L['Dc']['Dc'] for L in hd])
    ax.plot(e, v, 'o-', color=WONG['reddishpurple'], lw=2, ms=8)
    for L in hd:   # R^2 label under each point, clipped inside the axes
        ax.annotate(f"$R^2$={L['Dc']['R2']:.2f}", (L['eps'], L['Dc']['Dc']), textcoords='offset points',
                    xytext=(0, -12), ha='center', va='top', fontsize=10, color='0.35',
                    annotation_clip=True)
    ax.margins(x=0.12, y=0.15)
    ax.axhline(cfg.DC_SLOW_REF, color='0.4', ls=':', lw=1.5, label=f'slow reference $D_c$ = {cfg.DC_SLOW_REF}')
    ax.set_xlabel(r'applied strain  $\varepsilon$')
    ax.set_ylabel(r'$D_c$  ($\sigma^2/\tau$)')
    ax.set_title(r'(a) $D_c$ vs applied strain', fontsize=15)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=11)
    zff = np.linspace(0, 1, 300)
    for k, L in enumerate(hd):
        ax = axes[k + 1]
        F = L['Dc']
        i_lvl = levels.index(L)
        cmap = plt.cm.viridis
        norm = Normalize(vmin=F['ts'][F['early'][0]], vmax=F['ts'][F['early'][-1]])
        u0_b = F['u_IC'](F['zf'])
        for i in F['early']:
            c = cmap(norm(F['ts'][i]))
            ax.plot(F['zf'], u0_b + F['uhat'][i][F['idx']], 'o', color=c, ms=2.5, alpha=0.35)
            ax.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=1.5)
        ax.plot([0, 1], [0, -F['DL'] / F['L']], 'k:', lw=1.6)
        ax.set(xlabel=r'$\zeta$', ylabel=r'$u_z/L$', xlim=(0, 1))
        ax.set_title(fr"({'bcdefgh'[k]}) $\varepsilon={L['eps']:.2f}$:  $D_c={F['Dc']:.2e}$, $R^2={F['R2']:.3f}$",
                     fontsize=14, color=level_color(i_lvl))
        ax.grid(alpha=0.3)
    return _save(fig, cfg, 'sweep_Dc')


def fig_kappa_sweep(cfg, R, levels):
    hk = [L for L in levels if L.get('kappa')]
    if not hk:
        print('kappa figure skipped (no D_c / M)')
        return None
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    for key, lab, col, mk in (('net', r'$D_c/M_\mathrm{network}$', WONG['blue'], 'o'),
                              ('pist', r'$D_c/M_\mathrm{piston}$', WONG['vermillion'], 's')):
        Ls = [L for L in hk if key in L['kappa']]
        if not Ls:
            continue
        e = np.array([L['eps'] for L in Ls])
        k = np.array([L['kappa'][key]['k'] for L in Ls])
        ax.errorbar(e, k, yerr=[k - [L['kappa'][key]['lo'] for L in Ls], [L['kappa'][key]['hi'] for L in Ls] - k],
                    fmt=mk + '-', ms=10, lw=2, color=col, capsize=6, label=lab)
    ax.set_xlabel(r'applied strain  $\varepsilon$')
    ax.set_ylabel(r'$\kappa = D_c/M$  (LJ)')
    ax.set_title(r'Hydraulic permeability $\kappa = D_c/M$ vs strain', fontsize=15)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=12, loc='best')
    return _save(fig, cfg, 'sweep_kappa')


# ===========================================================================
#  11. SUMMARIES
# ===========================================================================
def print_hold_check(cfg, levels):
    """7b of the original notebook: was each level held long enough?"""
    hd = [L for L in levels if L.get('Dc') is not None]
    if not hd:
        return
    print(f'\nHOLD-ADEQUACY CHECK  (tau_1 = L^2/(pi^2 D_c); residual = mean excess stress over the last '
          f'{cfg.plateau_frac:.0%} of the hold, i.e. the window M is read from)')
    for L in hd:
        F = L['Dc']
        print(f"  level _c{L['lvl']}:  L = {F['L']:.1f} sigma   hold T = {F['hold_T']:.0f} tau "
              f"= {F['hold_T'] / cfg.dt_lj / 1e6:.2f}M steps")
        for tag, h in F['hold_check'].items():
            flag = '' if h['ok'] else '   <-- TOO SHORT'
            print(f"     {tag:<4s} D_c={h['Dc']:.3f}:  tau_1 = {h['tau1']:.0f} tau = {h['tau1'] / cfg.dt_lj / 1e6:.2f}M steps"
                  f"  |  held {F['hold_T'] / h['tau1']:.2f} tau_1  ->  residual {h['avg'] * 100:5.2f}% "
                  f"(end {h['end'] * 100:5.2f}%)  |  {cfg.DC_TARGET_RESID:.0%} needs {h['need']:.1f} tau_1 "
                  f"= {h['need'] * h['tau1'] / cfg.dt_lj / 1e6:.1f}M steps{flag}")
        if 'M_pist' in L:
            print(f"     observed M_piston/M_network - 1 = {L['M_pist'] / L['M_net'] - 1:+.1%}   vs   "
                  f"predicted unrelaxed excess (slow D_c) = {F['hold_check']['slow']['avg']:+.1%}")
    print(f'  (triaxial_compression.lmp sizes each hold as n_tau_hold * tau_1 from the live compressed '
          f'BB thickness and Dc_est = {cfg.DC_SLOW_REF:.2f}; tau_1 ~ (1-eps)^2.)')


def print_summary(cfg, levels):
    """One line per level with the headline numbers, then the hold check."""
    ci = int(cfg.ci_level * 100)
    print(f'\nSUMMARY  ({cfg.sim_name}; {ci}% CIs)')
    print(f"{'eps':>6s} {'M_net':>18s} {'M_pist':>18s} {'G_x':>18s} {'G_y':>18s} {'D_c':>10s} {'kappa_net':>10s} {'kappa_pist':>10s}")
    for L in levels:
        def ci_(v, lo, hi):
            return f'{v:.3f} [{lo:.3f},{hi:.3f}]'
        mp = ci_(L['M_pist'], L['M_pist_lo'], L['M_pist_hi']) if 'M_pist' in L else 'n/a'
        gx = ci_(L['G']['xx']['G'], L['G']['xx']['lo'], L['G']['xx']['hi']) if 'xx' in L['G'] else 'n/a'
        gy = ci_(L['G']['yy']['G'], L['G']['yy']['lo'], L['G']['yy']['hi']) if 'yy' in L['G'] else 'n/a'
        dc = f"{L['Dc']['Dc']:.3e}" if L.get('Dc') else 'n/a'
        kn = f"{L['kappa']['net']['k']:.3e}" if L.get('kappa', {}).get('net') else 'n/a'
        kp = f"{L['kappa']['pist']['k']:.3e}" if L.get('kappa', {}).get('pist') else 'n/a'
        print(f"{L['eps']:6.3f} {ci_(L['M_net'], L['M_net_lo'], L['M_net_hi']):>18s} {mp:>18s} {gx:>18s} {gy:>18s} "
              f"{dc:>10s} {kn:>10s} {kp:>10s}")
    print_hold_check(cfg, levels)
