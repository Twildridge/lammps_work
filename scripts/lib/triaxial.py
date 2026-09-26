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
   12. two-piston (2026-09-16)      mode='permeation' loader + figures, wet-piston
                                    panels for two-piston compression runs

Physics conventions (see the Notes section at the end of either notebook):
  * total stress sigma^t = sigma_p + sigma_s (group stress/atom, kinetic term included)
  * Terzaghi: sigma' = sigma^t - p_pore, p_pore read per curve from the flat
    reservoir at z/Lz ~ baseline_zf (one scalar per stress component)
  * uniaxial strain (fixed lateral box):  sigma'_zz = M eps,  sigma'_xx = lambda eps,
    lambda = M - 2G  ->  sigma'_zz/sigma'_xx = M/(M - 2G)  ->  G = (sigma'_zz - sigma'_xx)/(2 eps)
  * M_network = (<sigma'_zz>_interior,plateau - sigma'_zz,ref) / eps_applied
    M_piston  = (<F_z/A>_plateau - P_ref) / eps_applied          (M_SUBTRACT_REF; each estimator
    subtracts its own eps = 0 reading -- the piston preload and the profile bias differ by ~0.002)
  * D_c from the two-sided consolidation fit of the polymer displacement u_z(z,t)
  * kappa = D_c / M   (Darcy permeability over viscosity, k/eta, in LJ units)
"""
from __future__ import annotations

import hashlib
import re
import sys
import stat
import time
import warnings
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
import volfrac  # noqa: E402  (scripts/lib/volfrac.py -- Voronoi + lambda calibration)
import psd      # noqa: E402  (scripts/lib/psd.py -- geometric porosity + pore-size distribution, 2026-09-24)


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
    COMP_LEVELS: list = field(default_factory=list)   # applied-strain targets, as STRINGS ("0.10"), = STRAIN_TARGETS
                                          # in the .batch (compression); leave empty for mode='permeation'
    NSTEPS: object = None                 # the <steps> tag in the file names.  None (default) -> resolved PER LEVEL
                                          # from the files present (since 2026-09-03 the .lmp tags every level with its
                                          # own auto-sized hold length, and the reference files with the first level's);
                                          # an int forces one tag everywhere (old flat-hold runs); a dict {level: tag}
                                          # (key 'ref' for the reference files) pins them by hand.
    base_dir: str = '../../flow_data_local'
    # ---- which sequence (2026-09-16) ------------------------------------
    mode: str = 'compression'             # 'compression' | 'permeation': data folder flow_data_local/<mode>/<RUN_ID>,
                                          # Expanse folder triaxial_<mode>, and which loaders/figures apply
    two_pist: bool = False                # two-piston (feed/permeate NPT-piston) run: Expanse folder
                                          # triaxial_<mode>_two_pist; wet-piston files are synced and analysed.
                                          # The file readers auto-detect multi-piston columns regardless.
    DP_PISTON: object = None              # permeation: applied dP; None -> read from the piston_pressure file
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
    res_wall_margin: float = 3.0          # TWO-PISTON runs (2026-09-22): a reservoir bin enters the pore baseline only
                                          # if it lies ENTIRELY >= this many sigma from the nearest wet-piston plane.
                                          # The 1-sigma piston sheet perturbs the solvent (depletion, then layering)
                                          # out to ~3 sigma and its pair virial is split with the piston atoms, which
                                          # the solvent profile does not count: the bin touching the feed piston reads
                                          # sigma_zz ~1.43 for a 1.50 bath and shifted every sigma' by +0.016.
    res_gel_gap_bins: int = 2             # bins skipped between the gel edge and the reservoir window
    roll_win: int = 21                    # rolling-mean window (samples) for the piston plots
    q_win_steps: int = 150000             # permeation: length (steps) of the independent windows whose
                                          # N_permeate / z_perm slopes give Q and its standard error
    # ---- volume fractions (lib/volfrac.py) -------------------------------
    VOR_ENABLE: bool = True               # False -> mass-fraction only (no trajectory pass)
    VOR_MOBILE_ONLY: bool = True          # tessellate types 1,2,3 only
    VOR_NORM: str = 'bin'                 # 'bin' (absolute) | 'mobile' (saturating) for the raw phi^vor
    P_CAL: float = 1.5                    # P_local for lambda(phi_p, P); reference-state P_target
    P_BARO: float = 1.5                   # bath / barostat pressure: the total-stress panels are drawn as sigma^t / P_BARO
    PHI_FLOOR: float = 0.02
    REF_VOR_FRAMES: int = 3               # reference frames tessellated (~20 s each)
    VOR_MAX_FRAMES: int = 4               # plateau / steady-window frames tessellated per level or run (~20 s each)
    VOR_EVO_FRAMES: int = 5               # permeation: frames tessellated evenly over the drive (evolution panels)
    P_CAL_MODE: str = 'const'             # the P handed to lambda(phi_p, P) in every z-bin (2026-09-23):
                                          #   'const'  -> P_CAL everywhere: drained equilibrium, p_pore = P_target
                                          #               throughout the gel (compression)
                                          #   'pore'   -> the pore-pressure ramp across the membrane from the measured
                                          #               feed-reservoir baseline to the measured permeate baseline
                                          #               (permeation: p_pore drops feed -> permeate; Darcy, uniform k);
                                          #               reservoirs take their own baseline
                                          #   'thermo' -> the local P_th = -1/3 tr sigma^t of the nearest stress
                                          #               snapshot (includes the network's share; noisier)
                                          # Non-finite bins fall back to the ramp; every P is clipped to the
                                          # calibration sweep's range (volfrac.lambda_of clamps there anyway).
    GEL_SHADE_TO_PISTON: bool = True      # permeation figures (2026-09-24): the grey membrane shading runs up to the
                                          # FINAL feed-piston plane instead of ending at the polymer-stress edge
                                          # (which sits ~4 bins below the piston: the remaining feed reservoir)
    # ---- pore-size distribution (lib/psd.py, 2026-09-24) ---------------------
    PSD_ENABLE: bool = True               # geometric porosity + PSD on the tessellated frames (needs VOR_ENABLE)
    PSD_R_PROBE: float = 0.5              # probe radius (sigma): void = grid points >= this far from a bead surface
    PSD_GRID: float = 0.5                 # grid spacing (sigma): 0.5 ~ 7 s/frame, 0.25 ~ 2 min/frame (the covering step)
    PSD_DMAX: float = 8.0                 # pore-diameter histogram range (2 r_probe .. PSD_DMAX) ...
    PSD_DBIN: float = 0.25                # ... and bin width (= the resolution of the covering step)
    # ---- M and G as increments from the eps = 0 reference ------------------
    M_SUBTRACT_REF: bool = True           # M = (stress - its own eps = 0 reading) / eps for BOTH estimators
                                          # (2026-09-12: the seated piston carries a real preload ~+0.0014 and the
                                          # profile method a constant ~-0.002 bias; absolute M inherits that
                                          # ~0.002 gap.  Piston reference needs piston_force_avg_ref.)
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
    RUNS_ROOT: str = None                 # None -> /home/dpollard/Documents/lammps_runs/triaxial_<mode>[_two_pist]
    TRAJ_ROOT: str = '/expanse/lustre/scratch/dpollard/temp_project/lammps_trajectories'
    # ---- derived (filled by __post_init__) -------------------------------
    DATA_DIR: Path = field(init=False)
    PLOT_DIR: Path = field(init=False)
    TRAJ_DIR: Path = field(init=False)
    sim_name: str = field(init=False)

    def __post_init__(self):
        assert self.mode in ('compression', 'permeation'), "mode must be 'compression' or 'permeation'"
        if self.mode == 'compression':
            assert self.COMP_LEVELS, 'COMP_LEVELS is empty -- list at least one level, e.g. ["0.10"]'
        self.COMP_LEVELS = [str(l) for l in self.COMP_LEVELS]
        if self.RUNS_ROOT is None:
            self.RUNS_ROOT = f'/home/dpollard/Documents/lammps_runs/triaxial_{self.mode}' + ('_two_pist' if self.two_pist else '')
        base = Path(self.base_dir)
        self.DATA_DIR = base / self.mode / self.RUN_ID
        self.PLOT_DIR = base / 'plots' / self.mode / self.RUN_ID
        self.TRAJ_DIR = base / 'traj_files.nosync'
        for d in (self.DATA_DIR, self.PLOT_DIR, self.TRAJ_DIR):
            d.mkdir(parents=True, exist_ok=True)
        self._tags = {}
        if isinstance(self.NSTEPS, dict):
            self._tags = {str(k): str(v) for k, v in self.NSTEPS.items()}
        elif self.NSTEPS is not None:
            self.NSTEPS = int(self.NSTEPS)
        # sim_name (titles, sweep plot names) carries the reference-file tag when it is known
        base = f'{self.DATANAME}_{self.INTERACTION}'
        rt = self.tag_for(None)
        self.sim_name = f'{base}_{rt}' if rt is not None else base
        found = {('ref' if l is None else l): self.tag_for(l) for l in [None] + list(self.COMP_LEVELS)}
        print(f'Config: {self.sim_name}  |  levels {self.COMP_LEVELS}\n'
              f'  data  {self.DATA_DIR}\n  plots {self.PLOT_DIR}\n  traj  {self.TRAJ_DIR}\n'
              f'  file tags (hold length): ' + ', '.join(f'{k}: {v if v else "not local yet"}' for k, v in found.items()))

    # ---- <steps> tag resolution ---------------------------------------------
    def _glob_tag(self, lvl):
        """largest <tag> among sigmazz_polymer[_ref]_<DATANAME>_<INTERACTION>_<tag>[_c<lvl>].dat
        present in DATA_DIR (non-empty files only), else None."""
        suf = '' if lvl is None else '_c' + re.escape(str(lvl))
        pat = re.compile(rf'^sigmazz_polymer{"_ref" if lvl is None else ""}_{re.escape(self.DATANAME)}_'
                         rf'{re.escape(self.INTERACTION)}_(\d+){suf}\.dat$')
        hits = sorted({int(m.group(1)) for f in self.DATA_DIR.glob('sigmazz_polymer*.dat')
                       if (m := pat.match(f.name)) and f.stat().st_size > 0})
        return str(hits[-1]) if hits else None

    def tag_for(self, lvl=None):
        """the <steps> tag of level `lvl`'s files (None -> the shared reference files),
        or None when the files are not on disk yet (sync first)."""
        key = 'ref' if lvl is None else str(lvl)
        if key in self._tags:
            return self._tags[key]
        t = str(self.NSTEPS) if isinstance(self.NSTEPS, int) else self._glob_tag(lvl)
        if t is not None:
            self._tags[key] = t
        return t

    # ---- path builders: *_ref files carry no level suffix, production files do ----
    def csuf(self, lvl):
        return '' if lvl is None else f'_c{lvl}'

    def stem(self, lvl=None):
        """<DATANAME>_<INTERACTION>_<tag>[_c<lvl>]; '*' stands in for a tag not resolved yet
        (so the name doubles as the sync glob pattern)."""
        t = self.tag_for(lvl)
        return f'{self.DATANAME}_{self.INTERACTION}_{t if t is not None else "*"}{self.csuf(lvl)}'

    def path(self, name, lvl=None, ext='dat'):
        return self.DATA_DIR / f'{name}_{self.stem(lvl)}.{ext}'

    def traj(self, name, lvl=None):
        return self.TRAJ_DIR / f'{name}_{self.stem(lvl)}.lammpstrj'

    def plot(self, stem, lvl=None):
        t = self.tag_for(lvl) or 'untagged'
        return self.PLOT_DIR / f'{stem}_{self.DATANAME}_{self.INTERACTION}_{t}{self.csuf(lvl)}.png'


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


def read_piston_table(filepath):
    """Multi-column fix print / fix ave/time file -> (names, array).  Names come
    from the '# col col ...' header the two-piston decks write (2026-09-16); None
    for the one-piston 'step value' files.  Missing file -> (None, empty)."""
    filepath = Path(filepath)
    if not filepath.exists():
        return None, np.empty((0, 0))
    names, rows = None, []
    with open(filepath) as f:
        for line in f:
            t = line.strip()
            if not t:
                continue
            if t.startswith('#'):
                toks = t.lstrip('#').split()
                if toks and not toks[0].replace('.', '').replace('-', '').isdigit():
                    names = toks
                continue
            try:
                rows.append([float(v) for v in t.split()])
            except ValueError:
                continue
    if not rows:
        return None, np.empty((0, 0))
    ncol = min(len(r) for r in rows)
    arr = np.array([r[:ncol] for r in rows])
    if names is not None and len(names) < ncol:
        names = None
    return (names[:ncol] if names else None), arr


def table_col(names, arr, prefix, default=None):
    """column whose header name starts with `prefix` (case-insensitive), or column
    `default`; None when absent."""
    if arr.size == 0:
        return None
    if names is not None:
        for j, nm in enumerate(names):
            if nm.lower().startswith(prefix.lower()):
                return arr[:, j] if j < arr.shape[1] else None
    if default is not None and default < arr.shape[1]:
        return arr[:, default]
    return None


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


def windowed_slopes(steps, y, win_steps, dt=1.0, min_pts=3):
    """Least-squares slope dy/dt (per tau) of y(steps) in consecutive NON-overlapping
    windows of win_steps, from the start of the series.  Returns (slopes, centers)
    with centers in steps; a trailing partial window is dropped."""
    steps = np.asarray(steps, float)
    y = np.asarray(y, float)
    if len(steps) < min_pts:
        return np.array([]), np.array([])
    edges = np.arange(steps[0], steps[-1] + 1e-9, win_steps)
    slopes, centers = [], []
    for a in edges:
        sel = (steps >= a) & (steps < a + win_steps)
        if sel.sum() < min_pts or steps[sel][-1] - steps[sel][0] < 0.5 * win_steps:
            continue
        slopes.append(np.polyfit(steps[sel] * dt, y[sel], 1)[0])
        centers.append(0.5 * (steps[sel][0] + steps[sel][-1]))
    return np.array(slopes), np.array(centers)


def slope_estimate(steps, y, win_steps, dt=1.0, ci=0.95):
    """Slope of y(t) as the MEAN of the slopes of independent win_steps windows, with
    the standard error over the windows (t-interval).  Returns dict(mean, se, lo, hi,
    n_win, slopes, centers, slope_all) -- slope_all is the single fit over the whole
    span (the two agree when the drift is linear)."""
    steps = np.asarray(steps, float)
    y = np.asarray(y, float)
    sl, ce = windowed_slopes(steps, y, win_steps, dt)
    out = dict(slopes=sl, centers=ce, n_win=int(len(sl)), win_steps=int(win_steps))
    out['slope_all'] = float(np.polyfit(steps * dt, y, 1)[0]) if len(steps) >= 3 else np.nan
    if len(sl) == 0:
        out.update(mean=out['slope_all'], se=np.nan, lo=np.nan, hi=np.nan)
        return out
    m = float(sl.mean())
    if len(sl) >= 2:
        se = float(sl.std(ddof=1) / np.sqrt(len(sl)))
        tcrit = float(stats.t.ppf(0.5 * (1 + ci), len(sl) - 1))
        out.update(mean=m, se=se, lo=m - tcrit * se, hi=m + tcrit * se)
    else:
        out.update(mean=m, se=np.nan, lo=np.nan, hi=np.nan)
    return out


def plateau_window_slopes(steps, y, win_steps, dt=1.0, plateau_frac_auto=0.45, ci=0.95, fracs=None, min_win=4):
    """Steady window for a SLOPE (flux) estimate: the longest trailing fraction of
    (steps, y) whose per-window slopes show no trend -- the mean slope of the first
    half of the windows agrees with that of the second half within their combined
    standard errors (2 sigma).  Returns the slope_estimate dict of that window plus
    step0, frac [, warn]."""
    steps = np.asarray(steps, float)
    y = np.asarray(y, float)
    if fracs is None:
        fracs = plateau_frac_auto * np.arange(9, 1, -1) / 9.0
    t0, t1 = float(steps[0]), float(steps[-1])
    last = None
    for f in fracs:
        sel = steps >= t1 - f * (t1 - t0)
        E = slope_estimate(steps[sel], y[sel], win_steps, dt, ci)
        if E['n_win'] < min_win:
            continue
        last = (E, f, sel)
        h = E['n_win'] // 2
        a, b = E['slopes'][:h], E['slopes'][h:]
        se_ab = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) if min(len(a), len(b)) >= 2 else np.inf
        if abs(a.mean() - b.mean()) <= 2.0 * se_ab:
            E.update(step0=float(steps[sel][0]), frac=float(f))
            return E
    if last is None:
        E = slope_estimate(steps, y, win_steps, dt, ci)
        E.update(step0=t0, frac=1.0, warn=True)
        return E
    E, f, sel = last
    E.update(step0=float(steps[sel][0]), frac=float(f), warn=True)
    return E


def sig(x, n=3):
    """x rounded to n significant figures, written without an exponent for
    1e-4 <= |x| < 1e6 (figure text: never more than 3 sig figs)."""
    x = float(x)
    if not np.isfinite(x):
        return 'n/a'
    if x == 0:
        return '0'
    e = int(np.floor(np.log10(abs(x))))
    if -4 <= e < 6:
        dec = max(0, n - 1 - e)
        return f'{round(x, n - 1 - e):.{dec}f}'
    return f'{x:.{n - 1}e}'


def fmt_step(t):
    """timestep as 535k / 5.35M (3 sig figs)."""
    t = float(t)
    if abs(t) >= 1e6:
        return sig(t / 1e6) + 'M'
    if abs(t) >= 1e3:
        return sig(t / 1e3) + 'k'
    return sig(t)


def fmt_val_unc(v, u):
    """value (3 sig figs) +/- uncertainty (2 sig figs)."""
    if np.isfinite(u) and u > 0:
        return f'{sig(v, 3)} ± {sig(u, 2)}'
    return sig(v, 3)


def fmt_ci(v, lo, hi):
    return f'{sig(v)} [{sig(lo)}, {sig(hi)}]'


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


def reservoir_mask(z, cfg, side, z_pist, z_gel, relax=True):
    """Bins of ONE solvent reservoir usable as a pore-pressure baseline at ONE instant
    (two-piston runs).  side='feed': between the gel top z_gel (+ res_gel_gap_bins
    bins) and the feed-piston plane z_pist; side='perm': between the permeate-piston
    plane and the gel bottom.  A bin counts only when it lies ENTIRELY >=
    cfg.res_wall_margin sigma from the piston plane (see the Config note: the bin
    touching the piston sheet is depleted/layered and half its wall virial sits on
    the piston atoms).  The support is solvent-permeable and needs no margin.  With
    relax=True the margin is halved when the window would otherwise be empty (a thin
    reservoir late in a permeation run); an empty result means 'no usable bin'."""
    z = np.asarray(z, float)
    h = 0.5 * cfg.binWidth
    gap = cfg.res_gel_gap_bins * cfg.binWidth
    margins = (cfg.res_wall_margin, 0.5 * cfg.res_wall_margin) if relax else (cfg.res_wall_margin,)
    for margin in margins:
        if side == 'feed':
            m = (z >= z_gel + gap) & (z + h <= z_pist - margin)
        else:
            m = (z - h >= z_pist + margin) & (z <= z_gel - gap)
        if m.sum() >= 1:
            return m
    return np.zeros(len(z), bool)


def baseline_mask(z, cfg, R=None):
    """Flat far-reservoir bins used for the pore-pressure baseline (drops the
    half-empty extreme-edge bin).  Two-piston runs (R carries finite z_feed):
    the box top is VACUUM, so the baseline is the interior of the FEED reservoir
    (reservoir_mask with the REFERENCE piston plane) instead of z/Lz ~ baseline_zf."""
    z = np.asarray(z, float)
    if R is not None and np.isfinite(R.get('z_feed', np.nan)):
        m = reservoir_mask(z, cfg, 'feed', R['z_feed'], R['z_gel_hi'])
        if m.sum() >= 1:
            return m
    zf = (z - z.min()) / (z.max() - z.min())
    return ((zf >= cfg.baseline_zf - cfg.baseline_zf_half)
            & (zf <= cfg.baseline_zf + cfg.baseline_zf_half) & (zf < 0.995))


def baseline_masks_at(z, cfg, R, ts, z_feed_at=None, z_gel_hi=None):
    """Per-snapshot pore-baseline masks, shape (n_snap, n_bins).  Two-piston runs
    rebuild the feed-reservoir window at every snapshot from the MEASURED feed-piston
    plane z_feed_at(t) (the wet pistons drift as solvent is expelled during a
    compression hold, and travel during a permeation run); the gel top is the larger
    of the reference and the level's own.  One-piston runs tile the fixed window.
    A snapshot whose window would be empty falls back to the reference window."""
    z = np.asarray(z, float)
    ts = np.asarray(ts, float)
    if np.isfinite(R.get('z_feed', np.nan)):
        zg = R['z_gel_hi'] if z_gel_hi is None else max(float(z_gel_hi), R['z_gel_hi'])
        zf = (lambda t: R['z_feed']) if z_feed_at is None else z_feed_at
        rows = [reservoir_mask(z, cfg, 'feed', zf(t), zg) for t in ts]
        return np.array([r if r.any() else R['bw'] for r in rows], bool)
    return np.tile(np.asarray(R['bw'], bool), (len(ts), 1))


def mask_span(z, m):
    """'z in [lo, hi] (n bins)' for a boolean bin mask."""
    z = np.asarray(z, float)
    m = np.asarray(m, bool)
    return f'z in [{z[m].min():.1f}, {z[m].max():.1f}] ({int(m.sum())} bins)' if m.any() else '(no bins)'


def terzaghi_split(tot_stack, bw, sd_bin, ci=0.95):
    """sigma' = sigma^t - p_pore per snapshot; p_pore = mean over the baseline bins.
    `bw` is one mask for every snapshot (1-D) or one mask PER snapshot (2-D, the
    piston-tracking windows of baseline_masks_at).
    Returns (net_stack, pore_val, pore_half, net_half)."""
    tot = np.asarray(tot_stack, float)
    bw = np.asarray(bw, bool)
    net = np.zeros_like(tot)
    pore = np.zeros(len(tot))
    pore_h = np.zeros(len(tot))
    for i in range(len(tot)):
        v = tot[i][bw[i] if bw.ndim == 2 else bw]
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
    wz = wall_z_first_frame(cfg.traj('traj_ref'), (4, 5, 6, 7))
    R['z_support'] = wz.get(4, np.nan)
    z5, z6, z7 = wz.get(5, np.nan), wz.get(6, np.nan), wz.get(7, np.nan)
    # two-piston runs (2026-09-16): types 5/6 are the wet feed/permeate pistons and the
    # loading plate is the load piston (type 7).  One-piston runs: type 5 IS the piston,
    # and z_feed/z_perm stay NaN so every one-piston code path (pore baseline at
    # z/Lz ~ baseline_zf, interior mask, autoscaling) is exactly the pre-2026-09-16 one.
    R['two_pist'] = bool(np.isfinite(z5) and np.isfinite(z6))
    R['z_feed'], R['z_perm'] = (z5, z6) if R['two_pist'] else (np.nan, np.nan)
    R['z_load'] = z7
    R['z_piston'] = z7 if np.isfinite(z7) else z5
    print(f'box (fixed): Lx={box["lx"]:.2f} Ly={box["ly"]:.2f} Lz={box["lz"]:.2f}  |  '
          f'A={R["AREA"]:.2f}  V_bin={R["V_BIN"]:.2f}')
    if R['two_pist']:
        print(f'TWO-PISTON run: support z = {R["z_support"]:.2f}  |  permeate piston (6) z = {R["z_perm"]:.2f}  |  '
              f'feed piston (5) z = {R["z_feed"]:.2f}' + (f'  |  load piston (7) z = {R["z_load"]:.2f}' if np.isfinite(R['z_load']) else ''))
    else:
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
    # upper bound: the loading plate (one-piston / two-piston compression); in permeation
    # mode the only plate above the gel is the feed piston far up in the reservoir, so the
    # gel's own top bounds the interior
    z_top = R['z_piston'] if (np.isfinite(R['z_piston']) and cfg.mode == 'compression') else R['z_gel_hi']
    R['interior'] = (R['in_gel'] & (z >= R['z_gel_lo'] + cfg.wall_margin) & (z <= z_top - cfg.wall_margin))
    print(f'gel interior (reference): z in [{R["z_gel_lo"]:.1f}, {R["z_gel_hi"]:.1f}]  ({R["in_gel"].sum()} bins)')

    # ---- reference totals, Terzaghi split per component ----------------
    R['bw'] = baseline_mask(z, cfg, R)
    if R['two_pist']:
        print(f'pore baseline window (feed reservoir, bins >= {cfg.res_wall_margin:g} sigma clear of the feed piston): '
              f'{mask_span(z, R["bw"])}  [levels re-derive it per snapshot from the measured piston plane]')
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
        _, ilo, ihi = mean_ci(S['net_m'][R['interior']], cfg.ci_level) if R['interior'].sum() > 1 else (np.nan, np.nan, np.nan)
        S['net_interior_half'] = float(0.5 * (ihi - ilo)) if np.isfinite(ihi) else 0.0
    print('reference pore baseline p_pore  ' + '  '.join(
        f'{c}: {S["pore"]:.4f}' for c, S in R['stress'].items()))
    print('reference network stress in gel interior  ' + '  '.join(
        f"sigma'_{c}: {S['net_interior']:+.4f}" for c, S in R['stress'].items()))

    # ---- reference piston load (runs since 2026-09-05 write piston_force_avg_ref) ----
    # the preload the eps = 0 state carries; must match the reference sigma'_zz step
    # above if the profile method is unbiased at zero strain
    fp = cfg.path('piston_force_avg_ref')
    if fp.exists():
        pfa = read_print_file(fp, ['step', 'Fz'])          # first value column = the loading piston (load / one-piston)
        names, tab = read_piston_table(fp)
        # two-piston: the wet pistons' zero-flux baselines P = +/- F_fluid/A
        for lab, sign in (('feed', 1.0), ('perm', -1.0)):
            col = table_col(names, tab, f'F_fluid_{lab}')
            if col is not None:
                m, lo, hi, *_ = block_bootstrap_ci(sign * col / R['AREA'], cfg.ci_level)
                R[f'P_ref_{lab}'], R[f'P_ref_{lab}_lo'], R[f'P_ref_{lab}_hi'] = float(m), float(lo), float(hi)
        if cfg.mode == 'permeation':
            R['P_ref'] = np.nan
            print('reference (zero-flux) wet-piston pressures: ' + '  '.join(
                f'{lab}: {R[f"P_ref_{lab}"]:.4f} [{R[f"P_ref_{lab}_lo"]:.4f}, {R[f"P_ref_{lab}_hi"]:.4f}]'
                for lab in ('feed', 'perm') if f'P_ref_{lab}' in R))
        else:
            P = pfa['Fz'] / R['AREA']
            m, lo, hi, blk, tau = block_bootstrap_ci(P, cfg.ci_level)
            R.update(P_ref=float(m), P_ref_lo=float(lo), P_ref_hi=float(hi), P_ref_step=pfa['step'], P_ref_P=P)
            d = R['stress']['zz']['net_interior'] - m
            print(f'reference piston preload P_ref = <F_z>/A = {m:.4f} [{lo:.4f}, {hi:.4f}] (n={len(P)}, block={blk});  '
                  f"profile sigma'_zz(ref) - P_ref = {d:+.4f}  (zero-strain bias of the profile method)")
            if 'P_ref_feed' in R:
                print('reference wet-piston (bath) pressures: ' + '  '.join(
                    f'{lab}: {R[f"P_ref_{lab}"]:.4f} [{R[f"P_ref_{lab}_lo"]:.4f}, {R[f"P_ref_{lab}_hi"]:.4f}]'
                    for lab in ('feed', 'perm') if f'P_ref_{lab}' in R))
    else:
        R['P_ref'] = np.nan
        print('reference piston preload: no piston_force_avg_ref file (runs before 2026-09-05) -> preload unknown')

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


def _support_z_func(cfg, R, lvl):
    """z(t) of the support plate over one level.  Runs since the 2026-09-18 SYMMETRIC
    drive (piston down + support up, half the closure each) write support_position_*;
    older top-only runs have none, and the support is then its fixed reference plane."""
    f = cfg.path('support_position', lvl)
    if f.exists():
        sp = np.atleast_2d(np.loadtxt(f, comments='#'))
        if sp.size and sp.shape[1] >= 2:
            st, sz = sp[:, 0], sp[:, 1]
            return (lambda t: float(np.interp(t, st, sz))), (st, sz)
    return (lambda t: R['z_support']), None


def _wet_piston_z_funcs(cfg, R, lvl=None):
    """z(t) of the feed and permeate pistons from piston_position[_c<lvl>] (columns
    z_feed / z_perm of the two-piston decks); the reference planes when the file or
    the columns are absent (one-piston runs).  Returns {feed_at, perm_at, feed_track,
    perm_track}; a track is (step, z) or None."""
    out = {}
    names, tab = read_piston_table(cfg.path('piston_position', lvl))
    for lab in ('feed', 'perm'):
        col = table_col(names, tab, f'z_{lab}') if tab.size else None
        z0 = float(R.get(f'z_{lab}', np.nan))
        if col is not None and tab.shape[0] >= 1:
            st = tab[:, 0]
            out[f'{lab}_at'] = (lambda t, st=st, c=col: float(np.interp(t, st, c)))
            out[f'{lab}_track'] = (st, col)
        else:
            out[f'{lab}_at'] = (lambda t, z0=z0: z0)
            out[f'{lab}_track'] = None
    return out


def plate_tracks(R, L):
    """{plate: (z_start, z_end)} over the level's own position files (ramp + hold):
    the support, the loading piston and -- two-piston runs -- both wet pistons.  Only
    the plates a run does not drive are expected to sit still; under the symmetric
    drive the support moves too, so nothing is assumed: every plate is reported."""
    def ends(track, const):
        return (float(track[1][0]), float(track[1][-1])) if track is not None else (float(const), float(const))
    T = {'support': ends(L.get('support_pos'), R['z_support']),
         ('load piston' if R.get('two_pist') else 'piston'): ends(L.get('piston_pos'), R['z_piston'])}
    if R.get('two_pist') and L.get('wetz'):
        T['feed piston'] = ends(L['wetz']['feed_track'], R['z_feed'])
        T['permeate piston'] = ends(L['wetz']['perm_track'], R['z_perm'])
    return T


def fmt_plates(T):
    return ' | '.join(f'{k} {a:.2f} -> {b:.2f} ({b - a:+.2f})' for k, (a, b) in T.items())


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
    # held SUPPORT plane: moves up under the symmetric drive (support_position file);
    # top-only runs (no file) keep the reference plane
    d['z_supp_held'] = _support_z_func(cfg, R, lvl)[0](d['ts'][-1])
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
    (see the D_c notes in the notebooks).  Returns a dict or None.

    The domain is bounded by the HELD plate planes.  DL is the cumulative closure of
    the plate gap since the seated reference; f_sup is the share of it the support
    took (0 for the top-only drive of runs before 2026-09-18, 1/2 for the symmetric
    drive since).  Only the affine end state / plotted IC depend on f_sup: the fit
    itself models the hold-referenced u_dat with both faces pinned, so the odd-mode
    (centre-symmetric) expansion is the same in both cases -- the symmetric drive
    just makes the hold-onset state actually symmetric about the gel centre."""
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

    # absolute u_z/L (referenced to the pre-drive state): affine end state
    # (DL/L)(f_sup - zeta) -- the support face moved UP by DL_sup, the piston face
    # DOWN by DL_pist -- plus the fitted transient
    u_model = lambda zh, t: (DL / L) * (f_sup - zh) + T(zh, t)
    u_IC = lambda zh: u_model(zh, 0.0)
    beta = np.nan if cfg.DC_FREE_AMPS else float(A[0])
    hold_T = float(t_lj[-1])
    return dict(Dc=Dc, A=A, beta=beta, R2=R2, L=L, DL=DL, DL_pist=DL_pist, DL_sup=DL_sup, f_sup=f_sup,
                gap=gap, z_perm=z_perm, z_feed=z_feed, z_sup=z_sup, z_pist=disp['z_pist_held'],
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
    say(f'\n=== level _c{lvl}  (applied strain {float(lvl):.3f};  files tagged _{cfg.tag_for(lvl)} = hold length) ===')

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
    L['support_z_at'], L['support_pos'] = _support_z_func(cfg, R, lvl)
    L['z_supp'] = L['support_z_at'](t1)                 # = R['z_support'] for top-only runs
    spf = np.abs(zz['p'][-1])
    thr = cfg.gel_thresh * float(np.nanmax(spf))
    mem = spf > thr
    L['z_mem_lo'] = float(z[mem].min()) if mem.any() else R['z_gel_lo']
    L['z_mem_hi'] = float(z[mem].max()) if mem.any() else R['z_gel_hi']
    L['in_mem'] = (z >= L['z_mem_lo']) & (z <= L['z_mem_hi'])
    L['interior'] = (L['in_mem'] & (z >= L['z_mem_lo'] + cfg.wall_margin)
                     & (z <= L['z_pist'] - cfg.wall_margin))
    say(f'  membrane (from final polymer stress): z in [{L["z_mem_lo"]:.1f}, {L["z_mem_hi"]:.1f}] '
        f'({L["in_mem"].sum()} bins);  held piston z = {L["z_pist"]:.2f}, held support z = {L["z_supp"]:.2f}'
        + ('' if L['support_pos'] is None else
           f' (moved {L["z_supp"] - R["z_support"]:+.2f} from the reference: symmetric drive)'))

    # ---- plates: what moved over this level; wet-piston tracks (2026-09-22) ----
    L['wetz'] = _wet_piston_z_funcs(cfg, R, lvl)
    L['plates'] = plate_tracks(R, L)
    say('  plates over the level (start -> end): ' + fmt_plates(L['plates']))

    # ---- Terzaghi split per component ----------------------------------
    # pore baseline window per snapshot: two-piston runs follow the measured feed-piston
    # plane (it rises as solvent is expelled), one-piston runs use the fixed window
    bw = baseline_masks_at(z, cfg, R, ts, L['wetz']['feed_at'], L['z_mem_hi'])
    L['bw'] = bw
    if R['two_pist']:
        say(f'  pore baseline window (feed reservoir, tracking the feed piston): first snapshot {mask_span(z, bw[0])}'
            f' -> last {mask_span(z, bw[-1])}')
    for comp, S in L['stress'].items():
        Rs = R['stress'].get(comp)
        sd_bin = Rs['sd_bin'] if Rs is not None else None
        if sd_bin is None or Rs['t'].shape[0] < 2 or not np.any(sd_bin > 0):
            sd_bin = np.full(len(z), float(np.nanstd(S['t'][-1][bw[-1]])))
        S['net'], S['pore'], S['pore_half'], S['net_half'] = terzaghi_split(S['t'], bw, sd_bin, cfg.ci_level)
        S['ref_net'] = (Rs['net_interior'] if (cfg.G_SUBTRACT_REF and Rs is not None) else 0.0)
        S['net_mem'] = np.array([np.nanmean(S['net'][i][L['in_mem']]) for i in range(len(ts))])
    # plateau-averaged profiles (mean over the snapshots of the last plateau_frac of the
    # hold): M, G and the figure annotations read from these, not from the last snapshot
    # alone (2026-09-05: a single snapshot averages only ~1/num_stress_curves of the hold and
    # its membrane mean scatters by ~0.002 in stress = ~0.02 in M from snapshot to snapshot)
    L['plat'] = ts >= L['halt_ts']
    if L['plat'].sum() < 1:
        L['plat'][-1] = True
    for comp, S in L['stress'].items():
        S['net_plat'] = np.nanmean(S['net'][L['plat']], axis=0)
        S['t_plat'] = np.nanmean(S['t'][L['plat']], axis=0)
    say(f"  pore pressure (zz baseline): {zz['pore'][0]:.4f} -> {zz['pore'][-1]:.4f};  "
        f"network sigma'_zz in membrane: {zz['net_mem'][0]:.4f} -> {zz['net_mem'][-1]:.4f}  "
        f"(plateau mean over {int(L['plat'].sum())} snapshots: membrane {np.nanmean(zz['net_plat'][L['in_mem']]):.4f}, "
        f"interior {np.nanmean(zz['net_plat'][L['interior']]):.4f})")

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
    # network: plateau-averaged sigma'_zz over the WALL-TRIMMED interior bins (2026-09-05;
    # was the last snapshot over the whole membrane).  The two bins that contain the wall
    # planes are missing half of the wall-polymer virial (stress/atom hands it to the
    # piston/support atoms) and read ~0.01-0.03 low, so they are excluded, as for G.
    #
    # INCREMENT (2026-09-12, M_SUBTRACT_REF): each estimator subtracts ITS OWN eps = 0 reading,
    #   M_network = (<sigma'_zz>_int,plateau - sigma'_zz,ref) / eps
    #   M_piston  = (<P>_plateau - P_ref) / eps
    # M is a slope, and the two eps = 0 readings differ: the seated piston carries a real
    # thermal-contact preload (~+0.0014 at contact_gap 1.12, ~0.5 % pre-strain) while the
    # profile method reads ~0.002 low at every strain (the total sigma_zz in the gel interior
    # sits ~0.0017 below p_res + P both at eps = 0 and at eps = 0.10).  Absolute stress / eps
    # therefore inherits a constant ~0.002 offset between the estimators (6 % of M at 0.29);
    # the increments agree within their CIs.  The absolute values are kept as *_abs.
    # The piston increment needs piston_force_avg_ref (runs since 2026-09-05); without it the
    # piston M stays absolute and L['M_pist_ref'] says so.
    eps = L['eps']
    im = L['interior'] if L['interior'].sum() >= 3 else L['in_mem']
    L['M_net_mask'] = 'interior' if im is L['interior'] else 'membrane'
    Rzz = R['stress']['zz']
    sub = bool(cfg.M_SUBTRACT_REF)
    L['M_net_ref'] = float(Rzz['net_interior']) if (sub and np.isfinite(Rzz['net_interior'])) else 0.0
    mn_abs = zz['net_plat'][im] / eps
    mn_abs = mn_abs[np.isfinite(mn_abs)]
    L['M_net_abs'], L['M_net_abs_lo'], L['M_net_abs_hi'] = mean_ci(mn_abs, cfg.ci_level)
    # increment per bin (the reference is one scalar, so it shifts the mean, not the bin scatter);
    # the reference's own CI half-width is added in quadrature to the bin-scatter interval
    ref_h = float(Rzz.get('net_interior_half', 0.0)) / eps if sub else 0.0
    m, lo, hi = mean_ci(mn_abs - L['M_net_ref'] / eps, cfg.ci_level)
    half = np.sqrt(((hi - lo) / 2) ** 2 + ref_h ** 2)
    L['M_net'], L['M_net_lo'], L['M_net_hi'] = float(m), float(m - half), float(m + half)
    L['M_net_nbins'] = len(mn_abs)
    L['M_net_final'] = float(np.nanmean(zz['net'][-1][L['in_mem']]) / eps)   # the pre-2026-09-05 estimator (absolute), for reference
    if 'PF' in L:
        p = L['PF']
        L['M_pist_abs'], L['M_pist_abs_lo'], L['M_pist_abs_hi'] = p['mean'] / eps, p['lo'] / eps, p['hi'] / eps
        L['P_final'] = p['mean']
        have_pref = sub and np.isfinite(R.get('P_ref', np.nan))
        L['M_pist_ref'] = 'measured' if have_pref else ('absent' if sub else 'off')
        L['P_ref'] = float(R['P_ref']) if have_pref else 0.0
        if have_pref:
            # CI: piston plateau bootstrap half-width (+) reference preload bootstrap half-width, in quadrature
            ph = (p['hi'] - p['lo']) / 2
            rh = (R['P_ref_hi'] - R['P_ref_lo']) / 2
            half = np.sqrt(ph ** 2 + rh ** 2) / eps
            L['M_pist'] = (p['mean'] - R['P_ref']) / eps
            L['M_pist_lo'], L['M_pist_hi'] = L['M_pist'] - half, L['M_pist'] + half
        else:
            L['M_pist'], L['M_pist_lo'], L['M_pist_hi'] = L['M_pist_abs'], L['M_pist_abs_lo'], L['M_pist_abs_hi']
    how = 'increment from eps = 0' if sub else 'absolute'
    say(f"  M_network = {L['M_net']:.4f} [{L['M_net_lo']:.4f}, {L['M_net_hi']:.4f}] "
        f"({how}; plateau mean, {L['M_net_mask']}, {L['M_net_nbins']} bins"
        + (f"; ref sigma'_zz subtracted: {L['M_net_ref']:+.4f}; absolute {L['M_net_abs']:.4f}" if sub else '')
        + f"; last-snapshot/membrane estimator: {L['M_net_final']:.4f})")
    if 'M_pist' in L:
        note = {'measured': f"ref P_ref subtracted: {L['P_ref']:+.4f}; absolute {L['M_pist_abs']:.4f}",
                'absent': 'NO piston_force_avg_ref -> piston M is ABSOLUTE (pre-2026-09-05 run)',
                'off': 'absolute'}[L['M_pist_ref']]
        say(f"  M_piston  = {L['M_pist']:.4f} [{L['M_pist_lo']:.4f}, {L['M_pist_hi']:.4f}] ({note})"
            f"   ratio M_piston/M_network = {L['M_pist'] / L['M_net']:.4f}")

    # ---- G from the lateral network stress -------------------------------
    # uniaxial strain, fixed lateral box:  sigma'_zz = M eps,  sigma'_xx = (M - 2G) eps
    #   -> G = (sigma'_zz - sigma'_xx) / (2 eps),  sigma'_zz/sigma'_xx = M/(M - 2G)
    # increments relative to the eps = 0 reference when G_SUBTRACT_REF (the
    # lateral network stress need not vanish in the reference state).
    L['G'] = {}
    dzz = zz['net'] - zz['ref_net']
    dzz_p = zz['net_plat'] - zz['ref_net']            # plateau-averaged (final G, lambda, ratio)
    for comp in ('xx', 'yy'):
        S = L['stress'].get(comp)
        if S is None:
            continue
        # The lateral total stress JUMPS at the gel boundary (only sigma_zz is continuous
        # there), so the 2-3 edge bins of the membrane carry no anisotropy information and
        # inflate the bin scatter ~10x: G and the ratio use the wall_margin-trimmed interior.
        im = L['interior'] if L['interior'].sum() >= 3 else L['in_mem']
        dxx = S['net'] - S['ref_net']
        dxx_p = S['net_plat'] - S['ref_net']
        g_bins = (dzz_p - dxx_p)[im] / (2.0 * eps)
        g_bins = g_bins[np.isfinite(g_bins)]
        Gm, Glo, Ghi = mean_ci(g_bins, cfg.ci_level)
        lam = float(np.nanmean(dxx_p[im])) / eps
        with np.errstate(invalid='ignore', divide='ignore'):
            ratio = np.array([np.nanmean(dzz[i][im]) / np.nanmean(dxx[i][im]) for i in range(len(ts))])
        # ---- ratio with propagated uncertainty (item: error bars on sigma'_zz/sigma'_ii) ----
        # each membrane mean carries: the t-interval of its bin scatter (as M_network),
        # the pore baseline uncertainty of that snapshot (common to all bins, so it does
        # not average out), and the reference-state increment uncertainty; combined in
        # quadrature, then  d(a/b)/(a/b) = sqrt((da/a)^2 + (db/b)^2).
        def _mem_mean_err(D, Sx, Rx):
            m = np.array([np.nanmean(D[i][im]) for i in range(len(ts))])
            half = np.array([0.5 * (mean_ci(D[i][im], cfg.ci_level)[2]
                                    - mean_ci(D[i][im], cfg.ci_level)[1]) for i in range(len(ts))])
            ref_h = (Rx['net_interior_half'] if (cfg.G_SUBTRACT_REF and Rx is not None) else 0.0)
            return m, np.sqrt(half ** 2 + np.asarray(Sx['pore_half']) ** 2 + ref_h ** 2)
        a, da = _mem_mean_err(dzz, zz, R['stress'].get('zz'))
        b, db = _mem_mean_err(dxx, S, R['stress'].get(comp))
        with np.errstate(invalid='ignore', divide='ignore'):
            ratio_err = np.abs(ratio) * np.sqrt((da / a) ** 2 + (db / b) ** 2)
        # final ratio from the plateau-averaged profiles; its error from the plateau snapshots' errors
        with np.errstate(invalid='ignore', divide='ignore'):
            ratio_plat = float(np.nanmean(dzz_p[im]) / np.nanmean(dxx_p[im]))
        ratio_plat_err = float(np.sqrt(np.nanmean(ratio_err[L['plat']] ** 2) / max(int(L['plat'].sum()), 1)))
        L['G'][comp] = dict(G=Gm, lo=Glo, hi=Ghi, lam=lam, ratio=ratio, ratio_err=ratio_err,
                            nbins=len(g_bins), ratio_final=ratio_plat, ratio_final_err=ratio_plat_err)
        say(f"  G from {comp}: {Gm:.4f} [{Glo:.4f}, {Ghi:.4f}]   lambda_{comp} = {lam:.4f}   "
            f"sigma'_zz/sigma'_{comp} (plateau) = {ratio_plat:.3f} ± {ratio_plat_err:.3f}"
            + (f"   (ref sigma'_{comp} subtracted: {S['ref_net']:+.4f})" if cfg.G_SUBTRACT_REF else ''))

    # ---- two-piston: wet-piston bath check + solvent expelled (2026-09-16) ----
    L['wet'] = load_wet_pistons(cfg, R, lvl, plat_from=L['halt_ts'])
    if L['wet'] is not None:
        W = L['wet']
        say('  bath (wet pistons, plateau means): ' + '  '.join(
            f"{k}: {W['plat'][k]:.4f}" for k in ('P_feed_meas', 'P_perm_meas', 'P_load_meas') if k in W['plat'])
            + (f"   applied {W['plat'].get('P_feed_app', np.nan):.3f}" if 'P_feed_app' in W['plat'] else '')
            + (f"   expelled dV_total = {W['dV_total'][-1]:.1f} sigma^3 (= {W['dV_total'][-1] / R['AREA']:.2f} sigma of feed rise)"
               if W.get('dV_total') is not None else ''))

    # ---- D_c and kappa ---------------------------------------------------
    L['Dc'] = fit_Dc(cfg, R, load_disp(cfg, R, lvl))
    if L['Dc'] is None:
        say('  NOTE: no D_c (missing disp_z_polymer / piston_position, no support plane, or too few bins)')
    else:
        F = L['Dc']
        say(f"  D_c = {F['Dc']:.4e} sigma^2/tau  (R^2 = {F['R2']:.3f};  L = {F['L']:.2f}, "
            f"DL/L = {F['DL'] / F['L']:.4f} [support share {F['f_sup']:.2f}], hold = {F['hold_T']:.0f} tau)")
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
    """One streamed pass over `traj` -> dict(ts, phi_vor, phi_vor_mobile) on the
    R['z'] grid, cached on (file, frames, knobs).  The calibration is applied
    afterwards by _calibrate (so P_CAL / P_CAL_MODE never force a re-tessellation)."""
    ts_want = [int(t) for t in np.atleast_1d(ts_want)]
    key = (str(traj), tuple(ts_want), cfg.VOR_MOBILE_ONLY, cfg.VOR_NORM, cfg.binWidth)
    if key in _VF_CACHE:
        print(f'    {label}cached: {len(_VF_CACHE[key]["ts"])} frame(s), no re-read')
        return _VF_CACHE[key]
    # ---- on-disk cache (2026-09-26): <DATA_DIR>/vor_cache/<traj stem>__<hash>.npz keyed on the
    #      frames + knobs, so a kernel restart or a reload of this module never re-tessellates.
    #      Only phi profiles are stored; the polymer frames the PSD pass needs are re-streamed
    #      on demand (_frames_of).  Delete the folder to force a recompute.
    cdir = Path(cfg.DATA_DIR) / 'vor_cache'
    ck = hashlib.md5(repr((Path(traj).name, tuple(ts_want), cfg.VOR_MOBILE_ONLY, cfg.VOR_NORM, cfg.binWidth,
                           np.round(np.asarray(R['z'], float), 6).tolist())).encode()).hexdigest()[:12]
    cfile = cdir / f'{Path(traj).stem}__{ck}.npz'
    if cfile.exists():
        try:
            d = np.load(cfile, allow_pickle=False)
            out = {'ts': d['ts'], 'phi_vor': (d['phi_vor'] if d['phi_vor'].ndim == 2 else None),
                   'phi_vor_mobile': (d['phi_vor_mobile'] if d['phi_vor_mobile'].ndim == 2 else None),
                   'frames': {}, 'traj': str(traj)}
            print(f'    {label}disk-cached ({cfile.parent.name}/{cfile.name}): {len(out["ts"])} frame(s), no re-read')
            _VF_CACHE[key] = out
            return out
        except Exception as e:
            print(f'    {label}vor_cache unreadable ({type(e).__name__}) -> recomputing')
    print(f'    {label}streaming {Path(traj).name} for {len(ts_want)} frame(s): Voronoi (~20 s/frame) ...')
    fr = volfrac.stream_traj_frames(traj, ts_want)
    missing = [t for t in ts_want if t not in fr]
    if missing:
        print(f'    NOTE: no traj frame at {missing} -- dropped')
    ts_ok = [t for t in ts_want if t in fr]
    out = {'ts': np.array(ts_ok, float), 'phi_vor': None, 'phi_vor_mobile': None, 'frames': {}, 'traj': str(traj)}
    if not ts_ok:
        return out
    pv, pm, ok = [], [], []
    for t in ts_ok:
        box, typ, xyz = fr[t]
        keep = np.isin(typ, psd.POLYMER_TYPES)          # polymer beads only: the PSD pass (lib/psd.py) runs on
        out['frames'][t] = (box, typ[keep], xyz[keep])  # these without re-reading the dump
        try:
            zc, ph_b, ph_m = volfrac.phi_voronoi_frame(box, typ, xyz, cfg.binWidth,
                                                       mobile_only=cfg.VOR_MOBILE_ONLY, norm='both')
        except Exception as e:                       # e.g. a frame without mobile atoms
            print(f'      ts {t}: tessellation failed ({e}) -- frame dropped')
            continue
        ph = ph_b if cfg.VOR_NORM == 'bin' else ph_m
        pv.append(np.interp(R['z'], zc, ph, left=np.nan, right=np.nan))
        pm.append(np.interp(R['z'], zc, ph_m, left=np.nan, right=np.nan))
        ok.append(t)
        print(f'      ts {t}: phi^vor max = {np.nanmax(pv[-1]):.3f}')
    del fr
    if not ok:
        return out
    out['ts'] = np.array(ok, float)
    out['phi_vor'] = np.array(pv)
    out['phi_vor_mobile'] = np.array(pm)
    _VF_CACHE[key] = out
    try:
        cdir.mkdir(parents=True, exist_ok=True)
        np.savez(cfile, ts=out['ts'], phi_vor=out['phi_vor'], phi_vor_mobile=out['phi_vor_mobile'])
        print(f'    {label}saved to {cfile.parent.name}/{cfile.name}')
    except Exception as e:
        print(f'    {label}(vor_cache not written: {type(e).__name__}: {e})')
    return out


def _frames_of(vf, ts=None):
    """Polymer frames {ts: (box, types, xyz)} of a _volume_fractions result; re-streamed
    from vf['traj'] when the result came from the on-disk cache (frames are not stored)."""
    if vf is None:
        return {}
    ts = [int(t) for t in (vf['ts'] if ts is None else np.atleast_1d(ts))]
    if all(t in vf.get('frames', {}) for t in ts):
        return vf['frames']
    traj = vf.get('traj')
    if not traj or not Path(traj).exists():
        return vf.get('frames', {})
    print(f'    re-streaming {Path(traj).name} for the PSD pass ({len(ts)} frame(s); phi came from the disk cache)')
    fr = volfrac.stream_traj_frames(traj, ts)
    for t, (box, typ, xyz) in fr.items():
        keep = np.isin(typ, psd.POLYMER_TYPES)
        vf.setdefault('frames', {})[t] = (box, typ[keep], xyz[keep])
    return vf['frames']


def _calib_range(R):
    """(P_min, P_max) of the calibration sweep, or (-inf, inf) without an artifact."""
    c = R.get('CALIB')
    if c is None:
        return -np.inf, np.inf
    Ps = np.asarray(c['coeffs']['pressures'], float)
    return float(Ps.min()), float(Ps.max())


def _p_local_fn(cfg, R, D):
    """-> f(ts) = P_local(z) on R['z']: the pressure handed to lambda(phi_p, P) per
    z-bin for the state D (R itself, a level L or the permeation dict P) under
    cfg.P_CAL_MODE -- see Config.  'pore' builds the ramp from the measured
    reservoir baselines (sigma_zz of the reservoir interior, positive under
    compression, i.e. +P_res): the feed baseline above the membrane, the permeate
    baseline (when the run has one) below it, linear in between.  Every mode
    clips to the calibration's pressure range and fills non-finite bins from
    the ramp."""
    z = R['z']
    const = np.full(len(z), float(cfg.P_CAL))
    mode = cfg.P_CAL_MODE
    if mode not in ('const', 'pore', 'thermo'):
        raise ValueError(f"P_CAL_MODE must be 'const', 'pore' or 'thermo' (got {mode!r})")
    if mode == 'const' or D is None or 'stress' not in D or 'zz' not in D['stress']:
        return lambda t: const.copy()
    zz = D['stress']['zz']
    Pmin, Pmax = _calib_range(R)
    if D is R:                                          # reference: one baseline, uniform
        pf_at = pp_at = (lambda t, v=float(zz['pore']): v)
        zlo, zhi = R['z_gel_lo'], R['z_gel_hi']
    else:
        ts = np.asarray(D['ts'], float)
        pore = np.asarray(zz['pore'], float)
        pf_at = lambda t, ts=ts, pore=pore: float(np.interp(t, ts, pore))
        if zz.get('pore_perm') is not None:
            pp = np.asarray(zz['pore_perm'], float)
            pp_at = lambda t, ts=ts, pp=pp: float(np.interp(t, ts, pp))
        else:
            pp_at = pf_at
        zlo, zhi = D.get('z_mem_lo', R['z_gel_lo']), D.get('z_mem_hi', R['z_gel_hi'])

    def ramp(t):
        pf, pp = pf_at(t), pp_at(t)
        w = np.clip((z - zlo) / max(zhi - zlo, 1e-9), 0.0, 1.0)
        return pp + w * (pf - pp)

    f = ramp
    if mode == 'thermo':
        Pt = _tr3(D['stress'], 't')
        if Pt is None:
            print("NOTE: P_CAL_MODE='thermo' needs the xx and yy profiles -> falling back to the pore-pressure ramp")
        elif D is R:
            Pm = np.nanmean(Pt, axis=0)
            f = lambda t, Pm=Pm: Pm.copy()
        else:
            ts = np.asarray(D['ts'], float)
            f = lambda t, ts=ts, Pt=Pt: np.asarray(Pt[int(np.argmin(np.abs(ts - t)))], float).copy()

    def g(t):
        p = np.asarray(f(t), float).copy()
        bad = ~np.isfinite(p)
        if bad.any():
            p[bad] = ramp(t)[bad]
        return np.clip(p, Pmin, Pmax)
    return g


def _calibrate(cfg, R, vf, p_of_t):
    """Apply lambda(phi_p, P_local) to the frames in vf -> (phi_cal, P_local, lam)
    stacks [n_frames, nz], or (None, None, None) without a calibration artifact.
    Solvent primary (phi_p^vor = 1 - phi_f^vor, 'mobile' norm), polymer by complement."""
    if R.get('CALIB') is None or vf.get('phi_vor_mobile') is None:
        return None, None, None
    Pl = np.array([p_of_t(t) for t in vf['ts']])
    pm = vf['phi_vor_mobile']
    lam = volfrac.lambda_of(1.0 - pm, Pl, R['CALIB'])
    cal = volfrac.phi_calibrated(pm, Pl, R['CALIB'])
    return cal, Pl, lam


def _load_calib(R):
    if 'CALIB' not in R:
        try:
            R['CALIB'] = volfrac.load_calibration()
            print(f'calibration loaded ({volfrac.CALIBRATION_JSON.name}): lambda coeffs {R["CALIB"]["coeffs"]}')
        except FileNotFoundError:
            R['CALIB'] = None
            print('NOTE: no calibration artifact -> phi_cal skipped '
                  '(generate scripts/calibration/calibration_lambda.json with calibration_analysis.ipynb)')


def _reference_volume_fractions(cfg, R):
    """Voronoi + calibrated phi for the reference state: REF_VOR_FRAMES frames evenly
    over the reference window (mean + CI).  P_local from _p_local_fn(cfg, R, R)."""
    if R.get('phi_vor') is not None:
        return
    tr = cfg.traj('traj_ref')
    if tr.exists() and 'dens_ts' in R:
        want, _ = subsample(R['dens_ts'], R['dens_ts'][:, None], cfg.REF_VOR_FRAMES)
        vf = _volume_fractions(cfg, R, tr, want, 'reference: ')
        R['vf_ref'] = vf                                 # frames reused by add_perm_psd
        if vf['phi_vor'] is not None:
            R['phi_vor'] = mean_ci(vf['phi_vor'], cfg.ci_level)
            R['phi_vor_frames'] = vf['ts']
            cal, Pl, lam = _calibrate(cfg, R, vf, _p_local_fn(cfg, R, R))
            if cal is not None:
                R['phi_cal'] = mean_ci(cal, cfg.ci_level)
                R['P_cal'] = mean_ci(Pl, cfg.ci_level)
                R['lam'] = mean_ci(lam, cfg.ci_level)
    else:
        print(f'NOTE: {tr.name} or reference density missing -> reference Voronoi skipped')


def _vf_line(tag, D, mask):
    parts = []
    for k, lab in (('phi_mf', 'mass-frac'), ('phi_vor', 'Voronoi'), ('phi_cal', 'calibrated')):
        if D.get(k) is not None:
            parts.append(f'{lab} {fmt_mu(D[k][0][mask])}')
    print(f'  {tag:<22s} ' + '   '.join(parts))


def add_volume_fractions(cfg, R, levels=()):
    """Voronoi + lambda-calibrated solvent volume fractions for the reference
    state (R) and each level in `levels`.  Reference: REF_VOR_FRAMES frames
    evenly over the reference window (mean + CI).  Level: VOR_MAX_FRAMES frames
    inside the plateau window (plateau mean + CI).  Mass-fraction profiles were
    already set by load_reference / load_level.  The calibration pressure per bin
    follows cfg.P_CAL_MODE ('const' = P_CAL, the drained-equilibrium choice)."""
    _load_calib(R)
    if not cfg.VOR_ENABLE:
        print('VOR_ENABLE=False -> Voronoi / calibrated phi skipped')
        return
    _reference_volume_fractions(cfg, R)
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
            cal, Pl, lam = _calibrate(cfg, R, vf, _p_local_fn(cfg, R, L))
            if cal is not None:
                L['phi_cal'] = mean_ci(cal, cfg.ci_level)
                L['P_cal'] = mean_ci(Pl, cfg.ci_level)
                L['lam'] = mean_ci(lam, cfg.ci_level)
    # ---- printed in-gel means ------------------------------------------
    print(f'in-gel solvent volume fractions (interior, wall_margin trimmed; P_CAL_MODE = {cfg.P_CAL_MODE}):')
    _vf_line('reference (eps = 0)', R, R['interior'])
    for L in levels:
        if L is not None:
            _vf_line(f'compressed eps={L["eps"]:.2f}', L, L['interior'])


def traj_timesteps(traj_file):
    """The timesteps a lammpstrj holds (header scan only; the atom lines are skipped)."""
    out = []
    with open(traj_file) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith('ITEM: TIMESTEP'):
                continue
            ts = int(f.readline()); f.readline(); n = int(f.readline())
            out.append(ts)
            for _ in range(5 + n):            # BOX BOUNDS header, 3 bounds, ATOMS header, n atoms
                f.readline()
    return np.array(out, float)


def add_perm_volume_fractions(cfg, R, P):
    """Mass-fraction, Voronoi and lambda-calibrated solvent volume fractions for a
    PERMEATION run (2026-09-23).  Mass fraction from every density snapshot
    (steady mean = trailing plateau_frac window); Voronoi on VOR_EVO_FRAMES frames
    evenly over the drive plus VOR_MAX_FRAMES frames inside the steady window, one
    streamed pass over traj_stress (frames chosen from the ones the dump holds).
    The calibration pressure per bin follows cfg.P_CAL_MODE -- 'pore' (the pore-
    pressure ramp feed -> permeate across the membrane) is the choice under flow.
    Fills P['phi_mf' | 'phi_vor' | 'phi_cal' | 'P_cal' | 'lam'] (steady mean + CI)
    and P['vf'] (the per-frame stacks for the evolution figure)."""
    _load_calib(R)
    z = R['z']
    # ---- mass fraction ----------------------------------------------------
    P['phi_mf'] = None
    if 'dens_m' in P and np.isfinite(R.get('rho_s0', np.nan)):
        mf = np.array([np.interp(z, P['dens_z'], row) for row in P['dens_m']]) / R['rho_s0']
        P['mf_stack'], P['mf_ts'] = mf, np.asarray(P['dens_ts'], float)
        pl = P['mf_ts'] >= P['halt_ts']
        P['phi_mf'] = mean_ci(mf[pl] if pl.any() else mf[-1:], cfg.ci_level)
    else:
        print('NOTE: production density or reference rho_s,0 missing -> mass-fraction phi skipped')
    P['phi_vor'] = P['phi_cal'] = P['vf'] = None
    P['P_cal_mode'] = cfg.P_CAL_MODE
    if not cfg.VOR_ENABLE:
        print('VOR_ENABLE=False -> Voronoi / calibrated phi skipped')
        return
    _reference_volume_fractions(cfg, R)
    tp = cfg.traj('traj_stress')
    if not tp.exists():
        print(f'NOTE: {tp.name} missing -> production Voronoi skipped')
        return
    avail = traj_timesteps(tp)
    prod = avail[avail >= P['evol_from']]
    if not len(prod):
        prod = avail
    if not len(prod):
        print(f'NOTE: {tp.name} holds no frames -> production Voronoi skipped')
        return
    evo, _ = subsample(prod, prod[:, None], cfg.VOR_EVO_FRAMES)
    steady = prod[prod >= P['halt_ts']]
    if len(steady):
        st_sel, _ = subsample(steady, steady[:, None], cfg.VOR_MAX_FRAMES)
    else:
        st_sel = prod[-min(2, len(prod)):]
        print(f'NOTE: no traj frame inside the steady window (step >= {P["halt_ts"]}) -> '
              f'the last {len(st_sel)} frame(s) stand in for the steady state')
    want = np.unique(np.concatenate([evo, st_sel]))
    vf = dict(_volume_fractions(cfg, R, tp, want, 'permeation: '))
    if vf['phi_vor'] is None:
        return
    cal, Pl, lam = _calibrate(cfg, R, vf, _p_local_fn(cfg, R, P))
    vf.update(phi_cal=cal, P_local=Pl, lam=lam)
    P['vf'] = vf
    sm = np.isin(vf['ts'], st_sel)
    if not sm.any():
        sm[-1] = True
    P['vor_steady_ts'] = vf['ts'][sm]
    P['phi_vor'] = mean_ci(vf['phi_vor'][sm], cfg.ci_level)
    if cal is not None:
        P['phi_cal'] = mean_ci(cal[sm], cfg.ci_level)
        P['P_cal'] = mean_ci(Pl[sm], cfg.ci_level)
        P['lam'] = mean_ci(lam[sm], cfg.ci_level)
    # ---- printed in-gel means --------------------------------------------
    print(f'in-gel solvent volume fractions (interior, wall_margin trimmed; P_CAL_MODE = {cfg.P_CAL_MODE}):')
    _vf_line('reference (zero flux)', R, R['interior'])
    _vf_line('steady permeation', P, P['interior'])
    if P.get('P_cal') is not None:
        print(f'  calibration pressure over the membrane (steady): {fmt_mu(P["P_cal"][0][P["in_mem"]])}   '
              f'lambda in the interior: {fmt_mu(P["lam"][0][P["interior"]])}   '
              f'(frames {", ".join(fmt_step(t) for t in P["vor_steady_ts"])})')


# ---------------------------------------------------------------------------
#  6b. geometric porosity + pore-size distribution (lib/psd.py, 2026-09-24)
# ---------------------------------------------------------------------------
def _psd_state(cfg, R, frames, ts, z_lo, z_hi, interior, label=''):
    """psd.psd_frame on the frames `ts` of `frames` ({ts: (box, types, xyz)} of polymer
    beads, kept by _volume_fractions) -> dict(ts, por, d_mean, d_med [mean_ci on
    R['z']], d_edges, regions{name: (D, mean, lo, hi)}, D_reg{name: (mean, lo, hi)});
    None without frames.  The grid spans the membrane plus two bins on either side
    (porosity -> 1 in the reservoirs); regions = the interior and its feed-side /
    permeate-side halves."""
    ts = [int(t) for t in np.atleast_1d(ts) if int(t) in frames]
    if not ts:
        return None
    z, bw = R['z'], cfg.binWidth
    edges = np.concatenate([z - 0.5 * bw, [z[-1] + 0.5 * bw]])
    d_edges = np.arange(2.0 * cfg.PSD_R_PROBE, cfg.PSD_DMAX + 1e-9, cfg.PSD_DBIN)
    por, dm, dmed, hists, nvoid = [], [], [], [], []
    for t in ts:
        box, typ, xyz = frames[t]
        t0 = time.time()
        out = psd.psd_frame(box, typ, xyz, z_lo - 2 * bw, z_hi + 2 * bw, edges, h=cfg.PSD_GRID,
                            r_probe=cfg.PSD_R_PROBE, d_edges=d_edges)
        print(f'      {label}ts {t}: interior porosity {np.nanmean(out["por"][interior]):.3f}, '
              f'<D> = {np.nanmean(out["d_mean"][interior]):.2f} sigma  ({time.time() - t0:.0f} s)')
        por.append(out['por']); dm.append(out['d_mean']); dmed.append(out['d_med'])
        hists.append(out['hist']); nvoid.append(out['n_void'])
    S = dict(ts=np.array(ts, float), por=mean_ci(por, cfg.ci_level), d_mean=mean_ci(dm, cfg.ci_level),
             d_med=mean_ci(dmed, cfg.ci_level), d_edges=d_edges, regions={}, D_reg={})
    zmid = 0.5 * (z_lo + z_hi)
    Dc = 0.5 * (d_edges[:-1] + d_edges[1:])
    for name, mask in (('interior', interior), ('feed half', interior & (z >= zmid)), ('permeate half', interior & (z < zmid))):
        dens = np.array([psd.region_density(h, nv, mask, d_edges)[1] for h, nv in zip(hists, nvoid)])
        S['regions'][name] = (Dc,) + tuple(mean_ci(dens, cfg.ci_level))
        cnt = np.array([h[mask].sum(axis=0) for h in hists])
        with np.errstate(invalid='ignore', divide='ignore'):
            Dm = (cnt * Dc).sum(axis=1) / cnt.sum(axis=1)
        S['D_reg'][name] = tuple(float(v) for v in mean_ci(Dm[:, None], cfg.ci_level))
    return S


def add_perm_psd(cfg, R, P):
    """Geometric porosity and pore-size distribution (lib/psd.py) of the zero-flux
    reference and the steady permeation state, on the frames the Voronoi pass already
    streamed (R['vf_ref']: the reference frames; P['vf']: the steady-window frames), so
    the dump is not re-read.  Fills R['psd'] and P['psd'] (None when unavailable) and
    prints the interior means next to the hydraulic mesh size xi = sqrt(kappa)
    (Brinkman: kappa ~ xi^2).  Needs add_perm_volume_fractions with VOR_ENABLE."""
    P['psd'] = None
    if not (cfg.PSD_ENABLE and cfg.VOR_ENABLE):
        print('PSD_ENABLE/VOR_ENABLE False -> pore-size distribution skipped')
        return
    print(f'pore-size distribution (grid {cfg.PSD_GRID} sigma, r_probe = {cfg.PSD_R_PROBE}, largest-included-sphere covering):')
    vr = R.get('vf_ref')
    if R.get('psd') is None and vr is not None and _frames_of(vr):
        R['psd'] = _psd_state(cfg, R, vr['frames'], vr['ts'], R['z_gel_lo'], R['z_gel_hi'], R['interior'], 'reference: ')
    vf = P.get('vf')
    if vf is None or not _frames_of(vf, P.get('vor_steady_ts', vf['ts']) if vf else None):
        print('  NOTE: no tessellated production frames (traj_stress missing?) -> steady-state PSD skipped')
        return
    P['psd'] = _psd_state(cfg, R, vf['frames'], P.get('vor_steady_ts', vf['ts']), P['z_mem_lo'], P['z_mem_hi'],
                          P['interior'], 'steady: ')
    for tag, D, S in (('reference (zero flux)', R, R.get('psd')), ('steady permeation', P, P['psd'])):
        if S is None:
            continue
        print(f"  {tag:<22s} interior porosity {fmt_mu(S['por'][0][D['interior']])}   <D_pore> {fmt_mu(S['d_mean'][0][D['interior']])} sigma"
              f"   (feed half {S['D_reg']['feed half'][0]:.2f}, permeate half {S['D_reg']['permeate half'][0]:.2f});  frames "
              + ', '.join(fmt_step(t) for t in S['ts']))
    k = (P.get('flux') or {}).get('k') or {}
    kk = 'N_measured' if 'N_measured' in k else (next(iter(k)) if k else None)
    if kk is not None:
        xi = np.sqrt(max(k[kk]['k'], 0.0))
        print(f"  hydraulic mesh size xi = sqrt(kappa[{kk}]) = {xi:.2f} sigma  (Brinkman: kappa ~ xi^2)  vs geometric "
              f"<D_pore> = {np.nanmean(P['psd']['d_mean'][0][P['interior']]):.2f} sigma")


def fig_perm_psd(cfg, R, P):
    """(a) geometric porosity profile (reference dashed vs steady permeation, 95 % bands
    over frames) with phi_s^cal overlaid for comparison -- the porosity is the volume a
    solvent-sized probe can reach, NOT the thermodynamic solvent fraction; (b) mean pore
    diameter profile with 95 % bands; (c) the volume-weighted PSD of the membrane
    interior, reference vs steady, plus the feed-side and permeate-side halves.
    Needs add_perm_psd."""
    S, S0 = P.get('psd'), R.get('psd')
    if S is None and S0 is None:
        print('PSD figure skipped (run add_perm_psd first; needs tessellated frames)')
        return None
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    fig.suptitle(f'Geometric porosity and pore-size distribution of the network (largest included sphere; grid '
                 f'{cfg.PSD_GRID} $\\sigma$, $r_{{\\rm probe}}={cfg.PSD_R_PROBE}$)  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    zx = zn(R, R['z'])
    n_fr = lambda X: f' ({len(X["ts"])} frames)'
    # (a) porosity vs phi_s^cal
    ax = axes[0]
    if S0 is not None:
        m, lo, hi = S0['por']
        ax.fill_between(zx, lo, hi, color='0.5', alpha=0.25, lw=0, zorder=1)
        ax.plot(zx, m, '--', color='k', lw=2.0, alpha=0.9, zorder=2, label=r'$\epsilon_g$ reference (zero flux)' + n_fr(S0))
    if S is not None:
        m, lo, hi = S['por']
        ax.fill_between(zx, lo, hi, color=WONG['blue'], alpha=0.25, lw=0, zorder=3)
        ax.plot(zx, m, '-', color=WONG['blue'], lw=2.8, zorder=4, label=r'$\epsilon_g$ steady permeation' + n_fr(S))
    if R.get('phi_cal') is not None:
        ax.plot(zx, R['phi_cal'][0], '--', color=WONG['green'], lw=1.4, alpha=0.8, zorder=3, label=r'$\phi_s^{\rm cal}$ reference (thermodynamic)')
    if P.get('phi_cal') is not None:
        ax.plot(zx, P['phi_cal'][0], '-', color=WONG['green'], lw=1.6, alpha=0.9, zorder=4, label=r'$\phi_s^{\rm cal}$ steady (thermodynamic)')
    ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6)
    shade_gel(ax, R, P)
    mark_walls(ax, R, P)
    finish_axes(ax, r'geometric porosity $\epsilon_g$', f'(a) porosity ($r_{{\\rm probe}}={cfg.PSD_R_PROBE}$) vs $\\phi_s^{{\\rm cal}}$')
    ax.set_ylim(0, 1.15)
    smart_legend(ax, fontsize=11)
    txt = []
    if S0 is not None:
        txt.append(f"reference: {fmt_mu(S0['por'][0][R['interior']])}")
    if S is not None:
        txt.append(f"steady: {fmt_mu(S['por'][0][P['interior']])}")
        if S0 is not None:
            txt.append(f"$\\Delta\\epsilon_g$ = {np.nanmean(S['por'][0][P['interior']]) - np.nanmean(S0['por'][0][R['interior']]):+.2g}")
    annotate_box(ax, r'$\epsilon_g$ in gel interior' + '\n' + '\n'.join(txt), loc='lower left', fontsize=12)
    # (b) mean pore diameter
    ax = axes[1]
    curves = []
    if S0 is not None:
        m, lo, hi = S0['d_mean']
        ax.fill_between(zx, lo, hi, color='0.5', alpha=0.25, lw=0, zorder=1)
        ax.plot(zx, m, '--', color='k', lw=2.0, alpha=0.9, zorder=2, label='reference (zero flux)')
        curves.append(np.where(R['interior'], m, np.nan))
    if S is not None:
        m, lo, hi = S['d_mean']
        ax.fill_between(zx, lo, hi, color=WONG['blue'], alpha=0.25, lw=0, zorder=3)
        ax.plot(zx, m, '-', color=WONG['blue'], lw=2.8, zorder=4, label='steady permeation')
        curves.append(np.where(P['interior'], m, np.nan))
    shade_gel(ax, R, P)
    mark_walls(ax, R, P)
    finish_axes(ax, r'$\langle D_{\rm pore}\rangle$  ($\sigma$)', '(b) mean pore diameter per bin')
    robust_ylim(ax, curves, pad=0.3, qlo=0, qhi=100, include_zero=True)
    smart_legend(ax, fontsize=11)
    txt = []
    if S0 is not None:
        txt.append(f"reference: {fmt_mu(S0['d_mean'][0][R['interior']])}")
    if S is not None:
        txt.append(f"steady: {fmt_mu(S['d_mean'][0][P['interior']])}")
    annotate_box(ax, r'$\langle D_{\rm pore}\rangle$ in gel interior ($\sigma$)' + '\n' + '\n'.join(txt), loc='lower left', fontsize=12)
    # (c) PSD of the interior
    ax = axes[2]
    if S0 is not None:
        Dc, m, lo, hi = S0['regions']['interior']
        ax.fill_between(Dc, lo, hi, color='0.5', alpha=0.25, lw=0, zorder=1)
        ax.plot(Dc, m, '--', color='k', lw=2.0, alpha=0.9, zorder=2, label=f"reference interior  $\\langle D\\rangle$ = {sig(S0['D_reg']['interior'][0])}")
    if S is not None:
        Dc, m, lo, hi = S['regions']['interior']
        ax.fill_between(Dc, lo, hi, color=WONG['blue'], alpha=0.25, lw=0, zorder=3)
        ax.plot(Dc, m, '-', color=WONG['blue'], lw=2.8, zorder=4, label=f"steady interior  $\\langle D\\rangle$ = {sig(S['D_reg']['interior'][0])}")
        for name, col in (('feed half', WONG['vermillion']), ('permeate half', WONG['green'])):
            Dc, m, lo, hi = S['regions'][name]
            ax.plot(Dc, m, '-', color=col, lw=1.8, alpha=0.9, zorder=4, label=f"steady, {name}  $\\langle D\\rangle$ = {sig(S['D_reg'][name][0])}")
    ax.set_xlabel(r'pore diameter $D$  ($\sigma$)')
    ax.set_ylabel('probability density')
    ax.set_title('(c) PSD of the membrane interior (volume-weighted)')
    ax.set_xlim(2.0 * cfg.PSD_R_PROBE, cfg.PSD_DMAX)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=11)
    return _save(fig, cfg, 'perm_psd')


# ===========================================================================
#  7. EXPANSE SYNC
# ===========================================================================
_PROD_DAT = ('sigmazz_polymer', 'sigmazz_solvent', 'sigmaxx_polymer', 'sigmaxx_solvent',
             'sigmayy_polymer', 'sigmayy_solvent', 'solvent_density_z', 'disp_z_polymer',
             'strain_zz', 'strain_piston', 'piston_position', 'support_position', 'piston_force', 'piston_force_avg',
             'box_dimensions', 'gel_dimensions_bb', 'gel_dimensions_rg', 'polymer_com', 'gel_edges')
# support_position: written since the 2026-09-18 symmetric drive (support driven up with the
# piston); older top-only runs have none and every reader falls back to the reference plane.
_REF_DAT = ('sigmazz_polymer_ref', 'sigmazz_solvent_ref', 'sigmaxx_polymer_ref', 'sigmaxx_solvent_ref',
            'sigmayy_polymer_ref', 'sigmayy_solvent_ref', 'solvent_density_z_ref')
# written only by runs since 2026-09-05: staged when a login happens anyway, but their
# absence never triggers one (older runs never wrote them)
_REF_DAT_OPT = ('piston_force_avg_ref',)
_REQUIRED_PROD = ('sigmazz_polymer', 'sigmazz_solvent', 'solvent_density_z', 'strain_zz',
                  'piston_force', 'box_dimensions', 'gel_dimensions_bb', 'disp_z_polymer')
# two-piston compression extras (per level) and permeation-mode lists (no level tag), 2026-09-16
_TWO_PIST_DAT = ('piston_pressure', 'permeation', 'pressure_reservoirs')
_PERM_DAT = ('sigmazz_polymer', 'sigmazz_solvent', 'sigmaxx_polymer', 'sigmaxx_solvent',
             'sigmayy_polymer', 'sigmayy_solvent', 'solvent_density_z', 'disp_z_polymer', 'strain_zz',
             'piston_position', 'piston_velocity', 'piston_force', 'piston_force_avg', 'piston_pressure',
             'permeation', 'permeate_count', 'pressure_feed', 'pressure_permeate',
             'box_dimensions', 'gel_dimensions_bb', 'gel_dimensions_rg', 'polymer_com', 'stress_aniso')
_PERM_REQUIRED = ('sigmazz_polymer', 'sigmazz_solvent', 'solvent_density_z', 'piston_position', 'permeation')


def _present(p):
    """a Path exists, or (pattern with '*' standing in for an unresolved <steps> tag)
    a file whose tag is all digits matches it -- so a permeation pattern never
    accepts a compression level's `<tag>_c<lvl>` file (2026-09-24)."""
    if '*' in p.name:
        rx = re.compile('^' + re.escape(p.name).replace(r'\*', '[0-9]+') + '$')
        return p.parent.exists() and any(rx.match(f.name) for f in p.parent.iterdir())
    return p.exists()


def sync_files(cfg, levels=None):
    """(data_files, traj_files, required) the notebooks read for `levels`
    (default: every level in cfg.COMP_LEVELS)."""
    levels = cfg.COMP_LEVELS if levels is None else [str(l) for l in levels]
    data = [cfg.path(n) for n in _REF_DAT + _REF_DAT_OPT]
    traj = [cfg.traj('traj_ref')]
    req = [cfg.path('sigmazz_polymer_ref'), cfg.path('sigmazz_solvent_ref'), cfg.path('solvent_density_z_ref')]
    if cfg.mode == 'permeation':            # one continuous run, no level tags
        data += [cfg.path(n) for n in _PERM_DAT]
        traj += [cfg.traj('traj_stress')]
        req += [cfg.path(n) for n in _PERM_REQUIRED]
        return data, traj, req
    for l in levels:
        data += [cfg.path(n, l) for n in _PROD_DAT + (_TWO_PIST_DAT if cfg.two_pist else ())]
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
    data_files, traj_files, required = sync_files(cfg, levels)
    optional = {cfg.path(n).name for n in _REF_DAT_OPT}
    optional |= {cfg.path(n, l).name for n in _TWO_PIST_DAT for l in ([None] + list(cfg.COMP_LEVELS))}
    sync_pull(cfg, data_files, traj_files, required, optional, force,
              refresh=lambda: sync_files(cfg, levels))


def sync_pull(cfg, data_files, traj_files, required, optional, force=False, refresh=None):
    """The Expanse login + staging + SFTP pull behind sync_from_expanse, for any
    file list (2026-09-23: shared with lib/shear.py).  `optional` = basenames whose
    absence never triggers a login; `refresh()` -> (data_files, traj_files, required)
    is called after the pull to re-resolve file tags."""
    import paramiko
    import getpass
    import json
    absent_f = cfg.DATA_DIR / '.sync_absent.json'
    absent = set()
    if absent_f.exists() and not force:
        try:
            absent = set(json.loads(absent_f.read_text()))
        except Exception:
            absent = set()
    missing_dat = [f for f in data_files if not _present(f)]
    missing_traj = [f for f in traj_files if not _present(f)]
    # a missing trajectory logs in too (2026-09-24: before, only .dat files did, so a
    # traj_stress the first sync mis-staged was never fetched)
    missing_new = [f for f in missing_dat + missing_traj if f.name not in absent and f.name not in optional]
    missing_all = missing_dat + missing_traj
    missing_req = [f for f in required if not _present(f)]
    if not force and not missing_new:
        n_opt = sum(1 for f in missing_dat if f.name in optional)
        known = len(missing_all) - len(missing_new) - n_opt
        print(f'All data files present locally' + (f' except {known} known absent on the cluster' if known else '')
              + (f' (+{n_opt} optional file(s) this run never wrote, e.g. piston_force_avg_ref)' if n_opt else '')
              + f' ({len(missing_all)} target file(s) missing in total) -- skipping Expanse login.  '
              f'FORCE_SYNC=True re-checks everything.')
        if missing_req:
            print('  WARNING: required files missing: ' + ', '.join(f.name for f in missing_req))
        return
    why = 'force=True' if (force and not missing_new) else f'{len(missing_new)} data file(s) missing'
    print(f'Syncing from Expanse ({why}); {len(missing_all)} of {len(data_files) + len(traj_files)} target files missing locally.')
    if missing_new:
        print('  missing: ' + ', '.join(f.name for f in missing_new[:6]) + (' ...' if len(missing_new) > 6 else ''))
    # an unresolved <steps> tag ('*') must match DIGITS ONLY on the cluster (extglob), so a
    # permeation request never stages a compression level's <tag>_c<lvl> file (2026-09-24)
    bn = lambda paths: ' '.join('"' + p.name.replace('*', '+([0-9])') + '"' for p in paths)
    stage = f'{cfg.RUNS_ROOT}/triaxial_stage'
    script = r"""
set -u
shopt -s extglob
RUNS="__RUNS__"; TRAJ="__TRAJ__"; STAGE="__STAGE__"
rm -rf "$STAGE"; mkdir -p "$STAGE/data" "$STAGE/traj"
# newest-first lists of every candidate file; each requested name may be a glob
# pattern (the <steps> tag is '*' until the level's files exist locally), and the
# NEWEST match wins.
mapfile -t ALL_D < <(find "$RUNS" -name '*.dat' -not -path '*/triaxial_stage/*' -size +0 -printf '%T@ %p\n' 2>/dev/null | sort -rn | cut -d' ' -f2-)
for B in __DATA_BN__; do
  for p in "${ALL_D[@]}"; do b=${p##*/}; if [[ "$b" == $B ]]; then cp -p "$p" "$STAGE/data/" 2>/dev/null || true; break; fi; done
done
mapfile -t ALL_T < <(find "$TRAJ" "$RUNS" -name '*.lammpstrj' -not -path '*/triaxial_stage/*' -printf '%T@ %p\n' 2>/dev/null | sort -rn | cut -d' ' -f2-)
for B in __TRAJ_BN__; do
  for p in "${ALL_T[@]}"; do b=${p##*/}; if [[ "$b" == $B ]]; then cp -p "$p" "$STAGE/traj/" 2>/dev/null || true; break; fi; done
done
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
    if not isinstance(cfg.NSTEPS, (int, dict)):
        cfg._tags = {}                                     # re-resolve the tags from the new files
    if refresh is not None:
        data_files, traj_files, required = refresh()
    still = [f.name for f in data_files + traj_files if not _present(f)]
    absent_f.write_text(json.dumps(sorted(set(still)), indent=1))
    if still:
        print(f'Sync complete; {len(still)} file(s) not found on the cluster (remembered in {absent_f.name}): '
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
    """Grey membrane band: the reference gel, or level L's membrane bounds.  A
    permeation dict carries `shade_hi` (the final feed-piston plane when
    cfg.GEL_SHADE_TO_PISTON, 2026-09-24), which replaces the polymer-stress edge."""
    lo, hi = (L['z_mem_lo'], L['z_mem_hi']) if L is not None else (R['z_gel_lo'], R['z_gel_hi'])
    if L is not None and L.get('shade_hi') is not None:
        hi = L['shade_hi']
    ax.axvspan(zn(R, lo), zn(R, hi), **GEL_SHADE)


def _level_support_z(R, L, t=None):
    """Support plane of level L (held, or at time t); the reference plane when L carries
    no support track (top-only runs, permeation)."""
    if L is None:
        return R['z_support']
    if t is not None and callable(L.get('support_z_at')):
        return L['support_z_at'](t)
    return L.get('z_supp', R['z_support'])


def mark_walls(ax, R, L=None, ts=None):
    """Support (solid) + piston (dash-dot): reference positions when L is None, else
    the level's HELD positions -- per evolution timestep when ts is given (final dark,
    earlier faded).  Both plates move under the symmetric drive (2026-09-18)."""
    if L is None or ts is None:
        walls = [(_level_support_z(R, L), R['z_piston'] if L is None else L['z_pist'], True)]
    else:
        ts = np.asarray(ts)
        walls = [(_level_support_z(R, L, t), L['piston_z_at'](t), i == len(ts) - 1) for i, t in enumerate(ts)]
    for sz, pz, last in walls:
        if np.isfinite(sz):
            ax.axvline(zn(R, sz), color=('k' if last else '0.7'), ls='-',
                       lw=(1.5 if last else 1.0), alpha=(0.85 if last else 0.35), zorder=(4 if last else 3))
        if np.isfinite(pz):
            ax.axvline(zn(R, pz), color=('0.15' if last else '0.7'), ls='-.',
                       lw=(1.6 if last else 1.0), alpha=(0.9 if last else 0.35), zorder=(4 if last else 3))


def mark_level_walls(ax, R, levels):
    """Sweep panels: reference support (black solid) and, per level in its colour, the
    held piston (dash-dot) and -- when it moved (symmetric drive) -- the held support (solid)."""
    if np.isfinite(R['z_support']):
        ax.axvline(zn(R, R['z_support']), color='k', ls='-', lw=1.5, alpha=0.85, zorder=4)
    for i, L in enumerate(levels):
        zs = _level_support_z(R, L)
        if np.isfinite(zs) and abs(zs - R['z_support']) > 1e-6:
            ax.axvline(zn(R, zs), color=level_color(i), ls='-', lw=1.2, alpha=0.6, zorder=4)
        ax.axvline(zn(R, L['z_pist']), color=level_color(i), ls='-.', lw=1.2, alpha=0.6, zorder=4)


def finish_axes(ax, ylabel, title):
    ax.set_xlabel(r'$z/L$')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.grid(alpha=0.3)


def _data_obstacles(ax, renderer):
    """Display-space sample points of everything drawn on ax (lines, fills,
    markers) plus the boxes of its text annotations and existing legend."""
    from matplotlib.collections import PathCollection

    def dense(P, step=4.0, cap=4000):
        """resample a display-space polyline every `step` px (so a curve crossing
        a box counts however few vertices it has), capped at `cap` points."""
        if len(P) < 2:
            return P
        seg = np.hypot(*np.diff(P, axis=0).T)
        n_new = np.minimum(np.maximum((seg / step).astype(int), 1), 50)
        if n_new.sum() > cap:
            n_new = np.maximum((n_new * cap / n_new.sum()).astype(int), 1)
        out = [P[:1]]
        for i, k in enumerate(n_new):
            t = np.linspace(0, 1, k + 1)[1:, None]
            out.append(P[i] + t * (P[i + 1] - P[i]))
        return np.concatenate(out)

    pts, boxes = [], []
    for ln in ax.lines:
        if ln.get_transform() is not ax.transData:      # axhline / axvline guides are not data
            continue
        xy = ln.get_xydata()
        xy = xy[np.isfinite(xy).all(axis=1)]
        if len(xy):
            pts.append(dense(ln.get_transform().transform(xy)))
    for c in ax.collections:
        try:
            if isinstance(c, PathCollection):           # scatter: marker positions
                off = np.asarray(c.get_offsets(), float)
                off = off[np.isfinite(off).all(axis=1)]
                if len(off):
                    pts.append(c.get_offset_transform().transform(off))
                continue
            for path in c.get_paths():                  # fills, error bars: outline
                v = path.vertices
                v = v[np.isfinite(v).all(axis=1)]
                if len(v):
                    pts.append(dense(c.get_transform().transform(v)))
        except Exception:
            pass
    for t in ax.texts:
        try:
            boxes.append(t.get_window_extent(renderer))
        except Exception:
            pass
    leg = ax.get_legend()
    if leg is not None:
        try:
            boxes.append(leg.get_window_extent(renderer))
        except Exception:
            pass
    P = np.concatenate(pts) if pts else np.zeros((0, 2))
    return P, boxes


def _overlap_score(box, P, boxes):
    """number of sampled data points inside `box` (+ a large penalty per unit
    area-fraction of overlap with existing text / legend boxes)."""
    x0, y0, x1, y1 = box
    n_in = int(((P[:, 0] >= x0) & (P[:, 0] <= x1) & (P[:, 1] >= y0) & (P[:, 1] <= y1)).sum()) if len(P) else 0
    area = max((x1 - x0) * (y1 - y0), 1e-9)
    ov = 0.0
    for b in boxes:
        w = min(x1, b.x1) - max(x0, b.x0)
        h = min(y1, b.y1) - max(y0, b.y0)
        if w > 0 and h > 0:
            ov += w * h / area
    return n_in + 1000.0 * ov


_LOC_ORDER = ('upper right', 'upper left', 'lower right', 'lower left', 'center right',
              'center left', 'upper center', 'lower center', 'center')


def _candidate_boxes(ax, w, h, pad, locs=_LOC_ORDER):
    ab = ax.bbox
    out = {}
    for loc in locs:
        if 'left' in loc:
            x0 = ab.x0 + pad
        elif 'right' in loc:
            x0 = ab.x1 - pad - w
        else:
            x0 = 0.5 * (ab.x0 + ab.x1) - 0.5 * w
        if 'lower' in loc:
            y0 = ab.y0 + pad
        elif 'upper' in loc:
            y0 = ab.y1 - pad - h
        else:
            y0 = 0.5 * (ab.y0 + ab.y1) - 0.5 * h
        out[loc] = (x0, y0, x0 + w, y0 + h)
    return out


def smart_legend(ax, *args, outside='auto', clear_tol=3, **kw):
    """smart_legend(ax, ...) placed where it covers no data: scores the nine standard
    positions against every plotted point / fill / annotation and takes the
    clearest; if none is clear the legend goes OUTSIDE the axes (right, or below
    when the axes carries a colorbar).  outside=False keeps it inside."""
    kw.pop('loc', None)
    leg = ax.legend(*args, loc='upper right', **kw)
    try:
        ax.figure.canvas.draw()          # final autoscaled limits + layout before measuring anything
        renderer = ax.figure.canvas.get_renderer()
        bb = leg.get_window_extent(renderer)
        fs = leg._fontsize if hasattr(leg, '_fontsize') else 12.0
        pad = leg.borderaxespad * fs * ax.figure.dpi / 72.0
        leg.remove()                     # so the probe legend is not its own obstacle
        P, boxes = _data_obstacles(ax, renderer)
        cands = _candidate_boxes(ax, bb.width, bb.height, pad)
        scores = {loc: _overlap_score(box, P, boxes) for loc, box in cands.items()}
        best = min(_LOC_ORDER, key=lambda l: scores[l])
    except Exception as e:
        print(f'  (smart_legend fell back to loc="best": {type(e).__name__}: {e})')
        if leg.axes is not None:
            leg.remove()
        return ax.legend(*args, loc='best', **kw)
    if scores[best] <= clear_tol or outside is False:
        return ax.legend(*args, loc=best, **kw)
    if getattr(ax, '_tri_has_colorbar', False):     # a colorbar sits to the right -> go below
        labels = [t.get_text() for t in leg.get_texts()]
        ncol = 2 if max((len(l) for l in labels), default=0) > 24 else min(max(1, len(labels)), 4)
        return ax.legend(*args, loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=ncol, **kw)
    return ax.legend(*args, loc='upper left', bbox_to_anchor=(1.02, 1.0), **kw)


def annotate_box(ax, text, loc='lower left', fontsize=14, color='k'):
    """Text box in an axes corner.  `loc` is the preferred corner; if it would
    cover data or another box and a clearer corner exists, that one is used."""
    corners = {'lower left': (0.02, 0.03, 'left', 'bottom'), 'lower right': (0.98, 0.03, 'right', 'bottom'),
               'upper left': (0.02, 0.97, 'left', 'top'), 'upper right': (0.98, 0.97, 'right', 'top')}
    x, y, ha, va = corners[loc]
    t = ax.text(x, y, text, transform=ax.transAxes, va=va, ha=ha, fontsize=fontsize, color=color,
                bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.85))
    try:
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
        bb = t.get_window_extent(renderer)
        P, boxes = _data_obstacles(ax, renderer)
        boxes = [b for b in boxes if not (abs(b.x0 - bb.x0) < 1 and abs(b.y0 - bb.y0) < 1)]
        pad = 0.02 * ax.bbox.width
        cands = _candidate_boxes(ax, bb.width, bb.height, pad, locs=tuple(corners))
        scores = {c: _overlap_score(cands[c], P, boxes) for c in corners}
        if scores[loc] > 3:
            best = min(corners, key=lambda c: (scores[c], c != loc))
            if scores[best] < scores[loc]:
                x, y, ha, va = corners[best]
                t.set_position((x, y))
                t.set_ha(ha)
                t.set_va(va)
    except Exception as e:
        print(f'  (annotate_box placement check skipped: {type(e).__name__}: {e})')
    return t


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
                   annotate=True, legend=True):
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
        ax._tri_has_colorbar = True
    if ref is not None and legend:
        smart_legend(ax, fontsize=12)
    if annotate:
        interior = L['interior']
        if flat_inside(stack[-1], z, interior, cfg.flat_tol):
            annotate_box(ax, 'final (plateau)\nmean in gel = ' + fmt_mu(stack[-1][interior]))
        else:
            ax.text(0.02, 0.03, 'final not flat inside\n(means omitted)', transform=ax.transAxes,
                    va='bottom', ha='left', fontsize=13, color='0.35')


def overlay_levels(ax, R, levels, get_z, get_ts, get_stack, cfg, ref=None, autoscale_mask=None,
                   ylabel='', title='', include_zero=True, pad=0.15, qlo=2):
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
    mark_level_walls(ax, R, levels)
    finish_axes(ax, ylabel, title)
    if finals:
        robust_ylim(ax, finals, zmask=autoscale_mask, pad=pad, qlo=qlo, qhi=100, include_zero=include_zero)
    return finals


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
        ax.annotate(f'plateau: Rg {sig(L.get("eps_rg", np.nan))}  BB {sig(L.get("eps_bb", np.nan))}',
                    (st[-1], eps_rg[-1]), textcoords='offset points', xytext=(-6, 9), ha='right',
                    va='bottom', fontsize=10, color=col,
                    bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.8))
    ax.set_xlabel('time step')
    ax.set_ylabel(r'compression strain  $\varepsilon=(L_0-L)/L_0$')
    ax.set_title('Strain diagnostic: solid $\\varepsilon_{Rg}$, dashed $\\varepsilon_{BB}$,\n'
                 'dotted = applied target, shaded = plateau window', fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=10)
    if not any_:
        plt.close(fig)
        print('strain diagnostic skipped (no strain_zz files)')
        return None
    # one level -> tag the file with it (keeps clear of the long-form notebook's untagged file)
    return _save(fig, cfg, stem, levels[0]['lvl'] if len(levels) == 1 else None)


def _phi_panel(ax, cfg, R, D, L=None, title='', bands=True, delta_from=None):
    """One phi_s profile panel.  delta_from=<state dict> (2026-09-24) replaces the
    in-gel means box by the in-gel CHANGE of each estimator relative to that state
    (steady - reference), two significant figures."""
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
        drawn.append((key, lab, m))
    ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6)
    ax.axhline(cfg.PHI_FLOOR, color='r', ls=':', lw=1.0, alpha=0.5)
    shade_gel(ax, R, L)
    mark_walls(ax, R, L)
    finish_axes(ax, r'$\phi_s$', title)
    ax.set_ylim(0, 1.15)
    smart_legend(ax, fontsize=12)
    mask = L['interior'] if L is not None else R['interior']
    short = lambda key: key.split('_')[1]
    if delta_from is not None:
        txt = '\n'.join(f'{short(key)}: {np.nanmean(m[mask]) - np.nanmean(delta_from[key][0][delta_from["interior"]]):+.2g}'
                        for key, lab, m in drawn if delta_from.get(key) is not None)
        if txt:
            annotate_box(ax, r'$\Delta\phi_s$ (steady $-$ reference), in gel' + '\n' + txt, loc='lower right', fontsize=12)
        return
    txt = '\n'.join(f'{short(key)}: {fmt_mu(m[mask])}' for key, lab, m in drawn)
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
        if kind == 't':      # normalise by the bath pressure (permeation: the permeate pressure P_perm, 2026-09-24)
            Pb = _p_norm(cfg, L)
            ev = ev / Pb
            ref = None if ref is None else tuple(np.asarray(r) / Pb for r in ref)
            ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6, zorder=1)
            lab = lab + (r'$/P_{\rm perm}$' if _is_perm(L) else r'$/P_{\rm bath}$')
        plot_evolution(ax, cfg, R, L, R['z'], ts, ev, lab + '$(z,t)$', title, ref=ref, band=band,
                       annotate=False, legend=False)
        if kind == 't':      # zoom on the band around P_bath (edge bin excluded), keep 1 well inside
            robust_ylim(ax, list(ev) + ([ref[0]] if ref is not None else []), zmask=_scale_mask(cfg, R, L),
                        pad=0.45, qlo=2, qhi=100, include_zero=False)
            lo, hi = ax.get_ylim()
            ax.set_ylim(min(lo, 1 - 0.3 * (hi - lo)), max(hi, 1 + 0.3 * (hi - lo)))
            note = ('steady' if _is_perm(L) else 'plateau') + ' mean in gel interior = ' + fmt_mu(S['t_plat'][L['interior']] / Pb)
        if kind == 'net':
            mask = (R['z'] >= L['z_mem_lo'] + cfg.wall_margin)
            robust_ylim(ax, list(ev), zmask=mask, pad=0.15)
            note = 'plateau mean in gel interior = ' + fmt_mu(S['net_plat'][L['interior']])
        h, lab = ax.get_legend_handles_labels()
        h.append(Patch(alpha=0, label=note))          # the number rides in the legend, never on data
        lab.append(note)
        smart_legend(ax, handles=h, labels=lab, fontsize=12)
    return _save(fig, cfg, stem, L['lvl'])


def _is_perm(L):
    return L is not None and L.get('mode') == 'permeation'


def _p_norm(cfg, L):
    """Pressure the total-stress panels are normalised by: the applied permeate pressure
    of a permeation run (P_perm = P_target; 2026-09-24), else cfg.P_BARO."""
    W = L.get('wet') if L is not None else None
    if _is_perm(L) and W is not None and 'P_perm_app' in W.get('plat', {}) and np.isfinite(W['plat']['P_perm_app']):
        return float(W['plat']['P_perm_app'])
    return float(cfg.P_BARO)


def _pending_panels(cfg, R, L, stem, suptitle, titles, ylabels, note):
    """Placeholder figure (permeation, 2026-09-24): the axes, membrane shading and plate
    planes are drawn but no data -- for quantities that need the pore pressure from the
    solvent chemical potential, which the deck does not measure yet."""
    fig, axes = plt.subplots(1, len(titles), figsize=(25 if len(titles) == 3 else 13, 6.5),
                             constrained_layout=True, squeeze=False)
    fig.suptitle(suptitle, fontsize=13, fontweight='bold')
    for ax, t, yl in zip(axes[0], titles, ylabels):
        shade_gel(ax, R, L)
        mark_walls(ax, R, L)
        finish_axes(ax, yl, t)
        ax.set_ylim(-1, 1)
        ax.text(0.5, 0.5, note, ha='center', va='center', transform=ax.transAxes, fontsize=14, color='0.35')
    return _save(fig, cfg, stem, L['lvl'])


_PORE_NOTE = ('pending: needs the pore pressure from the solvent\nchemical-potential profile '
              '($p_{\\rm pore}=\\mu_s/\\bar v_s$), which the\npermeation deck does not measure yet')


def fig_total_stress(cfg, R, L):
    if _is_perm(L):
        return _stress_evo_panels(cfg, R, L, 't', 'total_stress_evolution',
                                  f'Total stress / permeate pressure ($P_{{\\rm perm}}={sig(_p_norm(cfg, L))}$), reference $\\rightarrow$ '
                                  f'permeation drive (production from step {fmt_step(L["evol_from"])}; steady window from '
                                  f'{fmt_step(L["halt_ts"])})  |  {cfg.sim_name}')
    return _stress_evo_panels(cfg, R, L, 't', 'total_stress_evolution',
                              f'Total stress / bath pressure ($P_{{\\rm bath}}={sig(cfg.P_BARO)}$), reference $\\rightarrow$ compressed '
                              f'(hold from step {fmt_step(L["evol_from"])}; plateau from {fmt_step(L["halt_ts"])})  |  {cfg.sim_name}')


def _scale_mask(cfg, R, L=None, levels=None):
    """z-bins used to autoscale profile figures: from wall_margin inside the gel bottom
    up to the feed piston (two-piston runs; the bins beyond the wet pistons are vacuum,
    where sigma^t = 0 and sigma' = -p_pore would otherwise dominate the axis)."""
    z = R['z']
    lo = (min(L['z_mem_lo'] for L in levels) if levels else L['z_mem_lo']) + cfg.wall_margin
    m = z >= lo
    if np.isfinite(R.get('z_feed', np.nan)):
        m &= z + 0.5 * cfg.binWidth <= R['z_feed'] - cfg.res_wall_margin
    return m


def _final_state(S, plat, ci):
    """plateau-averaged network profile with its 95 % band (t-interval across the
    plateau snapshots; each snapshot already carries the baseline noise)."""
    return mean_ci(S['net'][plat], ci)


def _net_panel(ax, cfg, R, L, comp, ylabel, title, band_color=WONG['blue'], line_color=WONG['blue'], label='final (plateau)'):
    """Reference (dashed black, grey band) + final plateau state (solid, band) of one
    network-stress component.  Returns the final mean profile or None."""
    S = L['stress'].get(comp)
    Rs = R['stress'].get(comp)
    zx = zn(R, R['z'])
    if Rs is not None:
        ax.fill_between(zx, Rs['net_lo'], Rs['net_hi'], color='0.5', alpha=0.25, lw=0, zorder=1)
        ax.plot(zx, Rs['net_m'], '--', color='k', lw=2.0, alpha=0.9, zorder=2, label=r'reference ($\varepsilon=0$)')
    m = None
    if S is not None:
        m, lo, hi = _final_state(S, L['plat'], cfg.ci_level)
        ax.fill_between(zx, lo, hi, color=band_color, alpha=0.25, lw=0, zorder=3)
        ax.plot(zx, m, '-', color=line_color, lw=2.8, zorder=4, label=label + f' ({int(L["plat"].sum())} snapshots)')
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    shade_gel(ax, R, L)
    mark_walls(ax, R, L)
    finish_axes(ax, ylabel, title)
    return m


def fig_network_stress(cfg, R, L):
    """Network stress sigma'_ii = sigma^t_ii - p_pore (Terzaghi) for zz, xx, yy: the
    eps = 0 REFERENCE (dashed, grey band; ~0 since the gel starts at equilibrium
    swelling) and the FINAL plateau-averaged state (solid blue, 95 % band) only.
    The evolution curves were dropped on 2026-09-17: while the gel consolidates the
    pore pressure is NOT uniform (solvent is still diffusing through the network),
    so subtracting one reservoir value is only valid in the equilibrated state.
    PERMEATION (2026-09-24): drawn blank -- under flow p_pore(z) must come from the
    solvent chemical-potential profile (p_pore = mu_s / v_s), not from a reservoir baseline."""
    if _is_perm(L):
        return _pending_panels(cfg, R, L, 'network_stress_final',
                               f"Network stress $\\sigma'=\\sigma^t-p_{{\\rm pore}}(z)$ under permeation  |  {cfg.sim_name}",
                               [f'({"abc"[i]}) ' + r"network $\sigma'_{%s}$" % c for i, c in enumerate(COMPONENTS)],
                               [r"$\sigma'_{%s}(z)$" % c for c in COMPONENTS], _PORE_NOTE)
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    fig.suptitle(f"Network stress (Terzaghi $\\sigma'=\\sigma^t-p_{{\\rm pore}}$): reference and final equilibrated state  |  "
                 f"p_pore(final) = {fmt_val_unc(L['stress']['zz']['pore'][-1], L['stress']['zz']['pore_half'][-1])}  |  {cfg.sim_name}",
                 fontsize=13, fontweight='bold')
    for ax, comp in zip(axes, COMPONENTS):
        title = f'({"abc"[COMPONENTS.index(comp)]}) ' + r"network $\sigma'_{%s}$" % comp
        if comp not in L['stress']:
            ax.text(0.5, 0.5, f'sigma{comp} files\nnot found', ha='center', va='center', transform=ax.transAxes)
            finish_axes(ax, r"$\sigma'_{%s}$" % comp, title)
            continue
        m = _net_panel(ax, cfg, R, L, comp, r"$\sigma'_{%s}(z)$" % comp, title)
        robust_ylim(ax, [m] + ([R['stress'][comp]['net_m']] if comp in R['stress'] else []), zmask=_scale_mask(cfg, R, L), pad=0.2)
        h, lab = ax.get_legend_handles_labels()
        note = 'final mean in gel interior = ' + fmt_mu(m[L['interior']])
        h.append(Patch(alpha=0, label=note))
        lab.append(note)
        smart_legend(ax, handles=h, labels=lab, fontsize=12)
    return _save(fig, cfg, 'network_stress_final', L['lvl'])


def _tr3(D, key):
    """(xx + yy + zz)/3 of the per-snapshot stacks D[comp][key]; None if a component is missing."""
    if not all(c in D for c in COMPONENTS):
        return None
    return (D['xx'][key] + D['yy'][key] + D['zz'][key]) / 3.0


def fig_thermo_pressure(cfg, R, L):
    """Thermodynamic pressure P_th = -(1/3) tr(sigma^t) (sign convention of the
    profiles: positive under compression), time evolution over the hold (cividis,
    final bold) with the eps = 0 reference dashed.  Needs the xx and yy profiles."""
    Pt = _tr3(L['stress'], 't')
    if Pt is None:
        print('thermodynamic-pressure figure skipped (sigmaxx / sigmayy files missing)')
        return None
    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    ts, ev = post_halt(cfg, L, L['ts'], Pt)
    ref = None
    Rt = _tr3(R['stress'], 't')
    if Rt is not None:
        ref = mean_ci(Rt, cfg.ci_level)
    plot_evolution(ax, cfg, R, L, R['z'], ts, ev, r'$P_{th}(z,t)$  (LJ)',
                   r'Thermodynamic pressure $P_{th}=-\frac{1}{3}\,\mathrm{tr}(\mathbf{\sigma}^t)$: reference $\rightarrow$ '
                   + ('evolution' if _is_perm(L) else r'hold $\rightarrow$ plateau'),
                   ref=ref, annotate=False, legend=False)
    ax.axhline(cfg.P_BARO, color='k', ls=':', lw=1.2, alpha=0.6, zorder=1)
    robust_ylim(ax, list(ev) + ([ref[0]] if ref is not None else []), zmask=_scale_mask(cfg, R, L), pad=0.45, qlo=2, qhi=100, include_zero=False)
    lo, hi = ax.get_ylim()
    ax.set_ylim(min(lo, cfg.P_BARO - 0.3 * (hi - lo)), max(hi, cfg.P_BARO + 0.3 * (hi - lo)))
    h, lab = ax.get_legend_handles_labels()
    if _is_perm(L):
        note = f'$P_{{th}}$ in gel interior: $\\approx$ {sig(np.nanmean(np.nanmean(Pt[L["plat"]], axis=0)[L["interior"]]))}'
    else:
        note = f'plateau mean in gel interior = {fmt_mu(np.nanmean(Pt[L["plat"]], axis=0)[L["interior"]])}   (dotted: $P_{{\\rm bath}}={sig(cfg.P_BARO)}$)'
    h.append(Patch(alpha=0, label=note))
    lab.append(note)
    smart_legend(ax, handles=h, labels=lab, fontsize=12)
    lab_eps = f'($\\varepsilon={L["eps"]:.2f}$)' if np.isfinite(L['eps']) else '(permeation drive)'
    fig.suptitle(f'Thermodynamic pressure evolution {lab_eps}  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    return _save(fig, cfg, 'thermo_pressure_evolution', L['lvl'])


def fig_osmotic_pressure(cfg, R, L):
    """Osmotic pressure Pi = -(1/3) tr(sigma') (positive under compression), reference
    (dashed, ~0) and final plateau state (solid blue) with 95 % bands.  Like the
    network stress, only the equilibrated state is meaningful.  PERMEATION
    (2026-09-24): drawn blank until the pore pressure comes from the solvent
    chemical potential (see fig_network_stress)."""
    if _is_perm(L):
        return _pending_panels(cfg, R, L, 'osmotic_pressure_final',
                               f'Osmotic pressure $\\mathit{{\\Pi}}=-\\frac{{1}}{{3}}\\,\\mathrm{{tr}}(\\mathbf{{\\sigma}}\')$ under permeation  |  {cfg.sim_name}',
                               [r"osmotic pressure $\mathit{\Pi}=-\frac{1}{3}\,\mathrm{tr}(\mathbf{\sigma}')$"],
                               [r'$\mathit{\Pi}(z)$  (LJ)'], _PORE_NOTE)
    Pn = _tr3(L['stress'], 'net')
    if Pn is None:
        print('osmotic-pressure figure skipped (sigmaxx / sigmayy files missing)')
        return None
    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    zx = zn(R, R['z'])
    Rn = _tr3(R['stress'], 'net')
    if Rn is not None:
        m, lo, hi = mean_ci(Rn, cfg.ci_level)
        ax.fill_between(zx, lo, hi, color='0.5', alpha=0.25, lw=0, zorder=1)
        ax.plot(zx, m, '--', color='k', lw=2.0, alpha=0.9, zorder=2, label=r'reference ($\varepsilon=0$)')
    m, lo, hi = mean_ci(Pn[L['plat']], cfg.ci_level)
    ax.fill_between(zx, lo, hi, color=WONG['blue'], alpha=0.25, lw=0, zorder=3)
    ax.plot(zx, m, '-', color=WONG['blue'], lw=2.8, zorder=4, label=f'final (plateau, {int(L["plat"].sum())} snapshots)')
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    shade_gel(ax, R, L)
    mark_walls(ax, R, L)
    finish_axes(ax, r'$\mathit{\Pi}(z)$  (LJ)',
                r"Osmotic pressure $\mathit{\Pi}=-\frac{1}{3}\,\mathrm{tr}(\mathbf{\sigma}')$: reference and final equilibrated state")
    robust_ylim(ax, [m] + ([mean_ci(Rn, cfg.ci_level)[0]] if Rn is not None else []), zmask=_scale_mask(cfg, R, L), pad=0.2)
    h, lab = ax.get_legend_handles_labels()
    note = 'final mean in gel interior = ' + fmt_mu(m[L['interior']])
    h.append(Patch(alpha=0, label=note))
    lab.append(note)
    smart_legend(ax, handles=h, labels=lab, fontsize=12)
    lab_eps = f'($\\varepsilon={L["eps"]:.2f}$)' if np.isfinite(L['eps']) else '(steady permeation; feed baseline -- p_pore is not uniform under flow)'
    fig.suptitle(f'Osmotic pressure {lab_eps}  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    return _save(fig, cfg, 'osmotic_pressure_final', L['lvl'])


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
    if not _is_perm(L):
        handles.append(Line2D([0], [0], color='0.4', lw=1.2, alpha=0.6, label='hold (faint = early)'))
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    shade_gel(ax, R, L)
    mark_walls(ax, R, L)
    finish_axes(ax, r'$\sigma_{zz}(z,t)$  (LJ)',
                'Partial and total normal stresses: reference $\\rightarrow$ evolution' if _is_perm(L) else
                f'Partial and total $\\sigma_{{zz}}$: reference $\\rightarrow$ compressed ($\\varepsilon={L["eps"]:.2f}$)')
    smart_legend(ax, handles=handles, fontsize=12)
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
    fig.suptitle(f'Piston pressure history ($A = {sig(L["area"])}\\,\\sigma^2$, $F_z = P\\,A$)  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    for ax in (axL, axG):
        ax.axvline(peak, color=WONG['blue'], ls='--', lw=1.6, alpha=0.7, label=f'peak $P$ (step {fmt_step(peak)})')
        if 'PF' in L:
            ax.axvspan(L['PF']['step0'], float(steps[-1]), color=WONG['green'], alpha=0.10,
                       label='plateau window (auto)')
        ax.set_xlabel('time step')
        ax.grid(alpha=0.3)
    axL.plot(steps, P, '-', color=WONG['vermillion'], lw=1.0, alpha=0.30, label=r'$P=F_z/A$ (raw)')
    axL.plot(steps, P_roll, '-', color=WONG['vermillion'], lw=2.6, alpha=0.95, label=f'rolling mean ({cfg.roll_win})')
    if 'pfa_P' in L:
        axL.plot(L['pfa_step'], L['pfa_P'], '-', color='k', lw=1.8, alpha=0.8, label='LMP block-avg')
    axL.axhline(0, color='k', ls='--', lw=0.8, alpha=0.4)
    axL.set_ylabel(r'$P = F_z/A$  (LJ / $\sigma^2$)')
    axL.set_title('(a) piston pressure')
    smart_legend(axL, fontsize=12)
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
    smart_legend(axG, fontsize=12)
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
        ax.errorbar(L['ts'], G['ratio'], yerr=G['ratio_err'], ls=ls, color=col, lw=2.0, marker='o', ms=5,
                    capsize=4, elinewidth=1.2,
                    label=fr"$\sigma'_{{zz}}/\sigma'_{{{comp}}}$   final = {fmt_val_unc(G['ratio_final'], G['ratio_final_err'])}")
    ax.axvspan(L['halt_ts'], float(L['ts'][-1]), color=WONG['green'], alpha=0.10, label='plateau window')
    ax.axhline(1.0, color='k', ls=':', lw=1.0, alpha=0.6)
    robust_ylim(ax, [G['ratio'] for G in L['G'].values()], pad=0.35, qlo=0, qhi=100)   # error bars may run off
    ax.set_xlabel('time step')
    ax.set_ylabel(r"$\langle\sigma'_{zz}\rangle_{\rm int}/\langle\sigma'_{ii}\rangle_{\rm int}$")
    ax.set_title(r"Network stress anisotropy $\sigma'_{zz}/\sigma'_{ii} = M/(M-2G)$"
                 + ('  (relative to $\\varepsilon=0$)' if cfg.G_SUBTRACT_REF else ''), fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=13)
    return _save(fig, cfg, 'network_stress_ratio', L['lvl'])


def fig_M(cfg, R, L):
    """Longitudinal modulus, two panes.
    (a) diagnostic: the reported (increment) network and piston M, filled, next to the
        absolute stress / eps values, hollow, with the eps = 0 readings each estimator
        subtracts -- shows where the old ~0.002 offset between the estimators came from.
    (b) presentation: the two increment M values alone, with their CIs."""
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    ci = int(cfg.ci_level * 100)
    sub = bool(cfg.M_SUBTRACT_REF)
    incr = ' (increment from $\\varepsilon=0$)' if sub else ''
    fig.suptitle('Longitudinal modulus' + incr + '   |   ' + cfg.RUN_ID + '   |   $\\varepsilon = ' + str(L['lvl']) + '$',
                 fontsize=14, fontweight='bold')
    has_p = 'M_pist' in L
    pist_abs = has_p and (L['M_pist_ref'] != 'measured')

    def _pt(ax, x, key, marker, color, ms, label, mfc=None, alpha=1.0, lw=2.5, cap=8):
        ax.errorbar([x], [L[key]], yerr=[[L[key] - L[key + '_lo']], [L[key + '_hi'] - L[key]]],
                    fmt=marker, ms=ms, color=color, capsize=cap, lw=lw, alpha=alpha, label=label,
                    **({'mfc': mfc} if mfc else {}))

    # ---- (a) diagnostic pane ------------------------------------------------
    _pt(axA, 0, 'M_net', 'o', WONG['blue'], 13,
        f"network  $M = {sig(L['M_net'])}$\n{ci}% CI [{sig(L['M_net_lo'])}, {sig(L['M_net_hi'])}]")
    axA.axhline(L['M_net'], color=WONG['blue'], ls='--', lw=1.2, alpha=0.5)
    if sub:
        _pt(axA, 0.15, 'M_net_abs', 'o', WONG['blue'], 10,
            f"absolute $\\sigma'_{{zz}}/\\varepsilon = {sig(L['M_net_abs'])}$\n(ref $\\sigma'_{{zz}} = {L['M_net_ref']:+.4f}$ subtracted)",
            mfc='none', alpha=0.7, lw=1.5, cap=5)
    if has_p:
        _pt(axA, 1, 'M_pist', 's', WONG['vermillion'], 13,
            f"piston  $M = {sig(L['M_pist'])}$\n{ci}% CI [{sig(L['M_pist_lo'])}, {sig(L['M_pist_hi'])}]"
            + ('\n(absolute: no $P_{\\rm ref}$ file)' if pist_abs and sub else ''))
        axA.axhline(L['M_pist'], color=WONG['vermillion'], ls='--', lw=1.2, alpha=0.5)
        if sub and not pist_abs:
            _pt(axA, 1.15, 'M_pist_abs', 's', WONG['vermillion'], 10,
                f"absolute $P/\\varepsilon = {sig(L['M_pist_abs'])}$\n(ref $P = {L['P_ref']:+.4f}$ subtracted)",
                mfc='none', alpha=0.7, lw=1.5, cap=5)
    axA.set_xticks([0, 1])
    if sub:
        pist_lab = (r'piston' + '\n' + r'$(P-P_{\rm ref})/\varepsilon$' if not pist_abs
                    else r'piston' + '\n' + r'$P/\varepsilon$  (no $P_{\rm ref}$ file)')
        axA.set_xticklabels([r'network' + '\n' + r"$(\langle\sigma'_{zz}\rangle_{\rm int}-\sigma'_{zz,\rm ref})/\varepsilon$",
                             pist_lab], fontsize=13)
    else:
        axA.set_xticklabels([r"network ($\langle\sigma'_{zz}\rangle_{\rm int}/\varepsilon$)", r'piston ($P/\varepsilon$)'], fontsize=14)
    axA.set_ylabel(r'$M$  (LJ units)')
    axA.set_title('(a) with the absolute values (hollow) they replace', fontsize=13)
    axA.set_xlim(-0.5, 1.6)
    axA.grid(axis='y', alpha=0.3)
    smart_legend(axA, fontsize=11)

    # ---- (b) presentation pane: the two M values only -------------------------
    _pt(axB, 0, 'M_net', 'o', WONG['blue'], 14,
        f"network  $M = {sig(L['M_net'])}$\n{ci}% CI [{sig(L['M_net_lo'])}, {sig(L['M_net_hi'])}]")
    axB.axhline(L['M_net'], color=WONG['blue'], ls='--', lw=1.2, alpha=0.5)
    if has_p:
        _pt(axB, 1, 'M_pist', 's', WONG['vermillion'], 14,
            f"piston  $M = {sig(L['M_pist'])}$\n{ci}% CI [{sig(L['M_pist_lo'])}, {sig(L['M_pist_hi'])}]")
        axB.axhline(L['M_pist'], color=WONG['vermillion'], ls='--', lw=1.2, alpha=0.5)
    axB.set_xticks([0, 1])
    axB.set_xticklabels([r'network' + '\n' + r"$\Delta\langle\sigma'_{zz}\rangle_{\rm int}/\varepsilon$",
                         r'piston' + '\n' + r'$\Delta P/\varepsilon$'] if sub else
                        [r"network ($\langle\sigma'_{zz}\rangle_{\rm int}/\varepsilon$)", r'piston ($P/\varepsilon$)'], fontsize=16)
    axB.set_ylabel(r'$M$  (LJ units)')
    axB.set_title('(b) longitudinal modulus, two independent estimates', fontsize=13)
    axB.set_xlim(-0.5, 1.5)
    axB.grid(axis='y', alpha=0.3)
    smart_legend(axB, fontsize=13)
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
                    capsize=8, lw=2.5, label=f"from {comp}:  $G = {sig(G['G'])}$\n{ci}% CI [{sig(G['lo'])}, {sig(G['hi'])}]")
        xs.append(k)
        vals.append(G['G'])
    if vals:
        ax.axhline(np.mean(vals), color='0.3', ls='--', lw=1.2, alpha=0.6,
                   label=f'mean of the two = {sig(np.mean(vals))}')
    ax.set_xticks([0, 1])
    ax.set_xticklabels([r"$G_x=(\sigma'_{zz}-\sigma'_{xx})/2\varepsilon$", r"$G_y=(\sigma'_{zz}-\sigma'_{yy})/2\varepsilon$"],
                       fontsize=14)
    ax.set_ylabel(r'$G$  (LJ units)')
    ax.set_title(f'Shear modulus from network-stress anisotropy\n{cfg.RUN_ID}  |  $\\varepsilon = {L["lvl"]}$',
                 fontsize=15)
    ax.set_xlim(-0.5, 1.5)
    ax.grid(axis='y', alpha=0.3)
    smart_legend(ax, fontsize=11)
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
        ax.set(xlabel=r'$\zeta=(z-z_\mathrm{perm})/L$', ylabel=r'$u_z/L$', xlim=(0, 1))
        ax.grid(alpha=0.3)
        smart_legend(ax, fontsize=11)
    axl.set_title(r'$u_z(\zeta,t)/L$ -- hold snapshots (data + fitted IC offset)', fontsize=15)
    lbl = '' if cfg.DC_FREE_AMPS else rf"$\beta=p_0/M={sig(F['beta'])}$, "
    axr.set_title(rf"Consolidation fit: $D_c={sig(F['Dc'])}\ \sigma^2/\tau$, " + lbl + rf"$R^2={sig(F['R2'])}$", fontsize=15)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=[axl, axr], fraction=0.015, pad=0.04).set_label('timestep')
    fig.suptitle(f'Cooperative diffusivity fit (level _c{L["lvl"]})  |  {cfg.sim_name}  |  '
                 f'plate-gap closure $\\Delta L={sig(F["DL"])}\\,\\sigma$: support {fs:.0%} / piston {1 - fs:.0%}'
                 + ('  (symmetric drive)' if abs(fs - 0.5) < 0.05 else '  (top-only drive)' if fs < 0.05 else ''),
                 fontsize=12, fontweight='bold')
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
                    capsize=8, lw=2.5, label=f"{lab} = {sig(v['k'])}")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['network $M$', 'piston $M$'], fontsize=16)
    ax.set_ylabel(r'$\kappa = D_c/M$  (LJ: $\sigma^5/(\epsilon\,\tau)$)')
    ax.set_title(f'Hydraulic permeability $\\kappa = D_c/M$\n$D_c = {sig(L["Dc"]["Dc"])}$  |  '
                 f'$\\varepsilon = {L["lvl"]}$', fontsize=15)
    ax.set_xlim(-0.5, 1.5)
    ax.grid(axis='y', alpha=0.3)
    smart_legend(ax, fontsize=13)
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
            zs = _level_support_z(R, L)              # held support (moved: symmetric drive)
            if np.isfinite(zs) and abs(zs - R['z_support']) > 1e-6:
                ax.axvline(zn(R, zs), color=level_color(i), ls='-', lw=1.2, alpha=0.6)
            txt.append(f'{L["eps"]:.2f}: {fmt_mu(m[L["interior"]])}')
        ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6)
        if np.isfinite(R['z_support']):
            ax.axvline(zn(R, R['z_support']), color='k', lw=1.5, alpha=0.85)
        finish_axes(ax, r'$\phi_s$', title)
        ax.set_ylim(0, 1.15)
        smart_legend(ax, handles=level_handles(levels, ref=True), fontsize=11)
        if txt:
            annotate_box(ax, 'in-gel means\n' + '\n'.join(txt), loc='lower right', fontsize=11)
        else:
            ax.text(0.5, 0.5, 'unavailable', ha='center', va='center', transform=ax.transAxes)
    return _save(fig, cfg, 'sweep_volfrac_profiles')


def _sweep_stress_panels(cfg, R, levels, kind, stem, suptitle):
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    fig.suptitle(suptitle, fontsize=13, fontweight='bold')
    mask = _scale_mask(cfg, R, levels=levels)
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
        scale = cfg.P_BARO if kind == 't' else 1.0
        if kind == 't':
            lab = lab + r'$/P_{\rm bath}$'
            ref = None if ref is None else tuple(np.asarray(r) / scale for r in ref)
            ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6, zorder=1)
        overlay_levels(ax, R, levels, lambda L: L['z'], lambda L: L['ts'],
                       lambda L, c=comp, sc=scale: (L['stress'][c][kind] / sc if c in L['stress'] else None), cfg,
                       ref=ref, autoscale_mask=mask, ylabel=lab, title=title,
                       include_zero=(kind != 't'), pad=(0.45 if kind == 't' else 0.15))
        if kind == 't':
            lo, hi = ax.get_ylim()
            ax.set_ylim(min(lo, 1 - 0.3 * (hi - lo)), max(hi, 1 + 0.3 * (hi - lo)))
        smart_legend(ax, handles=level_handles(levels, ref=Rs is not None), fontsize=11)
    return _save(fig, cfg, stem)


def fig_total_stress_sweep(cfg, R, levels):
    return _sweep_stress_panels(cfg, R, levels, 't', 'sweep_total_stress_evolution',
                                f'Total stress / bath pressure ($P_{{\\rm bath}}={sig(cfg.P_BARO)}$), all levels '
                                f'(faint = early hold, bold = plateau)  |  {cfg.sim_name}')


def _matter_mask(cfg, D, comp='zz', thr=0.05):
    """Bins that hold matter in the plateau state of D (a level or R): |sigma^t_comp| > thr.
    Beyond the wet pistons and below the support the box is vacuum, where sigma^t = 0 and
    the Terzaghi sigma' = -p_pore would draw as a -P_bath shelf; those bins are blanked in
    the network-stress / osmotic-pressure overlays (2026-09-26)."""
    S = D['stress'].get(comp) if isinstance(D.get('stress'), dict) else None
    if S is None:
        return None
    t = S.get('t_plat')
    if t is None:
        t = S.get('t_m')
    if t is None:
        return None
    return np.abs(np.asarray(t, float)) > thr


def _final_overlay(ax, cfg, R, levels, get_stack, ref_stack, ylabel, title, matter_only=False):
    """Sweep overlay of FINAL plateau states (colour = level, 95 % bands) + the
    reference (dashed black).  get_stack(L) -> per-snapshot stack or None.
    matter_only=True blanks the vacuum bins (see _matter_mask)."""
    zx = zn(R, R['z'])
    finals = []

    def _blank(m, lo, hi, D):
        if not matter_only:
            return m, lo, hi
        mk = _matter_mask(cfg, D)
        if mk is None:
            return m, lo, hi
        return tuple(np.where(mk, np.asarray(v, float), np.nan) for v in (m, lo, hi))

    if ref_stack is not None:
        m, lo, hi = _blank(*mean_ci(ref_stack, cfg.ci_level), R)
        ax.fill_between(zx, lo, hi, color='0.5', alpha=0.2, lw=0, zorder=1)
        ax.plot(zx, m, '--', color='k', lw=2.0, alpha=0.9, zorder=2)
        finals.append(m)
    for i, L in enumerate(levels):
        st = get_stack(L)
        if st is None:
            continue
        m, lo, hi = _blank(*mean_ci(st[L['plat']], cfg.ci_level), L)
        ax.fill_between(zx, lo, hi, color=level_color(i), alpha=0.18, lw=0, zorder=3)
        ax.plot(zx, m, '-', color=level_color(i), lw=2.6, zorder=4)
        finals.append(m)
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.4)
    mark_level_walls(ax, R, levels)
    finish_axes(ax, ylabel, title)
    if finals:
        robust_ylim(ax, finals, zmask=_scale_mask(cfg, R, levels=levels), pad=0.2)
    smart_legend(ax, handles=level_handles(levels, ref=ref_stack is not None), fontsize=11)


def fig_network_stress_sweep(cfg, R, levels):
    """Network stress sigma'_ii, FINAL plateau state of every level (+ reference dashed).
    No evolution curves (2026-09-17): p_pore is only uniform once consolidation is over."""
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    fig.suptitle(f'Network stress (Terzaghi), final equilibrated state of every level  |  {cfg.sim_name}\n'
                 r"(drawn where there is matter: beyond the wet pistons and below the support $\sigma^t=0$, so $\sigma'=-p_{\rm pore}$ there is blanked)",
                 fontsize=13, fontweight='bold')
    for ax, comp in zip(axes, COMPONENTS):
        title = f'({"abc"[COMPONENTS.index(comp)]}) ' + r"network $\sigma'_{%s}$" % comp
        if not any(comp in L['stress'] for L in levels):
            ax.text(0.5, 0.5, f'sigma{comp} files\nnot found', ha='center', va='center', transform=ax.transAxes)
            finish_axes(ax, r"$\sigma'_{%s}$" % comp, title)
            continue
        Rs = R['stress'].get(comp)
        _final_overlay(ax, cfg, R, levels, lambda L, c=comp: (L['stress'][c]['net'] if c in L['stress'] else None),
                       Rs['net'] if Rs is not None else None, r"$\sigma'_{%s}(z)$" % comp, title, matter_only=True)
    return _save(fig, cfg, 'sweep_network_stress_final')


def fig_thermo_pressure_sweep(cfg, R, levels):
    """P_th = -(1/3) tr(sigma^t) evolution, all levels overlaid (faint -> bold), reference dashed."""
    if not any(_tr3(L['stress'], 't') is not None for L in levels):
        print('thermodynamic-pressure sweep figure skipped (sigmaxx / sigmayy files missing)')
        return None
    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    Rt = _tr3(R['stress'], 't')
    ref = mean_ci(Rt, cfg.ci_level) if Rt is not None else None
    ax.axhline(cfg.P_BARO, color='k', ls=':', lw=1.2, alpha=0.6, zorder=1)
    overlay_levels(ax, R, levels, lambda L: L['z'], lambda L: L['ts'], lambda L: _tr3(L['stress'], 't'), cfg, ref=ref,
                   autoscale_mask=_scale_mask(cfg, R, levels=levels), ylabel=r'$P_{th}(z,t)$  (LJ)',
                   title=r'Thermodynamic pressure $P_{th}=-\frac{1}{3}\,\mathrm{tr}(\mathbf{\sigma}^t)$, all levels (faint = early hold, bold = plateau)',
                   include_zero=False, pad=0.45)
    lo, hi = ax.get_ylim()
    ax.set_ylim(min(lo, cfg.P_BARO - 0.3 * (hi - lo)), max(hi, cfg.P_BARO + 0.3 * (hi - lo)))
    smart_legend(ax, handles=level_handles(levels, ref=ref is not None), fontsize=11)
    return _save(fig, cfg, 'sweep_thermo_pressure_evolution')


def fig_osmotic_pressure_sweep(cfg, R, levels):
    """Pi = -(1/3) tr(sigma'), final plateau state of every level (+ reference dashed)."""
    if not any(_tr3(L['stress'], 'net') is not None for L in levels):
        print('osmotic-pressure sweep figure skipped (sigmaxx / sigmayy files missing)')
        return None
    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    _final_overlay(ax, cfg, R, levels, lambda L: _tr3(L['stress'], 'net'), _tr3(R['stress'], 'net'),
                   r'$\mathit{\Pi}(z)$  (LJ)',
                   r"Osmotic pressure $\mathit{\Pi}=-\frac{1}{3}\,\mathrm{tr}(\mathbf{\sigma}')$, final equilibrated state of every level",
                   matter_only=True)
    return _save(fig, cfg, 'sweep_osmotic_pressure_final')


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
        smart_legend(ax, handles=level_handles(levels, ref=True), fontsize=11)
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
    axP.set_xlabel('time step')
    axP.set_ylabel(r'$P = F_z/A$  (LJ / $\sigma^2$)')
    axP.set_title('(a) piston pressure (rolling mean; dotted = plateau)')
    axG.set_xlabel('time step')
    axG.set_ylabel(r'$\ln P$')
    axG.set_title('(b) log piston pressure (relaxation view)')
    for ax in (axP, axG):
        ax.grid(alpha=0.3)
        smart_legend(ax, handles=level_handles(levels), fontsize=11)
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
            ax.errorbar(L['ts'], G['ratio'], yerr=G['ratio_err'], ls=ls, color=level_color(i), lw=2.0,
                        marker='o', ms=4, capsize=3, elinewidth=1.0)
    ax.axhline(1.0, color='k', ls=':', lw=1.0, alpha=0.6)
    robust_ylim(ax, [G['ratio'] for L in levels for G in L['G'].values()], pad=0.35, qlo=0, qhi=100)
    h = level_handles(levels) + [Line2D([0], [0], color='0.3', ls='-', lw=2, label=r"$\sigma'_{zz}/\sigma'_{xx}$"),
                                 Line2D([0], [0], color='0.3', ls='--', lw=2, label=r"$\sigma'_{zz}/\sigma'_{yy}$")]
    ax.set_xlabel('time step')
    ax.set_ylabel(r"$\langle\sigma'_{zz}\rangle_{\rm int}/\langle\sigma'_{ii}\rangle_{\rm int}$")
    ax.set_title(r"Network stress anisotropy $\sigma'_{zz}/\sigma'_{ii}=M/(M-2G)$, all levels", fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, handles=h, fontsize=11)
    return _save(fig, cfg, 'sweep_network_stress_ratio')


def fig_M_sweep(cfg, R, levels):
    """Longitudinal modulus vs applied strain, two panes (same layout as fig_M):
    (a) diagnostic: the reported (increment) network and piston M per level, filled,
        next to the absolute stress / eps values, hollow;
    (b) presentation: the two increment M estimates per level alone, with CIs.
    The stress-strain curve and its least-squares slopes live in fig_stress_strain_sweep."""
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    sub = bool(cfg.M_SUBTRACT_REF)
    incr = ' (increment from $\\varepsilon=0$)' if sub else ''
    fig.suptitle('Longitudinal modulus across the sweep' + incr + '   |   ' + cfg.sim_name, fontsize=13, fontweight='bold')
    eps = np.array([L['eps'] for L in levels])
    hp = [L for L in levels if 'M_pist' in L]
    ep = np.array([L['eps'] for L in hp])
    n_abs = sum(L['M_pist_ref'] != 'measured' for L in hp)
    net_lab = (r"network  $(\langle\sigma'_{zz}\rangle_\mathrm{int}-\sigma'_{zz,\rm ref})/\varepsilon$" if sub
               else r"network  $\langle\sigma'_{zz}\rangle_\mathrm{int}/\varepsilon$")
    pist_lab = ((r'piston  $(P-P_{\rm ref})/\varepsilon$' if (sub and n_abs == 0) else r'piston  $P/\varepsilon$')
                + ('  (absolute: no $P_{\\rm ref}$ file)' if (sub and n_abs) else ''))

    def _ser(ax, e, key, Ls, marker, color, label, hollow=False, dx=0.0):
        v = np.array([L[key] for L in Ls])
        lo = np.array([L[key + '_lo'] for L in Ls]); hi = np.array([L[key + '_hi'] for L in Ls])
        kw = dict(mfc='none', alpha=0.7, lw=1.5, ms=8, capsize=4, ls=':') if hollow else dict(lw=2, ms=10, capsize=6, ls='-')
        ax.errorbar(e + dx, v, yerr=[v - lo, hi - v], fmt=marker, color=color, label=label, **kw)

    # ---- (a) diagnostic ---------------------------------------------------------
    _ser(axA, eps, 'M_net', levels, 'o', WONG['blue'], net_lab)
    if sub:
        _ser(axA, eps, 'M_net_abs', levels, 'o', WONG['blue'], r"network, absolute $\langle\sigma'_{zz}\rangle_\mathrm{int}/\varepsilon$",
             hollow=True, dx=0.003)
    if hp:
        _ser(axA, ep, 'M_pist', hp, 's', WONG['vermillion'], pist_lab)
        if sub and n_abs == 0:
            _ser(axA, ep, 'M_pist_abs', hp, 's', WONG['vermillion'], r'piston, absolute $P/\varepsilon$', hollow=True, dx=-0.003)
    if sub:
        Rzz = R['stress']['zz']
        txt = f"$\\varepsilon=0$ readings subtracted:\n  $\\sigma'_{{zz,\\rm ref}} = {float(Rzz['net_interior']):+.4f}$"
        if hp and n_abs == 0:
            txt += f"\n  $P_{{\\rm ref}} = {float(R['P_ref']):+.4f}$"
        annotate_box(axA, txt, loc='lower right', fontsize=11)
    axA.set_xlabel(r'applied strain  $\varepsilon$')
    axA.set_ylabel(r'$M$  (LJ units)')
    axA.set_title('(a) with the absolute values (hollow) they replace', fontsize=13)
    axA.grid(alpha=0.3)
    smart_legend(axA, fontsize=11)

    # ---- (b) presentation ----------------------------------------------------
    _ser(axB, eps, 'M_net', levels, 'o', WONG['blue'],
         r"network  $\Delta\langle\sigma'_{zz}\rangle_\mathrm{int}/\varepsilon$" if sub else net_lab)
    if hp:
        _ser(axB, ep, 'M_pist', hp, 's', WONG['vermillion'],
             (r'piston  $\Delta P/\varepsilon$' if (sub and n_abs == 0) else pist_lab))
    axB.set_xlabel(r'applied strain  $\varepsilon$')
    axB.set_ylabel(r'$M$  (LJ units)')
    axB.set_title('(b) longitudinal modulus per level, two independent estimates', fontsize=13)
    axB.grid(alpha=0.3)
    smart_legend(axB, fontsize=13)
    return _save(fig, cfg, 'sweep_modulus')


def fig_stress_strain_sweep(cfg, R, levels):
    """The stress-strain curve itself: plateau piston stress and interior network stress vs
    applied strain, INCLUDING each estimator's eps = 0 reading (hollow), with a least-squares
    line through each series.  The slopes are a separate estimate of M that uses only
    differences BETWEEN levels, so any constant offset an estimator carries (piston preload,
    profile bias) cancels whether or not M_SUBTRACT_REF is on; curvature of the points is the
    stress-strain nonlinearity.  (Was panel (b) of fig_M_sweep until 2026-09-14.)"""
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    eps = np.array([L['eps'] for L in levels])
    hp = [L for L in levels if 'M_pist' in L]
    Rzz = R['stress']['zz']
    fits = []

    def _series(e, s, err, color, marker, name, e0=None, s0=None, err0=None):
        ax.errorbar(e, s, yerr=err, fmt=marker + '-', lw=2, ms=9, color=color, capsize=6, label=name)
        ee, ss = list(e), list(s)
        if e0 is not None and s0 is not None and np.isfinite(s0):
            ax.errorbar([e0], [s0], yerr=err0, fmt=marker, ms=9, mfc='none', color=color, capsize=6)
            ee, ss = [e0] + ee, [s0] + ss
        if len(ee) >= 2:
            slope, icpt = np.polyfit(ee, ss, 1)
            xs = np.linspace(0, max(ee) * 1.05, 20)
            ax.plot(xs, slope * xs + icpt, ls=':', lw=1.5, color=color, alpha=0.8)
            fits.append((name.split()[0], slope, len(ee)))

    sn = np.array([L['M_net_abs'] * L['eps'] for L in levels])
    sn_err = [np.abs(np.array([L['M_net_abs_lo'] * L['eps'] for L in levels]) - sn),
              np.abs(np.array([L['M_net_abs_hi'] * L['eps'] for L in levels]) - sn)]
    _series(eps, sn, sn_err, WONG['blue'], 'o', "network $\\langle\\sigma'_{zz}\\rangle_{\\rm int}$ (plateau)",
            e0=0.0, s0=float(Rzz['net_interior']), err0=float(Rzz.get('net_interior_half', 0.0)))
    if hp:
        ep = np.array([L['eps'] for L in hp])
        Pp = np.array([L['P_final'] for L in hp])
        pref = R.get('P_ref', np.nan)
        _series(ep, Pp, [Pp - [L['PF']['lo'] for L in hp], [L['PF']['hi'] for L in hp] - Pp],
                WONG['vermillion'], 's', 'piston $P=\\langle F_z\\rangle/A$ (plateau)',
                e0=0.0, s0=pref, err0=(float(R['P_ref_hi'] - R['P_ref_lo']) / 2 if np.isfinite(pref) else None))
    ax.plot([], [], 'o', mfc='none', color='0.4', label=r'hollow = $\varepsilon=0$ reading')
    ax.set_title('Stress vs strain, least-squares slopes incl. $\\varepsilon=0$:  ' +
                 ',  '.join(f"$M_{{\\rm {n[:4]}}}\\approx{sig(s)}$ ({k} pts)" for n, s, k in fits) + '\n' + cfg.sim_name,
                 fontsize=12)
    ax.set_xlabel(r'applied strain  $\varepsilon$')
    ax.set_ylabel(r'plateau stress  (LJ)')
    ax.grid(alpha=0.3)
    ax.set_xlim(left=-0.01)
    smart_legend(ax, fontsize=12)
    return _save(fig, cfg, 'sweep_stress_strain')


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
    ax.set_xlabel(r'applied strain  $\varepsilon$')
    ax.set_ylabel(r'$G$  (LJ units)')
    ax.set_title('Shear modulus from network-stress anisotropy vs strain', fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=11)
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
        ax.annotate(f"$R^2$={sig(L['Dc']['R2'])}", (L['eps'], L['Dc']['Dc']), textcoords='offset points',
                    xytext=(0, -12), ha='center', va='top', fontsize=10, color='0.35',
                    annotation_clip=True)
    ax.margins(x=0.12, y=0.15)
    ax.axhline(cfg.DC_SLOW_REF, color='0.4', ls=':', lw=1.5, label=f'slow reference $D_c$ = {sig(cfg.DC_SLOW_REF)}')
    ax.set_xlabel(r'applied strain  $\varepsilon$')
    ax.set_ylabel(r'$D_c$  ($\sigma^2/\tau$)')
    ax.set_title(r'(a) $D_c$ vs applied strain', fontsize=15)
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
        for i in F['early']:
            c = cmap(norm(F['ts'][i]))
            ax.plot(F['zf'], u0_b + F['uhat'][i][F['idx']], 'o', color=c, ms=2.5, alpha=0.35)
            ax.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=1.5)
        dlL, fs = F['DL'] / F['L'], F['f_sup']
        ax.plot([0, 1], [fs * dlL, (fs - 1.0) * dlL], 'k:', lw=1.6)
        ax.set(xlabel=r'$\zeta$', ylabel=r'$u_z/L$', xlim=(0, 1))
        ax.set_title(fr"({'bcdefgh'[k]}) $\varepsilon={L['lvl']}$:  $D_c={sig(F['Dc'])}$, $R^2={sig(F['R2'])}$",
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
    smart_legend(ax, fontsize=12)
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
            mixed = (L.get('M_pist_ref') == 'absent')
            mp, mn = (L['M_pist_abs'], L['M_net_abs']) if mixed else (L['M_pist'], L['M_net'])
            print(f"     observed M_piston/M_network - 1 = {mp / mn - 1:+.1%}   vs   "
                  f"predicted unrelaxed excess (slow D_c) = {F['hold_check']['slow']['avg']:+.1%}"
                  + ('   [both ABSOLUTE here: no piston_force_avg_ref, so the increment ratio is not like-for-like]' if mixed else ''))
    print(f'  (triaxial_compression.lmp sizes each hold as n_tau_hold * tau_1 from the live compressed '
          f'BB thickness and Dc_est = {cfg.DC_SLOW_REF:.2f}; tau_1 ~ (1-eps)^2.)')


def print_summary(cfg, levels):
    """One line per level with the headline numbers, then the hold check."""
    ci = int(cfg.ci_level * 100)
    how = 'M = increment from the eps = 0 reference' if cfg.M_SUBTRACT_REF else 'M = absolute stress / eps'
    if cfg.M_SUBTRACT_REF and any(L.get('M_pist_ref') == 'absent' for L in levels):
        how += '; piston M ABSOLUTE (no piston_force_avg_ref file)'
    print(f'\nSUMMARY  ({cfg.sim_name}; {ci}% CIs; {how})')
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


# ===========================================================================
#  12. TWO-PISTON (feed / permeate NPT-piston) RUNS  -- 2026-09-16
# ===========================================================================
# Files written by triaxial_{compression,permeation}_two_pist.lmp carry ONE column
# set per piston with a '# ...' header naming them ([load |] feed | perm).  The
# compression loaders above keep reading the first value column (= the load
# piston, the network load), so every one-piston figure works unchanged; the
# extras below add the wet-piston bath check and the permeation analysis.
def load_wet_pistons(cfg, R, lvl=None, plat_from=None):
    """piston_pressure[_c<lvl>] (+ permeation[_c<lvl>] for compression: solvent expelled)
    -> dict(step, P_*_meas/app series, plat means, dV_*) or None when the run is one-piston.
    The MEASURED pressures come from piston_force_avg (fix ave/time block means over
    each volume_freq window, 2026-09-23) when that file exists, divided by l_x l_y;
    the instantaneous fix-print samples of piston_pressure are the fallback.  The
    applied pressures are interpolated onto the block-mean steps."""
    names, tab = read_piston_table(cfg.path('piston_pressure', lvl))
    if tab.size == 0 or names is None:
        return None
    W = dict(step=tab[:, 0], names=names, source='piston_pressure (instantaneous samples)')
    for nm in names[1:]:
        col = table_col(names, tab, nm)
        if col is not None:
            W[nm] = col
    fn, ft = read_piston_table(cfg.path('piston_force_avg', lvl))
    if ft.size and fn is not None and np.isfinite(R.get('AREA', np.nan)) and R['AREA'] > 0:
        fmap = {'F_load': ('P_load_meas', 1.0), 'F_fluid_feed': ('P_feed_meas', 1.0), 'F_fluid_perm': ('P_perm_meas', -1.0)}
        st_avg = ft[:, 0]
        got = {}
        for fcol, (pname, sgn) in fmap.items():
            c = table_col(fn, ft, fcol)
            if c is not None and pname in W:
                got[pname] = sgn * c / R['AREA']
        if got:
            st_pp = W['step']
            for nm in names[1:]:
                if nm in got:
                    W[nm] = got[nm]
                elif nm in W:
                    W[nm] = np.interp(st_avg, st_pp, W[nm])
            W['step'] = st_avg
            W['source'] = 'piston_force_avg (block means per volume_freq window) / A'
    plat = W['step'] >= (plat_from if plat_from is not None else W['step'][0] + 0.75 * (W['step'][-1] - W['step'][0]))
    if plat.sum() < 1:
        plat[-1] = True
    W['plat_mask'] = plat
    W['plat'] = {nm: float(np.mean(W[nm][plat])) for nm in names[1:] if nm in W}
    W['plat_ci'] = {}
    for nm in names[1:]:
        if nm in W and nm.endswith('_meas') and plat.sum() >= 4:
            m, lo, hi, *_ = block_bootstrap_ci(W[nm][plat], cfg.ci_level)
            W['plat_ci'][nm] = (m, lo, hi)
    en, et = read_piston_table(cfg.path('permeation', lvl))
    if et.size and en is not None and any(n.lower() == 'dv_total' for n in en):
        W['exp_step'] = et[:, 0]
        for key in ('dV_feed', 'dV_perm', 'dV_total', 'z_feed', 'z_perm'):
            W[key] = table_col(en, et, key)
    return W


def _wet_panel(ax, cfg, R, L, col=None, label_prefix='', applied=True, shade=True):
    W = L.get('wet')
    if W is None:
        return False
    st = W['step']
    pal = {'P_feed_meas': (WONG['vermillion'], '-'), 'P_perm_meas': (WONG['blue'], '-'),
           'P_load_meas': (WONG['orange'], '-'), 'P_feed_app': (WONG['vermillion'], '--'), 'P_perm_app': (WONG['blue'], '--')}
    for nm in W['names'][1:]:
        if nm not in W or (nm.endswith('_app') and not applied):
            continue
        c, ls = pal.get(nm, ('k', '-'))
        if col is not None:
            c = col
        y = W[nm]
        lab = label_prefix + nm.replace('_meas', ' measured').replace('_app', ' applied')
        if nm.endswith('_meas'):
            ax.plot(st, y, ls, color=c, lw=0.9, alpha=0.3)
            ax.plot(st, rolling_mean(y, cfg.roll_win), ls, color=c, lw=2.2, alpha=0.95,
                    label=lab + f'  (plateau {sig(W["plat"][nm])})')
            if 'P_load' not in nm and W.get('source', '').startswith('piston_force_avg'):
                lab_src = 'block-averaged $F_{\\rm fluid}/A$'
                if not any(h.get_label() == lab_src for h in ax.get_lines()):
                    ax.plot([], [], '-', color='0.5', lw=2.2, label=lab_src)
        else:
            ax.plot(st, y, ls, color=c, lw=1.4, alpha=0.8, label=lab)
    if shade:
        ax.axvspan(float(st[W['plat_mask']][0]), float(st[-1]), color=WONG['green'], alpha=0.08)
    ax.axhline(cfg.P_BARO, color='k', ls=':', lw=1.0, alpha=0.6)
    return True


def fig_wet_pistons(cfg, R, L):
    """Two-piston compression: bath pressure check over the level -- P_feed and
    P_perm = F_fluid/(lx ly) on the wet pistons (rolling mean) vs the applied P_target,
    plus the load-piston load P_load.  Shaded = plateau window."""
    if L.get('wet') is None:
        print('wet-piston figure skipped (one-piston run: no piston_pressure file)')
        return None
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    _wet_panel(ax, cfg, R, L)
    ax.set_xlabel('time step')
    ax.set_ylabel(r'$P = F_{\rm fluid}/(l_x l_y)$  (LJ)')
    ax.set_title(f'Wet pistons hold the bath at $P_{{\\rm target}}={sig(cfg.P_BARO)}$;  '
                 f'load piston = network load   ($\\varepsilon={L["lvl"]}$)', fontsize=14)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=11)
    return _save(fig, cfg, 'wet_piston_pressures', L['lvl'])


def fig_solvent_expelled(cfg, R, L):
    """Two-piston compression: solvent expelled from the gel = A * (feed-piston rise
    + permeate-piston descent) vs step, next to the applied strain * L0 * A guide."""
    W = L.get('wet')
    if W is None or W.get('dV_total') is None:
        print('solvent-expelled figure skipped (no permeation_<level> file)')
        return None
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    st = W['exp_step']
    ax.plot(st, W['dV_total'], '-', color=WONG['green'], lw=2.4, label=r'total  $A\,(\Delta z_{\rm feed} - \Delta z_{\rm perm})$')
    ax.plot(st, W['dV_feed'], '--', color=WONG['vermillion'], lw=1.6, label='into the feed reservoir')
    ax.plot(st, W['dV_perm'], '--', color=WONG['blue'], lw=1.6, label='into the permeate reservoir')
    if 'L_bb' in L.get('Dc', {}) if L.get('Dc') else False:
        pass
    ax.axhline(0, color='k', ls=':', lw=1, alpha=0.5)
    ax.set_xlabel('time step')
    ax.set_ylabel(r'solvent expelled  ($\sigma^3$)')
    ax.set_title(f'Solvent expelled since seating (wet-piston displacement), $\\varepsilon={L["lvl"]}$'
                 f'  |  final {sig(W["dV_total"][-1])} $\\sigma^3$ = {sig(W["dV_total"][-1] / R["AREA"])} $\\sigma$ of feed rise', fontsize=13)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=12)
    return _save(fig, cfg, 'solvent_expelled', L['lvl'])


def fig_wet_pistons_sweep(cfg, R, levels):
    hw = [L for L in levels if L.get('wet') is not None]
    if not hw:
        print('wet-piston sweep figure skipped (one-piston run)')
        return None
    fig, (axP, axE) = plt.subplots(1, 2, figsize=(18, 6), constrained_layout=True)
    for i, L in enumerate(hw):
        W = L['wet']
        st = W['step']
        for nm, ls in (('P_feed_meas', '-'), ('P_perm_meas', '--')):
            if nm in W:
                axP.plot(st, rolling_mean(W[nm], cfg.roll_win), ls, color=level_color(levels.index(L)), lw=2.0, alpha=0.9)
        if W.get('dV_total') is not None:
            axE.plot(W['exp_step'], W['dV_total'], '-', color=level_color(levels.index(L)), lw=2.0)
    axP.axhline(cfg.P_BARO, color='k', ls=':', lw=1.0, alpha=0.6, label=f'$P_{{\\rm target}}={sig(cfg.P_BARO)}$')
    axP.set_xlabel('time step')
    axP.set_ylabel(r'$P_{\rm bath}$ on the wet pistons  (LJ)')
    axP.set_title('(a) bath check per level: feed (solid), permeate (dashed), rolling means', fontsize=13)
    axP.grid(alpha=0.3)
    h = level_handles(levels) + [Line2D([0], [0], color='0.3', ls='-', lw=2, label='feed'),
                                 Line2D([0], [0], color='0.3', ls='--', lw=2, label='permeate')]
    smart_legend(axP, handles=h, fontsize=11)
    axE.set_xlabel('time step')
    axE.set_ylabel(r'solvent expelled  ($\sigma^3$)')
    axE.set_title('(b) solvent expelled since seating (cumulative over the sweep)', fontsize=13)
    axE.grid(alpha=0.3)
    smart_legend(axE, handles=level_handles(levels), fontsize=11)
    return _save(fig, cfg, 'sweep_wet_pistons')


# ---------------------------------------------------------------------------
#  permeation mode: loader
# ---------------------------------------------------------------------------
def load_permeation(cfg, R, verbose=True):
    """Everything for a two-piston PERMEATION run (one continuous drive, no levels):
    production stress stacks + Terzaghi split (feed-reservoir baseline), density
    evolution, both wet pistons (z, v, F, P measured vs applied), the permeate flux
    Q_perm(t) from the permeate-piston displacement with its steady-state
    (block-bootstrap) plateau and linear fit, the bead-count cross-check, and the
    permeability k = Q L/(A dP) (applied and measured dP).  Returns a dict P shaped
    like a level dict so the shared evolution figures (fig_total_stress,
    fig_network_stress, fig_partial_stress) work on it."""
    say = print if verbose else (lambda *a, **k: None)
    z = R['z']
    P = dict(lvl=None, eps=np.nan, z=z, mode='permeation')
    P['stress'] = {}
    for comp in COMPONENTS:
        S = _component_stacks(cfg, comp, lvl=None)
        if S is None:
            if comp == 'zz':
                raise FileNotFoundError('production sigmazz_{polymer,solvent} files missing -- run the sync cell')
            say(f'  NOTE: no sigma{comp} production files -> {comp} panels skipped')
            continue
        P['stress'][comp] = S
    zz = P['stress']['zz']
    ts = zz['ts']
    P['ts'] = ts
    t0, t1 = float(ts[0]), float(ts[-1])
    P['halt_ts'] = int(t1 - cfg.plateau_frac * (t1 - t0))
    P['evol_from'] = int(t0)
    P['plat'] = ts >= P['halt_ts']
    say(f'production: {len(ts)} stress snapshots, steps {int(t0)} -> {int(t1)}; steady window: steps >= {P["halt_ts"]}')

    # ---- membrane bounds from the final polymer stress; interior ----
    spf = np.abs(zz['p'][-1])
    mem = spf > cfg.gel_thresh * float(np.nanmax(spf))
    P['z_mem_lo'] = float(z[mem].min()) if mem.any() else R['z_gel_lo']
    P['z_mem_hi'] = float(z[mem].max()) if mem.any() else R['z_gel_hi']
    P['in_mem'] = (z >= P['z_mem_lo']) & (z <= P['z_mem_hi'])
    P['interior'] = P['in_mem'] & (z >= P['z_mem_lo'] + cfg.wall_margin) & (z <= P['z_mem_hi'] - cfg.wall_margin)

    # ---- pistons: position / velocity / force / pressure ------------------
    pn, pt = read_piston_table(cfg.path('piston_position'))
    P['pist'] = {}
    if pt.size:
        P['pist']['step'] = pt[:, 0]
        for lab in ('feed', 'perm'):
            c = table_col(pn, pt, f'z_{lab}')
            if c is not None:
                P['pist'][f'z_{lab}'] = c
    vn, vt = read_piston_table(cfg.path('piston_velocity'))
    if vt.size:
        for lab in ('feed', 'perm'):
            c = table_col(vn, vt, f'vz_{lab}')
            if c is not None:
                P['pist'][f'vz_{lab}'] = c
        P['pist']['vstep'] = vt[:, 0]
    zf = P['pist'].get('z_feed')
    if zf is not None:
        st = P['pist']['step']
        P['piston_z_at'] = lambda t, st=st, zf=zf: float(np.interp(t, st, zf))
        P['z_pist'] = float(zf[-1])
    else:
        P['piston_z_at'] = lambda t: R['z_feed']
        P['z_pist'] = R['z_feed']
    P['shade_hi'] = P['z_pist'] if (cfg.GEL_SHADE_TO_PISTON and np.isfinite(P['z_pist'])) else None
    P['wet'] = load_wet_pistons(cfg, R, None, plat_from=P['halt_ts'])

    # ---- Terzaghi split with the feed-reservoir baseline (+ the permeate side) ----
    # both windows are rebuilt per snapshot from the MEASURED piston planes: in a
    # permeation run the feed piston descends and the permeate piston rises all along
    # (2026-09-22); each bin stays >= res_wall_margin clear of its piston sheet
    P['wetz'] = _wet_piston_z_funcs(cfg, R, None)
    P['plates'] = {'support': (R['z_support'], R['z_support'])}
    P['plates'].update({k: v for k, v in plate_tracks(R, dict(wetz=P['wetz'], piston_pos=None, support_pos=None)).items()
                        if k in ('feed piston', 'permeate piston')})
    say('plates over the run (start -> end): ' + fmt_plates(P['plates']))
    bw = baseline_masks_at(z, cfg, R, ts, P['wetz']['feed_at'], P['z_mem_hi'])
    P['bw'] = bw
    bw_perm = None
    if np.isfinite(R.get('z_perm', np.nan)):
        z_bot = min(P['z_mem_lo'], R['z_gel_lo'])
        rows = [reservoir_mask(z, cfg, 'perm', P['wetz']['perm_at'](t), z_bot) for t in ts]
        bw_perm = np.array(rows, bool) if all(r.any() for r in rows) else None
    P['bw_perm'] = bw_perm
    say(f'pore baseline windows (tracking the pistons): feed {mask_span(z, bw[0])} -> {mask_span(z, bw[-1])}'
        + (f';  permeate {mask_span(z, bw_perm[0])} -> {mask_span(z, bw_perm[-1])}' if bw_perm is not None else ''))
    for comp, S in P['stress'].items():
        Rs = R['stress'].get(comp)
        sd_bin = Rs['sd_bin'] if (Rs is not None and Rs['t'].shape[0] >= 2 and np.any(Rs['sd_bin'] > 0)) \
            else np.full(len(z), float(np.nanstd(S['t'][-1][bw[-1]])))
        S['net'], S['pore'], S['pore_half'], S['net_half'] = terzaghi_split(S['t'], bw, sd_bin, cfg.ci_level)
        S['ref_net'] = Rs['net_interior'] if Rs is not None else 0.0
        S['net_plat'] = np.nanmean(S['net'][P['plat']], axis=0)
        S['t_plat'] = np.nanmean(S['t'][P['plat']], axis=0)
        S['net_mem'] = np.array([np.nanmean(S['net'][i][P['in_mem']]) for i in range(len(ts))])
        if bw_perm is not None:
            S['pore_perm'] = np.array([float(np.nanmean(S['t'][i][bw_perm[i]])) for i in range(len(ts))])
    say(f"reservoir sigma_zz baselines (profiles): feed {zz['pore'][0]:.4f} -> {zz['pore'][-1]:.4f}"
        + (f";  permeate {zz['pore_perm'][0]:.4f} -> {zz['pore_perm'][-1]:.4f}" if 'pore_perm' in zz else ''))

    # ---- density evolution ----
    fd = cfg.path('solvent_density_z')
    if fd.exists():
        d_ts, d_z, _, d_m = load_density(fd)
        P.update(dens_ts=d_ts, dens_z=d_z, dens_m=d_m)

    # ---- gel thickness (Rg) over the run ----
    G = load2c(cfg.path('gel_dimensions_rg'), 4)
    if G is not None:
        P['L_step'], P['L_rg'] = G[:, 0], G[:, 3]
        w = P['L_step'] >= P['halt_ts']
        P['L_gel'] = float(np.mean(P['L_rg'][w])) if w.any() else float(P['L_rg'][-1])
    else:
        P['L_gel'] = float(P['z_mem_hi'] - P['z_mem_lo'])

    # ---- flux: permeation file ----
    en, et = read_piston_table(cfg.path('permeation'))
    area = R['AREA']
    P['area'] = area
    P['flux'] = None
    if et.size and en is not None:
        F = dict(step=et[:, 0])
        for key, pre in (('z_feed', 'z_feed'), ('z_perm', 'z_perm'), ('P_feed_meas', 'P_feed_meas'),
                         ('P_perm_meas', 'P_perm_meas'), ('Q', 'Q_perm'), ('N', 'N_permeate')):
            F[key] = table_col(en, et, pre)
        st = F['step']
        prod = st >= P['evol_from']                      # the file also spans the ramp
        F['prod'] = prod
        rho0 = R.get('rho_s0', np.nan)
        # ---- PRIMARY flux (2026-09-23): slope of the permeate bead count N_permeate(t) ----
        # N is an integral of the flow, so its slope over long windows averages the
        # thermal jitter of the piston away; the block-averaged piston velocity does not
        # (v_thermal ~ v_drift for a 2.3e7-mass sheet).  Uncertainty = standard error of
        # the slopes of independent q_win_steps windows (t-interval), NOT a block std.
        cnt = load2c(cfg.path('permeate_count'), 2)
        if cnt is not None and (cnt[:, 0] >= P['evol_from']).sum() >= 8:
            cw = cnt[:, 0] >= P['evol_from']
            n_st, n_val, F['count'] = cnt[cw, 0], cnt[cw, 1], cnt
        elif F['N'] is not None and prod.sum() >= 8:
            n_st, n_val = st[prod], F['N'][prod]
        else:
            n_st = n_val = None
        if n_st is not None and np.isfinite(rho0) and rho0 > 0:
            NS = plateau_window_slopes(n_st, n_val, cfg.q_win_steps, cfg.dt_lj, cfg.plateau_frac_auto, cfg.ci_level)
            for k in ('mean', 'se', 'lo', 'hi', 'slope_all'):
                NS[k] = NS[k] / rho0
            NS['slopes'] = NS['slopes'] / rho0
            NS['rho0'] = float(rho0)
            F['Q_N'] = NS
            win = st >= NS['step0']
            F['win'] = win
            F['win_step0'] = NS['step0']
            # z_perm slope over the same windows: Q_z = -A dz_perm/dt (piston displacement, integrated)
            if F['z_perm'] is not None and win.sum() >= 3:
                ZS = slope_estimate(st[win], F['z_perm'][win], cfg.q_win_steps, cfg.dt_lj, cfg.ci_level)
                for k in ('mean', 'se', 'lo', 'hi', 'slope_all'):
                    ZS[k] = -area * ZS[k]
                ZS['se'] = abs(ZS['se'])                 # the sign flip above must not touch the SE (fixed 2026-09-24)
                ZS['lo'], ZS['hi'] = min(ZS['lo'], ZS['hi']), max(ZS['lo'], ZS['hi'])
                ZS['slopes'] = -area * ZS['slopes']
                F['Q_z'] = ZS
                F['Q_fit'] = ZS['slope_all']
        # ---- legacy: block-averaged piston velocity, drift-tested block-bootstrap plateau ----
        if F['Q'] is not None and prod.sum() >= 4:
            PW = plateau_window(st[prod], F['Q'][prod], cfg.plateau_frac_auto, cfg.ci_level)
            F['Q_ss'] = PW
            if 'win' not in F:
                F['win'] = st >= PW['step0']
                F['win_step0'] = PW['step0']
        if 'win' in F:
            win = F['win']
            # measured dP over the steady window (block means per volume_freq window)
            if F['P_feed_meas'] is not None and F['P_perm_meas'] is not None and win.sum() >= 4:
                dpm = F['P_feed_meas'][win] - F['P_perm_meas'][win]
                m, lo, hi, *_ = block_bootstrap_ci(dpm, cfg.ci_level)
                F['dP_meas'], F['dP_meas_lo'], F['dP_meas_hi'] = float(m), float(lo), float(hi)
            # applied dP
            W = P['wet']
            if cfg.DP_PISTON is not None:
                F['dP_app'] = float(cfg.DP_PISTON)
            elif W is not None and 'P_feed_app' in W and 'P_perm_app' in W:
                F['dP_app'] = float(np.mean((W['P_feed_app'] - W['P_perm_app'])[W['plat_mask']]))
            else:
                F['dP_app'] = np.nan
            # permeability k = Q L / (A dP)   (LJ: sigma^5/(eps tau) -> kappa = k/eta convention of the notebooks)
            Lg = P['L_gel']
            def _k(Q, Qlo, Qhi, dP, dPlo=None, dPhi=None):
                if not (np.isfinite(Q) and np.isfinite(dP)) or dP == 0:
                    return None
                k = Q * Lg / (area * dP)
                rq = 0.5 * (Qhi - Qlo) / abs(Q) if np.isfinite(Qlo) and Q != 0 else 0.0
                rp = 0.5 * (dPhi - dPlo) / abs(dP) if (dPlo is not None and np.isfinite(dPlo)) else 0.0
                h = abs(k) * np.sqrt(rq ** 2 + rp ** 2)
                return dict(k=float(k), lo=float(k - h), hi=float(k + h))
            F['k'] = {}
            if 'Q_N' in F:
                NS = F['Q_N']
                F['k']['N_applied'] = _k(NS['mean'], NS['lo'], NS['hi'], F['dP_app'])
                if 'dP_meas' in F:
                    F['k']['N_measured'] = _k(NS['mean'], NS['lo'], NS['hi'], F['dP_meas'], F['dP_meas_lo'], F['dP_meas_hi'])
            if 'Q_z' in F:
                ZS = F['Q_z']
                F['k']['z_applied'] = _k(ZS['mean'], ZS['lo'], ZS['hi'], F['dP_app'])
            F['k'] = {kk: v for kk, v in F['k'].items() if v is not None}
            if 'Q_N' in F:
                NS = F['Q_N']
                say(f"Q_perm PRIMARY (N_permeate slope): {NS['mean']:.4e} +/- {NS['se']:.2e} (SE over {NS['n_win']} windows of "
                    f"{NS['win_steps']} steps; {cfg.ci_level:.0%} CI [{NS['lo']:.4e}, {NS['hi']:.4e}]) sigma^3/tau; "
                    f"window last {NS['frac']:.0%} (step >= {int(NS['step0'])}); single fit {NS['slope_all']:.4e}"
                    + ('  DRIFT WARNING' if NS.get('warn') else ''))
            if 'Q_z' in F:
                ZS = F['Q_z']
                say(f"Q_perm from the z_perm slope (same windows): {ZS['mean']:.4e} +/- {ZS['se']:.2e}")
            say(f"dP applied = {F['dP_app']:.4f};  dP measured (wet pistons) = "
                + (f"{F['dP_meas']:.4f} [{F['dP_meas_lo']:.4f}, {F['dP_meas_hi']:.4f}]" if 'dP_meas' in F else 'n/a')
                + f";  L_gel (Rg) = {Lg:.2f};  A = {area:.1f}")
            say('permeability kappa = Q L/(A dP_ext):  ' + '   '.join(f"{kk}: {v['k']:.4e} [{v['lo']:.4e}, {v['hi']:.4e}]" for kk, v in F['k'].items()))
        P['flux'] = F
    else:
        say('NOTE: permeation_<stem>.dat missing -> flux / permeability skipped')
    return P


# ---------------------------------------------------------------------------
#  permeation mode: figures
# ---------------------------------------------------------------------------
def fig_perm_density(cfg, R, P):
    """Solvent mass density evolution (cividis, bold final) with the zero-flux reference
    (dropped from the notebook 2026-09-24: the mass-fraction panel carries the same field)."""
    if 'dens_m' not in P:
        print('density figure skipped (no solvent_density_z file)')
        return None
    fig, ax = plt.subplots(figsize=(11, 6.5), constrained_layout=True)
    ts, ev = post_halt(cfg, P, P['dens_ts'], P['dens_m'])
    ref = None
    if 'dens_m' in R:
        m, lo, hi = mean_ci(np.array([np.interp(P['dens_z'], R['dens_z'], row) for row in R['dens_m']]), cfg.ci_level)
        ref = (m, lo, hi)
    plot_evolution(ax, cfg, R, P, P['dens_z'], ts, ev, r'$\rho_s(z,t)$', 'Solvent mass density: reference -> permeation drive',
                   ref=ref, annotate=False)
    if np.isfinite(R.get('rho_s0', np.nan)):
        ax.axhline(R['rho_s0'], color='k', ls=':', lw=1.2, alpha=0.6)
    return _save(fig, cfg, 'perm_solvent_density_evolution')


def fig_perm_volfrac(cfg, R, P):
    """Solvent volume fraction under permeation: (a) the zero-flux reference and (b) the
    steady state (trailing-window means of the mass-fraction, Voronoi and lambda-
    calibrated estimators); (c) the calibration pressure P_local(z) handed to
    lambda(phi_p, P) -- under P_CAL_MODE='pore' the pore-pressure ramp from the feed
    to the permeate reservoir baseline across the membrane; (d) the resulting
    lambda(z).  Needs add_perm_volume_fractions."""
    if P.get('phi_mf') is None and P.get('phi_vor') is None:
        print('volume-fraction figure skipped (run add_perm_volume_fractions first)')
        return None
    fig, axes = plt.subplots(2, 2, figsize=(18, 12.5), constrained_layout=True)
    fig.suptitle(f'Solvent volume fraction under permeation  (P_CAL_MODE = {cfg.P_CAL_MODE!r}, '
                 f'$P_{{\\rm CAL}}={cfg.P_CAL}$)  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    _phi_panel(axes[0, 0], cfg, R, R, None, '(a) zero-flux reference (both reservoirs at $P_{\\rm target}$)')
    _phi_panel(axes[0, 1], cfg, R, P, P, '(b) steady permeation (trailing-window mean)', delta_from=R)
    zx = zn(R, R['z'])
    zz = P['stress']['zz']
    # (c) the calibration pressure per bin
    ax = axes[1, 0]
    if P.get('P_cal') is not None:
        m, lo, hi = P['P_cal']
        ax.fill_between(zx, lo, hi, color=WONG['orange'], alpha=0.22, lw=0, zorder=2)
        ax.plot(zx, m, '-', lw=2.4, color=WONG['orange'], label=r'$P_{\rm local}(z)$, steady permeation', zorder=3)
        if R.get('P_cal') is not None:
            ax.plot(zx, R['P_cal'][0], '--', lw=2.0, color='k', alpha=0.8, label=r'$P_{\rm local}(z)$, reference', zorder=3)
        pf = float(np.nanmean(zz['pore'][P['plat']]))
        ax.axhline(pf, color=WONG['blue'], ls=':', lw=1.4, alpha=0.8, label=f'feed baseline {pf:.3f}')
        if 'pore_perm' in zz:
            pp = float(np.nanmean(zz['pore_perm'][P['plat']]))
            ax.axhline(pp, color=WONG['vermillion'], ls=':', lw=1.4, alpha=0.8, label=f'permeate baseline {pp:.3f}')
        shade_gel(ax, R, P)
        mark_walls(ax, R, P)
        finish_axes(ax, r'$P_{\rm local}$  (LJ)', r'(c) calibration pressure handed to $\lambda(\phi_p, P)$')
        vals = np.concatenate([m[np.isfinite(m)], [pf, cfg.P_BARO]])
        span = max(float(vals.max() - vals.min()), 0.05)
        ax.set_ylim(vals.min() - 0.6 * span, vals.max() + 0.6 * span)
        Pmin, Pmax = _calib_range(R)
        smart_legend(ax, fontsize=11)
        annotate_box(ax, f'membrane mean = {fmt_mu(m[P["in_mem"]])}\n'
                         f'calibration sweep: $P$ = {Pmin:g} ... {Pmax:g}', loc='lower left', fontsize=12)
    else:
        ax.text(0.5, 0.5, 'unavailable (no calibration / no tessellated frames)', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(r'(c) calibration pressure handed to $\lambda(\phi_p, P)$')
    # (d) lambda(z)
    ax = axes[1, 1]
    if P.get('lam') is not None:
        m, lo, hi = P['lam']
        ax.fill_between(zx, lo, hi, color=WONG['green'], alpha=0.22, lw=0, zorder=2)
        ax.plot(zx, m, '-', lw=2.4, color=WONG['green'], label=r'$\lambda(z)$, steady permeation', zorder=3)
        if R.get('lam') is not None:
            ax.plot(zx, R['lam'][0], '--', lw=2.0, color='k', alpha=0.8, label=r'$\lambda(z)$, reference', zorder=3)
        ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6)
        shade_gel(ax, R, P)
        mark_walls(ax, R, P)
        finish_axes(ax, r'$\lambda$', r'(d) $\lambda(\phi_p^{\rm vor}, P_{\rm local})$ per bin  (reservoir bins: exactly 1)')
        fin = m[np.isfinite(m)]
        if fin.size:
            span = max(float(fin.max() - fin.min()), 0.02)
            ax.set_ylim(min(fin.min(), 1.0) - 0.5 * span, fin.max() + 0.5 * span)
        rt = R['CALIB'].get('raw_table', {}) if R.get('CALIB') else {}
        win = ''
        if rt.get('phi_p_vor'):
            pv = np.asarray(rt['phi_p_vor'], float)
            win = f'\ncalibration data window: $\\phi_p^{{\\rm vor}}$ = {pv.min():.2f} ... {pv.max():.2f}'
        smart_legend(ax, fontsize=11)
        annotate_box(ax, f'interior mean = {fmt_mu(m[P["interior"]])}' + win, loc='lower left', fontsize=12)
    else:
        ax.text(0.5, 0.5, 'unavailable (no calibration / no tessellated frames)', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(r'(d) $\lambda(z)$')
    return _save(fig, cfg, 'perm_volfrac_profiles')


def fig_perm_volfrac_evolution(cfg, R, P):
    """phi_s evolution over the drive (cividis, final bold; zero-flux reference dashed):
    (a) mass fraction from every density snapshot, (b) the lambda-calibrated Voronoi
    fraction on the tessellated frames, each with its own colour bar (2026-09-24; the
    raw Voronoi panel is kept only when no calibration artifact exists).  The box gives
    the in-gel change from the reference.  Needs add_perm_volume_fractions."""
    panels = []
    if P.get('mf_stack') is not None:
        panels.append(('phi_mf', P['mf_ts'], P['mf_stack'], r'mass fraction $\phi_s^{\rm mf}=\rho_s/\rho_{s,0}$'))
    vf = P.get('vf')
    if vf is not None and vf.get('phi_cal') is not None:
        panels.append(('phi_cal', vf['ts'], vf['phi_cal'],
                       r'$\lambda$-calibrated Voronoi $\phi_s^{\rm cal}$' + f'  (P_CAL_MODE = {cfg.P_CAL_MODE!r})'))
    elif vf is not None and vf.get('phi_vor') is not None:
        panels.append(('phi_vor', vf['ts'], vf['phi_vor'], r'Voronoi $\phi_s^{\rm vor}$ (no calibration artifact)'))
    if not panels:
        print('volume-fraction evolution skipped (run add_perm_volume_fractions first)')
        return None
    fig, axes = plt.subplots(1, len(panels), figsize=(8.5 * len(panels), 6.5), constrained_layout=True, squeeze=False)
    fig.suptitle(f'Solvent volume fraction evolution: zero-flux reference $\\rightarrow$ permeation drive  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    for i, (ax, (key, ts, stack, title)) in enumerate(zip(axes[0], panels)):
        if key == 'phi_mf':
            ts, ev = post_halt(cfg, P, ts, stack)
        else:
            ts, ev = np.asarray(ts, float), np.asarray(stack, float)
        plot_evolution(ax, cfg, R, P, R['z'], ts, ev, r'$\phi_s$', f'({"abc"[i]}) ' + title, ref=R.get(key), annotate=False, legend=False)
        ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6)
        ax.set_ylim(0, 1.15)
        h, lab = ax.get_legend_handles_labels()
        if h:
            smart_legend(ax, handles=h, labels=lab, fontsize=11)
        if P.get(key) is not None and R.get(key) is not None:
            d = np.nanmean(P[key][0][P['interior']]) - np.nanmean(R[key][0][R['interior']])
            annotate_box(ax, r'$\Delta\phi_s$ (steady $-$ reference), in gel: ' + f'{d:+.2g}', loc='lower right', fontsize=12)
    return _save(fig, cfg, 'perm_volfrac_evolution')


def fig_perm_pistons(cfg, R, P):
    """(a) both wet-piston positions vs step, (b) their measured vs applied pressures
    (+ the reservoir virial pressures pressure_feed / pressure_permeate, dotted)."""
    if not P['pist']:
        print('piston figure skipped (no piston_position file)')
        return None
    fig, (axZ, axP) = plt.subplots(1, 2, figsize=(18, 6), constrained_layout=True)
    st = P['pist']['step']
    for lab, col in (('feed', WONG['vermillion']), ('perm', WONG['blue'])):
        zc = P['pist'].get(f'z_{lab}')
        if zc is not None:
            axZ.plot(st, zc - zc[0], '-', color=col, lw=2.2, label=f'{lab} piston  ($z_0={sig(zc[0])}$)')
    axZ.axvspan(P['evol_from'], float(st[-1]), color=WONG['green'], alpha=0.06, label='production')
    axZ.axhline(0, color='k', ls=':', lw=1, alpha=0.5)
    axZ.set_xlabel('time step')
    axZ.set_ylabel(r'$z - z_0$  ($\sigma$)')
    axZ.set_title('(a) wet-piston displacement (feed descends, permeate retreats)', fontsize=13)
    axZ.grid(alpha=0.3)
    smart_legend(axZ, fontsize=11)
    if _wet_panel(axP, cfg, R, P):
        # reservoir virial pressures: block means at the pf cadence since 2026-09-23
        # (one point per volume_freq window; older runs have ~10 points at the stress cadence)
        for name, col in (('pressure_feed', WONG['vermillion']), ('pressure_permeate', WONG['blue'])):
            f = cfg.path(name)
            if f.exists():
                a = load2c(f, 2)
                if a is not None:
                    dense = a.shape[0] > 3 * cfg.roll_win
                    if dense:
                        axP.plot(a[:, 0], a[:, 1], ':', color=col, lw=0.8, alpha=0.25)
                    axP.plot(a[:, 0], rolling_mean(a[:, 1], cfg.roll_win) if dense else a[:, 1], ':', color=col, lw=2.0, alpha=0.95,
                             label=f'{name.split("_")[1]} reservoir virial $P_{{\\rm res}}$' + (' (rolling mean)' if dense else ''))
    src = (P.get('wet') or {}).get('source', '')
    axP.set_xlabel('time step')
    axP.set_ylabel('pressure  (LJ)')
    axP.set_title('(b) piston pressure: measured ' + ('block-averaged ' if src.startswith('piston_force_avg') else '')
                  + '$F_{\\rm fluid}/A$ vs applied (dashed); reservoir virial (dotted)', fontsize=12)
    axP.grid(alpha=0.3)
    smart_legend(axP, fontsize=10)
    return _save(fig, cfg, 'perm_pistons')


def fig_perm_flux(cfg, R, P):
    """(a) Q_perm(t) from the block-averaged permeate-piston velocity (context) with the
    steady window and the two estimates: N_perm slope (primary, +/- SE over independent
    windows) and the z_perm slope (the legacy block mean was dropped 2026-09-24);
    (b) N_perm(t) with the per-window fits that give the primary value."""
    F = P.get('flux')
    if F is None or (F.get('Q') is None and F.get('Q_N') is None):
        print('flux figure skipped (no permeation file)')
        return None
    fig, (axQ, axN) = plt.subplots(1, 2, figsize=(18, 6), constrained_layout=True)
    st, Q = F['step'], F['Q']
    if Q is not None:
        axQ.plot(st, Q, '-', color=WONG['green'], lw=1.0, alpha=0.35, label=r'$Q=A\,dz_{\rm perm}/dt$ (block-avg piston velocity)')
        axQ.plot(st, rolling_mean(Q, cfg.roll_win), '-', color=WONG['green'], lw=2.4, label=f'rolling mean ({cfg.roll_win})')
    if 'win_step0' in F:
        axQ.axvspan(F['win_step0'], float(st[-1]), color=WONG['green'], alpha=0.10, label='steady window (auto, slope drift test)')
    if 'Q_N' in F:
        NS = F['Q_N']
        axQ.axhline(NS['mean'], color='k', ls='--', lw=1.8,
                    label=f"$N_{{\\rm perm}}$ slope: $Q$ = {fmt_val_unc(NS['mean'], NS['se'])}  (SE, {NS['n_win']} windows)")
        axQ.axhspan(NS['lo'], NS['hi'], color='k', alpha=0.08)
    if 'Q_z' in F:
        ZS = F['Q_z']
        axQ.axhline(ZS['mean'], color=WONG['reddishpurple'], ls=':', lw=1.6,
                    label=f"$z_{{\\rm perm}}$ slope: {fmt_val_unc(ZS['mean'], ZS['se'])}")
    axQ.axvline(P['evol_from'], color='k', ls=':', lw=1, alpha=0.5)
    axQ.set_xlabel('time step')
    axQ.set_ylabel(r'$Q_{\rm perm}$  ($\sigma^3/\tau$)')
    axQ.set_title('(a) permeate flux: piston-velocity trace (context) and the slope estimates', fontsize=13)
    axQ.grid(alpha=0.3)
    smart_legend(axQ, fontsize=10)
    cnt = F.get('count')
    if cnt is not None:
        n_st, n_val = cnt[:, 0], cnt[:, 1]
    elif F.get('N') is not None:
        n_st, n_val = st, F['N']
    else:
        n_st = n_val = None
    if n_st is not None:
        axN.plot(n_st, n_val - n_val[0], '-', color='0.3', lw=1.8, label=r'$\Delta N_{\rm perm}$ (crossed below the support)')
        if 'Q_N' in F:
            NS = F['Q_N']
            cw = n_st >= NS['step0']
            if cw.sum() >= 3:
                sl, ic = np.polyfit(n_st[cw], n_val[cw] - n_val[0], 1)
                axN.plot(n_st[cw], sl * n_st[cw] + ic, '--', color=WONG['vermillion'], lw=2,
                         label=f"single fit $\\rightarrow$ $Q$ = {sig(NS['slope_all'])}")
            # the independent windows: short segments at the fitted slopes
            w = NS['win_steps']
            for j, (c, q) in enumerate(zip(NS['centers'], NS['slopes'])):
                sel = (n_st >= c - 0.5 * w) & (n_st <= c + 0.5 * w)
                if sel.sum() < 2:
                    continue
                m = np.polyfit(n_st[sel], n_val[sel] - n_val[0], 1)
                axN.plot(n_st[sel], np.polyval(m, n_st[sel]), '-', color=WONG['blue'], lw=2.6, alpha=0.85,
                         label=(f'{NS["n_win"]} independent windows of {w} steps $\\rightarrow$ $Q$ = {fmt_val_unc(NS["mean"], NS["se"])}' if j == 0 else None))
            axN.axvspan(NS['step0'], float(n_st[-1]), color=WONG['green'], alpha=0.10)
            axN.set_title(f"(b) permeate bead count: window slopes / $\\rho_{{s,0}}$ = {sig(NS['rho0'])}   (primary $Q$)", fontsize=13)
        else:
            axN.set_title('(b) permeate bead count', fontsize=13)
    axN.set_xlabel('time step')
    axN.set_ylabel(r'$\Delta N_{\rm perm}$  (beads)')
    axN.grid(alpha=0.3)
    smart_legend(axN, fontsize=10)
    return _save(fig, cfg, 'perm_flux')


def fig_perm_permeability(cfg, R, P):
    """kappa = Q L/(A dP_ext): the N_perm-slope Q with the MEASURED (wet-piston) dP_ext --
    the primary value -- and the z_perm-slope Q with the applied dP_ext; CI = the slope's
    SE-based interval (and dP_meas's block bootstrap) in quadrature.  (2026-09-24: the
    N-slope/applied-dP and the legacy block-mean entries are no longer drawn; both stay
    in P['flux']['k'] and the printed summary.)"""
    F = P.get('flux')
    if F is None or not F.get('k'):
        print('permeability figure skipped (no flux)')
        return None
    labs = {'N_measured': r'$N_{\rm perm}$ slope, measured $\Delta P_{\mathrm{ext}}$',
            'z_applied': r'$z_{\rm perm}$ slope, applied $\Delta P_{\mathrm{ext}}$'}
    cols = {'N_measured': WONG['vermillion'], 'z_applied': WONG['reddishpurple']}
    show = [kk for kk in labs if kk in F['k']]
    if not show:
        print('permeability figure skipped (neither the measured-dP nor the z_perm estimate is available)')
        return None
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    for i, kk in enumerate(show):
        v = F['k'][kk]
        ax.errorbar([i], [v['k']], yerr=[[v['k'] - v['lo']], [v['hi'] - v['k']]], fmt='o', ms=12, color=cols[kk],
                    capsize=8, lw=2.5, label=f"{labs[kk]}:  $\\kappa = {sig(v['k'])}$")
    ax.set_xticks(range(len(show)))
    ax.set_xticklabels([labs[kk] for kk in show], fontsize=11)
    ax.set_ylabel(r'$\kappa = Q L/(A\,\Delta P_{\mathrm{ext}})$  (LJ: $\sigma^5/(\epsilon\,\tau)$)')
    ax.set_title(f"Permeability  ($L={sig(P['L_gel'])}\\,\\sigma$, $A={sig(P['area'])}\\,\\sigma^2$, "
                 f"$\\Delta P_{{\\mathrm{{ext}}}}={sig(F['dP_app'])}$ applied)\n{cfg.sim_name}", fontsize=12)
    ax.set_xlim(-0.6, len(show) - 0.4)
    ax.grid(axis='y', alpha=0.3)
    smart_legend(ax, fontsize=11)
    return _save(fig, cfg, 'perm_permeability')


def print_perm_summary(cfg, R, P):
    F = P.get('flux') or {}
    print(f'\nPERMEATION SUMMARY  ({cfg.sim_name})')
    print(f"  gel L (Rg, steady window) = {P['L_gel']:.2f} sigma;  A = {P['area']:.1f};  membrane bins z in "
          f"[{P['z_mem_lo']:.1f}, {P['z_mem_hi']:.1f}]")
    W = P.get('wet')
    if W is not None:
        print('  wet pistons (steady means): ' + '  '.join(f"{k}: {v:.4f}" for k, v in W['plat'].items()))
    if 'Q_N' in F:
        NS = F['Q_N']
        print(f"  Q_perm (N_permeate slope) = {NS['mean']:.4e} +/- {NS['se']:.2e} sigma^3/tau  "
              f"(SE over {NS['n_win']} windows of {NS['win_steps']} steps; Darcy velocity Q/A = {NS['mean'] / P['area']:.3e})")
    if 'Q_z' in F:
        print(f"  Q_perm (z_perm slope)     = {F['Q_z']['mean']:.4e} +/- {F['Q_z']['se']:.2e}")
    if F.get('k'):
        for kk, v in F['k'].items():
            print(f"  kappa ({kk:10s}) = {v['k']:.4e} [{v['lo']:.4e}, {v['hi']:.4e}]")
    S = P.get('psd')
    if S is not None:
        print(f"  geometric porosity (r_probe {cfg.PSD_R_PROBE}), interior: {fmt_mu(S['por'][0][P['interior']])};  "
              f"mean pore diameter {fmt_mu(S['d_mean'][0][P['interior']])} sigma")
    elif 'Q_N' not in F:
        print('  no steady flux (permeation / permeate_count files missing or too short)')
