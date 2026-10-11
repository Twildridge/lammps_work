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
                                    panels for two-piston compression runs; D_c and M from
                                    the polymer displacement under permeation (2026-09-28)

Figure conventions (2026-10-06): cfg.FLIP_Z draws every z-profile with the feed / load piston on
the left (finish_axes / _zlim / flip_z_axis -- plotting only); cfg.PARTIAL_NORM = 'share' adds,
under each total-stress, partial-stress and thermodynamic-pressure profile figure, the version
with the back pressure removed and divided by the driving pressure (partial_norm, section 9).

Physics conventions (see the Notes section at the end of either notebook):
  * total stress sigma^t = sigma_p + sigma_s (group stress/atom, kinetic term included)
  * Terzaghi: sigma' = sigma^t - p_pore, p_pore read per curve from the flat
    reservoir at z/Lz ~ baseline_zf (one scalar per stress component)
  * uniaxial strain (fixed lateral box):  sigma'_zz = M eps,  sigma'_xx = lambda eps,
    lambda = M - 2G  ->  sigma'_zz/sigma'_xx = M/(M - 2G)  ->  G = (sigma'_zz - sigma'_xx)/(2 eps)
  * M_network = (<sigma'_zz>_interior,plateau - sigma'_zz,ref) / eps_applied
    M_piston  = (<F_z/A>_plateau - P_ref) / eps_applied          (M_SUBTRACT_REF; each estimator
    subtracts its own eps = 0 reading -- the piston preload and the profile bias differ by ~0.002)
  * D_c from the consolidation fit of the polymer displacement u_z(z,t) during a hold between
    two drained plates: modes sin(2 m pi zeta) and cos(2 m pi zeta) - 1 decaying as
    exp(-4 m^2 pi^2 D_c t/L^2), tau_1 = L^2/(4 pi^2 D_c)  (2026-09-28; the earlier
    (2 zeta - 1) + cos(k pi zeta) shapes at k^2 pi^2 D_c/L^2 overstated D_c ~5x)
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
from scipy.optimize import minimize_scalar, curve_fit
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter, FormatStrFormatter
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
    # ---- profile-figure conventions (2026-10-06) ---------------------------
    PARTIAL_NORM: str = 'share'           # 'raw' | 'share'.  'share': under each total-stress, partial-stress and
                                          # thermodynamic-pressure profile figure a second figure (same call, file stem
                                          # <stem>_norm) draws the profiles with the back pressure removed and divided
                                          # by the driving pressure -- the form of Marioni et al.'s partial P_zz figure,
                                          # whose permeate is at ~0 where ours is at P_perm = 1.5.  Per bin and snapshot:
                                          #   w_s = sigma_s/sigma_t,   sigma_s* = (sigma_s - w_s P_ref)/dP,
                                          #   sigma_p* = (sigma_p - (1 - w_s) P_ref)/dP,   sigma_t* = (sigma_t - P_ref)/dP
                                          # (partial_norm; the shares come from the stresses alone).  The original
                                          # figures and their PNGs are untouched.  'raw': the original figures only.
    FLIP_Z: bool = True                   # every z-profile figure is drawn with the feed / load piston on the LEFT and
                                          # the permeate / support on the RIGHT, as in Marioni et al.: the z/L axes
                                          # show 1 - z/L running 0 -> 1 and say which side is which; the sigma- and
                                          # zeta-unit axes are inverted and say so.  Plotting only -- no z array, mask,
                                          # window or fit sees it (zn, finish_axes, _zlim, flip_z_axis).  False = z
                                          # increasing to the right (the orientation of every figure before 2026-10-06)
    # ---- pore-size distribution (lib/psd.py, 2026-09-24) ---------------------
    PSD_ENABLE: bool = True               # geometric porosity + PSD on the tessellated frames (needs VOR_ENABLE)
    PSD_R_PROBE: float = 0.5              # probe radius (sigma): void = grid points >= this far from a bead surface
    PSD_GRID: float = 0.5                 # grid spacing (sigma): 0.5 ~ 7 s/frame, 0.25 ~ 2 min/frame (the covering step)
    PSD_DMAX: float = 8.0                 # pore-diameter histogram range (2 r_probe .. PSD_DMAX) ...
    PSD_DBIN: float = 0.25                # ... and bin width (= the resolution of the covering step)
    # ---- permeation: D_c and M from the polymer displacement (2026-09-28) ----------
    PERM_GAP: float = 1.0                 # sigma: the pinned polymer face z^P = z_support + PERM_GAP (the plate-bead
                                          # exclusion; the compression fit's plate gap is the same ~1 sigma).  zeta = (z - z^P)/L_0
    PERM_FORCING: str = 'applied'         # the dP(t) history driving the consolidation model (2026-10-07): 'applied' = the
                                          # piston_pressure ramp P_feed_app - P_perm_app; 'measured' = the reservoir pressure
                                          # difference P_res_feed - P_res_perm (pressure_feed / pressure_permeate, pf cadence,
                                          # smoothed over PERM_FORCING_SMOOTH samples; the applied ramp start stays t = 0 and
                                          # the ramp is assumed to carry no load until the first measured sample).  The
                                          # undamped pistons (damp_prod=0) ring after the ramp: on the quarter gel the
                                          # measured dP overshot to 2x the applied for ~0.3M steps and the gel followed it
                                          # (tau_1 ~0.35M steps); the full gel (tau_1 ~2M) filtered it.
    PERM_FORCING_SMOOTH: int = 5          # rolling-mean window (samples) of the measured dP before differencing
    PERM_L0_MODE: str = 'face'            # the L_0 of the permeation fit (2026-10-07): 'face' = reference polymer-stress edge
                                          # minus Z^P (the pinned face's REFERENCE position z^P - u_0): the thickness of the
                                          # dense gel the material coordinate spans, so zeta reaches 1 at the real top face.
                                          # 'bb' (old): the bounding-box thickness, which also counts the ~3 sigma the face
                                          # floats above the plate at zero flux (perm_2/3: 2.4 % of L_0; the quarter gel:
                                          # 10 %, the fit domain stopped at zeta ~0.8 and the modes ran wild above it).
                                          # The bounding-box L_0 still drives the thickness trace (L_bb(t) - L_bb(0)).
    PERM_L0_STEPS: int = 100000           # steps before the production onset averaged for L_0 (bounding-box thickness; lies
                                          # inside the zero-flux reference window)
    PERM_DC_WINDOW_AVG: bool = True       # the fitted modes are averaged over each ave/chunk window: the deck block-averages
                                          # u_z over the whole nfreq interval and tags the snapshot with the window END
    PERM_DISP_RESET: str = 'auto'         # where the deck reset displace/atom (u_z = 0): 'ramp_start' (decks since
                                          # 2026-09-28: the profile fix is live through the dP ramp) | 'ramp_end' (older
                                          # decks: reset at the production onset) | 'auto' (ramp_start when the first
                                          # snapshot precedes the end of the ramp).  The applied-dP history itself
                                          # (piston_pressure) drives the consolidation model either way.
    PERM_TRACE_SKIP: int = 0              # steps after the ramp start left out of the thickness-trace fit (the zero-IC
                                          # series still starts at the ramp start): use it to skip the piston ringing of
                                          # a short-ramp run (perm_2: ~2M steps)
    PERM_DC_BOUNDS: tuple = (1e-6, 50.0)  # D_c search interval of the permeation fits (sigma^2/tau); wide so that a
                                          # fit driven by a fast non-consolidation transient (the piston ringing of a
                                          # short ramp) shows up as a large D_c rather than pegging at DC_BOUNDS
    PERM_RIGID: str = 'contact'           # rigid shift removed before the consolidation fit (2026-09-28 pm): 'contact'
                                          # subtracts u_0(t), the displacement of the contact layer (the PERM_PIN_NBINS
                                          # lowest populated bins).  At zero flux a solvent layer separates the plate and
                                          # the polymer face (perm_2: the face sat ~3 sigma above z^P); the drag closes it
                                          # and the whole network arrives at the plate displaced by u_0 with no strain
                                          # attached -- the pinned modes cannot carry that offset.  'none' = old behaviour
    PERM_PIN_NBINS: int = 2               # bins of the contact layer averaged for u_0(t)
    PERM_COORDS: str = 'lagrangian'       # 'lagrangian': zeta = (z - u_z - Z^P)/L_0 -- each bin's atoms placed where they
                                          # were at the reset, Z^P = z^P - u_0 the reference position of the pinned face;
                                          # the top bins, which empty as the face comes down, then reach zeta -> 1.
                                          # 'eulerian': the bin centre, zeta = (z - z^P)/L_0 (old behaviour)
    PERM_DC_N_MODES: int = 0              # modes of the permeation profile fit (0 -> DC_N_MODES)
    PERM_M_PRIMARY: str = 'lag'           # the M that feeds kappa = D_c/M and the q(t) check: 'lag' (default since
                                          # 2026-10-03: the steady frames paired per atom with the zero-flux reference,
                                          # u_F = u(top material bin) - u(bottom material bin) -- the rigid drop cancels
                                          # and the contact layer's compaction counts; it agrees with the compression sweep
                                          # to 8 % on perm_3 where the parabola was 20 % stiff; needs local traj_ref +
                                          # traj_stress, else falls back to 'prof') | 'prof' (steady-profile
                                          # parabola, the default since 2026-09-28 pm: it reads the deformation of the
                                          # network alone) | 'bb' (bounding-box thickness change: its ends are the extreme
                                          # beads, which carry ~2 sigma of the rigid drop / tail into |u_F| in perm_2) | 'trace'
    # ---- M and G as increments from the eps = 0 reference ------------------
    M_SUBTRACT_REF: bool = True           # M = (stress - its own eps = 0 reading) / eps for BOTH estimators
                                          # (2026-09-12: the seated piston carries a real preload ~+0.0014 and the
                                          # profile method a constant ~-0.002 bias; absolute M inherits that
                                          # ~0.002 gap.  Piston reference needs piston_force_avg_ref.)
    G_SUBTRACT_REF: bool = True           # use increments relative to the eps = 0 reference state
    RELAX_SYS: bool = True                # 2026-09-27: unrelaxed-hold systematic on M.  An exponential tail is fitted to
                                          # the block-averaged piston P over the last RELAX_TAIL_FRAC of the hold; the
                                          # plateau mean minus the fitted asymptote, / eps, is delta_sys and is added to the
                                          # LOWER bound of BOTH M estimators (the true M is lower), hence to the UPPER bound
                                          # of kappa = D_c/M.  0 when the tail is flat or the fit is not credible.
    RELAX_TAIL_FRAC: float = 0.6          # trailing fraction of the hold the tail is fitted over
    SS_FIT_NPTS: int = 3                  # fig_stress_strain_sweep (2026-09-27): the dotted line of each estimator is
                                          # anchored at its eps = 0 reading and its slope is the through-anchor
                                          # least-squares fit of the first SS_FIT_NPTS levels (the linear regime);
                                          # the higher levels are NOT in the fit, so their stiffening shows as
                                          # departure from the line instead of a negative intercept
    # ---- strain in the denominator of M, G and kappa (2026-10-03) -----------
    G_REF: object = None                  # shear modulus from the shear notebooks (e.g. 0.2): drawn as G_REF / M on the
                                          # reservoir normal-stress sweep figure next to (1 - sigma'_lat/sigma'_zz)/2, and
                                          # (2026-10-07) as the guide dP_th = (4/3)(G/M) dP on the normalised
                                          # thermodynamic-pressure figure (fig_thermo_pressure_norm)
    G_COMP_REF: object = None             # the compression-mode G of the same gel (the lateral network stress of the
                                          # triaxial holds, ~0.05): a second dP_th guide on that figure, so the two G
                                          # estimates can be told apart there (M is trusted, G is not; 2026-10-07)
    M_REF: object = None                  # the M of those guides; None -> the run's own M (permeation: the primary
                                          # estimate of add_perm_displacement; compression: M_net of the level)
    M_STRAIN: str = 'disp'                # 'disp': the slope of the STEADY displacement profile u_z(Z) of the polymer
                                          #   (plateau frames vs the eps = 0 reference, per atom, binned by the reference
                                          #   position Z; u_z = u_0 - eps (Z - Z^P) exactly under uniform stress, and the
                                          #   slope is blind to the ~0.02 of rigid travel -- the gel falling onto the
                                          #   support -- that the applied strain counts; fit_disp_profile).  Falls back to
                                          #   'rg' then 'applied' when the trajectories / disp_z_polymer_cum are missing.
                                          # 'rg': the Rg thickness strain of strain_zz (plateau mean).
                                          # 'applied': the piston-travel target (the pre-2026-10-03 behaviour).
                                          # L['eps'] stays the applied level everywhere (identity, x axes of the D_c and
                                          # kappa sweeps); L['eps_M'] / L['eps_M_src'] is what M, G and kappa divide by.
    DISP_TRIM_BINS: int = 3               # populated bins dropped at each face before the line fit (the diffuse
                                          # faces and the contact layers are not in the uniform-stress interior)
    DISP_MAX_FRAMES: int = 0              # plateau frames paired with the reference (0 = all of them)
    # ---- permeation vs compression consistency check (2026-10-03) ------------
    M_RECORD: object = None               # path of the compression stress-strain record the permeation notebook checks
                                          # against (save_M_record writes it from the sweep notebook); None -> the newest,
                                          # flow_data_local/compression/M_record_latest.json
    PERM_M_TOL: float = 0.10              # relative tolerance of the check (|u_F measured / predicted - 1|) before it flags
    # ---- D_c consolidation fit ------------------------------------------
    DC_N_MODES: int = 5
    DC_KMAX_IC: int = 199
    DC_FRAC_EARLY: float = 1.0
    DC_TRIM_BINS: int = 2
    DC_FREE_AMPS: bool = True
    DC_USE_FINE: bool = True              # read disp_z_polymer_fine (the deck's fine-cadence hold profile, DISP_FINE_NFREQ
                                          # in the batch, 2026-10-05) for the D_c fit when a level has one; the coarse
                                          # disp_z_polymer otherwise.  DC_TRIM_BINS stays in COARSE bins (binWidth) either way.
    DC_PLOT_MAX: int = 25                 # at most this many snapshots drawn in the D_c fit panels (the fit uses all)
    DC_WINDOW_AVG: bool = True            # fit_Dc averages each mode's decay over the ave/chunk window, as the permeation fit
                                          # does (PERM_DC_WINDOW_AVG): the deck block-averages u_z over the whole nfreq
                                          # interval (~1/20 of the hold) and tags the snapshot with the window END.  Read as
                                          # instantaneous, the snapshots returned 0.85x (window = 0.3 tau_1) to 0.74x
                                          # (0.6 tau_1) the true D_c (tests/dc_convention_test.py, 2026-10-03)
    DC_BOUNDS: tuple = (1e-6, 1.0)
    DC_SLOW_REF: float = 0.17             # the deck's hold-sizing constant Dc_est, in the deck's tau_1 = L^2/(pi^2 Dc_est)
                                          # formula.  The held slab relaxes as L^2/(4 pi^2 D_c), so the hold the deck
                                          # prescribes is that of D_c = 0.17/4 = 0.0425 (2026-09-28; the factor 4 is the
                                          # held slab's L/2 drainage path, not a second D_c convention -- 2026-10-03)
    DC_TARGET_RESID: float = 0.01
    # ---- D_c: fit window, stability scan, stress-trace D_c (2026-10-09) ----
    DC_FIT_TAU1: float = 6.0              # fit_Dc fits the FIRST DC_FIT_TAU1 x tau_1 of the hold, tau_1 = L^2/(4 pi^2 D_c) of the
                                          # fit itself (iterated to self-consistency), and keeps the whole-hold fit as F['Dc_all'].
                                          # 0 = fit the DC_FRAC_EARLY fraction of the hold (the pre-2026-10-09 behaviour).  Why:
                                          # the quarter gel's hold (comp_3 level 0.30, fine file) is a poroelastic transient at
                                          # D_c = 0.109 over its first ~6 tau_1 and then a SLOWER process (+0.012 extra interior
                                          # strain = 20 % of the transient, -8 % load-piston stress, from ~1000 tau); fitted whole
                                          # with poroelastic modes only, D_c halves to 0.054 and the fit becomes unstable in the
                                          # number of modes (N = 8: 0.021) and the trim.  Run 7's coarse file hides the same split
                                          # inside its first 2120-tau block (profile 0.053 vs its stress trace 0.073 at eps 0.10).
    DC_STABILITY: bool = True             # refit at 1, 2, 4, 10 tau_1 and the whole hold, and with 3 and 8 modes at the chosen
                                          # window -> F['stab'] (one line in load_level); |Dc_all/Dc - 1| > DC_STAB_FLAG is flagged
    DC_STAB_FLAG: float = 0.20
    DC_STRESS: bool = True                # D_c from the load-piston stress trace as well (fit_Dc_stress; piston_force_avg at
                                          # the volume_freq cadence, 25 tau): the held-slab series with a Gaussian-skin hold-onset
                                          # state from DC_STRESS_TMIN, and a single exponential on the tail from
                                          # DC_STRESS_TAIL_TAU1 x tau_1.  Independent of the profile binning, trim and window
                                          # averaging, and the only resolved transient of runs without a fine displacement file.
    DC_STRESS_TMIN: float = 250.0         # tau after the hold onset left out of the series fit (ramp ringing; modes faster than a block)
    DC_STRESS_TAIL_TAU1: float = 1.0      # the tail exponential starts this many tau_1 (of the series fit) into the hold
    # ---- unload / free re-swelling (2026-10-09; deck Phase 2c, -var unload 1) ----
    DC_UNLOAD: bool = True                # load the level's _u<lvl> files when they exist and fit the re-swelling D_c
                                          # (fit_Dc_unload): the skin-free control of the hold's D_c, same gel, same strain
    UNLOAD_L: str = 'hold'                # the L in tau_1 = L^2/(pi^2 D_c) of the unload fits: 'hold' (the compressed thickness,
                                          # the hold's own L -- so hold and unload D_c compare like for like), 'ref' (the swollen
                                          # reference thickness; D_c then reads (L_ref/L_hold)^2 larger) or 'mean'
    UNLOAD_TRACE_SKIP: float = 0.0        # tau after the retraction start left out of the thickness-trace fit (0 = none: the
                                          # faces detach from the retracting plates within a few tau)
    UNLOAD_N_MODES: int = 0               # modes of the unload profile fit (0 -> DC_N_MODES)
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
        assert self.PARTIAL_NORM in ('raw', 'share'), "PARTIAL_NORM must be 'raw' or 'share'"
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

    def upath(self, name, lvl, ext='dat'):
        """the level's UNLOAD file (deck Phase 2c, 2026-10-09): <name>_<DATANAME>_<INTERACTION>_<tag>_u<lvl>.<ext>,
        the tag being the level's hold tag."""
        t = self.tag_for(lvl)
        return self.DATA_DIR / f'{name}_{self.DATANAME}_{self.INTERACTION}_{t if t is not None else "*"}_u{lvl}.{ext}'

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


def relax_tail(steps, x, plateau_mean, dt=1.0, frac=0.6):
    """Unrelaxed excess of a plateau reading (2026-09-27).  Fits x(t) = x_inf + A exp(-t/tau)
    over the trailing `frac` of the series and returns dict(x_inf, A, tau, T_fit, excess, ok,
    why) with excess = plateau_mean - x_inf (>= 0).  The fit is credible (ok) only when it
    converged, relaxes DOWNWARD (A > 0, compression-positive), and tau is shorter than the
    fitted window (otherwise the tail is indistinguishable from a line and the asymptote is
    unconstrained); excess is 0 whenever not ok or the asymptote lies above the plateau."""
    steps = np.asarray(steps, float)
    x = np.asarray(x, float)
    out = dict(x_inf=np.nan, A=np.nan, tau=np.nan, T_fit=np.nan, excess=0.0, ok=False, why='')
    if len(steps) < 20:
        out['why'] = 'too few points'
        return out
    t = (steps - steps[0]) * dt
    sel = t >= (1.0 - frac) * t[-1]
    tw, xw = t[sel] - t[sel][0], x[sel]
    T = float(tw[-1])
    out['T_fit'] = T
    if sel.sum() < 20 or T <= 0:
        out['why'] = 'window too short'
        return out
    f = lambda tt, xi, A, tau: xi + A * np.exp(-tt / tau)
    try:
        p, _ = curve_fit(f, tw, xw, p0=[xw[-1], xw[0] - xw[-1], 0.3 * T], maxfev=20000)
    except Exception:
        out['why'] = 'fit did not converge'
        return out
    xi, A, tau = (float(v) for v in p)
    out.update(x_inf=xi, A=A, tau=tau)
    if not np.isfinite(tau) or tau <= 0 or tau > T:
        out['why'] = 'tau not shorter than the window (asymptote unconstrained)'
        return out
    if A <= 0:
        out['why'] = 'tail not relaxing downward'
        return out
    out['ok'] = True
    out['excess'] = float(max(plateau_mean - xi, 0.0))
    return out


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


def reservoir_mask(z, cfg, side, z_pist, z_gel, relax=True, gap_bins=None):
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
    gap = (cfg.res_gel_gap_bins if gap_bins is None else gap_bins) * cfg.binWidth
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
    of the reference and the level's own (scalar z_gel_hi), or the snapshot's OWN gel top
    (z_gel_hi an array, one per snapshot: a permeation run, where the feed face and the
    feed piston both come down below the reference gel top -- with the larger-of rule the
    window was empty and fell back to the reference window, by then in the vacuum above
    the piston: perm_3 read a feed baseline of 0; 2026-10-03).  One-piston runs tile the
    fixed window.  A snapshot whose window would be empty falls back to the reference window."""
    z = np.asarray(z, float)
    ts = np.asarray(ts, float)
    if np.isfinite(R.get('z_feed', np.nan)):
        if z_gel_hi is None:
            zg = np.full(len(ts), R['z_gel_hi'])
        elif np.ndim(z_gel_hi) == 0:
            zg = np.full(len(ts), max(float(z_gel_hi), R['z_gel_hi']))
        else:
            zg = np.asarray(z_gel_hi, float)
        zf = (lambda t: R['z_feed']) if z_feed_at is None else z_feed_at
        rows = [reservoir_mask(z, cfg, 'feed', zf(t), zg[i]) for i, t in enumerate(ts)]
        if z_gel_hi is not None and np.ndim(z_gel_hi) > 0:
            # a nearly drained feed reservoir (late in a permeation run): close the gel-edge gap bin by bin
            for i, t in enumerate(ts):
                for gb in range(cfg.res_gel_gap_bins - 1, -1, -1):
                    if rows[i].any():
                        break
                    rows[i] = reservoir_mask(z, cfg, 'feed', zf(t), zg[i], gap_bins=gb)
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
    R = {'flip_z': bool(cfg.FLIP_Z),       # orientation of the z-profile figures (zn / finish_axes / _zlim; plotting only)
         'flip_ends': (r'feed $\rightarrow$ permeate' if cfg.mode == 'permeation' else r'load piston $\rightarrow$ support')}
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
    ff = cfg.path('disp_z_polymer_fine', lvl)
    fine = bool(cfg.DC_USE_FINE and ff.exists())
    if fine:
        f = ff
    if not f.exists():
        return None
    snaps = read_ave_chunk_file(f)
    if not snaps:
        return None
    d = dict(fine=fine, ts=np.array([s[0] for s in snaps], float), z=snaps[0][1][:, 1],
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


# ---------------------------------------------------------------------------
#  Steady displacement profile -> the network strain (2026-10-03)
# ---------------------------------------------------------------------------
# Under a uniform network stress sigma' (the drained plateau of a compression hold) the
# polymer displacement from the eps = 0 reference is linear in the REFERENCE position,
#   u_z(Z) = u_0 - eps (Z - Z^P),   eps = sigma'/M,
# so the slope of u_z(Z) is the network strain itself, whatever rigid travel (the gel
# falling the ~3 sigma onto the support, the piston closing on the diffuse face) the
# applied strain counts.  The deck's disp_z_polymer is reset at every hold onset (the D_c
# fit wants the relaxation alone), so the steady profile is built from the per-atom
# trajectories: the plateau frames of traj_stress paired by atom id with the frames of
# traj_ref, binned by the reference z (Lagrangian).  Runs since 2026-10-03 also write
# disp_z_polymer_cum (displace/atom never reset after the reference, binned by CURRENT z),
# the fallback when the trajectories are not local: there du/dz = -eps/(1 - eps).
def traj_timesteps(traj):
    """The timesteps a LAMMPS dump holds (one sequential pass, no parsing)."""
    ts = []
    with open(traj) as f:
        line = f.readline()
        while line:
            if line.startswith('ITEM: TIMESTEP'):
                ts.append(int(f.readline()))
                f.readline(); n = int(f.readline())
                for _ in range(5 + n):
                    f.readline()
            line = f.readline()
    return ts


def read_traj_id_z(traj, want_ts=None, types=psd.POLYMER_TYPES):
    """{ts: (ids, z)} of the atoms of `types`, sorted by id, for the requested frames
    (all frames when want_ts is None).  pandas-parsed per frame (~0.3 s per 300k atoms)."""
    import pandas as pd
    want = None if want_ts is None else {int(t) for t in np.atleast_1d(want_ts)}
    out = {}
    with open(traj) as f:
        line = f.readline()
        while line:
            if line.startswith('ITEM: TIMESTEP'):
                ts = int(f.readline()); f.readline(); n = int(f.readline()); f.readline()
                for _ in range(3):
                    f.readline()
                cols = f.readline().split()[2:]
                if want is not None and ts not in want:
                    for _ in range(n):
                        f.readline()
                else:
                    df = pd.read_csv(f, sep=r'\s+', header=None, names=cols, nrows=n, engine='c')
                    df = df[df['type'].isin(types)].sort_values('id')
                    out[ts] = (df['id'].to_numpy(), df['z'].to_numpy(float))
            line = f.readline()
    return out


def load_ref_positions(cfg, R):
    """Mean reference z of every polymer atom over the traj_ref frames -> dict(ids, z, ts, sd)
    cached on R['zref']; None without the trajectory."""
    if 'zref' in R:
        return R['zref']
    tr = cfg.traj('traj_ref')
    if not tr.exists():
        R['zref'] = None
        return None
    fr = read_traj_id_z(tr)
    if not fr:
        R['zref'] = None
        return None
    ts = sorted(fr)
    ids = fr[ts[0]][0]
    Z = np.array([fr[t][1] for t in ts if np.array_equal(fr[t][0], ids)])
    R['zref'] = dict(ids=ids, z=Z.mean(axis=0), sd=Z.std(axis=0), ts=ts)
    return R['zref']


def _fit_disp_line(cfg, zc, um, n_atoms, coord, per_frame=None):
    """Line through the interior of a binned displacement profile.  zc: bin centres (reference z
    for 'lagrangian', current z for 'eulerian'); um: mean displacement per bin; n_atoms: atoms per
    bin per frame; per_frame: optional (n_frames, n_bins) stack for the frame-to-frame CI."""
    idx = np.where(n_atoms >= cfg.Ncount_min)[0]
    trim = int(max(cfg.DISP_TRIM_BINS, 0))
    inner = idx[trim:len(idx) - trim] if len(idx) > 2 * trim + 3 else idx
    if len(inner) < 3:
        return None
    coef, cov = np.polyfit(zc[inner], um[inner], 1, cov=True)
    slope, icpt = float(coef[0]), float(coef[1])
    fit = slope * zc + icpt
    resid = um - fit
    ss = np.sum((um[inner] - um[inner].mean()) ** 2)
    R2 = float(1.0 - np.sum(resid[inner] ** 2) / ss) if ss > 0 else np.nan
    to_eps = (lambda s: -s) if coord == 'lagrangian' else (lambda s: -s / (1.0 - s))
    eps = float(to_eps(slope))
    tq = stats.t.ppf(0.5 + cfg.ci_level / 2, max(len(inner) - 2, 1))
    half_lsq = float(tq * np.sqrt(max(cov[0, 0], 0.0)) * abs(to_eps(slope + 1e-9) - to_eps(slope - 1e-9)) / 2e-9)
    pf_eps, half_frames = [], 0.0
    if per_frame is not None and len(per_frame) >= 2:
        for row in per_frame:
            c = np.polyfit(zc[inner], row[inner], 1)
            pf_eps.append(float(to_eps(c[0])))
        pf = np.array(pf_eps)
        half_frames = float(stats.t.ppf(0.5 + cfg.ci_level / 2, len(pf) - 1) * pf.std(ddof=1) / np.sqrt(len(pf)))
    half = max(half_lsq, half_frames)
    return dict(zc=zc, u=um, n_atoms=n_atoms, idx=idx, inner=inner, slope=slope, icpt=icpt, fit=fit, resid=resid,
                R2=R2, eps=eps, eps_lo=eps - half, eps_hi=eps + half, eps_half_lsq=half_lsq,
                eps_half_frames=half_frames, per_frame_eps=pf_eps, coord=coord,
                u_bot=float(um[idx[0]]), u_top=float(um[idx[-1]]))


def fit_disp_profile(cfg, R, L):
    """Steady displacement profile of one compression level and the network strain it gives
    (see the section comment).  Trajectories first (plateau frames of traj_stress vs traj_ref,
    per atom, binned by the reference z); else disp_z_polymer_cum (plateau mean, current z).
    Returns a dict (keys of _fit_disp_line + src, ts, n_frames, z_P, z_T) or None."""
    tp = cfg.traj('traj_stress', L['lvl'])
    ref = load_ref_positions(cfg, R) if tp.exists() else None
    if ref is not None:
        ts_all = traj_timesteps(tp)
        plat = [t for t in ts_all if t >= L['halt_ts']] or ts_all[-1:]
        if cfg.DISP_MAX_FRAMES > 0 and len(plat) > cfg.DISP_MAX_FRAMES:
            plat = list(subsample(np.array(plat), np.array(plat)[:, None], cfg.DISP_MAX_FRAMES)[0])
        fr = read_traj_id_z(tp, plat)
        U, Z = [], None
        for t in sorted(fr):
            ids, z = fr[t]
            if np.array_equal(ids, ref['ids']):
                Zt, u = ref['z'], z - ref['z']
            else:
                _, ia, ib = np.intersect1d(ids, ref['ids'], return_indices=True)
                Zt, u = ref['z'][ib], z[ia] - ref['z'][ib]
            if Z is None or len(Zt) == len(Z):
                Z = Zt
                U.append(u)
        if U:
            U = np.array(U)
            edges = np.arange(np.floor(Z.min()), Z.max() + cfg.binWidth, cfg.binWidth)
            k = np.clip(np.digitize(Z, edges) - 1, 0, len(edges) - 2)
            nb = len(edges) - 1
            cnt = np.bincount(k, minlength=nb)
            per = np.array([np.bincount(k, weights=u, minlength=nb) / np.maximum(cnt, 1) for u in U])
            out = _fit_disp_line(cfg, 0.5 * (edges[:-1] + edges[1:]), per.mean(axis=0), cnt.astype(float), 'lagrangian', per)
            if out is not None:
                out.update(src='trajectories', ts=[int(t) for t in sorted(fr)], n_frames=len(U),
                           ref_ts=ref['ts'], z_P=float(R['z_support']), z_T=float(R['z_piston']))
                return out
    fc = cfg.path('disp_z_polymer_cum', L['lvl'])
    if fc.exists():
        snaps = read_ave_chunk_file(fc)
        if snaps:
            ts = np.array([s[0] for s in snaps], float)
            pl = ts >= L['halt_ts']
            if not pl.any():
                pl[-1] = True
            z = snaps[0][1][:, 1]
            Nc = np.array([s[1][:, 2] for s in snaps])[pl]
            uz = np.array([s[1][:, 3] for s in snaps])[pl]
            w = Nc.sum(axis=0)
            um = np.where(w > 0, (uz * Nc).sum(axis=0) / np.maximum(w, 1), np.nan)
            out = _fit_disp_line(cfg, z, np.nan_to_num(um), Nc.mean(axis=0), 'eulerian', uz)
            if out is not None:
                out.update(src='disp_z_polymer_cum', ts=[int(t) for t in ts[pl]], n_frames=int(pl.sum()),
                           z_P=float(L.get('z_supp', R['z_support'])), z_T=float(L.get('z_pist', R['z_piston'])))
                return out
    return None


def _strain_for_M(cfg, L):
    """(eps, source, CI half-width) that M, G and kappa divide by, per cfg.M_STRAIN with its
    fallback chain disp -> rg -> applied."""
    pref = cfg.M_STRAIN
    chain = {'disp': ('disp', 'rg', 'applied'), 'rg': ('rg', 'applied'), 'applied': ('applied',)}.get(pref)
    if chain is None:
        raise ValueError(f"M_STRAIN must be 'disp', 'rg' or 'applied' (got {pref!r})")
    for src in chain:
        if src == 'disp' and np.isfinite(L.get('eps_disp', np.nan)):
            return float(L['eps_disp']), 'disp', float(L['eps_disp_hi'] - L['eps_disp_lo']) / 2
        if src == 'rg' and np.isfinite(L.get('eps_rg', np.nan)):
            return float(L['eps_rg']), 'rg', 0.0
        if src == 'applied':
            return float(L['eps']), 'applied', 0.0


EPS_LABEL = {'disp': r'network strain  $\varepsilon_{\rm disp}$  (steady displacement-profile slope)',
             'rg': r'measured strain  $\varepsilon_{Rg}$', 'applied': r'applied strain  $\varepsilon$'}


def eps_axis_label(levels):
    """x-axis label for figures plotted against L['eps_M']."""
    srcs = {L.get('eps_M_src', 'applied') for L in levels}
    return EPS_LABEL[srcs.pop()] if len(srcs) == 1 else r'strain  $\varepsilon$  (mixed sources: ' + ', '.join(sorted(srcs)) + ')'


# ---------------------------------------------------------------------------
#  D_c: one equation, two boundary-value problems (2026-10-03)
#
#  Both consolidation fits -- fit_Dc here (compression hold) and fit_perm_Dc (permeation,
#  section 12) -- solve  du_z/dt = q(t) + D_c d2u_z/dz2,  D_c = kappa M, with L the FULL
#  thickness between the two faces (never a half-thickness).  The strain eps = -du_z/dz
#  diffuses with D_c in both; only its boundary conditions differ:
#    held slab   (both faces pinned, same bath pressure): eps(0) = eps(L), eps'(0) = eps'(L)
#                -> cos/sin(2 m pi z/L) at 4 m^2 pi^2 D_c/L^2,  tau_1 = L^2/(4 pi^2 D_c)
#    permeation  (support carries dP, feed face free):    eps(0) = dP/M,  eps(L) = 0
#                -> sin(k pi z/L)       at   k^2 pi^2 D_c/L^2,  tau_1 = L^2/(  pi^2 D_c)
#  The 4 between the two tau_1 is the drainage path (L/2 to the nearer plate against L
#  across the membrane), so equal D_c means the hold relaxes 4x FASTER than the permeation
#  transient of the same slab.  A prescribed-FLUX membrane (eps'(0) fixed instead of eps(0))
#  would have tau_1 = 4 L^2/(pi^2 D_c); the pistons prescribe the pressure, not the flux.
#  tests/dc_convention_test.py integrates the equation by finite differences with each
#  problem's physical boundary conditions (no mode series) and runs BOTH fitters on the
#  result: each returns the solver's D_c within 4 % (the residue is the Eulerian binning).
#  So a difference between the two fitted D_c on real runs is not a convention factor.
# ---------------------------------------------------------------------------
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
            return float(np.sum((Xd @ np.linalg.lstsq(Xd, y, rcond=None)[0] - y) ** 2))

        # coarse log grid first (short windows have secondary minima), then bounded Brent in the best cell
        grid = np.geomspace(cfg.DC_BOUNDS[0], cfg.DC_BOUNDS[1], 61)
        j = int(np.argmin([resid(g) for g in grid]))
        Dc = float(minimize_scalar(resid, bounds=(grid[max(j - 1, 0)], grid[min(j + 1, len(grid) - 1)]), method='bounded').x)
        Xd = X(Dc)
        A = np.linalg.lstsq(Xd, y, rcond=None)[0]
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


def _decay_mat(lam, t, W):
    """window_decay for an array of times: (len(t), len(lam))."""
    lam = np.asarray(lam, float)[None, :]
    t = np.atleast_1d(np.asarray(t, float))[:, None]
    w = np.minimum(W, t)
    inst = np.exp(-lam * t)
    with np.errstate(divide='ignore', invalid='ignore'):
        avg = (np.exp(-lam * (t - w)) - inst) / (lam * w)
    return np.where(w > 0, avg, inst)


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


# ---------------------------------------------------------------------------
#  Unload / free re-swelling (deck Phase 2c, -var unload 1; 2026-10-09): the skin-free control
#
#  After a level's hold the plates go back to the seated gap and the gel re-swells with BOTH
#  faces free: zero network stress there = zero displacement slope = DIRICHLET strain,
#  modes sin(k pi zeta) at k^2 pi^2 D_c/L^2, tau_1 = L^2/(pi^2 D_c) -- FOUR times the hold's
#  tau_1 at the same D_c (the hold only redistributes solvent between the ramp's compressed
#  skins and the interior with no flux through the plates; here solvent has to come in from
#  the reservoirs and fill the half thickness).  The initial state is the UNIFORM equilibrated
#  compression eps_0 = -dL_inf/L, so the solution is known in full:
#      u(zeta,t) - u(zeta,0) = sum_k B_k cos(k pi zeta) [exp(-lam_k t) - 1]       (free B_k in the fit)
#      dL(t)/dL_inf          = 1 - sum_{k odd} 8/(k pi)^2 exp(-k^2 pi^2 D_c t/L^2)   (no free amplitude)
#  and with dL_inf = L_ref - L_hold known from the reference, the thickness trace is a
#  ONE-parameter fit for D_c.  The faces detach from the retracting plates within a few tau
#  (the free face moves as 2 eps_0 sqrt(D t/pi), faster than the plates at first), so t = 0 of
#  the trace is the retraction start; the displacement reference is reset when the plates stop.
# ---------------------------------------------------------------------------
def load_unload(cfg, R, lvl, F_hold):
    """The level's _u<lvl> files -> dict, or None when there is no unload (no displacement file).
    F_hold: the level's fit_Dc result (its L, face planes and t_hold are the unload's reference)."""
    if F_hold is None:
        return None
    ff, fc = cfg.upath('disp_z_polymer_fine', lvl), cfg.upath('disp_z_polymer', lvl)
    fine = bool(cfg.DC_USE_FINE and ff.exists())
    f = ff if fine else fc
    if not f.exists():
        return None
    snaps = read_ave_chunk_file(f)
    if len(snaps) < 2:
        return None
    U = dict(fine=fine, ts=np.array([s[0] for s in snaps], float), z=snaps[0][1][:, 1],
             Nc=np.array([s[1][:, 2] for s in snaps]), uz=np.array([s[1][:, 3] for s in snaps]))
    fp = cfg.upath('piston_position', lvl)
    if not fp.exists():
        return None
    pp = np.atleast_2d(np.loadtxt(fp, comments='#'))
    U['t_retract'] = float(pp[0, 0])                                   # the loggers start with the retraction
    U['z_pist_free'] = float(pp[-1, 1])
    U['t_free'] = float(pp[np.argmax(np.isclose(pp[:, 1], pp[-1, 1], atol=1e-6)), 0])   # plates stopped: the reset
    fs = cfg.upath('support_position', lvl)
    if fs.exists():
        sp = np.atleast_2d(np.loadtxt(fs, comments='#'))
        U['z_supp_free'] = float(sp[-1, 1])
    fb = cfg.upath('gel_dimensions_bb', lvl)
    if fb.exists():
        bb = np.atleast_2d(np.loadtxt(fb, comments='#'))
        U['bb_step'], U['bb_L'] = bb[:, 0], bb[:, 3]
    fcm = cfg.upath('polymer_com', lvl)
    if fcm.exists():
        cm = np.atleast_2d(np.loadtxt(fcm, comments='#'))
        U['com_step'], U['com_z'] = cm[:, 0], cm[:, 3]
    fpa = cfg.upath('piston_force_avg', lvl)
    if fpa.exists():
        pfa = read_print_file(fpa, ['step', 'Fz'])
        U['pfa_step'], U['pfa_F'] = pfa['step'], pfa['Fz']
    # reference (swollen) thickness: the level's own bb file starts at the seated reference state (the
    # unload design compresses every level from the reference), else the reference profile edges
    fl = cfg.path('gel_dimensions_bb', lvl)
    if fl.exists():
        bl = np.atleast_2d(np.loadtxt(fl, comments='#'))
        U['L_ref'], U['L_ref_src'] = float(np.median(bl[:3, 3])), "level bb file, first rows (seated reference)"
    else:
        U['L_ref'], U['L_ref_src'] = float(R['z_gel_hi'] - R['z_gel_lo'] + cfg.binWidth), 'reference stress-profile edges'
    U['L_hold'] = float(F_hold['L'])
    U['Z_bot'], U['Z_top'] = float(F_hold['z_perm']), float(F_hold['z_feed'])     # face planes at the end of the hold
    return U


def _unload_series(t, Dc, L, n=60):
    """dL(t)/dL_inf of the free re-swelling from a uniform strain: 1 - sum_odd 8/(k pi)^2 exp(-k^2 pi^2 Dc t/L^2)."""
    k = 2.0 * np.arange(1, n + 1) - 1.0
    t = np.atleast_1d(np.asarray(t, float))[:, None]
    return 1.0 - (8.0 / (k * np.pi) ** 2 * np.exp(-(k * np.pi / L) ** 2 * Dc * np.clip(t, 0.0, None))).sum(axis=1)


def fit_Dc_unload(cfg, R, U, F_hold):
    """D_c of the free re-swelling (see the section comment): (a) the thickness trace L_bb(t) against the
    exact uniform-IC series, dL_inf FIXED at L_ref - L_hold (one parameter) and free (two); (b) the
    displacement profiles since the reset in the Lagrangian frame of the compressed gel, modes
    cos(k pi zeta) with free amplitudes, the COM drift removed, the same self-consistent window and
    stability scan as fit_Dc but with tau_1 = L^2/(pi^2 D_c).  Returns a dict or None."""
    if U is None:
        return None
    Lh, Lr = U['L_hold'], U['L_ref']
    L = {'hold': Lh, 'ref': Lr, 'mean': 0.5 * (Lh + Lr)}.get(cfg.UNLOAD_L, Lh)
    dL_inf0 = Lr - Lh
    out = dict(L=L, L_hold=Lh, L_ref=Lr, L_ref_src=U['L_ref_src'], dL_inf_ref=dL_inf0, t_retract=U['t_retract'],
               t_free=U['t_free'], retract_T=(U['t_free'] - U['t_retract']) * cfg.dt_lj, fine=U['fine'], trace=None, prof=None)
    # ---- (a) the thickness trace ----
    if 'bb_L' in U and dL_inf0 > 0:
        tt = (U['bb_step'] - U['t_retract']) * cfg.dt_lj
        dL = U['bb_L'] - Lh
        w = tt >= cfg.UNLOAD_TRACE_SKIP
        if w.sum() >= 20:
            T = dict(t=tt, dL=dL, fit=w)
            try:
                p1, c1 = curve_fit(lambda t_, Dc: dL_inf0 * _unload_series(t_, Dc, L), tt[w], dL[w],
                                   p0=(cfg.DC_SLOW_REF / 4.0,), bounds=((1e-5,), (10.0,)), maxfev=20000)
                r1 = dL[w] - dL_inf0 * _unload_series(tt[w], p1[0], L)
                T.update(Dc=float(p1[0]), Dc_se=float(np.sqrt(c1[0, 0]) * np.sqrt(max(autocorr_time(r1), 1.0))),
                         rms=float(np.sqrt(np.mean(r1 ** 2))), tau1=L ** 2 / (np.pi ** 2 * float(p1[0])),
                         R2=float(1.0 - np.sum(r1 ** 2) / max(np.sum((dL[w] - dL[w].mean()) ** 2), 1e-30)),
                         model=lambda t_, D=float(p1[0]): dL_inf0 * _unload_series(t_, D, L))
                p2, c2 = curve_fit(lambda t_, A, Dc: A * _unload_series(t_, Dc, L), tt[w], dL[w],
                                   p0=(dL_inf0, float(p1[0])), bounds=((0.0, 1e-5), (np.inf, 10.0)), maxfev=20000)
                r2 = dL[w] - p2[0] * _unload_series(tt[w], p2[1], L)
                s2 = np.sqrt(np.abs(np.diag(c2))) * np.sqrt(max(autocorr_time(r2), 1.0))
                T.update(Dc_free=float(p2[1]), Dc_free_se=float(s2[1]), dL_inf_free=float(p2[0]), dL_inf_free_se=float(s2[0]),
                         model_free=lambda t_, A=float(p2[0]), D=float(p2[1]): A * _unload_series(t_, D, L),
                         relaxed=float(_unload_series(tt[-1], float(p1[0]), L)[0]))
                out['trace'] = T
            except (RuntimeError, ValueError) as e:
                print(f'  NOTE: unload thickness-trace fit failed ({e})')
    # ---- (b) the profiles: Lagrangian zeta of the compressed gel, COM drift removed ----
    ts, z, Nc, uz = U['ts'], U['z'], U['Nc'], U['uz']
    t_lj = (ts - U['t_free']) * cfg.dt_lj
    keep = t_lj > 0
    ts, Nc, uz, t_lj = ts[keep], Nc[keep], uz[keep], t_lj[keep]
    if len(ts) >= 2:
        if 'com_z' in U:
            com = np.interp(ts, U['com_step'], U['com_z']) - np.interp(U['t_free'], U['com_step'], U['com_z'])
        else:
            com = np.zeros(len(ts))
        W = float(np.min(np.diff(ts))) * cfg.dt_lj if (cfg.DC_WINDOW_AVG and len(ts) > 1) else 0.0
        n_trim = int(round(cfg.DC_TRIM_BINS * cfg.binWidth / abs(float(z[1] - z[0])))) if len(z) > 1 else cfg.DC_TRIM_BINS
        N = int(cfg.UNLOAD_N_MODES) or int(cfg.DC_N_MODES)
        k = np.arange(1, N + 1, dtype=float)
        zl, yl = [], []
        for i in range(len(ts)):
            p_ = np.where(Nc[i] > cfg.Ncount_min)[0]
            if n_trim and len(p_) > 2 * n_trim + 4:
                p_ = p_[n_trim:-n_trim]
            zi = (z[p_] - uz[i][p_] - U['Z_bot']) / Lh                 # material coordinate of the bin's atoms
            yi = (uz[i][p_] - com[i]) / L
            ok = (zi > 0) & (zi < 1)
            zl.append(zi[ok]); yl.append(yi[ok])
        if min(len(a) for a in zl) >= 4:
            def _fit(sel, k_):
                y = np.concatenate([yl[i] for i in sel])
                Wm = [np.cos(np.pi * np.outer(zl[i], k_)) for i in sel]

                def X(Dc):
                    d = _decay_mat((np.pi * k_ / L) ** 2 * Dc, t_lj[sel], W) - 1.0
                    return np.vstack([Wm[n] * d[n][None, :] for n in range(len(sel))])

                def resid(Dc):
                    Xd = X(Dc)
                    return float(np.sum((Xd @ np.linalg.lstsq(Xd, y, rcond=None)[0] - y) ** 2))

                grid = np.geomspace(cfg.DC_BOUNDS[0], cfg.DC_BOUNDS[1], 61)
                j = int(np.argmin([resid(g) for g in grid]))
                Dc = float(minimize_scalar(resid, bounds=(grid[max(j - 1, 0)], grid[min(j + 1, 60)]), method='bounded').x)
                Xd = X(Dc)
                A = np.linalg.lstsq(Xd, y, rcond=None)[0]
                ss = np.sum((y - y.mean()) ** 2)
                return Dc, A, float(1.0 - np.sum((Xd @ A - y) ** 2) / ss) if ss > 1e-30 else np.nan

            all_ = np.arange(len(ts))

            def _window(Dc, n_tau):
                sel = all_[t_lj <= n_tau * L ** 2 / (np.pi ** 2 * Dc)]
                return all_[:3] if len(sel) < 3 else sel

            Dc_all, A_all, R2_all = _fit(all_, k)
            early, Dc, A, R2 = all_, Dc_all, A_all, R2_all
            n_tau = float(cfg.DC_FIT_TAU1)
            if n_tau > 0:
                seen = {len(all_)}
                for _ in range(12):
                    sel = _window(Dc, n_tau)
                    if len(sel) == len(early) and np.array_equal(sel, early):
                        break
                    if len(sel) in seen:
                        break
                    seen.add(len(sel))
                    early, (Dc, A, R2) = sel, _fit(sel, k)
            tau1 = L ** 2 / (np.pi ** 2 * Dc)
            stab = []
            if cfg.DC_STABILITY:
                for nt in (0.5, 1.0, 2.0, 4.0):
                    sel = _window(Dc, nt)
                    d, _, r = _fit(sel, k)
                    stab.append(dict(what=f'{nt:g} tau_1', n=int(len(sel)), Dc=d, R2=r))
                stab.append(dict(what='whole unload', n=int(len(all_)), Dc=Dc_all, R2=R2_all))
            out['prof'] = dict(Dc=Dc, A=A, R2=R2, Dc_all=Dc_all, R2_all=R2_all, early=early, n_all=int(len(all_)), tau1=tau1,
                               n_tau1_fit=float(t_lj[early[-1]] / tau1), stab=stab, t_lj=t_lj, ts=ts, zl=zl, yl=yl, k=k, W=W, com=com,
                               stab_flag=bool(n_tau > 0 and len(early) < len(all_) and abs(Dc_all / Dc - 1.0) > cfg.DC_STAB_FLAG),
                               u_model=lambda zh, t, A_=A, D_=Dc: np.cos(np.pi * np.outer(np.atleast_1d(zh), k)) @ (A_ * (_decay_mat((np.pi * k / L) ** 2 * D_, t, W)[0] - 1.0)),
                               shown=early[::max(1, int(np.ceil(len(early) / max(cfg.DC_PLOT_MAX, 1))))])
    if out['trace'] is None and out['prof'] is None:
        return None
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
    say(f"  strain: applied {L['eps']:.4f} | measured Rg {L.get('eps_rg', np.nan):.4f}  "
        f"BB {L.get('eps_bb', np.nan):.4f}  boundary(diag) {L.get('eps_boundary', np.nan):.4f}")
    # ---- steady displacement profile -> the network strain (2026-10-03) ----
    L['disp_prof'] = D = fit_disp_profile(cfg, R, L)
    if D is not None:
        L['eps_disp'], L['eps_disp_lo'], L['eps_disp_hi'] = D['eps'], D['eps_lo'], D['eps_hi']
        say(f"  steady displacement profile ({D['src']}, {D['n_frames']} plateau frame(s), {D['coord']} bins, "
            f"{len(D['inner'])} interior bins): slope strain eps_disp = {D['eps']:.4f} [{D['eps_lo']:.4f}, {D['eps_hi']:.4f}]"
            f"  R^2 = {D['R2']:.4f};  u at the faces {D['u_bot']:+.2f} / {D['u_top']:+.2f} sigma")
    L['eps_M'], L['eps_M_src'], L['eps_M_half'] = _strain_for_M(cfg, L)
    say(f"  M, G, kappa divide by eps_{L['eps_M_src']} = {L['eps_M']:.4f} (M_STRAIN = {cfg.M_STRAIN!r})"
        + ('' if L['eps_M_src'] == cfg.M_STRAIN else f"  -- FALLBACK: no eps_{cfg.M_STRAIN} for this level"))

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
    eps = L['eps_M']                                   # the network strain (M_STRAIN, 2026-10-03; was the applied level)
    eps_h = L['eps_M_half']                            # its CI half-width, folded into the M CIs below
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
    half = np.sqrt(((hi - lo) / 2) ** 2 + ref_h ** 2 + (m * eps_h / eps) ** 2)
    L['M_net'], L['M_net_lo'], L['M_net_hi'] = float(m), float(m - half), float(m + half)
    L['M_net_nbins'] = len(mn_abs)
    L['dsig_net'] = float(m * eps)                     # the stress increment itself (strain-definition free)
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
            L['M_pist'] = (p['mean'] - R['P_ref']) / eps
            half = np.sqrt((np.sqrt(ph ** 2 + rh ** 2) / eps) ** 2 + (L['M_pist'] * eps_h / eps) ** 2)
            L['M_pist_lo'], L['M_pist_hi'] = L['M_pist'] - half, L['M_pist'] + half
        else:
            L['M_pist'], L['M_pist_lo'], L['M_pist_hi'] = L['M_pist_abs'], L['M_pist_abs_lo'], L['M_pist_abs_hi']
        L['dP_pist'] = float(L['M_pist'] * eps)
    # ---- unrelaxed-hold systematic (2026-09-27, RELAX_SYS): the piston tail's excess over its
    #      fitted asymptote, / eps, widens the LOWER bound of BOTH M estimators (same relaxation
    #      drives both); kappa = D_c/M below inherits it on its upper bound.  The block-bootstrap
    #      CIs see only the scatter inside the plateau window, not the drift that is still in it.
    L['M_sys'] = 0.0
    if cfg.RELAX_SYS and 'PF' in L and 'pfa_P' in L:
        rt = relax_tail(L['pfa_step'], L['pfa_P'], L['PF']['mean'], cfg.dt_lj, cfg.RELAX_TAIL_FRAC)
        L['relax'] = rt
        L['M_sys'] = rt['excess'] / eps
        L['M_net_lo'] -= L['M_sys']
        if 'M_pist' in L:
            L['M_pist_lo'] -= L['M_sys']
        if rt['ok']:
            say(f"  unrelaxed-hold systematic: piston tail -> P_inf = {rt['x_inf']:.4f} (tau ~ {rt['tau']:.0f} tau over the last "
                f"{cfg.RELAX_TAIL_FRAC:.0%}), plateau excess {rt['excess']:+.4f} = {rt['excess'] / max(L['PF']['mean'], 1e-30):+.1%}"
                f"  ->  delta_sys = {L['M_sys']:.4f} added to the LOWER bound of M_network and M_piston")
        else:
            say(f"  unrelaxed-hold systematic: 0 ({rt['why']})")
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
        # unrelaxed-hold systematic (RELAX_SYS, 2026-09-27): the unrelaxed excess is undrained pore
        # pressure, i.e. ISOTROPIC -- on run 7 sigma'_xx and sigma'_yy drift by the same amount as
        # sigma'_zz (0.142 vs 0.143 at eps = 0.5) -- so it cancels in sigma'_zz - sigma'_xx and G
        # carries only the anisotropic residual.  That residual is measured directly: the drop of the
        # interior anisotropy over the second half of the hold (3-snapshot means), / 2 eps, >= 0.
        G_sys = 0.0
        if cfg.RELAX_SYS and len(ts) >= 8:
            an = np.array([np.nanmean((dzz[i] - dxx[i])[im]) for i in range(len(ts))])
            h = len(ts) // 2
            G_sys = float(max(np.nanmean(an[h:h + 3]) - np.nanmean(an[-3:]), 0.0) / (2.0 * eps))
        Glo -= G_sys
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
        L['G'][comp] = dict(G=Gm, lo=Glo, hi=Ghi, sys=G_sys, lam=lam, ratio=ratio, ratio_err=ratio_err,
                            nbins=len(g_bins), ratio_final=ratio_plat, ratio_final_err=ratio_plat_err)
        say(f"  G from {comp}: {Gm:.4f} [{Glo:.4f}, {Ghi:.4f}]" + (f" (lower bound incl. anisotropy drift {G_sys:.4f})" if G_sys > 0 else '')
            + f"   lambda_{comp} = {lam:.4f}   "
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
            f"DL/L = {F['DL'] / F['L']:.4f} [support share {F['f_sup']:.2f}], hold = {F['hold_T']:.0f} tau = {F['hold_T'] / F['tau1']:.1f} tau_1; "
            f"fitted: first {F['n_tau1_fit']:.1f} tau_1 = {len(F['early'])} of {len(F['early_all'])} {'FINE' if F['fine'] else 'coarse'} snapshots, "
            f"first at {F['t_lj'][F['early'][0]]:.0f} tau, block {F['W']:.0f} tau)")
        if F['fit_tau1'] > 0 and len(F['early']) < len(F['early_all']):
            say(f"      whole hold: D_c = {F['Dc_all']:.4e} (R^2 = {F['R2_all']:.3f})"
                + (f"   <-- {F['Dc'] / F['Dc_all']:.2f}x the windowed value: a slower process after the poroelastic transient "
                   f"(fig_Dc panel c); kappa below uses the windowed D_c" if F['stab_flag'] else '   (consistent)'))
        if F.get('stab'):
            say('      stability: ' + '   '.join(f"{s['what']} ({s['n']} snaps): {s['Dc']:.3e}" for s in F['stab']))
        L['kappa'] = {}
        for key, Mk in (('net', 'M_net'), ('pist', 'M_pist')):
            if Mk in L:
                L['kappa'][key] = dict(k=F['Dc'] / L[Mk], lo=F['Dc'] / L[Mk + '_hi'], hi=F['Dc'] / L[Mk + '_lo'])
        say('  kappa = D_c/M:  ' + '   '.join(f"{k}: {v['k']:.4e} [{v['lo']:.4e}, {v['hi']:.4e}]"
                                             for k, v in L['kappa'].items()))
        # ---- the stress trace: an independent D_c (2026-10-09) ----
        L['Dc_stress'] = None
        if cfg.DC_STRESS and 'pfa_P' in L:
            S = fit_Dc_stress(cfg, L['pfa_step'], L['pfa_P'], F['L'], F['t_hold'])
            L['Dc_stress'] = S
            if S['ok']:
                say(f"  D_c (stress trace) = {S['Dc']:.4e} +/- {S['Dc_se']:.1e}  (held-slab series from {S['t_min']:.0f} tau, "
                    f"{S['n']} blocks of {S['block']:.0f} tau, skin {S['delta']:.1f} sigma, rms {S['rms']:.4f}; tau_1 = {S['tau1']:.0f} tau)"
                    + (f";  tail exponential from {S['t_tail']:.0f} tau: tau = {S['tau_tail']:.0f} +/- {S['tau_tail_se']:.0f} "
                       f"-> D_c = {S['Dc_tail']:.4e}" if S['tail_ok'] else ';  tail exponential not credible')
                    + f"   [profile/stress = {F['Dc'] / S['Dc']:.2f}]")
                L['kappa_stress'] = {key: S['Dc'] / L[Mk] for key, Mk in (('net', 'M_net'), ('pist', 'M_pist')) if Mk in L}
                say('      kappa from the stress-trace D_c:  ' + '   '.join(f"{k}: {v:.4e}" for k, v in L['kappa_stress'].items()))
            else:
                say(f"  D_c (stress trace): not resolved -- {S['why']}")
        # ---- the unload / free re-swelling (deck Phase 2c), the skin-free control (2026-10-09) ----
        L['unload'] = None
        if cfg.DC_UNLOAD:
            Uf = fit_Dc_unload(cfg, R, load_unload(cfg, R, lvl, F), F)
            L['unload'] = Uf
            if Uf is not None:
                T, Pf = Uf['trace'], Uf['prof']
                say(f"  UNLOAD (free re-swelling after the hold; retraction {Uf['retract_T']:.0f} tau; L = {Uf['L']:.2f} = {cfg.UNLOAD_L!r} "
                    f"[hold {Uf['L_hold']:.2f}, reference {Uf['L_ref']:.2f}: {Uf['L_ref_src']}]; tau_1 = L^2/(pi^2 D_c), faces free)")
                if T is not None:
                    say(f"    D_c (thickness trace, dL_inf fixed = {Uf['dL_inf_ref']:.2f} sigma) = {T['Dc']:.4e} +/- {T['Dc_se']:.1e}  "
                        f"(R^2 {T['R2']:.3f}, rms {T['rms']:.3f} sigma; tau_1 = {T['tau1']:.0f} tau, {T['t'][-1] / T['tau1']:.1f} tau_1 recorded, "
                        f"{100 * T['relaxed']:.1f}% relaxed);  free asymptote: D_c = {T['Dc_free']:.4e} +/- {T['Dc_free_se']:.1e}, "
                        f"dL_inf = {T['dL_inf_free']:.2f} +/- {T['dL_inf_free_se']:.2f}")
                if Pf is not None:
                    say(f"    D_c (profiles, first {Pf['n_tau1_fit']:.1f} tau_1 = {len(Pf['early'])} of {Pf['n_all']} {'FINE' if Uf['fine'] else 'coarse'} snapshots) = "
                        f"{Pf['Dc']:.4e} (R^2 {Pf['R2']:.3f});  whole unload {Pf['Dc_all']:.4e}"
                        + ('   <-- differ' if Pf['stab_flag'] else '') + ('   stability: ' + '  '.join(f"{x['what']} ({x['n']}): {x['Dc']:.3e}" for x in Pf['stab']) if Pf['stab'] else ''))
                ref = T['Dc'] if T is not None else Pf['Dc']
                say(f"    vs the hold:  windowed profile {F['Dc'] / ref:.2f}x,  whole hold {F['Dc_all'] / ref:.2f}x"
                    + (f",  stress trace {L['Dc_stress']['Dc'] / ref:.2f}x" if (L.get('Dc_stress') or {}).get('ok') else '')
                    + "  the unload trace D_c  (1.0 = the hold's transient is the same D_c; < 1 = the hold reads slower, e.g. skin-throttled)")
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
        S['D_reg'][name] = tuple(float(np.ravel(v)[0]) for v in mean_ci(Dm[:, None], cfg.ci_level))   # 1-element arrays (numpy >= 2)
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
    finish_axes(ax, r'geometric porosity $\epsilon_g$', f'(a) porosity ($r_{{\\rm probe}}={cfg.PSD_R_PROBE}$) vs $\\phi_s^{{\\rm cal}}$', R=R)
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
    finish_axes(ax, r'$\langle D_{\rm pore}\rangle$  ($\sigma$)', '(b) mean pore diameter per bin', R=R)
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
             'sigmayy_polymer', 'sigmayy_solvent', 'solvent_density_z', 'disp_z_polymer', 'disp_z_polymer_cum',
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
# fine-cadence hold displacement profile (DISP_FINE_NFREQ > 0 in the batch, 2026-10-05): optional, most runs have none
_TWO_PIST_OPT = ('disp_z_polymer_fine',)
# unload / free re-swelling files (deck Phase 2c, -var unload 1, 2026-10-09): optional, tagged _<tag>_u<lvl>
_UNLOAD_DAT = ('disp_z_polymer', 'disp_z_polymer_fine', 'gel_dimensions_bb', 'gel_dimensions_rg', 'piston_position',
               'support_position', 'piston_force_avg', 'piston_pressure', 'permeation', 'strain_piston', 'polymer_com',
               'gel_edges', 'pressure_reservoirs', 'sigmazz_polymer', 'sigmazz_solvent')
_PERM_DAT = ('sigmazz_polymer', 'sigmazz_solvent', 'sigmaxx_polymer', 'sigmaxx_solvent',
             'sigmayy_polymer', 'sigmayy_solvent', 'solvent_density_z', 'disp_z_polymer', 'strain_zz',
             'piston_position', 'piston_velocity', 'piston_force', 'piston_force_avg', 'piston_pressure',
             'permeation', 'permeate_count', 'pressure_feed', 'pressure_permeate', 'support_force',
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
        data += [cfg.path(n, l) for n in _PROD_DAT + ((_TWO_PIST_DAT + _TWO_PIST_OPT) if cfg.two_pist else ())]
        if cfg.two_pist:
            data += [cfg.upath(n, l) for n in _UNLOAD_DAT]
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
    optional |= {cfg.path(n, l).name for n in _TWO_PIST_DAT + _TWO_PIST_OPT for l in ([None] + list(cfg.COMP_LEVELS))}
    optional |= {cfg.upath(n, l).name for n in _UNLOAD_DAT for l in cfg.COMP_LEVELS}
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
    """Fractional box height z/Lz -- or 1 - z/Lz under cfg.FLIP_Z (R['flip_z'], set by load_reference,
    2026-10-07), so that the feed / load piston end (high z) draws at 0 on the left and the permeate /
    support end at 1 on the right.  Plotting coordinate only: every mask, window and fit uses z itself."""
    x = (np.asarray(z, float) - R['Z_LO']) / R['LZ']
    return 1.0 - x if R.get('flip_z') else x


def _flipped(R):
    return bool(R is not None and R.get('flip_z'))


def _zx_label(R):
    """x label of a z/L profile axis: z/L, or 1 - z/L with the two ends named under FLIP_Z."""
    return rf'$1 - z/L$   ({R["flip_ends"]})' if _flipped(R) else r'$z/L$'


def flip_hint(R, label):
    """Append 'left -> right' orientation to the x label of a sigma- or zeta-unit axis that is
    drawn inverted under FLIP_Z (flip_z_axis / _zlim); the label itself when not flipped."""
    return f'{label}   [{R["flip_ends"]}, left to right]' if _flipped(R) else label


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


def _zlim(R, lo=0.0, hi=1.0):
    """x-limits of a zeta-unit profile axis (zeta = 0 at the support).  Under cfg.FLIP_Z (R['flip_z'],
    2026-10-06) the axis is inverted -- (hi, lo) -- so the feed / load piston end is on the left; the
    data keep their coordinates and the label says so (flip_hint).  R = None or a dict without the
    key (other libraries) -> (lo, hi)."""
    return (hi, lo) if _flipped(R) else (lo, hi)


def flip_z_axis(ax, R):
    """FLIP_Z for a z-profile axis with autoscaled limits (z or Z in sigma): invert it BEFORE
    plotting, so the autoscaled limits stay inverted and the legend / box placement sees the
    final layout; label it with flip_hint.  No-op without R['flip_z']."""
    if _flipped(R):
        ax.xaxis.set_inverted(True)


def finish_axes(ax, ylabel, title, R=None):
    """x label, limits and grid of a z/L profile axis (zn coordinates: 0 -> 1 either way)."""
    ax.set_xlabel(_zx_label(R))
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
    finish_axes(ax, ylabel, title, R=R)
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
    finish_axes(ax, ylabel, title, R=R)
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
    finish_axes(ax, ylabel, title, R=R)
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
    finish_axes(ax, r'$\phi_s$', title, R=R)
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
    """1 x 3 evolution panels (zz, xx, yy) for kind='t' (total) or 'net' (network).  With
    cfg.PARTIAL_NORM = 'share' the total-stress figure is followed by (sigma^t - P_ref)/dP
    (_fig_total_norm, saved as <stem>_norm; 2026-10-06)."""
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
            finish_axes(ax, lab, title, R=R)
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
    fig = _save(fig, cfg, stem, L['lvl'])
    if kind == 't' and cfg.PARTIAL_NORM == 'share':        # normalised companion, drawn right under the original
        _fig_total_norm(cfg, R, L, stem + '_norm')
    return fig


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
        finish_axes(ax, yl, t, R=R)
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
    finish_axes(ax, ylabel, title, R=R)
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
            finish_axes(ax, r"$\sigma'_{%s}$" % comp, title, R=R)
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
    final bold) with the eps = 0 reference dashed.  Needs the xx and yy profiles.
    With cfg.PARTIAL_NORM = 'share' a second figure follows (thermo_pressure_evolution_norm,
    2026-10-06): the solvent, polymer and total traces with the back pressure removed, / dP
    (partial_norm on the trace; the total keeps its gradient across the gel)."""
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
    h, lab = ax.get_legend_handles_labels()
    if _is_perm(L):
        # PERMEATION (2026-10-02): the profile is not a plateau in general (it slopes as G grows), so no
        # interior mean is quoted and the bold final curve is just 'final'; the y-axis starts at the bath
        # pressure P_BARO so the P_th(z) structure above it is readable.
        ax.set_ylim(cfg.P_BARO, max(hi, cfg.P_BARO + 0.3 * (hi - lo)))
        lab = ['final' if s == 'final (plateau)' else s for s in lab]
        note = f'dotted: $P_{{\\rm bath}}={sig(cfg.P_BARO)}$'
    else:
        ax.set_ylim(min(lo, cfg.P_BARO - 0.3 * (hi - lo)), max(hi, cfg.P_BARO + 0.3 * (hi - lo)))
        note = f'plateau mean in gel interior = {fmt_mu(np.nanmean(Pt[L["plat"]], axis=0)[L["interior"]])}   (dotted: $P_{{\\rm bath}}={sig(cfg.P_BARO)}$)'
    h.append(Patch(alpha=0, label=note))
    lab.append(note)
    smart_legend(ax, handles=h, labels=lab, fontsize=12)
    lab_eps = f'($\\varepsilon={L["eps"]:.2f}$)' if np.isfinite(L['eps']) else '(permeation drive)'
    fig.suptitle(f'Thermodynamic pressure evolution {lab_eps}  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    fig = _save(fig, cfg, 'thermo_pressure_evolution', L['lvl'])
    if cfg.PARTIAL_NORM == 'share':
        _fig_partial_norm(cfg, R, L, 'tr', 'thermo_pressure_evolution_norm')
    return fig


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


def _partial_evo(ax, cfg, R, L, fam, stack_of, ref_of, zero=True):
    """Solvent / polymer / total profile evolutions on one axis: one colour family per entry
    of fam = ((key, label, cmap name), ...), the reference dashed, the final curve bold.
    stack_of(key) -> the per-snapshot stack, ref_of(key) -> the reference mean profile (None:
    not drawn).  Returns the legend handles and every curve drawn (for the y-range)."""
    zx = zn(R, R['z'])
    handles, curves = [], []
    for key, lab, cmap_name in fam:
        cmap = plt.get_cmap(cmap_name)
        ts, ev = post_halt(cfg, L, L['ts'], stack_of(key))
        ref = ref_of(key)
        base = cmap(0.85)
        if ref is not None:
            ax.plot(zx, ref, '--', color=base, lw=1.8, alpha=0.9, zorder=2)
            curves.append(ref)
        n = len(ts)
        for i in range(n):
            last = (i == n - 1)
            ax.plot(zx, ev[i], '-', color=(base if last else cmap(0.3 + 0.5 * i / max(n - 1, 1))),
                    lw=(3.2 if last else 1.2), alpha=(1.0 if last else 0.6), zorder=(5 if last else 3))
        curves.extend(ev)
        handles.append(Line2D([0], [0], color=base, lw=3, label=lab))
    handles.append(Line2D([0], [0], color='0.4', ls='--', lw=2, label=r'reference ($\varepsilon=0$)'))
    if not _is_perm(L):
        handles.append(Line2D([0], [0], color='0.4', lw=1.2, alpha=0.6, label='hold (faint = early)'))
    if zero:
        ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    shade_gel(ax, R, L)
    mark_walls(ax, R, L)
    return handles, curves


def fig_partial_stress(cfg, R, L):
    """One axis: solvent partial (blues), polymer partial (oranges) and total
    (greys) sigma_zz evolutions, reference dashed, final curves bold.  With
    cfg.PARTIAL_NORM = 'share' a second figure follows (partial_stress_evolution_norm,
    2026-10-06): the same three curves with the back pressure removed, / dP (partial_norm)."""
    zz, Rz = L['stress']['zz'], R['stress']['zz']
    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    fam = (('s', r'solvent $\sigma_{s,zz}$', 'Blues'), ('p', r'polymer $\sigma_{p,zz}$', 'Oranges'),
           ('t', r'total $\sigma^{t}_{zz}$', 'Greys'))
    handles, _ = _partial_evo(ax, cfg, R, L, fam, lambda k: zz[k], lambda k: Rz[k + '_m'])
    finish_axes(ax, r'$\sigma_{zz}(z,t)$  (LJ)',
                'Partial and total normal stresses: reference $\\rightarrow$ evolution' if _is_perm(L) else
                f'Partial and total $\\sigma_{{zz}}$: reference $\\rightarrow$ compressed ($\\varepsilon={L["eps"]:.2f}$)', R=R)
    smart_legend(ax, handles=handles, fontsize=12)
    fig = _save(fig, cfg, 'partial_stress_evolution', L['lvl'])
    if cfg.PARTIAL_NORM == 'share':
        _fig_partial_norm(cfg, R, L, 'zz', 'partial_stress_evolution_norm')
        if _is_perm(L):
            _fig_partial_norm_phi(cfg, R, L, 'zz', 'partial_stress_evolution_norm_phi')   # TEST (2026-10-07), see below
    return fig


# ---------------------------------------------------------------------------
#  Back-pressure-subtracted partial stresses (cfg.PARTIAL_NORM = 'share', 2026-10-06)
# ---------------------------------------------------------------------------
# Marioni et al. plot the raw partial P_zz of water and polyamide across the membrane with
# the permeate at ~0, so their partial / total IS each species' share of the driving pressure.
# Here the permeate (or the bath) sits at P_ref = 1.5 and the signal (dP ~ 0.1) rides on top
# of it, so the same figure needs each species' share of that back pressure taken off first.
# The share comes from the stresses themselves, per z-bin and per snapshot -- no mass or
# volume fractions, no calibration against another run:
#     w_s      = sigma_s / sigma_t                      (same snapshot, same component)
#     sigma_s* = (sigma_s - w_s P_ref) / dP             = w_s sigma_t*
#     sigma_p* = (sigma_p - (1 - w_s) P_ref) / dP       = (1 - w_s) sigma_t*
#     sigma_t* = (sigma_t - P_ref) / dP
# i.e. the excess of the total over the back pressure, split between the species in proportion
# to their current partial stresses; sigma_s* + sigma_p* = sigma_t* by construction (asserted).
# In a reservoir sigma_p = 0, so w_s = 1: the solvent curve reads 1 in the feed and 0 in the
# permeate, the polymer curve 0.  This is the literal analogue of Marioni's figure but it does
# NOT reproduce his end points (water 0.3 -> 0, polyamide 0.7 -> 1): the solvent's share of
# the back pressure cannot be told apart from the transferred load with the local stresses
# alone, so the solvent curve does not reach 0 at the permeate face.  That is accepted.
#     P_ref : permeation -> the applied permeate pressure (_p_norm); compression -> cfg.P_BARO
#     dP    : permeation -> the applied feed - permeate pressure; compression -> the level's
#             plateau load-piston pressure increment (L['dP_pist'])
# The trace (thermodynamic pressure) is treated the same way with w_s = tr(sigma_s)/tr(sigma_t).
# Only the constant P_ref is taken off the total, so the total trace keeps the gradient it has
# across the gel.  Bins that hold no matter (beyond the pistons, sigma_t ~ 0) are masked, and
# so are the bins within res_wall_margin of a wet-piston plane (the same clearance as the
# reservoir windows: the sheet depletes and layers the solvent there and takes half the wall
# virial, an offset of ~0.07 that the division by dP would turn into a spike of order 1).
_NORM_VAC = 0.05          # a bin is empty when |sigma_t| <= _NORM_VAC * P_ref


def _clear_of_pistons(cfg, R, L, n):
    """(n_snap, n_bins) mask of the z-bins lying ENTIRELY between the two wet-piston planes
    and >= cfg.res_wall_margin clear of them, per stress snapshot of level / run L (measured
    piston tracks; the extreme plane over the snapshot's averaging window, since the pistons
    travel) or of the reference (L = None: the planes at rest).  All True for one-piston runs."""
    z, h = np.asarray(R['z'], float), 0.5 * cfg.binWidth
    m = np.ones((n, len(z)), bool)
    wz = L.get('wetz') if L is not None else None
    for i in range(n):
        if wz:
            t1 = float(L['ts'][i])
            t0 = float(L['ts'][i - 1]) if i > 0 else t1 - (float(L['ts'][1] - L['ts'][0]) if n > 1 else 0.0)
            zf = min(wz['feed_at'](t0), wz['feed_at'](t1))
            zp = max(wz['perm_at'](t0), wz['perm_at'](t1))
        else:
            zf, zp = R.get('z_feed', np.nan), R.get('z_perm', np.nan)
        if np.isfinite(zf):
            m[i] &= z + h <= zf - cfg.res_wall_margin
        if np.isfinite(zp):
            m[i] &= z - h >= zp + cfg.res_wall_margin
    return m


def _dp_norm(cfg, L):
    """Driving pressure the normalised profiles are divided by (see the block comment above);
    nan when the run does not provide it."""
    if _is_perm(L):
        if cfg.DP_PISTON is not None:
            return float(cfg.DP_PISTON)
        pl = (L.get('wet') or {}).get('plat', {})
        if 'P_feed_app' in pl and 'P_perm_app' in pl:
            return float(pl['P_feed_app'] - pl['P_perm_app'])
        return np.nan
    return float(L.get('dP_pist', np.nan))


def _norm_symbols(L):
    """(P_ref, dP) as mathtext, for labels."""
    return (r'P_{\rm perm}', r'\Delta P_{\rm ext}') if _is_perm(L) else (r'P_{\rm bath}', r'\Delta P_{\rm pist}')


def _norm_skip(cfg, L, what):
    print(f'normalised {what} figure skipped (PARTIAL_NORM = {cfg.PARTIAL_NORM!r}): no driving pressure -- '
          + ('piston_pressure carries no applied feed / permeate pressure' if _is_perm(L) else 'no load-piston plateau (piston_force file)'))


def _spt(D, comp):
    """(solvent, polymer, total) per-snapshot stacks of D = L['stress'] or R['stress'] for
    comp in COMPONENTS, or comp = 'tr' for the trace / 3; None when a file is missing."""
    if comp == 'tr':
        return tuple(_tr3(D, k) for k in ('s', 'p', 't')) if all(c in D for c in COMPONENTS) else None
    S = D.get(comp)
    return None if S is None else (S['s'], S['p'], S['t'])


def partial_norm(cfg, R, L, comp='zz', ref=False):
    """Back-pressure-subtracted partial stresses of level / run L (see the block comment
    above): dict(s, p, t = the normalised solvent, polymer and total stacks (snapshot x z-bin),
    w = the solvent share sigma_s/sigma_t, P_ref, dP).  comp = 'zz' | 'xx' | 'yy' | 'tr' (trace).
    ref=True normalises the eps = 0 reference stacks of R with the SAME P_ref and dP: its total
    equals P_ref, so all three come out ~0 -- the sanity check drawn dashed.  None when the
    component or the driving pressure is missing."""
    P_ref, dP = _p_norm(cfg, L), _dp_norm(cfg, L)
    spt = _spt(R['stress'] if ref else L['stress'], comp)
    if spt is None or not (np.isfinite(dP) and dP > 0):
        return None
    s, p, t = (np.asarray(a, float) for a in spt)
    ok = (np.abs(t) > _NORM_VAC * P_ref) & _clear_of_pistons(cfg, R, None if ref else L, len(t))
    with np.errstate(invalid='ignore', divide='ignore'):
        w = np.where(ok, s / t, np.nan)                                  # nan = masked bin, in all three
    out = dict(w=w, P_ref=P_ref, dP=dP, s=(s - w * P_ref) / dP, p=(p - (1.0 - w) * P_ref) / dP,
               t=np.where(np.isfinite(w), (t - P_ref) / dP, np.nan))
    assert np.allclose(out['s'] + out['p'], out['t'], rtol=0, atol=1e-9, equal_nan=True), 'partial_norm: s* + p* != t*'
    return out


def _norm_ylim(ax, curves, mask, pad=0.2):
    """y-range of a normalised panel: robust over the scale bins, always showing 0 and 1."""
    robust_ylim(ax, curves, zmask=mask, pad=pad, qlo=1, qhi=99)
    lo, hi = ax.get_ylim()
    ax.set_ylim(min(lo, -0.15), max(hi, 1.15))


def _fig_total_norm(cfg, R, L, stem):
    """Normalised companion of fig_total_stress: (sigma^t_c - P_ref)/dP for c = zz, xx, yy on
    the same evolution panels (reference dashed with its 95 % band, final bold), dotted guide
    at 1.  Under permeation the zz panel reads 1 in the feed and across the membrane, 0 in the
    permeate; in a compression hold 1 in the gel and 0 in both reservoirs."""
    dP = _dp_norm(cfg, L)
    if not (np.isfinite(dP) and dP > 0):
        _norm_skip(cfg, L, 'total-stress')
        return None
    pr, dp = _norm_symbols(L)
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    fig.suptitle(f'Total stress with the back pressure removed, $(\\sigma^t-{pr})/{dp}$  ($' + pr + f'={sig(_p_norm(cfg, L))}$, $'
                 + dp + f'={sig(dP)}$)  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    for ax, comp in zip(axes, COMPONENTS):
        lab = rf'$(\sigma^{{t}}_{{{comp}}}-{pr})/{dp}$'
        title = f'({"abc"[COMPONENTS.index(comp)]}) ' + r'total $\sigma^{t}_{%s}$, back pressure removed' % comp
        N = partial_norm(cfg, R, L, comp)
        if N is None:
            ax.text(0.5, 0.5, f'sigma{comp} files\nnot found', ha='center', va='center', transform=ax.transAxes)
            finish_axes(ax, lab, title, R=R)
            continue
        Nr = partial_norm(cfg, R, L, comp, ref=True)
        ref = mean_ci(Nr['t'], cfg.ci_level) if Nr is not None else None
        ts, ev = post_halt(cfg, L, L['ts'], N['t'])
        ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6, zorder=1)
        plot_evolution(ax, cfg, R, L, R['z'], ts, ev, lab + '$(z,t)$', title, ref=ref, annotate=False, legend=False)
        _norm_ylim(ax, list(ev) + ([ref[0]] if ref is not None else []), _scale_mask(cfg, R, L), pad=0.3)
        note = (('steady' if _is_perm(L) else 'plateau') + ' mean in gel interior = '
                + fmt_mu(mean_ci(N['t'][L['plat']], cfg.ci_level)[0][L['interior']]))
        h, lb = ax.get_legend_handles_labels()
        h.append(Patch(alpha=0, label=note))
        lb.append(note)
        smart_legend(ax, handles=h, labels=lb, fontsize=12)
    return _save(fig, cfg, stem, L['lvl'])


def _norm_family(comp):
    """((key, label, cmap name), ...) and the y symbol of the normalised partial figures."""
    if comp == 'tr':
        return ((('s', r'solvent $P^*_{th,s}$', 'Blues'), ('p', r'polymer $P^*_{th,p}$', 'Oranges'),
                 ('t', r'total $P^*_{th}$', 'Greys')), r'P_{th}')
    return ((('s', r'solvent $\sigma^*_{s,%s}$' % comp, 'Blues'), ('p', r'polymer $\sigma^*_{p,%s}$' % comp, 'Oranges'),
             ('t', r'total $\sigma^{t*}_{%s}$' % comp, 'Greys')), r'\sigma_{%s}' % comp)


def _gm_guide(cfg, L):
    """(M, source of M, [(G, name, colour), ...]) for the dP_th guides of the normalised
    thermodynamic-pressure figure (2026-10-07): G = cfg.G_REF (the shear notebooks' modulus) and
    cfg.G_COMP_REF (the compression-mode estimate), whichever are set; M = cfg.M_REF when set,
    else the run's own -- the primary permeation estimate (add_perm_displacement) or the level's
    M_net.  None when no G or no M is available."""
    Gs = [(float(G), name, col) for G, name, col in ((cfg.G_REF, r'G_{\rm shear}', WONG['green']),
                                                     (cfg.G_COMP_REF, r'G_{\rm comp}', WONG['reddishpurple']))
          if G is not None and np.isfinite(G)]
    if not Gs:
        return None
    if cfg.M_REF is not None:
        return float(cfg.M_REF), 'M_REF', Gs
    if _is_perm(L):
        Mp = L.get('M_perm')
        if Mp and Mp.get('primary'):
            return float(Mp[Mp['primary']]['M']), f"$M_{{\\rm {Mp['primary']}}}$", Gs
        return None
    return (float(L['M_net']), r'$M_{\rm net}$', Gs) if np.isfinite(L.get('M_net', np.nan)) else None


def _draw_dpth_guide(ax, cfg, R, L, N):
    """The dP_th = (4/3)(G/M) dP guides on the normalised trace figure (2026-10-07), one per G
    estimate (G_REF: shear notebooks; G_COMP_REF: the triaxial holds' lateral stress).  Under
    uniaxial strain the lateral stresses carry (M - 2G)/M of the axial one, so the trace of the
    load is (1 - 4G/3M) of it: across a membrane that passes dP_ext the total P_th falls from the
    feed reservoir's value at the feed face to dP_th = (4/3)(G/M) dP_ext below it at the permeate
    side (drawn as the straight line between the membrane faces, 1 -> 1 - 4G/3M in normalised
    units); in a compression hold it sits at 1 - 4G/3M across the gel.  Returns the legend handles."""
    gm = _gm_guide(cfg, L)
    if gm is None:
        return []
    M, src, Gs = gm
    pr, dp = _norm_symbols(L)
    handles = []
    for G, name, col in Gs:
        r = 4.0 * G / (3.0 * M)
        if _is_perm(L):
            x0, x1 = zn(R, L['z_mem_hi']), zn(R, L['z_mem_lo'])            # feed face -> support face
            h, = ax.plot([x0, x1], [1.0, 1.0 - r], '-.', color=col, lw=2.4, zorder=6,
                         label=(rf'$\Delta P_{{th}}=\frac{{4}}{{3}}\frac{{{name}}}{{M}}\,{dp}$ = {r:.2f} ${dp}$  (${name}$ = {G:g}, {src} = {M:.3f}):' + '\n'
                                rf'feed-reservoir $P_{{th}}$ at the feed face $\rightarrow$ {1.0 - r:.2f} at the permeate side'))
        else:
            x0, x1 = zn(R, L['z_mem_lo']), zn(R, L['z_mem_hi'])
            h, = ax.plot([x0, x1], [1.0 - r, 1.0 - r], '-.', color=col, lw=2.4, zorder=6,
                         label=rf'$1-\frac{{4}}{{3}}\frac{{{name}}}{{M}}$ = {1.0 - r:.2f}  (${name}$ = {G:g}, {src} = {M:.3f})')
        handles.append(h)
    return handles


def fig_thermo_pressure_norm(cfg, R, L):
    """The normalised thermodynamic-pressure figure on its own (2026-10-07; fig_thermo_pressure
    draws it under the raw one when cfg.PARTIAL_NORM = 'share'): solvent, polymer and total
    traces with the back pressure removed, / dP, plus the dP_th = (4/3)(G/M) dP guides when
    cfg.G_REF / cfg.G_COMP_REF are set."""
    return _fig_partial_norm(cfg, R, L, 'tr', 'thermo_pressure_evolution_norm')


def _fig_partial_norm(cfg, R, L, comp, stem):
    """Normalised companion of fig_partial_stress (comp = 'zz') and of fig_thermo_pressure
    (comp = 'tr', the partial and total traces): solvent* (blues), polymer* (oranges) and total*
    (greys) evolutions of partial_norm on one axis, the reference dashed (~0: the sanity check),
    final curves bold, dotted guides at 0 and 1.  The trace figure also carries the
    dP_th = (4/3)(G/M) dP guide (_draw_dpth_guide) when cfg.G_REF is set."""
    what = 'thermodynamic-pressure' if comp == 'tr' else 'partial-stress'
    N = partial_norm(cfg, R, L, comp)
    if N is None:
        if np.isfinite(_dp_norm(cfg, L)):
            print(f'normalised {what} figure skipped (sigmaxx / sigmayy files missing)')
        else:
            _norm_skip(cfg, L, what)
        return None
    Nr = partial_norm(cfg, R, L, comp, ref=True)
    pr, dp = _norm_symbols(L)
    fam, sym = _norm_family(comp)
    fig, ax = plt.subplots(figsize=((16, 7) if comp == 'tr' and _gm_guide(cfg, L) else (13, 7)), constrained_layout=True)
    for y in (0.0, 1.0):
        ax.axhline(y, color='k', ls=':', lw=1.2, alpha=0.6, zorder=1)
    handles, curves = _partial_evo(ax, cfg, R, L, fam, lambda k: N[k],
                                   lambda k: (None if Nr is None else mean_ci(Nr[k], cfg.ci_level)[0]), zero=False)
    name = (r'Partial and total $P_{th}=-\frac{1}{3}\,\mathrm{tr}(\mathbf{\sigma})$' if comp == 'tr'
            else r'Partial and total $\sigma_{%s}$' % comp)
    finish_axes(ax, rf'$({sym}-w\,{pr})/{dp}$', name + ', back pressure removed: reference $\\rightarrow$ '
                + ('evolution' if _is_perm(L) else f'compressed ($\\varepsilon={L["eps"]:.2f}$)'), R=R)
    _norm_ylim(ax, curves, _scale_mask(cfg, R, L))
    if comp == 'tr':
        hs = _draw_dpth_guide(ax, cfg, R, L, N)
        handles.extend(hs)
        if not hs and cfg.G_REF is None and cfg.G_COMP_REF is None:
            print('  (no dP_th = 4/3 G/M dP guide: set cfg.G_REF / cfg.G_COMP_REF, the shear modulus estimates)')
    smart_legend(ax, handles=handles, fontsize=12)
    num = rf'\mathrm{{tr}}\,\sigma_s/\mathrm{{tr}}\,\sigma^t' if comp == 'tr' else rf'\sigma_{{s,{comp}}}/\sigma^t_{{{comp}}}'
    fig.suptitle(f'Each species\' share $w$ of the back pressure removed, per bin and snapshot:  $w_s={num}$, $w_p=1-w_s$, $w=1$ for the total\n'
                 f'(${pr}={sig(N["P_ref"])}$, ${dp}={sig(N["dP"])}$)  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    return _save(fig, cfg, stem, L['lvl'])


def _phi_cal_on_snapshots(L, ts):
    """lambda-calibrated solvent volume fraction phi_s^cal (snapshot x z-bin) on the stress
    snapshot times ts: the per-frame Voronoi stacks of L['vf'] interpolated linearly in time,
    bin by bin (held constant beyond the first / last tessellated frame); without the per-
    frame stacks the steady mean L['phi_cal'] is used for every snapshot.  None without a
    calibrated fraction."""
    vf = L.get('vf')
    if vf is not None and vf.get('phi_cal') is not None and len(vf['ts']):
        ft, fc = np.asarray(vf['ts'], float), np.asarray(vf['phi_cal'], float)
        o = np.argsort(ft)
        ft, fc = ft[o], fc[o]
        ts = np.asarray(ts, float)
        return np.column_stack([np.interp(ts, ft, fc[:, j]) for j in range(fc.shape[1])])
    if L.get('phi_cal') is not None:
        return np.broadcast_to(np.asarray(L['phi_cal'][0], float), (len(ts), len(L['phi_cal'][0]))).copy()
    return None


def _fig_partial_norm_phi(cfg, R, L, comp, stem):
    """TEST figure (2026-10-07, permeation only; may be dropped): the normalised total sigma^t*
    of _fig_partial_norm next to the normalised SOLVENT partial divided by the calibrated solvent
    volume fraction, sigma_s* / phi_s^cal -- the solvent's back-pressure-free excess per unit
    solvent volume (phi_s^cal = 1 in the reservoirs, so the two curves coincide there).  No
    polymer curve.  phi_s^cal per snapshot from _phi_cal_on_snapshots; the reference uses
    R['phi_cal'].  Bins with phi_s^cal <= cfg.PHI_FLOOR are masked."""
    N = partial_norm(cfg, R, L, comp)
    if N is None:
        return None
    phi = _phi_cal_on_snapshots(L, L['ts'])
    if phi is None:
        print('normalised partial-stress / phi_s^cal figure skipped (no calibrated volume fraction: run add_perm_volume_fractions with VOR_ENABLE)')
        return None
    Nr = partial_norm(cfg, R, L, comp, ref=True)
    phi_r = np.asarray(R['phi_cal'][0], float) if R.get('phi_cal') is not None else None
    with np.errstate(invalid='ignore', divide='ignore'):
        s_phi = np.where(phi > cfg.PHI_FLOOR, N['s'] / phi, np.nan)
        s_phi_r = (None if (Nr is None or phi_r is None)
                   else np.where(phi_r > cfg.PHI_FLOOR, mean_ci(Nr['s'], cfg.ci_level)[0] / phi_r, np.nan))
    pr, dp = _norm_symbols(L)
    fam = (('s', r'solvent $\sigma^*_{s,%s}/\phi_s^{\rm cal}$' % comp, 'Blues'),
           ('t', r'total $\sigma^{t*}_{%s}$' % comp, 'Greys'))
    stacks = {'s': s_phi, 't': N['t']}
    refs = {'s': s_phi_r, 't': (None if Nr is None else mean_ci(Nr['t'], cfg.ci_level)[0])}
    fig, ax = plt.subplots(figsize=(13, 7), constrained_layout=True)
    for y in (0.0, 1.0):
        ax.axhline(y, color='k', ls=':', lw=1.2, alpha=0.6, zorder=1)
    handles, curves = _partial_evo(ax, cfg, R, L, fam, lambda k: stacks[k], lambda k: refs[k], zero=False)
    finish_axes(ax, rf'$(\sigma_{{{comp}}}-w\,{pr})/{dp}$   [solvent: $/\,\phi_s^{{\rm cal}}$]',
                r'Total $\sigma^{t*}_{%s}$ and solvent $\sigma^*_{s,%s}$ per unit solvent volume: reference $\rightarrow$ evolution' % (comp, comp), R=R)
    _norm_ylim(ax, curves, _scale_mask(cfg, R, L))
    vf = L.get('vf')
    src = (rf'$\phi_s^{{\rm cal}}$ interpolated in time from {len(vf["ts"])} tessellated frames' if vf is not None and vf.get('phi_cal') is not None
           else r'$\phi_s^{\rm cal}$ = steady mean for every snapshot')
    note = 'steady mean in gel interior: solvent ' + fmt_mu(mean_ci(s_phi[L['plat']], cfg.ci_level)[0][L['interior']]) \
           + ', total ' + fmt_mu(mean_ci(N['t'][L['plat']], cfg.ci_level)[0][L['interior']])
    handles.append(Patch(alpha=0, label=note))
    smart_legend(ax, handles=handles, fontsize=12)
    fig.suptitle(f'TEST -- solvent partial with the back pressure removed, divided by the calibrated solvent volume fraction '
                 f'(P_CAL_MODE = {cfg.P_CAL_MODE!r}; {src})\n'
                 f'(${pr}={sig(N["P_ref"])}$, ${dp}={sig(N["dP"])}$)  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    return _save(fig, cfg, stem, L['lvl'])


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


def _unload_panels(cfg, axe, axf, F, Uf, S):
    """The unload row of fig_Dc (2026-10-09): (e) the thickness trace with the one-parameter series fit, the
    free-asymptote fit and the series drawn at the hold's D_c values; (f) the re-swelling profiles and the
    modal fit in the compressed gel's material coordinate."""
    T, Pf, L = Uf['trace'], Uf['prof'], Uf['L']
    if T is not None:
        tt, dL = T['t'], T['dL']
        td = np.geomspace(max(tt[tt > 0][0], 1.0), tt[-1], 400)
        axe.semilogx(tt[tt > 0], dL[tt > 0], '.', color='0.6', ms=2, label=r'$L_\mathrm{bb}(t) - L_\mathrm{hold}$ (every print)')
        axe.semilogx(td, T['model'](td), '-', color=WONG['vermillion'], lw=2.6,
                     label=rf"series, $\Delta L_\infty$ fixed = {Uf['dL_inf_ref']:.2f} $\sigma$: $D_c={sig(T['Dc'])}\pm{sig(T['Dc_se'], 1)}$, $\tau_1={T['tau1']:.0f}\,\tau$")
        axe.semilogx(td, T['model_free'](td), '--', color=WONG['blue'], lw=2.0,
                     label=rf"series, $\Delta L_\infty$ free = {T['dL_inf_free']:.2f}: $D_c={sig(T['Dc_free'])}$")
        for D_, lab, ls in ((F['Dc'], "hold, windowed profile", ':'), (F['Dc_all'], "hold, whole-hold profile", '-.')):
            axe.semilogx(td, Uf['dL_inf_ref'] * _unload_series(td, D_, L), ls, color='0.3', lw=1.6, label=rf"series at the {lab} $D_c={sig(D_)}$")
        if S is not None and S.get('ok'):
            axe.semilogx(td, Uf['dL_inf_ref'] * _unload_series(td, S['Dc'], L), ls='--', color=WONG['green'], lw=1.6,
                         label=rf"series at the hold's stress-trace $D_c={sig(S['Dc'])}$")
        axe.axhline(Uf['dL_inf_ref'], color='k', lw=0.8, alpha=0.5)
        axe.axvline(Uf['retract_T'], color='k', lw=0.8, ls=':', alpha=0.6)
        axe.text(Uf['retract_T'] * 1.1, 0.93 * Uf['dL_inf_ref'], 'plates stopped', fontsize=9, color='0.3')
        axe.set(xlabel=r'time since retraction start  ($\tau$)', ylabel=r'$\Delta L = L_\mathrm{bb}(t) - L_\mathrm{hold}$  ($\sigma$)')
        axe.xaxis.set_minor_formatter(NullFormatter())
        axe.set_title(rf"(e) unload re-swelling ($\tau_1=L^2/\pi^2 D_c$): trace $D_c={sig(T['Dc'])}$ vs hold {sig(F['Dc'])}, ratio {F['Dc'] / T['Dc']:.2f}", fontsize=14)
        axe.grid(alpha=0.3)
        axe.legend(loc='lower right', fontsize=9)
    else:
        axe.text(0.5, 0.5, 'no thickness trace (gel_dimensions_bb_*_u<lvl>)', ha='center', va='center', transform=axe.transAxes)
    if Pf is not None:
        zff = np.linspace(0.0, 1.0, 300)
        norm = Normalize(vmin=Pf['ts'][Pf['early'][0]], vmax=Pf['ts'][Pf['early'][-1]])
        cmap = plt.cm.viridis
        for i in Pf['shown']:
            c = cmap(norm(Pf['ts'][i]))
            axf.plot(Pf['zl'][i], Pf['yl'][i], 'o', color=c, ms=2.5, alpha=0.4)
            axf.plot(zff, Pf['u_model'](zff, Pf['t_lj'][i]), '-', color=c, lw=1.6)
        axf.set(xlabel=r'$\zeta$ (material coordinate of the compressed gel)',
                ylabel=r'$u_z/L$ since the plates stopped')
        axf.set_title(rf"(f) unload profiles ($\cos k\pi\zeta$, COM removed): $D_c={sig(Pf['Dc'])}$ over {Pf['n_tau1_fit']:.1f} $\tau_1$, whole {sig(Pf['Dc_all'])}",
                      fontsize=13, color=WONG['vermillion'] if Pf['stab_flag'] else 'k')
        axf.grid(alpha=0.3)
        axf._tri_has_colorbar = True
    else:
        axf.text(0.5, 0.5, 'no unload displacement file', ha='center', va='center', transform=axf.transAxes)


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
        fig, ((axl, axr), (axc, axd), (axe, axf)) = plt.subplots(3, 2, figsize=(18, 19), constrained_layout=True)
        _unload_panels(cfg, axe, axf, F, Uf, L.get('Dc_stress'))
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
        ax.set(xlabel=flip_hint(R, r'$\zeta=(z-z_\mathrm{perm})/L$'), ylabel=r'$u_z/L$', xlim=_zlim(R))
        ax.grid(alpha=0.3)
        smart_legend(ax, fontsize=11)
    axl.set_title(r'(a) $u_z(\zeta,t)/L$ -- fitted hold snapshots (data + fitted IC offset)', fontsize=15)
    lbl = '' if cfg.DC_FREE_AMPS else rf"$\beta=p_0/M={sig(F['beta'])}$, "
    axr.set_title(rf"(b) consolidation fit, first {F['n_tau1_fit']:.1f} $\tau_1$: $D_c={sig(F['Dc'])}\ \sigma^2/\tau$, " + lbl + rf"$R^2={sig(F['R2'])}$", fontsize=15)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=[axl, axr], fraction=0.015, pad=0.04).set_label('timestep')
    fig.suptitle(f'Cooperative diffusivity fit (level _c{L["lvl"]})  |  {cfg.sim_name}  |  '
                 f'plate-gap closure $\\Delta L={sig(F["DL"])}\\,\\sigma$: support {fs:.0%} / piston {1 - fs:.0%}'
                 + ('  (symmetric drive)' if abs(fs - 0.5) < 0.05 else '  (top-only drive)' if fs < 0.05 else ''),
                 fontsize=12, fontweight='bold')
    return _save(fig, cfg, 'Dc_consolidation_fit', L['lvl'])


def fig_v0_check(cfg, R, L, disp=None):
    """TEMPORARY diagnostic (2026-10-09; delete once the D_c discrepancy is resolved).  Does the barycentric velocity
    v0(t) = phi_s v_s + phi_p v_p = q_s + q_p bias the hold's D_c?  The displacement equation du/dt = v0(t) + D_c u''
    becomes the strain equation d eps/dt = D_c eps'' with v0 gone (uniform in z); the displacement fit (fit_Dc) keeps
    u and lets its pinned-face modes fix v0(t) -- the solvent shuttle of the antisymmetric block.  Here the same fit
    is run with the per-snapshot rigid translation projected out (= the strain-PDE fit) and the two D_c are compared
    over every fit window; (b) shows the mean displacement over the fit bins, whose rate IS v0(t) in the hold
    (d<u>/dt = v0 + D_c [eps(1) - eps(0)]/L with equal face strains), data against the displacement model."""
    F = L.get('Dc')
    if F is None:
        print('v0 check skipped (no D_c fit)')
        return None
    disp = load_disp(cfg, R, L['lvl']) if disp is None else disp
    F2 = fit_Dc(cfg, R, disp, free_offset=True)
    if F2 is None:
        print('v0 check skipped (strain fit failed)')
        return None
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(17, 6.5), constrained_layout=True)
    # (a) D_c over the fit windows, both fits
    labs = [x['what'] for x in F['stab'] if 'modes' not in x['what']] + [f"fit ({F['n_tau1_fit']:.1f} tau_1)"]
    v1 = [x['Dc'] for x in F['stab'] if 'modes' not in x['what']] + [F['Dc']]
    v2 = [x['Dc'] for x in F2['stab'] if 'modes' not in x['what']] + [F2['Dc']]
    xi = np.arange(len(labs))
    axa.plot(xi, v1, 'o-', color=WONG['blue'], ms=10, lw=2, label=r'displacement fit: $v^0(t)$ from the pinned-face modes (fit_Dc)')
    axa.plot(xi, v2, 's--', color=WONG['vermillion'], ms=9, lw=2, mfc='none', mew=2, label=r'strain-PDE fit: rigid translation projected out ($v^0$ drops out)')
    S = L.get('Dc_stress')
    if S is not None and S.get('ok'):
        axa.axhline(S['Dc'], color=WONG['green'], ls=':', lw=2, label=f"load-piston stress trace: {sig(S['Dc'])}")
    axa.set_xticks(xi)
    axa.set_xticklabels([l.replace(' tau_1', r' $\tau_1$') for l in labs], rotation=20, fontsize=11)
    axa.set_ylabel(r'$D_c$  ($\sigma^2/\tau$)')
    axa.set_title(rf"(a) level _c{L['lvl']}: $D_c$ with and without $v^0$ -- {sig(F['Dc'])} vs {sig(F2['Dc'])} (ratio {F2['Dc'] / F['Dc']:.3f})", fontsize=14)
    axa.grid(alpha=0.3)
    smart_legend(axa, fontsize=10)
    # (b) the rigid part: <u> over the fit bins, data vs the displacement model
    t, zf, Lh = F['t_lj'], F['zf'], F['L']
    um_data = F['u_mean'] * Lh
    um_model = np.array([F['u_tr'](zf, ti).mean() for ti in t]) * Lh
    axb.semilogx(t, um_data, '.', color='k', ms=4 if len(t) < 60 else 2.5, label=r'data: $\langle u_z\rangle$ over the fit bins')
    axb.semilogx(t, um_model, '-', color=WONG['blue'], lw=2.2, label='displacement model (its own rigid part = the modal solvent shuttle)')
    axb.semilogx(t, um_data - um_model, '-', color=WONG['vermillion'], lw=1.8, alpha=0.9, label=r'difference = translation the modes cannot represent, $\int v^0_\mathrm{extra}\,dt$')
    axb.axhline(0, color='k', lw=0.8)
    amp = float(np.nanmax(np.abs(F['uhat'][:, F['idx']]))) * Lh
    axb.set(xlabel=r'hold time  ($\tau$)', ylabel=r'mean displacement  ($\sigma$)')
    axb.xaxis.set_minor_formatter(NullFormatter())
    axb.set_title(rf"(b) rigid translation during the hold  (transient amplitude for scale: {amp:.2f} $\sigma$)", fontsize=14)
    axb.grid(alpha=0.3)
    smart_legend(axb, fontsize=10)
    fig.suptitle(f'TEMPORARY v0 check  |  {cfg.sim_name}  |  d u/dt = v0(t) + D_c u\'\' vs d eps/dt = D_c eps\'\'', fontsize=12, fontweight='bold')
    print(f"  v0 check (level {L['lvl']}): D_c displacement fit {F['Dc']:.4e} (whole hold {F['Dc_all']:.4e})  |  strain-PDE fit {F2['Dc']:.4e} "
          f"(whole hold {F2['Dc_all']:.4e})  ->  ratio {F2['Dc'] / F['Dc']:.3f};  rms rigid translation data - model = {np.sqrt(np.mean((um_data - um_model) ** 2)):.3f} sigma "
          f"(transient amplitude {amp:.2f} sigma)")
    return _save(fig, cfg, 'v0_check', L['lvl'])


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
        finish_axes(ax, r'$\phi_s$', title, R=R)
        ax.set_ylim(0, 1.15)
        smart_legend(ax, handles=level_handles(levels, ref=True), fontsize=11)
        if txt:
            annotate_box(ax, 'in-gel means\n' + '\n'.join(txt), loc='lower right', fontsize=11)
        else:
            ax.text(0.5, 0.5, 'unavailable', ha='center', va='center', transform=ax.transAxes)
    return _save(fig, cfg, 'sweep_volfrac_profiles')


def _sweep_stress_panels(cfg, R, levels, kind, stem, suptitle):
    """1 x 3 sweep overlays (zz, xx, yy) for kind='t' (total) or 'net' (network).  With
    cfg.PARTIAL_NORM = 'share' the total-stress figure is followed by its normalised version
    (_sweep_norm_panels, saved as <stem>_norm; 2026-10-06)."""
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
            finish_axes(ax, lab, title, R=R)
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
    fig = _save(fig, cfg, stem)
    if kind == 't' and cfg.PARTIAL_NORM == 'share':        # normalised companion, drawn right under the original
        pr, dp = _norm_symbols(None)
        _sweep_norm_panels(cfg, R, levels, [(c, 't', rf'$(\sigma^{{t}}_{{{c}}}-{pr})/{dp}$', f'({"abc"[i]}) '
                                             + r'total $\sigma^{t}_{%s}$, back pressure removed' % c) for i, c in enumerate(COMPONENTS)],
                           stem + '_norm', f'Total stress with the back pressure removed, $(\\sigma^t-{pr})/{dp}$, all levels', 'total-stress')
    return fig


def _norm_ref_level(cfg, levels):
    """The level whose driving pressure normalises the eps = 0 reference in the sweep overlays:
    the smallest one (the strictest reading of 'the reference comes out ~0')."""
    ok = [L for L in levels if np.isfinite(_dp_norm(cfg, L)) and _dp_norm(cfg, L) > 0]
    return min(ok, key=lambda L: _dp_norm(cfg, L)) if ok else None


def _sweep_norm_panels(cfg, R, levels, panels, stem, suptitle, what):
    """1 x 3 sweep overlays of partial_norm profiles (2026-10-06): panels = ((comp, key, ylabel,
    title), ...) with key 's' | 'p' | 't'.  Each level is normalised by its OWN driving pressure
    (plateau load-piston increment); the reference (dashed, 95 % band) by the smallest of them.
    Dotted guide at 1."""
    L0 = _norm_ref_level(cfg, levels)
    if L0 is None:
        _norm_skip(cfg, None, what + ' sweep')
        return None
    pr, dp = _norm_symbols(L0)
    fig, axes = plt.subplots(1, 3, figsize=(25, 6.5), constrained_layout=True)
    dps = ', '.join(f'{sig(_dp_norm(cfg, L))}' for L in levels)
    fig.suptitle(f'{suptitle}  (${pr}={sig(_p_norm(cfg, L0))}$; ${dp}$ = {dps} per level, reference / {sig(_dp_norm(cfg, L0))})\n{cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    mask = _scale_mask(cfg, R, levels=levels)
    for ax, (comp, key, ylabel, title) in zip(axes, panels):
        if not any(_spt(L['stress'], comp) is not None for L in levels):
            ax.text(0.5, 0.5, 'sigmaxx / sigmayy files\nnot found', ha='center', va='center', transform=ax.transAxes)
            finish_axes(ax, ylabel, title, R=R)
            continue
        Nr = partial_norm(cfg, R, L0, comp, ref=True)
        ref = mean_ci(Nr[key], cfg.ci_level) if Nr is not None else None
        ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6, zorder=1)

        def stack(L, c=comp, k=key):
            N = partial_norm(cfg, R, L, c)
            return None if N is None else N[k]
        finals = overlay_levels(ax, R, levels, lambda L: L['z'], lambda L: L['ts'], stack, cfg, ref=ref,
                                autoscale_mask=mask, ylabel=ylabel, title=title)
        _norm_ylim(ax, finals, mask, pad=0.3)
        smart_legend(ax, handles=level_handles(levels, ref=ref is not None), fontsize=11)
    return _save(fig, cfg, stem)


def _sweep_partial_norm(cfg, R, levels, comp, stem):
    """Normalised companion of fig_partial_stress_sweep (comp = 'zz') and of
    fig_thermo_pressure_sweep (comp = 'tr'): (a) solvent*, (b) polymer*, (c) total*, all levels."""
    fam, sym = _norm_family(comp)
    pr, dp = _norm_symbols(None)
    name = r'Partial and total $P_{th}$' if comp == 'tr' else r'Partial and total $\sigma_{%s}$' % comp
    return _sweep_norm_panels(cfg, R, levels,
                              [(comp, k, rf'$({sym}-w\,{pr})/{dp}$', f'({"abc"[i]}) ' + lab) for i, (k, lab, _) in enumerate(fam)],
                              stem, name + r" with each species' share $w$ of the back pressure removed, all levels",
                              'thermodynamic-pressure' if comp == 'tr' else 'partial-stress')


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
    finish_axes(ax, ylabel, title, R=R)
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
            finish_axes(ax, r"$\sigma'_{%s}$" % comp, title, R=R)
            continue
        Rs = R['stress'].get(comp)
        _final_overlay(ax, cfg, R, levels, lambda L, c=comp: (L['stress'][c]['net'] if c in L['stress'] else None),
                       Rs['net'] if Rs is not None else None, r"$\sigma'_{%s}(z)$" % comp, title, matter_only=True)
    return _save(fig, cfg, 'sweep_network_stress_final')


def fig_thermo_pressure_sweep(cfg, R, levels):
    """P_th = -(1/3) tr(sigma^t) evolution, all levels overlaid (faint -> bold), reference dashed.
    With cfg.PARTIAL_NORM = 'share' the normalised solvent | polymer | total traces follow
    (sweep_thermo_pressure_evolution_norm, 2026-10-06)."""
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
    fig = _save(fig, cfg, 'sweep_thermo_pressure_evolution')
    if cfg.PARTIAL_NORM == 'share':
        _sweep_partial_norm(cfg, R, levels, 'tr', 'sweep_thermo_pressure_evolution_norm')
    return fig


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
    """(a) solvent partial, (b) polymer partial, (c) total sigma_zz evolutions, all levels.  With
    cfg.PARTIAL_NORM = 'share' the normalised version follows (sweep_partial_stress_evolution_norm,
    2026-10-06)."""
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
    fig = _save(fig, cfg, 'sweep_partial_stress_evolution')
    if cfg.PARTIAL_NORM == 'share':
        _sweep_partial_norm(cfg, R, levels, 'zz', 'sweep_partial_stress_evolution_norm')
    return fig


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
    eps = np.array([L['eps_M'] for L in levels])          # the strain M divides by (M_STRAIN)
    hp = [L for L in levels if 'M_pist' in L]
    ep = np.array([L['eps_M'] for L in hp])
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
    axA.set_xlabel(eps_axis_label(levels))
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
    axB.set_xlabel(eps_axis_label(levels))
    axB.set_ylabel(r'$M$  (LJ units)')
    axB.set_title('(b) longitudinal modulus per level, two independent estimates', fontsize=13)
    axB.grid(alpha=0.3)
    smart_legend(axB, fontsize=13)
    return _save(fig, cfg, 'sweep_modulus')


def fig_stress_strain_sweep(cfg, R, levels):
    """The stress-strain curve itself: plateau piston stress and interior network stress vs
    applied strain, INCLUDING each estimator's eps = 0 reading (hollow), with a dotted line
    through each series.  The line is ANCHORED at the estimator's own eps = 0 reading and its
    slope is the through-anchor least-squares fit of the first cfg.SS_FIT_NPTS levels only
    (2026-09-27; was a free-intercept fit of every point, whose intercept went negative as
    the stiffening levels pulled it).  The slope is the small-strain M of that estimator,
    referenced to its own zero like M_SUBTRACT_REF; the higher levels' departure from the
    line is the stress-strain nonlinearity.  (Was panel (b) of fig_M_sweep until 2026-09-14.)"""
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    eps = np.array([L['eps_M'] for L in levels])          # the strain M divides by (M_STRAIN), not the applied level
    hp = [L for L in levels if 'M_pist' in L]
    Rzz = R['stress']['zz']
    n_fit = int(max(1, cfg.SS_FIT_NPTS))
    fits = []

    def _series(e, s, err, color, marker, name, e0=None, s0=None, err0=None):
        ax.errorbar(e, s, yerr=err, fmt=marker + '-', lw=2, ms=9, color=color, capsize=6, label=name)
        e, s = np.asarray(e, float), np.asarray(s, float)
        anchor = float(s0) if (e0 is not None and s0 is not None and np.isfinite(s0)) else 0.0
        if e0 is not None and s0 is not None and np.isfinite(s0):
            ax.errorbar([e0], [s0], yerr=err0, fmt=marker, ms=9, mfc='none', color=color, capsize=6)
        order = np.argsort(e)
        ef, sf = e[order][:n_fit], s[order][:n_fit]
        if len(ef) >= 1:
            # least squares through the anchor (0, s0):  slope = sum eps (s - s0) / sum eps^2
            slope = float(np.sum(ef * (sf - anchor)) / np.sum(ef ** 2))
            xs = np.linspace(0, max(e) * 1.05, 20)
            ax.plot(xs, anchor + slope * xs, ls=':', lw=1.5, color=color, alpha=0.8)
            fits.append((name.split()[0], slope, len(ef)))

    # unrelaxed-hold systematic (RELAX_SYS): the plateau stress is over by the piston tail's excess,
    # so both series carry it on their LOWER bar (as M does through delta_sys)
    exc = np.array([L.get('relax', {}).get('excess', 0.0) for L in levels])
    sn = np.array([L['M_net_abs'] * L['eps_M'] for L in levels])
    sn_err = [np.abs(np.array([L['M_net_abs_lo'] * L['eps_M'] for L in levels]) - sn) + exc,
              np.abs(np.array([L['M_net_abs_hi'] * L['eps_M'] for L in levels]) - sn)]
    _series(eps, sn, sn_err, WONG['blue'], 'o', "network $\\langle\\sigma'_{zz}\\rangle_{\\rm int}$ (plateau)",
            e0=0.0, s0=float(Rzz['net_interior']), err0=float(Rzz.get('net_interior_half', 0.0)))
    if hp:
        ep = np.array([L['eps_M'] for L in hp])
        Pp = np.array([L['P_final'] for L in hp])
        pref = R.get('P_ref', np.nan)
        exc_p = np.array([L.get('relax', {}).get('excess', 0.0) for L in hp])
        _series(ep, Pp, [Pp - [L['PF']['lo'] for L in hp] + exc_p, [L['PF']['hi'] for L in hp] - Pp],
                WONG['vermillion'], 's', 'piston $P=\\langle F_z\\rangle/A$ (plateau)',
                e0=0.0, s0=pref, err0=(float(R['P_ref_hi'] - R['P_ref_lo']) / 2 if np.isfinite(pref) else None))
    ax.plot([], [], 'o', mfc='none', color='0.4', label=r'hollow = $\varepsilon=0$ reading')
    ax.plot([], [], ls=':', lw=1.5, color='0.4', label=f'dotted = through the $\\varepsilon=0$ reading, slope from the first {n_fit} levels')
    if np.any(exc > 0):
        ax.plot([], [], ' ', label='lower bars include the unrelaxed-hold excess')
    ax.set_title('Stress vs strain, small-strain slopes through $\\varepsilon=0$:  ' +
                 ',  '.join(f"$M_{{\\rm {n[:4]}}}\\approx{sig(s)}$ ({k} levels)" for n, s, k in fits) + '\n' + cfg.sim_name,
                 fontsize=12)
    ax.set_xlabel(eps_axis_label(levels))
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
        e = np.array([L['eps_M'] for L in Ls])
        g = np.array([L['G'][comp]['G'] for L in Ls])
        ax.errorbar(e, g, yerr=[g - [L['G'][comp]['lo'] for L in Ls], [L['G'][comp]['hi'] for L in Ls] - g],
                    fmt=mk + '-', ms=10, lw=2, color=col, capsize=6,
                    label=fr"$G$ from {comp}:  $(\sigma'_{{zz}}-\sigma'_{{{comp}}})/2\varepsilon$")
    ax.set_xlabel(eps_axis_label(hg))
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
    ax.plot(e, v, 'o-', color=WONG['reddishpurple'], lw=2, ms=8,
            label=rf"profile fit, first {cfg.DC_FIT_TAU1:g} $\tau_1$ of the hold" if cfg.DC_FIT_TAU1 > 0 else 'profile fit')
    for L in hd:   # R^2 label under each point, clipped inside the axes
        ax.annotate(f"$R^2$={sig(L['Dc']['R2'])}", (L['eps'], L['Dc']['Dc']), textcoords='offset points',
                    xytext=(0, -12), ha='center', va='top', fontsize=10, color='0.35',
                    annotation_clip=True)
    if cfg.DC_FIT_TAU1 > 0:            # 2026-10-09: the whole-hold fit and the stress-trace D_c alongside
        ax.plot(e, [L['Dc']['Dc_all'] for L in hd], 'o:', mfc='none', color=WONG['reddishpurple'], lw=1.5, ms=8,
                label='profile fit, whole hold')
    hs = [L for L in hd if (L.get('Dc_stress') or {}).get('ok')]
    if hs:
        ax.errorbar([L['eps'] for L in hs], [L['Dc_stress']['Dc'] for L in hs], yerr=[L['Dc_stress']['Dc_se'] for L in hs],
                    fmt='^--', color=WONG['green'], lw=1.8, ms=9, capsize=5, label='load-piston stress trace (held-slab series)')
    ax.margins(x=0.12, y=0.15)
    ax.axhline(cfg.DC_SLOW_REF / 4.0, color='0.4', ls=':', lw=1.5,
               label=f"deck hold-sizing $D_c$ = {sig(cfg.DC_SLOW_REF / 4)}  ($Dc_{{est}}$ = {sig(cfg.DC_SLOW_REF)} in the deck's $L^2/\\pi^2 D_c$ formula)")
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
        for i in F['shown']:
            c = cmap(norm(F['ts'][i]))
            ax.plot(F['zf'], u0_b + F['uhat'][i][F['idx']], 'o', color=c, ms=2.5, alpha=0.35)
            ax.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=1.5)
        dlL, fs = F['DL'] / F['L'], F['f_sup']
        ax.plot([0, 1], [fs * dlL, (fs - 1.0) * dlL], 'k:', lw=1.6)
        ax.set(xlabel=flip_hint(R, r'$\zeta$'), ylabel=r'$u_z/L$', xlim=_zlim(R))
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
        lo = k - np.array([L['kappa'][key]['lo'] for L in Ls])
        hi = np.array([L['kappa'][key]['hi'] for L in Ls]) - k
        ax.errorbar(e, k, yerr=[np.clip(lo, 0, None), np.clip(hi, 0, None)],      # an M CI through zero makes hi < 0
                    fmt=mk + '-', ms=10, lw=2, color=col, capsize=6, label=lab)
    hs = [L for L in hk if 'net' in (L.get('kappa_stress') or {})]
    if hs:
        ax.plot([L['eps'] for L in hs], [L['kappa_stress']['net'] for L in hs], '^--', ms=10, lw=1.8, color=WONG['green'],
                label=r'$D_c^\mathrm{stress\ trace}/M_\mathrm{network}$')
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


def _ref_thickness(R):
    """Reference gel thickness (sigma) for turning a strain offset into a length: the first of
    the thickness keys the reference loader sets, else the seated plate gap."""
    for k in ('L0', 'L0_bb', 'L_bb', 'L_gel', 'L_rg'):
        if np.isfinite(R.get(k, np.nan)):
            return float(R[k])
    return float(R['z_piston'] - R['z_support'])


def _disp_panel(ax, cfg, L, col, label=None, show_fit=True, ms=5):
    """One steady displacement profile with its line fit on `ax` (shared by the single and
    sweep figures).  Filled = interior bins in the fit, hollow = trimmed / tail bins."""
    D = L['disp_prof']
    zc, u, idx, inner = D['zc'], D['u'], D['idx'], D['inner']
    out = np.setdiff1d(idx, inner)
    ax.plot(zc[inner], u[inner], 'o', ms=ms, color=col, label=label)
    ax.plot(zc[out], u[out], 'o', ms=ms, color=col, mfc='none', alpha=0.7)
    if show_fit:
        zz = zc[idx]
        ax.plot(zz, D['slope'] * zz + D['icpt'], '-', lw=1.4, color=col, alpha=0.9)


# ---------------------------------------------------------------------------
#  Reservoir normal-stress check (2026-10-06)
# ---------------------------------------------------------------------------
_NCOL = {'xx': WONG['blue'], 'yy': WONG['orange'], 'zz': WONG['vermillion']}


def _perm_masks(cfg, R, L, ts):
    """Per-snapshot permeate-reservoir masks (between the permeate piston and the support,
    res_wall_margin clear of the piston), the mirror of L['bw'] for the lower reservoir;
    all-False rows when the run has no permeate piston."""
    z = R['z']
    if not np.isfinite(R.get('z_perm', np.nan)):
        return np.zeros((len(ts), len(z)), bool)
    pa = L['wetz']['perm_at'] if L is not None and L.get('wetz') else (lambda t: R['z_perm'])
    rows = []
    for t in ts:
        zs = _level_support_z(R, L, t)
        rows.append(reservoir_mask(z, cfg, 'perm', float(pa(t)), float(zs) if np.isfinite(zs) else R['z_gel_lo']))
    return np.array(rows, bool)


def _masked_series(stack, masks):
    """mean of stack[i] over masks[i] per snapshot (nan where the mask is empty)."""
    out = np.full(len(stack), np.nan)
    for i in range(len(stack)):
        v = stack[i][masks[i]]
        v = v[np.isfinite(v)]
        if len(v):
            out[i] = float(np.mean(v))
    return out


def reservoir_normal_stress(cfg, R, L):
    """Total normal stresses sigma^t_xx, yy, zz during a compression hold, in the gel interior
    and in the two solvent reservoirs, as INCREMENTS from the eps = 0 reference of the same
    region (plateau means, CI = t-interval over the plateau snapshots).  The reservoirs are
    pure fluid: their increments must be 0 in every component (the bath stays at P_bath);
    the gel's lateral increments are the sigma'_xx that sets G_comp = (M/2)(1 - sigma'_xx/
    sigma'_zz).  Also returns the per-snapshot series and the plateau profiles.  None when
    a component is missing."""
    if not all(c in L['stress'] and c in R['stress'] for c in COMPONENTS):
        return None
    ts = L['ts']
    pm = _perm_masks(cfg, R, L, ts)
    nref = len(R['stress']['zz']['t'])
    pm_ref = np.tile(_perm_masks(cfg, R, None, [0.0])[0], (nref, 1))      # reference pistons / support at rest
    regions = {'gel': ('interior', np.tile(L['interior'], (len(ts), 1)), np.tile(R['interior'], (nref, 1))),
               'feed': ('feed reservoir', L['bw'], np.tile(np.asarray(R['bw'], bool), (nref, 1))),
               'perm': ('permeate reservoir', pm, pm_ref)}
    out = dict(ts=ts, plat=L['plat'], series={}, ref={}, inc={}, prof={}, prof_ref={}, has_perm=bool(pm.any()))
    for reg, (lab, m, mref) in regions.items():
        out['series'][reg], out['ref'][reg], out['inc'][reg] = {}, {}, {}
        for c in COMPONENTS:
            ser = _masked_series(L['stress'][c]['t'], m)
            rser = _masked_series(R['stress'][c]['t'], mref)
            out['series'][reg][c] = ser
            rm = float(np.nanmean(rser)) if np.isfinite(rser).any() else np.nan
            out['ref'][reg][c] = rm
            v = ser[L['plat']] - rm
            v = v[np.isfinite(v)]
            if len(v) >= 2:
                h = stats.t.ppf(0.5 + cfg.ci_level / 2, len(v) - 1) * np.std(v, ddof=1) / np.sqrt(len(v))
                out['inc'][reg][c] = (float(np.mean(v)), float(h))
            elif len(v) == 1:
                out['inc'][reg][c] = (float(v[0]), np.nan)
            else:
                out['inc'][reg][c] = (np.nan, np.nan)
    for c in COMPONENTS:
        out['prof'][c] = np.nanmean(L['stress'][c]['t'][L['plat']], axis=0)
        out['prof_ref'][c] = np.nanmean(R['stress'][c]['t'], axis=0)
    g, z = out['inc']['gel'], out['inc']['gel']['zz'][0]
    lat = 0.5 * (g['xx'][0] + g['yy'][0])
    out['ratio'] = lat / z if z else np.nan                     # sigma'_xx / sigma'_zz (pore pressure cancels: same bath)
    out['G_comp_M'] = 0.5 * (1.0 - out['ratio'])                # G_comp / M implied by the lateral stress
    return out


def fig_reservoir_normal_stress(cfg, R, L):
    """(a) sigma^t_xx, yy, zz over the hold in the gel interior (solid) and in the solvent
    reservoirs (feed dashed, permeate dotted), with P_bath; (b) plateau-mean z-profiles of the
    three components (solid) over the reference (thin dashed).  The reservoirs are the in-situ
    control: a fluid must read P_bath in all three components through the whole hold, so a
    lateral drift there is a box / normalisation artefact, while a lateral rise confined to the
    gel is network stress -- the sigma'_xx of G_comp = (M/2)(1 - sigma'_xx/sigma'_zz)."""
    if _is_perm(L):
        print('reservoir normal-stress figure: compression holds only')
        return None
    Q = reservoir_normal_stress(cfg, R, L)
    if Q is None:
        print('reservoir normal-stress figure skipped (sigmaxx / sigmayy files missing)')
        return None
    L['res_check'] = Q
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(20, 6.8), constrained_layout=True)
    lab_eps = f'$\\varepsilon={L["eps"]:.2f}$' if np.isfinite(L['eps']) else ''
    fig.suptitle(f'Normal stresses in the gel vs the solvent reservoirs during the hold  |  {lab_eps}  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    ts = Q['ts']
    styles = {'gel': ('-', 2.4, 1.0, 'gel interior'), 'feed': ('--', 1.6, 0.9, 'feed reservoir'), 'perm': (':', 1.8, 0.9, 'permeate reservoir')}
    for reg, (ls, lw, al, name) in styles.items():
        if reg == 'perm' and not Q['has_perm']:
            continue
        for c in COMPONENTS:
            axA.plot(ts, Q['series'][reg][c], ls, color=_NCOL[c], lw=lw, alpha=al,
                     label=(f'$\\sigma^t_{{{c}}}$ {name}' if reg == 'gel' else None))
    for reg, (ls, lw, al, name) in styles.items():
        if reg == 'gel' or (reg == 'perm' and not Q['has_perm']):
            continue
        axA.plot([], [], ls, color='0.3', lw=lw, label=name + ' (all three components)')
    axA.axhline(cfg.P_BARO, color='k', ls=':', lw=1.2, alpha=0.7, label=f'$P_{{\\rm bath}}={sig(cfg.P_BARO)}$')
    axA.axvspan(L['halt_ts'], float(ts[-1]), color=WONG['green'], alpha=0.10, label='plateau window')
    axA.set_xlabel('time step'); axA.set_ylabel(r'$\sigma^t$  (region mean, compression-positive, LJ)')
    axA.set_title('(a) gel interior vs reservoirs over the hold'); axA.grid(alpha=0.3)
    # y-range from the plateau window (the early-hold consolidation transient would otherwise set the scale)
    pl = Q['plat']
    curves = [Q['series'][r][c][pl] for r in styles for c in COMPONENTS if not (r == 'perm' and not Q['has_perm'])]
    robust_ylim(axA, curves + [[cfg.P_BARO]], pad=0.8, qlo=0, qhi=100, include_zero=False)
    lines = ['plateau increments from the $\\varepsilon=0$ reference (region means):']
    for reg, (_, _, _, name) in styles.items():
        if reg == 'perm' and not Q['has_perm']:
            continue
        lines.append(name + ':  ' + '   '.join(f'$\\Delta\\sigma^t_{{{c}}}$ = {fmt_val_unc(*Q["inc"][reg][c])}' for c in COMPONENTS))
    lines.append(f"gel $\\sigma'_{{\\rm lat}}/\\sigma'_{{zz}}$ = {Q['ratio']:.3f}  $\\rightarrow$  $G_{{\\rm comp}}/M = \\frac{{1}}{{2}}(1 - \\sigma'_{{xx}}/\\sigma'_{{zz}})$ = {Q['G_comp_M']:.3f}")
    annotate_box(axA, '\n'.join(lines), loc='lower left', fontsize=10.5)
    smart_legend(axA, fontsize=10.5)
    zx = zn(R, R['z'])
    for c in COMPONENTS:
        axB.plot(zx, Q['prof_ref'][c], '--', color=_NCOL[c], lw=1.2, alpha=0.6, zorder=2)
        axB.plot(zx, Q['prof'][c], '-', color=_NCOL[c], lw=2.4, alpha=0.95, zorder=3, label=f'$\\sigma^t_{{{c}}}$ plateau')
    axB.plot([], [], '--', color='0.3', lw=1.2, label=r'reference ($\varepsilon=0$), same colours')
    axB.axhline(cfg.P_BARO, color='k', ls=':', lw=1.2, alpha=0.7, label=f'$P_{{\\rm bath}}={sig(cfg.P_BARO)}$')
    shade_gel(axB, R, L)
    mark_walls(axB, R, L)
    finish_axes(axB, r'$\sigma^t(z)$  (compression-positive, LJ)', '(b) plateau profiles: reservoirs must read $P_{\\rm bath}$ in all three', R=R)
    robust_ylim(axB, [Q['prof'][c] for c in COMPONENTS] + [Q['prof_ref'][c] for c in COMPONENTS] + [np.full(len(zx), cfg.P_BARO)],
                zmask=_scale_mask(cfg, R, L), pad=0.5, qlo=1, qhi=99, include_zero=False)
    smart_legend(axB, fontsize=11)
    return _save(fig, cfg, 'reservoir_normal_stress', L['lvl'])


def fig_reservoir_normal_stress_sweep(cfg, R, levels):
    """Plateau increments of sigma^t_xx, yy, zz from the reference, normalised by the gel's
    Delta sigma^t_zz, vs strain: gel interior (filled: = sigma'_lat/sigma'_zz) and the solvent
    reservoirs (hollow: feed squares, permeate triangles -- must sit at 0), plus the implied
    G_comp/M = (1 - sigma'_lat/sigma'_zz)/2 per level (right panel) against G_REF / M_net of
    the lowest level when cfg.G_REF is set."""
    Qs = [(L, L.get('res_check') or reservoir_normal_stress(cfg, R, L)) for L in levels if not _is_perm(L)]
    Qs = [(L, Q) for L, Q in Qs if Q is not None]
    if not Qs:
        print('reservoir normal-stress sweep figure skipped (no level with the three total-stress components)')
        return None
    for L, Q in Qs:
        L['res_check'] = Q
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(19, 6.5), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.35, 1]))
    fig.suptitle(f'Normal-stress increments vs strain: gel interior vs solvent reservoirs (plateau means)  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    eps = np.array([L['eps'] for L, Q in Qs])
    off = {'xx': -0.004, 'yy': 0.0, 'zz': 0.004}
    zz_gel = np.array([Q['inc']['gel']['zz'][0] for L, Q in Qs])
    for c in COMPONENTS:
        if c != 'zz':
            g = np.array([Q['inc']['gel'][c] for L, Q in Qs])
            axA.errorbar(eps + off[c], g[:, 0] / zz_gel, yerr=g[:, 1] / zz_gel, fmt='o-', color=_NCOL[c], ms=7, lw=1.6, capsize=3,
                         label=f'gel interior $\\Delta\\sigma^t_{{{c}}} / \\Delta\\sigma^t_{{zz}}$  ($= \\sigma\'_{{{c}}}/\\sigma\'_{{zz}}$)')
        f = np.array([Q['inc']['feed'][c] for L, Q in Qs])
        axA.errorbar(eps + off[c], f[:, 0] / zz_gel, yerr=f[:, 1] / zz_gel, fmt='s', mfc='none', color=_NCOL[c], ms=7, lw=1.2, capsize=3,
                     label=f'feed reservoir $\\Delta\\sigma^t_{{{c}}} / \\Delta\\sigma^t_{{zz,\\rm gel}}$')
        if any(Q['has_perm'] for L, Q in Qs):
            pr = np.array([Q['inc']['perm'][c] for L, Q in Qs])
            axA.errorbar(eps + off[c], pr[:, 0] / zz_gel, yerr=pr[:, 1] / zz_gel, fmt='^', mfc='none', color=_NCOL[c], ms=7, lw=1.2, capsize=3,
                         label=f'permeate reservoir $\\Delta\\sigma^t_{{{c}}} / \\Delta\\sigma^t_{{zz,\\rm gel}}$')
    axA.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    axA.axhline(1, color='k', ls=':', lw=1, alpha=0.5)
    axA.set_xlabel(eps_axis_label([L for L, Q in Qs]))
    axA.set_ylabel(r'$\Delta\sigma^t / \Delta\sigma^t_{zz,\rm gel}$  (plateau increments)')
    axA.set_title("(a) lateral / axial increment: filled = gel interior, hollow = reservoirs (must be 0)"); axA.grid(alpha=0.3)
    axA.set_ylim(-0.1, 1.1)
    smart_legend(axA, fontsize=10)
    ratio = np.array([Q['G_comp_M'] for L, Q in Qs])
    axB.plot(eps, ratio, 'o-', color=WONG['blue'], ms=8, lw=2, label=r"$\frac{1}{2}\left(1-\sigma'_{\rm lat}/\sigma'_{zz}\right)$ from the hold")
    Ms = np.array([L.get('M_net', np.nan) for L, Q in Qs], float)
    if cfg.G_REF is not None and np.isfinite(Ms).any():
        i0 = int(np.nanargmin(np.where(np.isfinite(Ms), eps, np.inf)))      # small-strain M: the lowest level
        M0 = float(Ms[i0])
        axB.axhline(cfg.G_REF / M0, color=WONG['vermillion'], ls='--', lw=1.8,
                    label=f'$G_{{\\rm shear}}/M$ = {cfg.G_REF / M0:.3f}  (G_REF {cfg.G_REF:g}, $M_{{\\rm net}}$ {M0:.3f} at $\\varepsilon$ = {eps[i0]:.2f})')
    axB.axhline(0.5, color='k', ls=':', lw=1, alpha=0.6, label=r'$\nu = 0$ ($\sigma_{\rm lat}^\prime = 0$)')
    axB.axhline(0.0, color='k', ls='--', lw=1, alpha=0.4, label=r'$\nu = 0.5$ ($\sigma_{\rm lat}^\prime = \sigma_{zz}^\prime$)')
    axB.set_xlabel(eps_axis_label([L for L, Q in Qs])); axB.set_ylabel(r'$G/M$')
    axB.set_title(r'(b) $G_{\rm comp}/M$ implied by the lateral network stress'); axB.grid(alpha=0.3)
    axB.set_ylim(-0.05, 0.6)
    smart_legend(axB, fontsize=11)
    return _save(fig, cfg, 'sweep_reservoir_normal_stress')


def fig_disp_profile(cfg, R, L):
    """Steady displacement profile of one level and the network strain it gives (2026-10-03):
    (a) u_z of the polymer, plateau frames vs the eps = 0 reference, binned by the reference
        position Z (or by the current z for the disp_z_polymer_cum fallback), with the line fit
        over the interior bins -> eps_disp = -slope; (b) the residuals.  The slope is the strain
        of the network alone: the applied strain also counts the rigid travel (the gel falling
        onto the support), which the title quantifies as the offset."""
    D = L.get('disp_prof')
    if D is None:
        print('displacement-profile figure skipped (no local traj_ref / traj_stress and no disp_z_polymer_cum)')
        return None
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(11, 9), sharex=True, constrained_layout=True,
                                   gridspec_kw={'height_ratios': [3, 1.3]})
    flip_z_axis(axA, R)                                    # shared x: both panels
    coord = 'reference position $Z$' if D['coord'] == 'lagrangian' else 'current position $z$'
    fig.suptitle(f"Steady displacement profile of the network   |   {cfg.sim_name}   |   level $\\varepsilon = {L['lvl']}$\n"
                 f"{D['src']}: {D['n_frames']} plateau frame(s) ({D['ts'][0]:.0f}-{D['ts'][-1]:.0f}) vs the $\\varepsilon=0$ reference, "
                 f"{cfg.binWidth:g}$\\sigma$ bins by the {coord}", fontsize=12, fontweight='bold')
    _disp_panel(axA, cfg, L, WONG['blue'], ms=6,
                label=f"polymer $u_z$ per bin (filled = {len(D['inner'])} interior bins in the fit, hollow = trimmed {cfg.DISP_TRIM_BINS} per face)")
    zz = D['zc'][D['idx']]
    axA.plot(zz, D['slope'] * zz + D['icpt'], '-', lw=1.8, color=WONG['vermillion'],
             label=(f"line fit:  $u_z = u_0 - \\varepsilon_{{\\rm disp}}\\,(Z - Z^P)$,  "
                    f"$\\varepsilon_{{\\rm disp}} = {D['eps']:.4f}$ [{D['eps_lo']:.4f}, {D['eps_hi']:.4f}],  $R^2 = {D['R2']:.4f}$"))
    for zp, lab in ((D['z_P'], 'support plane (reference)'), (D['z_T'], 'load-piston plane (reference)')):
        axA.axvline(zp, color='0.4', ls='--', lw=1.0, label=lab)
    axA.axhline(0, color='k', lw=0.8, alpha=0.5)
    axA.set_ylabel(r'$u_z$  ($\sigma$)')
    ci = int(cfg.ci_level * 100)
    txt = (f"$\\varepsilon_{{\\rm disp}}$ = {D['eps']:.4f}  ({ci}% CI half-width {max(D['eps_hi'] - D['eps'], 0):.4f}: "
           f"LSQ {D['eps_half_lsq']:.4f}, frame-to-frame {D['eps_half_frames']:.4f})\n"
           f"$\\varepsilon_{{Rg}}$ = {L.get('eps_rg', np.nan):.4f}      applied $\\varepsilon$ = {L['eps']:.4f}  "
           f"(rigid travel counted by the applied strain: {L['eps'] - D['eps']:+.4f} = {(L['eps'] - D['eps']) * _ref_thickness(R):+.2f}$\\sigma$)\n"
           f"with $\\varepsilon_{{\\rm {L['eps_M_src']}}}$ (M_STRAIN):  $M_{{\\rm pist}}$ = {L.get('M_pist', np.nan):.3f},  "
           f"$M_{{\\rm net}}$ = {L['M_net']:.3f};  with the applied strain they would be "
           f"{L.get('dP_pist', np.nan) / L['eps']:.3f}, {L['dsig_net'] / L['eps']:.3f}")
    annotate_box(axA, txt, loc='upper right', fontsize=10.5)
    axA.set_title('(a) $u_z$ of the polymer in the plateau, relative to the $\\varepsilon=0$ reference', fontsize=13)
    axA.grid(alpha=0.3)
    smart_legend(axA, fontsize=10)
    axB.plot(D['zc'][D['inner']], D['resid'][D['inner']], 'o-', ms=5, lw=1, color=WONG['blue'])
    out = np.setdiff1d(D['idx'], D['inner'])
    axB.plot(D['zc'][out], D['resid'][out], 'o', ms=5, color=WONG['blue'], mfc='none', alpha=0.7)
    axB.axhline(0, color='k', lw=0.8)
    axB.set_ylabel(r'$u_z - $ fit  ($\sigma$)')
    axB.set_xlabel(flip_hint(R, 'reference position $Z$  ($\\sigma$)' if D['coord'] == 'lagrangian' else 'current position $z$  ($\\sigma$)'))
    axB.set_title('(b) residuals from the line (uniform-stress interior should be flat; the faces and contact layers are hollow)', fontsize=12)
    axB.grid(alpha=0.3)
    return _save(fig, cfg, 'disp_profile', L['lvl'])


def fig_disp_profile_sweep(cfg, R, levels):
    """Steady displacement profiles across the sweep (2026-10-03):
    (a) every level's profile with its line fit; (b) the slope strain eps_disp and eps_Rg vs
    the applied strain (the constant offset is the rigid travel the applied strain counts);
    (c) M_pist and M_net per level under the three strain definitions -- the stress increments
    are the same, only the denominator changes -- with the small-strain slope of each."""
    hd = [L for L in levels if L.get('disp_prof') is not None]
    if not hd:
        print('displacement-profile sweep figure skipped (no level has a profile)')
        return None
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(21, 6.5), constrained_layout=True)
    flip_z_axis(axA, R)
    fig.suptitle('Steady displacement profiles -> the network strain that divides $M$   |   ' + cfg.sim_name,
                 fontsize=13, fontweight='bold')
    for i, L in enumerate(hd):
        D = L['disp_prof']
        _disp_panel(axA, cfg, L, level_color(levels.index(L)),
                    label=f"$\\varepsilon$ = {L['lvl']}:  $\\varepsilon_{{\\rm disp}}$ = {D['eps']:.4f}  ($R^2$ = {D['R2']:.4f})")
    D0 = hd[0]['disp_prof']
    for zp in (D0['z_P'], D0['z_T']):
        axA.axvline(zp, color='0.4', ls='--', lw=1.0)
    axA.axhline(0, color='k', lw=0.8, alpha=0.5)
    axA.set_xlabel(flip_hint(R, 'reference position $Z$  ($\\sigma$)' if D0['coord'] == 'lagrangian' else 'current position $z$  ($\\sigma$)'))
    axA.set_ylabel(r'$u_z$  ($\sigma$)')
    axA.set_title('(a) plateau $u_z$ vs the $\\varepsilon=0$ reference, line fits over the interior (filled)', fontsize=12)
    axA.grid(alpha=0.3)
    smart_legend(axA, fontsize=10)
    # ---- (b) strains ----
    ea = np.array([L['eps'] for L in hd])
    ed = np.array([L['disp_prof']['eps'] for L in hd])
    ed_err = [ed - [L['disp_prof']['eps_lo'] for L in hd], [L['disp_prof']['eps_hi'] for L in hd] - ed]
    er = np.array([L.get('eps_rg', np.nan) for L in hd])
    axB.plot([0, ea.max() * 1.05], [0, ea.max() * 1.05], ':', color='0.4', lw=1.2, label='1:1')
    axB.errorbar(ea, ed, yerr=ed_err, fmt='o-', ms=9, lw=2, color=WONG['blue'], capsize=5,
                 label=r'$\varepsilon_{\rm disp}$ (profile slope)')
    if np.isfinite(er).any():
        axB.plot(ea, er, 's--', ms=7, lw=1.5, color=WONG['vermillion'], mfc='none', label=r'$\varepsilon_{Rg}$ (strain_zz plateau)')
    off = ea - ed
    axB.set_xlabel(r'applied strain  $\varepsilon$')
    axB.set_ylabel('measured network strain')
    axB.set_title(f"(b) applied - $\\varepsilon_{{\\rm disp}}$ = {off.mean():+.4f} ± {off.std(ddof=1) if len(off) > 1 else 0:.4f}\n"
                  f"(= {off.mean() * _ref_thickness(R):+.2f}$\\sigma$ of rigid travel, the same at every level)", fontsize=12)
    axB.grid(alpha=0.3)
    smart_legend(axB, fontsize=11)
    # ---- (c) M under each strain definition ----
    n_fit = int(max(1, cfg.SS_FIT_NPTS))
    hp = [L for L in hd if 'dP_pist' in L]
    rows = []
    for key, Ls, col, mk, name in (('dsig_net', hd, WONG['blue'], 'o', 'network'), ('dP_pist', hp, WONG['vermillion'], 's', 'piston')):
        if not Ls:
            continue
        ds = np.array([L[key] for L in Ls])
        for src, style, lab in (('applied', dict(mfc='none', ls=':', alpha=0.8, lw=1.3, ms=8), 'applied $\\varepsilon$'),
                                ('rg', dict(mfc='none', ls='--', alpha=0.8, lw=1.3, ms=8, marker='x'), '$\\varepsilon_{Rg}$'),
                                ('disp', dict(ls='-', lw=2.2, ms=10), '$\\varepsilon_{\\rm disp}$')):
            e = np.array([{'applied': L['eps'], 'rg': L.get('eps_rg', np.nan), 'disp': L['disp_prof']['eps']}[src] for L in Ls])
            if not np.isfinite(e).any():
                continue
            style = dict(style)
            marker = style.pop('marker', mk)
            x = np.array([L['eps'] for L in Ls])
            M = ds / e
            o = np.argsort(x)[:n_fit]
            slope = float(np.sum(e[o] * ds[o]) / np.sum(e[o] ** 2))
            axC.plot(x, M, marker=marker, color=col, label=f"{name} / {lab}:  small-strain $M$ = {slope:.3f}", **style)
            rows.append((name, src, slope))
    axC.set_xlabel(r'applied strain  $\varepsilon$  (level)')
    axC.set_ylabel(r'$M = \Delta\sigma / \varepsilon_{\rm x}$  (LJ)')
    axC.set_title(f"(c) the same stress increments / three strains\n(slope = through-origin fit of the first {n_fit} levels)", fontsize=12)
    axC.grid(alpha=0.3)
    smart_legend(axC, fontsize=10)
    return _save(fig, cfg, 'sweep_disp_profile')


# ---------------------------------------------------------------------------
#  Compression record -> permeation consistency check (2026-10-03)
# ---------------------------------------------------------------------------
# The permeation drive and the compression sweep measure the same modulus two ways, and on
# 2026-10-02/03 they disagreed by 10-20 % after the strain definitions were made consistent.
# So the sweep notebook SAVES its stress-strain curve (save_M_record) and the permeation
# notebook CHECKS against it every time it runs (fig_perm_vs_compression): the steady
# permeation frames are paired per atom with the zero-flux reference frames (Lagrangian, so
# the rigid drop onto the support cancels and the contact layer is included), and the
# measured network deformation is compared with the one the compression curve predicts for
# the same drag load, sigma'(Z) = dP (Z_F - Z)/L_0.  A verdict line is printed and put in
# the figure title; the check also refuses to stay silent when the record is missing.
def save_M_record(cfg, R, levels, path=None):
    """Write the sweep's stress-strain record (per level: the strain M divides by, the network
    and piston stress increments with CIs, the M values) to <DATA_DIR>/M_record.json and to
    <base>/compression/M_record_latest.json (the one the permeation check reads by default)."""
    import json, datetime
    rows = []
    for L in sorted(levels, key=lambda L: L['eps']):
        rows.append(dict(level=str(L['lvl']), eps_applied=float(L['eps']), eps=float(L['eps_M']), eps_src=L['eps_M_src'],
                         eps_half=float(L.get('eps_M_half', 0.0)), eps_rg=float(L.get('eps_rg', np.nan)),
                         dsig_net=float(L['dsig_net']), dsig_net_half=float((L['M_net_hi'] - L['M_net_lo']) / 2 * L['eps_M']),
                         dP_pist=float(L.get('dP_pist', np.nan)),
                         dP_pist_half=float((L['M_pist_hi'] - L['M_pist_lo']) / 2 * L['eps_M']) if 'M_pist' in L else np.nan,
                         M_net=float(L['M_net']), M_pist=float(L.get('M_pist', np.nan)), M_sys=float(L.get('M_sys', 0.0))))
    n_fit = int(max(1, cfg.SS_FIT_NPTS))
    e = np.array([r['eps'] for r in rows[:n_fit]])
    slopes = {k: float(np.sum(e * np.array([r[k] for r in rows[:n_fit]])) / np.sum(e ** 2)) for k in ('dsig_net', 'dP_pist')}
    rec = dict(run_id=cfg.RUN_ID, sim_name=cfg.sim_name, saved=datetime.datetime.now().isoformat(timespec='minutes'),
               M_STRAIN=cfg.M_STRAIN, M_SUBTRACT_REF=bool(cfg.M_SUBTRACT_REF), L0=_ref_thickness(R),
               z_support=float(R['z_support']), z_piston=float(R['z_piston']),
               M_small_net=slopes['dsig_net'], M_small_pist=slopes['dP_pist'], n_fit=n_fit, levels=rows)
    outs = [Path(path)] if path else [cfg.DATA_DIR / 'M_record.json', cfg.DATA_DIR.parent / 'M_record_latest.json']
    for o in outs:
        o.parent.mkdir(parents=True, exist_ok=True)
        o.write_text(json.dumps(rec, indent=1, default=lambda x: None if (isinstance(x, float) and not np.isfinite(x)) else float(x)))
    print(f"compression record saved ({len(rows)} levels, strain = eps_{cfg.M_STRAIN}, small-strain M net {slopes['dsig_net']:.3f} / "
          f"pist {slopes['dP_pist']:.3f}):\n  " + '\n  '.join(str(o) for o in outs))
    return rec


def load_M_record(cfg):
    """The compression record for the permeation check (cfg.M_RECORD or the newest), or None."""
    import json
    p = Path(cfg.M_RECORD) if cfg.M_RECORD else Path(cfg.base_dir) / 'compression' / 'M_record_latest.json'
    if not p.exists():
        return None
    rec = json.loads(p.read_text())
    rec['path'] = str(p)
    return rec


def _record_curve(rec, key):
    """eps(sigma') of the compression record (through the origin, linear beyond the last level)."""
    rows = [r for r in rec['levels'] if r.get(key) is not None and np.isfinite(r[key])]
    sg = np.array([0.0] + [r[key] for r in rows]); ep = np.array([0.0] + [r['eps'] for r in rows])
    o = np.argsort(sg); sg, ep = sg[o], ep[o]

    def f(x):
        x = np.asarray(x, float)
        y = np.interp(x, sg, ep)
        hi = x > sg[-1]
        if hi.any() and len(sg) > 2:
            y[hi] = ep[-1] + (x[hi] - sg[-1]) * (ep[-1] - ep[-2]) / (sg[-1] - sg[-2])
        return y
    return f


def perm_lagrangian_profile(cfg, R, P):
    """Steady permeation displacement in the MATERIAL frame: the steady frames of traj_stress
    paired per atom with the zero-flux traj_ref frames, binned by the reference position Z.
    -> dict(zc, u, per_frame, n_atoms, idx, eps_loc (= -du/dZ, central differences), ts) or None."""
    tp = cfg.traj('traj_stress')
    if not tp.exists():
        return None
    ref = load_ref_positions(cfg, R)
    if ref is None:
        return None
    ts_all = traj_timesteps(tp)
    plat = [t for t in ts_all if t >= P['halt_ts']] or ts_all[-1:]
    if cfg.DISP_MAX_FRAMES > 0 and len(plat) > cfg.DISP_MAX_FRAMES:
        plat = list(subsample(np.array(plat), np.array(plat)[:, None], cfg.DISP_MAX_FRAMES)[0])
    fr = read_traj_id_z(tp, plat)
    U, Z = [], None
    for t in sorted(fr):
        ids, z = fr[t]
        if np.array_equal(ids, ref['ids']):
            Zt, u = ref['z'], z - ref['z']
        else:
            _, ia, ib = np.intersect1d(ids, ref['ids'], return_indices=True)
            Zt, u = ref['z'][ib], z[ia] - ref['z'][ib]
        if Z is None or len(Zt) == len(Z):
            Z, U = Zt, U + [u]
    if not U:
        return None
    U = np.array(U)
    edges = np.arange(np.floor(Z.min()), Z.max() + cfg.binWidth, cfg.binWidth)
    k = np.clip(np.digitize(Z, edges) - 1, 0, len(edges) - 2)
    nb = len(edges) - 1
    cnt = np.bincount(k, minlength=nb).astype(float)
    per = np.array([np.bincount(k, weights=u, minlength=nb) / np.maximum(cnt, 1) for u in U])
    zc = 0.5 * (edges[:-1] + edges[1:])
    um = per.mean(axis=0)
    idx = np.where(cnt >= cfg.Ncount_min)[0]
    eps_loc = np.full(nb, np.nan)
    if len(idx) >= 3:
        zi, ui = zc[idx], um[idx]
        g = np.gradient(ui, zi)
        eps_loc[idx] = -g
    return dict(zc=zc, u=um, per_frame=per, n_atoms=cnt, idx=idx, eps_loc=eps_loc, ts=[int(t) for t in sorted(fr)],
                n_frames=len(U), ref_ts=ref['ts'])


def perm_vs_compression(cfg, R, P, verbose=True):
    """The consistency check (see the section comment).  Returns a dict with the measured and
    predicted deformation, the ratio, the verdict and the pieces the figure draws; None when
    the permeation M is not available.  Prints the verdict."""
    say = print if verbose else (lambda *a, **k: None)
    Mp = P.get('M_perm'); F = P.get('Dc')
    if not Mp or F is None:
        say('CHECK vs compression: skipped (no permeation M)')
        return None
    rec = load_M_record(cfg)
    out = dict(rec=rec, dP=Mp['dP'], L0=F['L0'], z_P=F['z_P'], Z_P=F['Z_P'], Z_F=F['Z_F'], M_perm=Mp)
    lag = P.get('lag') if P.get('lag') is not None else perm_lagrangian_profile(cfg, R, P)
    out['lag'] = lag
    # ---- measured deformation of the network, rigid drop excluded ----
    #   (i) the fig-12 parabola: c L_0 (contact-layer u_0 subtracted, Eulerian bins)
    #   (ii) the Lagrangian profile: u(top material bin) - u(bottom material bin), which cancels
    #        any rigid translation and INCLUDES the contact layer's own compaction
    out['uF_prof'] = Mp['prof']['uF'] if 'prof' in Mp else np.nan
    if lag is not None and len(lag['idx']) >= 4:
        i0, i1 = lag['idx'][0], lag['idx'][-1]
        out['uF_lag'] = float(lag['u'][i1] - lag['u'][i0])
        out['u_contact'] = float(lag['u'][i0])                 # the rigid drop + whatever the lowest layer did
        out['Z_bot'], out['Z_top'] = float(lag['zc'][i0]), float(lag['zc'][i1])
        pf = lag['per_frame']
        if len(pf) >= 2:
            d = pf[:, i1] - pf[:, i0]
            out['uF_lag_half'] = float(stats.t.ppf(0.5 + cfg.ci_level / 2, len(d) - 1) * d.std(ddof=1) / np.sqrt(len(d)))
        else:
            out['uF_lag_half'] = np.nan
    else:
        out['uF_lag'] = np.nan
    # ---- support force (deck output since 2026-10-03): the network load at the contact ----
    fs = cfg.path('support_force')
    out['F_supp'] = None
    if fs.exists():
        try:
            t = read_print_file(fs, ['step', 'F_poly', 'F_solv'])
            st = np.where(t['step'] >= P['halt_ts'])[0]
            st = st if len(st) else np.arange(len(t['step']))
            A = R['AREA']
            out['F_supp'] = dict(P_poly=float(np.mean(-t['F_poly'][st]) / A), P_solv=float(np.mean(-t['F_solv'][st]) / A),
                                 P_poly_half=float(stats.t.ppf(0.5 + cfg.ci_level / 2, max(len(st) - 1, 1))
                                                   * np.std(t['F_poly'][st], ddof=1) / np.sqrt(len(st)) / A) if len(st) > 1 else np.nan)
        except Exception as e:                                  # noqa: BLE001
            say(f'  support_force file present but unreadable: {e}')
    if rec is None:
        out['verdict'] = 'NO COMPRESSION RECORD'
        say('CHECK vs compression: NO RECORD -- run the compression sweep notebook (it saves flow_data_local/compression/'
            'M_record_latest.json via tri.save_M_record) and rerun this cell.  The permeation M is unchecked.')
        return out
    # ---- prediction from the compression curve for the same drag load ----
    Z = np.linspace(F['Z_P'], F['Z_F'], 400)
    sig = Mp['dP'] * (F['Z_F'] - Z) / F['L0']                   # drag load, linear in the material coordinate
    preds = {}
    for key, name in (('dsig_net', 'network'), ('dP_pist', 'piston')):
        f = _record_curve(rec, key)
        eps = f(sig)
        u = -np.concatenate([[0.0], np.cumsum(0.5 * (eps[1:] + eps[:-1]) * np.diff(Z))])
        preds[name] = dict(Z=Z, eps=eps, u=u, uF=float(u[-1]), eps_mean=float(-u[-1] / F['L0']),
                           M_eff=float(Mp['dP'] * F['L0'] / (2 * abs(u[-1]))), M_small=rec['M_small_net' if key == 'dsig_net' else 'M_small_pist'])
    out['pred'] = preds
    uF_pred = float(np.mean([p['uF'] for p in preds.values()]))
    out['uF_pred'] = uF_pred
    # the same material span as the Lagrangian measurement
    if np.isfinite(out['uF_lag']):
        up = {n: float(np.interp(out['Z_top'], Z, p['u']) - np.interp(out['Z_bot'], Z, p['u'])) for n, p in preds.items()}
        out['uF_pred_span'] = float(np.mean(list(up.values())))
        out['ratio_lag'] = out['uF_lag'] / out['uF_pred_span']
    else:
        out['uF_pred_span'] = out['ratio_lag'] = np.nan
    out['ratio_prof'] = out['uF_prof'] / uF_pred if np.isfinite(out['uF_prof']) else np.nan
    r = out['ratio_lag'] if np.isfinite(out['ratio_lag']) else out['ratio_prof']
    out['ratio'] = r
    out['ratio_src'] = 'Lagrangian profile' if np.isfinite(out['ratio_lag']) else 'fig-12 parabola'
    tol = float(cfg.PERM_M_TOL)
    out['ok'] = bool(np.isfinite(r) and abs(r - 1) <= tol)
    out['verdict'] = ('CONSISTENT' if out['ok'] else 'DISCREPANCY') + f': deformation {r:.3f} x the compression prediction ({out["ratio_src"]}; tolerance {tol:.0%})'
    say(f"CHECK vs compression ({rec['run_id']}, saved {rec['saved']}, strain eps_{rec['M_STRAIN']}):  {out['verdict']}\n"
        f"  measured |u_F| (rigid drop excluded): parabola {abs(out['uF_prof']):.2f} sigma, Lagrangian {abs(out['uF_lag']):.2f} sigma"
        + (f" (span {out['Z_bot']:.0f}-{out['Z_top']:.0f}, contact layer u = {out['u_contact']:+.2f} sigma)" if np.isfinite(out['uF_lag']) else '')
        + f";  predicted from the compression curve for dP = {Mp['dP']:.4f}: {abs(uF_pred):.2f} sigma (span-matched {abs(out['uF_pred_span']):.2f})\n"
        f"  M: permeation {Mp[Mp['primary']]['M']:.3f} ({Mp['primary']}, L = {Mp['L_src']}; L_0 basis {Mp[Mp['primary']].get('M_L0', np.nan):.3f})"
        f"  vs compression small-strain {rec['M_small_net']:.3f} (net) / {rec['M_small_pist']:.3f} (pist), "
        f"effective at this load {preds['network']['M_eff']:.3f} / {preds['piston']['M_eff']:.3f}"
        + (f"\n  network load at the support (polymer-support pair force / A): {out['F_supp']['P_poly']:.4f} vs dP_ext {Mp['dP']:.4f}"
           f"  (solvent-support {out['F_supp']['P_solv']:.4f}, should be 0: transparent plate)" if out['F_supp'] else
           '\n  (no support_force file: deck output since 2026-10-03 -- the network load at the support is not checked)'))
    return out


def fig_perm_vs_compression(cfg, R, P):
    """The permeation-vs-compression check as a figure (verdict in the title):
    (a) the steady displacement in the material frame (per-atom pairing with the zero-flux
        reference) against the deformation the compression curve predicts for the drag load,
        shifted by the measured contact-layer displacement; the fig-12 parabola for reference;
    (b) the local strain -du/dZ against the predicted eps(sigma'(Z));
    (c) M: the permeation estimates vs the compression small-strain slopes and the effective
        modulus the compression curve gives at this load."""
    C = perm_vs_compression(cfg, R, P, verbose=False)
    if C is None:
        print('check figure skipped (no permeation M)')
        return None
    perm_vs_compression(cfg, R, P, verbose=True)            # the printed verdict, always
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(21, 6.5), constrained_layout=True)
    flip_z_axis(axA, R)
    flip_z_axis(axB, R)
    ok = C.get('ok', False)
    col_v = WONG['green'] if ok else WONG['vermillion']
    fig.suptitle(f"CHECK -- permeation vs compression sweep:  {C['verdict']}\n{cfg.sim_name}", fontsize=13,
                 fontweight='bold', color=col_v)
    Mp, F, lag = C['M_perm'], P['Dc'], C['lag']
    # ---- (a) ----
    if lag is not None:
        i = lag['idx']
        axA.plot(lag['zc'][i], lag['u'][i], 'o', ms=5, color=WONG['blue'],
                 label=f"measured: steady frames ({lag['n_frames']}) vs zero-flux reference, per atom, by reference $Z$")
    zeta = np.linspace(0, 1, 200)
    if 'prof' in Mp:
        axA.plot(F['Z_P'] + zeta * F['L0'], C.get('u_contact', F['u0_ss']) + Mp['prof']['uF'] * zeta * (2 - zeta), '--', color='0.3', lw=1.3,
                 label=f"fig-12 parabola (+ contact-layer $u$): $u_F$ = {Mp['prof']['uF']:.2f}$\\sigma$")
    if C.get('pred'):
        for name, p, col in (('network', C['pred']['network'], WONG['vermillion']), ('piston', C['pred']['piston'], WONG['orange'])):
            axA.plot(p['Z'], C.get('u_contact', F['u0_ss']) + p['u'], '-', lw=1.8, color=col, alpha=0.9,
                     label=f"predicted from the compression {name} curve: $u_F$ = {p['uF']:.2f}$\\sigma$")
    for zp, lab in ((F['Z_P'], 'contact plane $Z^P$'), (F['Z_F'], 'free face $Z^F$')):
        axA.axvline(zp, color='0.4', ls=':', lw=1)
    axA.axhline(0, color='k', lw=0.8, alpha=0.5)
    axA.set_xlabel(flip_hint(R, 'reference (material) position $Z$  ($\\sigma$)'))
    axA.set_ylabel('$u_z$ since the zero-flux reference  ($\\sigma$)')
    axA.set_title('(a) steady displacement in the material frame\nvs the compression-curve prediction', fontsize=12)
    axA.grid(alpha=0.3)
    smart_legend(axA, fontsize=10)
    # ---- (b) ----
    if lag is not None:
        i = lag['idx']
        axB.plot(lag['zc'][i], lag['eps_loc'][i], 'o-', ms=4, lw=1, color=WONG['blue'], label=r'measured  $-\mathrm{d}u/\mathrm{d}Z$')
    if C.get('pred'):
        for name, p, col in (('network', C['pred']['network'], WONG['vermillion']), ('piston', C['pred']['piston'], WONG['orange'])):
            axB.plot(p['Z'], p['eps'], '-', lw=1.8, color=col, label=f"predicted $\\varepsilon(\\sigma'(Z))$, {name} curve")
    if 'prof' in Mp:
        axB.plot(F['Z_P'] + zeta * F['L0'], -2 * Mp['prof']['uF'] / F['L0'] * (1 - zeta), '--', color='0.3', lw=1.3, label='fig-12 parabola')
    axB.axhline(0, color='k', lw=0.8, alpha=0.5)
    axB.set_xlabel(flip_hint(R, 'reference (material) position $Z$  ($\\sigma$)'))
    axB.set_ylabel('local compressive strain')
    axB.set_title("(b) local strain $-\\mathrm{d}u/\\mathrm{d}Z$ vs the compression curve at the\nlocal drag load $\\sigma'(Z) = \\Delta P\\,(Z^F - Z)/L_0$ (material-linear)", fontsize=12)
    axB.grid(alpha=0.3)
    smart_legend(axB, fontsize=10)
    # ---- (c) ----
    bars = []
    for k in ('lag', 'prof', 'bb', 'trace'):
        if k in Mp:
            bars.append((f"permeation\n{k}" + (' (primary)' if k == Mp['primary'] else ''), Mp[k]['M'], Mp[k]['lo'], Mp[k]['hi'],
                         WONG['skyblue'] if k == 'lag' else WONG['blue']))
    if C.get('rec'):
        rec = C['rec']
        bars.append(('compression\nsmall-strain (net)', rec['M_small_net'], np.nan, np.nan, WONG['vermillion']))
        bars.append(('compression\nsmall-strain (pist)', rec['M_small_pist'], np.nan, np.nan, WONG['orange']))
        bars.append(('compression curve\neffective at this load', C['pred']['network']['M_eff'], np.nan, np.nan, WONG['green']))
    for j, (lab, m, lo, hi, col) in enumerate(bars):
        axC.bar(j, m, color=col, alpha=0.85)
        if np.isfinite(lo):
            axC.errorbar(j, m, yerr=[[m - lo], [hi - m]], color='k', capsize=5, lw=1.2)
        axC.text(j, m, f'{m:.3f}', ha='center', va='bottom', fontsize=10)
    axC.set_xticks(range(len(bars)))
    axC.set_xticklabels([b[0] for b in bars], fontsize=8)
    axC.set_ylabel('$M$  (LJ)')
    axC.set_title(f"(c) $M$: permeation (L = {Mp['L_src']}) vs the compression sweep\n" + (f"[{C['rec']['run_id']}]" if C.get('rec') else '[NO RECORD]'), fontsize=12)
    axC.grid(axis='y', alpha=0.3)
    txt = (f"support load: polymer {C['F_supp']['P_poly']:.4f} vs $\\Delta P_{{\\rm ext}}$ {Mp['dP']:.4f}" if C.get('F_supp')
           else 'no support_force file\n(deck output since 2026-10-03)')
    annotate_box(axC, txt, loc='upper left', fontsize=10)
    return _save(fig, cfg, 'perm_vs_compression')


def print_summary(cfg, levels):
    """One line per level with the headline numbers, then the hold check."""
    ci = int(cfg.ci_level * 100)
    how = 'M = increment from the eps = 0 reference' if cfg.M_SUBTRACT_REF else 'M = absolute stress / eps'
    if cfg.M_SUBTRACT_REF and any(L.get('M_pist_ref') == 'absent' for L in levels):
        how += '; piston M ABSOLUTE (no piston_force_avg_ref file)'
    print(f'\nSUMMARY  ({cfg.sim_name}; {ci}% CIs; {how})')
    if cfg.RELAX_SYS:
        how += '; M lower bounds include delta_sys (unrelaxed-hold excess, last column)'
    print(f"{'eps':>6s} {'eps_M':>9s} {'M_net':>18s} {'M_pist':>18s} {'G_x':>18s} {'G_y':>18s} {'D_c':>10s} {'D_c_all':>10s} {'D_c_stress':>10s} {'D_c_unload':>10s} {'kappa_net':>10s} {'kappa_pist':>10s} {'delta_sys':>10s}")
    for L in levels:
        def ci_(v, lo, hi):
            return f'{v:.3f} [{lo:.3f},{hi:.3f}]'
        em = f"{L['eps_M']:.4f}{ {'disp': 'd', 'rg': 'r', 'applied': 'a'}[L['eps_M_src']] }"
        ds = L.get('M_sys', 0.0)
        ds_s = f"{ds:.3f}" + (' *' if ds > 0.5 * (L['M_net_hi'] - L['M_net']) else '') if cfg.RELAX_SYS else ''
        mp = ci_(L['M_pist'], L['M_pist_lo'], L['M_pist_hi']) if 'M_pist' in L else 'n/a'
        gx = ci_(L['G']['xx']['G'], L['G']['xx']['lo'], L['G']['xx']['hi']) if 'xx' in L['G'] else 'n/a'
        gy = ci_(L['G']['yy']['G'], L['G']['yy']['lo'], L['G']['yy']['hi']) if 'yy' in L['G'] else 'n/a'
        dc = f"{L['Dc']['Dc']:.3e}" + (' !' if L['Dc'].get('stab_flag') else '') if L.get('Dc') else 'n/a'
        da = f"{L['Dc']['Dc_all']:.3e}" if L.get('Dc') else 'n/a'
        S = L.get('Dc_stress')
        dst = f"{S['Dc']:.3e}" if (S and S.get('ok')) else 'n/a'
        Uf = L.get('unload')
        dun = f"{Uf['trace']['Dc']:.3e}" if (Uf and Uf.get('trace')) else (f"{Uf['prof']['Dc']:.3e}p" if (Uf and Uf.get('prof')) else 'n/a')
        kn = f"{L['kappa']['net']['k']:.3e}" if L.get('kappa', {}).get('net') else 'n/a'
        kp = f"{L['kappa']['pist']['k']:.3e}" if L.get('kappa', {}).get('pist') else 'n/a'
        print(f"{L['eps']:6.3f} {em:>9s} {ci_(L['M_net'], L['M_net_lo'], L['M_net_hi']):>18s} {mp:>18s} {gx:>18s} {gy:>18s} "
              f"{dc:>10s} {da:>10s} {dst:>10s} {dun:>10s} {kn:>10s} {kp:>10s} {ds_s:>10s}")
    print('  eps_M = the strain M, G and kappa divide by (M_STRAIN): d = steady displacement-profile slope, r = Rg, a = applied')
    print(f'  D_c = profile fit over the first {cfg.DC_FIT_TAU1:g} tau_1 of the hold (kappa uses it; ! = the whole-hold fit D_c_all differs by '
          f'more than {cfg.DC_STAB_FLAG:.0%}: a slower process after the transient);  D_c_stress = the load-piston trace, held-slab series;  '
          f'D_c_unload = the free re-swelling after the hold, thickness trace (p: profile fit), tau_1 = L^2/(pi^2 D_c)')
    if cfg.RELAX_SYS:
        print('  delta_sys = (piston plateau - fitted tail asymptote) / eps;  * = larger than the M_net bootstrap half-width')
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
    # the feed face comes down with the piston: each snapshot's own gel top (its polymer-stress edge)
    top = np.array([float(z[np.abs(q) > cfg.gel_thresh * float(np.nanmax(np.abs(q)))].max()) for q in zz['p']])
    bw = baseline_masks_at(z, cfg, R, ts, P['wetz']['feed_at'], top)
    P['bw'] = bw
    off = [int(t) for i, t in enumerate(ts) if z[bw[i]].max() + 0.5 * cfg.binWidth > P['wetz']['feed_at'](t)]
    if off:
        say(f'  WARNING: no usable feed-reservoir bin at step(s) {off} -- the pore baseline window there lies beyond the feed piston')
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
        finish_axes(ax, r'$P_{\rm local}$  (LJ)', r'(c) calibration pressure handed to $\lambda(\phi_p, P)$', R=R)
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
        finish_axes(ax, r'$\lambda$', r'(d) $\lambda(\phi_p^{\rm vor}, P_{\rm local})$ per bin  (reservoir bins: exactly 1)', R=R)
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


# ---------------------------------------------------------------------------
#  permeation mode: D_c and M from the polymer displacement (2026-09-28)
#
#  Linear 1-D poroelasticity of the membrane z^P <= z <= z^F between two reservoirs whose
#  pressures the pistons PRESCRIBE (P_perm at the support, P_feed = P_perm + dP at the free feed
#  face).  Mixture incompressibility + Darcy + quasi-static force balance give
#      du_z/dt = q(t) + D_c d2u_z/dz2,   q(t) = phi_s v_s + phi_p v_p  (the TOTAL flux: uniform in z
#  but a function of time -- the permeate flux that overshoots 25x in perm_2; = Q_perm/A only in
#  the steady state).  The total stress is uniform, so the network stress at the support is -dP
#  from t = 0+ while sigma' = 0 at the free face: the STRAIN eps = -du_z/dz obeys the diffusion
#  equation with DIRICHLET ends, eps(z^P) = dP/M, eps(z^F) = 0, and relaxes with
#      tau_1 = L^2 / (pi^2 D_c)     (4x the held slab's L^2/(4 pi^2 D_c) of the compression fit for the
#                                    same D_c: see "D_c: one equation, two boundary-value problems" above fit_Dc).
#  With zeta = (z - z^P)/L and mu_k = k pi the displacement modes are phi_k = 1 - cos(mu_k zeta)
#  (u(0) = 0 and u'(1) = 0 built in), decaying as exp(-mu_k^2 D_c t/L^2); the steady state is
#      u_ss(zeta) = -(dP L / 2M) [zeta^2 - 2 zeta],   u_F = u_ss(1) = -dP L/(2M),
#      M = dP L / (2 |u_F|) = dP / (2 |c|),   u_ss(zeta)/L = c zeta (2 - zeta).
#  NOTE (2026-09-28, second pass): the first version took q constant and expanded u in
#  sin((k - 1/2) pi zeta), the modes of a prescribed FLUX (tau_1 = 4 L^2/pi^2 D_c); with the
#  pressure drop prescribed the decay is 4x faster for the same D_c, so that reading overstated
#  D_c by 4 (1.37 instead of ~0.34 for perm_2).
#  The drive is the applied-dP HISTORY (piston_pressure: increments w_j at t_j, sum w_j = 1;
#  a 5 k-step ramp is one increment, the 1M-step ramp of decks since 2026-09-28 ~200).  Each
#  mode responds with  H_k(t) = sum_j w_j [1 - exp(-mu_k^2 D_c (t - t_j)/L^2)]_(t > t_j),  so
#      u_z(zeta, t) = sum_k a_k phi_k(zeta) H_k(t)   (a_k fixed by the zero IC: only odd k, and
#      u_F(t) = u_F sum_k b_k H_k(t),  b_k = 8/(k pi)^2 for odd k, sum b_k = 1: the feed-face trace).
#  The deck resets displace/atom at the ramp start (since 2026-09-28) or at the ramp end
#  (older decks, PERM_DISP_RESET): the data are u(t) - u(t_reset), so the profile fit uses
#  H_k(t) - H_k(t_reset) and the amplitudes stay free.
#
#  Rigid drop (2026-09-28 pm, perm_2): at zero flux the polymer face sits ~3 sigma above the
#  contact plane z^P (a solvent layer between the plate and the network, see the reference
#  polymer profile); the drag closes it and the whole network arrives at the plate displaced by
#  u_0 ~ -3 sigma with no strain attached.  The contact-layer bins then read u_0 for the rest of
#  the run.  Left in the data, u_0 makes the pinned modes ring (a step at zeta = 0 built from 5
#  cosines) and biases the steady parabola; PERM_RIGID = 'contact' subtracts it (u_0(t) from the
#  PERM_PIN_NBINS lowest populated bins) and PERM_COORDS = 'lagrangian' places each bin at its
#  material coordinate zeta = (z - u_z - Z^P)/L_0, Z^P = z^P - u_0 (the emptying top bins then
#  reach zeta -> 1 instead of stopping at the last bin populated in EVERY snapshot).
#  Three readings of the displacement data:
#    * profile fit  (fig 12 a-b): (u_z - u_0)(zeta, t)/L of the disp_z_polymer snapshots with FREE
#      amplitudes (the compression fit's idiom, DC_FREE_AMPS) -> D_c; u = 0 at t_reset whatever
#      the amplitudes.
#    * thickness trace (fig 12 c): u_F(t) = L_bb(t) - L_0 (5 k-step cadence) against the
#      zero-IC series -> (u_F(inf), D_c); this resolves the early transient that the
#      stress-cadence snapshots (first window >= 1 nfreq after the reset) missed in perm_2.
#    * M (fig 13): M = dP L_ss/(2|u_F|) (L_ss the steady thickness; 2026-10-02, was L_0) from the
#      steady thickness change L_0 - L_ss (the definition), from the
#      steady-profile parabola c, and from the trace asymptote u_F(inf).
# ---------------------------------------------------------------------------
_PERM_TRACE_MODES = 60        # terms of the exact zero-IC series (coefficients fall as 1/mu^3)


def _perm_mu(n):
    """wavenumbers of the strain modes sin(mu_k zeta), mu_k = k pi (Dirichlet strain at both faces);
    lam_k = mu_k^2 D_c/L^2.  The displacement modes are phi_k = 1 - cos(mu_k zeta)."""
    return np.arange(1, n + 1) * np.pi


def _perm_phi(zh, mu):
    """displacement modes phi_k(zeta) = 1 - cos(mu_k zeta): (len(zh), len(mu))"""
    return 1.0 - np.cos(np.outer(np.atleast_1d(np.asarray(zh, float)), np.asarray(mu, float)))


def _perm_H(mu, Dc, L, t, forcing, W=0.0, nsub=12):
    """Modal response to the applied-dP history: H_k(t) = sum_j w_j [1 - exp(-lam_k (t - t_j))]
    over the increments already applied (t > t_j), lam_k = mu_k^2 Dc/L^2, sum_j w_j = 1 (so
    H_k -> 1; one increment at t = 0 is the step response 1 - exp(-lam_k t)).  W > 0 averages
    over the ave/chunk window [t - W, t] (nsub sub-samples; the deck block-averages u_z over
    the whole nfreq interval and tags the snapshot with the window END).  t: LJ time from
    the ramp start, scalar or array; returns (len(t), n_modes)."""
    mu = np.atleast_1d(np.asarray(mu, float))
    lam = (mu / L) ** 2 * Dc
    t = np.atleast_1d(np.asarray(t, float))
    tj, wj = forcing
    tt = (t[:, None] - W * (np.arange(nsub) + 0.5)[None, :] / nsub).ravel() if W > 0 else t
    dt = tt[:, None] - tj[None, :]
    on = dt > 0
    dtp = np.where(on, dt, 0.0)
    out = np.empty((len(tt), len(lam)))
    for k, lk in enumerate(lam):
        out[:, k] = np.sum(wj[None, :] * on * (1.0 - np.exp(-lk * dtp)), axis=1)
    if W > 0:
        out = out.reshape(len(t), nsub, len(lam)).mean(axis=1)
    return out


def _perm_trace_series(t_lj, Dc, L, forcing, n=_PERM_TRACE_MODES):
    """u_F(t)/u_F(inf) for the zero initial condition under the applied-dP history:
    sum_k b_k H_k(t) over the odd modes, b_k = 8/(k pi)^2 (sum b_k = 1)."""
    mu = (2.0 * np.arange(1, n + 1) - 1.0) * np.pi
    b = 8.0 / mu ** 2
    return _perm_H(mu, Dc, L, t_lj, forcing) @ b


def perm_forcing(cfg):
    """The applied-dP history from piston_pressure (P_feed_app - P_perm_app, one sample per
    volume_freq) -> dict(t0 = ramp start step (the last zero sample), t_end = step at which the
    full dP is first reached, full, n = number of increments, tj = increment times (LJ time from
    t0; each increment sits mid-way through its sampling interval), wj = weights (sum 1), steps,
    dp).  None when the file or the columns are missing."""
    names, tab = read_piston_table(cfg.path('piston_pressure'))
    if tab.size == 0 or names is None:
        return None
    fa, pa = table_col(names, tab, 'P_feed_app'), table_col(names, tab, 'P_perm_app')
    if fa is None or pa is None:
        return None
    st, dp = tab[:, 0], fa - pa
    full = float(np.nanmax(dp))
    if not full > 0:
        return None
    i0 = int(np.where(dp > 0.005 * full)[0][0])
    t0 = float(st[i0 - 1]) if i0 > 0 else float(st[i0])
    i_end = int(np.argmax(dp >= 0.999 * full))
    dsamp = float(np.median(np.diff(st))) if len(st) > 1 else 0.0
    a = max(i0 - 1, 0)
    inc = np.diff(dp[a:i_end + 1])
    tinc = st[a + 1:i_end + 1] - 0.5 * dsamp
    keep = inc > 0
    inc, tinc = inc[keep], tinc[keep]
    if inc.size == 0:
        inc, tinc = np.array([full]), np.array([t0])
    fz = dict(t0=t0, t_end=float(st[i_end]), full=full, n=int(inc.size), steps=st, dp=dp,
              tj=(tinc - t0) * cfg.dt_lj, wj=inc / inc.sum(), forcing='applied')
    if cfg.PERM_FORCING == 'measured':
        pf, pp = load2c(cfg.path('pressure_feed'), 2), load2c(cfg.path('pressure_permeate'), 2)
        if pf is not None and pp is not None and len(pf) >= 8:
            stm = pf[:, 0]
            dpm = pf[:, 1] - np.interp(stm, pp[:, 0], pp[:, 1])
            if cfg.PERM_FORCING_SMOOTH > 1:
                dpm = rolling_mean(dpm, int(cfg.PERM_FORCING_SMOOTH))
            after = stm > t0
            n_ss = max(int(0.25 * after.sum()), 4)
            dp_ss = float(np.mean(dpm[after][-n_ss:]))            # the steady measured dP normalises the weights
            if dp_ss > 0 and after.sum() >= 4:
                # the deck logs the reservoir pressures from the production start, i.e. after the ramp, so the
                # reservoir dP during the ramp is unobserved: it is taken to rise linearly from 0 at the ramp start
                # to the first measured sample (the quarter runs: the reservoir lags the applied ramp, 0.003 of
                # 0.1 at the end of a 200k-step ramp, 0.036 of 0.02 after 400k)
                t_first = float(stm[after][0])
                dsm = float(np.median(np.diff(stm)))
                n_pre = max(int(round((t_first - t0) / dsm)), 1)
                st_pre = t0 + (np.arange(1, n_pre + 1)) * (t_first - t0) / n_pre
                dp_pre = dpm[after][0] * np.arange(1, n_pre + 1) / n_pre
                st2 = np.concatenate([st_pre, stm[after][1:]])
                dp2 = np.concatenate([dp_pre, dpm[after][1:]])
                tinc_m = st2 - 0.5 * np.diff(np.concatenate([[t0], st2]))
                inc_m = np.diff(np.concatenate([[0.0], dp2]))
                tj_app = st_pre
                fz.update(tj=(tinc_m - t0) * cfg.dt_lj, wj=inc_m / dp_ss, full=dp_ss, n=int(inc_m.size), forcing='measured',
                          dp_meas=dpm, steps_meas=stm, dp_applied_full=full, n_ramp_prefix=int(n_pre),
                          dp_meas_max=float(np.max(dpm[after])), step_meas_max=float(stm[after][np.argmax(dpm[after])]))
    return fz


def load_perm_disp(cfg, R, P):
    """disp_z_polymer (fix ave/chunk of the per-atom z displacement since the deck's reset,
    binned by CURRENT z) + the applied-dP history + the reset convention + the polymer
    bounding-box thickness trace -> dict, or None when the displacement file is missing /
    has < 2 snapshots."""
    f = cfg.path('disp_z_polymer')
    if not f.exists():
        return None
    snaps = read_ave_chunk_file(f)
    if len(snaps) < 2:
        return None
    d = dict(ts=np.array([s[0] for s in snaps], float), z=snaps[0][1][:, 1],
             Nc=np.array([s[1][:, 2] for s in snaps]), uz=np.array([s[1][:, 3] for s in snaps]))
    d['stride'] = float(np.min(np.diff(d['ts'])))       # = the ave/chunk window (nfreq = nevery x nrepeat in the deck)
    fz = perm_forcing(cfg)
    if fz is None:
        t0 = float(d['ts'][0] - d['stride'])              # the first window must lie entirely after the reset
        fz = dict(t0=t0, t_end=t0, full=np.nan, n=1, tj=np.array([0.0]), wj=np.array([1.0]),
                  source='assumed: a step at the first displacement window start (no usable piston_pressure)')
    elif fz.get('forcing') == 'measured':
        fz['source'] = (f"MEASURED reservoir dP (pressure_feed - pressure_permeate, {cfg.PERM_FORCING_SMOOTH}-sample smoothing; "
                        f"linear rise over {fz['n_ramp_prefix']} increments to the first sample): "
                        f"steady {fz['full']:.4f}, peak {fz['dp_meas_max']:.4f} at step {int(fz['step_meas_max'])} "
                        f"(applied {fz['dp_applied_full']:.4f})")
    else:
        fz['source'] = 'piston_pressure (applied P_feed_app - P_perm_app)'
    d['forcing'] = fz
    mode = cfg.PERM_DISP_RESET
    if mode == 'auto':
        # the profile fix is created at the reset, so its FIRST WINDOW starts within one stride of
        # it: a window starting before the end of the ramp means a ramp-start reset (the snapshot
        # itself can land after a ramp shorter than two strides -- the local smoke test)
        mode = 'ramp_start' if d['ts'][0] - d['stride'] < fz['t_end'] else 'ramp_end'
    d['reset_mode'] = mode
    d['t_reset'] = fz['t0'] if mode == 'ramp_start' else fz['t_end']
    bb = load2c(cfg.path('gel_dimensions_bb'), 4)
    if bb is not None:
        d['bb_step'], d['bb_L'] = bb[:, 0], bb[:, 3]
    return d


def fit_perm_Dc(cfg, R, P, disp, free_offset=False):
    """One-sided consolidation fit of a permeation run (see the section comment):
    D_c from the displacement profiles (free amplitudes) and from the feed-face
    thickness trace (exact zero-IC series), M from the steady thickness change,
    the steady-profile parabola and the trace asymptote.  Returns a dict or None."""
    if disp is None or not np.isfinite(R.get('z_support', np.nan)):
        return None
    ts, z, Nc, uz = disp['ts'], disp['z'], disp['Nc'], disp['uz']
    fz = disp['forcing']
    forcing = (fz['tj'], fz['wj'])
    t_on = fz['t0']                                        # the ramp start = t = 0 of the model
    t_reset_lj = (disp['t_reset'] - t_on) * cfg.dt_lj
    z_P = float(R['z_support']) + cfg.PERM_GAP
    # ---- L_0 (bounding box): thickness over the PERM_L0_STEPS before the ramp start -- the reference of the
    #      thickness trace; the fit's own L_0 is chosen below (PERM_L0_MODE), once u_0 is known ----
    L0_bb_ci = (np.nan, np.nan)
    if 'bb_L' in disp:
        w0 = (disp['bb_step'] >= t_on - cfg.PERM_L0_STEPS) & (disp['bb_step'] < t_on)
        if w0.sum() >= 4:
            L0_bb, lo, hi, *_ = block_bootstrap_ci(disp['bb_L'][w0], cfg.ci_level)
            L0_bb_ci = (lo, hi)
        else:
            L0_bb = float(disp['bb_L'][disp['bb_step'] <= t_on][-1]) if (disp['bb_step'] <= t_on).any() else float(disp['bb_L'][0])
        L0_bb_source = 'gel_dimensions_bb'
    else:
        L0_bb = float(R['z_gel_hi'] + 0.5 * cfg.binWidth - z_P)
        L0_bb_source = 'reference polymer-stress edge'
    if not L0_bb > 0:
        return None
    n_snap = len(ts)
    steady = np.where(ts >= P['halt_ts'])[0]
    if len(steady) == 0:
        steady = np.array([n_snap - 1])
    # ---- populated bins per snapshot: the top bins empty as the face comes down, so the old
    #      all-snapshot intersection threw the early snapshots' top ~10 % away ----
    pop = [np.where(Nc[i] > cfg.Ncount_min)[0] for i in range(n_snap)]
    nb = max(int(cfg.PERM_PIN_NBINS), 1)
    if min(len(p) for p in pop) < 4 + 2 * cfg.DC_TRIM_BINS + nb:
        return None
    # ---- rigid shift u_0(t) = displacement of the contact layer (Nc-weighted mean of the nb lowest
    #      populated bins).  At zero flux a solvent layer separates the plate and the polymer face
    #      (perm_2: the face sat ~3 sigma above z^P = z_support + PERM_GAP); the drag closes it and
    #      the whole network arrives at the plate displaced by u_0 with no strain attached.  The
    #      pinned modes 1 - cos(mu zeta) vanish at zeta = 0 and cannot carry it: left in, u_0 shows
    #      up as ringing of the fitted modes (a step at zeta = 0 built from 5 cosines) and as a
    #      parabola c biased by the small-zeta bins (u/[zeta(2 - zeta)] -> large). ----
    if cfg.PERM_RIGID == 'contact':
        u0 = np.array([np.average(uz[i][p[:nb]], weights=Nc[i][p[:nb]]) for i, p in enumerate(pop)])
        u0_sd = np.array([np.sqrt(np.average((uz[i][p[:nb]] - u0[i]) ** 2, weights=Nc[i][p[:nb]]))
                          for i, p in enumerate(pop)])
    else:
        u0 = np.zeros(n_snap)
        u0_sd = np.zeros(n_snap)
    u0_ss = float(np.median(u0[steady]))
    Z_P = z_P - u0_ss                                      # reference (material) position of the pinned face
    # ---- the fit's L_0: dense-gel thickness Z_top - Z^P ('face') or the bounding box ('bb') ----
    if cfg.PERM_L0_MODE == 'face' and np.isfinite(R.get('z_gel_hi', np.nan)):
        L0 = float(R['z_gel_hi'] + 0.5 * cfg.binWidth - Z_P)
        L0_ci = (np.nan, np.nan)
        L0_source = f"reference polymer-stress edge {R['z_gel_hi'] + 0.5 * cfg.binWidth:.1f} - Z^P (PERM_L0_MODE 'face'; bounding box {L0_bb:.2f})"
    else:
        L0, L0_ci, L0_source = L0_bb, L0_bb_ci, L0_bb_source
    if not L0 > 0:
        return None
    z_F = z_P + L0
    Z_F = Z_P + L0
    # ---- coordinates and the deformation, per snapshot ----
    lagr = cfg.PERM_COORDS == 'lagrangian'
    zl, yl, il = [], [], []
    for i, p in enumerate(pop):
        if cfg.DC_TRIM_BINS:
            p = p[cfg.DC_TRIM_BINS:len(p) - cfg.DC_TRIM_BINS]
        zi = (z[p] - uz[i][p] - Z_P) / L0 if lagr else (z[p] - z_P) / L0
        yi = (uz[i][p] - u0[i]) / L0
        ok = (zi > 0) & (zi < 1)
        zl.append(zi[ok])
        yl.append(yi[ok])
        il.append(p[ok])
    if min(len(zz) for zz in zl) < 4:
        return None
    zeta = (z - z_P) / L0                                  # Eulerian bin coordinate (bookkeeping)
    uhat = uz / L0
    t_lj = (ts - t_on) * cfg.dt_lj
    W = disp['stride'] * cfg.dt_lj if cfg.PERM_DC_WINDOW_AVG else 0.0
    early = np.where(t_lj <= cfg.DC_FRAC_EARLY * t_lj[-1])[0]
    if len(early) < 2:
        return None
    n_modes = int(cfg.PERM_DC_N_MODES) or int(cfg.DC_N_MODES)
    mu = _perm_mu(n_modes)
    modes_l = [_perm_phi(zl[i], mu) for i in range(n_snap)]          # (points_i, modes)
    if free_offset:                                                   # strain-PDE fit (fig_perm_v0_check): the per-snapshot
        modes_l = [m - m.mean(axis=0, keepdims=True) for m in modes_l]    # rigid translation projected out of data and modes
        yl_fit = [y - y.mean() for y in yl]
    else:
        yl_fit = yl
    y_all = np.concatenate([yl_fit[i] for i in early])

    def dH(Dc):
        """H_k(t_i) - H_k(t_reset) for the fitted snapshots: (n_early, modes)"""
        H = _perm_H(mu, Dc, L0, t_lj[early], forcing, W)
        return H - _perm_H(mu, Dc, L0, t_reset_lj, forcing)[0][None, :]

    def amps(Dc, D=None):
        D = dH(Dc) if D is None else D
        X = np.vstack([modes_l[i] * D[n][None, :] for n, i in enumerate(early)])
        return np.linalg.lstsq(X, y_all, rcond=None)[0]

    def predict(Dc, A=None):
        D = dH(Dc)
        A = amps(Dc, D) if A is None else A
        return np.concatenate([(modes_l[i] * D[n][None, :]) @ A for n, i in enumerate(early)]), A

    def resid(Dc):
        p, _ = predict(Dc)
        return float(np.sum((p - y_all) ** 2))

    # coarse log-spaced scan first: the residual has secondary minima once the transient is
    # under-sampled (perm_2: 6 late snapshots), and bounded Brent alone landed on one at D_c ~ 30
    grid = np.geomspace(cfg.PERM_DC_BOUNDS[0], cfg.PERM_DC_BOUNDS[1], 61)
    rg = np.array([resid(g) for g in grid])
    j = int(np.nanargmin(rg))
    lo_b, hi_b = grid[max(j - 1, 0)], grid[min(j + 1, len(grid) - 1)]
    Dc = float(minimize_scalar(resid, bounds=(lo_b, hi_b), method='bounded').x)
    p_all, A = predict(Dc)
    ss_t = np.sum((y_all - np.mean(y_all)) ** 2)
    R2 = float(1.0 - np.sum((y_all - p_all) ** 2) / ss_t) if ss_t > 1e-30 else np.nan
    H_reset = _perm_H(mu, Dc, L0, t_reset_lj, forcing)[0]

    def u_model(zh, t, Wm=W):
        """deformation (u_z - u_0)/L_0 since the RESET at zeta = zh (array) and time t (LJ from the
        ramp start, scalar); window-averaged like the data."""
        zh = np.atleast_1d(np.asarray(zh, float))
        return _perm_phi(zh, mu) @ (A * (_perm_H(mu, Dc, L0, t, forcing, Wm)[0] - H_reset))

    def u_phys(zh, t):
        """deformation since the ramp start (instantaneous; = u_model + the fitted offset at the reset)."""
        zh = np.atleast_1d(np.asarray(zh, float))
        return _perm_phi(zh, mu) @ (A * _perm_H(mu, Dc, L0, t, forcing)[0])

    def u_inf(zh):
        """the reset-referenced steady state the snapshots approach"""
        zh = np.atleast_1d(np.asarray(zh, float))
        return _perm_phi(zh, mu) @ (A * (1.0 - H_reset))

    tau1 = L0 ** 2 / (np.pi ** 2 * Dc)                         # first-mode time (mu_1 = pi: Dirichlet strain at both faces)
    T_prod = float(t_lj[-1])
    F = dict(Dc=Dc, A=A, R2=R2, L0=L0, L0_ci=L0_ci, L0_source=L0_source, L0_bb=L0_bb, L0_bb_ci=L0_bb_ci, free_offset=bool(free_offset),
             z_P=z_P, z_F=z_F, Z_P=Z_P, Z_F=Z_F,
             u0=u0, u0_sd=u0_sd, u0_ss=u0_ss, rigid=cfg.PERM_RIGID, pin_nbins=nb, lagr=lagr, n_modes=n_modes,
             t_onset=t_on, forcing=fz, t_reset=disp['t_reset'], t_reset_lj=t_reset_lj, reset_mode=disp['reset_mode'],
             H_reset=H_reset, stride=disp['stride'], W=W, zeta=zeta, uhat=uhat, pop=pop, zl=zl, yl=yl, il=il,
             zf_lo=float(min(zz[0] for zz in zl)), zf_hi=float(max(zz[-1] for zz in zl)),
             n_pts=int(sum(len(zz) for zz in zl)), steady=steady,
             early=early, t_lj=t_lj, ts=ts, mu=mu, u_model=u_model, u_phys=u_phys, u_inf=u_inf, tau1=tau1,
             T_prod=T_prod, n_tau=T_prod / tau1, unrelaxed=float(8.0 / np.pi ** 2 * np.exp(-T_prod / tau1)))

    # ---- steady profile: parabola c zeta (2 - zeta) by least squares over the steady snapshots'
    #      points (the old per-bin mean of u/[zeta(2 - zeta)] let the small-zeta bins dominate) ----
    def _c(zz, yy):
        sh = zz * (2.0 - zz)
        return float(np.sum(yy * sh) / np.sum(sh * sh))

    c_snap = np.array([_c(zl[i], yl[i]) for i in steady])
    zs_ = np.concatenate([zl[i] for i in steady])
    ys_ = np.concatenate([yl[i] for i in steady])
    c = _c(zs_, ys_)
    sh_ = zs_ * (2.0 - zs_)
    res_c = ys_ - c * sh_
    se_lsq = float(np.sqrt(np.sum(res_c ** 2) / max(len(ys_) - 1, 1) / np.sum(sh_ * sh_)))
    h = _Z95 * se_lsq
    if len(steady) >= 2:
        h = max(h, _Z95 * float(np.std(c_snap, ddof=1)) / np.sqrt(len(steady)))
    ss_p = np.sum((ys_ - np.mean(ys_)) ** 2)
    F.update(c=c, c_lo=c - h, c_hi=c + h, c_npts=int(len(ys_)), c_snap=c_snap, u_ss=ys_, zeta_ss=zs_,
             parab_R2=float(1.0 - np.sum(res_c ** 2) / ss_p) if ss_p > 1e-30 else np.nan)

    # ---- feed-face thickness trace: exact zero-IC series, (u_F(inf), D_c) by least squares ----
    F['trace'] = None
    if 'bb_L' in disp:
        wt = disp['bb_step'] >= t_on
        st, uF = disp['bb_step'][wt], disp['bb_L'][wt] - L0_bb      # the trace is a thickness CHANGE: bounding box minus bounding box
        tt = (st - t_on) * cfg.dt_lj
        fit = (tt <= cfg.DC_FRAC_EARLY * tt[-1]) & (st >= t_on + cfg.PERM_TRACE_SKIP)
        if fit.sum() >= 8:
            def model(t, u_inf_, Dc_):
                return u_inf_ * _perm_trace_series(t, Dc_, L0, forcing)
            try:
                p0 = (float(np.mean(uF[-max(3, len(uF) // 10):])), float(cfg.DC_SLOW_REF))
                p0 = (p0[0], float(np.clip(p0[1], *cfg.PERM_DC_BOUNDS)))
                popt, pcov = curve_fit(model, tt[fit], uF[fit], p0=p0,
                                       bounds=((-np.inf, cfg.PERM_DC_BOUNDS[0]), (np.inf, cfg.PERM_DC_BOUNDS[1])))
                res = uF[fit] - model(tt[fit], *popt)
                tau_ac = autocorr_time(res)                # samples; the 5 k-step points are correlated
                se = np.sqrt(np.diag(pcov)) * np.sqrt(tau_ac)
                ss = np.sum((uF[fit] - np.mean(uF[fit])) ** 2)
                F['trace'] = dict(step=st, t_lj=tt, uF=uF, fit=fit, u_inf=float(popt[0]), u_inf_se=float(se[0]),
                                  Dc=float(popt[1]), Dc_se=float(se[1]), tau_ac=float(tau_ac),
                                  R2=float(1.0 - np.sum(res ** 2) / ss) if ss > 1e-30 else np.nan,
                                  model=lambda t, p=popt: model(t, *p),
                                  tau1=L0 ** 2 / (np.pi ** 2 * float(popt[1])))
            except (RuntimeError, ValueError) as e:
                print(f'  NOTE: thickness-trace fit failed ({e})')
        # steady thickness (the steady window of the stress snapshots) for M from dL
        ws = disp['bb_step'] >= P['halt_ts']
        if ws.sum() >= 4:
            Lss, lo, hi, *_ = block_bootstrap_ci(disp['bb_L'][ws], cfg.ci_level)
            F['L_ss'], F['L_ss_ci'] = float(Lss), (float(lo), float(hi))
        elif ws.any():
            F['L_ss'], F['L_ss_ci'] = float(np.mean(disp['bb_L'][ws])), (np.nan, np.nan)
    return F


def perm_M_estimates(cfg, P, F):
    """M = dP_ext L_ss/(2 |u_F|) from (bb) the steady thickness change L_0 - L_ss, (prof) the
    steady-profile parabola c = u_F/L_0 and (trace) the fitted asymptote u_F(inf).  dP_ext is
    the MEASURED piston difference (block-bootstrap CI) when the flux loader has it, else the
    applied one; every CI combines the relative half-widths in quadrature.  L_ss (the steady
    compressed thickness, the state the modulus describes) is the primary L (2026-10-02; L_0
    before); the L_0 variants are stored as M_L0.  Falls back to L_0 when L_ss is missing."""
    flux = P.get('flux') or {}
    if 'dP_meas' in flux:
        dP, dP_h, dP_src = flux['dP_meas'], 0.5 * (flux['dP_meas_hi'] - flux['dP_meas_lo']), 'measured'
    elif np.isfinite(flux.get('dP_app', np.nan)):
        dP, dP_h, dP_src = flux['dP_app'], 0.0, 'applied'
    else:
        return None
    if not dP > 0:
        return None
    L0 = F['L0']
    L_ss = F.get('L_ss', np.nan)
    Lm = float(L_ss) if np.isfinite(L_ss) else float(L0)              # the L in M = dP L/(2|u_F|)
    out = dict(dP=float(dP), dP_h=float(dP_h), dP_src=dP_src, L0=L0, L_ss=L_ss, L=Lm,
               L_src=('L_ss' if np.isfinite(L_ss) else 'L_0'))
    rdp = dP_h / dP

    def _m(uF, uF_h, tag, how):
        if not (np.isfinite(uF) and uF != 0):
            return None
        M = dP * Lm / (2.0 * abs(uF))
        rel = np.sqrt((uF_h / abs(uF)) ** 2 + rdp ** 2) if np.isfinite(uF_h) else rdp
        d = dict(M=float(M), lo=float(M * (1 - rel)), hi=float(M * (1 + rel)), uF=float(uF), uF_h=float(uF_h),
                 eps=float(abs(uF) / L0), how=how)
        if np.isfinite(L_ss):
            d['M_L0'] = float(M * L0 / Lm)
        return d

    L0_bb = F.get('L0_bb', L0)
    if np.isfinite(F.get('L_ss', np.nan)):
        h0 = 0.5 * (F['L0_bb_ci'][1] - F['L0_bb_ci'][0]) if np.isfinite(F.get('L0_bb_ci', (np.nan,))[0]) else 0.0
        hs = 0.5 * (F['L_ss_ci'][1] - F['L_ss_ci'][0]) if np.isfinite(F['L_ss_ci'][0]) else 0.0
        out['bb'] = _m(F['L_ss'] - L0_bb, np.sqrt(h0 ** 2 + hs ** 2), 'bb', 'steady thickness change L_0 - L_ss (bounding box, both)')
    out['prof'] = _m(F['c'] * L0, 0.5 * (F['c_hi'] - F['c_lo']) * L0, 'prof',
                     f"steady-profile parabola c zeta(2 - zeta), least squares over {F['c_npts']} points of {len(F['steady'])} steady snapshot(s)")
    T = F.get('trace')
    if T is not None:
        out['trace'] = _m(T['u_inf'], T['u_inf_se'] * _Z95, 'trace', 'thickness-trace asymptote u_F(inf), zero-IC series')
    lag = P.get('lag')
    if lag is not None and len(lag['idx']) >= 4:
        i0, i1 = lag['idx'][0], lag['idx'][-1]
        pf = lag['per_frame'][:, i1] - lag['per_frame'][:, i0]
        h = float(stats.t.ppf(0.5 + cfg.ci_level / 2, len(pf) - 1) * pf.std(ddof=1) / np.sqrt(len(pf))) if len(pf) >= 2 else np.nan
        out['lag'] = _m(float(lag['u'][i1] - lag['u'][i0]), h, 'lag',
                        f"per-atom pairing of {lag['n_frames']} steady frame(s) with the zero-flux reference, binned by the reference Z: "
                        f"u(Z = {lag['zc'][i1]:.0f}) - u(Z = {lag['zc'][i0]:.0f}) (rigid drop cancels, contact layer included)")
    out = {k: v for k, v in out.items() if v is not None}
    order = [cfg.PERM_M_PRIMARY] + [k for k in ('lag', 'prof', 'bb', 'trace') if k != cfg.PERM_M_PRIMARY]
    out['primary'] = next((k for k in order if k in out), None)
    return out


def add_perm_displacement(cfg, R, P, verbose=True):
    """P['disp'], P['Dc'] (fit_perm_Dc), P['M_perm'] (perm_M_estimates) and P['kappa_Dc']
    (D_c/M against the Darcy kappa of figure 9).  Missing files -> the keys are None and
    figures 12-13 are skipped."""
    say = print if verbose else (lambda *a, **k: None)
    P['disp'] = load_perm_disp(cfg, R, P)
    P['Dc'] = P['M_perm'] = P['kappa_Dc'] = None
    if P['disp'] is None:
        say('  NOTE: no disp_z_polymer file (or < 2 snapshots) -> D_c / M from the displacement skipped')
        return P
    F = fit_perm_Dc(cfg, R, P, P['disp'])
    if F is None:
        say('  NOTE: D_c fit not possible (no support plane, L_0 <= 0 or too few polymer bins)')
        return P
    P['Dc'] = F
    d = P['disp']
    fz = F['forcing']
    say(f"displacement: {len(d['ts'])} snapshots (windows of {int(d['stride'])} steps, first {int(d['ts'][0])});  "
        f"dP ramp: start step {int(fz['t0'])} -> full dP at {int(fz['t_end'])} ({fz['n']} increment(s); {fz['source']});  "
        f"displacement reference = {F['reset_mode']} (step {int(F['t_reset'])})\n"
        f"  L_0 = {F['L0']:.2f} sigma ({F['L0_source']}), pinned face z^P = {F['z_P']:.2f} (contact plane), "
        f"Z^P = {F['Z_P']:.2f} (reference position = z^P - u_0), z^F = {F['z_F']:.2f}; "
        f"fit domain {F['n_pts']} points ({'material' if F['lagr'] else 'Eulerian'} zeta in [{F['zf_lo']:.3f}, {F['zf_hi']:.3f}]), "
        f"{F['n_modes']} modes")
    if F['rigid'] == 'contact':
        say(f"  rigid shift u_0 (contact layer, {F['pin_nbins']} lowest bins): "
            + ', '.join(f"{v:+.2f}" for v in F['u0']) + f" sigma (bin spread {np.max(F['u0_sd']):.2f}); "
            f"steady u_0 = {F['u0_ss']:+.2f} sigma = the network's drop onto the plate -- subtracted before the fit")
    say(f"  D_c (profile fit, {len(F['early'])} snapshots{', window-averaged modes' if F['W'] > 0 else ''}) = "
        f"{F['Dc']:.4e} sigma^2/tau  (R^2 = {F['R2']:.3f});  tau_1 = L_0^2/(pi^2 D_c) = {F['tau1']:.0f} tau = "
        f"{F['tau1'] / cfg.dt_lj / 1e6:.2f}M steps;  production = {F['n_tau']:.2f} tau_1 "
        f"-> {100 * F['unrelaxed']:.1f}% of u_F still unrelaxed at the end (first mode)")
    T = F.get('trace')
    if T is not None:
        say(f"  D_c (thickness trace, zero-IC series) = {T['Dc']:.4e} +/- {T['Dc_se']:.1e} sigma^2/tau  "
            f"(u_F(inf) = {T['u_inf']:.2f} +/- {T['u_inf_se']:.2f} sigma, R^2 = {T['R2']:.3f}, tau_ac = {T['tau_ac']:.0f} samples);  "
            f"tau_1 = {T['tau1']:.0f} tau = {T['tau1'] / cfg.dt_lj / 1e6:.2f}M steps")
    say(f"  steady profile: c = u_F/L_0 = {F['c']:.4f} [{F['c_lo']:.4f}, {F['c_hi']:.4f}] (parabola R^2 = {F['parab_R2']:.3f})"
        + (f";  L_ss = {F['L_ss']:.2f} sigma (bounding box; dL = {F['L0_bb'] - F['L_ss']:.2f} from the bounding-box L_0 {F['L0_bb']:.2f}, eps_F = {(F['L0_bb'] - F['L_ss']) / F['L0']:.4f})" if 'L_ss' in F else ''))
    P['lag'] = perm_lagrangian_profile(cfg, R, P)             # per-atom steady profile (None without local trajectories)
    if P['lag'] is None:
        say('  NOTE: no local traj_ref / traj_stress -> no Lagrangian u_F; the parabola is the primary M')
    Mp = perm_M_estimates(cfg, P, F)
    P['M_perm'] = Mp
    if Mp:
        say(f"  M = dP_ext L_ss/(2 |u_F|)  with dP_ext = {Mp['dP']:.4f} ({Mp['dP_src']}), L_ss = {Mp['L']:.2f} ({Mp['L_src']}; L_0 = {Mp['L0']:.2f}):")
        for key in ('lag', 'prof', 'bb', 'trace'):
            v = Mp.get(key)
            if v:
                say(f"    {key:5s}{' *' if key == Mp.get('primary') else '  '} M = {v['M']:.4f} [{v['lo']:.4f}, {v['hi']:.4f}]  (u_F = {v['uF']:+.2f} sigma, eps_F = {v['eps']:.4f}"
                    + (f"; with L_0 instead -> {v['M_L0']:.4f}" if 'M_L0' in v else '') + f")  {v['how']}")
        Mref = Mp.get(Mp['primary']) if Mp.get('primary') else None
        if Mref:
            K = {'profile': dict(k=F['Dc'] / Mref['M'], lo=F['Dc'] / Mref['hi'], hi=F['Dc'] / Mref['lo'])}
            if T is not None:
                K['trace'] = dict(k=T['Dc'] / Mref['M'], lo=(T['Dc'] - _Z95 * T['Dc_se']) / Mref['hi'],
                                  hi=(T['Dc'] + _Z95 * T['Dc_se']) / Mref['lo'])
            P['kappa_Dc'] = K
            kd = (P.get('flux') or {}).get('k', {}).get('N_measured')
            say(f"  kappa = D_c/M (M = {Mp['primary']}, {Mref['M']:.4f}):  "
                + '   '.join(f"{k}: {v['k']:.4e} [{v['lo']:.4e}, {v['hi']:.4e}]" for k, v in K.items())
                + (f"   |  Darcy kappa = Q L/(A dP) (fig 9): {kd['k']:.4e}  ->  D_c(Darcy) = kappa M = {kd['k'] * Mref['M']:.4e}" if kd else ''))
    return P


def fig_perm_Dc(cfg, R, P):
    """Figure 12: (a) the deformation (u_z - u_0)(zeta, t)/L_0 of the snapshots, (b) data + the
    one-sided consolidation model (D_c), (c) the feed-face displacement u_F(t) from the
    thickness trace with the exact zero-IC series and the profile fit's own u(1, t)."""
    F = P.get('Dc')
    if F is None:
        print('D_c figure skipped (no fit)')
        return None
    zff = np.linspace(0.0, 1.0, 400)
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85])
    axl, axr, axt = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
    norm = Normalize(vmin=F['ts'][F['early'][0]], vmax=F['ts'][F['early'][-1]])
    cmap = plt.cm.viridis
    for i in F['early']:
        c = cmap(norm(F['ts'][i]))
        axl.plot(F['zl'][i], F['yl'][i], 'o-', color=c, ms=3, alpha=0.6)
        axr.plot(F['zl'][i], F['yl'][i], 'o', color=c, ms=3, alpha=0.35)
        axr.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=2.0)
    u_end = float(F['u_inf'](1.0)[0])
    contact = F['rigid'] == 'contact'
    ylab = r'$(u_z-u_0)/L_0$' if contact else r'$u_z/L_0$'
    xlab = (r'$\zeta=(z-u_z-Z^P)/L_0$  (material coordinate)' if F['lagr'] else r'$\zeta=(z-z^P)/L_0$')
    for ax in (axl, axr):
        ax._tri_has_colorbar = True
        ax.axhline(0.0, color='0.45', lw=1.6, label=(r'IC: $u_z(\zeta,0)=0$ at the ramp start' if F['reset_mode'] == 'ramp_start'
                                                    else r'reference: $u_z=0$ at the ramp end'))
        ax.plot(zff, F['u_inf'](zff), 'k:', lw=1.8, label=r'steady state ($t\to\infty$, fitted modes)')
        ax.plot(zff, F['c'] * zff * (2.0 - zff), '--', color='0.35', lw=1.4,
                label=(rf"steady-window parabola $c\,\zeta(2-\zeta)$, $c={sig(F['c'])}$  ($\to M$)" if ax is axl
                       else 'steady-window parabola'))
        ax.plot(0.0, 0.0, 's', color=WONG['blue'], ms=9, zorder=5, label=r'BC: $u_z(0,t)=u_0$ (support, pinned)')
        ax.plot([0.93, 1.0], [u_end, u_end], '-', color=WONG['vermillion'], lw=4, solid_capstyle='butt', zorder=5,
                label=r'BC: $\partial u_z/\partial\zeta\,|_{\zeta=1}=0$ (free feed face)')
        if contact and ax is axl:
            ax.plot([], [], ' ', label=rf"rigid drop $u_0={F['u0_ss']:+.2f}\,\sigma$ subtracted (contact layer)")
        ax.set(xlabel=flip_hint(R, xlab), ylabel=ylab, xlim=_zlim(R))
        ax.grid(alpha=0.3)
        smart_legend(ax, fontsize=10)
    axl.set_title(r'(a) ' + ylab + r' -- production snapshots (since the reset)'
                  + ('\n' + rf"face $Z^P={F['Z_P']:.1f}$ at zero flux $\to z^P={F['z_P']:.1f}$ under flow" if contact else ''), fontsize=13)
    axr.set_title(rf"(b) consolidation fit: $D_c={sig(F['Dc'])}\ \sigma^2/\tau$, $R^2={sig(F['R2'])}$, {F['n_modes']} modes"
                  + (r' (window-averaged)' if F['W'] > 0 else ''), fontsize=13)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=[axl, axr], fraction=0.015, pad=0.02).set_label('timestep')
    # ---- (c) feed-face displacement vs time ----
    T = F.get('trace')
    if T is not None:
        axt.plot(T['step'], T['uF'] / F['L0'], '-', color='0.6', lw=0.7, alpha=0.6,
                 label=r'$u_F(t)=L_{bb}(t)-L_0$  (bounding box, 5k-step samples)')
        axt.plot(T['step'], rolling_mean(T['uF'], cfg.roll_win) / F['L0'], '-', color='0.25', lw=1.8,
                 label=f'rolling mean ({cfg.roll_win})')
        tt = np.linspace(0.0, T['t_lj'][-1], 400)
        axt.plot(F['t_onset'] + tt / cfg.dt_lj, T['model'](tt) / F['L0'], '-', color=WONG['vermillion'], lw=2.4,
                 label=(rf"zero-IC series fit: $D_c={sig(T['Dc'])}$, $u_F(\infty)/L_0={sig(T['u_inf'] / F['L0'])}$, "
                        rf"$R^2={sig(T['R2'])}$"))
        if T['fit'].sum() < len(T['fit']):
            axt.axvline(T['step'][T['fit']][-1], color=WONG['vermillion'], ls=':', lw=1.2, alpha=0.7)
        if cfg.PERM_TRACE_SKIP > 0:
            axt.axvline(F['t_onset'] + cfg.PERM_TRACE_SKIP, color=WONG['vermillion'], ls='-.', lw=1.2, alpha=0.7,
                        label=f'trace fit starts (PERM_TRACE_SKIP = {cfg.PERM_TRACE_SKIP / 1e6:.1f}M steps)')
    tt = np.linspace(0.0, F['t_lj'][-1], 400)
    axt.plot(F['t_onset'] + tt / cfg.dt_lj, [float(F['u_phys'](1.0, t)[0]) for t in tt], '--',
             color=WONG['blue'], lw=2.0, label=rf"profile fit at $\zeta=1$ (instantaneous): $D_c={sig(F['Dc'])}$")
    # the top fitted point of each snapshot: the deformation there (u_z - u_0), at the window centre,
    # plus the fitted deformation already present at the reset (zero for a ramp-start reset)
    ztop = np.array([F['zl'][i][-1] for i in range(len(F['ts']))])
    ytop = np.array([F['yl'][i][-1] for i in range(len(F['ts']))])
    off = np.array([float(F['u_phys'](zt, F['t_reset_lj'])[0]) for zt in ztop])
    axt.plot(F['ts'] - 0.5 * F['stride'], ytop + off, 'o', color=WONG['blue'], ms=7, zorder=5,
             label=(rf"top fitted point of each snapshot ($\zeta={ztop.min():.2f}$–${ztop.max():.2f}$, "
                    + (r'$u_z-u_0$, ' if contact else '') + 'window centre'
                    + (rf"; $+{np.max(np.abs(off)):.4f}$ fitted offset at the reset)" if np.max(np.abs(off)) > 1e-4 else ')')))
    fz = F['forcing']
    if fz['t_end'] - fz['t0'] > 0.02 * (float(F['ts'][-1]) - fz['t0']):
        axt.axvspan(fz['t0'], fz['t_end'], color='0.5', alpha=0.12, label=rf"$\Delta P$ ramp ({fz['n']} increments)")
    axt.axvspan(P['halt_ts'], float(F['ts'][-1]), color=WONG['green'], alpha=0.10, label='steady window')
    axt.axvline(F['t_onset'], color='k', ls=':', lw=1, alpha=0.6)
    axt.axhline(0.0, color='k', ls=':', lw=1, alpha=0.4)
    axt.set_xlabel('time step')
    axt.set_ylabel(r'$u_F/L_0$')
    axt.set_title(r'(c) feed-face displacement: thickness trace vs the two fits', fontsize=13)
    axt.grid(alpha=0.3)
    smart_legend(axt, fontsize=9)
    fig.suptitle(f'Cooperative diffusivity, permeation drive  |  {cfg.sim_name}  |  '
                 rf"$L_0={sig(F['L0'])}\,\sigma$, $z^P={F['z_P']:.1f}$"
                 + (rf" ($u_0={F['u0_ss']:+.2f}\,\sigma$)" if contact else '')
                 + rf", ramp {(fz['t_end'] - fz['t0']) / 1e6:.2f}M steps, "
                 rf"drive $={F['n_tau']:.2f}\,\tau_1$ "
                 rf"($\tau_1=L_0^2/\pi^2 D_c={F['tau1']:.0f}\,\tau$, profile fit)", fontsize=12, fontweight='bold')
    return _save(fig, cfg, 'perm_Dc_consolidation_fit')


def fig_perm_M(cfg, R, P):
    """Figure 13: M = dP_ext L_ss/(2 |u_F|) -- one point per reading of the feed-face
    displacement (steady-profile parabola, thickness change, trace asymptote), with CIs;
    the primary (PERM_M_PRIMARY) is marked, and the caption below the axes says why the
    bounding-box readings differ from the profile."""
    Mp = P.get('M_perm')
    if not Mp:
        print('M figure skipped (no displacement fit or no dP)')
        return None
    F = P.get('Dc') or {}
    prim = Mp.get('primary')
    order = (('lag', 'per-atom (Lagrangian)' + '\n' + r'$u(Z^F_{\rm mat}) - u(Z^P_{\rm mat})$, zero-flux ref', WONG['reddishpurple'], '^'),
             ('prof', 'steady profile' + '\n' + r'parabola $c\,\zeta(2-\zeta)$', WONG['vermillion'], 's'),
             ('bb', 'thickness change' + '\n' + r'$\Delta L = L_0 - L_{ss}$ (bounding box)', WONG['blue'], 'o'),
             ('trace', 'thickness trace' + '\n' + r'asymptote $u_F(\infty)$ (bounding box)', WONG['green'], 'D'))
    show = [o for o in order if o[0] in Mp]
    fig, ax = plt.subplots(figsize=(10, 7.4), constrained_layout=True)
    ci = int(cfg.ci_level * 100)
    for k, (key, lab, col, mk) in enumerate(show):
        v = Mp[key]
        ax.errorbar([k], [v['M']], yerr=[[v['M'] - v['lo']], [v['hi'] - v['M']]], fmt=mk, ms=13, color=col,
                    capsize=8, lw=2.5, mfc=(col if key == prim else 'white'), mew=2.5,
                    label=(('PRIMARY  ' if key == prim else '') + f"{lab.splitlines()[0]}:  $M = {sig(v['M'])}$\n{ci}% CI [{sig(v['lo'])}, {sig(v['hi'])}]"
                           rf"  ($u_F = {v['uF']:+.1f}\,\sigma$, $\varepsilon_F = {v['eps']:.3f}$)"))
        ax.axhline(v['M'], color=col, ls='--', lw=1.0, alpha=0.4)
    ax.set_xticks(range(len(show)))
    ax.set_xticklabels([o[1] + ('\n(primary)' if o[0] == prim else '') for o in show], fontsize=12)
    ax.set_ylabel(r'$M$  (LJ units)')
    ax.set_xlim(-0.6, len(show) - 0.4)
    ax.grid(axis='y', alpha=0.3)
    if np.isfinite(Mp.get('L_ss', np.nan)):
        lline = (rf"$L_{{ss}}={sig(Mp['L_ss'])}\,\sigma$ (steady), $L_0={sig(Mp['L0'])}\,\sigma$ (onset; "
                 rf"with $L_0$ instead every $M$ is $\times{Mp['L0'] / Mp['L_ss']:.3f}$)")
    else:
        lline = rf"$L_0={sig(Mp['L0'])}\,\sigma$ (onset; no steady window, used in place of $L_{{ss}}$)"
    ax.set_title(r'Longitudinal modulus from the permeation drive:  $M = \Delta P_{\mathrm{ext}}\,L_{ss}\,/\,(2\,|u_F|)$'
                 + '\n' + rf"$\Delta P_{{\mathrm{{ext}}}}={sig(Mp['dP'])}$ ({Mp['dP_src']}), " + lline + f'\n{cfg.sim_name}', fontsize=12)
    smart_legend(ax, fontsize=11)
    # ---- caption: why the two families of u_F differ (so that the choice of primary is not forgotten) ----
    if 'prof' in Mp and 'bb' in Mp:
        up, ub = Mp['prof']['uF'], Mp['bb']['uF']
        u0 = F.get('u0_ss', np.nan)
        dM = 100 * (Mp['bb']['M'] / Mp['prof']['M'] - 1)
        cap = (f"Two readings of u_F.  The parabola reads the deformation of the network alone: the binned displacement with the rigid drop\n"
               f"u_0 = {u0:+.2f} \u03c3 subtracted (the network's fall onto the plate when the flow closes the zero-flux gap), in the material coordinate,\n"
               f"pinned at the contact plane  ->  u_F = {up:+.1f} \u03c3.   The bounding-box thickness runs from the lowest to the highest polymer bead:\n"
               f"its bottom end is a tail bead that drops with the network onto the plate and its top end the extreme bead a few \u03c3 above the\n"
               f"face, so \u0394L_bb (u_F = {ub:+.1f} \u03c3) carries {abs(ub - up):.1f} \u03c3 of the drop / tail statistics into |u_F| and the bounding-box M "
               f"({'bb and trace' if 'trace' in Mp else 'bb'}) come out\n"
               f"{dM:+.0f}% relative to the parabola.   Primary = {prim} (PERM_M_PRIMARY): it feeds \u03ba = D_c/M and the q(t) check of figure 14."
               + (f"\nThe per-atom (Lagrangian) reading pairs the steady frames with the zero-flux reference atom by atom, so the rigid drop cancels and the\n"
                  f"contact layer's own compaction counts: u_F = {Mp['lag']['uF']:+.1f} \u03c3.  It is the one the compression sweep reproduces (figure 15)." if 'lag' in Mp else ''))
        fig.text(0.01, -0.02, cap, ha='left', va='top', fontsize=10, family='serif',
                 bbox=dict(boxstyle='round', fc='0.97', ec='0.75'))
    return _save(fig, cfg, 'perm_M')


def _perm_Hdot(mu, Dc, L, t, forcing):
    """dH_k/dt = sum_j w_j lam_k exp(-lam_k (t - t_j)) over the increments already applied: (len(t), n_modes)."""
    mu = np.atleast_1d(np.asarray(mu, float))
    lam = (mu / L) ** 2 * Dc
    t = np.atleast_1d(np.asarray(t, float))
    tj, wj = forcing
    dt = t[:, None] - tj[None, :]
    on = dt > 0
    dtp = np.where(on, dt, 0.0)
    out = np.empty((len(t), len(lam)))
    for k, lk in enumerate(lam):
        out[:, k] = lk * np.sum(wj[None, :] * on * np.exp(-lk * dtp), axis=1)
    return out


def perm_flux_model(cfg, R, P):
    """The q(t) CHECK (2026-09-28): the permeate flux the displacement fits imply.  Integrating
    du/dt = q + D_c u'' over the membrane with u'(z^F) = 0 and u'(z^P) = -dP(t)/M gives
        q(t) L = d/dt int u dz - D_c dP(t)/M   ->   Q(t) = -A q = A [ (D_c/M) dP(t)/L_0 - d<u>/dt ],
    <u> the membrane-mean displacement.  The steady level A (D_c/M) dP/L_0 is figure 9's Darcy
    check in flux form (kappa = D_c/M); the transient is the solvent the compaction expels.
    Profile fit: <u> = L_0 sum_k a_k H_k(t) (int phi_k dzeta = 1); trace fit (zero-IC shape):
    <u> = u_F sum_k 4/(k pi)^2 H_k(t).  dP(t) is the applied staircase; M the primary M (PERM_M_PRIMARY).
    Nothing here is fitted to the flux.  Returns a dict or None."""
    F, Mp = P.get('Dc'), P.get('M_perm')
    if F is None or not Mp or not Mp.get('primary'):
        return None
    fz = F['forcing']
    forcing = (fz['tj'], fz['wj'])
    L0, A, M = F['L0'], P['area'], Mp[Mp['primary']]['M']
    dP_full = fz['full'] if np.isfinite(fz.get('full', np.nan)) else Mp['dP']
    T_end = max(float(F['t_lj'][-1]), float(F['trace']['t_lj'][-1]) if F.get('trace') else 0.0)
    tt = np.linspace(0.0, T_end, 1500)
    dP_t = dP_full * np.array([fz['wj'][fz['tj'] < t].sum() for t in tt])
    out = dict(t_lj=tt, step=fz['t0'] + tt / cfg.dt_lj, dP=dP_t, M=M, L0=L0, A=A, dP_full=dP_full)
    dudt = L0 * (_perm_Hdot(F['mu'], F['Dc'], L0, tt, forcing) @ F['A'])
    out['Q_prof'] = A * ((F['Dc'] / M) * dP_t / L0 - dudt)
    out['Q_prof_ss'] = float(A * (F['Dc'] / M) * dP_full / L0)
    T = F.get('trace')
    if T is not None:
        mu_all = np.arange(1, _PERM_TRACE_MODES + 1) * np.pi
        dudt = T['u_inf'] * (_perm_Hdot(mu_all, T['Dc'], L0, tt, forcing) @ (4.0 / mu_all ** 2))
        out['Q_trace'] = A * ((T['Dc'] / M) * dP_t / L0 - dudt)
        out['Q_trace_ss'] = float(A * (T['Dc'] / M) * dP_full / L0)
    return out


def fig_perm_v0_check(cfg, R, P):
    """TEMPORARY diagnostic (2026-10-09; delete once the D_c discrepancy is resolved).  The permeation counterpart of
    fig_v0_check: (a) D_c from the displacement profile fit (v0(t) fixed by the modes vanishing at the support), from
    the same fit with the per-snapshot rigid translation projected out (the strain-PDE fit), the thickness-trace fit
    and the Darcy kappa M; (b) the rigid motion the standard fit removes by measurement (the contact layer u_0(t))
    against what the strain fit would add to it; (c) v0 at the support: the measured permeate flux Q/A against the
    flux the displacement fit implies (perm_flux_model) -- v0(t) is NOT constant during a permeation transient, and
    this is the direct test of whether the fit's v0 is the real one."""
    F, disp = P.get('Dc'), P.get('disp')
    if F is None or disp is None:
        print('permeation v0 check skipped (no displacement fit)')
        return None
    F2 = fit_perm_Dc(cfg, R, P, disp, free_offset=True)
    if F2 is None:
        print('permeation v0 check skipped (strain fit failed)')
        return None
    fig, (axa, axb, axc) = plt.subplots(1, 3, figsize=(21, 6.5), constrained_layout=True)
    # (a) the D_c readings
    names = [r'profile fit' + '\n' + r'($v^0$ from the modes)', r'strain-PDE fit' + '\n' + r'($v^0$ projected out)']
    vals, cols = [F['Dc'], F2['Dc']], [WONG['blue'], WONG['vermillion']]
    T = F.get('trace')
    if T is not None:
        names.append('thickness trace\n(zero-IC series)'); vals.append(T['Dc']); cols.append(WONG['skyblue'])
    kd = (P.get('flux') or {}).get('k', {}).get('N_measured') or (P.get('flux') or {}).get('k', {}).get('N_applied')
    Mp = P.get('M_perm')
    if kd and Mp and Mp.get('primary'):
        names.append(r'Darcy $\kappa\,M$' + f"\n(M = {Mp[Mp['primary']]['M']:.3f})"); vals.append(kd['k'] * Mp[Mp['primary']]['M']); cols.append(WONG['green'])
    axa.bar(np.arange(len(vals)), vals, color=cols, alpha=0.85)
    for i, v in enumerate(vals):
        axa.text(i, v, f' {sig(v)}', ha='center', va='bottom', fontsize=12)
    axa.set_xticks(np.arange(len(vals)))
    axa.set_xticklabels(names, fontsize=10)
    axa.set_ylabel(r'$D_c$  ($\sigma^2/\tau$)')
    axa.set_title(rf"(a) $D_c$: with $v^0$ {sig(F['Dc'])} vs without {sig(F2['Dc'])} (ratio {F2['Dc'] / F['Dc']:.3f})", fontsize=14)
    axa.grid(alpha=0.3, axis='y')
    # (b) rigid motion: the measured contact-layer u_0(t), and <u> over the fit points (after u_0) -- data vs the
    #     displacement model, whose own rigid part is what the support-pinned modes allow
    t = F['t_lj']
    forcing = (F['forcing']['tj'], F['forcing']['wj'])
    um_data = np.array([np.mean(F['yl'][i]) for i in range(len(t))]) * F['L0']
    um_model = np.array([np.mean(_perm_phi(F['zl'][i], F['mu']) @ (F['A'] * (_perm_H(F['mu'], F['Dc'], F['L0'], t[i], forcing, F['W'])[0] - F['H_reset'])))
                         for i in range(len(t))]) * F['L0']
    off = um_data - um_model
    axb.plot(t, F['u0'], '-', color='k', lw=2, label=r'contact layer $u_0(t)$ (measured, subtracted before both fits)')
    axb.plot(t, um_data, '.', color='0.45', ms=4, label=r'data: $\langle u_z\rangle$ over the fit points (after $u_0$)')
    axb.plot(t, um_model, '-', color=WONG['blue'], lw=2.2, label='displacement model (its own rigid part)')
    axb.plot(t, off, '-', color=WONG['vermillion'], lw=1.8, label='difference = translation the modes cannot represent')
    axb.axhline(0, color='k', lw=0.8)
    axb.axvline((F['t_reset'] - F['t_onset']) * cfg.dt_lj, color='0.5', ls=':', lw=1)
    axb.set(xlabel=r'time from the ramp start  ($\tau$)', ylabel=r'rigid displacement  ($\sigma$)')
    axb.set_title('(b) rigid motion of the network', fontsize=14)
    axb.grid(alpha=0.3)
    smart_legend(axb, fontsize=10)
    # (c) v0 at the support: measured permeate flux / A vs the flux the displacement fit implies
    Q = perm_flux_model(cfg, R, P)
    Fx = P.get('flux') or {}
    A = P.get('area', np.nan)
    if Fx.get('Q') is not None and Fx.get('step') is not None:
        tq = (Fx['step'] - F['t_onset']) * cfg.dt_lj
        ok = tq > 0
        axc.plot(tq[ok], Fx['Q'][ok] / A, '.', color='0.6', ms=3, label=r'measured: permeate-piston $Q/A$ (block averages)')
        qs = Fx.get('Q_N')
        if qs:
            axc.axhline(qs['mean'] / A, color='k', ls='--', lw=1.5, label=f"measured steady $Q/A$ ({qs['mean'] / A:.2e})")
    if Q is not None:
        axc.plot(Q['t_lj'], Q['Q_prof'] / A, '-', color=WONG['blue'], lw=2.2, label=r'implied by the displacement fit: $v^0 = (D_c/M)\,\Delta P/L_0 - d\langle u\rangle/dt$')
        if 'Q_trace' in Q:
            axc.plot(Q['t_lj'], Q['Q_trace'] / A, '-', color=WONG['skyblue'], lw=1.6, label='implied by the trace fit')
    axc.set(xlabel=r'time from the ramp start  ($\tau$)', ylabel=r'$v^0 = Q/A$  ($\sigma/\tau$)')
    axc.set_title(r'(c) $v^0(t)$ at the support: measured vs implied', fontsize=14)
    axc.grid(alpha=0.3)
    if axc.get_legend_handles_labels()[0]:
        smart_legend(axc, fontsize=10)
    else:
        axc.text(0.5, 0.5, 'no flux data / no M (perm_flux_model)', ha='center', va='center', transform=axc.transAxes)
    fig.suptitle(f'TEMPORARY v0 check (permeation)  |  {cfg.sim_name}', fontsize=12, fontweight='bold')
    print(f"  v0 check (permeation): D_c profile fit {F['Dc']:.4e}  |  strain-PDE fit {F2['Dc']:.4e}  ->  ratio {F2['Dc'] / F['Dc']:.3f}"
          + (f"  |  trace {T['Dc']:.4e}" if T is not None else '') + f";  rms translation data - displacement model {np.sqrt(np.mean(off ** 2)):.3f} sigma")
    return _save(fig, cfg, 'perm_v0_check')


def fig_perm_flux_check(cfg, R, P):
    """Figure 14 (check): (a) the model-implied Q(t) of perm_flux_model against the measured
    piston-velocity trace and the N_perm-slope Q; (b) its time integral against the measured
    permeate bead count."""
    Q = perm_flux_model(cfg, R, P)
    Fx = P.get('flux') or {}
    if Q is None or Fx.get('Q') is None and Fx.get('N') is None and Fx.get('count') is None:
        print('flux-check figure skipped (no displacement fit or no flux data)')
        return None
    fig, (axQ, axN) = plt.subplots(1, 2, figsize=(18, 6), constrained_layout=True)
    st = Fx.get('step')
    if Fx.get('Q') is not None:
        axQ.plot(st, Fx['Q'], '-', color=WONG['green'], lw=0.9, alpha=0.3, label=r'measured: $A\,dz_{\rm perm}/dt$ (block-averaged)')
        axQ.plot(st, rolling_mean(Fx['Q'], cfg.roll_win), '-', color=WONG['green'], lw=2.2, label=f'measured: rolling mean ({cfg.roll_win})')
    if 'Q_N' in Fx:
        NS = Fx['Q_N']
        axQ.axhline(NS['mean'], color='k', ls='--', lw=1.8, label=rf"measured steady $Q$ ($N_{{\rm perm}}$ slope) = {fmt_val_unc(NS['mean'], NS['se'])}")
        axQ.axhspan(NS['lo'], NS['hi'], color='k', alpha=0.08)
    for key, col, lab in (('Q_prof', WONG['blue'], 'profile fit'), ('Q_trace', WONG['vermillion'], 'trace fit')):
        if key in Q:
            axQ.plot(Q['step'], Q[key], '-', color=col, lw=2.2,
                     label=rf"model, {lab}: $A[(D_c/M)\,\Delta P(t)/L_0 - d\langle u_z\rangle/dt]$;  steady $= {sig(Q[key + '_ss'])}$")
            axQ.axhline(Q[key + '_ss'], color=col, ls=':', lw=1.2, alpha=0.7)
    qref = abs(Fx['Q_N']['mean']) if 'Q_N' in Fx else abs(Q['Q_prof_ss'])
    axQ.set_ylim(-1.0 * qref, 4.0 * qref)                  # the model's 1/sqrt(t) onset spike is clipped on purpose
    axQ.axvline(Q['step'][0], color='k', ls=':', lw=1, alpha=0.6)
    axQ.set_xlabel('time step')
    axQ.set_ylabel(r'$Q_{\rm perm}$  ($\sigma^3/\tau$)')
    axQ.set_title(r'(a) permeate flux implied by the displacement fits vs measured (nothing fitted to the flux)', fontsize=12)
    axQ.grid(alpha=0.3)
    smart_legend(axQ, fontsize=9)
    # ---- (b) cumulative beads ----
    cnt = Fx.get('count')
    if cnt is not None:
        n_st, n_val = cnt[:, 0], cnt[:, 1]
    elif Fx.get('N') is not None:
        n_st, n_val = st, Fx['N']
    else:
        n_st = n_val = None
    rho0 = R.get('rho_s0', np.nan)
    if n_st is not None:
        axN.plot(n_st, n_val - n_val[0], '-', color='0.3', lw=1.8, label=r'measured $\Delta N_{\rm perm}$')
    if np.isfinite(rho0):
        from scipy.integrate import cumulative_trapezoid
        for key, col, lab in (('Q_prof', WONG['blue'], 'profile fit'), ('Q_trace', WONG['vermillion'], 'trace fit')):
            if key in Q:
                Nm = rho0 * cumulative_trapezoid(Q[key], Q['t_lj'], initial=0.0)
                if n_st is not None:                        # both start at zero at the first count sample
                    Nm = Nm - np.interp(n_st[0], Q['step'], Nm)
                axN.plot(Q['step'], Nm, '-', color=col, lw=2.2, label=rf'model, {lab}: $\rho_{{s,0}}\int Q\,dt$')
    axN.set_xlabel('time step')
    axN.set_ylabel(r'$\Delta N_{\rm perm}$  (beads)')
    axN.set_title(r'(b) cumulative permeate: measured bead count vs the integrated model flux', fontsize=12)
    axN.grid(alpha=0.3)
    smart_legend(axN, fontsize=10)
    fig.suptitle(f"$q(t)$ check  |  {cfg.sim_name}  |  $M = {sig(Q['M'])}$ (thickness change), "
                 rf"$L_0 = {sig(Q['L0'])}\,\sigma$, $\Delta P = {sig(Q['dP_full'])}$", fontsize=12, fontweight='bold')
    return _save(fig, cfg, 'perm_flux_check')


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
    F2, Mp = P.get('Dc'), P.get('M_perm')
    if F2 is not None:
        T = F2.get('trace')
        print(f"  D_c (profile fit) = {F2['Dc']:.4e} sigma^2/tau (R^2 {F2['R2']:.3f}; production {F2['n_tau']:.2f} tau_1)"
              + (f";  D_c (thickness trace) = {T['Dc']:.4e} +/- {T['Dc_se']:.1e} (R^2 {T['R2']:.3f})" if T else ''))
    if Mp:
        print('  M = dP L_ss/(2|u_F|): ' + '  '.join(f"{k}{'*' if k == Mp.get('primary') else ''}: {v['M']:.4f} [{v['lo']:.4f}, {v['hi']:.4f}]"
                                                 for k, v in Mp.items() if isinstance(v, dict)) + '   (* = primary)')
    K = P.get('kappa_Dc')
    if K:
        kd = F.get('k', {}).get('N_measured')
        print('  kappa = D_c/M: ' + '  '.join(f"{k}: {v['k']:.4e} [{v['lo']:.4e}, {v['hi']:.4e}]" for k, v in K.items())
              + (f"   vs Darcy kappa (N slope, measured dP) {kd['k']:.4e} [{kd['lo']:.4e}, {kd['hi']:.4e}]" if kd else ''))
    S = P.get('psd')
    if S is not None:
        print(f"  geometric porosity (r_probe {cfg.PSD_R_PROBE}), interior: {fmt_mu(S['por'][0][P['interior']])};  "
              f"mean pore diameter {fmt_mu(S['d_mean'][0][P['interior']])} sigma")
    elif 'Q_N' not in F:
        print('  no steady flux (permeation / permeate_count files missing or too short)')
