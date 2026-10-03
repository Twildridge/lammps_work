"""bulk.py -- shared analysis library for the compress_slab (bulk modulus) notebooks (2026-09-29).

Used by  scripts/bulk_modulus_analysis_single.ipynb  (one volumetric-strain level)  and
         scripts/bulk_modulus_analysis_sweep.ipynb   (every level of a sweep).

The bulk counterpart of lib/triaxial.py, laid out section for section like it so the two
notebooks read like the compression ones and -- more important -- so K is estimated the way
the compression notebooks estimate M:

  * geometry: an ISOLATED gel loaded on its six faces by solvent-transparent plates, the bath
    held at P_target by two wet pistons (compress_slab.lmp since 2026-09-29);
  * every profile is taken in a CORE COLUMN through the gel (bath | gel | bath), along z, x and
    y; the deck divides by the live column cross-section, so the files hold true stresses;
  * Terzaghi: sigma' = sigma^t - p_pore, p_pore = the total stress in the BATH bins of the same
    column (one scalar per component and snapshot), on both sides of the gel;
  * osmotic pressure Pi = (1/3)(sigma'_xx + sigma'_yy + sigma'_zz) in the profiles' sign
    (compression positive, so this is -(1/3) tr sigma' in the mechanical sign);
  * M_network <-> K_network = (<Pi>_interior,plateau - Pi_ref) / eps_vol, one estimate per
    profile axis; the headline value is their mean;
  * M_piston  <-> K_plates  = (<P_plates>_plateau - P_plates,ref) / eps_vol, P_plates = the mean
    over the six plates of the compressive load / face area, read over the longest drift-free
    trailing window with the same circular block bootstrap;
  * every K is an INCREMENT from the measured eps_vol = 0 reference window (K_SUBTRACT_REF),
    i.e. a slope, like M;
  * eps_vol is the MEASURED plate strain 1 - prod(1 - closure_i / L0_i), L0_i = seated plate gap
    - 2 contact_gap (the edge of gel the seated plates enclose), averaged over the hold (the
    plates are frozen there); the network's own Rg and BB volumetric strains are reported next
    to it and K_STRAIN selects the denominator.  The plates load the network through its
    sparse outer shell, so the network strain can fall short of the plate strain: K_network
    is printed for both (K_net and K_net(Rg)) and a gap between them is that shortfall;
  * unrelaxed-hold systematic (RELAX_SYS), isotropy of the network stress (the deviator must
    vanish under isotropic loading: there is no G in this measurement);
  * D_c from the THREE-DIMENSIONAL consolidation fit (fit_Dc_cube): the held-cube version of
    the triaxial held-slab fit -- every face drained at the bath pressure (the feed BC of the
    compression fit), zero displacement at the centre by symmetry -- so u_a along axis a
    relaxes in sin(2 l pi zeta_a) times the transverse cos(2 m pi zeta_b) cos(2 n pi zeta_c)
    modes, all at 4 pi^2 D_c (l^2/L_a^2 + m^2/L_b^2 + n^2/L_c^2); one D_c is fitted to the
    three column profiles at once; kappa = D_c / M; the wet-piston bath check and the
    solvent expelled -- as in the triaxial notebooks;
  * network anisotropy: sigma'_ii / Pi in the gel interior over every hold (fig_ratio) and
    the whole-run polymer partial-stress ratios sig_p,xx / sig_p,zz, sig_p,yy / sig_p,zz of
    stress_aniso (fig_anisotropy_run; the same file slab_with_support writes, so the free
    swelling can be drawn next to the compression);
  * solvent volume fractions in the core columns along x, y and z: mass fraction from the
    density profiles, Voronoi and lambda-calibrated Voronoi from the trajectories
    (add_volume_fractions, one tessellation per frame serves the three axes);
  * closure: with M from the triaxial run (M_REF), G = (3/4)(M - K) and kappa = D_c / M.

Shared machinery (style, readers, block bootstrap, plateau window, tail fit, Terzaghi split,
the consolidation fit, smart_legend, the Expanse puller and the PROFILE FIGURES) is
triaxial.py's; nothing statistical is re-implemented here.  The profile figures are
triaxial's own functions, handed a per-axis "view" of the data (see view()).

Layout of this file
    0. style                         (triaxial's)
    1. Config                        knobs, path builders (per-axis file names)
    2. file readers                  plate tables, per-axis profile stacks
    3. reference state               load_reference   (eps_vol = 0, shared by every level)
    4. one strain level              load_level       (stresses, Terzaghi, plates, K, D_c)
    5. Expanse sync                  sync_from_expanse (via triaxial.sync_pull)
    6. per-axis views                view, axis context for triaxial's figures
    7. figures, single level         fig_strain ... fig_closure
    8. figures, sweep                fig_*_sweep
    9. summaries                     print_summary, print_hold_check
   10. D_c of the cube               fit_Dc_cube (3-D held-cube modes, one D_c for the three axes)
   11. anisotropy                    fig_ratio, fig_anisotropy_run (+ sweep)
   12. volume fractions              add_volume_fractions (column Voronoi + calibration), fig_volfrac_evolution
"""
from __future__ import annotations

import sys
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import hashlib

import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import triaxial as tri  # noqa: E402
import volfrac       # noqa: E402  (scripts/lib/volfrac.py -- Voronoi + lambda calibration)

# ===========================================================================
#  0. STYLE  (triaxial's palette and rcParams)
# ===========================================================================
WONG, EVO_CMAP = tri.WONG, tri.EVO_CMAP
level_color, setup_style = tri.level_color, tri.setup_style
sig, fmt_step, fmt_val_unc, fmt_ci, fmt_mu = tri.sig, tri.fmt_step, tri.fmt_val_unc, tri.fmt_ci, tri.fmt_mu
mean_ci, block_bootstrap_ci, plateau_window, relax_tail = (tri.mean_ci, tri.block_bootstrap_ci,
                                                           tri.plateau_window, tri.relax_tail)
smart_legend, annotate_box, robust_ylim, rolling_mean = (tri.smart_legend, tri.annotate_box,
                                                         tri.robust_ylim, tri.rolling_mean)
COMPONENTS = tri.COMPONENTS                       # ('zz', 'xx', 'yy')
AXI = {'x': 0, 'y': 1, 'z': 2}
PLATES = ('zlo', 'zhi', 'xlo', 'xhi', 'ylo', 'yhi')   # column order of every plate_* file
AX_COLOR = {'x': WONG['blue'], 'y': WONG['vermillion'], 'z': WONG['green']}
PLATE_STYLE = {'zlo': (WONG['green'], '--'), 'zhi': (WONG['green'], '-'), 'xlo': (WONG['blue'], '--'),
               'xhi': (WONG['blue'], '-'), 'ylo': (WONG['vermillion'], '--'), 'yhi': (WONG['vermillion'], '-')}


# ===========================================================================
#  1. CONFIG
# ===========================================================================
@dataclass
class Config(tri.Config):
    """Every knob of the bulk analysis.  Inherits triaxial.Config (same names, same meaning:
    plateau_frac, wall_margin, res_wall_margin, RELAX_SYS, SS_FIT_NPTS, DC_*, ...); COMP_LEVELS
    holds the cumulative VOLUMETRIC strain targets as strings, = STRAIN_TARGETS in the .batch."""
    AXES: tuple = ('z', 'x', 'y')         # profile axes analysed; the first is the one the profile figures draw by default
    contact_gap: float = 1.12             # compress_slab.lmp contact_gap: seated plate plane to the gel face
    core_margin: float = 6.0              # compress_slab.lmp core_margin: the core column stops this far inside the plates
    K_SUBTRACT_REF: bool = True           # K = (stress - its own eps_vol = 0 reading) / eps_vol for BOTH estimators
    K_STRAIN: str = 'plate'               # eps_vol in the denominator: 'plate' (measured plate closure, the applied
                                          # strain), 'rg' / 'bb' (the network's own Rg / bounding-box volume), or
                                          # 'nominal' (the level's target)
    PLATE_AREA: str = 'gap'               # plate pressure = load / face area.  'gap': the face enclosed by the other two
                                          # plate pairs, (gap_a - 2 contact_gap)(gap_b - 2 contact_gap), constant over a
                                          # hold; 'bb': the live bounding-box face of the network (follows its extreme
                                          # beads, so it fluctuates)
    M_REF: object = None                  # longitudinal modulus M of the SAME gel from the triaxial notebooks: gives
                                          # G = (3/4)(M - K) and kappa = D_c / M.  None -> both skipped
    G_REF: object = None                  # shear modulus to compare that G with (shear or triaxial notebook); optional
    SYNC_TRAJ: bool = True                # pull traj_ref / traj_stress too (large): add_volume_fractions tessellates them
    DC_N_T: int = 0                       # transverse modes 0..DC_N_T per transverse direction in the held-cube fit (DC_N_MODES
                                          # = the longitudinal modes l along the profile axis, as for the slab).  0 = the pure
                                          # modes (l, 0, 0) only: identifiable exactly like the slab fit (one rate per shape).
                                          # 1 adds the mixed modes (l, m, n) that carry the transverse drainage into the core
                                          # column; they share the shape sin(2 l pi zeta) with the pure ones, so with free
                                          # amplitudes D_c is then pinned only through the multi-l structure (a single-rate
                                          # signal is degenerate: D_c and D_c/2 fit alike) -- a cross-check, not the default
    HOST: str = 'expanse'                 # 'expanse' | 'pod' | 'bridges': which cluster folder RUNS_ROOT defaults to
                                          # (the sync cell logs in to Expanse only; copy Pod / Bridges output by hand)

    def __post_init__(self):
        assert self.COMP_LEVELS, 'COMP_LEVELS is empty -- list at least one level, e.g. ["0.10"]'
        assert self.K_STRAIN in ('plate', 'rg', 'bb', 'nominal'), "K_STRAIN must be 'plate', 'rg', 'bb' or 'nominal'"
        assert self.PLATE_AREA in ('gap', 'bb'), "PLATE_AREA must be 'gap' or 'bb'"
        self.AXES = tuple(a for a in self.AXES if a in AXI)
        assert self.AXES, "AXES must name at least one of 'x', 'y', 'z'"
        self.mode, self.two_pist = 'compression', True
        self.COMP_LEVELS = [str(l) for l in self.COMP_LEVELS]
        if self.RUNS_ROOT is None:
            self.RUNS_ROOT = '/home/dpollard/Documents/lammps_runs/compress_slab'
        base = Path(self.base_dir)
        self.DATA_DIR = base / 'bulk' / self.RUN_ID
        self.PLOT_DIR = base / 'plots' / 'bulk' / self.RUN_ID
        self.TRAJ_DIR = base / 'traj_files.nosync'
        for d in (self.DATA_DIR, self.PLOT_DIR):
            d.mkdir(parents=True, exist_ok=True)
        self._axis = 'z'
        self._tags = {}
        if isinstance(self.NSTEPS, dict):
            self._tags = {str(k): str(v) for k, v in self.NSTEPS.items()}
        elif self.NSTEPS is not None:
            self.NSTEPS = int(self.NSTEPS)
        stem = f'{self.DATANAME}_{self.INTERACTION}'
        rt = self.tag_for(None)
        self.sim_name = f'{stem}_{rt}' if rt is not None else stem
        found = {('ref' if l is None else l): self.tag_for(l) for l in [None] + list(self.COMP_LEVELS)}
        print(f'Config: {self.sim_name}  |  eps_vol levels {self.COMP_LEVELS}  |  profile axes {list(self.AXES)}\n'
              f'  data  {self.DATA_DIR}\n  plots {self.PLOT_DIR}\n'
              f'  file tags (hold length): ' + ', '.join(f'{k}: {v if v else "not local yet"}' for k, v in found.items()))

    # ---- per-axis profile file names ----------------------------------------------------
    @staticmethod
    def pname(kind, axis, ref=False, comp=None, species=None):
        """base name of a profile file: kind = 'sigma' (needs comp, species), 'density', 'disp'.
        Along z the files carry the two-piston names; along x / y they are tagged _along<axis>."""
        r = '_ref' if ref else ''
        if kind == 'sigma':
            return f'sigma{comp}_{species}' + ('' if axis == 'z' else f'_along{axis}') + r
        if kind == 'density':
            return ('solvent_density_z' if axis == 'z' else f'solvent_density_along{axis}') + r
        if kind == 'disp':
            return 'disp_z_polymer' if axis == 'z' else f'disp_along{axis}_polymer'
        raise ValueError(kind)

    def plot(self, stem, lvl=None):
        t = self.tag_for(lvl) or 'untagged'
        ax = '' if self._axis == 'z' else f'_along{self._axis}'
        return self.PLOT_DIR / f'{stem}{ax}_{self.DATANAME}_{self.INTERACTION}_{t}{self.csuf(lvl)}.png'


# ===========================================================================
#  2. FILE READERS
# ===========================================================================
def read_table(path):
    """'# step a b c' fix print / fix ave/time table -> dict(name -> column) with 'step'
    (the first column, whatever its header says); {} when the file is missing or empty."""
    names, tab = tri.read_piston_table(path)
    if tab.size == 0:
        return {}
    if names is None:
        names = ['step'] + [f'col_{i}' for i in range(1, tab.shape[1])]
    out = {nm: tab[:, j] for j, nm in enumerate(names[:tab.shape[1]])}
    out['step'] = tab[:, 0]
    return out


def _stacks(cfg, axis, lvl=None, ref=False):
    """{comp: dict(ts, bins, p, s, t)} for one profile axis; a component whose files are
    missing is left out."""
    out = {}
    for comp in COMPONENTS:
        fp = cfg.path(cfg.pname('sigma', axis, ref, comp, 'polymer'), None if ref else lvl)
        fs = cfg.path(cfg.pname('sigma', axis, ref, comp, 'solvent'), None if ref else lvl)
        if not (fp.exists() and fs.exists()):
            continue
        sp, ss = tri.read_ave_time_file(fp), tri.read_ave_time_file(fs)
        if not sp or not ss:
            continue
        n = min(len(sp), len(ss))
        P = np.array([x[2] for x in sp[:n]])
        S = np.array([x[2] for x in ss[:n]])
        out[comp] = dict(ts=np.array([x[0] for x in sp[:n]], float), bins=sp[0][1], p=P, s=S, t=P + S)
    return out


def _density(cfg, axis, lvl=None, ref=False):
    f = cfg.path(cfg.pname('density', axis, ref), None if ref else lvl)
    if not f.exists():
        return None
    snaps = tri.read_ave_chunk_file(f)
    if not snaps:
        return None
    return dict(ts=np.array([s[0] for s in snaps], float), z=snaps[0][1][:, 1],
                n=np.array([s[1][:, 3] for s in snaps]))


def _track(T, key):
    """(step, value) of one column of a plate/piston table, or None."""
    return (T['step'], T[key]) if (T and key in T) else None


def _at(track, const):
    if track is None:
        return lambda t: float(const)
    st, v = track
    return lambda t: float(np.interp(t, st, v))


def face_areas(cfg, gaps):
    """{axis: area of the gel face normal to it} from the three plate gaps (PLATE_AREA='gap')."""
    e = {a: np.asarray(gaps[a], float) - 2.0 * cfg.contact_gap for a in 'xyz'}
    return {'x': e['y'] * e['z'], 'y': e['x'] * e['z'], 'z': e['x'] * e['y']}


def plate_pressures(cfg, F, P):
    """Plate loads -> pressures.  F: table with N_<face> (+ A_x A_y A_z, the live BB faces);
    P: plate-position table (for the gap areas).  Returns dict(step, <face>: P(t), 'x'|'y'|'z':
    pair means, 'mean': all six, area: {axis: A(t)})."""
    st = F['step']
    if cfg.PLATE_AREA == 'gap' and P:
        gaps = {a: np.interp(st, P['step'], P[a + 'hi'] - P[a + 'lo']) for a in 'xyz'}
        A = face_areas(cfg, gaps)
    else:
        A = {a: F[f'A_{a}'] for a in 'xyz'}
    out = dict(step=st, area=A, area_mode=('gap' if (cfg.PLATE_AREA == 'gap' and P) else 'bb'))
    for f in PLATES:
        out[f] = F[f'N_{f}'] / A[f[0]]
    for a in 'xyz':
        out[a] = 0.5 * (out[a + 'lo'] + out[a + 'hi'])
    out['mean'] = (out['x'] + out['y'] + out['z']) / 3.0
    return out


# ===========================================================================
#  3. REFERENCE STATE  (eps_vol = 0; shared by every level)
# ===========================================================================
def bath_mask(cfg, z, plate_lo, plate_hi, lo_lim, hi_lim, relax=True):
    """Bins of a core column that hold BATH only: entirely outside the plate planes by
    res_gel_gap_bins bins and inside [lo_lim, hi_lim] (the box faces along x and y; the wet
    piston planes -/+ res_wall_margin along z).  With relax=True the gap is reduced, then
    dropped, when the window would otherwise be empty (a thin bath)."""
    z = np.asarray(z, float)
    h = 0.5 * cfg.binWidth
    gaps = [cfg.res_gel_gap_bins * cfg.binWidth]
    if relax:
        gaps += [cfg.binWidth, 0.0]
    for g in gaps:
        m = ((z - h >= plate_hi + g) | (z + h <= plate_lo - g)) & (z - h >= lo_lim - 1e-9) & (z + h <= hi_lim + 1e-9)
        if m.sum() >= 1:
            return m
    return np.zeros(len(z), bool)


def _limits(cfg, axis, G, z_perm, z_feed):
    """[lo, hi] a bath bin must lie inside along `axis`."""
    if axis == 'z' and np.isfinite(z_perm) and np.isfinite(z_feed):
        return z_perm + cfg.res_wall_margin, z_feed - cfg.res_wall_margin
    return G['lo'], G['hi']


def _axis_geometry(cfg, axis, nbins, dens, box):
    """bin centres and box span along one profile axis.  The bin origin is the box face
    (compute chunk/atom bin/1d ... lower); it is read from the density file's Coord1."""
    bw = cfg.binWidth
    lo = float(dens['z'][0] - 0.5 * bw) if dens is not None else 0.0
    L = float(box[AXI[axis]]) if box is not None else nbins * bw
    z = lo + (np.arange(nbins) + 0.5) * bw
    full = z + 0.5 * bw <= lo + L + 1e-6            # the last bin is cut by the box face when L is no multiple of bw
    return dict(lo=lo, hi=lo + L, L=L, z=z, full=full)


def _box(cfg):
    for lvl in cfg.COMP_LEVELS:
        B = tri.load2c(cfg.path('box_dimensions', lvl), 4)
        if B is not None:
            return B[-1, 1:4]
    return None


def load_reference(cfg, verbose=True):
    """Everything that belongs to the uncompressed (eps_vol = 0) state: plate and wet-piston
    planes, reference volumes, the plate preload, and per profile axis the reference stress
    profiles with their Terzaghi split, gel bounds and solvent density.  Returns a dict R;
    R['ax'][axis] is the per-axis part."""
    say = print if verbose else (lambda *a, **k: None)
    R = dict(ax={}, box=_box(cfg))
    if R['box'] is None:
        raise FileNotFoundError('no box_dimensions_<stem>_c<level> file found -- run the sync cell')
    R['AREA'] = float(R['box'][0] * R['box'][1])
    say(f'box: {R["box"][0]:.2f} x {R["box"][1]:.2f} x {R["box"][2]:.2f}   (x, y periodic bath; wet pistons close z)')

    # ---- plate planes, gaps, reference volumes -----------------------------------
    P = read_table(cfg.path('plate_position_ref'))
    if not P:
        raise FileNotFoundError(f'{cfg.path("plate_position_ref").name} missing -- run the sync cell')
    R['plate_tab'] = P
    R['plane'] = {f: float(np.mean(P[f])) for f in PLATES}
    R['gap'] = {a: R['plane'][a + 'hi'] - R['plane'][a + 'lo'] for a in 'xyz'}
    R['A_gap'] = {a: float(v) for a, v in face_areas(cfg, R['gap']).items()}
    R['L0'] = {a: R['gap'][a] - 2.0 * cfg.contact_gap for a in 'xyz'}      # reference edges = the deck's L0_x, L0_y, L0_z
    R['V0_gap'] = float(np.prod([R['L0'][a] for a in 'xyz']))
    say('seated plates: ' + '   '.join(f'{a}: {R["plane"][a + "lo"]:.2f} | {R["plane"][a + "hi"]:.2f} (gap {R["gap"][a]:.2f})' for a in 'xyz'))
    S = read_table(cfg.path('strain_vol_ref'))
    R['strain_tab'] = S
    R['V0_bb'] = float(np.mean(S['V_bb'])) if S else np.nan
    R['V0_rg'] = float(np.mean(S['V_rg'])) if S else np.nan
    B = tri.load2c(cfg.path('gel_dimensions_bb_ref'), 4)
    R['L0_bb'] = {a: (float(np.mean(B[:, 1 + AXI[a]])) if B is not None else np.nan) for a in 'xyz'}
    say('reference edges L0 (plate gap - 2 contact_gap) = ' + ' x '.join(f'{R["L0"][a]:.2f}' for a in 'xyz')
        + f'   V0 = {R["V0_gap"]:.1f}')
    say('reference network: L_bb = ' + ' x '.join(f'{R["L0_bb"][a]:.2f}' for a in 'xyz')
        + f'   V_bb = {R["V0_bb"]:.1f}   V_rg = {R["V0_rg"]:.1f}   (V_rg / V0 = {R["V0_rg"] / R["V0_gap"]:.3f})')
    W = read_table(cfg.path('piston_position_ref'))
    R['z_feed'] = float(np.mean(W['z_feed'])) if 'z_feed' in W else np.nan
    R['z_perm'] = float(np.mean(W['z_perm'])) if 'z_perm' in W else np.nan
    say(f'wet pistons: permeate z = {R["z_perm"]:.2f}   feed z = {R["z_feed"]:.2f}')

    # ---- plate preload ---------------------------------------------------------
    F = read_table(cfg.path('plate_force_avg_ref'))
    if F:
        PP = plate_pressures(cfg, F, P)
        R['plates'] = PP
        m, lo, hi, blk, _ = block_bootstrap_ci(PP['mean'], cfg.ci_level)
        R.update(P_ref=float(m), P_ref_lo=float(lo), P_ref_hi=float(hi))
        R['P_ref_face'] = {f: float(np.mean(PP[f])) for f in PLATES}
        R['P_ref_ax'] = {a: float(np.mean(PP[a])) for a in 'xyz'}
        say(f'reference plate preload <P_plates> = {m:.4f} [{lo:.4f}, {hi:.4f}] (n={len(PP["mean"])}, block={blk}, '
            f'area = {PP["area_mode"]});  per axis ' + '  '.join(f'{a}: {R["P_ref_ax"][a]:+.4f}' for a in 'xyz'))
    else:
        R['P_ref'] = np.nan
        say('reference plate preload: no plate_force_avg_ref file -> the plate K stays absolute')
    Wf = read_table(cfg.path('piston_force_avg_ref'))
    for lab, sgn in (('feed', 1.0), ('perm', -1.0)):
        if f'F_fluid_{lab}' in Wf:
            m, lo, hi, *_ = block_bootstrap_ci(sgn * Wf[f'F_fluid_{lab}'] / R['AREA'], cfg.ci_level)
            R[f'P_ref_{lab}'], R[f'P_ref_{lab}_lo'], R[f'P_ref_{lab}_hi'] = float(m), float(lo), float(hi)
    if 'P_ref_feed' in R:
        say('reference wet-piston (bath) pressures: ' + '  '.join(
            f'{lab}: {R[f"P_ref_{lab}"]:.4f} [{R[f"P_ref_{lab}_lo"]:.4f}, {R[f"P_ref_{lab}_hi"]:.4f}]'
            for lab in ('feed', 'perm') if f'P_ref_{lab}' in R))

    # ---- per-axis reference profiles ------------------------------------------------
    for a in cfg.AXES:
        st = _stacks(cfg, a, ref=True)
        if not all(c in st for c in COMPONENTS):
            if a == cfg.AXES[0]:
                raise FileNotFoundError(f'reference profile files along {a} missing -- run the sync cell')
            say(f'NOTE: reference profiles along {a} incomplete -> axis skipped')
            continue
        dens = _density(cfg, a, ref=True)
        nb = st['zz']['t'].shape[1]
        G = _axis_geometry(cfg, a, nb, dens, R['box'])
        z = G['z']
        A = dict(axis=a, stress=st, z=z, full=G['full'], lo=G['lo'], hi=G['hi'], L=G['L'],
                 plate_lo=R['plane'][a + 'lo'], plate_hi=R['plane'][a + 'hi'], n_ref=len(st['zz']['ts']))
        ptr = np.abs(np.nanmean((st['xx']['p'] + st['yy']['p'] + st['zz']['p']) / 3.0, axis=0))
        gel = (ptr > cfg.gel_thresh * float(ptr.max())) if ptr.max() > 0 else np.zeros(nb, bool)
        A['z_gel_lo'] = float(z[gel].min()) if gel.any() else A['plate_lo']
        A['z_gel_hi'] = float(z[gel].max()) if gel.any() else A['plate_hi']
        A['in_gel'] = (z >= A['z_gel_lo']) & (z <= A['z_gel_hi'])
        A['interior'] = A['in_gel'] & (z >= A['plate_lo'] + cfg.wall_margin) & (z <= A['plate_hi'] - cfg.wall_margin)
        lim = _limits(cfg, a, G, R['z_perm'], R['z_feed'])
        A['bw'] = bath_mask(cfg, z, A['plate_lo'], A['plate_hi'], *lim) & G['full']
        if not A['bw'].any():
            A['bw'] = bath_mask(cfg, z, A['plate_lo'], A['plate_hi'], *lim)
        for comp, Sx in st.items():
            Sx['t_m'], Sx['t_lo'], Sx['t_hi'] = mean_ci(Sx['t'], cfg.ci_level)
            Sx['p_m'], Sx['s_m'] = np.nanmean(Sx['p'], axis=0), np.nanmean(Sx['s'], axis=0)
            Sx['sd_bin'] = np.nanstd(Sx['t'], axis=0)
            net, pore, _, _ = tri.terzaghi_split(Sx['t'], A['bw'], Sx['sd_bin'], cfg.ci_level)
            Sx['net'] = net
            Sx['net_m'], Sx['net_lo'], Sx['net_hi'] = mean_ci(net, cfg.ci_level)
            Sx['pore'] = float(np.mean(pore))
            im = A['interior']
            Sx['net_interior'] = float(np.nanmean(Sx['net_m'][im])) if im.any() else np.nan
            _, ilo, ihi = mean_ci(Sx['net_m'][im], cfg.ci_level) if im.sum() > 1 else (np.nan, np.nan, np.nan)
            Sx['net_interior_half'] = float(0.5 * (ihi - ilo)) if np.isfinite(ihi) else 0.0
        Pn = tri._tr3(st, 'net')
        Pm = np.nanmean(Pn, axis=0)
        im = A['interior']
        A['Pi_ref'] = float(np.nanmean(Pm[im])) if im.any() else np.nan
        _, ilo, ihi = mean_ci(Pm[im], cfg.ci_level) if im.sum() > 1 else (np.nan, np.nan, np.nan)
        A['Pi_ref_half'] = float(0.5 * (ihi - ilo)) if np.isfinite(ihi) else 0.0
        if dens is not None:
            rho_m = np.nanmean(dens['n'], axis=0)
            A['rho_s0'] = float(np.nanmean(rho_m[A['bw']])) if A['bw'].any() else float(np.nanmax(rho_m))
            A['phi_mf'] = mean_ci(dens['n'] / A['rho_s0'], cfg.ci_level)
        else:
            A['rho_s0'], A['phi_mf'] = np.nan, None
        R['ax'][a] = A
        say(f'along {a}: {A["n_ref"]} snapshots on {nb} bins [{z.min():.1f}, {z.max():.1f}];  gel {A["z_gel_lo"]:.1f}..{A["z_gel_hi"]:.1f} '
            f'({int(A["in_gel"].sum())} bins, interior {int(im.sum())});  bath window {tri.mask_span(z, A["bw"])}')
        say(f'          p_pore  ' + '  '.join(f'{c}: {st[c]["pore"]:.4f}' for c in COMPONENTS)
            + f"   |   sigma'(interior)  " + '  '.join(f"{c}: {st[c]['net_interior']:+.4f}" for c in COMPONENTS)
            + f'   |   Pi_ref = {A["Pi_ref"]:+.4f} ± {A["Pi_ref_half"]:.4f}   rho_bath = {A["rho_s0"]:.4f}')
    fin = [R['ax'][a]['Pi_ref'] for a in R['ax'] if np.isfinite(R['ax'][a]['Pi_ref'])]
    R['Pi_ref'] = float(np.mean(fin)) if fin else np.nan
    R['Pi_ref_half'] = float(np.mean([R['ax'][a]['Pi_ref_half'] for a in R['ax']])) if fin else 0.0
    if np.isfinite(R.get('P_ref', np.nan)) and fin:
        say(f'zero-strain readings: profiles Pi_ref = {R["Pi_ref"]:+.4f}, plates P_ref = {R["P_ref"]:+.4f}  '
            f'(difference {R["Pi_ref"] - R["P_ref"]:+.4f}: each estimator subtracts its own)')
    return R


# ===========================================================================
#  4. ONE STRAIN LEVEL
# ===========================================================================
def _load_disp(cfg, axis, lvl, Pt, R):
    """displacement profile along `axis` + the held plate planes, in triaxial.fit_Dc's format."""
    f = cfg.path(cfg.pname('disp', axis), lvl)
    if not f.exists() or not Pt:
        return None
    snaps = tri.read_ave_chunk_file(f)
    if not snaps:
        return None
    hi, lo = Pt[axis + 'hi'], Pt[axis + 'lo']
    d = dict(ts=np.array([s[0] for s in snaps], float), z=snaps[0][1][:, 1],
             Nc=np.array([s[1][:, 2] for s in snaps]), uz=np.array([s[1][:, 3] for s in snaps]),
             z_pist_held=float(hi[-1]), z_supp_held=float(lo[-1]),
             t_hold=float(Pt['step'][np.argmax(np.isclose(hi, hi[-1]))]))
    B = tri.load2c(cfg.path('gel_dimensions_bb', lvl), 4)
    if B is not None:
        h = B[:, 0] >= d['t_hold'] + 0.1 * (B[-1, 0] - d['t_hold'])
        rows = B[h] if h.any() else B[-1:]
        d['L_bb'] = float(np.median(rows[:, 1 + AXI[axis]]))
    return d


def load_level(cfg, R, lvl, verbose=True):
    """Everything for ONE volumetric-strain level `lvl` (string, e.g. "0.10"): measured strains,
    per-axis stress stacks with their Terzaghi split, the plate histories and plateau, K
    (network per axis + mean, plates), the unrelaxed-hold systematic, the isotropy check,
    the D_c fit per axis, the wet-piston bath check.  Returns a dict L, or None when the
    core files are missing."""
    say = print if verbose else (lambda *a, **k: None)
    a0 = cfg.AXES[0]
    st0 = _stacks(cfg, a0, lvl=lvl) if a0 in R['ax'] else {}
    if 'zz' not in st0:
        say(f'  level {lvl}: profile files along {a0} missing -- skipping level')
        return None
    L = dict(lvl=lvl, eps_nominal=float(lvl), ax={})
    ts = st0['zz']['ts']
    t0, t1 = float(ts[0]), float(ts[-1])
    L['ts'] = ts
    L['halt_ts'] = int(t1 - cfg.plateau_frac * (t1 - t0))
    L['evol_from'] = int(t0)
    say(f'\n=== level _c{lvl}  (target eps_vol {float(lvl):.3f};  files tagged _{cfg.tag_for(lvl)} = hold length) ===')
    say(f'  production: {len(ts)} snapshots, steps {int(t0)} -> {int(t1)};  plateau window: steps >= {L["halt_ts"]} '
        f'(last {cfg.plateau_frac:.0%})')

    # ---- plates, strains ---------------------------------------------------------
    Pt = read_table(cfg.path('plate_position', lvl))
    L['plate_tab'] = Pt
    L['plane'] = {f: (float(Pt[f][-1]) if Pt else R['plane'][f]) for f in PLATES}
    L['gap'] = {a: L['plane'][a + 'hi'] - L['plane'][a + 'lo'] for a in 'xyz'}
    L['closure'] = {a: R['gap'][a] - L['gap'][a] for a in 'xyz'}
    L['eps_lin'] = {a: L['closure'][a] / R['L0'][a] for a in 'xyz'}
    eps_plate_geo = 1.0 - float(np.prod([1.0 - L['eps_lin'][a] for a in 'xyz']))
    S = read_table(cfg.path('strain_vol', lvl))
    L['strain_tab'] = S
    L['strain'] = dict(nominal=float(lvl), plate=eps_plate_geo, plate_geo=eps_plate_geo)
    if S:
        pl = S['step'] >= L['halt_ts']
        pl = pl if pl.any() else np.ones(len(S['step']), bool)
        L['strain']['plate'] = float(np.mean(S['eps_vol_plate'][pl]))     # the deck's own reading (its L0)
        for key, V0 in (('bb', R['V0_bb']), ('rg', R['V0_rg'])):
            if np.isfinite(V0) and f'V_{key}' in S:
                v = 1.0 - S[f'V_{key}'][pl] / V0
                L['strain'][key] = float(np.mean(v))
                L['strain'][key + '_sd'] = float(np.std(v))
    L['eps'] = L['strain'].get(cfg.K_STRAIN, np.nan)
    if not np.isfinite(L['eps']) or L['eps'] <= 0:
        say(f"  WARNING: K_STRAIN='{cfg.K_STRAIN}' strain is {L['eps']} -- falling back to the plate strain")
        L['eps'] = L['strain']['plate']
    eps = L['eps']
    say(f"  plates (held): " + '   '.join(f'{a}: gap {L["gap"][a]:.3f} (closed {L["closure"][a]:.3f} = {L["eps_lin"][a]:.4f} L0)' for a in 'xyz'))
    say(f"  eps_vol: target {float(lvl):.4f} | plates {L['strain']['plate']:.5f} (K denominator: {cfg.K_STRAIN} = {eps:.5f}) | "
        f"network Rg {L['strain'].get('rg', np.nan):.4f} ± {L['strain'].get('rg_sd', np.nan):.4f}   "
        f"BB {L['strain'].get('bb', np.nan):.4f} ± {L['strain'].get('bb_sd', np.nan):.4f}")
    if abs(L['strain']['plate'] - float(lvl)) > 0.02 * float(lvl):
        say(f"  NOTE: the plate strain differs from the target by {L['strain']['plate'] / float(lvl) - 1:+.1%} "
            '(drive stopped early / overshot): the measured value is used')

    # ---- wet pistons -------------------------------------------------------------
    Wp = read_table(cfg.path('piston_position', lvl))
    feed_tr, perm_tr = _track(Wp, 'z_feed'), _track(Wp, 'z_perm')
    feed_at, perm_at = _at(feed_tr, R['z_feed']), _at(perm_tr, R['z_perm'])
    L['z_feed'], L['z_perm'] = feed_at(t1), perm_at(t1)
    L['wetz'] = dict(feed_at=feed_at, perm_at=perm_at, feed_track=feed_tr, perm_track=perm_tr)

    # ---- per-axis stress, Terzaghi, Pi -------------------------------------------
    L['plat'] = ts >= L['halt_ts']
    if L['plat'].sum() < 1:
        L['plat'][-1] = True
    for a in cfg.AXES:
        if a not in R['ax']:
            continue
        st = st0 if a == a0 else _stacks(cfg, a, lvl=lvl)
        if not all(c in st for c in COMPONENTS):
            say(f'  NOTE: profiles along {a} incomplete for level {lvl} -> axis skipped')
            continue
        Ra = R['ax'][a]
        z = Ra['z']
        n = min(len(ts), len(st['zz']['ts']))
        A = dict(axis=a, stress=st, z=z, plate_lo=L['plane'][a + 'lo'], plate_hi=L['plane'][a + 'hi'],
                 lo_at=_at(_track(Pt, a + 'lo'), R['plane'][a + 'lo']), hi_at=_at(_track(Pt, a + 'hi'), R['plane'][a + 'hi']))
        ptr = np.abs((st['xx']['p'][-1] + st['yy']['p'][-1] + st['zz']['p'][-1]) / 3.0)
        mem = ptr > cfg.gel_thresh * float(np.nanmax(ptr))
        A['z_mem_lo'] = float(z[mem].min()) if mem.any() else Ra['z_gel_lo']
        A['z_mem_hi'] = float(z[mem].max()) if mem.any() else Ra['z_gel_hi']
        A['in_mem'] = (z >= A['z_mem_lo']) & (z <= A['z_mem_hi'])
        A['interior'] = A['in_mem'] & (z >= A['plate_lo'] + cfg.wall_margin) & (z <= A['plate_hi'] - cfg.wall_margin)
        G = dict(lo=Ra['lo'], hi=Ra['hi'])
        rows = []
        for t in ts[:n]:
            lim = _limits(cfg, a, G, perm_at(t), feed_at(t))
            m = bath_mask(cfg, z, A['lo_at'](t), A['hi_at'](t), *lim) & Ra['full']
            rows.append(m if m.any() else Ra['bw'])
        A['bw'] = np.array(rows, bool)
        for comp, Sx in st.items():
            Rs = Ra['stress'][comp]
            sd_bin = Rs['sd_bin']
            if Rs['t'].shape[0] < 2 or not np.any(sd_bin > 0):
                sd_bin = np.full(len(z), float(np.nanstd(Sx['t'][-1][A['bw'][-1]])))
            for k in ('p', 's', 't', 'ts'):
                Sx[k] = Sx[k][:n]
            Sx['net'], Sx['pore'], Sx['pore_half'], Sx['net_half'] = tri.terzaghi_split(Sx['t'], A['bw'], sd_bin, cfg.ci_level)
            Sx['ref_net'] = Rs['net_interior'] if cfg.K_SUBTRACT_REF else 0.0
            Sx['net_plat'] = np.nanmean(Sx['net'][L['plat'][:n]], axis=0)
            Sx['t_plat'] = np.nanmean(Sx['t'][L['plat'][:n]], axis=0)
        im = A['interior'] if A['interior'].sum() >= 3 else A['in_mem']
        A['K_mask'] = 'interior' if im is A['interior'] else 'gel'
        Pn = tri._tr3(st, 'net')
        A['Pi'] = Pn                                        # (n_snap, n_bins)
        A['Pi_plat'] = np.nanmean(Pn[L['plat'][:n]], axis=0)
        A['Pi_series'] = np.array([np.nanmean(Pn[i][im]) for i in range(n)])
        pib = A['Pi_plat'][im]
        pib = pib[np.isfinite(pib)]
        A['Pi_int'], plo, phi = mean_ci(pib, cfg.ci_level)
        A['Pi_int'], A['Pi_half'] = float(A['Pi_int']), float(0.5 * (phi - plo))
        ref = Ra['Pi_ref'] if (cfg.K_SUBTRACT_REF and np.isfinite(Ra['Pi_ref'])) else 0.0
        ref_h = Ra['Pi_ref_half'] if cfg.K_SUBTRACT_REF else 0.0
        A['K_ref'] = ref
        A['K_abs'], A['K_abs_lo'], A['K_abs_hi'] = (A['Pi_int'] / eps, (A['Pi_int'] - A['Pi_half']) / eps,
                                                    (A['Pi_int'] + A['Pi_half']) / eps)
        half = np.sqrt(A['Pi_half'] ** 2 + ref_h ** 2) / eps
        A['K'] = (A['Pi_int'] - ref) / eps
        A['K_lo'], A['K_hi'], A['K_nbins'] = A['K'] - half, A['K'] + half, int(len(pib))
        # isotropy: each normal component's increment against the mean (the deviator)
        A['dsig'] = {c: float(np.nanmean(st[c]['net_plat'][im]) - st[c]['ref_net']) for c in COMPONENTS}
        A['dev'] = {c: A['dsig'][c] - float(np.mean(list(A['dsig'].values()))) for c in COMPONENTS}
        dens = _density(cfg, a, lvl=lvl)
        if dens is not None and np.isfinite(Ra['rho_s0']):
            pl = dens['ts'] >= L['halt_ts']
            mf = dens['n'] / Ra['rho_s0']
            A['phi_mf'] = mean_ci(mf[pl] if pl.any() else mf[-1:], cfg.ci_level)
        else:
            A['phi_mf'] = None
        L['ax'][a] = A
        say(f"  along {a}: p_pore(zz) {st['zz']['pore'][0]:.4f} -> {st['zz']['pore'][-1]:.4f};  Pi in the {A['K_mask']} "
            f"({A['K_nbins']} bins, plateau) = {A['Pi_int']:.4f} ± {A['Pi_half']:.4f}  (ref {Ra['Pi_ref']:+.4f})   "
            f"K_network = {A['K']:.4f} [{A['K_lo']:.4f}, {A['K_hi']:.4f}]")

    # ---- plate loads + plateau ------------------------------------------------
    F = read_table(cfg.path('plate_force_avg', lvl))
    Fr = read_table(cfg.path('plate_force', lvl))
    if F:
        L['plates'] = plate_pressures(cfg, F, Pt)
    if Fr:
        L['plates_raw'] = plate_pressures(cfg, Fr, Pt)
    if 'plates' in L or 'plates_raw' in L:
        PP = L.get('plates', L.get('plates_raw'))
        L['PF'] = plateau_window(PP['step'], PP['mean'], cfg.plateau_frac_auto, cfg.ci_level)
        L['PF']['src'] = 'LMP block-avg' if 'plates' in L else 'raw'
        sel = PP['step'] >= L['PF']['step0']
        L['P_face'] = {f: float(np.mean(PP[f][sel])) for f in PLATES}
        L['P_ax'] = {a: float(np.mean(PP[a][sel])) for a in 'xyz'}
        p = L['PF']
        say(f"  plate plateau <P_plates> = {p['mean']:.4f} [{p['lo']:.4f}, {p['hi']:.4f}]  ({p['src']}; area = {PP['area_mode']}; "
            f"auto window last {p['frac']:.0%}, n={p['n']}, block={p['block']}, tau~{p['tau']:.1f})"
            + ('  DRIFT WARNING: no drift-free window -- extend the hold' if p.get('warn') else ''))
        say('  per axis: ' + '  '.join(f'{a}: {L["P_ax"][a]:.4f}' for a in 'xyz')
            + '   per plate: ' + '  '.join(f'{f}: {L["P_face"][f]:.4f}' for f in PLATES))
    else:
        say('  NOTE: no plate_force file -> K_plates skipped')

    # ---- K: network (mean over the profile axes) and plates --------------------------
    sub = bool(cfg.K_SUBTRACT_REF)
    axs = [a for a in cfg.AXES if a in L['ax']]
    L['K_net'] = float(np.mean([L['ax'][a]['K'] for a in axs]))
    h = float(np.mean([0.5 * (L['ax'][a]['K_hi'] - L['ax'][a]['K_lo']) for a in axs]))   # overlapping columns: not independent
    L['K_net_lo'], L['K_net_hi'] = L['K_net'] - h, L['K_net'] + h
    L['K_net_abs'] = float(np.mean([L['ax'][a]['K_abs'] for a in axs]))
    ha = float(np.mean([0.5 * (L['ax'][a]['K_abs_hi'] - L['ax'][a]['K_abs_lo']) for a in axs]))
    L['K_net_abs_lo'], L['K_net_abs_hi'] = L['K_net_abs'] - ha, L['K_net_abs'] + ha
    L['K_net_ref'] = float(np.mean([L['ax'][a]['K_ref'] for a in axs]))
    L['Pi_int'] = float(np.mean([L['ax'][a]['Pi_int'] for a in axs]))
    L['Pi_half'] = float(np.mean([L['ax'][a]['Pi_half'] for a in axs]))
    L['K_axes_spread'] = float(np.std([L['ax'][a]['K'] for a in axs])) if len(axs) > 1 else 0.0
    # the same stress increment over the network's OWN volumetric strain (Rg)
    e_rg = L['strain'].get('rg', np.nan)
    L['K_net_rg'] = float(L['K_net'] * eps / e_rg) if (np.isfinite(e_rg) and e_rg > 0) else np.nan
    if 'PF' in L:
        p = L['PF']
        L['K_pl_abs'], L['K_pl_abs_lo'], L['K_pl_abs_hi'] = p['mean'] / eps, p['lo'] / eps, p['hi'] / eps
        L['P_final'] = p['mean']
        have = sub and np.isfinite(R.get('P_ref', np.nan))
        L['K_pl_ref'] = 'measured' if have else ('absent' if sub else 'off')
        L['P_ref'] = float(R['P_ref']) if have else 0.0
        if have:
            half = np.sqrt(((p['hi'] - p['lo']) / 2) ** 2 + ((R['P_ref_hi'] - R['P_ref_lo']) / 2) ** 2) / eps
            L['K_pl'] = (p['mean'] - R['P_ref']) / eps
            L['K_pl_lo'], L['K_pl_hi'] = L['K_pl'] - half, L['K_pl'] + half
        else:
            L['K_pl'], L['K_pl_lo'], L['K_pl_hi'] = L['K_pl_abs'], L['K_pl_abs_lo'], L['K_pl_abs_hi']
    # ---- unrelaxed-hold systematic (RELAX_SYS): the plate tail's excess over its fitted
    #      asymptote, / eps_vol, widens the LOWER bound of both K estimators
    L['K_sys'] = 0.0
    if cfg.RELAX_SYS and 'PF' in L and 'plates' in L:
        rt = relax_tail(L['plates']['step'], L['plates']['mean'], L['PF']['mean'], cfg.dt_lj, cfg.RELAX_TAIL_FRAC)
        L['relax'] = rt
        L['K_sys'] = rt['excess'] / eps
        L['K_net_lo'] -= L['K_sys']
        for a in axs:
            L['ax'][a]['K_lo'] -= L['K_sys']
        if 'K_pl' in L:
            L['K_pl_lo'] -= L['K_sys']
        if rt['ok']:
            say(f"  unrelaxed-hold systematic: plate tail -> P_inf = {rt['x_inf']:.4f} (tau ~ {rt['tau']:.0f} tau over the last "
                f"{cfg.RELAX_TAIL_FRAC:.0%}), plateau excess {rt['excess']:+.4f}  ->  delta_sys = {L['K_sys']:.4f} "
                'added to the LOWER bound of K_network and K_plates')
        else:
            say(f"  unrelaxed-hold systematic: 0 ({rt['why']})")
    how = 'increment from eps_vol = 0' if sub else 'absolute'
    say(f"  K_network = {L['K_net']:.4f} [{L['K_net_lo']:.4f}, {L['K_net_hi']:.4f}]  ({how}; mean of the profiles along "
        f"{', '.join(axs)}; spread between them {L['K_axes_spread']:.4f}"
        + (f"; ref Pi subtracted {L['K_net_ref']:+.4f}; absolute {L['K_net_abs']:.4f}" if sub else '') + ')')
    if np.isfinite(L['K_net_rg']) and cfg.K_STRAIN != 'rg':
        say(f"  K_network over the network's own strain (Rg, {e_rg:.4f}) = {L['K_net_rg']:.4f}   "
            f"[network strain / {cfg.K_STRAIN} strain = {e_rg / eps:.3f}]")
    if 'K_pl' in L:
        note = {'measured': f"ref P subtracted {L['P_ref']:+.4f}; absolute {L['K_pl_abs']:.4f}",
                'absent': 'NO plate_force_avg_ref -> ABSOLUTE', 'off': 'absolute'}[L['K_pl_ref']]
        say(f"  K_plates  = {L['K_pl']:.4f} [{L['K_pl_lo']:.4f}, {L['K_pl_hi']:.4f}]  ({note})"
            f"   ratio K_plates/K_network = {L['K_pl'] / L['K_net']:.4f}")

    # ---- isotropy --------------------------------------------------------------
    dv = [abs(v) for a in axs for v in L['ax'][a]['dev'].values()]
    dP = L['Pi_int'] - (L['K_net_ref'] if sub else 0.0)
    L['dev_max'] = float(max(dv)) if dv else np.nan
    L['dev_rel'] = float(L['dev_max'] / abs(dP)) if dP else np.nan
    say(f"  isotropy: largest |sigma'_ii increment - mean| over the axes = {L['dev_max']:.4f} = {L['dev_rel']:.1%} of the Pi increment"
        + ('' if 'P_ax' not in L else '   plate pairs vs their mean: ' + '  '.join(
            f'{a}: {L["P_ax"][a] - R.get("P_ref_ax", {}).get(a, 0.0) - (L["PF"]["mean"] - L.get("P_ref", 0.0)):+.4f}' for a in 'xyz')))

    # ---- wet pistons: bath check + solvent expelled --------------------------------
    L['wet'] = tri.load_wet_pistons(cfg, R, lvl, plat_from=L['halt_ts'])
    if L['wet'] is not None:
        W = L['wet']
        say('  bath (wet pistons, plateau means): ' + '  '.join(
            f"{k}: {W['plat'][k]:.4f}" for k in ('P_feed_meas', 'P_perm_meas') if k in W['plat'])
            + (f"   applied {W['plat'].get('P_feed_app', np.nan):.3f}" if 'P_feed_app' in W['plat'] else ''))
        if W.get('dV_total') is not None:
            L['eps_expelled'] = float(W['dV_total'][-1] / R['V0_gap'])
            say(f"  solvent expelled dV_total = {W['dV_total'][-1]:.1f} sigma^3 = {L['eps_expelled']:.4f} of the volume the seated "
                f"plates enclosed (plate strain {L['strain']['plate']:.4f})")

    # ---- D_c: the held-cube fit (one D_c, three axes) + the 1-D per-axis fits as a check ----
    L['disp'] = {a: _load_disp(cfg, a, lvl, Pt, R) for a in axs}
    L['Dc_axis'] = {}
    for a in axs:
        Ra = R['ax'][a]
        F_ = tri.fit_Dc(cfg, dict(z_support=Ra['plate_lo'], z_piston=Ra['plate_hi']), L['disp'][a])
        if F_ is not None:
            L['Dc_axis'][a] = F_
    L['Dc'] = fit_Dc_cube(cfg, R, L)
    if L['Dc'] is not None:
        F = L['Dc']
        L['Dc_mean'] = F['Dc']
        say(f"  D_c (held-cube fit, {len(F['ax'])} axes jointly, {F['n_amp']} amplitudes) = {F['Dc']:.4e} sigma^2/tau  "
            f"(R^2 = {F['R2']:.3f}; L = " + ' x '.join(f"{F['ax'][a]['L']:.1f}" for a in F['ax']) + f"; hold = {F['hold_T']:.0f} tau)")
        if L['Dc_axis']:
            say('    1-D held-slab fits per axis (apparent): ' + '   '.join(
                f"{a}: {F_['Dc']:.4e} (R^2 {F_['R2']:.3f})" for a, F_ in L['Dc_axis'].items()))
    else:
        say('  NOTE: no D_c (missing disp files / plate_position, or too few polymer atoms per bin: see Ncount_min)')

    # ---- closure with the triaxial M ------------------------------------------------
    if cfg.M_REF is not None:
        M = float(cfg.M_REF)
        L['G_from_MK'] = 0.75 * (M - L['K_net'])
        L['G_from_MK_lo'], L['G_from_MK_hi'] = 0.75 * (M - L['K_net_hi']), 0.75 * (M - L['K_net_lo'])
        say(f"  closure M = K + 4G/3 with M_REF = {M:.4f}:  G = {L['G_from_MK']:.4f} [{L['G_from_MK_lo']:.4f}, {L['G_from_MK_hi']:.4f}]"
            + (f"   (G_REF = {float(cfg.G_REF):.4f})" if cfg.G_REF is not None else ''))
        if L['Dc'] is not None:
            L['kappa'] = L['Dc_mean'] / M
            say(f"  kappa = D_c / M_REF = {L['kappa']:.4e}")
    return L


# ===========================================================================
#  5. EXPANSE SYNC
# ===========================================================================
_PROD = ('strain_vol', 'strain_zz', 'strain_piston', 'plate_position', 'plate_force', 'plate_force_avg',
         'piston_position', 'piston_force', 'piston_force_avg', 'piston_pressure', 'permeation', 'pressure_reservoirs',
         'support_position', 'box_dimensions', 'gel_dimensions_bb', 'gel_dimensions_rg', 'polymer_com', 'gel_edges')
_REF = ('strain_vol_ref', 'plate_position_ref', 'plate_force_avg_ref', 'piston_position_ref', 'piston_force_avg_ref',
        'gel_dimensions_bb_ref', 'pressure_reservoirs_ref')
_REQUIRED = ('strain_vol', 'plate_position', 'plate_force_avg', 'box_dimensions')


def _profile_names(cfg, ref):
    out = []
    for a in cfg.AXES:
        out += [cfg.pname('sigma', a, ref, c, s) for c in COMPONENTS for s in ('polymer', 'solvent')]
        out.append(cfg.pname('density', a, ref))
        if not ref:
            out.append(cfg.pname('disp', a))
    return out


def sync_files(cfg, levels=None):
    """(data_files, traj_files, required) the notebooks read for `levels`."""
    levels = cfg.COMP_LEVELS if levels is None else [str(l) for l in levels]
    data = [cfg.path(n) for n in tuple(_REF) + tuple(_profile_names(cfg, True))]
    data.append(cfg.DATA_DIR / f'stress_aniso_{cfg.DATANAME}_{cfg.INTERACTION}_*.dat')   # whole-run file, tagged with NSTEPS
    req = [cfg.path(n) for n in ('plate_position_ref', 'strain_vol_ref', cfg.pname('sigma', cfg.AXES[0], True, 'zz', 'polymer'))]
    traj = [cfg.traj('traj_ref')] if cfg.SYNC_TRAJ else []
    for l in levels:
        data += [cfg.path(n, l) for n in tuple(_PROD) + tuple(_profile_names(cfg, False))]
        req += [cfg.path(n, l) for n in _REQUIRED + (cfg.pname('sigma', cfg.AXES[0], False, 'zz', 'polymer'),)]
        if cfg.SYNC_TRAJ:
            traj.append(cfg.traj('traj_stress', l))
    return data, traj, req


def sync_from_expanse(cfg, levels=None, force=False):
    """Pull every file the notebooks read from Expanse in ONE login (triaxial.sync_pull).
    Runs made on Pod or Bridges-2 are not reachable this way: copy their output_files/*/
    contents into cfg.DATA_DIR by hand and leave SYNC off."""
    if cfg.HOST != 'expanse':
        print(f"HOST = '{cfg.HOST}': no automatic sync -- copy the run's output_files/*/*.dat into {cfg.DATA_DIR}")
        return
    data, traj, req = sync_files(cfg, levels)
    optional = {cfg.path(n).name for n in ('piston_force_avg_ref', 'pressure_reservoirs_ref', 'gel_dimensions_bb_ref')}
    optional.add(f'stress_aniso_{cfg.DATANAME}_{cfg.INTERACTION}_*.dat')
    optional |= {cfg.path(n, l).name for n in ('pressure_reservoirs', 'support_position', 'gel_edges', 'polymer_com',
                                               'strain_zz', 'strain_piston', 'piston_force', 'gel_dimensions_rg')
                 for l in cfg.COMP_LEVELS}
    tri.sync_pull(cfg, data, traj, req, optional, force, refresh=lambda: sync_files(cfg, levels))


# ===========================================================================
#  6. PER-AXIS VIEWS  (triaxial's profile figures, drawn along any axis)
# ===========================================================================
def view(cfg, R, L=None, axis=None):
    """(Rv, Lv): the data along one profile axis in the layout triaxial's figure functions
    read (R['z'], R['stress'][comp], L['interior'], L['z_pist'], ...).  The lower plate of the
    axis plays the support, the upper one the loading piston; the wet pistons exist along z only."""
    a = axis or cfg.AXES[0]
    Ra = R['ax'][a]
    Rv = dict(z=Ra['z'], Z_LO=Ra['lo'], Z_HI=Ra['hi'], LZ=Ra['L'], AREA=R['AREA'], stress=Ra['stress'],
              z_gel_lo=Ra['z_gel_lo'], z_gel_hi=Ra['z_gel_hi'], in_gel=Ra['in_gel'], interior=Ra['interior'],
              z_support=Ra['plate_lo'], z_piston=Ra['plate_hi'], z_load=Ra['plate_hi'], bw=Ra['bw'],
              z_feed=(R['z_feed'] if a == 'z' else np.nan), z_perm=(R['z_perm'] if a == 'z' else np.nan),
              two_pist=(a == 'z'), phi_mf=Ra['phi_mf'], phi_vor=None, phi_cal=None, rho_s0=Ra['rho_s0'], axis=a)
    if L is None:
        return Rv, None
    La = L['ax'][a]
    Lv = dict(lvl=L['lvl'], eps=L['eps'], z=La['z'], ts=La['stress']['zz']['ts'], stress=La['stress'],
              halt_ts=L['halt_ts'], evol_from=L['evol_from'], plat=L['plat'][:len(La['stress']['zz']['ts'])],
              in_mem=La['in_mem'], interior=La['interior'], z_mem_lo=La['z_mem_lo'], z_mem_hi=La['z_mem_hi'],
              z_pist=La['plate_hi'], z_supp=La['plate_lo'], piston_z_at=La['hi_at'], support_z_at=La['lo_at'],
              support_pos=_track(L['plate_tab'], a + 'lo'), piston_pos=_track(L['plate_tab'], a + 'hi'),
              bw=La['bw'], wet=L.get('wet'), wetz=L['wetz'], phi_mf=La['phi_mf'], phi_vor=None, phi_cal=None,
              G={}, Dc=None, axis=a)
    return Rv, Lv


@contextmanager
def _axis(cfg, a):
    """Draw one of triaxial's figures along axis `a`: the saved file is tagged _along<a> and
    the profile coordinate in every label reads <a> instead of z."""
    orig_save, prev = tri._save, cfg._axis

    def _save(fig, c, stem, lvl=None):
        if a != 'z':
            for ax in fig.axes:
                for get, put in ((ax.get_xlabel, ax.set_xlabel), (ax.get_ylabel, ax.set_ylabel), (ax.get_title, ax.set_title)):
                    s = get()
                    s2 = s.replace('$z/L$', f'${a}/L$').replace('(z,t)', f'({a},t)').replace('(z)', f'({a})')
                    if s2 != s:
                        put(s2)
        return orig_save(fig, c, stem, lvl)

    tri._save, cfg._axis = _save, a
    try:
        yield
    finally:
        tri._save, cfg._axis = orig_save, prev


def _each_axis(cfg, R, L, fn, axes=None):
    """call triaxial figure `fn(cfg, Rv, Lv)` along every axis in `axes` (default: the first)."""
    out = []
    for a in (axes or cfg.AXES[:1]):
        if a not in R['ax'] or (L is not None and a not in L['ax']):
            print(f'{fn.__name__} along {a} skipped (no profiles)')
            continue
        with _axis(cfg, a):
            out.append(fn(cfg, *view(cfg, R, L, a)))
    return out[0] if len(out) == 1 else out


def _each_axis_sweep(cfg, R, levels, fn, axes=None):
    out = []
    for a in (axes or cfg.AXES[:1]):
        Ls = [L for L in levels if a in L['ax']]
        if a not in R['ax'] or not Ls:
            print(f'{fn.__name__} along {a} skipped (no profiles)')
            continue
        with _axis(cfg, a):
            Rv = view(cfg, R, None, a)[0]
            out.append(fn(cfg, Rv, [view(cfg, R, L, a)[1] for L in Ls]))
    return out[0] if len(out) == 1 else out


# ---- triaxial's profile figures along any axis (axes=None -> cfg.AXES[0]; axes=cfg.AXES -> all) ----
def fig_total_stress(cfg, R, L, axes=None):
    return _each_axis(cfg, R, L, tri.fig_total_stress, axes)


def fig_partial_stress(cfg, R, L, axes=None):
    return _each_axis(cfg, R, L, tri.fig_partial_stress, axes)


def fig_network_stress(cfg, R, L, axes=None):
    return _each_axis(cfg, R, L, tri.fig_network_stress, axes)


def fig_thermo_pressure(cfg, R, L, axes=None):
    return _each_axis(cfg, R, L, tri.fig_thermo_pressure, axes)


def fig_osmotic_pressure(cfg, R, L, axes=None):
    return _each_axis(cfg, R, L, tri.fig_osmotic_pressure, axes)


def fig_volfrac(cfg, R, L, axes=None):
    """mass-fraction solvent volume fraction rho_s / rho_bath in the column (the Voronoi and
    lambda-calibrated estimators of the triaxial notebooks are not computed for this geometry)."""
    return _each_axis(cfg, R, L, tri.fig_volfrac, axes)


def fig_wet_pistons(cfg, R, L):
    return _each_axis(cfg, R, L, tri.fig_wet_pistons, ('z',))


def fig_total_stress_sweep(cfg, R, levels, axes=None):
    return _each_axis_sweep(cfg, R, levels, tri.fig_total_stress_sweep, axes)


def fig_partial_stress_sweep(cfg, R, levels, axes=None):
    return _each_axis_sweep(cfg, R, levels, tri.fig_partial_stress_sweep, axes)


def fig_network_stress_sweep(cfg, R, levels, axes=None):
    return _each_axis_sweep(cfg, R, levels, tri.fig_network_stress_sweep, axes)


def fig_thermo_pressure_sweep(cfg, R, levels, axes=None):
    return _each_axis_sweep(cfg, R, levels, tri.fig_thermo_pressure_sweep, axes)


def fig_osmotic_pressure_sweep(cfg, R, levels, axes=None):
    return _each_axis_sweep(cfg, R, levels, tri.fig_osmotic_pressure_sweep, axes)


def fig_volfrac_sweep(cfg, R, levels, axes=None):
    return _each_axis_sweep(cfg, R, levels, tri.fig_volfrac_sweep, axes)


def fig_wet_pistons_sweep(cfg, R, levels):
    return _each_axis_sweep(cfg, R, levels, tri.fig_wet_pistons_sweep, ('z',))


# ===========================================================================
#  7. FIGURES -- SINGLE LEVEL
# ===========================================================================
def _save(fig, cfg, stem, lvl=None):
    return tri._save(fig, cfg, stem, lvl)


def fig_strain(cfg, R, levels, stem='strain_diagnostic'):
    """Volumetric strain vs step: the applied plate strain (bold), the network's own Rg
    (solid) and bounding-box (dashed) volumetric strains; dotted = target, shaded = plateau
    window.  Works for one level or a sweep (colour = level)."""
    fig, ax = plt.subplots(figsize=(10, 6.5), constrained_layout=True)
    any_ = False
    for i, L in enumerate(levels):
        S = L.get('strain_tab')
        if not S:
            continue
        any_ = True
        col = level_color(i) if len(levels) > 1 else WONG['blue']
        st = S['step']
        ax.plot(st, S['eps_vol_plate'], '-', color=col, lw=3.2, alpha=0.30, label=fr'plates (target {L["lvl"]})')
        if np.isfinite(R['V0_rg']):
            ax.plot(st, 1.0 - S['V_rg'] / R['V0_rg'], '-', color=col, lw=2.2, label=r'network $1-V_{Rg}/V_{Rg,0}$')
        if np.isfinite(R['V0_bb']):
            ax.plot(st, rolling_mean(1.0 - S['V_bb'] / R['V0_bb'], cfg.roll_win), '--', color=col, lw=1.5, alpha=0.7,
                    label=r'network $1-V_{BB}/V_{BB,0}$ (rolling mean)')
        ax.axvspan(L['halt_ts'], st[-1], color=col, alpha=0.06)
        ax.axhline(L['eps_nominal'], color=col, ls=':', lw=1.0, alpha=0.6)
        ax.annotate(f"plateau: plates {sig(L['strain']['plate'])}  Rg {sig(L['strain'].get('rg', np.nan))}  "
                    f"BB {sig(L['strain'].get('bb', np.nan))}", (st[-1], S['eps_vol_plate'][-1]),
                    textcoords='offset points', xytext=(-6, 9), ha='right', va='bottom', fontsize=10, color=col,
                    bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.8))
    if not any_:
        plt.close(fig)
        print('strain diagnostic skipped (no strain_vol files)')
        return None
    ax.set_xlabel('time step')
    ax.set_ylabel(r'volumetric strain  $\varepsilon_{vol}=1-V/V_0$')
    ax.set_title('Strain diagnostic: bold = applied (plates), solid = network $V_{Rg}$, dashed = network $V_{BB}$,\n'
                 'dotted = target, shaded = plateau window', fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=10)
    return _save(fig, cfg, stem, levels[0]['lvl'] if len(levels) == 1 else None)


def fig_plates(cfg, R, L):
    """Plate pressures P = load / face area vs step: (a) the six plates (rolling means) and
    their mean with the auto plateau window, (b) log of the mean (relaxation view)."""
    PP = L.get('plates_raw', L.get('plates'))
    if PP is None:
        print('plate figure skipped (no plate_force file)')
        return None
    st = PP['step']
    fig, (axL, axG) = plt.subplots(1, 2, figsize=(17, 6), constrained_layout=True)
    fig.suptitle(f'Plate pressure histories ($P$ = compressive load / face area; area = {PP["area_mode"]})  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    for f in PLATES:
        c, ls = PLATE_STYLE[f]
        axL.plot(st, rolling_mean(PP[f], cfg.roll_win), ls, color=c, lw=1.4, alpha=0.75, label=f'plate {f[0]}-{f[1:]}')
    Pm = rolling_mean(PP['mean'], cfg.roll_win)
    axL.plot(st, PP['mean'], '-', color='k', lw=0.8, alpha=0.25)
    axL.plot(st, Pm, '-', color='k', lw=2.8, label=f'mean of the six (rolling {cfg.roll_win})')
    if 'plates' in L and 'plates_raw' in L:
        axL.plot(L['plates']['step'], L['plates']['mean'], '-', color=WONG['orange'], lw=1.6, alpha=0.9, label='mean, LMP block-avg')
    if np.isfinite(R.get('P_ref', np.nan)):
        axL.axhline(R['P_ref'], color='0.4', ls=':', lw=1.4, label=fr"reference preload $P_{{\rm ref}}={sig(R['P_ref'])}$")
    pos = Pm > 0
    axG.plot(st[pos], np.log(Pm[pos]), '-', color='k', lw=2.6, label=f'$\\ln$ mean (rolling {cfg.roll_win})')
    for ax in (axL, axG):
        if 'PF' in L:
            ax.axvspan(L['PF']['step0'], float(st[-1]), color=WONG['green'], alpha=0.10, label='plateau window (auto)')
        ax.set_xlabel('time step')
        ax.grid(alpha=0.3)
    axL.axhline(0, color='k', ls='--', lw=0.8, alpha=0.4)
    axL.set_ylabel(r'$P = N/A$  (LJ / $\sigma^2$)')
    axL.set_title('(a) plate pressures')
    axG.set_ylabel(r'$\ln P$')
    axG.set_title('(b) log mean plate pressure (relaxation view)')
    smart_legend(axL, fontsize=11)
    smart_legend(axG, fontsize=12)
    if 'PF' in L:
        p = L['PF']
        annotate_box(axL, f"plateau $\\langle P\\rangle$ = {fmt_val_unc(p['mean'], 0.5 * (p['hi'] - p['lo']))}\n"
                          f"last {p['frac']:.0%}, n={p['n']}, block={p['block']}", loc='upper right', fontsize=12)
    return _save(fig, cfg, 'plate_pressure_history', L['lvl'])


def _pt(ax, x, D, key, marker, color, ms, label, hollow=False):
    kw = dict(mfc='none', alpha=0.7, lw=1.5, capsize=5) if hollow else dict(lw=2.5, capsize=8)
    ax.errorbar([x], [D[key]], yerr=[[D[key] - D[key + '_lo']], [D[key + '_hi'] - D[key]]], fmt=marker, ms=ms,
                color=color, label=label, **kw)


def fig_K(cfg, R, L):
    """Bulk modulus, two panes (the layout of triaxial.fig_M).
    (a) diagnostic: the network K from each profile axis and their mean, the plate K, and the
        absolute stress / eps_vol values (hollow) with the eps_vol = 0 readings each subtracts;
    (b) presentation: the two K estimates alone, with their CIs."""
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
    ci = int(cfg.ci_level * 100)
    sub = bool(cfg.K_SUBTRACT_REF)
    fig.suptitle('Bulk modulus' + (' (increment from $\\varepsilon_{vol}=0$)' if sub else '') + '   |   ' + cfg.RUN_ID
                 + f"   |   $\\varepsilon_{{vol}} = {sig(L['eps'])}$ ({cfg.K_STRAIN})", fontsize=14, fontweight='bold')
    axs = [a for a in cfg.AXES if a in L['ax']]
    has_p = 'K_pl' in L
    # ---- (a) diagnostic ---------------------------------------------------------
    for k, a in enumerate(axs):
        _pt(axA, -0.30 + 0.15 * k, L['ax'][a], 'K', 'o', AX_COLOR[a], 9, f"profile along {a}:  $K = {sig(L['ax'][a]['K'])}$", hollow=False)
    _pt(axA, 0.0 + 0.15 * len(axs) - 0.30, L, 'K_net', 'o', 'k', 13,
        f"network (mean)  $K = {sig(L['K_net'])}$\n{ci}% CI [{sig(L['K_net_lo'])}, {sig(L['K_net_hi'])}]")
    axA.axhline(L['K_net'], color='k', ls='--', lw=1.2, alpha=0.5)
    if sub:
        _pt(axA, 0.45, L, 'K_net_abs', 'o', 'k', 10,
            f"absolute $\\mathit{{\\Pi}}/\\varepsilon = {sig(L['K_net_abs'])}$\n(ref $\\mathit{{\\Pi}} = {L['K_net_ref']:+.4f}$ subtracted)", hollow=True)
    if has_p:
        _pt(axA, 1.0, L, 'K_pl', 's', WONG['vermillion'], 13,
            f"plates  $K = {sig(L['K_pl'])}$\n{ci}% CI [{sig(L['K_pl_lo'])}, {sig(L['K_pl_hi'])}]"
            + ('\n(absolute: no $P_{\\rm ref}$ file)' if (sub and L['K_pl_ref'] != 'measured') else ''))
        axA.axhline(L['K_pl'], color=WONG['vermillion'], ls='--', lw=1.2, alpha=0.5)
        if sub and L['K_pl_ref'] == 'measured':
            _pt(axA, 1.15, L, 'K_pl_abs', 's', WONG['vermillion'], 10,
                f"absolute $P/\\varepsilon = {sig(L['K_pl_abs'])}$\n(ref $P = {L['P_ref']:+.4f}$ subtracted)", hollow=True)
    axA.set_xticks([0, 1])
    axA.set_xticklabels([r'network' + '\n' + r"$(\langle\mathit{\Pi}\rangle_{\rm int}-\mathit{\Pi}_{\rm ref})/\varepsilon_{vol}$",
                         r'plates' + '\n' + r'$(\bar P-\bar P_{\rm ref})/\varepsilon_{vol}$'] if sub else
                        [r"network ($\langle\mathit{\Pi}\rangle_{\rm int}/\varepsilon_{vol}$)", r'plates ($\bar P/\varepsilon_{vol}$)'],
                        fontsize=13)
    axA.set_ylabel(r'$K$  (LJ units)')
    axA.set_title('(a) per profile axis, and the absolute values (hollow)', fontsize=13)
    axA.set_xlim(-0.6, 1.6)
    axA.grid(axis='y', alpha=0.3)
    smart_legend(axA, fontsize=10)
    # ---- (b) presentation -------------------------------------------------------
    _pt(axB, 0, L, 'K_net', 'o', WONG['blue'], 14,
        f"network  $K = {sig(L['K_net'])}$\n{ci}% CI [{sig(L['K_net_lo'])}, {sig(L['K_net_hi'])}]")
    axB.axhline(L['K_net'], color=WONG['blue'], ls='--', lw=1.2, alpha=0.5)
    if has_p:
        _pt(axB, 1, L, 'K_pl', 's', WONG['vermillion'], 14,
            f"plates  $K = {sig(L['K_pl'])}$\n{ci}% CI [{sig(L['K_pl_lo'])}, {sig(L['K_pl_hi'])}]")
        axB.axhline(L['K_pl'], color=WONG['vermillion'], ls='--', lw=1.2, alpha=0.5)
    axB.set_xticks([0, 1])
    axB.set_xticklabels([r'network' + '\n' + r"$\Delta\langle\mathit{\Pi}\rangle_{\rm int}/\varepsilon_{vol}$",
                         r'plates' + '\n' + r'$\Delta\bar P/\varepsilon_{vol}$'] if sub else
                        [r"network ($\langle\mathit{\Pi}\rangle_{\rm int}/\varepsilon_{vol}$)", r'plates ($\bar P/\varepsilon_{vol}$)'],
                        fontsize=16)
    axB.set_ylabel(r'$K$  (LJ units)')
    axB.set_title('(b) bulk modulus, two independent estimates', fontsize=13)
    axB.set_xlim(-0.5, 1.5)
    axB.grid(axis='y', alpha=0.3)
    smart_legend(axB, fontsize=13)
    return _save(fig, cfg, 'K_comparison', L['lvl'])


def fig_isotropy(cfg, R, L):
    """Is the loading isotropic?  (a) the increment of each normal network-stress component
    from its reference, per profile axis (equal bars = no deviator; their mean is the Pi
    increment); (b) the six plate pressure increments.  Replaces the anisotropy-ratio and G
    figures of the triaxial notebooks: under isotropic compression the deviator is zero by
    symmetry, so there is no G in this measurement -- a deviator here is a loading defect."""
    axs = [a for a in cfg.AXES if a in L['ax']]
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    fig.suptitle(f"Isotropy of the response   |   {cfg.RUN_ID}   |   $\\varepsilon_{{vol}} = {sig(L['eps'])}$", fontsize=14, fontweight='bold')
    w = 0.8 / max(len(axs), 1)
    comps = ('xx', 'yy', 'zz')
    for k, a in enumerate(axs):
        v = [L['ax'][a]['dsig'][c] for c in comps]
        axA.bar(np.arange(3) + (k - 0.5 * (len(axs) - 1)) * w, v, w * 0.9, color=AX_COLOR[a], alpha=0.85, label=f'profile along {a}')
        m = float(np.mean(v))
        axA.plot([-0.5, 2.5], [m, m], ':', color=AX_COLOR[a], lw=1.4)
    axA.set_xticks(range(3))
    axA.set_xticklabels([r"$\Delta\sigma'_{%s}$" % c for c in comps], fontsize=16)
    axA.axhline(0, color='k', lw=0.8, alpha=0.5)
    axA.set_ylabel(r"$\langle\sigma'_{ii}\rangle_{\rm int}-\sigma'_{ii,\rm ref}$  (LJ)")
    axA.set_title(r"(a) network-stress increments (dotted: their mean $=\Delta\mathit{\Pi}$)", fontsize=13)
    axA.grid(axis='y', alpha=0.3)
    smart_legend(axA, fontsize=11)
    annotate_box(axA, f"largest deviator = {sig(L['dev_max'])}\n= {L['dev_rel']:.1%} of $\\Delta\\mathit{{\\Pi}}$", loc='lower right', fontsize=12)
    if 'P_face' in L:
        ref = R.get('P_ref_face', {f: 0.0 for f in PLATES}) if cfg.K_SUBTRACT_REF else {f: 0.0 for f in PLATES}
        v = [L['P_face'][f] - ref[f] for f in PLATES]
        axB.bar(range(6), v, 0.7, color=[PLATE_STYLE[f][0] for f in PLATES], alpha=0.85)
        m = float(np.mean(v))
        axB.axhline(m, color='k', ls=':', lw=1.4, label=fr'mean $=\Delta\bar P={sig(m)}$')
        axB.set_xticks(range(6))
        axB.set_xticklabels([f'{f[0]}-{f[1:]}' for f in PLATES], fontsize=14)
        axB.axhline(0, color='k', lw=0.8, alpha=0.5)
        smart_legend(axB, fontsize=12)
    else:
        axB.text(0.5, 0.5, 'no plate_force file', ha='center', va='center', transform=axB.transAxes)
    axB.set_ylabel(r'$P_{\rm plate}-P_{\rm plate,ref}$  (LJ)' if cfg.K_SUBTRACT_REF else r'$P_{\rm plate}$  (LJ)')
    axB.set_title('(b) plate pressure increments, plate by plate', fontsize=13)
    axB.grid(axis='y', alpha=0.3)
    return _save(fig, cfg, 'isotropy_check', L['lvl'])


def fig_solvent_expelled(cfg, R, L):
    """Solvent expelled (wet-piston displacement) next to the volume the plates closed."""
    W = L.get('wet')
    if W is None or W.get('dV_total') is None:
        print('solvent-expelled figure skipped (no permeation_<level> file)')
        return None
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    st = W['exp_step']
    ax.plot(st, W['dV_total'], '-', color=WONG['green'], lw=2.4, label=r'expelled  $A\,(\Delta z_{\rm feed}-\Delta z_{\rm perm})$')
    ax.plot(st, W['dV_feed'], '--', color=WONG['vermillion'], lw=1.6, label='taken up on the feed side')
    ax.plot(st, W['dV_perm'], '--', color=WONG['blue'], lw=1.6, label='taken up on the permeate side')
    S = L.get('strain_tab')
    if S:
        ax.plot(S['step'], S['eps_vol_plate'] * R['V0_gap'], ':', color='k', lw=2.0,
                label=r'volume closed by the plates  $\varepsilon_{vol}\,V_0$')
    ax.axhline(0, color='k', ls=':', lw=1, alpha=0.5)
    ax.set_xlabel('time step')
    ax.set_ylabel(r'volume  ($\sigma^3$)')
    ax.set_title(f"Solvent expelled since seating, $\\varepsilon_{{vol}}={L['lvl']}$  |  final {sig(W['dV_total'][-1])} $\\sigma^3$ "
                 f"= {sig(W['dV_total'][-1] / R['V0_gap'])} of $V_0$ (plates closed {sig(L['strain']['plate'])})", fontsize=13)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=12)
    return _save(fig, cfg, 'solvent_expelled', L['lvl'])


def _dc_panel(ax, F, a, title, legend=True):
    zff = np.linspace(0.0, 1.0, 300)
    cmap = plt.cm.viridis
    norm = Normalize(vmin=F['ts'][F['early'][0]], vmax=F['ts'][F['early'][-1]])
    u0_b = F['u_IC'](F['zf'])
    for i in F['early']:
        c = cmap(norm(F['ts'][i]))
        ax.plot(F['zf'], u0_b + F['uhat'][i][F['idx']], 'o', color=c, ms=3, alpha=0.35)
        ax.plot(zff, F['u_model'](zff, F['t_lj'][i]), '-', color=c, lw=1.8)
    dlL, fs = F['DL'] / F['L'], F['f_sup']
    ax.plot([0, 1], [fs * dlL, (fs - 1.0) * dlL], 'k:', lw=1.8, label=r'affine ($t\to\infty$)')
    ax.set(xlabel=rf'$\zeta=({a}-{a}_{{\rm lo}})/L$', ylabel=rf'$u_{a}/L$', xlim=(0, 1))
    ax.set_title(title, fontsize=14, color=AX_COLOR[a])
    ax.grid(alpha=0.3)
    ax._tri_has_colorbar = True
    if legend:
        smart_legend(ax, fontsize=11)
    return norm, cmap


def fig_Dc(cfg, R, L):
    """The held-cube consolidation fit (fit_Dc_cube): u_a(zeta_a, t)/L_a along each profile
    axis, data + the joint model, ONE D_c for the three panels.  Every face is drained at the
    bath pressure and the displacement vanishes at the centre, so only the antisymmetric
    sin(2 l pi zeta) modes appear along an axis; their decay carries the transverse
    cos(2 m pi zeta_b) cos(2 n pi zeta_c) drainage of the other two axes."""
    F = L.get('Dc')
    if F is None:
        print('D_c figure skipped (no fit)')
        return None
    axs = list(F['ax'])
    fig, axes = plt.subplots(1, len(axs), figsize=(8.5 * len(axs), 6.5), constrained_layout=True, squeeze=False)
    for ax, a in zip(axes[0], axs):
        Fa = F['ax'][a]
        one = L.get('Dc_axis', {}).get(a)
        norm, cmap = _dc_panel(ax, Fa, a, rf"along {a}:  $L={sig(Fa['L'])}$, $\Delta L={sig(Fa['DL'])}$"
                               + (rf"  (1-D fit alone: $D_c={sig(one['Dc'])}$)" if one else ''))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=list(axes[0]), fraction=0.015, pad=0.02).set_label('timestep')
    fig.suptitle(f'Held-cube consolidation fit (level _c{L["lvl"]})  |  {cfg.sim_name}  |  '
                 f"$D_c = {sig(F['Dc'])}\\ \\sigma^2/\\tau$,  $R^2 = {sig(F['R2'])}$  "
                 f"(modes: $l \\leq {cfg.DC_N_MODES}$ along the axis, $m, n \\leq {cfg.DC_N_T}$ across; zero at the centre)",
                 fontsize=12, fontweight='bold')
    return _save(fig, cfg, 'Dc_consolidation_fit', L['lvl'])


def fig_closure(cfg, R, L):
    """K next to the triaxial M (cfg.M_REF): G = (3/4)(M - K) and kappa = D_c / M."""
    if cfg.M_REF is None:
        print('closure figure skipped (set M_REF in the Config to the M of the triaxial notebook)')
        return None
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    M = float(cfg.M_REF)
    _pt(ax, 0, L, 'K_net', 'o', WONG['blue'], 13, f"$K$ (network) = {sig(L['K_net'])}")
    if 'K_pl' in L:
        _pt(ax, 0.2, L, 'K_pl', 's', WONG['vermillion'], 11, f"$K$ (plates) = {sig(L['K_pl'])}")
    ax.plot([1], [M], 'D', ms=12, color='k', label=f'$M$ (triaxial, M_REF) = {sig(M)}')
    _pt(ax, 2, L, 'G_from_MK', '^', WONG['green'], 13, fr"$G=\frac{{3}}{{4}}(M-K)$ = {sig(L['G_from_MK'])}")
    if cfg.G_REF is not None:
        ax.plot([2.2], [float(cfg.G_REF)], 'v', ms=11, mfc='none', color=WONG['green'], label=f'$G$ (G_REF) = {sig(float(cfg.G_REF))}')
    ax.axhline(0, color='k', lw=0.8, alpha=0.5)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['$K$', '$M$', '$G$'], fontsize=16)
    ax.set_ylabel('modulus  (LJ units)')
    ax.set_title(f"Closure $M = K + \\frac{{4}}{{3}}G$   |   {cfg.RUN_ID}   |   $\\varepsilon_{{vol}}={sig(L['eps'])}$"
                 + (f"\n$\\kappa=D_c/M$ = {sig(L['kappa'])}" if 'kappa' in L else ''), fontsize=14)
    ax.set_xlim(-0.5, 2.7)
    ax.grid(axis='y', alpha=0.3)
    smart_legend(ax, fontsize=12)
    return _save(fig, cfg, 'moduli_closure', L['lvl'])


# ===========================================================================
#  8. FIGURES -- SWEEP
# ===========================================================================
def level_handles(levels, ref=False):
    h = [Line2D([0], [0], color=level_color(i), lw=3, label=fr'$\varepsilon_{{vol}}={L["eps"]:.3f}$') for i, L in enumerate(levels)]
    if ref:
        h.insert(0, Line2D([0], [0], color='k', ls='--', lw=2, label=r'reference ($\varepsilon_{vol}=0$)'))
    return h


def fig_plates_sweep(cfg, R, levels):
    fig, (axP, axG) = plt.subplots(1, 2, figsize=(18, 6), constrained_layout=True)
    fig.suptitle(f'Mean plate pressure histories, all levels  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    for i, L in enumerate(levels):
        PP = L.get('plates_raw', L.get('plates'))
        if PP is None:
            continue
        col = level_color(i)
        st, Pm = PP['step'], rolling_mean(PP['mean'], cfg.roll_win)
        axP.plot(st, PP['mean'], '-', color=col, lw=0.8, alpha=0.20)
        axP.plot(st, Pm, '-', color=col, lw=2.4, alpha=0.95)
        pos = Pm > 0
        axG.plot(st[pos], np.log(Pm[pos]), '-', color=col, lw=2.4, alpha=0.95)
        if 'PF' in L:
            axP.axvspan(L['PF']['step0'], float(st[-1]), color=col, alpha=0.06)
            axP.axhline(L['PF']['mean'], color=col, ls=':', lw=1.2, alpha=0.7)
    axP.axhline(0, color='k', ls='--', lw=0.8, alpha=0.4)
    axP.set_xlabel('time step')
    axP.set_ylabel(r'$\bar P$ = mean of the six $N/A$  (LJ / $\sigma^2$)')
    axP.set_title('(a) mean plate pressure (rolling mean; dotted = plateau)')
    axG.set_xlabel('time step')
    axG.set_ylabel(r'$\ln \bar P$')
    axG.set_title('(b) log mean plate pressure (relaxation view)')
    for ax in (axP, axG):
        ax.grid(alpha=0.3)
        smart_legend(ax, handles=level_handles(levels), fontsize=11)
    return _save(fig, cfg, 'sweep_plate_pressure_history')


def _ser(ax, e, Ls, key, marker, color, label, hollow=False, dx=0.0):
    v = np.array([L[key] for L in Ls])
    lo = np.array([L[key + '_lo'] for L in Ls])
    hi = np.array([L[key + '_hi'] for L in Ls])
    kw = dict(mfc='none', alpha=0.7, lw=1.5, ms=8, capsize=4, ls=':') if hollow else dict(lw=2, ms=10, capsize=6, ls='-')
    ax.errorbar(np.asarray(e) + dx, v, yerr=[v - lo, hi - v], fmt=marker, color=color, label=label, **kw)


def fig_K_sweep(cfg, R, levels):
    """Bulk modulus vs volumetric strain, two panes (the layout of triaxial.fig_M_sweep):
    (a) diagnostic (per profile axis + the absolute values, hollow), (b) presentation."""
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    sub = bool(cfg.K_SUBTRACT_REF)
    fig.suptitle('Bulk modulus across the sweep' + (' (increment from $\\varepsilon_{vol}=0$)' if sub else '') + '   |   ' + cfg.sim_name,
                 fontsize=13, fontweight='bold')
    eps = np.array([L['eps'] for L in levels])
    hp = [L for L in levels if 'K_pl' in L]
    ep = np.array([L['eps'] for L in hp])
    for a in cfg.AXES:
        Ls = [L for L in levels if a in L['ax']]
        if Ls:
            e = np.array([L['eps'] for L in Ls])
            v = np.array([L['ax'][a]['K'] for L in Ls])
            axA.plot(e, v, 'o--', color=AX_COLOR[a], ms=6, lw=1.2, alpha=0.8, label=f'profile along {a}')
    _ser(axA, eps, levels, 'K_net', 'o', 'k', 'network (mean of the axes)')
    if sub:
        _ser(axA, eps, levels, 'K_net_abs', 'o', 'k', r'network, absolute $\langle\mathit{\Pi}\rangle_{\rm int}/\varepsilon_{vol}$', hollow=True, dx=0.002)
    if hp:
        _ser(axA, ep, hp, 'K_pl', 's', WONG['vermillion'], 'plates')
        if sub and all(L['K_pl_ref'] == 'measured' for L in hp):
            _ser(axA, ep, hp, 'K_pl_abs', 's', WONG['vermillion'], r'plates, absolute $\bar P/\varepsilon_{vol}$', hollow=True, dx=-0.002)
    if sub:
        txt = f"$\\varepsilon_{{vol}}=0$ readings subtracted:\n  $\\mathit{{\\Pi}}_{{\\rm ref}} = {R['Pi_ref']:+.4f}$"
        if np.isfinite(R.get('P_ref', np.nan)):
            txt += f"\n  $\\bar P_{{\\rm ref}} = {R['P_ref']:+.4f}$"
        annotate_box(axA, txt, loc='lower right', fontsize=11)
    axA.set_title('(a) per profile axis, and the absolute values (hollow)', fontsize=13)
    _ser(axB, eps, levels, 'K_net', 'o', WONG['blue'], r"network  $\Delta\langle\mathit{\Pi}\rangle_{\rm int}/\varepsilon_{vol}$")
    if hp:
        _ser(axB, ep, hp, 'K_pl', 's', WONG['vermillion'], r'plates  $\Delta\bar P/\varepsilon_{vol}$')
    axB.set_title('(b) bulk modulus per level, two independent estimates', fontsize=13)
    for ax in (axA, axB):
        ax.set_xlabel(r'volumetric strain  $\varepsilon_{vol}$')
        ax.set_ylabel(r'$K$  (LJ units)')
        ax.grid(alpha=0.3)
    smart_legend(axA, fontsize=10)
    smart_legend(axB, fontsize=13)
    return _save(fig, cfg, 'sweep_modulus')


def fig_stress_strain_sweep(cfg, R, levels):
    """Pi and the mean plate pressure vs volumetric strain, INCLUDING each estimator's
    eps_vol = 0 reading (hollow).  The dotted line of each series is anchored at its own
    eps_vol = 0 reading and its slope is the through-anchor least-squares fit of the first
    cfg.SS_FIT_NPTS levels (the linear regime): that slope is the small-strain K."""
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    n_fit = int(max(1, cfg.SS_FIT_NPTS))
    fits = []

    def _series(e, s, err, color, marker, name, s0, err0):
        ax.errorbar(e, s, yerr=err, fmt=marker + '-', lw=2, ms=9, color=color, capsize=6, label=name)
        e, s = np.asarray(e, float), np.asarray(s, float)
        anchor = float(s0) if np.isfinite(s0) else 0.0
        if np.isfinite(s0):
            ax.errorbar([0.0], [s0], yerr=err0, fmt=marker, ms=9, mfc='none', color=color, capsize=6)
        o = np.argsort(e)
        ef, sf = e[o][:n_fit], s[o][:n_fit]
        slope = float(np.sum(ef * (sf - anchor)) / np.sum(ef ** 2))
        xs = np.linspace(0, max(e) * 1.05, 20)
        ax.plot(xs, anchor + slope * xs, ls=':', lw=1.5, color=color, alpha=0.8)
        fits.append((name.split()[0], slope, len(ef)))
        return slope

    eps = np.array([L['eps'] for L in levels])
    exc = np.array([L.get('relax', {}).get('excess', 0.0) for L in levels])
    sn = np.array([L['Pi_int'] for L in levels])
    hn = np.array([L['Pi_half'] for L in levels])
    R['K_small_net'] = _series(eps, sn, [hn + exc, hn], WONG['blue'], 'o', 'network $\\langle\\mathit{\\Pi}\\rangle_{\\rm int}$ (plateau)',
                               R['Pi_ref'], R['Pi_ref_half'])
    hp = [L for L in levels if 'PF' in L]
    if hp:
        ep = np.array([L['eps'] for L in hp])
        Pp = np.array([L['PF']['mean'] for L in hp])
        ex = np.array([L.get('relax', {}).get('excess', 0.0) for L in hp])
        pref = R.get('P_ref', np.nan)
        R['K_small_pl'] = _series(ep, Pp, [Pp - [L['PF']['lo'] for L in hp] + ex, [L['PF']['hi'] for L in hp] - Pp],
                                  WONG['vermillion'], 's', 'plates $\\bar P$ (plateau)', pref,
                                  (float(R['P_ref_hi'] - R['P_ref_lo']) / 2 if np.isfinite(pref) else None))
    ax.plot([], [], 'o', mfc='none', color='0.4', label=r'hollow = $\varepsilon_{vol}=0$ reading')
    ax.plot([], [], ls=':', lw=1.5, color='0.4', label=f'dotted = through the $\\varepsilon_{{vol}}=0$ reading, slope from the first {n_fit} levels')
    if np.any(exc > 0):
        ax.plot([], [], ' ', label='lower bars include the unrelaxed-hold excess')
    ax.set_title('Stress vs volumetric strain, small-strain slopes through $\\varepsilon_{vol}=0$:  ' +
                 ',  '.join(f"$K_{{\\rm {n[:4]}}}\\approx{sig(s)}$ ({k} levels)" for n, s, k in fits) + '\n' + cfg.sim_name, fontsize=12)
    ax.set_xlabel(r'volumetric strain  $\varepsilon_{vol}$')
    ax.set_ylabel(r'plateau stress  (LJ)')
    ax.grid(alpha=0.3)
    ax.set_xlim(left=-0.01)
    smart_legend(ax, fontsize=12)
    return _save(fig, cfg, 'sweep_stress_strain')


def fig_isotropy_sweep(cfg, R, levels):
    """Deviator of the network-stress increment vs volumetric strain: each normal component
    minus the mean, per profile axis (0 = isotropic response)."""
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    eps = np.array([L['eps'] for L in levels])
    for a in cfg.AXES:
        Ls = [L for L in levels if a in L['ax']]
        if not Ls:
            continue
        e = np.array([L['eps'] for L in Ls])
        for c, ls in (('xx', '-'), ('yy', '--'), ('zz', ':')):
            ax.plot(e, [L['ax'][a]['dev'][c] for L in Ls], ls, marker='o', ms=5, lw=1.6, color=AX_COLOR[a])
    ax.plot(eps, [L['Pi_int'] - (L['K_net_ref'] if cfg.K_SUBTRACT_REF else 0.0) for L in levels], '-', color='0.6', lw=3, alpha=0.5)
    ax.axhline(0, color='k', lw=0.8, alpha=0.5)
    h = [Line2D([0], [0], color=AX_COLOR[a], lw=3, label=f'profile along {a}') for a in cfg.AXES]
    h += [Line2D([0], [0], color='0.3', ls=ls, lw=2, label=r"$\Delta\sigma'_{%s}-\Delta\mathit{\Pi}$" % c)
          for c, ls in (('xx', '-'), ('yy', '--'), ('zz', ':'))]
    h += [Line2D([0], [0], color='0.6', lw=3, alpha=0.5, label=r'$\Delta\mathit{\Pi}$ (scale)')]
    ax.set_xlabel(r'volumetric strain  $\varepsilon_{vol}$')
    ax.set_ylabel(r'network-stress deviator  (LJ)')
    ax.set_title('Isotropy of the response across the sweep (0 = no deviator)', fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, handles=h, fontsize=11)
    return _save(fig, cfg, 'sweep_isotropy')


def fig_Dc_sweep(cfg, R, levels):
    """(a) D_c vs volumetric strain (the held-cube fit; hollow dotted = the 1-D fits per axis);
    then the cube fit along cfg.AXES[0] of every level."""
    hd = [L for L in levels if L.get('Dc') is not None]
    if not hd:
        print('D_c figure skipped (no level produced a fit)')
        return None
    a0 = next(a for a in cfg.AXES if any(a in L['Dc']['ax'] for L in hd))
    hp = [L for L in hd if a0 in L['Dc']['ax']]
    fig, axes = plt.subplots(1, len(hp) + 1, figsize=(7 + 6.5 * len(hp), 6), constrained_layout=True)
    fig.suptitle(f'Cooperative diffusivity across the sweep  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    ax = axes[0]
    ax.plot([L['eps'] for L in hd], [L['Dc']['Dc'] for L in hd], 'o-', color='k', lw=2.2, ms=9, label='held-cube fit (3 axes jointly)')
    for a in cfg.AXES:
        Ls = [L for L in levels if a in L.get('Dc_axis', {})]
        if Ls:
            ax.plot([L['eps'] for L in Ls], [L['Dc_axis'][a]['Dc'] for L in Ls], 'o:', mfc='none', color=AX_COLOR[a], lw=1.2, ms=7,
                    alpha=0.8, label=f'1-D fit along {a} (apparent)')
    ax.axhline(cfg.DC_SLOW_REF / 4.0, color='0.4', ls=':', lw=1.5,
               label=f"deck hold-sizing $D_c$ = {sig(cfg.DC_SLOW_REF / 4)}  (= {sig(cfg.DC_SLOW_REF)} in the deck's $L^2/\\pi^2 D_c$ convention)")
    ax.margins(x=0.12, y=0.15)
    ax.set_xlabel(r'volumetric strain  $\varepsilon_{vol}$')
    ax.set_ylabel(r'$D_c$  ($\sigma^2/\tau$)')
    ax.set_title(r'(a) $D_c$ vs volumetric strain', fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=11)
    for k, L in enumerate(hp):
        F = L['Dc']
        _dc_panel(axes[k + 1], F['ax'][a0], a0, fr"({'bcdefgh'[k]}) $\varepsilon_{{vol}}={L['lvl']}$, along {a0}:  $D_c={sig(F['Dc'])}$, $R^2={sig(F['R2'])}$",
                  legend=False)
    return _save(fig, cfg, 'sweep_Dc')


def fig_closure_sweep(cfg, R, levels):
    """K, the triaxial M (M_REF) and G = (3/4)(M - K) vs volumetric strain."""
    if cfg.M_REF is None:
        print('closure figure skipped (set M_REF in the Config to the M of the triaxial notebook)')
        return None
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    eps = np.array([L['eps'] for L in levels])
    _ser(ax, eps, levels, 'K_net', 'o', WONG['blue'], '$K$ (network)')
    _ser(ax, eps, levels, 'G_from_MK', '^', WONG['green'], r'$G=\frac{3}{4}(M-K)$')
    ax.axhline(float(cfg.M_REF), color='k', ls='--', lw=1.5, label=f'$M$ (triaxial, M_REF) = {sig(float(cfg.M_REF))}')
    if cfg.G_REF is not None:
        ax.axhline(float(cfg.G_REF), color=WONG['green'], ls=':', lw=1.5, label=f'$G$ (G_REF) = {sig(float(cfg.G_REF))}')
    ax.set_xlabel(r'volumetric strain  $\varepsilon_{vol}$')
    ax.set_ylabel('modulus  (LJ units)')
    ax.set_title(r'Closure $M = K + \frac{4}{3}G$ across the sweep', fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, fontsize=12)
    return _save(fig, cfg, 'sweep_moduli_closure')


# ===========================================================================
#  9. SUMMARIES
# ===========================================================================
def print_hold_check(cfg, levels):
    """Was each level held long enough?  (triaxial.hold_adequacy on the held-cube fit.)"""
    hd = [L for L in levels if L.get('Dc') is not None]
    if not hd:
        return
    print(f'\nHOLD-ADEQUACY CHECK  (fit: tau_1 = L_max^2/(4 pi^2 D_c), the slowest held-cube mode (1,0,0) along the longest '
          f'edge; slow: the deck\'s sizing formula L^2/(pi^2 Dc_est); residual = mean excess stress over the last '
          f'{cfg.plateau_frac:.0%} of the hold)')
    for L in hd:
        F = L['Dc']
        print(f"  level _c{L['lvl']}:  L = " + ' x '.join(f"{F['ax'][a]['L']:.1f}" for a in F['ax'])
              + f" sigma   hold T = {F['hold_T']:.0f} tau = {F['hold_T'] / cfg.dt_lj / 1e6:.2f}M steps")
        for tag, h in F['hold_check'].items():
            flag = '' if h['ok'] else '   <-- TOO SHORT'
            print(f"     {tag:<4s} D_c={h['Dc']:.3f}:  tau_1 = {h['tau1']:.0f} tau = {h['tau1'] / cfg.dt_lj / 1e6:.2f}M steps"
                  f"  |  held {F['hold_T'] / h['tau1']:.2f} tau_1  ->  residual {h['avg'] * 100:5.2f}%  |  "
                  f"{cfg.DC_TARGET_RESID:.0%} needs {h['need']:.1f} tau_1{flag}")
        if 'K_pl' in L:
            print(f"     observed K_plates/K_network - 1 = {L['K_pl'] / L['K_net'] - 1:+.1%}")
    print('  (the direct measure of what is still unrelaxed is delta_sys, the plate tail fit)')


def print_summary(cfg, levels):
    """One line per level with the headline numbers, then the hold check."""
    ci = int(cfg.ci_level * 100)
    how = 'K = increment from the eps_vol = 0 reference' if cfg.K_SUBTRACT_REF else 'K = absolute stress / eps_vol'
    print(f'\nSUMMARY  ({cfg.sim_name}; {ci}% CIs; {how}; eps_vol = {cfg.K_STRAIN}'
          + ('; K lower bounds include delta_sys' if cfg.RELAX_SYS else '') + ')')
    print(f"{'target':>7s} {'eps_vol':>8s} {'eps_Rg':>8s} {'K_net':>22s} {'K_net(Rg)':>10s} {'K_plates':>22s} {'dev/dPi':>8s} {'D_c':>10s} {'delta_sys':>10s}"
          + (f" {'G=3(M-K)/4':>11s}" if cfg.M_REF is not None else ''))
    for L in levels:
        c = lambda v, lo, hi: f'{v:.4f} [{lo:.4f},{hi:.4f}]'      # noqa: E731
        kp = c(L['K_pl'], L['K_pl_lo'], L['K_pl_hi']) if 'K_pl' in L else 'n/a'
        dc = f"{L['Dc_mean']:.3e}" if L.get('Dc') is not None else 'n/a'
        ds = L.get('K_sys', 0.0)
        print(f"{L['eps_nominal']:7.3f} {L['eps']:8.4f} {L['strain'].get('rg', np.nan):8.4f} "
              f"{c(L['K_net'], L['K_net_lo'], L['K_net_hi']):>22s} {L.get('K_net_rg', np.nan):10.4f} {kp:>22s} {L['dev_rel']:8.1%} {dc:>10s} "
              f"{ds:10.4f}" + (' *' if ds > 0.5 * (L['K_net_hi'] - L['K_net']) else '')
              + (f" {L['G_from_MK']:11.4f}" if 'G_from_MK' in L else ''))
    print(f"  eps_vol = the K denominator ('{cfg.K_STRAIN}');  K_net(Rg) = the same stress increment over the network's own Rg strain")
    if cfg.RELAX_SYS:
        print('  delta_sys = (plate plateau - fitted tail asymptote) / eps_vol;  * = larger than the K_net half-width')
    print_hold_check(cfg, levels)


# ===========================================================================
#  10. D_c OF THE CUBE  (held-cube modes, one D_c for the three axes)
# ===========================================================================
def _cube_design(cfg, Fa, Dc, t, zh=None, decay_only=False):
    """Design matrix of one axis at hold time t: one column per (l, s) with l = 1..DC_N_MODES
    the longitudinal mode sin(2 l pi zeta) and s a distinct transverse rate m^2/L_b^2 + n^2/L_c^2
    (m, n <= DC_N_T); lambda = 4 pi^2 D_c (l^2/L_a^2 + s).  Hold-referenced (u = 0 at the hold
    onset), the time factor is exp(-lambda t) - 1; decay_only gives exp(-lambda t) (the
    transient itself).  The exponential is averaged over the ave/chunk window Fa['W'] like the
    block-averaged snapshots (tri.window_decay, 2026-10-03).  The transverse cos modes average
    over the core column to constants, so every (l, s) pair is one free amplitude."""
    zf = (Fa['zf'] if zh is None else np.asarray(zh, float))[:, None]
    cols = []
    for l in range(1, cfg.DC_N_MODES + 1):
        shape = np.sin(2.0 * np.pi * l * zf)
        for sv in Fa['s_vals']:
            lam = 4.0 * np.pi ** 2 * Dc * (l ** 2 / Fa['L'] ** 2 + sv)
            dec = float(tri.window_decay(lam, t, Fa.get('W', 0.0)))
            cols.append(shape * (dec if decay_only else dec - 1.0))
    return np.hstack(cols)


def fit_Dc_cube(cfg, R, L):
    """The three-dimensional counterpart of triaxial.fit_Dc for the held cube.
    Boundary conditions: every face is drained at the bath pressure -- the feed BC of the
    compression fit (u'' = -q/D_c and equal pore pressure on both faces of every axis, so the
    strain modes are the periodic cos(2 m pi zeta) family) -- and the displacement vanishes at
    the centre (the plates close symmetrically about it).  The dilatation therefore relaxes in
    cos(2 l pi zeta_x) cos(2 m pi zeta_y) cos(2 n pi zeta_z) at 4 pi^2 D_c (l^2/L_x^2 + m^2/L_y^2
    + n^2/L_z^2), and the core-column displacement along axis a in sin(2 l pi zeta_a) times the
    column average of the transverse factors.  One D_c is fitted to the hold-referenced
    u_a/L_a of every axis at once (free amplitudes, linear least squares at each D_c, bounded
    scalar search over D_c, as triaxial.fit_Dc).  Returns dict(Dc, R2, hold_T, n_amp, ax={axis:
    per-axis dict in triaxial's plotting layout}, hold_check) or None."""
    axs = [a for a in cfg.AXES if L.get('disp', {}).get(a) is not None and a in R['ax']]
    F_ax, Lhat = {}, {}
    for a in axs:
        d = L['disp'][a]
        Ra = R['ax'][a]
        z_sup, z_pist = float(d['z_supp_held']), float(d['z_pist_held'])
        span = z_pist - z_sup
        Lg = float(d.get('L_bb', span - 2.0))
        gap = 0.5 * (span - Lg)
        DL_pist, DL_sup = Ra['plate_hi'] - z_pist, z_sup - Ra['plate_lo']
        DL = DL_pist + DL_sup
        if not (Lg > 0 and DL > 0):
            continue
        Lhat[a] = Lg
        F_ax[a] = dict(L=Lg, DL=DL, DL_pist=DL_pist, DL_sup=DL_sup, f_sup=float(DL_sup / DL), gap=gap, z_sup=z_sup,
                       z_pist=z_pist, zeta=(d['z'] - (z_sup + gap)) / Lg, uhat=d['uz'] / Lg, ts=d['ts'],
                       t_lj=(d['ts'] - d['t_hold']) * cfg.dt_lj, Nc=d['Nc'],
                       W=float(np.min(np.diff(d['ts']))) * cfg.dt_lj if (cfg.DC_WINDOW_AVG and len(d['ts']) > 1) else 0.0)
    if not F_ax:
        return None
    for a, Fa in F_ax.items():
        idx = np.where((Fa['Nc'].min(axis=0) > cfg.Ncount_min) & (Fa['zeta'] > 0) & (Fa['zeta'] < 1))[0]
        if len(idx) < 4 + 2 * cfg.DC_TRIM_BINS:
            return None
        if cfg.DC_TRIM_BINS:
            idx = idx[cfg.DC_TRIM_BINS:-cfg.DC_TRIM_BINS]
        Fa['idx'], Fa['zf'] = idx, Fa['zeta'][idx]
        Fa['early'] = np.where(Fa['t_lj'] <= cfg.DC_FRAC_EARLY * Fa['t_lj'][-1])[0]
        others = [b for b in Lhat if b != a]
        others = (others + others)[:2] if others else [a, a]      # a lone axis: its own edge stands in for the other two
        sv = {round(m ** 2 / Lhat[others[0]] ** 2 + n ** 2 / Lhat[others[1]] ** 2, 12)
              for m in range(cfg.DC_N_T + 1) for n in range(cfg.DC_N_T + 1)}
        Fa['s_vals'] = np.array(sorted(sv))
        if np.nanmax(np.abs(Fa['uhat'][:, idx] * Lg)) > 0.5 * Fa['DL']:
            print(f'  WARNING: |u| ~ DL along {a} -- the disp file does not look hold-referenced.')
        if len(Fa['early']) < 2:
            return None

    def solve(Dc):
        res, amps = 0.0, {}
        for a, Fa in F_ax.items():
            X = np.vstack([_cube_design(cfg, Fa, Dc, Fa['t_lj'][i]) for i in Fa['early']])
            y = np.concatenate([Fa['uhat'][i][Fa['idx']] for i in Fa['early']])
            A = np.linalg.lstsq(X, y, rcond=None)[0]
            amps[a] = A
            res += float(np.sum((X @ A - y) ** 2))
        return res, amps

    Dc = float(minimize_scalar(lambda D: solve(D)[0], bounds=cfg.DC_BOUNDS, method='bounded').x)
    ss_res, amps = solve(Dc)
    y_all = np.concatenate([Fa['uhat'][i][Fa['idx']] for Fa in F_ax.values() for i in Fa['early']])
    ss_t = float(np.sum((y_all - y_all.mean()) ** 2))
    R2 = float(1.0 - ss_res / ss_t) if ss_t > 1e-30 else np.nan
    hold_T = float(max(Fa['t_lj'][-1] for Fa in F_ax.values()))
    for a, Fa in F_ax.items():
        A = amps[a]
        Fa.update(A=A, Dc=Dc, R2=R2, hold_T=hold_T, kk=2.0 * np.arange(1, cfg.DC_N_MODES + 1))

        def T(zh, t, Fa=Fa, A=A):
            return _cube_design(cfg, Fa, Dc, t, zh=zh, decay_only=True) @ A
        dlL, fs = Fa['DL'] / Fa['L'], Fa['f_sup']
        Fa['T'] = T
        Fa['u_model'] = lambda zh, t, T=T, dlL=dlL, fs=fs: dlL * (fs - np.asarray(zh, float)) + T(zh, t)
        Fa['u_IC'] = lambda zh, um=Fa['u_model']: um(zh, 0.0)
    return dict(Dc=Dc, R2=R2, hold_T=hold_T, ax=F_ax, n_amp=int(sum(len(A) for A in amps.values())),
                hold_check=tri.hold_adequacy(cfg, max(Lhat.values()), hold_T, Dc))


# ===========================================================================
#  11. NETWORK ANISOTROPY
# ===========================================================================
def _ratio_series(La, ref=False):
    """sigma'_ii(t) / Pi(t) in the interior of one axis view (ref=True: increments from the
    eps_vol = 0 reference over the Pi increment) -> {comp: series}."""
    im = La['interior'] if La['interior'].sum() >= 3 else La['in_mem']
    n = len(La['Pi'])
    Pi = np.array([np.nanmean(La['Pi'][i][im]) for i in range(n)])
    p0 = float(np.mean([La['stress'][c]['ref_net'] for c in COMPONENTS]))
    out = {}
    for c in COMPONENTS:
        S = La['stress'][c]
        sv = np.array([np.nanmean(S['net'][i][im]) for i in range(n)])
        with np.errstate(invalid='ignore', divide='ignore'):
            out[c] = (sv - S['ref_net']) / (Pi - p0) if ref else sv / Pi
    return out


def fig_ratio(cfg, R, L):
    """Network-stress anisotropy over the hold, per profile axis: (a) sigma'_ii / Pi in the gel
    interior (absolute; 1 = isotropic), (b) the increments from the eps_vol = 0 reference,
    Delta sigma'_ii / Delta Pi (the triaxial ratio figure's convention).  Shaded = plateau."""
    axs = [a for a in cfg.AXES if a in L['ax']]
    if not axs:
        print('anisotropy figure skipped (no profiles)')
        return None
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(18, 6), constrained_layout=True)
    fig.suptitle(f"Network-stress anisotropy over the hold   |   {cfg.RUN_ID}   |   $\\varepsilon_{{vol}}={sig(L['eps'])}$",
                 fontsize=14, fontweight='bold')
    ls = {'xx': '-', 'yy': '--', 'zz': ':'}
    for ax, ref in ((axA, False), (axB, True)):
        vals = []
        for a in axs:
            La = L['ax'][a]
            ts = La['stress']['zz']['ts']
            r = _ratio_series(La, ref)
            for c in COMPONENTS:
                ax.plot(ts, r[c], ls[c], color=AX_COLOR[a], lw=1.8, marker='o', ms=3.5, alpha=0.9)
                vals.append(r[c])
        ax.axvspan(L['halt_ts'], float(L['ts'][-1]), color=WONG['green'], alpha=0.10)
        ax.axhline(1.0, color='k', ls=':', lw=1.0, alpha=0.6)
        ax.set_xlabel('time step')
        ax.grid(alpha=0.3)
        h = [Line2D([0], [0], color=AX_COLOR[a], lw=3, label=f'profile along {a}') for a in axs]
        h += [Line2D([0], [0], color='0.3', ls=ls[c], lw=2,
                     label=(r"$\Delta\sigma'_{%s}/\Delta\mathit{\Pi}$" if ref else r"$\sigma'_{%s}/\mathit{\Pi}$") % c) for c in COMPONENTS]
        robust_ylim(ax, vals, pad=0.3, qlo=2, qhi=98, include_zero=False)
        smart_legend(ax, handles=h, fontsize=11)
    axA.set_ylabel(r"$\langle\sigma'_{ii}\rangle_{\rm int}/\langle\mathit{\Pi}\rangle_{\rm int}$")
    axA.set_title('(a) absolute (1 = isotropic network stress)', fontsize=13)
    axB.set_ylabel(r"$\Delta\langle\sigma'_{ii}\rangle_{\rm int}/\Delta\langle\mathit{\Pi}\rangle_{\rm int}$")
    axB.set_title(r'(b) increments from $\varepsilon_{vol}=0$ (the response to the plates)', fontsize=13)
    return _save(fig, cfg, 'network_stress_ratio', L['lvl'])


def fig_ratio_sweep(cfg, R, levels):
    """Delta sigma'_ii / Delta Pi over the hold along cfg.AXES[0], every level (colour = level)."""
    a0 = cfg.AXES[0]
    Ls = [L for L in levels if a0 in L['ax']]
    if not Ls:
        print('anisotropy sweep figure skipped (no profiles)')
        return None
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    ls = {'xx': '-', 'yy': '--', 'zz': ':'}
    vals = []
    for L in Ls:
        i = levels.index(L)
        r = _ratio_series(L['ax'][a0], True)
        ts = L['ax'][a0]['stress']['zz']['ts']
        for c in COMPONENTS:
            ax.plot(ts, r[c], ls[c], color=level_color(i), lw=1.8, marker='o', ms=3.5, alpha=0.9)
            vals.append(r[c])
    ax.axhline(1.0, color='k', ls=':', lw=1.0, alpha=0.6)
    robust_ylim(ax, vals, pad=0.3, qlo=2, qhi=98, include_zero=False)
    h = level_handles(levels) + [Line2D([0], [0], color='0.3', ls=ls[c], lw=2, label=r"$\Delta\sigma'_{%s}/\Delta\mathit{\Pi}$" % c)
                                 for c in COMPONENTS]
    ax.set_xlabel('time step')
    ax.set_ylabel(r"$\Delta\langle\sigma'_{ii}\rangle_{\rm int}/\Delta\langle\mathit{\Pi}\rangle_{\rm int}$")
    ax.set_title(f'Network-stress anisotropy (increments from $\\varepsilon_{{vol}}=0$) along {a0}, all levels', fontsize=15)
    ax.grid(alpha=0.3)
    smart_legend(ax, handles=h, fontsize=11)
    return _save(fig, cfg, 'sweep_network_stress_ratio')


def whole_run_file(cfg, name):
    """The stem-tagged whole-run file <name>_<DATANAME>_<INTERACTION>_<nsteps>.dat (its tag is
    the batch NSTEPS, not a hold length): the newest match in DATA_DIR, or None."""
    import re as _re
    pat = _re.compile(rf'^{_re.escape(name)}_{_re.escape(cfg.DATANAME)}_{_re.escape(cfg.INTERACTION)}_(\d+)\.dat$')
    hits = [f for f in cfg.DATA_DIR.glob(f'{name}_*.dat') if pat.match(f.name) and f.stat().st_size > 0]
    return max(hits, key=lambda f: f.stat().st_mtime) if hits else None


def read_stress_aniso(path):
    """stress_aniso file (slab_with_support, the triaxial decks, compress_slab): TimeStep
    sig_p_xx sig_p_yy sig_p_zz P_xx P_yy P_zz -> dict of columns, or None."""
    if path is None or not Path(path).exists():
        return None
    A = tri.load2c(path, 4)
    if A is None:
        return None
    out = dict(step=A[:, 0], xx=A[:, 1], yy=A[:, 2], zz=A[:, 3])
    with np.errstate(invalid='ignore', divide='ignore'):
        out['xx/zz'], out['yy/zz'] = out['xx'] / out['zz'], out['yy'] / out['zz']
    return out


def fig_anisotropy_run(cfg, R, levels=None, swell_file=None):
    """Whole-run polymer partial-stress anisotropy sig_p,xx/sig_p,zz and sig_p,yy/sig_p,zz from
    the stress_aniso file the deck writes from step 0 (settle, seating, reference, every level);
    the levels' holds are shaded.  `swell_file` = the stress_aniso file of the slab_with_support
    run that produced the input (free swelling of the isolated gel): drawn in a first panel so
    the swelling anisotropy and the compression anisotropy stand side by side."""
    run = read_stress_aniso(whole_run_file(cfg, 'stress_aniso'))
    sw = read_stress_aniso(swell_file) if swell_file else None
    if run is None and sw is None:
        print('anisotropy-run figure skipped (no stress_aniso file; sync it or set swell_file)')
        return None
    panels = ([('free swelling (slab_with_support): ' + Path(swell_file).name, sw, False)] if sw is not None else []) \
        + ([('compress_slab run (shaded: the holds)', run, True)] if run is not None else [])
    fig, axes = plt.subplots(1, len(panels), figsize=(10 * len(panels), 6), constrained_layout=True, squeeze=False)
    for k, (ax, (title, D, mark)) in enumerate(zip(axes[0], panels)):
        ax.plot(D['step'], D['xx/zz'], '-', marker='o', ms=4, color=WONG['blue'], lw=2, label=r'$\sigma_{p,xx}/\sigma_{p,zz}$')
        ax.plot(D['step'], D['yy/zz'], '--', marker='s', ms=4, color=WONG['vermillion'], lw=2, label=r'$\sigma_{p,yy}/\sigma_{p,zz}$')
        if mark:
            for i, L in enumerate(levels or []):
                ax.axvspan(float(L['ts'][0]), float(L['ts'][-1]), color=level_color(i), alpha=0.10)
                ax.text(float(L['ts'][0]), 1.0, f' $\\varepsilon_{{vol}}={L["lvl"]}$', fontsize=10, color=level_color(i), va='bottom')
        ax.axhline(1.0, color='k', ls=':', lw=1.0, alpha=0.6)
        ax.set_xlabel('time step')
        ax.set_ylabel(r'$\sigma_{p,ii}/\sigma_{p,zz}$  (whole-polymer partial stress)')
        ax.set_title(f'({"ab"[k]}) {title}', fontsize=11)
        ax.grid(alpha=0.3)
        robust_ylim(ax, [D['xx/zz'], D['yy/zz']], pad=0.3, qlo=1, qhi=99, include_zero=False)
        smart_legend(ax, fontsize=12)
    fig.suptitle(f'Polymer partial-stress anisotropy over the whole run  |  {cfg.sim_name}', fontsize=13, fontweight='bold')
    return _save(fig, cfg, 'anisotropy_whole_run')


# ===========================================================================
#  12. SOLVENT VOLUME FRACTIONS IN THE CORE COLUMNS  (Voronoi + calibration)
# ===========================================================================
def _column_masks(cfg, R, L, p):
    """({axis: bool mask of the tessellated points inside the axis' core column}, core bounds)
    for the plate planes of the state (R, or a level L)."""
    planes = (L or R)['plane']
    core = {a: (planes[a + 'lo'] + cfg.core_margin, planes[a + 'hi'] - cfg.core_margin) for a in 'xyz'}
    inside = {a: (p[:, AXI[a]] > core[a][0]) & (p[:, AXI[a]] < core[a][1]) for a in 'xyz'}
    return {a: np.logical_and.reduce([inside[b] for b in 'xyz' if b != a]) for a in 'xyz'}, core


def _phi_columns(cfg, R, L, box, typ, xyz):
    """One tessellation of a frame -> {axis: (phi_bin, phi_mobile)} on the axis' bin grid: the
    solvent Voronoi volume in the core column per bin over (a) the column bin volume and (b)
    the summed cell volume in the bin (exact saturation, the choice for the calibrated phi)."""
    p, t_k, vol, origin, Lbox = volfrac._tessellate(box, typ, xyz, mobile_only=cfg.VOR_MOBILE_ONLY)
    masks, core = _column_masks(cfg, R, L, p)
    out = {}
    for a in cfg.AXES:
        if a not in R['ax']:
            continue
        z = R['ax'][a]['z']
        bw = cfg.binWidth
        A_col = float(np.prod([core[b][1] - core[b][0] for b in 'xyz' if b != a]))
        m = masks[a]
        bi = np.clip(np.floor((p[m, AXI[a]] - (z[0] - 0.5 * bw)) / bw).astype(int), 0, len(z) - 1)
        s = t_k[m] == volfrac.SOLVENT_TYPE
        V_s = np.bincount(bi[s], weights=vol[m][s], minlength=len(z))[:len(z)]
        V_t = np.bincount(bi, weights=vol[m], minlength=len(z))[:len(z)]
        with np.errstate(invalid='ignore', divide='ignore'):
            out[a] = (V_s / (A_col * bw), np.where(V_t > 0, V_s / V_t, np.nan))
    return out


def _vf_columns(cfg, R, L, traj, ts_want, label=''):
    """Stream `traj` once for the frames ts_want -> {axis: dict(ts, phi_vor, phi_vor_mobile)}
    (stacks [n_frames, n_bins]); cached on disk in DATA_DIR/vor_cache."""
    ts_want = [int(t) for t in np.atleast_1d(ts_want)]
    planes = (L or R)['plane']
    key = repr((Path(traj).name, tuple(ts_want), cfg.VOR_MOBILE_ONLY, cfg.binWidth, tuple(cfg.AXES), cfg.core_margin,
                tuple(round(planes[f], 4) for f in PLATES)))
    cdir = Path(cfg.DATA_DIR) / 'vor_cache'
    cfile = cdir / f'{Path(traj).stem}__col_{hashlib.md5(key.encode()).hexdigest()[:12]}.npz'
    if cfile.exists():
        try:
            d = np.load(cfile, allow_pickle=False)
            out = {a: dict(ts=d['ts'], phi_vor=d[f'{a}_bin'], phi_vor_mobile=d[f'{a}_mob']) for a in cfg.AXES if f'{a}_bin' in d}
            print(f'    {label}disk-cached ({cfile.parent.name}/{cfile.name}): {len(d["ts"])} frame(s), no re-read')
            return out
        except Exception as e:
            print(f'    {label}vor_cache unreadable ({type(e).__name__}) -> recomputing')
    print(f'    {label}streaming {Path(traj).name} for {len(ts_want)} frame(s): Voronoi, three columns per frame ...')
    fr = volfrac.stream_traj_frames(traj, ts_want)
    ok, per = [], {a: ([], []) for a in cfg.AXES}
    for t in ts_want:
        if t not in fr:
            print(f'    NOTE: no traj frame at {t} -- dropped')
            continue
        try:
            res = _phi_columns(cfg, R, L, *fr[t])
        except Exception as e:
            print(f'      ts {t}: tessellation failed ({e}) -- frame dropped')
            continue
        for a in res:
            per[a][0].append(res[a][0])
            per[a][1].append(res[a][1])
        ok.append(t)
        print(f'      ts {t}: phi^vor in the gel interior ~ ' + '  '.join(
            f'{a}: {np.nanmean(res[a][0][R["ax"][a]["interior"]]):.3f}' for a in res))
    del fr
    if not ok:
        return {}
    out = {a: dict(ts=np.array(ok, float), phi_vor=np.array(v[0]), phi_vor_mobile=np.array(v[1])) for a, v in per.items() if v[0]}
    out = {a: v for a, v in out.items() if np.isfinite(v['phi_vor']).any()}
    if not out:
        print(f'    {label}no usable frame (tessellation empty in the columns) -> Voronoi skipped')
        return {}
    try:
        cdir.mkdir(parents=True, exist_ok=True)
        np.savez(cfile, ts=np.array(ok, float), **{f'{a}_bin': out[a]['phi_vor'] for a in out},
                 **{f'{a}_mob': out[a]['phi_vor_mobile'] for a in out})
    except Exception as e:
        print(f'    {label}(vor_cache not written: {type(e).__name__}: {e})')
    return out


def _p_local_axis(cfg, R, D, a):
    """P handed to lambda(phi_p, P) per bin of axis a: 'const' -> P_CAL; 'pore' -> the bath
    baseline of that axis (uniform: one connected bath); 'thermo' -> the local P_th."""
    Ra = R['ax'][a]
    z = Ra['z']
    Da = D['ax'][a]
    Pmin, Pmax = tri._calib_range(R)
    if cfg.P_CAL_MODE == 'const':
        return lambda t: np.clip(np.full(len(z), float(cfg.P_CAL)), Pmin, Pmax)
    ts = np.asarray(Da['stress']['zz']['ts'], float)
    if cfg.P_CAL_MODE == 'thermo':
        Pt = tri._tr3(Da['stress'], 't')
        return lambda t: np.clip(np.asarray(Pt[int(np.argmin(np.abs(ts - t)))], float), Pmin, Pmax)
    pore = np.asarray(Da['stress']['zz']['pore'], float)
    if pore.ndim == 0:
        return lambda t: np.clip(np.full(len(z), float(pore)), Pmin, Pmax)
    return lambda t: np.clip(np.full(len(z), float(np.interp(t, ts, pore))), Pmin, Pmax)


def _apply_vf(cfg, R, D, vf):
    """attach phi_vor / phi_cal (mean, lo, hi) to every axis view of D from the frame stacks vf."""
    for a, v in vf.items():
        Da = D['ax'][a]
        Da['phi_vor'] = mean_ci(v['phi_vor'], cfg.ci_level)
        Da['phi_vor_frames'] = v['ts']
        if R.get('CALIB') is not None:
            pfn = _p_local_axis(cfg, R, D, a)
            Pl = np.array([pfn(t) for t in v['ts']])
            Da['phi_cal'] = mean_ci(volfrac.phi_calibrated(v['phi_vor_mobile'], Pl, R['CALIB']), cfg.ci_level)


def add_volume_fractions(cfg, R, levels=()):
    """Voronoi + lambda-calibrated solvent volume fractions in the core columns along every
    axis, for the reference (REF_VOR_FRAMES frames of traj_ref) and each level (VOR_MAX_FRAMES
    frames inside the plateau of traj_stress, plus VOR_EVO_FRAMES frames over the whole hold
    for the evolution figure).  The mass-fraction profiles were set by load_reference /
    load_level.  P per bin follows cfg.P_CAL_MODE ('const' = P_CAL, the drained-equilibrium
    choice; 'pore' = the axis' bath baseline; 'thermo' = the local P_th)."""
    tri._load_calib(R)
    if not cfg.VOR_ENABLE:
        print('VOR_ENABLE=False -> Voronoi / calibrated phi skipped')
        return
    tr = cfg.traj('traj_ref')
    if tr.exists():
        rts = R['ax'][cfg.AXES[0]]['stress']['zz']['ts']
        want, _ = tri.subsample(rts, rts[:, None], cfg.REF_VOR_FRAMES)
        vf = _vf_columns(cfg, R, None, tr, want, 'reference: ')
        if vf:
            _apply_vf(cfg, R, R, vf)
    else:
        print(f'NOTE: {tr.name} missing -> reference Voronoi skipped')
    for L in levels:
        if L is None:
            continue
        tp = cfg.traj('traj_stress', L['lvl'])
        if not tp.exists():
            print(f'NOTE: {tp.name} missing -> Voronoi skipped for level {L["lvl"]}')
            continue
        ts = L['ts']
        plat = ts[ts >= L['halt_ts']]
        plat = plat if len(plat) else ts[-1:]
        want_p, _ = tri.subsample(plat, plat[:, None], cfg.VOR_MAX_FRAMES)
        want_e, _ = tri.subsample(ts, ts[:, None], cfg.VOR_EVO_FRAMES)
        vf = _vf_columns(cfg, R, L, tp, np.unique(np.concatenate([want_p, want_e])), f'level {L["lvl"]}: ')
        if not vf:
            continue
        a0 = next(iter(vf))
        sel = np.isin(vf[a0]['ts'], want_p)
        if not sel.any():
            sel[:] = True
        _apply_vf(cfg, R, L, {a: dict(ts=v['ts'][sel], phi_vor=v['phi_vor'][sel], phi_vor_mobile=v['phi_vor_mobile'][sel])
                              for a, v in vf.items()})
        for a, v in vf.items():                      # every frame of the hold: fig_volfrac_evolution
            cal = None
            if R.get('CALIB') is not None:
                pfn = _p_local_axis(cfg, R, L, a)
                cal = volfrac.phi_calibrated(v['phi_vor_mobile'], np.array([pfn(t) for t in v['ts']]), R['CALIB'])
            L['ax'][a]['vf_evo'] = dict(ts=v['ts'], phi_vor=v['phi_vor'], phi_cal=cal)
    print(f'in-gel solvent volume fractions (interior; P_CAL_MODE = {cfg.P_CAL_MODE}):')
    for a in cfg.AXES:
        if a in R['ax']:
            tri._vf_line(f'reference, along {a}', R['ax'][a], R['ax'][a]['interior'])
    for L in levels:
        if L is None:
            continue
        for a in cfg.AXES:
            if a in L['ax']:
                tri._vf_line(f'eps={L["eps"]:.3f}, along {a}', L['ax'][a], L['ax'][a]['interior'])


def fig_volfrac_evolution(cfg, R, L, key='phi_cal'):
    """phi_s(a, t) over the hold along x, y and z (cividis, final bold; reference dashed):
    the lambda-calibrated Voronoi fraction (key='phi_cal'; 'phi_vor' raw, 'phi_mf' the mass
    fraction from every density snapshot).  Needs add_volume_fractions (not for 'phi_mf')."""
    axs = [a for a in cfg.AXES if a in L['ax']]
    lab = {'phi_cal': r'$\lambda$-calibrated Voronoi $\phi_s^{\rm cal}$', 'phi_vor': r'Voronoi $\phi_s^{\rm vor}$',
           'phi_mf': r'mass fraction $\phi_s^{\rm mf}=\rho_s/\rho_{\rm bath}$'}[key]
    fig, axes = plt.subplots(1, len(axs), figsize=(8.5 * len(axs), 6.5), constrained_layout=True, squeeze=False)
    drawn = False
    for i, (ax, a) in enumerate(zip(axes[0], axs)):
        La, Ra = L['ax'][a], R['ax'][a]
        Rv, Lv = view(cfg, R, L, a)
        if key == 'phi_mf':
            dens = _density(cfg, a, lvl=L['lvl'])
            if dens is None or not np.isfinite(Ra['rho_s0']):
                continue
            ts, ev = tri.post_halt(cfg, Lv, dens['ts'], dens['n'] / Ra['rho_s0'])
        else:
            ve = La.get('vf_evo')
            if ve is None or ve.get(key) is None:
                continue
            ts, ev = np.asarray(ve['ts'], float), np.asarray(ve[key], float)
        drawn = True
        tri.plot_evolution(ax, cfg, Rv, Lv, Rv['z'], ts, ev, r'$\phi_s$', f'({"abc"[i]}) along {a}: ' + lab,
                           ref=Ra.get(key), annotate=False, legend=False)
        ax.set_xlabel(f'${a}/L$')
        ax.axhline(1.0, color='k', ls=':', lw=1.2, alpha=0.6)
        ax.set_ylim(0, 1.15)
        h, lb = ax.get_legend_handles_labels()
        if h:
            smart_legend(ax, handles=h, labels=lb, fontsize=11)
        if La.get(key) is not None and Ra.get(key) is not None:
            d = np.nanmean(La[key][0][La['interior']]) - np.nanmean(Ra[key][0][Ra['interior']])
            annotate_box(ax, r'$\Delta\phi_s$ (plateau $-$ reference), in gel: ' + f'{d:+.2g}', loc='lower right', fontsize=12)
    if not drawn:
        plt.close(fig)
        print(f'volume-fraction evolution skipped ({key}: run add_volume_fractions first, or no calibration artifact)')
        return None
    fig.suptitle(f'Solvent volume fraction evolution over the hold, $\\varepsilon_{{vol}}={L["lvl"]}$  |  {cfg.sim_name}',
                 fontsize=13, fontweight='bold')
    return _save(fig, cfg, f'volfrac_evolution_{key}', L['lvl'])
