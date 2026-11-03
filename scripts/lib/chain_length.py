#!/usr/bin/env python3
"""chain_length.py -- the beads-per-strand comparison (scripts/chain_length_comparison.ipynb, 2026-11-03).

Two kinds of input, both keyed by N = beads per strand:
  * the EQUILIBRATION runs of simulations/chain_length_testing (slab_with_support.lmp with the Rg-plateau
    halt): the volume_data files (gel_volume_rg, gel_dimensions, box_dimensions) and the log's RG_BLOCK /
    RG_PLATEAU_CONVERGED lines -> Rg(t), thickness and in-plane box per N, where the halt fired;
  * the COMPRESSION runs (triaxial_compression_two_pist on the converted final_config of each N), each loaded
    by lib/triaxial.py exactly as the compression notebooks load it (compare.load_deck) -> M, G, phi_s per N.

The figures: the equilibration traces overlaid (one line per N), and the observables against N on log axes with
a fitted power-law guide each.  Nothing here recomputes what triaxial.py computes; it only collects and draws.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# the sweep (scripts/lib/build_beads_sweep.py + simulations/chain_length_testing/chain_length_testing_sweep.batch)
# ---------------------------------------------------------------------------
A_REF, N_REF, NU = 4.82, 5, 0.59


def a_eq(n):
    """The lattice constant the N-bead random slab was built at (= the mesh the gel is expected to swell to)."""
    return A_REF * ((n - 1) / (N_REF - 1)) ** NU


SWEEP = {   # N: (equilibration DATANAME, step budget = the <totsteps> tag of its files, in-plane cells, z cells)
    5:  ('slab_support_periodic_5beads_tall_rho04_random',  14_000_000, 9, 23),
    10: ('slab_support_periodic_10beads_tall_rho04_random', 24_000_000, 6, 14),
    15: ('slab_support_periodic_15beads_tall_rho04_random', 30_000_000, 4, 11),
    20: ('slab_support_periodic_20beads_tall_rho04_random', 35_000_000, 4, 9),
    30: ('slab_support_periodic_30beads_tall_rho04_random', 45_000_000, 3, 7),
}
COLORS = {5: '#444444', 10: '#1f77b4', 15: '#2ca02c', 20: '#ff7f0e', 30: '#d62728'}


def eq_dataname(n):
    return SWEEP[n][0]


def eq_final_config(n, interaction='1.0_1.0'):
    """Stem of the equilibrated final configuration (the two-piston converter's input)."""
    return f'final_config_{eq_dataname(n)}_{interaction}_{SWEEP[n][1]}'


# ---------------------------------------------------------------------------
# equilibration runs
# ---------------------------------------------------------------------------
def _read_fix_print(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        try:
            rows.append([float(x) for x in line.split()])
        except ValueError:
            continue
    return np.array(rows) if rows else np.zeros((0, 2))


def find_eq_dir(n, roots, interaction='1.0_1.0'):
    """The newest run directory of N's equilibration under any of `roots` (each root is searched for
    <folder>_<DATANAME>_<interaction>_<timestamp> and for the bare files of a synced copy)."""
    dn = eq_dataname(n)
    hits = []
    for root in roots:
        root = Path(root).expanduser()
        if not root.exists():
            continue
        hits += [d for d in root.glob(f'*_{dn}_{interaction}_*') if d.is_dir()]
        if (root / f'gel_volume_rg_{dn}_{interaction}_{SWEEP[n][1]}.dat').exists():
            hits.append(root)
    return sorted(hits, key=lambda d: d.name)[-1] if hits else None


def load_equilibration(n, roots=('~/Documents/lammps_runs/chain_length_testing',
                                '~/Documents/lammps_runs/slab_with_support',
                                '../../flow_data_local/chain_length'),
                       interaction='1.0_1.0', dt=0.005):
    """Rg(t), gel thickness and box traces of N's equilibration run (+ its continuations), the Rg-halt
    bookkeeping from the log.  Returns None (with a note) when no run directory is found."""
    d = find_eq_dir(n, roots, interaction)
    if d is None:
        print(f'N = {n:2d}: no equilibration run found for {eq_dataname(n)} under {list(roots)}')
        return None
    tag = f'{eq_dataname(n)}_{interaction}_{SWEEP[n][1]}'
    vd = d / 'output_files' / 'volume_data' if (d / 'output_files').exists() else d
    parts = [vd] + sorted(p / 'output_files' / 'volume_data' for p in d.glob('continuation_*'))
    out = dict(N=n, dir=d, tag=tag, dt=dt)
    for key, fname, cols in (('rg', 'gel_volume_rg', (0, 1)), ('gel', 'gel_dimensions', (0, 1, 2, 3)),
                             ('box', 'box_dimensions', (0, 1, 2, 3))):
        arrs = [_read_fix_print(p / f'{fname}_{tag}.dat') for p in parts if (p / f'{fname}_{tag}.dat').exists()]
        arrs = [a for a in arrs if len(a)]
        if not arrs:
            out[key] = None
            continue
        a = np.vstack(arrs)
        a = a[np.argsort(a[:, 0], kind='stable')]
        out[key] = a
    if out['rg'] is not None:
        out['rg'] = np.column_stack([out['rg'][:, 0], out['rg'][:, 1] ** (1 / 3)])   # the file holds Rg^3
    # halt bookkeeping from the log(s)
    blocks, det, conv = [], None, None
    for log in [d / 'log.lammps'] + sorted(d.glob('continuation_*/log.lammps')):
        if not log.exists():
            continue
        for line in log.read_text(errors='replace').splitlines():
            m = re.match(r'>>> RG_BLOCK step (\d+) Rg_block ([\d.]+) Rg_prev ([\d.]+) rel_change ([\d.e+-]+) consecutive_under_tol_after (\d+)', line)
            if m:
                blocks.append(tuple(float(x) for x in m.groups()))
            m = re.match(r'>>> RG_PLATEAU_DETECTED at step (\d+) after (\d+) production steps', line)
            if m:
                det = (int(m.group(1)), int(m.group(2)))
            m = re.match(r'>>> RG_PLATEAU_CONVERGED at step (\d+) after (\d+) production steps', line)
            if m:
                conv = (int(m.group(1)), int(m.group(2)))
    out['blocks'] = np.array(blocks) if blocks else None
    out['detected'] = det
    out['converged'] = conv
    g = out['gel']
    out['thickness'] = float(np.mean(g[-10:, 3])) if g is not None and len(g) else np.nan
    b = out['box']
    out['lx'] = float(np.mean(b[-10:, 1])) if b is not None and len(b) else np.nan
    out['Rg'] = float(np.mean(out['rg'][-10:, 1])) if out['rg'] is not None and len(out['rg']) else np.nan
    print(f"N = {n:2d}: {d.name}  Rg {out['Rg']:.2f}  thickness {out['thickness']:.1f}  lx {out['lx']:.2f}  "
          f"halt: {'converged at step %d' % conv[0] if conv else 'not converged (budget or still running)'}  "
          f"({0 if out['blocks'] is None else len(out['blocks'])} blocks logged)")
    return out


def fig_equilibration(EQ, out=None):
    """(a) Rg / Rg(first sample), (b) gel thickness, (c) in-plane box lx vs production steps, one line per N; the
    block means of the Rg-halt as dots on (a) with a star where the halt fired."""
    EQ = [e for e in EQ if e is not None]
    if not EQ:
        print('no equilibration run loaded'); return None
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
    for e in EQ:
        n, c = e['N'], COLORS.get(e['N'], None)
        if e['rg'] is not None:
            s, r = e['rg'][:, 0], e['rg'][:, 1]
            axs[0].plot(s / 1e6, r / r[0], color=c, lw=1.2, label=f'N = {n}')
            if e['blocks'] is not None:
                axs[0].plot(e['blocks'][:, 0] / 1e6, e['blocks'][:, 1] / r[0], 'o', color=c, ms=3)
            if e['detected']:
                axs[0].plot(e['detected'][0] / 1e6, np.interp(e['detected'][0], s, r) / r[0], 'D', color=c, ms=6, mfc='none')
            if e['converged']:
                axs[0].plot(e['converged'][0] / 1e6, np.interp(e['converged'][0], s, r) / r[0], '*', color=c, ms=12)
        if e['gel'] is not None:
            axs[1].plot(e['gel'][:, 0] / 1e6, e['gel'][:, 3], color=c, lw=1.2, label=f'N = {n}')
        if e['box'] is not None:
            axs[2].plot(e['box'][:, 0] / 1e6, e['box'][:, 1], color=c, lw=1.2, label=f'N = {n}')
    axs[0].set_ylabel('$R_g$ / $R_g$(first sample)'); axs[0].set_title('gel radius of gyration (dots: block means, diamond: plateau detected, star: halt)')
    axs[1].set_ylabel('gel thickness  $L_z$ (bounding box, σ)'); axs[1].set_title('thickness')
    axs[2].set_ylabel('box  $l_x$ (σ)'); axs[2].set_title('in-plane box (aniso NPH)')
    for ax in axs:
        ax.set_xlabel('production steps / 10$^6$'); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle('Equilibration of the periodic random slab vs beads per strand', fontweight='bold')
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150, bbox_inches='tight'); print('saved', out)
    plt.show()
    return fig


# ---------------------------------------------------------------------------
# compression runs -> observables vs N
# ---------------------------------------------------------------------------
def _interior_mean(R, prof):
    if prof is None:
        return np.nan
    m = np.asarray(prof[0] if isinstance(prof, tuple) else prof, float)
    idx = R.get('interior')
    try:
        return float(np.nanmean(m[idx])) if idx is not None else float(np.nanmean(m))
    except Exception:
        return float(np.nanmean(m))


def collect_observables(DECKS, EQ=None):
    """One row per N from the loaded compression decks (compare.load_deck output, the first = lowest level):
    M_net, M_pist (with CIs), G_xx, G_yy, the reference and compressed interior solvent fractions phi_s, the
    level's strain; plus the equilibration geometry (Rg, thickness, lx, a_eq) when EQ is given."""
    eq = {e['N']: e for e in (EQ or []) if e is not None}
    rows = []
    for n, D in sorted(DECKS.items()):
        if D is None or not D.get('levels'):
            print(f'N = {n}: no compression level loaded'); continue
        L, R = D['levels'][0], D['R']
        G = L.get('G', {})
        row = dict(N=n, eps=L['eps'], eps_M=L.get('eps_M', np.nan),
                   M_net=L['M_net'], M_net_lo=L['M_net_lo'], M_net_hi=L['M_net_hi'],
                   M_pist=L.get('M_pist', np.nan), M_pist_lo=L.get('M_pist_lo', np.nan), M_pist_hi=L.get('M_pist_hi', np.nan),
                   G_xx=G.get('xx', {}).get('G', np.nan), G_xx_lo=G.get('xx', {}).get('lo', np.nan), G_xx_hi=G.get('xx', {}).get('hi', np.nan),
                   G_yy=G.get('yy', {}).get('G', np.nan),
                   phi_s_ref=_interior_mean(R, R.get('phi_cal')), phi_s=_interior_mean(R, L.get('phi_cal')),
                   a_eq=a_eq(n))
        e = eq.get(n)
        row.update(Rg=e['Rg'] if e else np.nan, thickness=e['thickness'] if e else np.nan, lx=e['lx'] if e else np.nan)
        rows.append(row)
    return rows


def print_table(rows):
    print(f"{'N':>3} {'a_eq':>6} {'eps':>6} {'M_net [lo,hi]':>24} {'M_pist':>8} {'G_xx [lo,hi]':>24} {'G_yy':>7} "
          f"{'phi_s ref':>9} {'phi_s':>7} {'G/M':>6} | {'Rg':>6} {'L_z':>6} {'l_x':>6}")
    for r in rows:
        print(f"{r['N']:>3d} {r['a_eq']:>6.2f} {r['eps']:>6.3f} {r['M_net']:>8.3f} [{r['M_net_lo']:.3f},{r['M_net_hi']:.3f}] "
              f"{r['M_pist']:>8.3f} {r['G_xx']:>8.3f} [{r['G_xx_lo']:.3f},{r['G_xx_hi']:.3f}] {r['G_yy']:>7.3f} "
              f"{r['phi_s_ref']:>9.3f} {r['phi_s']:>7.3f} {r['G_xx'] / r['M_net'] if r['M_net'] else np.nan:>6.3f} | "
              f"{r['Rg']:>6.2f} {r['thickness']:>6.1f} {r['lx']:>6.2f}")
    print('  M, G from the lowest compression level of each N (lib/triaxial.py: M_net = network stress increment / eps_M, '
          'G from the lateral network stress); phi_s = calibrated interior solvent volume fraction (reference = uncompressed)')


def _powerlaw(ax, x, y, color='k', label='fit'):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if ok.sum() < 2:
        return None
    p = np.polyfit(np.log(x[ok]), np.log(y[ok]), 1)
    xx = np.linspace(x[ok].min(), x[ok].max(), 50)
    ax.plot(xx, np.exp(p[1]) * xx ** p[0], '--', color=color, lw=1, label=f'{label}: ∝ N$^{{{p[0]:.2f}}}$')
    return p[0]


def fig_observables_vs_N(rows, out=None):
    """M_net / M_pist, G_xx / G_yy, phi_s (reference and compressed) and the polymer fraction 1 - phi_s, each vs N on
    log-log axes with a fitted power law.  Flory-Rehner guides: the modulus of an ideal network at fixed volume goes as
    the strand density ~ 1/N; in good-solvent swelling equilibrium (c* theorem) phi_p ~ N^(-4/5) and G ~ N^(-9/5)."""
    rows = [r for r in rows if np.isfinite(r.get('M_net', np.nan))]
    if not rows:
        print('nothing to draw'); return None
    N = np.array([r['N'] for r in rows], float)
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.2))
    ax = axs[0]
    ax.errorbar(N, [r['M_net'] for r in rows], yerr=[[r['M_net'] - r['M_net_lo'] for r in rows], [r['M_net_hi'] - r['M_net'] for r in rows]],
                fmt='o', color='#1f77b4', label='$M_{net}$', capsize=3)
    ax.errorbar(N, [r['M_pist'] for r in rows], yerr=[[r['M_pist'] - r['M_pist_lo'] for r in rows], [r['M_pist_hi'] - r['M_pist'] for r in rows]],
                fmt='s', color='#ff7f0e', label='$M_{pist}$', capsize=3, alpha=0.8)
    _powerlaw(ax, N, [r['M_net'] for r in rows], '#1f77b4', '$M_{net}$')
    ax.set_ylabel('longitudinal modulus  M  (ε/σ³)'); ax.set_title('M vs N')
    ax = axs[1]
    ax.errorbar(N, [r['G_xx'] for r in rows], yerr=[[r['G_xx'] - r['G_xx_lo'] for r in rows], [r['G_xx_hi'] - r['G_xx'] for r in rows]],
                fmt='o', color='#2ca02c', label='$G_{xx}$', capsize=3)
    ax.plot(N, [r['G_yy'] for r in rows], 'v', color='#2ca02c', alpha=0.6, label='$G_{yy}$')
    _powerlaw(ax, N, [r['G_xx'] for r in rows], '#2ca02c', '$G_{xx}$')
    g = N ** (-9 / 5); g *= np.nanmax([r['G_xx'] for r in rows]) / g[0] if np.isfinite(rows[0]['G_xx']) else 1
    ax.plot(N, g, ':', color='gray', label='c* guide ∝ N$^{-9/5}$')
    ax.set_ylabel('shear modulus  G  (ε/σ³)'); ax.set_title('G vs N')
    ax = axs[2]
    ax.plot(N, [r['phi_s_ref'] for r in rows], 'o-', color='#9467bd', label='reference (ε = 0)')
    ax.plot(N, [r['phi_s'] for r in rows], 's--', color='#9467bd', alpha=0.6, label='compressed level')
    ax.set_ylabel('interior solvent fraction  φ$_s$'); ax.set_title('φ$_s$ vs N'); ax.set_yscale('linear')
    ax = axs[3]
    pp = 1 - np.array([r['phi_s_ref'] for r in rows])
    ax.plot(N, pp, 'o', color='#8c564b', label='1 − φ$_s$ (reference)')
    _powerlaw(ax, N, pp, '#8c564b', 'polymer fraction')
    pg = N ** (-4 / 5); pg *= pp[0] / pg[0] if np.isfinite(pp[0]) and pp[0] > 0 else 1
    ax.plot(N, pg, ':', color='gray', label='c* guide ∝ N$^{-4/5}$')
    ax.set_ylabel('polymer fraction  1 − φ$_s$'); ax.set_title('swelling vs N')
    for ax in axs:
        ax.set_xscale('log'); ax.set_xlabel('beads per strand  N'); ax.grid(alpha=0.3, which='both'); ax.legend(fontsize=8)
        ax.set_xticks(N); ax.set_xticklabels([f'{int(v)}' for v in N])
    for ax in (axs[0], axs[1], axs[3]):
        ax.set_yscale('log')
    fig.suptitle('Compression observables vs beads per strand (lowest level of each N)', fontweight='bold')
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150, bbox_inches='tight'); print('saved', out)
    plt.show()
    return fig
