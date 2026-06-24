#!/usr/bin/env python3
"""Cluster-side quick-look for the shear_slab stress-strain sweep.

Reads the per-strain output files written by shear_slab.lmp (tagged _g<strain>)
and produces two PNGs in <folder>/output_plots/:

  1. <stem>_stress_strain_curve.png
       sigma_zx vs gamma for every strain in the sweep, with block-SEM error
       bars (scatter across the ~num_stress_curves nfreq blocks of each
       production hold) and a through-origin weighted fit for G.

  2. <stem>_sigma_zx_profiles.png  (only if stress_profile_x files are present)
       polymer sigma_zx(x) across the gap, one curve per strain.

This complements (and matches) shear_analysis.ipynb Steps 12-13; it is just an
automatic first look so you don't have to open the notebook to sanity-check a run.

Usage:
  python plot_shear_strain_sweep.py <folder> <stem> "<strain1 strain2 ...>"
    folder : run dir (".")   stem : DATANAME_INTERACTION_TOTSTEPS
    strains: space-separated list, exactly the literals used as _g tags
Only numpy + matplotlib are required.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ZX = 4  # index of sigma_zx among the 6 tensor columns: xx yy zz xy zx yz


def read_ave_time(path):
    """fix ave/time scalar file -> (timesteps[N], data[N, ncols]). [] if absent/empty."""
    if not os.path.exists(path):
        return None, None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            rows.append([float(v) for v in line.split()])
    if not rows:
        return None, None
    arr = np.array(rows)
    return arr[:, 0].astype(int), arr[:, 1:]


def read_print(path, ncol):
    """fix print file -> array[N, ncol]. None if absent/empty."""
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= ncol:
                rows.append([float(v) for v in parts[:ncol]])
    return np.array(rows) if rows else None


def read_ave_chunk(path):
    """fix ave/chunk file -> list of (timestep, array[nchunks, ncols]).
    Per-timestep header has 3 fields: timestep nchunks total-count."""
    if not os.path.exists(path):
        return []
    with open(path) as f:
        lines = [l for l in f if l.strip() and not l.startswith('#')]
    snaps, i = [], 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) == 3:
            ts, nch = int(float(parts[0])), int(parts[1])
            block = []
            for j in range(1, nch + 1):
                if i + j < len(lines):
                    block.append([float(v) for v in lines[i + j].split()])
            if block:
                snaps.append((ts, np.array(block)))
            i += nch + 1
        else:
            i += 1
    return snaps


def sem(x):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    n = len(x)
    return float(np.std(x, ddof=1) / np.sqrt(n)) if n > 1 else 0.0


def plot_stress_strain(folder, stem, strains, idx_start=0):
    sd = os.path.join(folder, 'output_files', 'stress_data')
    plot_dir = os.path.join(folder, 'output_plots')
    os.makedirs(plot_dir, exist_ok=True)

    g, mp, ep, mt, et = [], [], [], [], []
    for s in strains:
        tag = f"{stem}_g{s}"
        ts_p, tp = read_ave_time(os.path.join(sd, f"stress_tensor_polymer_{tag}.dat"))
        ts_s, tsv = read_ave_time(os.path.join(sd, f"stress_tensor_solvent_{tag}.dat"))
        ss = read_print(os.path.join(sd, f"shear_strain_{tag}.dat"), 4)
        if tp is None or ss is None:
            print(f"  [skip] strain {s}: missing tensor or shear_strain file")
            continue
        gamma = float(ss[-1, 3])
        bp = tp[idx_start:, ZX]
        if tsv is not None:
            n = min(len(bp), len(tsv[idx_start:, ZX]))
            bt = bp[:n] + tsv[idx_start:idx_start + n, ZX]
        else:
            bt = bp
        g.append(gamma)
        mp.append(float(np.mean(bp))); ep.append(sem(bp))
        mt.append(float(np.mean(bt))); et.append(sem(bt))
        print(f"  gamma={gamma:.4f}  <sig_p_zx>={np.mean(bp):+.5f} +/- {sem(bp):.5f} "
              f"(n_blocks={len(bp)})")

    if not g:
        print("  No usable per-strain stress data found — nothing to plot.")
        return None

    g = np.array(g); mp = np.array(mp); ep = np.array(ep)
    mt = np.array(mt); et = np.array(et)
    order = np.argsort(g)
    g, mp, ep, mt, et = g[order], mp[order], ep[order], mt[order], et[order]

    # Through-origin weighted fit  sigma = G * gamma
    if len(g) >= 2:
        w = 1.0 / np.clip(ep, 1e-12, None) ** 2
        G = float(np.sum(w * g * mp) / np.sum(w * g ** 2))
        G_err = float(np.sqrt(1.0 / np.sum(w * g ** 2)))
    elif len(g) == 1:
        G, G_err = mp[0] / g[0], ep[0] / g[0]
    else:
        G = G_err = np.nan
    print(f"\n  Through-origin shear modulus  G = {G:.4f} +/- {G_err:.4f}  (LJ, polymer)")

    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.errorbar(g, mp, yerr=ep, fmt='o', ms=9, color='steelblue', capsize=5, lw=2,
                label=r'$\sigma_{p,zx}$ polymer')
    ax.errorbar(g, mt, yerr=et, fmt='s', ms=8, color='purple', capsize=5, lw=2,
                alpha=0.85, label=r'$\sigma_{\mathrm{tot},zx}$ total')
    if np.isfinite(G):
        gg = np.linspace(0, g.max() * 1.05, 100)
        ax.plot(gg, G * gg, '--', color='steelblue', lw=1.6, alpha=0.8,
                label=fr'$G={G:.3f}\pm{G_err:.3f}$')
    ax.axhline(0, color='k', lw=0.8, alpha=0.3); ax.axvline(0, color='k', lw=0.8, alpha=0.3)
    ax.set_xlabel(r'$\gamma_{zx}$'); ax.set_ylabel(r'$\sigma_{zx}$ (LJ)')
    ax.set_title('Shear stress-strain sweep'); ax.legend(fontsize=13); ax.grid(alpha=0.3)
    out = os.path.join(plot_dir, f"{stem}_stress_strain_curve.png")
    fig.tight_layout(); fig.savefig(out, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f"  Saved: {out}")
    return G


def plot_sigma_zx_profiles(folder, stem, strains):
    """Overlay polymer sigma_zx(x) (last snapshot) for each strain."""
    sd = os.path.join(folder, 'output_files', 'stress_data')
    plot_dir = os.path.join(folder, 'output_plots')
    os.makedirs(plot_dir, exist_ok=True)

    gvals = [float(s) for s in strains]
    norm = Normalize(vmin=min(gvals), vmax=max(gvals)) if len(set(gvals)) > 1 else Normalize(0, 1)
    cmap = plt.cm.viridis
    fig, ax = plt.subplots(figsize=(9, 6))
    plotted = False
    for s in strains:
        snaps = read_ave_chunk(os.path.join(sd, f"stress_profile_x_polymer_{stem}_g{s}.dat"))
        if not snaps:
            continue
        _, arr = snaps[-1]                       # last production snapshot
        x = arr[:, 1]                            # Coord1 (reduced x)
        zx = arr[:, 2 + 1 + ZX]                  # cols: chunk,coord,Ncount, then 6 comps
        ax.plot(x, zx, '-o', ms=3, lw=1.5, color=cmap(norm(float(s))), label=fr'$\gamma={s}$')
        plotted = True
    if not plotted:
        plt.close(fig)
        print("  No stress_profile_x_polymer files found — skipping profile plot.")
        return
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.4)
    ax.set_xlabel(r'$\hat{x}$ (gap, reduced)'); ax.set_ylabel(r'$\sigma_{p,zx}(x)$ (LJ)')
    ax.set_title('Polymer shear-stress profile across the gap'); ax.legend(fontsize=12)
    ax.grid(alpha=0.3)
    out = os.path.join(plot_dir, f"{stem}_sigma_zx_profiles.png")
    fig.tight_layout(); fig.savefig(out, dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print('Usage: python plot_shear_strain_sweep.py <folder> <stem> "<s1 s2 ...>"')
        sys.exit(1)
    folder = sys.argv[1]
    stem = sys.argv[2]
    strains = sys.argv[3].split()
    print(f"Shear sweep post-processing: stem={stem}  strains={strains}")
    try:
        plot_stress_strain(folder, stem, strains)
        plot_sigma_zx_profiles(folder, stem, strains)
    except Exception as e:
        # Never let a quick-look plot fail the whole post-processing chain.
        print(f"  plot_shear_strain_sweep.py warning: {e}")
