#!/usr/bin/env python3
"""plot_nodescan.py -- results.csv of the allocation benchmark -> speed and cost against the node
count at one rank per core (128 x 1 on Expanse, 40 x 1 on Pod), one line per setting (2026-11-07,
after the nodescan tier).

    python plot_nodescan.py results.csv [--network rot45] [--out nodescan_rot45.png]

Left column: steps/s against nodes; right column: SU per 1 M steps.  One row per N.  Each point
is the median of the valid repeats; the bar through it spans the slowest and fastest repeat, so
a difference inside the bar is not a difference.  The three settings are the production baseline
(--map-by node, comm cutoff 5.0), the Slurm-only change (--map-by core, i.e. ppr:128:node, cutoff
5.0) and both changes (core + cutoff 2.5).  The text report lists the medians and spreads.
"""
import argparse
import csv
import statistics
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from collect_results import PRODUCTION_STEPS  # noqa: E402
from plot_results import SERIES, SURFACE, INK, INK2, style  # noqa: E402

SETTINGS = (                       # fixed order and colour, never cycled
    ('node', 5.0, 'map node, cutoff 5.0 (production today)', SERIES[0]),
    ('core', 5.0, 'map core, cutoff 5.0 (Slurm change only)', SERIES[1]),
    ('core', 2.5, 'map core, cutoff 2.5 (both changes)', SERIES[2]),
)


def load(path, network, cluster=None):
    """-> (groups, cores, cluster): the full-node x 1 runs of one cluster (the csv's only one, or --cluster)."""
    with open(path) as f:
        rows = [r for r in csv.DictReader(f)
                if r['role'] == 'bench' and r['valid'] == '1' and r['steps_per_s']
                and r.get('network', 'random') == network and float(r['scale']) == 1.0]
    clusters = sorted({r.get('cluster') or 'expanse' for r in rows})
    if cluster is None and len(clusters) > 1:
        raise SystemExit(f'{path} holds runs of {clusters}: pick one with --cluster')
    cluster = cluster or (clusters[0] if clusters else 'expanse')
    rows = [r for r in rows if (r.get('cluster') or 'expanse') == cluster]
    cores = max((int(r.get('cores_per_node') or 128) for r in rows), default=128)
    groups = defaultdict(list)
    for r in rows:
        if int(r['tasks_per_node']) != cores or int(r['cpus_per_task']) != 1:
            continue
        key = (int(r['nbeads']), r['map'], float(r['comm_cut']), int(r['nodes']))
        groups[key].append(float(r['steps_per_s']))
    return groups, cores, cluster


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv')
    ap.add_argument('--network', default='rot45')
    ap.add_argument('--out', default=None)
    ap.add_argument('--cluster', default=None, help='expanse | pod, when the csv holds both')
    a = ap.parse_args()
    groups, cores, cluster = load(a.csv, a.network, a.cluster)
    ns = sorted({k[0] for k in groups})
    if not ns:
        raise SystemExit(f'no {cores} x 1 runs of the {a.network} slab in {a.csv}')
    fig, axes = plt.subplots(len(ns), 2, figsize=(12, 4.6 * len(ns)), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    fig.suptitle(f'LAMMPS node-count scan on {cluster}: {a.network} slab, {cores} ranks x 1 thread per node '
                 f'(point = median of the repeats, bar = slowest to fastest repeat)', color=INK, x=0.01, ha='left')
    report = []
    for row, n in enumerate(ns):
        report.append(f'\n=== N = {n} ({a.network}): speed and cost against nodes, {cores} x 1 on {cluster} ===')
        report.append(f"{'setting':<42} {'nodes':>5} {'steps/s':>8} {'spread':>14} {'SU/M-step':>10} {'reps':>4}")
        for col, (quantity, ylabel) in enumerate((('speed', 'steps per second'), ('cost', 'SU per 1 M steps'))):
            ax = axes[row][col]
            style(ax)
            for mapping, cut, label, color in SETTINGS:
                pts = sorted((k[3], v) for k, v in groups.items() if k[:3] == (n, mapping, cut))
                if not pts:
                    continue
                xs = [p[0] for p in pts]
                med = [statistics.median(p[1]) for p in pts]
                lo = [min(p[1]) for p in pts]
                hi = [max(p[1]) for p in pts]
                if quantity == 'speed':
                    y, ylo, yhi = med, lo, hi
                else:   # cost falls when speed rises: swap the ends of the bar
                    conv = lambda s, nodes: 1e6 / s / 3600 * nodes * cores  # noqa: E731
                    y = [conv(s, x) for s, x in zip(med, xs)]
                    ylo = [conv(s, x) for s, x in zip(hi, xs)]
                    yhi = [conv(s, x) for s, x in zip(lo, xs)]
                err = [[a_ - b_ for a_, b_ in zip(y, ylo)], [b_ - a_ for a_, b_ in zip(y, yhi)]]
                ax.errorbar(xs, y, yerr=err, color=color, linewidth=2, marker='o', markersize=7,
                            markeredgecolor=SURFACE, markeredgewidth=1.2, capsize=3, label=label)
                ax.annotate(f'{y[-1]:.0f}', (xs[-1], y[-1]), textcoords='offset points', xytext=(7, 0),
                            color=INK2, fontsize=9, va='center')
                if quantity == 'speed':
                    for x, s, l, h in zip(xs, med, lo, hi):
                        report.append(f'{label:<42} {x:>5} {s:>8.1f} {l:>6.0f} - {h:<6.0f} '
                                      f'{1e6 / s / 3600 * x * cores:>10.0f} {len(groups[(n, mapping, cut, x)]):>4}')
            ax.set_xscale('log', base=2)
            xticks = sorted({k[3] for k in groups if k[0] == n})
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(x) for x in xticks])
            ax.minorticks_off()
            ax.set_xlabel(f'nodes ({cores} cores each)')
            ax.set_ylabel(ylabel)
            ax.set_ylim(bottom=0)
            ax.set_title(f'N = {n}: {quantity}', loc='left', color=INK)
            if row == 0 and col == 0:
                ax.legend(frameon=False, fontsize=9, loc='upper left')
        if not any(k[:3] == (n, 'core', 2.5) for k in groups):
            continue
        best = max(((k[3], statistics.median(v)) for k, v in groups.items() if k[:3] == (n, 'core', 2.5)),
                   key=lambda t: t[1])
        cheapest = min(((k[3], statistics.median(v)) for k, v in groups.items() if k[:3] == (n, 'core', 2.5)),
                       key=lambda t: 1e6 / t[1] * t[0])
        steps = PRODUCTION_STEPS.get(n, 1e6) if isinstance(PRODUCTION_STEPS, dict) else 1e6
        report.append(f'core + 2.5: fastest at {best[0]} nodes ({best[1]:.0f} steps/s, '
                      f'{steps / best[1] / 3600:.1f} h and {steps / best[1] / 3600 * best[0] * cores:.0f} SU per '
                      f'{steps / 1e6:g} M steps); cheapest at {cheapest[0]} node(s) ({cheapest[1]:.0f} steps/s, '
                      f'{steps / cheapest[1] / 3600:.1f} h, {steps / cheapest[1] / 3600 * cheapest[0] * cores:.0f} SU)')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = a.out or f'nodescan_{a.network}.png'
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    print(f'wrote {out}')
    print('\n'.join(report))


if __name__ == '__main__':
    main()
