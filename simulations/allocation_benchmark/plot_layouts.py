#!/usr/bin/env python3
"""plot_layouts.py -- results.csv of the allocation benchmark -> the MPI x OpenMP layout figure
of one slab family, plus a ranked table and the #SBATCH lines to copy (2026-11-07).

    python plot_layouts.py results.csv [--network rot45] [--out layouts_rot45.png]

For each N (one row of the figure): steps/s and SU per 1 M steps of every timed layout, bars
coloured by node count.  Only valid runs (collect_results.py) are used, repeats reduced to their
median.  The text report ranks the layouts, gives the gain over the production baseline
(4 nodes x 128 x 1, --map-by node, cutoff 5.0, or the 1-node 128 x 1 run when the baseline is
absent) and the wall hours and SUs of the full production budget (PRODUCTION_STEPS).
"""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from collect_results import PRODUCTION_STEPS  # noqa: E402
from plot_results import SERIES, SURFACE, INK, INK2, style  # noqa: E402

NODE_ORDER = (1, 2, 3, 4, 6, 8, 16)    # fixed colour per node count, never cycled


def load(path, network):
    groups = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f):
            if r['role'] != 'bench' or r['valid'] != '1' or not r['steps_per_s']:
                continue
            if r.get('network', 'random') != network or float(r['scale']) != 1.0:
                continue
            key = (int(r['nbeads']), int(r['nodes']), int(r['tasks_per_node']), int(r['cpus_per_task']),
                   r['map'], float(r['comm_cut']), r['partition'])
            groups[key].append(r)
    rows = []
    for key, g in groups.items():
        n, nodes, tasks, cpus, mapping, cut, part = key
        speed = statistics.median(float(x['steps_per_s']) for x in g)
        cores = int(g[0]['cores_charged'])
        rows.append(dict(nbeads=n, nodes=nodes, tasks=tasks, cpus=cpus, map=mapping, cut=cut, part=part, speed=speed,
                         cores=cores, natoms=int(g[0]['natoms'] or 0), nrep=len(g),
                         su_per_M=1e6 / speed / 3600 * cores, hours_per_M=1e6 / speed / 3600,
                         pct_comm=g[0].get('pct_comm'), pct_pair=g[0].get('pct_pair')))
    return rows


def layout_label(r):
    s = f"{r['tasks']}x{r['cpus']}"
    if r['map'] != 'node':
        s += f" {r['map']}"
    if r['cut'] != 5.0:
        s += f" cc{r['cut']:g}"
    if r['part'] == 'compute' and r['nodes'] <= 2:
        s += ' [compute]'      # the cross-check against the same layout on the debug partition
    return s


def sbatch_lines(r):
    return (f"#SBATCH --nodes={r['nodes']}\n#SBATCH --ntasks-per-node={r['tasks']}\n"
            f"#SBATCH --cpus-per-task={r['cpus']}   # OMP_NUM_THREADS={r['cpus']}; mpirun --map-by "
            + (f"ppr:{r['tasks']}:node:pe={r['cpus']}" if r['map'] == 'core' else f"node:pe={r['cpus']}")
            + (f"; comm_modify cutoff {r['cut']:g}" if r['cut'] != 5.0 else ''))


def panel(ax, rows, value, ylabel, title):
    labels = sorted({layout_label(r) for r in rows},
                    key=lambda s: (-int(s.split('x')[0]), 'core' in s, 'cc' in s))
    nodes = [n for n in NODE_ORDER if any(r['nodes'] == n for r in rows)]
    width = 0.8 / max(len(nodes), 1)
    x = np.arange(len(labels))
    for i, n in enumerate(nodes):
        vals = []
        for lab in labels:
            match = [r[value] for r in rows if r['nodes'] == n and layout_label(r) == lab]
            vals.append(match[0] if match else np.nan)
        pos = x - 0.4 + width * (i + 0.5)
        ax.bar(pos, vals, width=width * 0.92, color=SERIES[NODE_ORDER.index(n) % 8],
               edgecolor=SURFACE, linewidth=1, label=f"{n} node{'s' if n > 1 else ''}")
        for p, v in zip(pos, vals):
            if not np.isnan(v):
                ax.annotate(f'{v:.0f}', (p, v), xytext=(0, 2), textcoords='offset points',
                            ha='center', fontsize=7, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8, rotation=20, ha='right')
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    ax.grid(False, axis='x')
    ax.set_title(title, loc='left', fontsize=11)
    if len(nodes) > 1:
        ax.legend(frameon=False, fontsize=8, labelcolor=INK2)


def report(rows, network):
    lines = []
    for n in sorted({r['nbeads'] for r in rows}):
        rs = sorted((r for r in rows if r['nbeads'] == n), key=lambda r: -r['speed'])
        base = [r for r in rs if (r['nodes'], r['tasks'], r['cpus'], r['map'], r['cut']) == (4, 128, 1, 'node', 5.0)] \
            or [r for r in rs if (r['nodes'], r['tasks'], r['cpus'], r['map'], r['cut'], r['part']) == (1, 128, 1, 'node', 5.0, 'debug')] \
            or [r for r in rs if (r['nodes'], r['tasks'], r['cpus'], r['map'], r['cut']) == (1, 128, 1, 'node', 5.0)]
        base = base[0] if base else None
        budget = PRODUCTION_STEPS.get(n, 14e6)
        lines.append(f"\n=== N = {n} ({network}, {rs[0]['natoms']} atoms): {len(rs)} layouts, "
                     f"production budget {budget / 1e6:.0f} M steps ===")
        lines.append(f"{'nodes':>5} {'layout':<26} {'steps/s':>8} {'h/M-step':>9} {'SU/M-step':>10} "
                     f"{'full run h':>10} {'full run SU':>11} {'vs base':>8} {'comm %':>7}")
        for r in rs:
            gain = f"{r['speed'] / base['speed']:.2f}x" if base else '-'
            comm = f"{float(r['pct_comm']):.1f}" if r['pct_comm'] else '-'
            lines.append(f"{r['nodes']:>5} {layout_label(r):<26} {r['speed']:>8.1f} {r['hours_per_M']:>9.2f} "
                         f"{r['su_per_M']:>10.0f} {budget / r['speed'] / 3600:>10.1f} "
                         f"{budget / r['speed'] / 3600 * r['cores']:>11.0f} {gain:>8} {comm:>7}")
        fastest = rs[0]
        per_nodes = {}
        for r in rs:
            if r['nodes'] not in per_nodes:
                per_nodes[r['nodes']] = r      # rs is sorted by speed
        lines.append(f"\nfastest overall: {fastest['nodes']} node(s), {layout_label(fastest)}  "
                     f"({fastest['speed']:.1f} steps/s)")
        lines.append(sbatch_lines(fastest))
        for nodes, r in sorted(per_nodes.items()):
            if r is not fastest:
                lines.append(f"best on {nodes} node(s): {layout_label(r)}  ({r['speed']:.1f} steps/s, "
                             f"{r['su_per_M']:.0f} SU/M-step)")
        cheapest = min(rs, key=lambda r: r['su_per_M'])
        lines.append(f"cheapest per step: {cheapest['nodes']} node(s), {layout_label(cheapest)}  "
                     f"({cheapest['su_per_M']:.0f} SU/M-step, {cheapest['speed']:.1f} steps/s)")
        if base:
            lines.append(f"baseline ({base['nodes']} x 128 x 1 node cc5.0): {base['speed']:.1f} steps/s, "
                         f"{base['su_per_M']:.0f} SU/M-step")
    return '\n'.join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv')
    ap.add_argument('--network', default='rot45')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    rows = load(args.csv, args.network)
    if not rows:
        raise SystemExit(f'no valid timed {args.network} runs in {args.csv}')
    ns = sorted({r['nbeads'] for r in rows})
    fig, axes = plt.subplots(len(ns), 2, figsize=(12.5, 4.2 * len(ns)), facecolor=SURFACE, squeeze=False)
    for i, n in enumerate(ns):
        rs = [r for r in rows if r['nbeads'] == n]
        for ax in axes[i]:
            style(ax)
        panel(axes[i, 0], rs, 'speed', 'steps per second',
              f'N = {n}: speed  ({rs[0]["natoms"] / 1e3:.0f}k atoms)')
        panel(axes[i, 1], rs, 'su_per_M', 'SU per 1 M steps', f'N = {n}: cost')
    fig.suptitle(f'LAMMPS layout benchmark: {args.network} slab, tasks x cpus per node '
                 f'(map node unless marked; comm cutoff 5.0 unless marked)',
                 x=0.01, ha='left', fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = args.out or f'layouts_{args.network}.png'
    fig.savefig(out, dpi=170)
    print(f'wrote {Path(out).resolve()}')
    print(report(rows, args.network))


if __name__ == '__main__':
    main()
