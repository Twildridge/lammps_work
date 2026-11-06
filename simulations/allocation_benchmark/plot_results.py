#!/usr/bin/env python3
"""plot_results.py -- results.csv of the allocation benchmark -> figure + bead-cost fit (2026-11-05).

    python plot_results.py results.csv [--out allocation_benchmark.png] [--title "..."]

Only valid runs are used (collect_results.py), with repeats reduced to their median.

  A  hours per 1 M steps vs beads per chain, one line per node count (best layout at each N);
     dashed = production settings when the csv has both
  B  the same cost vs total beads, including the 0.5x and 2x slabs when present
  C  hours for the whole production run of each N: budget of chain_length_testing_sweep.batch /
     the fastest speed.  Longer strands need more steps even if each step costs the same.
  D  parallel scaling at N = 5 and N = 30: steps/s vs total cores, one line per tasks x cpus layout

The fit: seconds per step = a N_polymer + b N_other, one fit per layout with enough runs, to answer
whether polymer or solvent beads drive the cost.  N_other is solvent + frozen sheets: at one slab
size they rise and fall together (r = 0.995) and a fit cannot tell them apart.  They are reported
separately only for data that separates them (|r| < 0.95, which the 0.5x and 2x slabs may give).
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

SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e6e5e0'
PROBE = (5, 30)


def load(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            if r['role'] != 'bench' or r['valid'] != '1' or not r['steps_per_s']:
                continue
            for k in ('nbeads', 'nodes', 'tasks_per_node', 'cpus_per_task', 'natoms', 'n_polymer', 'n_solvent', 'n_frozen'):
                r[k] = int(r[k])
            for k in ('scale', 'comm_cut', 'steps_per_s'):
                r[k] = float(r[k])
            rows.append(r)
    groups = defaultdict(list)
    for r in rows:
        groups[(r['nbeads'], r['scale'], r['partition'], r['nodes'], r['tasks_per_node'], r['cpus_per_task'],
                r['map'], r['comm_cut'])].append(r)
    out = []
    for g in groups.values():
        r = dict(g[0])
        r['steps_per_s'] = statistics.median(x['steps_per_s'] for x in g)
        r['nrep'] = len(g)
        r['hours_per_M'] = 1e6 / r['steps_per_s'] / 3600
        out.append(r)
    return out


def settings_of(r):
    return f"map {r['map']}, cutoff {r['comm_cut']:g}" if r['partition'] != 'local' else f"cutoff {r['comm_cut']:g}"


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.title.set_color(INK)


def line(ax, xs, ys, color, dashed=False, label=None):
    order = np.argsort(xs)
    ax.plot(np.array(xs)[order], np.array(ys)[order], color=color, linewidth=2, linestyle='--' if dashed else '-',
            marker='o', markersize=6, markeredgecolor=SURFACE, markeredgewidth=1.2, label=label)


def group_key(rows):
    """Node count on Expanse; thread count on a Mac, where there is only ever one rank."""
    if len({r['nodes'] for r in rows}) > 1 or any(r['tasks_per_node'] > 1 for r in rows):
        return (lambda r: r['nodes']), (lambda v: f'{v} node' + ('s' if v > 1 else ''))
    return (lambda r: r['cpus_per_task']), (lambda v: f'{v} thread' + ('s' if v > 1 else ''))


def panel_a(ax, rows):
    base = [r for r in rows if r['scale'] == 1.0 and r['partition'] != 'shared']
    key, name = group_key(base)
    groups = sorted({key(r) for r in base})
    settings = sorted({settings_of(r) for r in base}, key=lambda s: ('2.5' not in s, s))
    for gi, g in enumerate(groups):
        for si, s in enumerate(settings):
            best = {}
            for r in base:
                if key(r) == g and settings_of(r) == s and (r['nbeads'] not in best or r['steps_per_s'] > best[r['nbeads']]['steps_per_s']):
                    best[r['nbeads']] = r
            if best:
                label = name(g) + (f', {s}' if len(settings) > 1 else '')
                line(ax, list(best), [best[n]['hours_per_M'] for n in best], SERIES[gi % 8], dashed=si > 0, label=label)
    atoms = {r['nbeads']: r['natoms'] for r in base}
    ax.set_xticks(sorted(atoms))
    ax.set_xticklabels([f'{n}\n{atoms[n] / 1e3:.0f}k' for n in sorted(atoms)])
    ax.set_xlabel('beads per chain  (total beads below)')
    ax.set_ylabel('hours per 1 M steps')
    ax.set_ylim(bottom=0)
    ax.set_title('A. Simulation time vs beads per chain', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, handlelength=3.5)


def fastest_per_n(rows, scale=1.0):
    best = {}
    for r in rows:
        if r['scale'] == scale and (r['nbeads'] not in best or r['steps_per_s'] > best[r['nbeads']]['steps_per_s']):
            best[r['nbeads']] = r
    return best


def panel_b(ax, rows):
    scales = sorted({r['scale'] for r in rows})
    for i, scale in enumerate(scales):
        best = fastest_per_n(rows, scale)
        ax.scatter([b['natoms'] / 1e3 for b in best.values()], [b['hours_per_M'] for b in best.values()], s=60,
                   color=SERIES[i], edgecolor=SURFACE, linewidth=1.2, label=f'{scale:g}x slab', zorder=3)
        # alternate above and below along x: neighbours in bead count sit close together
        ordered = sorted(best.items(), key=lambda kv: kv[1]['natoms'])
        for j, (n, b) in enumerate(ordered):
            ax.annotate(f'N={n}', (b['natoms'] / 1e3, b['hours_per_M']), xytext=(0, 8 if j % 2 else -15),
                        textcoords='offset points', ha='center', fontsize=8, color=INK2)
    ax.set_xlabel('total beads (thousands)')
    ax.set_ylabel('hours per 1 M steps, fastest layout')
    ax.set_xlim(0, max(r['natoms'] for r in rows) / 1e3 * 1.08)
    ax.set_ylim(bottom=0)
    ax.set_title('B. Simulation time vs total beads', loc='left', fontsize=11)
    if len(scales) > 1:
        ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc='lower right')


def panel_c(ax, rows):
    best = fastest_per_n(rows)
    ns = sorted(best)
    hours = [PRODUCTION_STEPS[n] / best[n]['steps_per_s'] / 3600 for n in ns]
    ax.bar(range(len(ns)), hours, width=0.6, color=SERIES[0])
    for i, n in enumerate(ns):
        ax.annotate(f'{hours[i]:.0f} h\n{PRODUCTION_STEPS[n] / 1e6:.0f} M steps', (i, hours[i]), xytext=(0, 3),
                    textcoords='offset points', ha='center', fontsize=8, color=INK2)
    ax.set_xticks(range(len(ns)))
    ax.set_xticklabels(ns)
    ax.set_xlabel('beads per chain')
    ax.set_ylabel('hours for the whole run, fastest layout')
    ax.set_ylim(0, max(hours) * 1.22)
    ax.grid(False, axis='x')
    ax.set_title('C. Full production run', loc='left', fontsize=11)


def panel_d(ax, rows):
    base = [r for r in rows if r['scale'] == 1.0 and r['nbeads'] in PROBE]
    settings = sorted({settings_of(r) for r in base}, key=lambda s: ('2.5' not in s, s))
    layouts = sorted({(r['tasks_per_node'], r['cpus_per_task']) for r in base}, reverse=True)
    used = 0
    for n, dashed in zip(PROBE, (False, True)):
        s = settings[0] if settings else None
        for t, c in layouts:
            pick = [r for r in base if r['nbeads'] == n and settings_of(r) == s and (r['tasks_per_node'], r['cpus_per_task']) == (t, c)]
            cores = defaultdict(list)
            for r in pick:
                cores[r['nodes'] * t * c].append(r['steps_per_s'])
            if cores:
                color = SERIES[layouts.index((t, c)) % 8]
                line(ax, list(cores), [max(v) for v in cores.values()], color, dashed=dashed, label=f'N={n}, {t} x {c}')
                used += 1
    if len({r['nodes'] for r in base}) == 1 and all(r['tasks_per_node'] == 1 for r in base):
        # A Mac: threads are the only variable, so draw them across layouts.
        ax.clear()
        style(ax)
        for n, dashed in zip(PROBE, (False, True)):
            pick = [r for r in base if r['nbeads'] == n and settings_of(r) == settings[0]]
            line(ax, [r['cpus_per_task'] for r in pick], [r['steps_per_s'] for r in pick], SERIES[0 if not dashed else 1],
                 dashed=dashed, label=f'N={n}')
    ax.set_xscale('log', base=2)
    ax.set_xlabel('total cores')
    ticks = sorted({int(x) for ln in ax.get_lines() for x in ln.get_xdata()})
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_ylabel('steps per second')
    ax.set_ylim(bottom=0)
    ax.set_title(f'D. Parallel scaling ({settings[0] if settings else ""})', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, ncol=2, handlelength=3.5)


def least_squares(X, y):
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    dof = max(len(y) - X.shape[1], 1)
    se = np.sqrt(np.diag(np.linalg.inv(X.T @ X)) * (resid @ resid) / dof)
    return coef * 1e6, se * 1e6, np.sqrt(np.mean((resid / y) ** 2)) * 100


def composition_fit(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[(r['partition'], r['nodes'], r['tasks_per_node'], r['cpus_per_task'], r['map'], r['comm_cut'])].append(r)
    print('\nCost per bead, microseconds per bead per step (fit: seconds per step = a N_polymer + b N_other)')
    for key, g in sorted(groups.items()):
        if len({r['nbeads'] for r in g}) < 5:
            continue
        pol, sol, fro = (np.array([r[k] for r in g], float) for k in ('n_polymer', 'n_solvent', 'n_frozen'))
        y = np.array([1 / r['steps_per_s'] for r in g])
        coef, se, rms = least_squares(np.c_[pol, sol + fro], y)
        print(f"  {key[1]} node x {key[2]} x {key[3]} {key[4]} cc{key[5]:g} ({len(y)} points, rms {rms:.1f} %): "
              f"polymer {coef[0]:.3f} +/- {se[0]:.3f}, solvent + frozen {coef[1]:.3f} +/- {se[1]:.3f}, "
              f"polymer / other = {coef[0] / coef[1]:.2f}")
        r_sf = np.corrcoef(sol, fro)[0, 1]
        if abs(r_sf) < 0.95:
            coef, se, rms = least_squares(np.c_[pol, sol, fro], y)
            print(f"      separated (solvent-frozen r = {r_sf:.2f}, rms {rms:.1f} %): polymer {coef[0]:.3f} +/- {se[0]:.3f}, "
                  f"solvent {coef[1]:.3f} +/- {se[1]:.3f}, frozen {coef[2]:.3f} +/- {se[2]:.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv')
    ap.add_argument('--out', default='allocation_benchmark.png')
    ap.add_argument('--title', default='LAMMPS allocation benchmark: x/y-periodic random slab')
    args = ap.parse_args()
    rows = load(args.csv)
    if not rows:
        raise SystemExit('no valid timed runs in ' + args.csv)
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), facecolor=SURFACE)
    for ax in axes.flat:
        style(ax)
    panel_a(axes[0, 0], rows)
    panel_b(axes[0, 1], rows)
    panel_c(axes[1, 0], rows)
    panel_d(axes[1, 1], rows)
    fig.suptitle(args.title, x=0.01, ha='left', fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(args.out, dpi=170)
    print(f'wrote {Path(args.out).resolve()}')
    composition_fit(rows)


if __name__ == '__main__':
    main()
