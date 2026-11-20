#!/usr/bin/env python3
"""summarize_node_screen.py -- rows.tsv of a node screen -> ranked table, suspect list, exclude.txt.

    python3 summarize_node_screen.py [--out DIR]      (default: the newest node_screen/<date> folder)

Reads DIR/rows.tsv (one line per finished node_screen.sh job) and DIR/nodes.tsv (every node at
submission, with the jobs it got).  Speed numbers are core-microseconds per atom-step, lower is
better; a node is compared with the median of all tested nodes.  Standard library only.
"""
import argparse
import csv
import statistics
from pathlib import Path

ROOT = Path.home() / 'Documents/lammps_runs/allocation_benchmark/node_screen'

# relative thresholds: a node is suspect when it is this much slower than the median
SLOW = 1.25
HOT_C = 85            # package temperature during a LAMMPS run
LOW_MHZ = 2200        # mean all-core clock during the 40-rank runs (Xeon Gold 6148: 2.4 base, ~3.1 all-core turbo)
STRAY_BUSY = 3.0      # % CPU busy on a node that --exclusive gave us empty
THROTTLE_NEW_MS = 1000


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def median(vals):
    vals = [v for v in vals if v is not None]
    return statistics.median(vals) if vals else None


def read_rows(path):
    if not path.is_file():
        return []
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', help='screen folder (default: newest under %s)' % ROOT)
    ap.add_argument('--best', type=int, default=12, help='how many nodes to list as best')
    args = ap.parse_args()
    out = Path(args.out) if args.out else max(p for p in ROOT.iterdir() if p.is_dir())

    rows = read_rows(out / 'rows.tsv')
    nodes_at_submit = {r['node']: r for r in read_rows(out / 'nodes.tsv')}
    per = {}                                   # node -> {'core': row, 'node': row}
    for r in sorted(rows, key=lambda r: r['date']):
        per.setdefault(r['node'], {})[r['mode']] = r   # the latest row of a mode wins

    def pick(d, key, modes=('node', 'core')):
        for m in modes:
            if m in d and num(d[m].get(key)) is not None:
                return num(d[m][key])
        return None

    # the fleet's medians
    med = {k: median([pick(d, k) for d in per.values()]) for k in ('lj1_us', 'lj40_us', 'deck_us', 'scratch_mbps')}
    thr_rate = {n: (pick(d, 'thr_pkg_ms') or 0) / max(pick(d, 'uptime_d') or 1, 0.1) for n, d in per.items()}
    med_thr = median(list(thr_rate.values()))
    kernels = [d[m]['kernel'] for d in per.values() for m in d]
    common_kernel = max(set(kernels), key=kernels.count) if kernels else ''

    table = []
    for node, d in per.items():
        flags = []
        lj1, lj40, deck = pick(d, 'lj1_us'), pick(d, 'lj40_us'), pick(d, 'deck_us')
        if med['lj1_us'] and lj1 and lj1 > SLOW * med['lj1_us']:
            flags.append('SLOW_CORE(%.0f%%)' % (100 * (lj1 / med['lj1_us'] - 1)))
        if med['lj40_us'] and lj40 and lj40 > SLOW * med['lj40_us']:
            flags.append('SLOW_LJ40(%.0f%%)' % (100 * (lj40 / med['lj40_us'] - 1)))
        if med['deck_us'] and deck and deck > SLOW * med['deck_us']:
            flags.append('SLOW_DECK(%.0f%%)' % (100 * (deck / med['deck_us'] - 1)))
        pkg = max([v for v in (pick(d, 'deck_pkg_c'), pick(d, 'lj40_pkg_c'), pick(d, 'lj1_pkg_c')) if v is not None] or [0])
        if pkg >= HOT_C:
            flags.append('HOT(%dC)' % pkg)
        mhz = pick(d, 'deck_mhz_mean', ('node',)) or pick(d, 'lj40_mhz_mean', ('node',))
        if mhz and mhz < LOW_MHZ:
            flags.append('LOW_CLOCK(%dMHz)' % mhz)
        new_ms = max([v for v in (pick(d, 'thr_pkg_ms_new', (m,)) for m in d) if v is not None] or [0])
        if new_ms > THROTTLE_NEW_MS:
            flags.append('THROTTLING_NOW(%.0fs)' % (new_ms / 1000))
        if med_thr and thr_rate[node] > 3 * med_thr and thr_rate[node] > 60000:
            flags.append('THROTTLE_HISTORY(%.0fmin/day)' % (thr_rate[node] / 60000))
        busy = pick(d, 'busy_idle_pct', ('node',))
        if busy is not None and busy > STRAY_BUSY:
            flags.append('STRAY_LOAD(%.0f%%)' % busy)
        for m in d:
            if d[m]['status'] != 'ok':
                flags.append('%s:%s' % (m, d[m]['status']))
            if d[m]['kernel'] != common_kernel:
                flags.append('KERNEL(%s)' % d[m]['kernel'])
                break
        if any(d[m].get('opa') not in ('ACTIVE', None) for m in d):
            flags.append('OPA(%s)' % d[next(iter(d))].get('opa'))
        if 'node' not in d:
            flags.append('core-only')
        serious = [f for f in flags if f.split('(')[0] in
                   ('SLOW_CORE', 'SLOW_LJ40', 'SLOW_DECK', 'HOT', 'LOW_CLOCK', 'THROTTLING_NOW', 'STRAY_LOAD', 'OPA')
                   or ':' in f]
        verdict = 'suspect' if serious else ('good' if 'node' in d else 'core-ok')
        table.append(dict(node=node, verdict=verdict, deck_sps=pick(d, 'deck_steps_s'), deck=deck, lj40=lj40, lj1=lj1,
                          mhz=mhz, pkg=pkg, thr=thr_rate[node] / 60000, busy=busy, scratch=pick(d, 'scratch_mbps'),
                          uptime=pick(d, 'uptime_d'), flags=flags))

    def sort_key(t):
        return (t['verdict'] == 'suspect', t['deck'] is None, t['deck'] or 0, t['lj1'] or 0)
    table.sort(key=sort_key)

    untested = [n for n in nodes_at_submit if n not in per]
    skipped = [n for n, r in nodes_at_submit.items() if r.get('jobs') == '-']
    pending = [n for n in untested if n not in skipped]

    lines = []
    p = lines.append
    p('Pod node screen: %s   (%d of %d nodes tested, %d suspect)' % (
        out.name, len(per), len(nodes_at_submit) or len(per), sum(t['verdict'] == 'suspect' for t in table)))
    p('core-us/atom-step, lower is better; medians: LJ 1 rank %s, LJ 40 ranks %s, deck N15 %s' % tuple(
        '%.3f' % med[k] if med[k] else 'NA' for k in ('lj1_us', 'lj40_us', 'deck_us')))
    p('')
    hdr = '%-7s %-8s %8s %7s %7s %7s %6s %5s %9s %5s %7s  %s' % (
        'node', 'verdict', 'deck_s/s', 'deck', 'lj40', 'lj1', 'MHz', 'pkgC', 'thr_min/d', 'busy%', 'scr_MBs', 'flags')
    p(hdr)
    p('-' * len(hdr))
    f = lambda v, fmt: (fmt % v) if v is not None else 'NA'
    for t in table:
        p('%-7s %-8s %8s %7s %7s %7s %6s %5s %9s %5s %7s  %s' % (
            t['node'], t['verdict'], f(t['deck_sps'], '%.2f'), f(t['deck'], '%.3f'), f(t['lj40'], '%.3f'),
            f(t['lj1'], '%.3f'), f(t['mhz'], '%.0f'), f(t['pkg'], '%d'), f(t['thr'], '%.1f'),
            f(t['busy'], '%.1f'), f(t['scratch'], '%.0f'), ' '.join(t['flags'])))
    p('')
    best = [t['node'] for t in table if t['verdict'] == 'good'][:args.best]
    p('best %d (whole-node test passed, fastest deck first): %s' % (len(best), ' '.join(best)))
    suspects = [t['node'] for t in table if t['verdict'] == 'suspect']
    p('suspect: %s' % (' '.join(suspects) or 'none'))
    p('skipped at submit (down/drained): %s' % (' '.join(skipped) or 'none'))
    p('not yet tested (jobs pending): %s' % (' '.join(pending) or 'none'))
    if suspects:
        p('')
        p('exclude line for the pod block of scripts/lib/site.sh (sbatch honours it):')
        p('    export SBATCH_EXCLUDE="%s"' % ','.join(suspects))
        (out / 'exclude.txt').write_text(','.join(suspects) + '\n')
    text = '\n'.join(lines)
    print(text)
    (out / 'summary.txt').write_text(text + '\n')


if __name__ == '__main__':
    main()
