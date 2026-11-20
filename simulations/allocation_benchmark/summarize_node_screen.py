#!/usr/bin/env python3
"""summarize_node_screen.py -- node screen results -> ranked table, verdicts, exclude list.

    python3 summarize_node_screen.py [--out DIR] [--best N]     one screen folder (default: the newest)
    python3 summarize_node_screen.py --current [--quiet]        every screen folder merged -> node_screen/current/

Sources per folder: rows*.tsv (one line per finished node_screen.sh job, any schema version), nodes.tsv
(every node at submission) and <node>/<mode>/report.txt, which fills in what a job that was killed before
writing its row still measured, and the busiest process that was not ours.  With --current the latest
measurement of each node and mode, over all folders, counts; a bad verdict therefore stands until a newer
measurement clears it.  current/ holds summary.txt, exclude.txt (suspect + retest, what scripts/lib/site.sh
exports as SBATCH_EXCLUDE on Pod), suspects.txt, retest.txt and last_screen.txt (when the last 1-core and
whole-node measurements were taken; submit.sh starts a new screen when they are old).

Speed numbers are core-microseconds per atom-step (Loop time x ranks / (steps x atoms)); lower is
better; a node is compared with the median of all tested nodes.  Standard library only.

Verdicts: suspect = thermal trouble (HOT, THROTTLING_NOW, or slow with a low clock) or a slow/failed
whole-node run on a node without the exporter problem; retest = the Prometheus node exporter was spinning
when the node was measured (CSC restarts it, 2026-11-20), so its performance flags are stale and the node
needs a fresh whole-node run; good = whole-node run passed; core-ok = only the 1-core job so far.
"""
import argparse
import csv
import datetime as dt
import re
import statistics
from pathlib import Path

ROOT = Path.home() / 'Documents/lammps_runs/allocation_benchmark/node_screen'

SLOW = 1.25           # slower than the median by this factor
HOT_C = 92            # package temperature during a LAMMPS run; 85..91 is reported as WARM only
LOW_MHZ = 2200        # mean all-core clock during the 40-rank runs (Xeon Gold 6148: 2.4 base, ~3.1 all-core turbo)
STRAY_BUSY = 1.5      # % CPU busy on a node that --exclusive gave us empty (one core of 40 = 2.5)
STRAY_PROC = 50       # % CPU of the busiest process that is not ours
COMM_FRAC = 0.40      # Comm share of the 40-rank LJ loop; a healthy node sits near 0.05-0.10
THROTTLE_NEW_MS = 1000
THERMAL = ('HOT', 'THROTTLING_NOW')                       # suspect whatever else is going on
PERF = ('SLOW_LJ40', 'SLOW_DECK', 'COMM_WAIT', 'STRAY_LOAD', 'FAILED', 'TIMEOUT', 'OPA', 'SLOW_CORE_X')
# SLOW_CORE counts only from the exclusive whole-node job (SLOW_CORE_X); on a shared node it is contention.
# LOW_CLOCK counts only together with a slow run: the mean clock also covers the idle data-reading phase.


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def median(vals):
    vals = [v for v in vals if v is not None]
    return statistics.median(vals) if vals else None


def read_tsv(path):
    if not path.is_file():
        return []
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def parse_report(path):
    """What a report.txt holds even when the job died: run results, status, the foreign process."""
    d = {}
    text = path.read_text(errors='replace')
    m = re.search(r'^node=(\S+) job=(\S+) mode=(\S+) date=(\d{4}-\d\d-\d\dT\d\d:\d\d)', text, re.M)
    if m:
        d['job'], d['date'] = m.group(2), m.group(4)
    loops = {}
    for m in re.finditer(r'^(lj1|lj40|lj39|deck): ranks=\d+ .*?loop_s=(\S+) steps=(\S+) atoms=(\S+) core_us_per_atom_step=(\S+)', text, re.M):
        d[m.group(1) + '_us'] = num(m.group(5)); loops[m.group(1)] = num(m.group(2))
        if m.group(1) == 'deck' and num(m.group(2)) and num(m.group(3)):
            d['deck_steps_s'] = num(m.group(3)) / num(m.group(2))
    for m in re.finditer(r'^(lj1|lj40|lj39|deck): avg_s .*?comm=(\S+)', text, re.M):
        if num(m.group(2)) is not None and loops.get(m.group(1)):
            d[m.group(1) + '_comm_frac'] = num(m.group(2)) / loops[m.group(1)]
    for m in re.finditer(r'^(lj40|deck): during run clock_MHz mean=(\S+) max=\S+ package_C_max=(\S+)', text, re.M):
        d[m.group(1) + '_mhz_mean'] = num(m.group(2)); d[m.group(1) + '_pkg_c'] = num(m.group(3))
    m = re.search(r'^lj1: during run clock_MHz mean=\S+ max=(\S+) package_C_max=(\S+)', text, re.M)
    if m:
        d['lj1_mhz_max'] = num(m.group(1)); d['lj1_pkg_c'] = num(m.group(2))
    # the top block: the busiest process that is not ours, and the node exporter's share on its own
    blk = re.search(r'== processes by CPU.*?\n(.*?)\n==', text, re.S)
    best, exporter = (0.0, 'none'), 0.0
    for line in (blk.group(1).splitlines() if blk else []):
        f = line.split()
        if len(f) >= 12 and f[0].isdigit() and not f[1].startswith('pollard') and num(f[8]) is not None:
            if num(f[8]) > best[0]:
                best = (num(f[8]), '%s (%s)' % (f[11], f[1]))
            if re.search(r'prometh|node_exporter', f[11] + f[1]):
                exporter = max(exporter, num(f[8]))
    d['foreign_pct'], d['foreign_cmd'] = best
    d['exporter_pct'] = exporter
    m = re.search(r'^job_wall_s=\d+ status=(\S+)', text, re.M)
    d['status'] = m.group(1) if m else ('timeout' if 'signal received' in text else 'killed')
    for key, pat in (('kernel', r'^kernel=(\S+)'), ('uptime_d', r'^uptime_days=(\S+)'), ('thr_pkg_ms', r'^throttle .*package_time_ms=(\d+)'),
                     ('busy_idle_pct', r'busy_pct_5s=(\S+)'), ('scratch_mbps', r'^write_MB_per_s=(\S+)'), ('opa', r'^state=\d+: (\S+)')):
        m = re.search(pat, text, re.M)
        if m:
            d[key] = m.group(1) if key in ('kernel', 'opa') else num(m.group(1))
    m = re.search(r'\(new during job: core=\d+ package=\d+ ms=(\d+)\)', text)
    if m:
        d['thr_pkg_ms_new'] = num(m.group(1))
    return d


def jobid(v):
    """559427, '559427' and 559427.0 are the same job."""
    return str(int(v)) if isinstance(v, float) else str(v)


def load(folders):
    """per[node][mode] = the latest measurement over the folders; nodes_at_submit = union of nodes.tsv."""
    per, nodes_at_submit = {}, {}
    for out in folders:
        for f in sorted(out.glob('rows*.tsv')):
            for r in read_tsv(f):
                d = {k: (num(v) if num(v) is not None else v) for k, v in r.items() if k}
                d['_source'] = 'row'
                have = per.setdefault(r['node'], {}).get(r['mode'])
                if have is None or str(d.get('date', '')) >= str(have.get('date', '')):
                    per[r['node']][r['mode']] = d
        for rep in sorted(out.glob('node*/*/report.txt')):
            node, mode = rep.parent.parent.name, rep.parent.name
            d = parse_report(rep)
            have = per.setdefault(node, {}).get(mode)
            if have is None or str(d.get('date', '')) > str(have.get('date', '')):
                d['_source'] = 'report'; per[node][mode] = d
            elif jobid(have.get('job')) == jobid(d.get('job')):   # the row is authoritative; the report adds what it lacks
                for k, v in d.items():
                    if have.get(k) in (None, 'NA', '') and v is not None:
                        have[k] = v
        for r in read_tsv(out / 'nodes.tsv'):
            nodes_at_submit[r['node']] = r
    return per, nodes_at_submit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', help='screen folder (default: newest under %s)' % ROOT)
    ap.add_argument('--current', action='store_true', help='merge every screen folder; write node_screen/current/')
    ap.add_argument('--quiet', action='store_true')
    ap.add_argument('--best', type=int, default=15, help='how many nodes to list as best')
    args = ap.parse_args()
    dated = sorted(p for p in ROOT.iterdir() if p.is_dir() and re.fullmatch(r'\d{8}', p.name)) if ROOT.is_dir() else []
    if args.current:
        folders, out = dated, ROOT / 'current'
        out.mkdir(parents=True, exist_ok=True)
    else:
        out = Path(args.out) if args.out else dated[-1]
        folders = [out]
    per, nodes_at_submit = load(folders)

    def val(d, key, modes=('node', 'core')):
        for m in modes:
            if m in d:
                v = d[m].get(key)
                if isinstance(v, float):
                    return v
        return None

    med = {k: median([val(d, k) for d in per.values()]) for k in ('lj1_us', 'lj40_us', 'lj39_us', 'deck_us')}
    thr_rate = {n: (val(d, 'thr_pkg_ms') or 0) / max(val(d, 'uptime_d') or 1, 0.1) for n, d in per.items()}
    med_thr = median(list(thr_rate.values()))
    kernels = [d[m].get('kernel') for d in per.values() for m in d if d[m].get('kernel')]
    common_kernel = max(set(kernels), key=kernels.count) if kernels else ''

    table = []
    for node, d in per.items():
        flags = []
        lj1, lj40, lj39, deck = (val(d, k) for k in ('lj1_us', 'lj40_us', 'lj39_us', 'deck_us'))
        for name, v, m in (('SLOW_CORE', lj1, med['lj1_us']), ('SLOW_LJ40', lj40, med['lj40_us']), ('SLOW_DECK', deck, med['deck_us'])):
            if m and v and v > SLOW * m:
                if name == 'SLOW_CORE' and 'node' in d and val(d, 'lj1_us', ('node',)) == v:
                    name = 'SLOW_CORE_X'
                flags.append('%s(%.0f%%)' % (name, 100 * (v / m - 1)))
        pkg = max([v for v in (val(d, 'deck_pkg_c'), val(d, 'lj40_pkg_c'), val(d, 'lj1_pkg_c')) if v is not None] or [0])
        if pkg >= HOT_C:
            flags.append('HOT(%dC)' % pkg)
        elif pkg >= 85:
            flags.append('WARM(%dC)' % pkg)
        mhz = val(d, 'deck_mhz_mean', ('node',)) or val(d, 'lj40_mhz_mean', ('node',))
        if mhz and mhz < LOW_MHZ:
            flags.append('LOW_CLOCK(%dMHz)' % mhz)
        new_ms = max([val(d, 'thr_pkg_ms_new', (m,)) or 0 for m in d] or [0])
        if new_ms > THROTTLE_NEW_MS:
            flags.append('THROTTLING_NOW(%.0fs)' % (new_ms / 1000))
        if med_thr and thr_rate[node] > 3 * med_thr and thr_rate[node] > 60000:
            flags.append('THROTTLE_HISTORY(%.0fmin/day)' % (thr_rate[node] / 60000))
        busy = val(d, 'busy_idle_pct', ('node',))
        if busy is not None and busy > STRAY_BUSY:
            flags.append('STRAY_LOAD(%.1f%%)' % busy)
        fpct = max([val(d, 'foreign_pct', (m,)) or 0 for m in d] or [0])
        fcmd = next((d[m].get('foreign_cmd') for m in ('node', 'core') if m in d and (val(d, 'foreign_pct', (m,)) or 0) == fpct), '')
        epct = max([val(d, 'exporter_pct', (m,)) or 0 for m in d] or [0])
        if epct >= STRAY_PROC:
            flags.append('EXPORTER_SPIN(%.0f%%)' % epct)
        cf = val(d, 'lj40_comm_frac', ('node',))
        if cf is not None and cf > COMM_FRAC:
            flags.append('COMM_WAIT(%.0f%%)' % (100 * cf))
        if lj40 and lj39 and lj39 < lj40 / 1.5:
            flags.append('39_RANKS_FIX(%.1fx)' % (lj40 / lj39))
        for m in d:
            st = str(d[m].get('status', 'ok'))
            if 'timeout' in st or st == 'killed':
                flags.append('TIMEOUT(%s)' % m)
            elif st != 'ok':
                flags.append('FAILED(%s:%s)' % (m, st))
        kern = next((d[m].get('kernel') for m in d if d[m].get('kernel')), None)
        if kern and kern != common_kernel:
            flags.append('KERNEL(%s)' % kern)
        opa = next((d[m].get('opa') for m in d if d[m].get('opa')), None)
        if opa and opa != 'ACTIVE':
            flags.append('OPA(%s)' % opa)
        if 'node' not in d:
            flags.append('core-only')
        kinds = [f.split('(')[0] for f in flags]
        slow = any(k in ('SLOW_LJ40', 'SLOW_DECK', 'SLOW_CORE_X') for k in kinds)
        exporter = 'EXPORTER_SPIN' in kinds
        thermal = any(k in THERMAL for k in kinds) or ('LOW_CLOCK' in kinds and slow and not exporter)
        perf = any(k in PERF for k in kinds)
        if thermal:
            verdict = 'suspect'
        elif exporter:
            verdict = 'retest'
        elif perf:
            verdict = 'suspect'
        else:
            verdict = 'good' if 'node' in d else 'core-ok'
        when = max(str(d[m].get('date', '')) for m in d)[:16]
        table.append(dict(node=node, verdict=verdict, deck_sps=val(d, 'deck_steps_s'), deck=deck, lj40=lj40, lj39=lj39, lj1=lj1,
                          mhz=mhz, pkg=pkg, thr=thr_rate[node] / 60000, busy=busy, fpct=fpct, fcmd=str(fcmd).split(' ')[0],
                          when=when, flags=flags))
    order = {'good': 0, 'core-ok': 1, 'retest': 2, 'suspect': 3}
    table.sort(key=lambda t: (order[t['verdict']], t['deck'] is None, t['deck'] or 0, t['lj1'] or 0))

    skipped = [n for n, r in nodes_at_submit.items() if r.get('jobs') == '-' and n not in per]
    pending = [n for n in nodes_at_submit if n not in per and n not in skipped]
    suspects = [t['node'] for t in table if t['verdict'] == 'suspect']
    retest = [t['node'] for t in table if t['verdict'] == 'retest']
    best = [t['node'] for t in table if t['verdict'] == 'good'][:args.best]

    L = []
    L.append('Pod node screen %s: %d of %d nodes with results (%d with a whole-node run), %d suspect, %d to retest' % (
        'current (%d folders)' % len(folders) if args.current else out.name, len(per), len(nodes_at_submit) or len(per),
        sum(t['lj40'] is not None for t in table), len(suspects), len(retest)))
    L.append('core-us per atom-step, lower is better.  medians: LJ 1 rank %s, LJ 40 ranks %s, LJ 39 ranks %s, deck N15 %s' % tuple(
        ('%.3f' % med[k]) if med[k] else 'NA' for k in ('lj1_us', 'lj40_us', 'lj39_us', 'deck_us')))
    L.append('')
    hdr = '%-7s %-8s %8s %7s %7s %7s %7s %5s %5s %7s %5s %-20s %-11s %s' % (
        'node', 'verdict', 'deck_s/s', 'deck', 'lj40', 'lj39', 'lj1', 'MHz', 'pkgC', 'thr_m/d', 'busy', 'top_foreign(%cpu)', 'measured', 'flags')
    L.append(hdr); L.append('-' * len(hdr))
    f = lambda v, fmt: (fmt % v) if v is not None else 'NA'
    for t in table:
        L.append('%-7s %-8s %8s %7s %7s %7s %7s %5s %5s %7s %5s %-20s %-11s %s' % (
            t['node'], t['verdict'], f(t['deck_sps'], '%.1f'), f(t['deck'], '%.3f'), f(t['lj40'], '%.3f'), f(t['lj39'], '%.3f'),
            f(t['lj1'], '%.3f'), f(t['mhz'], '%.0f'), f(t['pkg'], '%d'), f(t['thr'], '%.1f'), f(t['busy'], '%.1f'),
            ('%s(%.0f)' % (t['fcmd'][:13], t['fpct'])) if t['fpct'] else '-', t['when'][5:16], ' '.join(t['flags'])))
    L.append('')
    L.append('best (whole-node run passed, fastest deck first): %s' % (' '.join(best) or 'none yet'))
    L.append('suspect (thermal or slow, exclude): %s' % (' '.join(suspects) or 'none'))
    L.append('retest (node exporter was spinning when measured): %s' % (' '.join(retest) or 'none'))
    L.append('skipped at submit (down/drained): %s' % (' '.join(skipped) or 'none'))
    L.append('no result yet (jobs pending): %s' % (' '.join(pending) or 'none'))
    if suspects or retest:
        L.append('')
        L.append('excluded by lsbatch on Pod (#SBATCH --exclude line from batch_header in scripts/lib/site.sh):')
        L.append('    %s' % ','.join(suspects + retest))
        if retest:
            L.append('(suspects %s; the rest until retested)' % (','.join(suspects) or 'none'))
    text = '\n'.join(L)
    if not args.quiet:
        print(text)
    (out / 'summary.txt').write_text(text + '\n')
    (out / 'exclude.txt').write_text(','.join(suspects + retest) + '\n')
    (out / 'suspects.txt').write_text(','.join(suspects) + '\n')
    (out / 'retest.txt').write_text(' '.join(retest) + '\n')
    if args.current:
        lines = []
        for mode in ('core', 'node'):
            dates = [str(d[mode].get('date', '')) for d in per.values() if mode in d and d[mode].get('date')]
            newest = max(dates) if dates else ''
            try:
                epoch = int(dt.datetime.strptime(newest[:16], '%Y-%m-%dT%H:%M').timestamp())
            except ValueError:
                epoch = 0
            lines += ['%s_date=%s' % (mode, newest[:16]), '%s_epoch=%d' % (mode, epoch)]
        lines.append('nodes_measured=%d' % len(per))
        lines.append('written=%s' % dt.datetime.now().strftime('%Y-%m-%dT%H:%M'))
        (out / 'last_screen.txt').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
