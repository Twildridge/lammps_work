#!/usr/bin/env python3
"""collect_results.py -- logs of the allocation benchmark -> results.csv + verdict (2026-11-05).

    python collect_results.py [--logs DIR] [--out results.csv]
    python collect_results.py --su-total          # only print the SUs spent so far

Standard library only (Expanse python/3.8).  A timed run counts only if it passes:
  complete   LAMMPS reached "Total wall time"
  safe       Dangerous builds = 0 (a dangerous build can silently miss pair interactions)
  step 0     PE and pressure at step 0 match every other run of the same prepared state.  They
             start from identical coordinates, so a decomposition or a communication cutoff that
             loses an interaction shows up here, before the Langevin noise makes trajectories differ.
"""
import argparse
import csv
import re
import statistics
from collections import defaultdict
from pathlib import Path

STEP0_TOL = 1e-6      # relative
# Steps to a converged slab: the budgets of chain_length_testing_sweep.batch.  N = 25 is not in
# that table; it uses the same rule (14 M x mesh ratio = 14 M x 2.88).
PRODUCTION_STEPS = {5: 14e6, 10: 24e6, 15: 30e6, 20: 35e6, 25: 40e6, 30: 45e6}
CORES_PER_NODE = 128

COLUMNS = ('tag role nbeads scale nodes tasks_per_node cpus_per_task partition map comm_cut rep job '
           'cores_charged ranks job_wall_s complete mpirun_rc natoms n_polymer n_solvent n_frozen steps '
           'loop_s steps_per_s katom_step_s su_per_Mstep hours_per_Mstep '
           'pct_pair pct_bond pct_neigh pct_comm pct_output pct_modify pct_other imbalance_pair imbalance_neigh '
           'ave_ghost ave_neighs neigh_builds dangerous step0_pe step0_press mean_mobile_temp '
           'lx ly lz valid why_invalid').split()
TIMING_ROW = re.compile(r'^(Pair|Bond|Neigh|Comm|Output|Modify|Other)\s*\|(.*)$')


def num(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def parse_log(path):
    """Fields of the LAST run in the log, the timed production block."""
    row = {}
    lines = path.read_text(errors='replace').splitlines()
    for line in lines:
        m = re.match(r'\s*(\d+) atoms in group (polymer|solvent|support|piston)\b', line)
        if m:
            key = {'polymer': 'n_polymer', 'solvent': 'n_solvent'}.get(m.group(2))
            if key:
                row[key] = int(m.group(1))
            else:
                row['n_frozen'] = row.get('n_frozen', 0) + int(m.group(1))
    loop_i = [i for i, l in enumerate(lines) if l.startswith('Loop time of')]
    if not loop_i:
        return row
    i = loop_i[-1]
    m = re.match(r'Loop time of (\S+) on (\d+) procs for (\d+) steps with (\d+) atoms', lines[i])
    if m:
        row.update(loop_s=float(m.group(1)), steps=int(m.group(3)), natoms=int(m.group(4)))
    # thermo table of this run: walk back to the header
    h = i
    while h > 0 and not lines[h].lstrip().startswith('Step '):
        h -= 1
    head = lines[h].split()
    table = []
    for line in lines[h + 1:i]:
        vals = line.split()
        if len(vals) == len(head) and num(vals[0]) is not None:
            table.append(dict(zip(head, map(num, vals))))
    if table:
        row['step0_pe'] = table[0].get('PotEng')
        row['step0_press'] = table[0].get('Press')
        temps = [t['c_mobile_temp'] for t in table[1:] if t.get('c_mobile_temp') is not None]
        if temps:
            row['mean_mobile_temp'] = statistics.mean(temps)
        for box in ('Lx', 'Ly', 'Lz'):
            row[box.lower()] = table[-1].get(box)
    for line in lines[i:]:
        m = re.search(r'([\d.eE+-]+) timesteps/s', line)
        if line.startswith('Performance:') and m:
            row['steps_per_s'] = float(m.group(1))
            k = re.search(r'([\d.eE+-]+) katom-step/s', line)
            row['katom_step_s'] = float(k.group(1)) if k else None
        m = TIMING_ROW.match(line)
        if m:
            name = m.group(1).lower()
            f = [x.strip() for x in m.group(2).split('|')]
            if 'pct_' + name not in row:
                row['pct_' + name] = num(f[-1])
                # slowest rank / average rank; the "Other" row has no min/avg/max columns
                if name in ('pair', 'neigh') and len(f) >= 3:
                    avg, mx = num(f[1]), num(f[2])
                    if avg and mx:
                        row['imbalance_' + name] = mx / avg
        if line.startswith('Nghost:'):
            row['ave_ghost'] = num(re.search(r'ave\s*=\s*(\S+)|^Nghost:\s*(\S+)', line).group(1)
                                   or line.split()[1])
        if line.startswith('Ave neighs/atom'):
            row['ave_neighs'] = num(line.split('=')[1])
        if line.startswith('Neighbor list builds'):
            row['neigh_builds'] = int(line.split('=')[1])
        if line.startswith('Dangerous builds'):
            row['dangerous'] = int(line.split('=')[1])
    return row


def read_all(logs):
    rows = []
    for meta_file in sorted(logs.glob('*.meta')):
        row = dict(line.split('=', 1) for line in meta_file.read_text().splitlines() if '=' in line)
        log_file = meta_file.with_suffix('.log')
        if log_file.exists():
            row.update(parse_log(log_file))
        nodes, tasks, cpus = int(row['nodes']), int(row['tasks_per_node']), int(row['cpus_per_task'])
        row['ranks'] = nodes * tasks
        # compute charges whole nodes; shared (and a Mac) charges the cores asked for
        row['cores_charged'] = nodes * CORES_PER_NODE if row['partition'] == 'compute' else nodes * tasks * cpus
        rate = row.get('steps_per_s')
        if rate:
            row['hours_per_Mstep'] = 1e6 / rate / 3600
            row['su_per_Mstep'] = row['hours_per_Mstep'] * row['cores_charged']
        rows.append(row)
    validate(rows)
    return rows


def validate(rows):
    groups = defaultdict(list)
    for r in rows:
        if r['role'] == 'bench':
            groups[(r['nbeads'], r['scale'])].append(r)
    for r in rows:
        why = []
        if r.get('complete') != '1':
            why.append('incomplete')
        if r['role'] == 'bench' and r.get('dangerous') not in (0, None):
            why.append(f"dangerous={r['dangerous']}")
        r['why_invalid'] = why
    for group in groups.values():
        for key in ('step0_pe', 'step0_press'):
            vals = [r[key] for r in group if r.get(key) is not None]
            if not vals:
                continue
            ref = statistics.median(vals)
            for r in group:
                if r.get(key) is not None and abs(r[key] - ref) > STEP0_TOL * max(abs(ref), 1e-12):
                    r['why_invalid'].append(f'{key} {r[key]:.10g} vs {ref:.10g}')
    for r in rows:
        r['valid'] = int(not r['why_invalid'])
        r['why_invalid'] = '; '.join(r['why_invalid'])


def su_total(rows):
    return sum(float(r.get('job_wall_s') or 0) * r['cores_charged'] / 3600 for r in rows)


def label(r):
    tail = '' if r['partition'] == 'compute' else f" [{r['partition']}]"
    return f"{r['nodes']} node(s) x {r['tasks_per_node']} tasks x {r['cpus_per_task']} cpus{tail}"


def verdict(rows):
    timed = [r for r in rows if r['role'] == 'bench' and r['valid'] and r.get('steps_per_s')]
    bad = [r for r in rows if not r['valid']]
    if bad:
        print(f'\nEXCLUDED {len(bad)} run(s):')
        for r in bad:
            print(f"  {r['tag']}: {r['why_invalid']}")
    configs = defaultdict(list)
    for r in timed:
        configs[(int(r['nbeads']), float(r['scale']), r['map'], float(r['comm_cut']), r['partition'],
                 int(r['nodes']), int(r['tasks_per_node']), int(r['cpus_per_task']))].append(r)
    medians = defaultdict(list)
    for key, group in configs.items():
        speed = statistics.median(g['steps_per_s'] for g in group)
        medians[key[:4]].append((speed, group[0], len(group)))
    print('\nFor each N, scale, mapping and comm cutoff (median of repeats):')
    print('  fastest = highest steps/s;  value = fewest SU/M-step among runs within 10 % of the fastest')
    for key in sorted(medians):
        n, scale, mapping, cut = key
        entries = medians[key]
        top = max(e[0] for e in entries)
        fastest = max(entries, key=lambda e: e[0])
        value = min((e for e in entries if e[0] >= 0.9 * top), key=lambda e: e[1]['su_per_Mstep'])
        budget = PRODUCTION_STEPS[n]     # the budget of the scale 1.0 slab of this N
        print(f"\nN={n} scale={scale} map={mapping} comm_cutoff={cut}  ({fastest[1].get('natoms')} atoms, {len(entries)} layouts)")
        for name, (speed, r, reps) in (('fastest', fastest), ('value  ', value)):
            hours = budget / speed / 3600
            su = hours * r['cores_charged']
            print(f"  {name} {label(r):42s} {speed:8.1f} steps/s  {r['su_per_Mstep']:7.1f} SU/M-step  "
                  f"{budget / 1e6:.0f} M steps: {hours:6.1f} h, {su:8.0f} SU  (x{reps})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--logs', default=str(Path.home() / 'Documents/lammps_runs/allocation_benchmark/logs'))
    ap.add_argument('--out', default=None, help='default: results.csv beside the logs folder')
    ap.add_argument('--su-total', action='store_true')
    args = ap.parse_args()
    logs = Path(args.logs)
    rows = read_all(logs)
    if args.su_total:
        print(f'{su_total(rows):.0f}')
        return
    out = Path(args.out) if args.out else logs.parent / 'results.csv'
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    print(f'{len(rows)} runs -> {out}')
    print(f'SUs spent so far: {su_total(rows):.0f} (job wall time x cores charged, runs with logs only)')
    verdict(rows)


if __name__ == '__main__':
    main()
