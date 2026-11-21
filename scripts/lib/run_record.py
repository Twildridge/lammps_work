#!/usr/bin/env python3
"""run_record.py -- the run log behind the "Recent runs" table of the README.

Every job writes one small JSON record, runs/<cluster>_<jobid>.json, when it starts and again when
it ends (run_record_start / run_record_end in lib/run_record.sh, called by run_lammps.sh and
continue_sim.sh). On a cluster the records are kept in $RUNS_DIR/.run_records/ and published to
GitHub by lib/run_record.sh without touching the working tree. `lsync` (git_sync.sh) rebuilds the
tables from the records on every machine.

    run_record.py start   <store> <work_dir>            record of the job starting in work_dir
    run_record.py end     <store> <work_dir> <status>   close it (complete | partial | failed | incomplete)
    run_record.py refresh <store> <runs_dir> [days]     records for runs that have none (run_info.txt newer
                                                        than `days`, default 30) and close records whose
                                                        job is gone
    run_record.py table   [repo_root]                   rewrite the README section and runs/README.md

Standard library only, Python 3.6+ (it runs on the clusters' system python).
The laptop side: lib/uniaxial.py sync_pull writes runs/local/<deck>__<RUN_ID>.json naming the cluster
run folders a notebook's files came from; the table shows that RUN_ID next to those jobs.
"""
import glob
import json
import os
import re
import subprocess
import sys
import time

N_RECENT = 5
BEGIN, END = '<!-- RECENT_RUNS:BEGIN (written by scripts/lib/run_record.py; edits here are overwritten) -->', '<!-- RECENT_RUNS:END -->'
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ---------------------------------------------------------------------------------------------
#  records
# ---------------------------------------------------------------------------------------------
def _read(path):
    try:
        with open(path) as f:
            return f.read()
    except (OSError, UnicodeDecodeError):
        return ''


def parse_run_info(work_dir):
    """run_info.txt of a run folder -> dict (empty when there is none)."""
    out = {}
    for line in _read(os.path.join(work_dir, 'run_info.txt')).splitlines():
        k, sep, v = line.partition(':')
        if sep:
            out[k.strip()] = v.strip()
    return out


def log_facts(work_dir):
    """What log.lammps says: steps run (sum over the run commands), the timestep, the wall time."""
    steps, dt, wall, done = 0, None, None, False
    path = os.path.join(work_dir, 'log.lammps')
    try:
        with open(path, errors='replace') as f:
            for line in f:
                if line.startswith('Loop time of'):
                    m = re.search(r'for (\d+) steps', line)
                    if m:
                        steps += int(m.group(1))
                elif line.startswith('timestep'):
                    m = re.match(r'timestep\s+([0-9.eE+-]+)\s*$', line)
                    if m:
                        try:
                            dt = float(m.group(1))
                        except ValueError:
                            pass
                elif line.startswith('Total wall time'):
                    done = True
                    m = re.search(r'(\d+):(\d+):(\d+)', line)
                    if m:
                        wall = int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3))
    except OSError:
        pass
    return dict(steps=steps, dt=dt, wall=wall, done=done)


def final_data_files(work_dir):
    """The configurations the run wrote at its end (data_files/, or the run folder for old runs)."""
    names = set()
    for d in (os.path.join(work_dir, 'data_files'), work_dir):
        for pat in ('final_*.data', 'purepol_*.data', 'puresolv_*.data'):
            for p in glob.glob(os.path.join(d, pat)):
                if not os.path.islink(p) and '_snap' not in os.path.basename(p):
                    names.add(os.path.basename(p))
    return sorted(names)


def record_from_dir(work_dir):
    """A record built from what the run folder holds (run_info.txt is required)."""
    work_dir = os.path.abspath(work_dir)
    info = parse_run_info(work_dir)
    if not info:
        return None
    m = re.match(r'(\S+ \S+) on (\S+) \(cluster: ([^)]*)\)', info.get('started', ''))
    job = re.match(r'(\S+) \((.*)\)', info.get('slurm job', ''))
    vars_ = info.get('lammps vars', '')

    def var(name):
        mm = re.search(r'-var %s (\S+)' % name, vars_)
        return mm.group(1) if mm else None

    parent, base = os.path.split(work_dir)
    cont = base.startswith('continuation_')
    run_name = os.path.basename(parent) if cont else base
    folder = os.path.basename(os.path.dirname(parent) if cont else parent)
    deck_file = os.path.basename(info.get('deck', ''))
    return dict(
        cluster=(m.group(3) if m else '?'), job_id=(job.group(1) if job else 'none'), job_name=(job.group(2) if job else ''),
        deck=folder, deck_file=deck_file, run_name=run_name, continuation=(base if cont else None), run_dir=work_dir,
        dataname=var('dataname'), input_data=os.path.basename(info.get('input data', '')), interaction=var('interaction'),
        nsteps=var('nsteps'), allocation=info.get('allocation', ''), commit=info.get('lammps_work', ''), vars=vars_,
        started=(m.group(1) if m else ''), status='running', ended=None, steps_done=None, dt=None, wall_seconds=None, final_data=[])


def legacy_record_from_dir(work_dir, cluster):
    """A record for a run folder older than run_info.txt (before 2026-11-17), from its name
    <folder>_<dataname>_<interaction>_<YYYYMMDD>_<HHMMSS> and its log.lammps. No job ID is known."""
    work_dir = os.path.abspath(work_dir)
    parent, base = os.path.split(work_dir)
    folder = os.path.basename(parent)
    m = re.match(r'(.*)_(\d{4})(\d\d)(\d\d)_(\d\d)(\d\d)(\d\d)$', base)
    if not m or not os.path.isfile(os.path.join(work_dir, 'log.lammps')):
        return None
    rest = m.group(1)
    interaction = None
    mi = re.match(r'(.*)_(\d+(?:\.\d+)?_\d+(?:\.\d+)?)$', rest)
    if mi:
        rest, interaction = mi.group(1), mi.group(2)
    dataname = rest[len(folder) + 1:] if rest.startswith(folder + '_') else None
    if not dataname:                                      # folder renamed since: ask the log
        with open(os.path.join(work_dir, 'log.lammps'), errors='replace') as f:
            for i, line in enumerate(f):
                mm = re.match(r'read_data\s+(\S+)', line)
                if mm and '$' not in mm.group(1):
                    dataname = re.sub(r'\.data$', '', os.path.basename(mm.group(1)))
                    break
                if i > 5000:
                    break
    return dict(
        cluster=cluster or '?', job_id='none', job_name='', deck=folder, deck_file='', run_name=base, continuation=None,
        run_dir=work_dir, dataname=dataname or rest, input_data='', interaction=interaction, nsteps=None, allocation='',
        commit='', vars='', backfilled=True, started='%s-%s-%s %s:%s:%s' % m.groups()[1:],
        status='running', ended=None, steps_done=None, dt=None, wall_seconds=None, final_data=[])


def close_record(rec, work_dir, status):
    L = log_facts(work_dir)
    rec.update(status=status, steps_done=L['steps'] or None, dt=L['dt'], final_data=final_data_files(work_dir))
    try:
        end_t = os.path.getmtime(os.path.join(work_dir, 'log.lammps'))
    except OSError:
        end_t = time.time()
    if status != 'running':
        rec['ended'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_t if L['done'] or status == 'ended' else time.time()))
        if L['wall'] is not None:
            rec['wall_seconds'] = L['wall']
        elif rec.get('started'):
            try:
                rec['wall_seconds'] = int(time.mktime(time.strptime(rec['ended'], '%Y-%m-%d %H:%M:%S'))
                                          - time.mktime(time.strptime(rec['started'], '%Y-%m-%d %H:%M:%S')))
            except ValueError:
                pass
    return rec


def record_name(rec):
    jid = rec['job_id']
    if jid in ('none', '?', ''):                          # not a batch job: name it after the folder
        jid = re.sub(r'[^A-Za-z0-9_.-]', '_', rec['run_name'] + ('_' + rec['continuation'] if rec.get('continuation') else ''))
    return '%s_%s.json' % (rec['cluster'], jid)


def write_record(store, rec):
    os.makedirs(store, exist_ok=True)
    path = os.path.join(store, record_name(rec))
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(rec, f, indent=1, sort_keys=True)
        f.write('\n')
    os.replace(tmp, path)
    return path


def cmd_start(store, work_dir):
    rec = record_from_dir(work_dir)
    if rec is None:
        print('run_record: no run_info.txt in %s' % work_dir, file=sys.stderr)
        return 1
    print('run record: ' + write_record(store, rec))
    return 0


def cmd_end(store, work_dir, status):
    rec = record_from_dir(work_dir)
    if rec is None:
        return 1
    old = os.path.join(store, record_name(rec))
    if os.path.exists(old):
        try:
            with open(old) as f:
                rec.update({k: v for k, v in json.load(f).items() if k in ('started',) and v})
        except ValueError:
            pass
    if status == 'complete' and 'STATUS INCOMPLETE' in _read(os.path.join(work_dir, 'log.lammps'))[-200000:]:
        status = 'incomplete'                             # the deck stopped at its wall-clock budget (a continuation follows)
    print('run record: ' + write_record(store, close_record(rec, work_dir, status)) + ' (' + status + ')')
    return 0


def _job_alive(job_id):
    try:
        out = subprocess.run(['squeue', '-h', '-j', str(job_id), '-o', '%T'], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, timeout=30).stdout.decode()
        return bool(out.strip())
    except (OSError, subprocess.SubprocessError):
        return True                                       # cannot tell: leave the record alone


def cmd_refresh(store, runs_dir, days=30):
    """Records for runs that have none, and an end for records whose job is gone."""
    os.makedirs(store, exist_ok=True)
    have = {}
    for p in glob.glob(os.path.join(store, '*.json')):
        try:
            with open(p) as f:
                have[p] = json.load(f)
        except ValueError:
            pass
    # <deck>_stage folders are the notebooks' sync staging area (copies of run folders), not runs; a
    # record made from one is replaced by the real folder's. Two folders that give the same record
    # name would otherwise overwrite each other on alternate passes: the first one recorded stays.
    staged = lambda d: os.path.basename(os.path.dirname(d or '')).endswith('_stage')
    have = {p: r for p, r in have.items() if not staged(r.get('run_dir'))}
    known = {r.get('run_dir') for r in have.values()}
    names = {os.path.basename(p) for p in have}
    cutoff, n_new, n_closed = time.time() - float(days) * 86400, 0, 0
    for info in glob.glob(os.path.join(runs_dir, '*', '*', 'run_info.txt')) + glob.glob(os.path.join(runs_dir, '*', '*', 'continuation_*', 'run_info.txt')):
        d = os.path.dirname(os.path.abspath(info))
        if d in known or staged(d) or os.path.getmtime(info) < cutoff:
            continue
        rec = record_from_dir(d)
        if rec is None or record_name(rec) in names:
            continue
        names.add(record_name(rec))
        alive = rec['job_id'] not in ('none', '?', '') and _job_alive(rec['job_id'])
        L = log_facts(d)
        write_record(store, rec if alive else close_record(rec, d, 'complete' if L['done'] else 'ended'))
        n_new += 1
    # run folders older than run_info.txt: from the folder name and the log, when the log is no longer written
    cut_day = time.strftime('%Y%m%d', time.localtime(cutoff))
    for log in glob.glob(os.path.join(runs_dir, '*', '*', 'log.lammps')):
        d = os.path.dirname(os.path.abspath(log))
        m = re.search(r'_(\d{8})_\d{6}$', d)
        if d in known or staged(d) or not m or m.group(1) < cut_day or os.path.exists(os.path.join(d, 'run_info.txt')):
            continue
        if time.time() - os.path.getmtime(log) < 3600:    # may still be running
            continue
        rec = legacy_record_from_dir(d, os.environ.get('CLUSTER', ''))
        if rec is None or record_name(rec) in names:
            continue
        names.add(record_name(rec))
        write_record(store, close_record(rec, d, 'complete' if log_facts(d)['done'] else 'ended'))
        n_new += 1
    for p, rec in have.items():
        if rec.get('status') == 'running' and not (rec.get('job_id') not in ('none', '?', '') and _job_alive(rec['job_id'])):
            d = rec.get('run_dir', '')
            if os.path.isdir(d):
                write_record(store, close_record(rec, d, 'complete' if log_facts(d)['done'] else 'ended'))
                n_closed += 1
    if n_new or n_closed:
        print('run records: %d added, %d closed' % (n_new, n_closed))
    return 0


# ---------------------------------------------------------------------------------------------
#  tables
# ---------------------------------------------------------------------------------------------
def lattice(dataname):
    d = dataname or ''
    if 'puresolv' in d or 'solvent_box' in d or 'pure_solvent' in d:
        return 'no network (solvent)'
    net = 'random' if 'random' in d else 'diamond rot45' if 'rot45' in d else 'diamond'
    m = re.search(r'(\d+)beads', d)
    return net + (', %s beads/strand' % m.group(1) if m else '')


def _hms(sec):
    if sec is None:
        return ''
    sec = int(sec)
    return '%d h %02d min' % (sec // 3600, sec % 3600 // 60)


def _steps(rec):
    s = rec.get('steps_done')
    if not s:
        return ('planned %s (production)' % rec['nsteps']) if rec.get('nsteps') and rec.get('status') == 'running' else ''
    txt = ('%.2f M steps' % (s / 1e6)) if s >= 1e5 else '%d steps' % s
    if rec.get('dt'):
        txt += ' = %s τ' % format(int(round(s * rec['dt'])), ',')
    return txt


STATUS = {'running': 'running', 'complete': 'complete', 'partial': 'partial (wall limit)', 'failed': 'failed',
          'incomplete': 'stopped at wall budget, continued', 'ended': 'ended without finishing'}


def load_all(repo):
    recs = []
    for p in glob.glob(os.path.join(repo, 'runs', '*.json')):
        try:
            with open(p) as f:
                recs.append(json.load(f))
        except ValueError:
            pass
    links = {}                                            # cluster run folder -> laptop RUN_ID(s)
    for p in glob.glob(os.path.join(repo, 'runs', 'local', '*.json')):
        try:
            with open(p) as f:
                L = json.load(f)
        except ValueError:
            continue
        for rn in L.get('run_names', []):
            links.setdefault(rn, [])
            if L.get('run_id') and L['run_id'] not in links[rn]:
                links[rn].append(L['run_id'])
    recs.sort(key=lambda r: (r.get('started') or '', r.get('job_id') or '', r.get('cluster') or '', r.get('run_name') or ''), reverse=True)
    return recs, links


def cross_reference(repo):
    """Link the RUN_ID folders of ../flow_data_local/raw_data/<deck>/ to the recorded jobs without a
    cluster login, for data that is already on this computer. Writes runs/local/<deck>__<RUN_ID>.json
    (the same file lib/uniaxial.py sync_pull writes; a link written by a real pull is left alone).

    A folder named like a cluster run folder is that run. Otherwise a job matches when it is of the
    same deck and a file in the folder carries its tag, _<dataname>_<interaction>_<steps>; of several
    such jobs (the same data file run twice) each file is given the newest one that had started when
    the file was downloaded, which is the one the sync would have staged. The output depends only on
    the records and the files, so a second pass changes nothing."""
    raw = os.path.join(os.path.dirname(os.path.abspath(repo)), 'flow_data_local', 'raw_data')
    if not os.path.isdir(raw):
        return 0
    recs, _ = load_all(repo)
    if not recs:
        return 0
    by_name = {r.get('run_name'): r for r in recs}
    out_dir, n = os.path.join(repo, 'runs', 'local'), 0
    for deck_dir in sorted(glob.glob(os.path.join(raw, '*', ''))):
        deck = os.path.basename(deck_dir.rstrip(os.sep))
        same_deck = [r for r in recs if r.get('deck') == deck and r.get('dataname') and r.get('started')]
        for run_dir in sorted(glob.glob(os.path.join(deck_dir, '*', ''))):
            run_id = os.path.basename(run_dir.rstrip(os.sep))
            path = os.path.join(out_dir, '%s__%s.json' % (deck, run_id))
            try:
                with open(path) as f:
                    if 'pulled' in json.load(f):          # written by a real pull: exact, keep it
                        continue
            except (OSError, ValueError):
                pass
            names, dataname, cluster = set(), None, None
            if run_id in by_name:
                r = by_name[run_id]
                names, dataname, cluster = {run_id}, r.get('dataname'), r.get('cluster')
            elif same_deck:
                try:
                    files = [(e.name, e.stat().st_mtime) for e in os.scandir(run_dir) if e.is_file() and not e.name.startswith('.')]
                except OSError:
                    files = []
                for r in same_deck:
                    r['_tag'] = re.compile('_' + re.escape(r['dataname'])
                                           + ('_' + re.escape(r['interaction']) if r.get('interaction') else '') + r'_\d')
                    r['_t0'] = time.mktime(time.strptime(r['started'][:19], '%Y-%m-%d %H:%M:%S'))
                for fname, mtime in files:
                    cands = [r for r in same_deck if r['_t0'] <= mtime and r['_tag'].search(fname)]
                    if cands:
                        r = max(cands, key=lambda r: r['_t0'])
                        names.add(r['run_name'])
                        dataname, cluster = r['dataname'], r.get('cluster')
            if not names:
                continue
            os.makedirs(out_dir, exist_ok=True)
            text = json.dumps(dict(run_id=run_id, deck=deck, dataname=dataname, cluster=cluster, run_names=sorted(names),
                                   matched='from the file names and download times in flow_data_local (lib/run_record.py cross_reference)'),
                              indent=1, sort_keys=True) + '\n'
            if _read(path) != text:
                with open(path, 'w') as f:
                    f.write(text)
                n += 1
    if n:
        print('run log: %d RUN_ID folder(s) of flow_data_local linked to recorded jobs (runs/local/)' % n)
    return 0


def table(recs, links):
    if not recs:
        return '*No runs recorded yet.*'
    rows = ['| Started | Cluster | Job | Deck | Lattice | Initial data file | Final data file(s) | Simulated | Wall time | Status | Laptop RUN_ID |',
            '|---|---|---|---|---|---|---|---|---|---|---|']
    for r in recs:
        deck = r.get('deck', '') + (' (`%s`)' % r['deck_file'] if r.get('deck_file') and r['deck_file'] != r.get('deck', '') + '.lmp' else '')
        if r.get('continuation'):
            deck += ', continuation'
        init = r.get('input_data') or ((r.get('dataname') or '') + '.data')
        rows.append('| %s | %s | %s | %s | %s | `%s` | %s | %s | %s | %s | %s |' % (
            (r.get('started') or '')[:16], r.get('cluster', ''), '' if r.get('job_id') in ('none', '?') else r.get('job_id', ''), deck, lattice(r.get('dataname')), init,
            '<br>'.join('`%s`' % n for n in r.get('final_data') or []), _steps(r), _hms(r.get('wall_seconds')),
            STATUS.get(r.get('status'), r.get('status') or ''),
            ', '.join('`%s`' % x for x in links.get(r.get('run_name'), []))))
    return '\n'.join(rows)


def cmd_table(repo=REPO):
    try:
        cross_reference(repo)
    except Exception as e:                                # never stop the sync over the cross-reference
        print('run log: flow_data_local not cross-referenced (%s: %s)' % (type(e).__name__, e))
    recs, links = load_all(repo)
    os.makedirs(os.path.join(repo, 'runs'), exist_ok=True)
    full = ('# Run log\n\nEvery job that `scripts/run_lammps.sh` or `scripts/continue_sim.sh` started since 2026-11-20, newest first, '
            'on any cluster. One JSON record per job in this folder (`<cluster>_<jobid>.json`), written by the job when it starts and '
            'when it ends; `lsync` rebuilds this table and the five-row copy in the top-level README. "Laptop RUN_ID" is filled in '
            'once an analysis notebook has pulled that run\'s files (`runs/local/`). How it works: `docs/run_log.md`.\n\n'
            + table(recs, links) + '\n')
    _write_if_changed(os.path.join(repo, 'runs', 'README.md'), full)
    readme = os.path.join(repo, 'README.md')
    text = _read(readme)
    if BEGIN in text and END in text:
        block = (BEGIN + '\n' + table(recs[:N_RECENT], links) + '\n\n'
                 + 'The %d most recent of %d recorded jobs; all of them: [`runs/README.md`](runs/README.md).\n' % (min(N_RECENT, len(recs)), len(recs)) + END)
        _write_if_changed(readme, text[:text.index(BEGIN)] + block + text[text.index(END) + len(END):])
    return 0


def _write_if_changed(path, text):
    if _read(path) != text:
        with open(path, 'w') as f:
            f.write(text)


def main(argv):
    if len(argv) >= 4 and argv[1] == 'start':
        return cmd_start(argv[2], argv[3])
    if len(argv) >= 5 and argv[1] == 'end':
        return cmd_end(argv[2], argv[3], argv[4])
    if len(argv) >= 4 and argv[1] == 'refresh':
        return cmd_refresh(argv[2], argv[3], *(argv[4:5]))
    if len(argv) >= 2 and argv[1] == 'table':
        return cmd_table(*(argv[2:3]))
    sys.exit(__doc__)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
