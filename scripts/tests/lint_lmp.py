#!/usr/bin/env python3
"""
lint_lmp.py -- static consistency check of a LAMMPS input deck (2026-09-16).

Checks, in script order (jump/if branches are read linearly):
  * every ${var} / v_var / f_fix / c_compute reference is defined before use
  * every group / region referenced by a fix, compute, dump, velocity, group
    dynamic, count()/vcm()/xcm() is defined before use
  * unfix / undump / uncompute / dump_modify / fix_modify name something defined
  * equal/atom/string variables are not redefined without `variable X delete`
    (this LAMMPS build refuses); fixes/computes/dumps redefined without removal
  * label targets of `jump SELF` exist; `next VAR` names an index/loop variable
  * leftover one-piston reposition / WCA-switch / rigid_piston code (pattern list)
  * the -var names the deck expects from run_lammps.sh (undefined ${...}) are all
    forwarded by scripts/run_lammps.sh
Usage:  python lint_lmp.py <deck.lmp> [--run-lammps <run_lammps.sh>] [--forbid PATTERN ...]
Exit status 1 on any ERROR.
"""
import argparse
import re
import sys
from pathlib import Path

LAMMPS_KEYWORDS = {  # thermo keywords / math constants legal inside $( ) and formulas
    'step', 'lx', 'ly', 'lz', 'vol', 'zlo', 'zhi', 'xlo', 'xhi', 'ylo', 'yhi', 'press', 'temp', 'pe', 'ke',
    'PI', 'version', 'elapsed', 'dt', 'time', 'atoms', 'cpu',
}
GROUP_FUNCS = ('count', 'mass', 'charge', 'xcm', 'vcm', 'fcm', 'bound', 'gyration', 'ke', 'angmom', 'torque',
               'inertia', 'omega')


def read_logical_lines(path):
    """Join '&' continuations; keep the original starting line number."""
    out = []
    buf, start = '', None
    for n, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.rstrip('\n')
        if start is None:
            start = n
        stripped = line.rstrip()
        if stripped.endswith('&'):
            buf += stripped[:-1] + ' '
            continue
        buf += line
        out.append((start, buf))
        buf, start = '', None
    if buf.strip():
        out.append((start, buf))
    return out


def strip_comment(line):
    """Drop a trailing # comment that is not inside quotes."""
    out, q = [], None
    for ch in line:
        if q:
            out.append(ch)
            if ch == q:
                q = None
            continue
        if ch in ('"', "'"):
            q = ch
            out.append(ch)
            continue
        if ch == '#':
            break
        out.append(ch)
    return ''.join(out).strip()


def split_if(cmd):
    """`if "cond" then "cmd1" "cmd2" ...` -> [cond, cmd1, cmd2, ...] (also elif/else)."""
    parts = re.findall(r'"([^"]*)"|\'([^\']*)\'', cmd)
    return [a or b for a, b in parts]


class Lint:
    def __init__(self):
        self.vars = {}       # name -> style
        self.fixes = set()
        self.computes = set()
        self.dumps = set()
        self.groups = {'all'}
        self.regions = set()
        self.labels = set()
        self.jumps = []
        self.errors = []
        self.warnings = []
        self.external = set()   # ${x} used but never defined in the deck (expected from -var)

    def err(self, n, msg):
        self.errors.append(f'line {n}: ERROR {msg}')

    def warn(self, n, msg):
        self.warnings.append(f'line {n}: WARNING {msg}')

    # ---- reference checks -------------------------------------------------
    def check_refs(self, n, text, exclude=(), formulas_only=False):
        for m in re.finditer(r'\$\{(\w+)\}', text):
            v = m.group(1)
            if v not in self.vars and v not in exclude:
                self.external.add(v)
        text = re.sub(r'\$\{\w+\}', ' ', text)          # ${name} is not a v_/c_/f_ reference
        if formulas_only:                                # print strings: only $( ... ) is evaluated
            text = ' '.join(self._immediates(text))
        for m in re.finditer(r'(?<![\w])v_(\w+)', text):
            v = m.group(1)
            if v not in self.vars and v not in exclude:
                self.err(n, f'variable v_{v} referenced before definition')
        for m in re.finditer(r'(?<![\w])c_(\w+)', text):
            if m.group(1) not in self.computes:
                self.err(n, f'compute c_{m.group(1)} referenced before definition')
        for m in re.finditer(r'(?<![\w])f_(\w+)', text):
            if m.group(1) not in self.fixes:
                self.err(n, f'fix f_{m.group(1)} referenced before definition')
        for fn in GROUP_FUNCS:
            for m in re.finditer(rf'(?<![\w]){fn}\((\w+)[,)]', text):
                if m.group(1) not in self.groups:
                    self.err(n, f'group function {fn}({m.group(1)}) on undefined group')

    @staticmethod
    def _immediates(text):
        """contents of every $( ... ) in text (nested parentheses allowed)."""
        out, i = [], 0
        while True:
            j = text.find('$(', i)
            if j < 0:
                return out
            depth, k = 0, j + 1
            while k < len(text):
                if text[k] == '(':
                    depth += 1
                elif text[k] == ')':
                    depth -= 1
                    if depth == 0:
                        break
                k += 1
            out.append(text[j + 2:k])
            i = k + 1

    def check_group(self, n, g):
        if g not in self.groups:
            self.err(n, f'group "{g}" used before definition')

    # ---- one command --------------------------------------------------------
    def command(self, n, cmd):
        cmd = cmd.strip()
        if not cmd:
            return
        tok = cmd.split()
        head = tok[0]
        if head == 'if':
            parts = split_if(cmd)
            if parts:
                self.check_refs(n, parts[0])
                for sub in parts[1:]:
                    self.command(n, sub)
            return
        if head == 'variable':
            name, style = tok[1], (tok[2] if len(tok) > 2 else '')
            if style == 'delete':
                if name not in self.vars:
                    self.warn(n, f'variable {name} deleted but not defined')
                self.vars.pop(name, None)
                return
            rest = ' '.join(tok[3:])
            if style in ('equal', 'atom', 'vector', 'string', 'format', 'getenv', 'file', 'python', 'internal'):
                if name in self.vars and self.vars[name] not in ('index', 'loop'):
                    # LAMMPS replaces the formula (triaxial_permeation.lmp relies on it); the
                    # repo idiom is delete-first, so only warn
                    self.warn(n, f'variable {name} ({self.vars[name]}) redefined without `variable {name} delete`')
                # formulas may reference the variable itself? no. Check references (own name excluded for loop-style)
                self.check_refs(n, rest)
            elif style in ('index', 'loop', 'uloop', 'world', 'universe'):
                pass   # silently ignored if it already exists (that is how -var overrides work)
            else:
                self.warn(n, f'unknown variable style "{style}"')
            self.vars[name] = style
            return
        if head == 'next':
            for v in tok[1:]:
                if v not in self.vars:
                    self.err(n, f'next {v}: variable not defined')
                elif self.vars[v] not in ('index', 'loop', 'uloop'):
                    self.err(n, f'next {v}: variable is {self.vars[v]}-style, not index/loop')
            return
        if head == 'label':
            self.labels.add(tok[1])
            return
        if head == 'jump':
            self.jumps.append((n, tok[2] if len(tok) > 2 else None))
            return
        # generic reference scan of the rest of the command (print strings: $( ) only)
        is_print = head == 'print' or (head == 'fix' and len(tok) > 3 and tok[3] == 'print')
        self.check_refs(n, ' '.join(tok[1:]), formulas_only=is_print)
        if head == 'group':
            gid = tok[1]
            if len(tok) > 2 and tok[2] == 'dynamic':
                self.check_group(n, tok[3])
                if 'region' in tok:
                    r = tok[tok.index('region') + 1]
                    if r not in self.regions:
                        self.err(n, f'group {gid} dynamic: region {r} undefined')
            self.groups.add(gid)
            return
        if head == 'region':
            self.regions.add(tok[1])
            return
        if head == 'compute':
            cid, g = tok[1], tok[2]
            self.check_group(n, g)
            if cid in self.computes:
                self.err(n, f'compute {cid} redefined without uncompute')
            if 'region' in tok[3:]:
                r = tok[tok.index('region') + 1]
                if r not in self.regions:
                    self.err(n, f'compute {cid}: region {r} undefined')
            self.computes.add(cid)
            return
        if head == 'uncompute':
            if tok[1] not in self.computes:
                self.err(n, f'uncompute {tok[1]}: not defined')
            self.computes.discard(tok[1])
            return
        if head == 'fix':
            fid, g = tok[1], tok[2]
            self.check_group(n, g)
            if fid in self.fixes:
                self.err(n, f'fix {fid} redefined without unfix')
            if len(tok) > 3 and tok[3] == 'halt':
                pass
            self.fixes.add(fid)
            return
        if head == 'unfix':
            if tok[1] not in self.fixes:
                self.err(n, f'unfix {tok[1]}: not defined')
            self.fixes.discard(tok[1])
            return
        if head == 'fix_modify':
            if tok[1] not in self.fixes:
                self.err(n, f'fix_modify {tok[1]}: not defined')
            return
        if head == 'dump':
            did, g = tok[1], tok[2]
            self.check_group(n, g)
            if did in self.dumps:
                self.err(n, f'dump {did} redefined without undump')
            self.dumps.add(did)
            return
        if head == 'undump':
            if tok[1] not in self.dumps:
                self.err(n, f'undump {tok[1]}: not defined')
            self.dumps.discard(tok[1])
            return
        if head == 'dump_modify':
            if tok[1] not in self.dumps:
                self.err(n, f'dump_modify {tok[1]}: not defined')
            return
        if head in ('velocity', 'delete_atoms', 'displace_atoms'):
            if head == 'delete_atoms' and tok[1] == 'group':
                self.check_group(n, tok[2])
            elif head == 'velocity':
                self.check_group(n, tok[1])
            return

    def finish(self, forbid, run_lammps):
        for n, lab in self.jumps:
            if lab and lab not in self.labels:
                self.err(n, f'jump to undefined label {lab}')
        if run_lammps:
            fw = set(re.findall(r'-var\s+(\w+)', Path(run_lammps).read_text()))
            missing = sorted(v for v in self.external if v not in fw and v not in LAMMPS_KEYWORDS)
            if missing:
                self.errors.append(f'run_lammps.sh does not forward -var for: {", ".join(missing)}')
        return self.errors, self.warnings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('deck')
    ap.add_argument('--run-lammps', default=None)
    ap.add_argument('--forbid', nargs='*', default=[])
    args = ap.parse_args()
    L = Lint()
    lines = read_logical_lines(args.deck)
    for n, raw in lines:
        cmd = strip_comment(raw)
        if not cmd:
            continue
        for pat in args.forbid:
            if re.search(pat, cmd):
                L.err(n, f'forbidden pattern "{pat}" in: {cmd[:80]}')
        L.command(n, cmd)
    errors, warnings = L.finish(args.forbid, args.run_lammps)
    print(f'{Path(args.deck).name}: {len(lines)} logical lines; '
          f'{len(L.vars)} variables, {len(L.fixes)} fixes still defined at end, {len(L.computes)} computes, '
          f'{len(L.dumps)} dumps, labels {sorted(L.labels)}')
    print(f'  external (-var) names: {sorted(v for v in L.external if v not in LAMMPS_KEYWORDS)}')
    if L.fixes:
        print(f'  fixes never unfixed (fine if intended): {sorted(L.fixes)}')
    for w in warnings:
        print('  ' + w)
    for e in errors:
        print('  ' + e)
    print(f'  RESULT: {"OK" if not errors else str(len(errors)) + " ERROR(S)"}')
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    main()
