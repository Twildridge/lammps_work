#!/usr/bin/env python
"""screen_gate.py -- mechanical pass/fail check of a screen_isotropy smoke-test pair (2026-11-13).

The production jobs of the gel screen (screen_isotropy.batch with GATE=<screen dir>/smoke_passed)
refuse to run until that file exists.  This script decides whether to create it, from the two
smoke runs of one job: the same short deck run with -sf omp and with plain styles (..._sfoff).

Checks, for BOTH runs unless noted:
  1. log.lammps reached SCREEN_COMPLETE, no ERROR, no 'FENE bond too long', no lost/missing atoms
  2. all 13 stages logged; stress and dims series complete to the last step, finite
  3. every hold sits at its target strain (gamma = xy/ly, e_x = lx/lx0 - 1, ly = ly0/(1+e)) and
     returns to zero; box volume and lz unchanged in every hold
  4. mobile temperature within 3 % of 1.0, mean pressure between 1.3 and 1.7
  5. reference dump: expected number of frames, polymer only
  6. deformation dump: support and piston stay planes (z spread < 1e-4) and their atoms keep their
     fractional (lamda) x, y coordinates in every frame -> the plates deform affinely with the box
  7. (pair) step-0 polymer-group stress and box pressure identical with and without -sf omp
     (relative difference < 1e-9): the bond-virial under-tally of the hybrid/omp decks is absent
  8. (omp run) measured steps/s gives the largest production gel (354k atoms, 3.8 M steps) a wall
     time under 4.5 h of the 5 h job limit

Result: <screen dir>/gate_check.log always; <screen dir>/smoke_passed on success,
<screen dir>/smoke_FAILED.txt otherwise.

Called automatically by screen_isotropy.py after the plain-style twin of a smoke pair when
<screen dir>/gate_armed exists, or by hand:
    python screen_gate.py <omp_run_dir> <sfoff_run_dir>
"""
import glob
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

STAGES = ["ref", "ss_ramp_p", "ss_hold_p", "ss_ramp_m", "ss_hold_m", "ss_ramp_0", "ss_hold_0",
          "ps_ramp_p", "ps_hold_p", "ps_ramp_m", "ps_hold_m", "ps_ramp_0", "ps_hold_0"]
LARGEST_ATOMS = 353825
PROD_STEPS = 3.8e6
MAX_HOURS = 4.5


def check_run(run_dir, msgs):
    """Checks 1-6 on one run; returns (ok, info)."""
    import screen_isotropy as si
    import screen_voids as sv
    ok = True
    info = {}

    def bad(m):
        nonlocal ok
        ok = False
        msgs.append("  FAIL %s" % m)

    def good(m):
        msgs.append("  ok   %s" % m)

    msgs.append("run %s" % run_dir)
    log = open(os.path.join(run_dir, "log.lammps"), errors="replace").read()
    # 1
    if "SCREEN_COMPLETE" not in log:
        bad("log has no SCREEN_COMPLETE")
    nerr = len(re.findall(r"^ERROR", log, flags=re.M))
    nfene = log.count("FENE bond too long")
    nlost = len(re.findall(r"Lost atoms|Bond atoms .* missing|Out of range atoms", log))
    if nerr or nfene or nlost:
        bad("log: %d ERROR, %d FENE warnings, %d lost/missing-atom messages" % (nerr, nfene, nlost))
    else:
        good("log complete: no ERROR, no FENE warnings, no lost atoms")
    m = re.search(r"^\s*(\d+) atoms\s*$", log, flags=re.M)
    info["natoms"] = int(m.group(1)) if m else 0
    info["step0"] = [float(x) for ln in re.findall(r"^>>> STEP0 (.*)$", log, flags=re.M)
                     for x in ln.split()[1::2]]
    perf = [float(x) for x in re.findall(r"Performance:.*?([0-9.]+) timesteps/s", log)]
    info["steps_per_s"] = float(np.median(perf)) if perf else 0.0
    if not ok:
        return ok, info
    # 2-4 through the analysis module (raises on missing files)
    try:
        res = si.analyse(run_dir, log=lambda *_: None)
    except Exception as exc:                                    # noqa: BLE001
        bad("screen_isotropy.analyse raised %r" % (exc,))
        return ok, info
    st = res["stages"]
    missing = [s for s in STAGES if s not in st]
    if missing:
        bad("stages missing: %s" % missing)
        return ok, info
    last = st["ps_hold_0"]["last"]
    of = os.path.join(run_dir, "output_files")
    tag = res["tag"]
    P = np.loadtxt(os.path.join(of, "stress_data", "screen_stress_poly_%s.dat" % tag), ndmin=2)
    B = np.loadtxt(os.path.join(of, "stress_data", "screen_stress_box_%s.dat" % tag), ndmin=2)
    Dm = np.loadtxt(os.path.join(of, "volume_data", "screen_dims_%s.dat" % tag), ndmin=2)
    rec = int(P[1, 0] - P[0, 0])
    full = all(a[-1, 0] == last for a in (P, B, Dm)) and len(P) == last // rec and len(B) == last // rec
    finite = all(np.isfinite(a).all() for a in (P, B, Dm))
    if full and finite:
        good("13 stages; stress series complete (%d records of %d steps to step %d), finite" % (len(P), rec, last))
    else:
        bad("stress series incomplete or not finite (poly %d rows to %d, box %d rows to %d, expected step %d)"
            % (len(P), P[-1, 0], len(B), B[-1, 0], last))
    g = res["geometry"]
    sref = res["strain"]["ref"]
    worst = 0.0
    for s in res["holds_available"]:
        t = st[s]["target"]
        x = res["strain"][s]
        if s.startswith("ss") or s == "ref":
            dev = max(abs(x["gamma"] - (t if s.startswith("ss") else 0.0)), abs(x["e_x"]), abs(x["e_y"]))
        else:
            dev = max(abs(x["gamma"]), abs(x["e_x"] - t), abs(x["e_y"] - (1.0 / (1.0 + t) - 1.0)))
        dev = max(dev, abs(x["vol"] / sref["vol"] - 1.0))
        worst = max(worst, dev)
    mh = np.zeros(len(B), dtype=bool)
    for s in res["holds_available"]:
        mh |= (B[:, 0] > st[s]["first"]) & (B[:, 0] <= st[s]["last"])
    dlz = float(np.abs(B[:, 10] / g["lz0"] - 1.0).max())
    if len(res["holds_available"]) == 7 and worst < 1e-6 and dlz < 1e-9:
        good("7 holds at their target strains (worst deviation %.1e incl. volume), lz fixed (%.1e)" % (worst, dlz))
    else:
        bad("strain / volume off target: worst deviation %.2e over %d holds, lz change %.1e"
            % (worst, len(res["holds_available"]), dlz))
    th = res.get("thermo", {})
    if th and 0.97 < th["T_mobile_mean"] < 1.03 and 1.3 < th["press_mean"] < 1.7:
        good("T_mobile %.4f, <P> %.3f" % (th["T_mobile_mean"], th["press_mean"]))
    else:
        bad("thermo out of range: %s" % th)
    if g["wrapped_in_z"]:
        bad("gel straddles the z boundary (thickness %.1f vs lz %.1f)" % (g["thickness_polymer_bbox"], g["lz0"]))
    info["G"] = {k: res[k].get(q) for k in ("poly", "box") for q in ("G1",)} if "poly" in res else {}
    info["summary"] = "G1 poly %.4f box %.4f, G2 poly %.4f box %.4f (smoke values, not converged)" % (
        res["poly"].get("G1", float("nan")), res["box"].get("G1", float("nan")),
        res["poly"].get("G2", float("nan")), res["box"].get("G2", float("nan")))
    # 5 reference dump
    tr = glob.glob(os.path.join(run_dir, "traj_files", "screen_ref_polymer_*.lammpstrj"))
    offs, size = sv.frame_offsets(tr[0]) if tr else ([], 0)
    m = re.search(r"^dump ref_traj polymer custom (\d+)", log, flags=re.M)
    nexp = st["ref"]["last"] // int(m.group(1)) + 1 if m else -1
    if tr and (len(offs) == nexp or (nexp < 0 and len(offs) >= 2)):
        with open(tr[0], "rb") as fh:
            fr = sv.read_frame(fh, offs[-1], size, keep_types=(1, 2, 3, 4, 5))
        if fr is not None and set(np.unique(fr["type"]).tolist()) <= {1, 2}:
            good("reference dump: %d frames, %d polymer beads, last step %d" % (len(offs), len(fr["id"]), fr["step"]))
        else:
            bad("reference dump: last frame truncated or contains non-polymer types")
    else:
        bad("reference dump: %d frames, expected %d" % (len(offs), nexp))
    # 6 plates in the deformation dump
    td = glob.glob(os.path.join(run_dir, "traj_files", "screen_def_*.lammpstrj"))
    if not td:
        bad("no deformation dump")
        return ok, info
    offs, size = sv.frame_offsets(td[0])
    bounds = list(zip(offs, offs[1:] + [size]))
    s0 = None
    worst_s, worst_z, tilts = 0.0, 0.0, []
    with open(td[0], "rb") as fh:
        for a, b in bounds:
            fr = sv.read_frame(fh, a, b, keep_types=(4, 5))
            if fr is None:
                bad("deformation dump: truncated frame")
                return ok, info
            xy = fr["tilt"][0]
            tilts.append(xy / fr["L"][1])
            sy = (fr["x"][:, 1] - fr["lo"][1]) / fr["L"][1]
            sx = (fr["x"][:, 0] - fr["lo"][0] - xy * sy) / fr["L"][0]
            s = np.column_stack([sx, sy])
            for t in (4, 5):
                worst_z = max(worst_z, float(fr["x"][fr["type"] == t, 2].std()))
            if s0 is None:
                s0 = s
            else:
                d = s - s0
                d -= np.round(d)
                worst_s = max(worst_s, float(np.abs(d).max()))
    if worst_s < 2e-4 and worst_z < 1e-4 and max(tilts) > 0.9 * st["ss_hold_p"]["target"]:
        good("deformation dump: %d frames, plates planar (z spread %.1e) and affine (max fractional x,y change %.1e), tilt range %+.3f .. %+.3f"
             % (len(bounds), worst_z, worst_s, min(tilts), max(tilts)))
    else:
        bad("plates not affine/planar: fractional change %.2e, z spread %.2e, tilt range %+.3f .. %+.3f"
            % (worst_s, worst_z, min(tilts), max(tilts)))
    return ok, info


def find_twin(sfoff_dir):
    """The -sf omp run of the same job: newest run dir with the same name minus _sfoff, started earlier."""
    base = os.path.basename(sfoff_dir.rstrip("/"))
    m = re.match(r"(.*)_sfoff_(\d{8}_\d{6})$", base)
    if not m:
        return None
    cands = [d for d in sorted(glob.glob(os.path.join(os.path.dirname(sfoff_dir), m.group(1) + "_2*")))
             if re.match(re.escape(m.group(1)) + r"_\d{8}_\d{6}$", os.path.basename(d))
             and os.path.basename(d)[-15:] < m.group(2)]
    return cands[-1] if cands else None


def run(omp_dir, sfoff_dir, target_atoms=LARGEST_ATOMS, max_hours=MAX_HOURS, prod_steps=PROD_STEPS):
    """target_atoms = 0: the smoke run's own system (a production job's preflight on its own gel)."""
    screen = os.path.dirname(os.path.abspath(sfoff_dir.rstrip("/")))
    msgs = ["=== screen_gate %s (job %s) ===" % (time.strftime("%Y-%m-%d %H:%M:%S"), os.environ.get("SLURM_JOB_ID", "none"))]
    ok = True
    infos = []
    for d in (omp_dir, sfoff_dir):
        if d is None or not os.path.isdir(d):
            msgs.append("  FAIL run directory missing: %s" % d)
            ok = False
            infos.append({})
            continue
        try:
            o, info = check_run(os.path.abspath(d), msgs)
        except Exception as exc:                                # noqa: BLE001
            o, info = False, {}
            msgs.append("  FAIL exception %r" % (exc,))
        ok = ok and o
        infos.append(info)
    a, b = infos[0].get("step0", []), infos[1].get("step0", [])
    if len(a) == len(b) and len(a) >= 10:
        rel = max(abs(x - y) / max(abs(x), 1e-12) for x, y in zip(a, b))
        if rel < 1e-9:
            msgs.append("  ok   step-0 polymer-group stress and box pressure identical with -sf omp and plain styles (max rel. diff %.1e)" % rel)
        else:
            msgs.append("  FAIL step-0 stress differs between -sf omp and plain styles: max rel. diff %.3e" % rel)
            ok = False
    else:
        msgs.append("  FAIL STEP0 lines missing (%d vs %d values)" % (len(a), len(b)))
        ok = False
    rate, nat = infos[0].get("steps_per_s", 0.0), infos[0].get("natoms", 0)
    if rate > 0 and nat > 0:
        tgt = target_atoms if target_atoms > 0 else nat
        hours = prod_steps / (rate * nat / tgt) / 3600.0
        line = "%.0f steps/s for %d atoms with -sf omp -> %.2f h for %.2g steps of a %d-atom gel (limit %.2f h)" % (rate, nat, hours, prod_steps, tgt, max_hours)
        if hours <= max_hours:
            msgs.append("  ok   " + line)
        else:
            msgs.append("  FAIL " + line)
            ok = False
    else:
        msgs.append("  FAIL no Performance line in the omp log")
        ok = False
    for i in infos:
        if i.get("summary"):
            msgs.append("  info " + i["summary"])
    msgs.append("RESULT: %s" % ("PASS -> gate opened" if ok else "FAIL -> gate stays closed"))
    text = "\n".join(msgs) + "\n"
    with open(os.path.join(screen, "gate_check.log"), "a") as fh:
        fh.write(text)
    print(text, flush=True)
    target = os.path.join(screen, "smoke_passed" if ok else "smoke_FAILED.txt")
    with open(target, "a") as fh:
        fh.write(text)
    return ok


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("omp_dir")
    ap.add_argument("sfoff_dir", nargs="?")
    ap.add_argument("--target-atoms", type=int, default=LARGEST_ATOMS, help="0 = the smoke run's own atom count")
    ap.add_argument("--max-hours", type=float, default=MAX_HOURS)
    ap.add_argument("--prod-steps", type=float, default=PROD_STEPS)
    a = ap.parse_args(argv)
    omp, off = (a.omp_dir, a.sfoff_dir) if a.sfoff_dir else (find_twin(a.omp_dir), a.omp_dir)
    return 0 if run(omp, off, a.target_atoms, a.max_hours, a.prod_steps) else 1


if __name__ == "__main__":
    sys.exit(main())
