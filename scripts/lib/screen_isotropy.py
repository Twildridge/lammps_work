#!/usr/bin/env python
"""screen_isotropy.py -- G1, G2, r = G1/G2 and the Zener ratio from a screen_isotropy.lmp run (2026-11-13).

Part B of the random-network screen.  Reads, in <run_dir>/output_files/:
    screen_stages_<tag>.dat                 stage, first step, last step, kind, target strain
    stress_data/screen_stress_poly_<tag>.dat   step, S_ab = sum of compute stress/atom over the polymer
    stress_data/screen_stress_box_<tag>.dat    step, P_ab (whole box), vol, lx, ly, lz, xy
    volume_data/screen_dims_<tag>.dat          step, polymer z min/max, crosslinker z min/max

Sign convention.  LAMMPS: stress/atom = -(pressure x volume), P_ab = -stress_ab.  Everything below
is tension-positive STRESS of the gel, per gel volume V_gel = lx0 * ly0 * <thickness>, where the
thickness is the polymer bounding-box z extent averaged over the reference stage (the deck's
gel_dimensions convention; the crosslinker extent is stored next to it):
    polymer partial stress   sigma_ab = +S_ab / V_gel
    whole box                sigma_ab = -P_ab * V_box / V_gel

    G1 = [sigma_xy(+g) - sigma_xy(-g)] / (2 g)                       simple shear, g = xy / ly
    G2 = [(sigma_xx - sigma_yy)(+e) - (sigma_xx - sigma_yy)(-e)] / (8 e)   pure shear at fixed area
    r  = G1 / G2   (1 = isotropic);   A = 3 r / (4 - r)   (Zener ratio of a cubic lattice in the
    rot45 orientation: G1 = C44, G2 = (3 C' + C44) / 4)

Each hold is averaged over its last avg_frac (0.75: the last 300k of 400k steps); the error bar is
the standard error of the means of blocks of block_steps (50k) steps.

Usage:  python screen_isotropy.py <run_dir> [--label L] [--block-steps 50000] [--avg-frac 0.75]
Writes <run_dir>/screen_isotropy_<label>.json and prints a summary.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np


def _load(path):
    return np.loadtxt(path, comments="#", ndmin=2)


def read_stages(path):
    st = {}
    order = []
    with open(path) as fh:
        for ln in fh:
            w = ln.split()
            if not w or w[0].startswith("#"):
                continue
            st[w[0]] = dict(first=int(float(w[1])), last=int(float(w[2])), kind=w[3], target=float(w[4]))
            order.append(w[0])
    return st, order


def block_stats(step, y, first, last, avg_frac, block_steps):
    """Mean and block standard error of y over the last avg_frac of (first, last]."""
    t0 = last - avg_frac * (last - first)
    m = (step > t0 + 1e-9) & (step <= last)
    s, v = step[m], y[m]
    if len(v) == 0:
        return None
    nb = int(round((last - t0) / block_steps))
    nb = max(nb, 1)
    edges = t0 + (last - t0) * np.arange(nb + 1) / nb
    means = []
    for a, b in zip(edges[:-1], edges[1:]):
        mm = (s > a + 1e-9) & (s <= b + 1e-9)
        if mm.any():
            means.append(v[mm].mean())
    means = np.asarray(means)
    sem = float(means.std(ddof=1) / np.sqrt(len(means))) if len(means) > 1 else float("nan")
    half = len(v) // 2
    return dict(mean=float(v.mean()), sem=sem, n_blocks=int(len(means)), n_records=int(len(v)),
                block_means=[float(x) for x in means],
                half_diff=float(v[half:].mean() - v[:half].mean()) if half > 0 else 0.0,
                first=float(t0), last=float(last))


def ratio(a, da, b, db):
    r = a / b
    dr = abs(r) * np.sqrt((da / a) ** 2 + (db / b) ** 2)
    return float(r), float(dr)


def thermo_means(log_path):
    """Mean c_mobile_temp and press over the thermo rows of log.lammps (runs of >= 3 rows)."""
    if not os.path.isfile(log_path):
        return {}
    cols, rows, out_t, out_p = None, [], [], []
    with open(log_path, errors="replace") as fh:
        for ln in fh:
            w = ln.split()
            if w and w[0] == "Step":
                cols, rows = w, []
                continue
            if cols is None:
                continue
            if w and w[0] == "Loop":
                if len(rows) >= 3 and "c_mobile_temp" in cols and "Press" in cols:
                    a = np.asarray(rows)
                    out_t.extend(a[:, cols.index("c_mobile_temp")])
                    out_p.extend(a[:, cols.index("Press")])
                cols = None
                continue
            try:
                vals = [float(x) for x in w]
            except ValueError:
                continue
            if len(vals) == len(cols):
                rows.append(vals)
    if not out_t:
        return {}
    return dict(T_mobile_mean=float(np.mean(out_t)), T_mobile_std=float(np.std(out_t)),
                press_mean=float(np.mean(out_p)), n_rows=len(out_t))


def analyse(run_dir, label=None, block_steps=50000, avg_frac=0.75, log=print):
    of = os.path.join(run_dir, "output_files")
    stage_files = sorted(glob.glob(os.path.join(of, "screen_stages_*.dat")))
    if not stage_files:
        raise RuntimeError("no screen_stages_*.dat in %s" % of)
    tag = os.path.basename(stage_files[0])[len("screen_stages_"):-len(".dat")]
    label = label or tag
    stages, order = read_stages(stage_files[0])
    P = _load(os.path.join(of, "stress_data", "screen_stress_poly_%s.dat" % tag))
    B = _load(os.path.join(of, "stress_data", "screen_stress_box_%s.dat" % tag))
    Dm = _load(os.path.join(of, "volume_data", "screen_dims_%s.dat" % tag))
    last_step = min(P[-1, 0], B[-1, 0])
    res = dict(label=label, tag=tag, run_dir=run_dir, block_steps=block_steps, avg_frac=avg_frac,
               last_step_on_disk=int(last_step), stages=stages, stage_order=order)
    res["complete"] = bool(os.path.isfile(os.path.join(run_dir, "log.lammps")) and
                           "SCREEN_COMPLETE" in open(os.path.join(run_dir, "log.lammps"), errors="replace").read())
    res["thermo"] = thermo_means(os.path.join(run_dir, "log.lammps"))

    # reference geometry
    ref = stages["ref"]
    mref = (B[:, 0] > ref["first"]) & (B[:, 0] <= ref["last"])
    lx0, ly0, lz0 = B[mref, 8].mean(), B[mref, 9].mean(), B[mref, 10].mean()
    dref = (Dm[:, 0] > ref["first"]) & (Dm[:, 0] <= ref["last"])
    thick = float((Dm[dref, 2] - Dm[dref, 1]).mean())
    thick_xl = float((Dm[dref, 4] - Dm[dref, 3]).mean())
    vgel = lx0 * ly0 * thick
    res["geometry"] = dict(lx0=float(lx0), ly0=float(ly0), lz0=float(lz0), thickness_polymer_bbox=thick,
                           thickness_xlink=thick_xl, thickness_std=float((Dm[dref, 2] - Dm[dref, 1]).std()),
                           V_gel=float(vgel), V_box=float(B[mref, 7].mean()),
                           wrapped_in_z=bool(thick > lz0 - 2.0))
    if res["geometry"]["wrapped_in_z"]:
        log("  WARNING: polymer z extent %.2f ~ lz %.2f: the gel straddles the z boundary, V_gel is wrong" % (thick, lz0))

    series = {
        "poly": dict(step=P[:, 0], sxy=P[:, 4] / vgel, nd=(P[:, 1] - P[:, 2]) / vgel,
                     sxx=P[:, 1] / vgel, syy=P[:, 2] / vgel, szz=P[:, 3] / vgel),
        "box": dict(step=B[:, 0], sxy=-B[:, 4] * B[:, 7] / vgel, nd=-(B[:, 1] - B[:, 2]) * B[:, 7] / vgel,
                    sxx=-B[:, 1] * B[:, 7] / vgel, syy=-B[:, 2] * B[:, 7] / vgel, szz=-B[:, 3] * B[:, 7] / vgel),
    }
    holds = [s for s in order if stages[s]["kind"] == "hold" and stages[s]["last"] <= last_step]
    res["holds_available"] = holds

    # actual strains of each hold
    strain = {}
    for s in holds:
        m = (B[:, 0] > stages[s]["first"]) & (B[:, 0] <= stages[s]["last"])
        strain[s] = dict(gamma=float((B[m, 11] / B[m, 9]).mean()), e_x=float(B[m, 8].mean() / lx0 - 1.0),
                         e_y=float(B[m, 9].mean() / ly0 - 1.0), vol=float(B[m, 7].mean()))
    res["strain"] = strain

    for name, S in series.items():
        out = dict(holds={})
        for s in holds:
            st = stages[s]
            bs = block_steps if (st["last"] - st["first"]) * avg_frac >= 2 * block_steps else \
                max(1.0, (st["last"] - st["first"]) * avg_frac / 3.0)
            out["holds"][s] = {q: block_stats(S["step"], S[q], st["first"], st["last"], avg_frac, bs)
                               for q in ("sxy", "nd", "sxx", "syy", "szz")}
        H = out["holds"]
        if "ss_hold_p" in H and "ss_hold_m" in H:
            gp, gm = strain["ss_hold_p"]["gamma"], strain["ss_hold_m"]["gamma"]
            a, b = H["ss_hold_p"]["sxy"], H["ss_hold_m"]["sxy"]
            out["G1"] = float((a["mean"] - b["mean"]) / (gp - gm))
            out["G1_err"] = float(np.hypot(a["sem"], b["sem"]) / abs(gp - gm))
            z = H["ref"]["sxy"]
            out["G1_plus"] = float((a["mean"] - z["mean"]) / gp)
            out["G1_minus"] = float((z["mean"] - b["mean"]) / (-gm))
        if "ps_hold_p" in H and "ps_hold_m" in H:
            ep, em = strain["ps_hold_p"]["e_x"], strain["ps_hold_m"]["e_x"]
            a, b = H["ps_hold_p"]["nd"], H["ps_hold_m"]["nd"]
            out["G2"] = float((a["mean"] - b["mean"]) / (4.0 * (ep - em)))
            out["G2_err"] = float(np.hypot(a["sem"], b["sem"]) / (4.0 * abs(ep - em)))
            z = H["ref"]["nd"]
            out["G2_plus"] = float((a["mean"] - z["mean"]) / (4.0 * ep))
            out["G2_minus"] = float((z["mean"] - b["mean"]) / (-4.0 * em))
        if "G1" in out and "G2" in out:
            r, dr = ratio(out["G1"], out["G1_err"], out["G2"], out["G2_err"])
            out["r"], out["r_err"] = r, dr
            out["r_rel_err"] = float(dr / abs(r))
            out["A"] = float(3.0 * r / (4.0 - r))
            out["A_err"] = float(12.0 / (4.0 - r) ** 2 * dr)
        res[name] = out

    with open(os.path.join(run_dir, "screen_isotropy_%s.json" % label), "w") as fh:
        json.dump(res, fh, indent=1)

    g = res["geometry"]
    log("[%s] %s  complete=%s  last step %d" % (label, run_dir, res["complete"], last_step))
    log("  box %.3f x %.3f x %.3f, gel thickness %.2f (crosslinkers %.2f), V_gel %.0f, V_box %.0f; T_mobile %s"
        % (g["lx0"], g["ly0"], g["lz0"], g["thickness_polymer_bbox"], g["thickness_xlink"], g["V_gel"], g["V_box"],
           "%.4f" % res["thermo"]["T_mobile_mean"] if res["thermo"] else "n/a"))
    for name in ("poly", "box"):
        out = res[name]
        log("  --- %s stress ---" % ("polymer-group" if name == "poly" else "whole-box"))
        for s in holds:
            h = out["holds"][s]
            log("    %-10s gamma %+8.5f e_x %+8.5f | sigma_xy %+.3e +- %.1e | sxx-syy %+.3e +- %.1e (%d blocks)"
                % (s, strain[s]["gamma"], strain[s]["e_x"], h["sxy"]["mean"], h["sxy"]["sem"], h["nd"]["mean"],
                   h["nd"]["sem"], h["nd"]["n_blocks"]))
        if "G1" in out:
            log("    G1 = %.4f +- %.4f   (one-sided: + %.4f, - %.4f)" % (out["G1"], out["G1_err"], out["G1_plus"], out["G1_minus"]))
        if "G2" in out:
            log("    G2 = %.4f +- %.4f   (one-sided: + %.4f, - %.4f)" % (out["G2"], out["G2_err"], out["G2_plus"], out["G2_minus"]))
        if "r" in out:
            log("    r = G1/G2 = %.3f +- %.3f (%.0f %%)   A = %.3f +- %.3f" % (out["r"], out["r_err"], 100 * out["r_rel_err"], out["A"], out["A_err"]))
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--label")
    ap.add_argument("--block-steps", type=float, default=50000)
    ap.add_argument("--avg-frac", type=float, default=0.75)
    a = ap.parse_args(argv)
    analyse(os.path.abspath(a.run_dir), a.label, a.block_steps, a.avg_frac)
    return 0


if __name__ == "__main__":
    sys.exit(main())
