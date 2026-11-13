#!/usr/bin/env python
"""screen_voids.py -- void persistence Phi and network statics for the gel screen (2026-11-13).

Part A of the random-network screen (report: lammps_runs/chain_length_testing/screen_report.md).
Nothing here runs MD; it reads an existing polymer trajectory (.lammpstrj) and, optionally, the
gel's final configuration (.data, for the bond list).

Void persistence
----------------
For the last N frames of a trajectory (polymer beads = types 1, 2; type 1 = crosslinker):

1. Drift.  x, y: every frame is put in fractional coordinates of ITS OWN box; the polymer
   centre-of-mass shift between consecutive frames is the mean minimum-image fractional
   displacement of the beads (matched by id), accumulated and subtracted.  z: the gel does not
   span z, so the polymer COM z of each frame is the origin (the gel is first un-wrapped across
   the z boundary if it straddles it).
2. Region and grid.  Full periodic x, y cross-section; in z the middle half of the gel, from the
   crosslinker z extent (window mean of min/max z of type 1, measured from the COM).  A grid of
   spacing h = 0.5 sigma: the SAME number of points in every frame, at fixed fractional x, y
   (so the grid deforms affinely with the box under the barostat) and, in z, at fixed offsets
   from the COM scaled by Rg_z(frame) / <Rg_z> (affine in z as well; identity at fixed box).
3. A grid point is empty in a frame if no bead centre lies within the probe radius r_p
   (nearest-bead distance > r_p; cKDTree with explicit periodic images in x, y).
4. emp[frame, point] boolean, p = emp.mean(), pe = emp.mean(axis=0), nt frames:
       Phi = (pe.var() / (p (1 - p)) - 1 / nt) / (1 - 1 / nt)
   0 = every point empty equally often, 1 = frozen pores.  Skipped (nan) when p < 0.005 or
   p > 0.995.  Phi is given for the whole window and for the last w frames (w = 5, 10, ...), and
   the two-time overlap C(lag) = (<emp(t) emp(t+lag)> - p^2) / (p (1 - p)) for a set of lags:
   both decay towards the frozen fraction, a plateau is the frozen part.

Network statics (--data)
------------------------
From the final configuration and its bond list: strands = paths from crosslinker to crosslinker
through chain beads; R = end-to-end vector (sum of minimum-image bond vectors), n = bonds in the
strand, b = mean bond length; chord extension = |R| / (n b).  "Interior" = middle half of the
crosslinker z extent (strand midpoint for strands, bead position for the number density).

Usage
-----
  python screen_voids.py --traj T.lammpstrj --label NAME --out NAME.json [--last 100]
         [--data final_config.data] [--probes 1.0 1.5 2.5] [--rel-probe 0.4] [--h 0.5]
"""
import argparse
import io
import json
import mmap
import os
import sys
import time

import numpy as np

POLYMER_TYPES = (1, 2)
XLINK_TYPE = 1


# ------------------------------------------------------------------------------------------------
# trajectory reading
# ------------------------------------------------------------------------------------------------
def frame_offsets(path):
    """Byte offset of every 'ITEM: TIMESTEP' record, and the file size."""
    key = b"ITEM: TIMESTEP"
    offs = []
    with open(path, "rb") as fh:
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        pos = 0
        while True:
            i = mm.find(key, pos)
            if i < 0:
                break
            offs.append(i)
            pos = i + len(key)
        size = mm.size()
        mm.close()
    return offs, size


def read_frame(fh, a, b, keep_types=POLYMER_TYPES):
    """One frame from bytes [a, b): dict(step, lo, L, tilt, id, type, x) with the beads of
    keep_types sorted by id.  Returns None for a truncated frame."""
    import pandas as pd
    fh.seek(a)
    buf = fh.read(b - a)
    head = []
    pos = 0
    for _ in range(9):
        j = buf.find(b"\n", pos)
        if j < 0:
            return None
        head.append(buf[pos:j].decode())
        pos = j + 1
    step = int(head[1])
    natoms = int(head[3])
    tri = "xy" in head[4]
    bnd = np.array([[float(v) for v in head[5 + k].split()] for k in range(3)])
    if tri:
        xy, xz, yz = bnd[0, 2], bnd[1, 2], bnd[2, 2]
        lo = np.array([bnd[0, 0] - min(0.0, xy, xz, xy + xz), bnd[1, 0] - min(0.0, yz), bnd[2, 0]])
        hi = np.array([bnd[0, 1] - max(0.0, xy, xz, xy + xz), bnd[1, 1] - max(0.0, yz), bnd[2, 1]])
        tilt = (xy, xz, yz)
    else:
        lo, hi = bnd[:, 0].copy(), bnd[:, 1].copy()
        tilt = (0.0, 0.0, 0.0)
    cols = head[8].split()[2:]
    df = pd.read_csv(io.BytesIO(buf[pos:]), delim_whitespace=True, header=None, names=cols)
    if len(df) < natoms or df[cols[-1]].isna().any():
        return None
    df = df.iloc[:natoms]
    typ = df["type"].to_numpy().astype(np.int32)
    m = np.isin(typ, keep_types)
    ids = df["id"].to_numpy().astype(np.int64)[m]
    order = np.argsort(ids, kind="stable")
    x = df[["x", "y", "z"]].to_numpy(dtype=np.float64)[m][order]
    return dict(step=step, lo=lo, L=hi - lo, tilt=tilt, id=ids[order], type=typ[m][order], x=x)


def read_last_frames(path, n_last, stride=1, log=print):
    """The last n_last complete frames (every stride-th, counted back from the end)."""
    offs, size = frame_offsets(path)
    bounds = list(zip(offs, offs[1:] + [size]))
    frames = []
    with open(path, "rb") as fh:
        # the last record may be truncated: test it first
        last = read_frame(fh, *bounds[-1])
        if last is None:
            log("  last frame of %s is truncated: dropped" % os.path.basename(path))
            bounds = bounds[:-1]
        pick = bounds[::-1][::stride][:n_last][::-1]
        for a, b in pick:
            fr = read_frame(fh, a, b)
            if fr is None:
                raise RuntimeError("truncated frame inside %s at byte %d" % (path, a))
            frames.append(fr)
    return frames, len(bounds)


def unwrap_z(z, lo, Lz):
    """Un-wrap a slab that may straddle the periodic z boundary: cut at the largest empty gap."""
    zs = np.sort(z)
    gaps = np.diff(zs)
    wrap_gap = zs[0] + Lz - zs[-1]
    k = int(np.argmax(gaps)) if len(gaps) else 0
    if len(gaps) and gaps[k] > wrap_gap:
        cut = 0.5 * (zs[k] + zs[k + 1])
        z = np.where(z < cut, z + Lz, z)
    return z


# ------------------------------------------------------------------------------------------------
# the metric
# ------------------------------------------------------------------------------------------------
def phi_stat(emp):
    """(p, Phi) of a boolean [frame, point] array; Phi = nan outside 0.005 <= p <= 0.995."""
    nt = emp.shape[0]
    p = float(emp.mean())
    if nt < 2 or p < 0.005 or p > 0.995:
        return p, float("nan")
    pe = emp.mean(axis=0, dtype=np.float64)
    return p, float((pe.var() / (p * (1.0 - p)) - 1.0 / nt) / (1.0 - 1.0 / nt))


def nearest_distance_field(frames, h=0.5, rmax=3.0, workers=-1, log=print):
    """Nearest-bead distance on the co-moving grid of every frame.

    Returns D [nt, npts] float32 (capped at rmax + 1 where nothing is within rmax) and a dict of
    the per-frame geometry (box, gel extent, drift)."""
    from scipy.spatial import cKDTree
    nt = len(frames)
    ids0 = frames[0]["id"]
    S = np.zeros(2)
    prev = None
    geo = dict(step=[], lx=[], ly=[], lz=[], zc=[], rgz=[], xl_lo=[], xl_hi=[], pol_lo=[], pol_hi=[],
               drift_x=[], drift_y=[])
    sxy = []
    zun = []
    for fr in frames:
        if any(abs(t) > 1e-12 for t in fr["tilt"]):
            raise RuntimeError("void metric expects an orthogonal box (frame %d is tilted)" % fr["step"])
        if len(fr["id"]) != len(ids0) or not np.array_equal(fr["id"], ids0):
            raise RuntimeError("polymer ids differ between frames (step %d)" % fr["step"])
        L, lo = fr["L"], fr["lo"]
        s = (fr["x"][:, :2] - lo[:2]) / L[:2]
        if prev is not None:
            ds = s - prev
            ds -= np.round(ds)
            S += ds.mean(axis=0)
        prev = s
        sxy.append((s - S) % 1.0)
        z = unwrap_z(fr["x"][:, 2], lo[2], L[2])
        zun.append(z)
        xl = z[fr["type"] == XLINK_TYPE]
        zc = z.mean()
        geo["step"].append(int(fr["step"]))
        geo["lx"].append(float(L[0])); geo["ly"].append(float(L[1])); geo["lz"].append(float(L[2]))
        geo["zc"].append(float(zc)); geo["rgz"].append(float(z.std()))
        geo["xl_lo"].append(float(xl.min() - zc)); geo["xl_hi"].append(float(xl.max() - zc))
        geo["pol_lo"].append(float(z.min() - zc)); geo["pol_hi"].append(float(z.max() - zc))
        geo["drift_x"].append(float(S[0] * L[0])); geo["drift_y"].append(float(S[1] * L[1]))
    g = {k: np.asarray(v) for k, v in geo.items()}
    lo_off, hi_off = g["xl_lo"].mean(), g["xl_hi"].mean()
    T = hi_off - lo_off
    z0, z1 = lo_off + 0.25 * T, hi_off - 0.25 * T
    nx = max(1, int(round(g["lx"].mean() / h)))
    ny = max(1, int(round(g["ly"].mean() / h)))
    nz = max(1, int(np.floor((z1 - z0) / h)))
    fx = (np.arange(nx) + 0.5) / nx
    fy = (np.arange(ny) + 0.5) / ny
    zrel = z0 + (np.arange(nz) + 0.5) * h
    rgz_mean = g["rgz"].mean()
    npts = nx * ny * nz
    log("  grid %d x %d x %d = %d points; z window [%.2f, %.2f] about the COM (crosslinker extent %.2f)"
        % (nx, ny, nz, npts, z0, z1, T))
    D = np.empty((nt, npts), dtype=np.float32)
    cap = rmax + 1.0
    for t, fr in enumerate(frames):
        L = fr["L"]
        scale = g["rgz"][t] / rgz_mean
        P = np.column_stack([sxy[t][:, 0] * L[0], sxy[t][:, 1] * L[1], zun[t] - g["zc"][t]])
        zg = zrel * scale
        P = P[(P[:, 2] > zg[0] - rmax - 1.0) & (P[:, 2] < zg[-1] + rmax + 1.0)]
        # explicit periodic images in x, then y (covers the corners)
        for ax in (0, 1):
            lo_side = P[P[:, ax] < rmax + 0.5].copy()
            lo_side[:, ax] += L[ax]
            hi_side = P[P[:, ax] > L[ax] - rmax - 0.5].copy()
            hi_side[:, ax] -= L[ax]
            P = np.vstack([P, lo_side, hi_side])
        X, Y, Z = np.meshgrid(fx * L[0], fy * L[1], zg, indexing="ij")
        G = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
        d, _ = cKDTree(P).query(G, k=1, distance_upper_bound=rmax, workers=workers)
        d[~np.isfinite(d)] = cap
        D[t] = d
    out = {k: v.tolist() for k, v in g.items()}
    out.update(grid=[nx, ny, nz], z_window=[float(z0), float(z1)], xl_extent=float(T), h=h)
    return D, out


def void_persistence(frames, probes, h=0.5, windows=(5, 10, 20, 30, 50, 100),
                     lags=(1, 2, 3, 5, 10, 20, 30, 50, 75, 99), workers=-1, log=print):
    nt = len(frames)
    rmax = max(probes) + 0.01
    D, geo = nearest_distance_field(frames, h=h, rmax=rmax, workers=workers, log=log)
    steps = np.asarray(geo["step"])
    dstep = int(np.median(np.diff(steps))) if nt > 1 else 0
    res = dict(n_frames=nt, step_first=int(steps[0]), step_last=int(steps[-1]), frame_spacing=dstep,
               probes=[], geometry=geo)
    wins = sorted(set([w for w in windows if 2 <= w < nt] + [nt]))
    for rp in probes:
        emp = D > rp
        p, phi = phi_stat(emp)
        row = dict(r_p=float(rp), empty_fraction=p, phi=phi, windows=[], lags=[])
        for w in wins:
            pw, phw = phi_stat(emp[-w:])
            row["windows"].append(dict(frames=int(w), span_steps=int(steps[-1] - steps[-w]),
                                       empty_fraction=pw, phi=phw))
        if 0.005 <= p <= 0.995:
            for k in lags:
                if k >= nt:
                    continue
                acc = 0.0
                for t in range(nt - k):
                    acc += np.count_nonzero(emp[t] & emp[t + k])
                c = (acc / ((nt - k) * emp.shape[1]) - p * p) / (p * (1.0 - p))
                row["lags"].append(dict(lag_frames=int(k), lag_steps=int(k * dstep), overlap=float(c)))
        res["probes"].append(row)
        log("  r_p %.3f: empty fraction %.4f, Phi %s" % (rp, p, "%.4f" % phi if phi == phi else "skipped"))
    return res


def drift_summary(geo):
    """Box and gel-thickness drift over the frames used (first vs last tenth, and a linear fit)."""
    steps = np.asarray(geo["step"], dtype=float)
    out = {}
    n = len(steps)
    k = max(1, n // 10)
    series = dict(lx=np.asarray(geo["lx"]), ly=np.asarray(geo["ly"]),
                  thick_xl=np.asarray(geo["xl_hi"]) - np.asarray(geo["xl_lo"]),
                  thick_pol=np.asarray(geo["pol_hi"]) - np.asarray(geo["pol_lo"]),
                  thick_rg=np.sqrt(12.0) * np.asarray(geo["rgz"]))
    for name, v in series.items():
        first, last = float(v[:k].mean()), float(v[-k:].mean())
        slope = float(np.polyfit(steps, v, 1)[0]) if n > 2 and steps[-1] > steps[0] else 0.0
        out[name] = dict(mean=float(v.mean()), first=first, last=last, change=last - first,
                         change_pct=100.0 * (last - first) / first, std=float(v.std()),
                         fit_change=slope * (steps[-1] - steps[0]),
                         fit_change_pct=100.0 * slope * (steps[-1] - steps[0]) / float(v.mean()))
    out["com_drift_removed"] = dict(x=float(geo["drift_x"][-1]), y=float(geo["drift_y"][-1]))
    return out


# ------------------------------------------------------------------------------------------------
# network statics from a write_data file
# ------------------------------------------------------------------------------------------------
def read_data(path):
    """Minimal reader of a LAMMPS write_data file (atom_style molecular): box, atoms, bonds."""
    import pandas as pd
    with open(path, "rb") as fh:
        raw = fh.read()
    lines_end = raw.find(b"\nAtoms")
    head = raw[:lines_end].decode().splitlines()
    box = np.zeros((3, 2))
    step = None
    natoms = nbonds = 0
    tilt = (0.0, 0.0, 0.0)
    if "timestep" in head[0]:
        try:
            step = int(head[0].split("timestep =")[1].split(",")[0])
        except (IndexError, ValueError):
            step = None
    for ln in head:
        w = ln.split()
        if len(w) == 2 and w[1] == "atoms":
            natoms = int(w[0])
        elif len(w) == 2 and w[1] == "bonds":
            nbonds = int(w[0])
        elif len(w) == 4 and w[2] in ("xlo", "ylo", "zlo"):
            box["xyz".index(w[2][0])] = [float(w[0]), float(w[1])]
        elif len(w) == 6 and w[3] == "xy":
            tilt = (float(w[0]), float(w[1]), float(w[2]))

    def section(name, nrows):
        i = raw.find(b"\n" + name.encode())
        i = raw.find(b"\n", i + 1) + 1          # end of the header line
        i = raw.find(b"\n", i) + 1              # the blank line
        return pd.read_csv(io.BytesIO(raw[i:]), delim_whitespace=True, header=None, nrows=nrows,
                           comment="#")

    at = section("Atoms", natoms)
    bd = section("Bonds", nbonds)
    ids = at[0].to_numpy().astype(np.int64)
    typ = at[2].to_numpy().astype(np.int32)
    x = at[[3, 4, 5]].to_numpy(dtype=np.float64)
    bonds = bd[[2, 3]].to_numpy().astype(np.int64)
    return dict(step=step, lo=box[:, 0], L=box[:, 1] - box[:, 0], tilt=tilt, id=ids, type=typ, x=x,
                bonds=bonds, natoms=natoms)


def network_statics(path, log=print):
    d = read_data(path)
    L, lo = d["L"], d["lo"]
    if any(abs(t) > 1e-12 for t in d["tilt"]):
        raise RuntimeError("network_statics expects an orthogonal box")
    idx = np.full(int(d["id"].max()) + 1, -1, dtype=np.int64)
    idx[d["id"]] = np.arange(len(d["id"]))
    typ, x = d["type"], d["x"]
    pol = np.isin(typ, POLYMER_TYPES)
    bi, bj = idx[d["bonds"][:, 0]], idx[d["bonds"][:, 1]]
    dv = x[bj] - x[bi]
    dv -= L * np.round(dv / L)
    blen = np.linalg.norm(dv, axis=1)
    b = float(blen.mean())
    # adjacency with the bond vector pointing away from the atom
    nbr = {}
    for a, c, v in zip(bi.tolist(), bj.tolist(), dv):
        nbr.setdefault(a, []).append((c, v))
        nbr.setdefault(c, []).append((a, -v))
    z = x[:, 2].copy()
    zpol = unwrap_z(z[pol], lo[2], L[2])
    z[pol] = zpol
    xl = np.where(typ == XLINK_TYPE)[0]
    zx = z[xl]
    T = float(zx.max() - zx.min())
    zin0, zin1 = zx.min() + 0.25 * T, zx.max() - 0.25 * T
    n_int = int(np.count_nonzero((zpol > zin0) & (zpol < zin1)))
    rho_int = n_int / (L[0] * L[1] * (zin1 - zin0))
    rho_all = int(pol.sum()) / (L[0] * L[1] * float(zpol.max() - zpol.min()))
    deg = np.array([len(nbr.get(int(a), [])) for a in xl])
    seen = set()
    R, nb, zmid, dang, loops = [], [], [], 0, 0
    for a in xl.tolist():
        for (j, v) in nbr.get(a, []):
            prev, cur = a, j
            vec = v.copy()
            n = 1
            first_mid = last_mid = None
            ok = True
            while typ[cur] != XLINK_TYPE:
                if first_mid is None:
                    first_mid = cur
                last_mid = cur
                nxt = [(c, w) for (c, w) in nbr[cur] if c != prev]
                if len(nxt) != 1:
                    ok = False           # dangling end (0) or a branch on a chain bead (>1)
                    break
                prev, (cur, w) = cur, nxt[0]
                vec = vec + w
                n += 1
            if not ok:
                dang += 1
                continue
            key = (min(first_mid, last_mid), max(first_mid, last_mid)) if first_mid is not None \
                else ("x", min(a, cur), max(a, cur))
            if key in seen:
                continue
            seen.add(key)
            if cur == a:
                loops += 1
            R.append(vec)
            nb.append(n)
            zmid.append(z[a] + 0.5 * vec[2])
    R = np.asarray(R)
    nb = np.asarray(nb)
    zmid = np.asarray(zmid)
    Rn = np.linalg.norm(R, axis=1)
    ext = Rn / (nb * b)
    inner = (zmid > zin0) & (zmid < zin1)

    def stats(m):
        e = ext[m]
        return dict(n_strands=int(m.sum()), R_mean=float(Rn[m].mean()), R_rms=float(np.sqrt((Rn[m] ** 2).mean())),
                    R_std=float(Rn[m].std()), ext_mean=float(e.mean()), ext_cv=float(e.std() / e.mean()),
                    R2_xyz=[float(v) for v in (R[m] ** 2).mean(axis=0)])

    vals, counts = np.unique(nb, return_counts=True)
    dvals, dcounts = np.unique(deg, return_counts=True)
    out = dict(data_file=path, step=d["step"], natoms=int(d["natoms"]), n_polymer=int(pol.sum()),
               n_crosslinkers=int(len(xl)), box=[float(v) for v in L],
               bond_length_mean=b, bond_length_max=float(blen.max()),
               bonds_per_strand={int(k): int(c) for k, c in zip(vals, counts)},
               bonds_per_strand_mode=int(vals[np.argmax(counts)]),
               crosslinker_degree={int(k): int(c) for k, c in zip(dvals, dcounts)},
               dangling_walks=int(dang), self_loops=int(loops),
               xl_extent=T, polymer_extent=float(zpol.max() - zpol.min()),
               rho_polymer_interior=float(rho_int), rho_polymer_bbox=float(rho_all),
               interior=stats(inner), all=stats(np.ones(len(nb), dtype=bool)))
    log("  statics: %d strands (%d interior), n = %d bonds, b = %.4f, <R> = %.3f, ext = %.3f (CV %.3f), rho_int = %.4f"
        % (len(nb), int(inner.sum()), out["bonds_per_strand_mode"], b, out["interior"]["R_mean"],
           out["interior"]["ext_mean"], out["interior"]["ext_cv"], rho_int))
    return out


# ------------------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--traj", help="polymer trajectory (.lammpstrj)")
    ap.add_argument("--data", help="final configuration (.data) for the network statics")
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--last", type=int, default=100, help="number of frames from the end")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--probes", type=float, nargs="*", default=[1.0, 1.5, 2.5])
    ap.add_argument("--rel-probe", type=float, default=0.4, help="extra probe = this x <R_ee> (needs --data)")
    ap.add_argument("--h", type=float, default=0.5)
    ap.add_argument("--workers", type=int, default=-1)
    a = ap.parse_args(argv)

    def log(msg):
        print(msg, flush=True)

    t0 = time.time()
    res = dict(label=a.label)
    probes = list(a.probes)
    if a.data:
        log("[%s] statics from %s" % (a.label, a.data))
        res["statics"] = network_statics(a.data, log=log)
        if a.rel_probe > 0:
            rp = a.rel_probe * res["statics"]["interior"]["R_mean"]
            res["rel_probe"] = dict(factor=a.rel_probe, r_p=rp)
            probes.append(rp)
    if a.traj:
        log("[%s] trajectory %s" % (a.label, a.traj))
        frames, ntot = read_last_frames(a.traj, a.last, a.stride, log=log)
        log("  %d of %d frames read, steps %d .. %d (%.0f s)"
            % (len(frames), ntot, frames[0]["step"], frames[-1]["step"], time.time() - t0))
        res["traj"] = a.traj
        res["frames_in_file"] = ntot
        res["voids"] = void_persistence(frames, probes, h=a.h, workers=a.workers, log=log)
        res["drift"] = drift_summary(res["voids"]["geometry"])
        for k in ("lx", "ly", "thick_xl", "thick_rg"):
            dk = res["drift"][k]
            log("  %-9s first %.3f last %.3f change %+.3f (%+.2f %%)" % (k, dk["first"], dk["last"], dk["change"], dk["change_pct"]))
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    log("[%s] wrote %s (%.0f s)" % (a.label, a.out, time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
