#!/usr/bin/env python3
"""cavity_widom.py — Cavity-biased Widom test-particle insertion for LAMMPS trajectories.

Standard Widom insertion fails in dense polymer gels because every random
insertion point overlaps with polymer, driving exp(-ΔU/kT) → 0.  Cavity-biased
Widom restricts insertions to void positions and corrects for the bias:

    μ_ex = -kT ln( p_cav × ⟨exp(-ΔU/kT)⟩_cavity )

where p_cav = N_cavity / N_trial is the fraction of trial points that qualify
as cavities (no real atom within r_cavity of the insertion point).

Usage
-----
    python cavity_widom.py  --traj widom_run.lammpstrj  \\
                            --eps-sp 1.0  --eps-ss 1.0   \\
                            --n-bins 20  --n-trial 2000  \\
                            --r-cavity 0.5  --temperature 1.0

Outputs (written to --out-dir, named with --out-stem):
    mu_z_cavity_<stem>.dat          — one row per frame per bin
    mu_z_cavity_summary_<stem>.dat  — time-averaged μ_ex(z) with stderr and p_cav

Atom-type → epsilon mapping (solvent ghost = type 3):
    type 1  polymer backbone    → eps_SP  (--eps-sp)
    type 2  crosslinker         → eps_SP  (--eps-sp)
    type 3  solvent             → eps_SS  (--eps-ss)
    type 4  support             → 0.0     (no solvent–support interaction)
    type 5  piston              → 1.0     (WCA with solvent)

All units are LJ reduced units (sigma=1, mass=1, epsilon=1 sets the scale).

Performance notes
-----------------
For large systems (400 k atoms) the two main costs per frame are:
  1. Trajectory I/O     — mitigated by bulk np.fromstring reads
  2. WCA energy         — mitigated by KD-tree neighbor lookup (only ~3–5
                          atoms within the 1.122σ WCA cutoff contribute;
                          O(N) full-distance loops are avoided)
eps_arr and the atom coordinate array are pre-computed once per frame and
reused across all bins and all cavity insertions.
"""

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np

try:
    from scipy.spatial import cKDTree
except ImportError:
    sys.exit("scipy is required: pip install scipy")


# ---------------------------------------------------------------------------
# Trajectory parser
# ---------------------------------------------------------------------------

def iter_frames(traj_path):
    """Yield (timestep, box, atoms_xyz, atom_types) tuples from a LAMMPS dump file.

    Returns
    -------
    timestep   : int
    box        : dict  xlo xhi ylo yhi zlo zhi
    atoms_xyz  : (N, 3) float64 array of (x, y, z) positions
    atom_types : (N,)   int32   array of atom types
    """
    with open(traj_path) as fh:
        while True:
            # ITEM: TIMESTEP
            line = fh.readline()
            if not line:
                return
            if "TIMESTEP" not in line:
                continue
            timestep = int(fh.readline().strip())

            # ITEM: NUMBER OF ATOMS
            fh.readline()
            n_atoms = int(fh.readline().strip())

            # ITEM: BOX BOUNDS
            fh.readline()
            xlo, xhi = map(float, fh.readline().split())
            ylo, yhi = map(float, fh.readline().split())
            zlo, zhi = map(float, fh.readline().split())
            box = dict(xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi)

            # ITEM: ATOMS id type x y z
            header = fh.readline()
            cols = header.split()[2:]   # e.g. ['id', 'type', 'x', 'y', 'z']

            # Bulk-read all atom lines — much faster than per-line readline for
            # large N because it avoids 400k+ Python function calls.
            raw = ''.join(fh.readline() for _ in range(n_atoms))
            data = np.fromstring(raw, sep=' ').reshape(n_atoms, len(cols))

            col_idx = {c: i for i, c in enumerate(cols)}
            atoms_xyz  = data[:, [col_idx['x'], col_idx['y'], col_idx['z']]]
            atom_types = data[:, col_idx['type']].astype(np.int32)

            yield timestep, box, atoms_xyz, atom_types


# ---------------------------------------------------------------------------
# WCA energy calculation — fully vectorized, KD-tree-aware
# ---------------------------------------------------------------------------

def wca_energy_vectorized(r2, eps, sigma=1.0, cutoff=1.122):
    """Compute WCA pair energy for arrays of squared distances and epsilons.

    Returns array of per-pair energies (0.0 for r >= cutoff).
    """
    rc2 = cutoff ** 2
    mask = r2 < rc2
    U = np.zeros_like(r2)
    if mask.any():
        inv_r2  = (sigma * sigma) / r2[mask]
        inv_r6  = inv_r2 * inv_r2 * inv_r2
        inv_r12 = inv_r6 * inv_r6
        U[mask] = 4.0 * eps[mask] * (inv_r12 - inv_r6) + eps[mask]
    return U


# ---------------------------------------------------------------------------
# Per-frame processing (optimized)
# ---------------------------------------------------------------------------

def process_frame(timestep, box, atoms_xyz, atom_types, eps_arr,
                  n_bins, n_trial, r_cavity,
                  temperature, sigma=1.0, cutoff=1.122,
                  rng=None,
                  cavity_types=None):
    """Run cavity-biased Widom for one snapshot.

    Parameters
    ----------
    atoms_xyz    : (N, 3) float64  — atom positions (world coordinates)
    atom_types   : (N,)   int32    — atom types
    eps_arr      : (N,)   float64  — epsilon for each atom vs. ghost solvent,
                                     PRE-COMPUTED once per frame in main()
    cavity_types : set of int or None
        Atom types that count as cavity-blockers (default: all types with
        eps > 0).  Use this to exclude piston (type 5) and support (type 4)
        from void detection so that boundary beads do not spuriously reduce
        p_cav or drag μ_ex toward the piston's lattice structure.

        Energy is ALWAYS computed from all atoms whose eps > 0, regardless
        of cavity_types — the ghost particle still feels piston WCA repulsion
        if inserted near a piston bead, which is physically correct.

    Returns list of dicts, one per z-bin:
        z_lo, z_hi, z_center, n_trial, n_cavity, p_cav,
        mu_ex (NaN if no cavity insertions), beta_dU_mean,
        skipped (bool) — True if the bin is outside the active fluid region
    """
    if rng is None:
        rng = np.random.default_rng()

    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    boxsize = np.array([lx, ly, lz])

    # Wrap atom positions to [0, L) for PBC-aware KD-tree
    xyz_wrapped = (atoms_xyz - np.array([box['xlo'], box['ylo'], box['zlo']])) % boxsize

    # ── Separate cavity-detection atoms from energy atoms ──────────────────
    # Cavity detection uses only fluid atoms (types 1,2,3 by default) so that
    # piston (type 5) and support (type 4) beads do not falsely reduce the
    # void fraction.  Energy calculation uses all atoms with eps > 0, so the
    # ghost particle still feels WCA repulsion from piston if it's close.
    if cavity_types is None:
        # Default: any type with non-zero epsilon counts as a cavity blocker
        cav_mask = eps_arr > 0.0
    else:
        cav_mask = np.isin(atom_types, list(cavity_types))

    xyz_cav     = xyz_wrapped[cav_mask]   # (N_cav, 3)
    eps_nz_mask = eps_arr > 0.0           # for energy: all interacting atoms

    # Build KD-trees — one for cavity detection, one for energy lookup
    # (same tree if cavity_types covers all eps>0 atoms)
    if cav_mask.sum() == 0:
        tree_cav = None
    else:
        tree_cav = cKDTree(xyz_cav, boxsize=boxsize)

    if eps_nz_mask.sum() == 0:
        tree_en = None
    else:
        tree_en  = cKDTree(xyz_wrapped[eps_nz_mask], boxsize=boxsize)
        eps_en   = eps_arr[eps_nz_mask]
        xyz_en   = xyz_wrapped[eps_nz_mask]

    # ── Active z-range: skip bins that are inside piston or support ─────────
    # Piston (type 5): piston beads sit at the top of the box.
    #   → exclude bins where z_lo >= min z of piston atoms
    # Support (type 4): support beads sit at the bottom of the box.
    #   → exclude bins where z_hi <= max z of support atoms
    PISTON_TYPE  = 5
    SUPPORT_TYPE = 4
    piston_mask  = atom_types == PISTON_TYPE
    support_mask = atom_types == SUPPORT_TYPE
    z_piston_min = float(atoms_xyz[piston_mask,  2].min()) if piston_mask.any()  else box['zhi']
    z_support_max= float(atoms_xyz[support_mask, 2].max()) if support_mask.any() else box['zlo']

    kT   = temperature
    beta = 1.0 / kT

    bin_edges = np.linspace(box['zlo'], box['zhi'], n_bins + 1)
    results   = []

    for ib in range(n_bins):
        z_lo     = bin_edges[ib]
        z_hi     = bin_edges[ib + 1]
        z_center = 0.5 * (z_lo + z_hi)

        # Skip bins fully inside the piston or fully inside the support
        if z_lo >= z_piston_min or z_hi <= z_support_max:
            results.append(dict(
                z_lo=z_lo, z_hi=z_hi, z_center=z_center,
                n_trial=0, n_cavity=0, p_cav=np.nan,
                mu_ex=np.nan, beta_dU_mean=np.nan, skipped=True,
            ))
            continue

        # --- Generate trial positions uniformly in this z-slab ---------------
        xs = rng.uniform(box['xlo'], box['xhi'], n_trial)
        ys = rng.uniform(box['ylo'], box['yhi'], n_trial)
        zs = rng.uniform(z_lo, z_hi,            n_trial)

        # Wrap to [0, L) for the KD-tree query
        trial_pts = np.column_stack([
            (xs - box['xlo']) % lx,
            (ys - box['ylo']) % ly,
            (zs - box['zlo']) % lz,
        ])

        # --- Cavity detection: use fluid-only tree (excludes piston/support) -
        if tree_cav is None:
            cavity_mask = np.ones(n_trial, dtype=bool)
        else:
            nearby      = tree_cav.query_ball_point(trial_pts, r=r_cavity, workers=-1)
            cavity_mask = np.array([len(nb) == 0 for nb in nearby], dtype=bool)

        n_cav = int(cavity_mask.sum())
        p_cav = n_cav / n_trial

        if n_cav == 0:
            results.append(dict(
                z_lo=z_lo, z_hi=z_hi, z_center=z_center,
                n_trial=n_trial, n_cavity=0, p_cav=0.0,
                mu_ex=np.nan, beta_dU_mean=np.nan, skipped=False,
            ))
            continue

        # --- WCA energy at cavity positions only -----------------------------
        # Use KD-tree to find only the ~3–5 atoms within WCA cutoff.
        # Energy tree includes all atoms with eps > 0 (including piston type 5).
        cav_pts = trial_pts[cavity_mask]         # (n_cav, 3) in [0,L) space
        beta_dU = np.empty(n_cav)

        if tree_en is None:
            beta_dU[:] = 0.0
        else:
            nbr_lists = tree_en.query_ball_point(cav_pts, r=cutoff, workers=-1)
            for k, nbrs in enumerate(nbr_lists):
                if len(nbrs) == 0:
                    beta_dU[k] = 0.0
                    continue
                nbrs_arr = np.asarray(nbrs, dtype=np.intp)
                pos      = cav_pts[k]                    # (3,) in [0, L)
                d        = xyz_en[nbrs_arr] - pos        # (m, 3)
                d       -= np.round(d / boxsize) * boxsize
                r2       = (d * d).sum(axis=1)
                U        = wca_energy_vectorized(r2, eps_en[nbrs_arr], sigma, cutoff)
                beta_dU[k] = beta * float(U.sum())

        # Clamp to avoid overflow in exp (underflow → 0 is correct and fine)
        beta_dU_clipped = np.minimum(beta_dU, 700.0)
        exp_neg_bdU     = np.exp(-beta_dU_clipped)
        avg_exp         = float(exp_neg_bdU.mean())

        if avg_exp > 0.0 and p_cav > 0.0:
            mu_ex = -kT * np.log(p_cav * avg_exp)
        else:
            mu_ex = np.nan

        results.append(dict(
            z_lo=z_lo, z_hi=z_hi, z_center=z_center,
            n_trial=n_trial, n_cavity=n_cav, p_cav=p_cav,
            mu_ex=mu_ex, beta_dU_mean=float(beta_dU.mean()), skipped=False,
        ))

    return results


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

FRAME_HEADER = (
    "# Cavity-biased Widom excess chemical potential μ_ex(z)\n"
    "# step  z_center  z_lo  z_hi  p_cav  n_cavity  n_trial  mu_ex\n"
)

SUMMARY_HEADER = (
    "# Cavity-biased Widom — time-averaged summary\n"
    "# z_center  z_lo  z_hi  mu_ex_mean  mu_ex_stderr  p_cav_mean  n_frames\n"
)


def write_frame(fh, timestep, bin_results):
    for b in bin_results:
        if b.get('skipped', False):
            continue   # piston/support bins — don't write to frame file
        mu_str  = f"{b['mu_ex']:.6f}"  if not np.isnan(b['mu_ex'])  else "nan"
        pcav_str= f"{b['p_cav']:.5f}"  if not np.isnan(b['p_cav'])  else "nan"
        fh.write(
            f"{timestep}  {b['z_center']:.4f}  {b['z_lo']:.4f}  {b['z_hi']:.4f}"
            f"  {pcav_str}  {b['n_cavity']}  {b['n_trial']}  {mu_str}\n"
        )


def write_summary(path, all_results, n_bins):
    """Collect per-bin time series and write mean ± stderr."""
    mu_series = [[] for _ in range(n_bins)]
    pc_series = [[] for _ in range(n_bins)]
    z_info    = [None] * n_bins

    for frame_bins in all_results:
        for ib, b in enumerate(frame_bins):
            z_info[ib] = (b['z_center'], b['z_lo'], b['z_hi'])
            if not np.isnan(b['mu_ex']):
                mu_series[ib].append(b['mu_ex'])
            pc_series[ib].append(b['p_cav'])

    with open(path, 'w') as fh:
        fh.write(SUMMARY_HEADER)
        for ib in range(n_bins):
            zc, zlo, zhi = z_info[ib]
            mus  = np.array(mu_series[ib])
            pcs  = np.array(pc_series[ib])
            n_ok = len(mus)

            if n_ok >= 2:
                mu_mean   = float(mus.mean())
                mu_stderr = float(mus.std(ddof=1) / np.sqrt(n_ok))
            elif n_ok == 1:
                mu_mean, mu_stderr = float(mus[0]), np.nan
            else:
                mu_mean = mu_stderr = np.nan

            pc_mean = float(pcs.mean()) if len(pcs) > 0 else np.nan

            mu_str = f"{mu_mean:.6f}"   if not np.isnan(mu_mean)   else "nan"
            se_str = f"{mu_stderr:.6f}" if not np.isnan(mu_stderr) else "nan"
            pc_str = f"{pc_mean:.5f}"   if not np.isnan(pc_mean)   else "nan"

            fh.write(
                f"{zc:.4f}  {zlo:.4f}  {zhi:.4f}"
                f"  {mu_str}  {se_str}  {pc_str}  {n_ok}\n"
            )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_eps_map(eps_sp, eps_ss):
    """Return dict mapping atom type → epsilon for WCA with ghost solvent."""
    return {
        1: eps_sp,  # polymer backbone
        2: eps_sp,  # crosslinker
        3: eps_ss,  # solvent
        4: 0.0,     # support (no interaction with solvent)
        5: 1.0,     # piston (WCA with solvent)
    }


def main():
    parser = argparse.ArgumentParser(
        description="Cavity-biased Widom insertion for LAMMPS trajectories.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--traj",        required=True,
                        help="LAMMPS dump file (id type x y z)")
    parser.add_argument("--n-bins",      type=int,   default=20,
                        help="Number of z-bins")
    parser.add_argument("--n-trial",     type=int,   default=2000,
                        help="Trial insertions per bin per frame")
    parser.add_argument("--r-cavity",    type=float, default=0.5,
                        help="Cavity radius in sigma units")
    parser.add_argument("--temperature", type=float, default=1.0,
                        help="kT in LJ units")
    parser.add_argument("--eps-sp",      type=float, default=1.0,
                        help="Solvent–polymer epsilon")
    parser.add_argument("--eps-ss",      type=float, default=1.0,
                        help="Solvent–solvent epsilon")
    parser.add_argument("--sigma",       type=float, default=1.0,
                        help="LJ sigma")
    parser.add_argument("--cutoff",      type=float, default=1.122,
                        help="WCA cutoff (default = 2^(1/6) sigma)")
    parser.add_argument("--seed",        type=int,   default=42,
                        help="Random seed for reproducibility")
    parser.add_argument("--out-dir",     default=None,
                        help="Output directory (default: same as --traj)")
    parser.add_argument("--out-stem",    default=None,
                        help="Output filename stem, e.g. 'dataname_interaction_totsteps'. "
                             "Produces mu_z_cavity_{stem}.dat and "
                             "mu_z_cavity_summary_{stem}.dat.  "
                             "Defaults to the trajectory filename stem.")
    parser.add_argument("--cavity-types", default="1,2,3",
                        help="Comma-separated atom types that count as cavity blockers "
                             "(default: '1,2,3' = polymer + solvent only). "
                             "Piston (type 5) and support (type 4) are excluded by "
                             "default so their lattice structure does not distort "
                             "void-fraction measurements.")
    args = parser.parse_args()

    traj_path = Path(args.traj)
    if not traj_path.exists():
        sys.exit(f"Trajectory not found: {traj_path}")

    out_dir = Path(args.out_dir) if args.out_dir else traj_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    stem         = args.out_stem if args.out_stem else traj_path.stem
    frame_path   = out_dir / f"mu_z_cavity_{stem}.dat"
    summary_path = out_dir / f"mu_z_cavity_summary_{stem}.dat"

    eps_map      = build_eps_map(args.eps_sp, args.eps_ss)
    rng          = np.random.default_rng(args.seed)
    cavity_types = set(int(t) for t in args.cavity_types.split(','))

    print(f"Trajectory   : {traj_path}")
    print(f"z-bins       : {args.n_bins}")
    print(f"trials/bin   : {args.n_trial}")
    print(f"r_cavity     : {args.r_cavity} σ")
    print(f"kT           : {args.temperature}")
    print(f"eps_SP       : {args.eps_sp}   eps_SS : {args.eps_ss}")
    print(f"WCA cutoff   : {args.cutoff} σ")
    print(f"cavity_types : {sorted(cavity_types)}  (piston/support excluded from void detection)")
    print(f"Frame file   : {frame_path}")
    print(f"Summary      : {summary_path}")
    print()

    all_results = []

    with open(frame_path, 'w') as fh_frame:
        fh_frame.write(FRAME_HEADER)

        for i_frame, (timestep, box, atoms_xyz, atom_types) in \
                enumerate(iter_frames(str(traj_path))):

            print(f"  Frame {i_frame+1}  step {timestep}  "
                  f"n_atoms {len(atoms_xyz)} ...", end=" ", flush=True)

            # ── Pre-compute eps_arr ONCE per frame ──────────────────────────
            # This avoids a Python loop over N_atoms inside every insertion.
            eps_arr = np.array([eps_map.get(int(t), 0.0) for t in atom_types],
                               dtype=np.float64)

            bin_results = process_frame(
                timestep, box, atoms_xyz, atom_types, eps_arr,
                n_bins       = args.n_bins,
                n_trial      = args.n_trial,
                r_cavity     = args.r_cavity,
                temperature  = args.temperature,
                sigma        = args.sigma,
                cutoff       = args.cutoff,
                rng          = rng,
                cavity_types = cavity_types,
            )

            write_frame(fh_frame, timestep, bin_results)
            all_results.append(bin_results)

            n_nan    = sum(1 for b in bin_results if np.isnan(b['mu_ex']))
            mean_pc  = np.mean([b['p_cav'] for b in bin_results])
            print(f"done  (p_cav_mean={mean_pc:.3f}, {n_nan}/{args.n_bins} bins NaN)")

    write_summary(str(summary_path), all_results, args.n_bins)

    print(f"\nWrote {len(all_results)} frames → {frame_path.name}")
    print(f"Summary → {summary_path.name}")

    # Quick sanity check: print μ_ex range from summary
    try:
        data = np.genfromtxt(str(summary_path), comments='#')
        if data.ndim == 2 and data.shape[1] >= 4:
            mu    = data[:, 3]
            valid = mu[~np.isnan(mu)]
            if len(valid):
                print(f"\nμ_ex range: {valid.min():.4f} to {valid.max():.4f}  "
                      f"(Δμ_ex = {valid.max()-valid.min():.4f})")
                if (valid.max() - valid.min()) > 0.05:
                    print("  *** WARNING: Δμ_ex > 0.05 ε — significant chemical "
                          "potential gradient ***")
    except Exception:
        pass


if __name__ == "__main__":
    main()
