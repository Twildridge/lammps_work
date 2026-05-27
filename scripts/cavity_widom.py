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
    mu_z_cavity_<stem>.dat            — one row per frame per bin
    mu_z_cavity_summary_<stem>.dat    — time-averaged μ_ex(z) with stderr and p_cav
    mu_total_diagnostic_<stem>.png    — 4-panel chemical-equilibrium diagnostic
                                         (rho_s, mu_ex, mu_total, p_p/P_ext)

Diagnostic plot
---------------
For an inhomogeneous fluid mu_ex(z) need NOT be flat at equilibrium — what is
flat is mu_total(z) = mu_ex(z) + kT*ln(rho_s(z)).  The script computes this
automatically (rho_s is histogrammed from the trajectory itself in the same
z-bins) and saves a 4-panel diagnostic:

    panel 1   rho_s(z)                                 — solvent density
    panel 2   mu_ex(z)                                 — what cavity-Widom returns
    panel 3   mu_total = mu_ex + kT*ln(rho_s)           — should be FLAT
    panel 4   p_p(z)/P_ext = 1 + (mu_total - mu_total,res)*rho_s,res/P_ext
              — pore pressure, referenced so reservoir = P_ext by construction

P_ext is set via --p-ext (LAMMPS barostat target; default 1.5 eps/sigma^3).
Pass --no-plot to skip the plot if matplotlib isn't available.

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
# Solvent density binning — same z-grid as cavity-Widom
# ---------------------------------------------------------------------------

def bin_solvent_density(atoms_xyz, atom_types, box, n_bins, solvent_type=3):
    """Histogram solvent (type-3) atoms in z-bins and divide by bin volume.

    Uses the SAME z-bin grid as process_frame so the resulting rho_s(z)
    can be added directly to mu_ex(z) to form mu_total(z) without
    interpolation.

    Returns
    -------
    rho_z : (n_bins,) float64  — number density per bin in sigma^-3
    """
    sol_mask  = atom_types == solvent_type
    z_solv    = atoms_xyz[sol_mask, 2]
    bin_edges = np.linspace(box['zlo'], box['zhi'], n_bins + 1)
    counts, _ = np.histogram(z_solv, bins=bin_edges)
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    dz = (box['zhi'] - box['zlo']) / n_bins
    return counts.astype(np.float64) / (lx * ly * dz)


# ---------------------------------------------------------------------------
# Per-frame processing (optimized)
# ---------------------------------------------------------------------------

def process_frame(timestep, box, atoms_xyz, atom_types, eps_arr,
                  n_bins, n_trial, r_cavity,
                  temperature, sigma=1.0, cutoff=1.122,
                  rng=None,
                  cavity_types=None,
                  exclusion_buffer=None,
                  xy_extent=None,
                  xy_buffer=None,
                  wall_type=6):
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
    exclusion_buffer : float or None
        Buffer (in sigma) applied around the piston extent and support top.
        Bins whose z_center falls within this distance of the support top OR
        within the window [z_piston_min - buffer, z_piston_max + buffer] are
        marked skipped=True and returned as NaN.  Bins above z_piston_max +
        buffer (i.e. the reservoir above the piston in slab_with_flow) are
        kept active.  Defaults to r_cavity if None; set to at least one bin
        width (e.g. 2.0 sigma) for slab_with_flow to avoid interface spikes.
    xy_extent : (x_lo, x_hi, y_lo, y_hi) or None
        Lateral sampling box for trial-point insertions.  When None,
        auto-detected: if any atoms of type ``wall_type`` are present in this
        frame, xs/ys are restricted to the wall interior with an inward buffer
        of ``xy_buffer`` (so trial points cannot land in the empty vacuum
        OUTSIDE a lateral wall tube — that void otherwise inflates p_cav with
        zero-energy "cavities" and creates a spurious μ_ex floor at
        −kT ln(void_fraction)).  If no wall atoms are present (e.g.
        slab_with_support), the full box [xlo, xhi] × [ylo, yhi] is used and
        behavior is identical to pre-patch.
    xy_buffer : float or None
        Inward buffer (σ) subtracted from the wall extent when auto-detecting
        xy_extent.  Defaults to ``cutoff`` (1.122 σ for WCA), which keeps trial
        points at least one WCA cutoff away from the wall plane so a ghost
        solvent never overlaps a wall bead.  Ignored if xy_extent is given
        explicitly.
    wall_type : int
        Atom type used for the lateral wall (default 6, matching
        slab_with_flow.lmp).  Set to a value that does not appear in the
        trajectory to disable auto-detection.

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

    # ── Active z-range: skip bins inside/near the piston or support ──────────
    # Support (type 4): sits at the bottom; exclude bins ≤ support_top + excl_buf.
    # Piston  (type 5): can be mid-box (slab_with_flow) or at the ceiling
    #   (slab_with_support).  We exclude a window AROUND the piston extent
    #   [z_piston_min - excl_buf, z_piston_max + excl_buf] rather than
    #   everything above z_piston_min.  This is critical for slab_with_flow
    #   where a solvent reservoir sits ABOVE the piston: bins clearly above the
    #   piston top (z_center > z_piston_max + excl_buf) must remain active.
    PISTON_TYPE  = 5
    SUPPORT_TYPE = 4
    piston_mask  = atom_types == PISTON_TYPE
    support_mask = atom_types == SUPPORT_TYPE
    z_piston_min = float(atoms_xyz[piston_mask,  2].min()) if piston_mask.any()  else box['zhi']
    z_piston_max = float(atoms_xyz[piston_mask,  2].max()) if piston_mask.any()  else box['zhi']
    z_support_max= float(atoms_xyz[support_mask, 2].max()) if support_mask.any() else box['zlo']

    kT   = temperature
    beta = 1.0 / kT

    # exclusion_buffer is decoupled from r_cavity so it can be widened without
    # affecting the cavity-detection radius.  Default: same as r_cavity.
    excl_buf = exclusion_buffer if exclusion_buffer is not None else r_cavity

    # ── Lateral (x, y) sampling box for trial-point insertions ─────────────
    # Auto-detect wall confinement from atoms of `wall_type` (default 6).  This
    # is essential for slab_with_flow.lmp, where a hollow tube of wall beads
    # confines polymer + solvent to the interior; sampling trial points in the
    # full periodic box would place ~30 % of insertions in the empty vacuum
    # outside the tube, where cavity detection trivially passes and the WCA
    # energy is zero — producing a spurious μ_ex floor at −kT ln(void_fraction)
    # that wipes out the gel/reservoir contrast.
    #
    # slab_with_support.lmp has no wall atoms; auto-detection then falls back
    # to the full box, reproducing pre-patch behavior bit-for-bit.
    xy_buf_eff = xy_buffer if xy_buffer is not None else cutoff
    if xy_extent is not None:
        x_lo_s, x_hi_s, y_lo_s, y_hi_s = (float(v) for v in xy_extent)
    else:
        wall_mask = atom_types == wall_type
        if wall_mask.any():
            xw = atoms_xyz[wall_mask, 0]
            yw = atoms_xyz[wall_mask, 1]
            x_lo_s = float(xw.min()) + xy_buf_eff
            x_hi_s = float(xw.max()) - xy_buf_eff
            y_lo_s = float(yw.min()) + xy_buf_eff
            y_hi_s = float(yw.max()) - xy_buf_eff
            # Degenerate-buffer guard: if the user passed a buffer wider than
            # the tube, fall back to the raw wall extent (no inward shrink).
            if x_hi_s <= x_lo_s:
                x_lo_s, x_hi_s = float(xw.min()), float(xw.max())
            if y_hi_s <= y_lo_s:
                y_lo_s, y_hi_s = float(yw.min()), float(yw.max())
        else:
            x_lo_s, x_hi_s = box['xlo'], box['xhi']
            y_lo_s, y_hi_s = box['ylo'], box['yhi']

    bin_edges = np.linspace(box['zlo'], box['zhi'], n_bins + 1)
    results   = []

    for ib in range(n_bins):
        z_lo     = bin_edges[ib]
        z_hi     = bin_edges[ib + 1]
        z_center = 0.5 * (z_lo + z_hi)

        # Skip bins that are:
        #   (a) within excl_buf of the support top, OR
        #   (b) in the gel-side interface buffer [z_piston_min-excl_buf, z_piston_min), OR
        #   (c) in the reservoir-side interface buffer (z_piston_max, z_piston_max+excl_buf].
        #
        # Bins INSIDE the piston body [z_piston_min, z_piston_max] are kept active so
        # that μ_ex can be measured there.  Set --piston-eps 0.0 for compression
        # mode (solvent transparent to piston) so the ghost particle does not
        # artificially see WCA repulsion from piston beads in the energy calc.
        near_support      = z_center <= z_support_max + excl_buf
        gel_side_iface    = (z_piston_min - excl_buf <= z_center < z_piston_min)
        res_side_iface    = (z_piston_max < z_center <= z_piston_max + excl_buf)
        if near_support or gel_side_iface or res_side_iface:
            results.append(dict(
                z_lo=z_lo, z_hi=z_hi, z_center=z_center,
                n_trial=0, n_cavity=0, p_cav=np.nan,
                mu_ex=np.nan, beta_dU_mean=np.nan, skipped=True,
            ))
            continue

        # --- Generate trial positions uniformly in this z-slab ---------------
        # xs/ys are sampled inside the lateral wall tube (x_lo_s..x_hi_s, etc.)
        # when walls are present; otherwise the full box.  See the xy_extent
        # auto-detect block above for why this matters in slab_with_flow.
        xs = rng.uniform(x_lo_s, x_hi_s, n_trial)
        ys = rng.uniform(y_lo_s, y_hi_s, n_trial)
        zs = rng.uniform(z_lo,   z_hi,   n_trial)

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
# Diagnostic plot
# ---------------------------------------------------------------------------

def make_diagnostic_plot(all_results, rho_per_frame, steps, kT,
                         p_ext, out_path, title="",
                         res_thresh=0.85, gel_thresh=0.70):
    """4-panel chemical-equilibrium diagnostic.

    Panels (top→bottom):
        1. rho_s(z)              solvent number density
        2. mu_ex(z)              cavity-Widom excess chemical potential
        3. mu_total(z)           = mu_ex + kT*ln(rho_s)   — should be FLAT
        4. p_p(z)/P_ext          = mu_total * rho_s,res / P_ext

    Parameters
    ----------
    all_results : list[list[dict]]
        Per-frame bin results from process_frame.  Outer list length = n_frames,
        inner list length = n_bins.
    rho_per_frame : (n_frames, n_bins) float64
        Solvent number density per bin per frame — must use the SAME bin grid
        as all_results (use bin_solvent_density with the same n_bins).
    steps : list[int]
        LAMMPS timestep for each frame (for legend ordering).
    kT : float
        Reduced temperature.
    p_ext : float
        LAMMPS barostat target pressure (eps/sigma^3).  Used to normalize p_p.
    out_path : str | Path
        Where to save the PNG.
    title : str
        Plot title (typically the run stem).
    res_thresh, gel_thresh : float
        Reservoir is bins with rho >= res_thresh*rho_max.
        Gel interior  is bins with rho <  gel_thresh*rho_max.
    """
    # Force a non-interactive backend BEFORE importing pyplot.  On HPC nodes
    # (Expanse, Bridges, POD) there is no display, and the default backend
    # raises before plt.savefig() ever runs.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[diagnostic plot] matplotlib unavailable / backend error: {e}")
        print( "                  skipping plot (run again with --no-plot to suppress this)")
        return

    n_frames = len(all_results)
    if n_frames == 0:
        print("[diagnostic plot] no frames in all_results — skipping")
        return
    print(f"[diagnostic plot] building 4-panel figure: {n_frames} frames, "
          f"{len(all_results[0])} bins → {out_path}")

    # ── Extract per-frame arrays ──────────────────────────────────────────────
    z_centers = np.array([b['z_center'] for b in all_results[0]])
    z_lo      = np.array([b['z_lo']     for b in all_results[0]])
    z_hi      = np.array([b['z_hi']     for b in all_results[0]])

    mu_per_frame = np.array(
        [[b['mu_ex'] for b in frame] for frame in all_results],
        dtype=np.float64,
    )  # (n_frames, n_bins)  — NaN where bin was skipped or no cavity

    # ── Time-averaged μ_ex (NaN-safe) ─────────────────────────────────────────
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        n_ok    = np.sum(~np.isnan(mu_per_frame), axis=0)
        mu_mean = np.nanmean(mu_per_frame, axis=0)
        mu_std  = np.nanstd(mu_per_frame, axis=0, ddof=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        mu_se = np.where(n_ok >= 2, mu_std / np.sqrt(np.maximum(n_ok, 1)), np.nan)

    rho_mean = rho_per_frame.mean(axis=0)
    rho_se   = (rho_per_frame.std(axis=0, ddof=1) / np.sqrt(n_frames)
                if n_frames >= 2 else np.zeros_like(rho_mean))

    # ── Reservoir + gel regions ───────────────────────────────────────────────
    rho_max  = rho_mean.max() if rho_mean.size else np.nan
    res_mask = (rho_mean >= res_thresh * rho_max) & np.isfinite(mu_mean)
    gel_mask = (rho_mean <  gel_thresh * rho_max) & np.isfinite(mu_mean)

    rho_res = rho_mean[res_mask].mean() if res_mask.any() else np.nan
    rho_gel = rho_mean[gel_mask].mean() if gel_mask.any() else np.nan
    V_bar   = 1.0 / rho_res if rho_res and np.isfinite(rho_res) else np.nan

    # ── μ_total = μ_ex + kT*ln(ρ_s) ───────────────────────────────────────────
    with np.errstate(divide='ignore', invalid='ignore'):
        mu_total_mean      = mu_mean      + kT * np.log(rho_mean)
        mu_total_se        = mu_se   # ρ uncertainty is negligible relative to μ_ex
        mu_total_per_frame = mu_per_frame + kT * np.log(rho_per_frame)

    # ── Summary numbers for annotations ───────────────────────────────────────
    mu_ex_res    = np.nanmean(mu_mean[res_mask])       if res_mask.any() else np.nan
    mu_ex_gel    = np.nanmean(mu_mean[gel_mask])       if gel_mask.any() else np.nan
    mu_total_res = np.nanmean(mu_total_mean[res_mask]) if res_mask.any() else np.nan
    mu_total_gel = np.nanmean(mu_total_mean[gel_mask]) if gel_mask.any() else np.nan

    # ── Pore pressure (normalized, referenced to reservoir) ───────────────────
    # Cavity-Widom mu_ex has an implicit Lambda=1 convention, so the absolute
    # value of mu_total is convention-dependent.  Reference it to the reservoir
    # so that p_p,res = P_ext by construction:
    #   p_p(z) = P_ext + (mu_total(z) - mu_total_res) / V_bar
    #   p_p(z) / P_ext = 1 + (mu_total(z) - mu_total_res) * rho_s,res / P_ext
    pp_mean      = 1.0 + (mu_total_mean      - mu_total_res) * rho_res / p_ext
    pp_per_frame = 1.0 + (mu_total_per_frame - mu_total_res) * rho_res / p_ext
    pp_res       = np.nanmean(pp_mean[res_mask])       if res_mask.any() else np.nan
    pp_gel       = np.nanmean(pp_mean[gel_mask])       if gel_mask.any() else np.nan
    delta_mu_ex_meas = (mu_ex_gel - mu_ex_res) if np.isfinite(mu_ex_gel) and np.isfinite(mu_ex_res) else np.nan
    delta_mu_ex_pred = (kT * np.log(rho_res / rho_gel)
                        if np.isfinite(rho_res) and np.isfinite(rho_gel) and rho_gel > 0
                        else np.nan)
    delta_mu_total   = (mu_total_gel - mu_total_res) if np.isfinite(mu_total_gel) and np.isfinite(mu_total_res) else np.nan
    delta_pp         = (pp_gel - pp_res) if np.isfinite(pp_gel) and np.isfinite(pp_res) else np.nan

    # ── Helper: linear interpolation across NaN gaps ──────────────────────────
    def _interp_gap(z, y):
        """Return list of (z_seg, y_seg) linearly interpolated across each NaN
        gap in y.  Each segment includes the flanking non-NaN endpoints so the
        dashed bridge meets the solid line seamlessly."""
        finite = np.isfinite(y)
        segs = []
        n = len(z)
        i = 0
        while i < n:
            if not finite[i]:
                j = i
                while j < n and not finite[j]:
                    j += 1
                if i > 0 and j < n:
                    z_seg = np.concatenate([[z[i - 1]], z[i:j], [z[j]]])
                    y_seg = np.interp(z_seg, [z[i - 1], z[j]],
                                              [y[i - 1], y[j]])
                    segs.append((z_seg, y_seg))
                i = max(j, i + 1)
            else:
                i += 1
        return segs

    # ── Identify piston zone: NaN gap in mu_mean flanked by data on both sides ─
    # (excludes support region which is NaN only at one edge)
    piston_spans = []
    finite_mean = np.isfinite(mu_mean)
    _i = 0
    while _i < len(mu_mean):
        if not finite_mean[_i]:
            _j = _i
            while _j < len(mu_mean) and not finite_mean[_j]:
                _j += 1
            if _i > 0 and _j < len(mu_mean):   # internal NaN gap = piston zone
                piston_spans.append((z_lo[_i], z_hi[_j - 1]))
            _i = max(_j, _i + 1)
        else:
            _i += 1

    # ── Build the figure ──────────────────────────────────────────────────────
    fig, axes = plt.subplots(
        4, 1, figsize=(10, 14), sharex=True,
        gridspec_kw={'hspace': 0.10, 'height_ratios': [1, 1.2, 1.2, 1.2]},
    )
    cmap_t = plt.get_cmap('viridis')

    def _frame_color(k):
        return cmap_t(k / max(n_frames - 1, 1))

    # Panel 1: ρ_s(z)
    ax = axes[0]
    ax.plot(z_centers, rho_mean, 'o-', color='steelblue',
            markersize=3, lw=1.4, label=r"$\langle\rho_s(z)\rangle$")
    if n_frames >= 2:
        ax.fill_between(z_centers, rho_mean - rho_se, rho_mean + rho_se,
                        color='steelblue', alpha=0.25, lw=0)
    if np.isfinite(rho_res):
        ax.axhline(rho_res, color='steelblue', ls=':', lw=1, alpha=0.7,
                   label=fr"$\rho_{{s,\rm res}}={rho_res:.3g}\,\sigma^{{-3}}$")
    for k in np.where(res_mask)[0]:
        ax.axvspan(z_lo[k], z_hi[k], color='lightblue', alpha=0.15, lw=0)
    ax.set_ylabel(r"$\rho_s$  ($\sigma^{-3}$)")
    ax.set_title(
        (title + "  —  " if title else "") +
        r"cavity-Widom diagnostic  ($\mu_{\rm total}$ flat $\Rightarrow$ chemical equilibrium)",
        fontsize=11,
    )
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.25)

    # Panel 2: μ_ex(z)
    ax = axes[1]
    for zlo_p, zhi_p in piston_spans:
        ax.axvspan(zlo_p, zhi_p, color='lightgray', alpha=0.35, lw=0,
                   label='piston zone (interpolated)')
    for k in range(n_frames):
        is_last = (k == n_frames - 1)
        color = 'darkorange' if is_last else _frame_color(k)
        ax.plot(z_centers, mu_per_frame[k], '-',
                color=color,
                lw=2.0 if is_last else 0.7,
                alpha=1.0 if is_last else 0.45,
                zorder=4 if is_last else 1,
                label=f'step {steps[k]:,}' if n_frames <= 12 else None)
        for z_seg, y_seg in _interp_gap(z_centers, mu_per_frame[k]):
            ax.plot(z_seg, y_seg, '--', color=color,
                    lw=1.6 if is_last else 0.5,
                    alpha=0.75 if is_last else 0.25,
                    zorder=4 if is_last else 0)
    if np.isfinite(mu_ex_res):
        ax.axhline(mu_ex_res, color='steelblue', ls='--', lw=1, alpha=0.8,
                   label=fr"$\langle\mu_{{ex}}\rangle_{{\rm res}}={mu_ex_res:.3g}\,\varepsilon$")
    if np.isfinite(mu_ex_gel):
        ax.axhline(mu_ex_gel, color='firebrick', ls='--', lw=1, alpha=0.8,
                   label=fr"$\langle\mu_{{ex}}\rangle_{{\rm gel}}={mu_ex_gel:.3g}\,\varepsilon$")
    ax.text(0.02, 0.04,
            fr"$\Delta\mu_{{ex}}$(gel$-$res) measured = ${delta_mu_ex_meas:+.3g}\,\varepsilon$"
            "\n"
            fr"  equilibrium $kT\ln(\rho_{{\rm res}}/\rho_{{\rm gel}})$ = ${delta_mu_ex_pred:+.3g}\,\varepsilon$",
            transform=ax.transAxes, fontsize=9, va='bottom',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.85))
    ax.set_ylabel(r"$\mu_{ex}$  ($\varepsilon$)")
    ax.legend(loc='lower right', fontsize=8, framealpha=0.9, ncol=2)
    ax.grid(alpha=0.25)

    # Panel 3: μ_total(z)
    ax = axes[2]
    for zlo_p, zhi_p in piston_spans:
        ax.axvspan(zlo_p, zhi_p, color='lightgray', alpha=0.35, lw=0)
    for k in range(n_frames):
        is_last = (k == n_frames - 1)
        color = 'darkorange' if is_last else _frame_color(k)
        ax.plot(z_centers, mu_total_per_frame[k], '-',
                color=color,
                lw=2.0 if is_last else 0.7,
                alpha=1.0 if is_last else 0.45,
                zorder=4 if is_last else 1)
        for z_seg, y_seg in _interp_gap(z_centers, mu_total_per_frame[k]):
            ax.plot(z_seg, y_seg, '--', color=color,
                    lw=1.6 if is_last else 0.5,
                    alpha=0.75 if is_last else 0.25,
                    zorder=4 if is_last else 0)
    if np.isfinite(mu_total_res):
        ax.axhline(mu_total_res, color='steelblue', ls='--', lw=1.2, alpha=0.9,
                   label=fr"$\mu_{{\rm total,res}}={mu_total_res:.3g}\,\varepsilon$  (equilibrium ref)")
    if np.isfinite(mu_total_gel):
        ax.axhline(mu_total_gel, color='firebrick', ls='--', lw=1, alpha=0.7,
                   label=fr"$\langle\mu_{{\rm total}}\rangle_{{\rm gel}}={mu_total_gel:.3g}\,\varepsilon$")

    ax.text(0.02, 0.04,
            fr"$\Delta\mu_{{\rm total}}$(gel$-$res) $= {delta_mu_total:+.3g}\,\varepsilon$",
            transform=ax.transAxes, fontsize=9, va='bottom',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.85))
    ax.set_ylabel(r"$\mu_{\rm total} = \mu_{ex} + kT\ln\rho_s$" "\n" r"($\varepsilon$)")
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.25)

    # Panel 4: p_p / P_ext
    ax = axes[3]
    for zlo_p, zhi_p in piston_spans:
        ax.axvspan(zlo_p, zhi_p, color='lightgray', alpha=0.35, lw=0,
                   label='piston zone (interpolated)')
    for k in range(n_frames):
        is_last = (k == n_frames - 1)
        color = 'darkorange' if is_last else _frame_color(k)
        ax.plot(z_centers, pp_per_frame[k], '-',
                color=color,
                lw=2.0 if is_last else 0.7,
                alpha=1.0 if is_last else 0.45,
                zorder=4 if is_last else 1)
        for z_seg, y_seg in _interp_gap(z_centers, pp_per_frame[k]):
            ax.plot(z_seg, y_seg, '--', color=color,
                    lw=1.6 if is_last else 0.5,
                    alpha=0.75 if is_last else 0.25,
                    zorder=4 if is_last else 0)
    if np.isfinite(pp_res):
        ax.axhline(pp_res, color='steelblue', ls='--', lw=1.2, alpha=0.9,
                   label=fr"$\langle p_p/P_{{\rm ext}}\rangle_{{\rm res}}={pp_res:.3g}$")
    if np.isfinite(pp_gel):
        ax.axhline(pp_gel, color='firebrick', ls='--', lw=1, alpha=0.7,
                   label=fr"$\langle p_p/P_{{\rm ext}}\rangle_{{\rm gel}}={pp_gel:.3g}$")
    ax.text(0.02, 0.04,
            fr"$\bar V_s = 1/\rho_{{s,\rm res}} = {V_bar:.3g}\,\sigma^3$"
            "\n"
            fr"$P_{{\rm ext}} = {p_ext:.3g}\,\varepsilon/\sigma^3$  (barostat target)"
            "\n"
            fr"$\Delta(p_p/P_{{\rm ext}})$(gel$-$res) $= {delta_pp:+.3g}$",
            transform=ax.transAxes, fontsize=9, va='bottom',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.85))
    ax.set_ylabel(r"$p_p / P_{\rm ext}$" "\n"
                  r"$= 1 + (\mu_{\rm total}-\mu_{\rm total,res})\,\rho_{s,\rm res}/P_{\rm ext}$")
    ax.set_xlabel(r"$z$  ($\sigma$)")
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.25)

    fig.savefig(out_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"Diagnostic plot   → {out_path}")
    print(f"  Δμ_ex(gel-res)    measured  = {delta_mu_ex_meas:+.4f} eps   "
          f"(equilibrium kT·ln(ρ_res/ρ_gel) = {delta_mu_ex_pred:+.4f} eps)")
    print(f"  Δμ_total(gel-res) measured  = {delta_mu_total:+.4f} eps   "
          f"(0.0 at equilibrium; |Δ| < 0.05 ε is flat)")
    print(f"  Δ(p_p/P_ext)(gel-res)       = {delta_pp:+.4f}")


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

def build_eps_map(eps_sp, eps_ss, piston_eps=1.0):
    """Return dict mapping atom type → epsilon for WCA with ghost solvent.

    piston_eps controls how the ghost test particle interacts with piston beads:
      1.0  — piston repels ghost (use for slab_with_support or permeation mode)
      0.0  — piston is transparent to ghost (use for compression mode, where the
             real solvent does NOT interact with the piston via WCA)

    Type 6 (walls, slab_with_flow only): WCA repulsion with solvent (eps=1.0).
    Absent from slab_with_support trajectories — harmless if type 6 never appears.
    """
    return {
        1: eps_sp,     # polymer backbone
        2: eps_sp,     # crosslinker
        3: eps_ss,     # solvent
        4: 0.0,        # support (no interaction with solvent)
        5: piston_eps, # piston — see docstring
        6: 1.0,        # walls — slab_with_flow lateral enclosure (WCA with solvent)
    }


def main():
    parser = argparse.ArgumentParser(
        description="Cavity-biased Widom insertion for LAMMPS trajectories.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--traj",        required=True,
                        help="LAMMPS dump file (id type x y z)")
    parser.add_argument("--n-bins",      type=int,   default=40,
                        help="Number of z-bins")
    parser.add_argument("--n-trial",     type=int,   default=20000,
                        help="Trial insertions per bin per frame")
    parser.add_argument("--r-cavity",    type=float, default=0.5,
                        help="Cavity radius in sigma units")
    parser.add_argument("--exclusion-buffer", type=float, default=None,
                        dest="exclusion_buffer",
                        help="Buffer (sigma) excluded around the piston extent "
                             "[z_piston_min-buf, z_piston_max+buf] and around "
                             "the support top.  Bins clearly above the piston "
                             "(reservoir in slab_with_flow) are kept active. "
                             "Defaults to r_cavity if omitted; set to ~one bin "
                             "width (e.g. 2.0) for slab_with_flow.")
    parser.add_argument("--temperature", type=float, default=1.0,
                        help="kT in LJ units")
    parser.add_argument("--eps-sp",      type=float, default=1.0,
                        help="Solvent–polymer epsilon")
    parser.add_argument("--eps-ss",      type=float, default=1.0,
                        help="Solvent–solvent epsilon")
    parser.add_argument("--piston-eps",  type=float, default=1.0,
                        dest="piston_eps",
                        help="Epsilon for ghost-particle WCA with piston beads. "
                             "Use 1.0 for slab_with_support or permeation mode "
                             "(piston repels solvent). Use 0.0 for compression "
                             "mode (piston is transparent to solvent). With 0.0 "
                             "the piston interior bins are physically meaningful "
                             "and reveal whether μ_total is continuous across the "
                             "piston.")
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
    parser.add_argument("--cavity-types", default="1,2,3,6",
                        help="Comma-separated atom types that count as cavity blockers "
                             "(default: '1,2,3,6' = polymer + solvent + walls). "
                             "Piston (type 5) and support (type 4) are excluded by "
                             "default so their lattice structure does not distort "
                             "void-fraction measurements. Walls (type 6) are included "
                             "as a defensive backstop: with the wall-interior trial-"
                             "point sampling (see --xy-extent / --wall-type), trial "
                             "points should never come near a wall plane anyway, but "
                             "this guards against partial-buffer setups.  If the "
                             "trajectory has no type-6 atoms (slab_with_support) the "
                             "extra entry is harmless.")
    parser.add_argument("--xy-extent", default=None,
                        help="Lateral sampling box for trial-point insertions, given "
                             "as four floats 'x_lo x_hi y_lo y_hi' (space- or comma-"
                             "separated, σ units, in world coordinates).  Overrides "
                             "the auto-detection below.  Use this if the wall type is "
                             "not 6 or if you want to restrict sampling to a polymer-"
                             "COM-tracked sub-region (e.g. mirroring reg_x / reg_y in "
                             "the LAMMPS chunk computes).")
    parser.add_argument("--xy-buffer", type=float, default=None,
                        dest="xy_buffer",
                        help="Inward buffer (σ) subtracted from the wall extent when "
                             "auto-detecting the lateral sampling box.  Default = "
                             "--cutoff (one WCA cutoff), which keeps trial points "
                             "from overlapping a wall bead.  Has no effect when "
                             "--xy-extent is specified or when no wall atoms exist.")
    parser.add_argument("--wall-type", type=int, default=6,
                        dest="wall_type",
                        help="Atom type used for the lateral wall tube (default 6, "
                             "matching slab_with_flow.lmp).  Set to a type that does "
                             "not appear in the trajectory to disable wall auto-"
                             "detection (full-box sampling), regardless of whether "
                             "type 6 is present.")
    parser.add_argument("--p-ext",       type=float, default=1.5,
                        help="LAMMPS barostat target pressure (eps/sigma^3). "
                             "Used to normalize the pore-pressure panel "
                             "(p_p / P_ext).  Default matches slab_with_support.lmp.")
    parser.add_argument("--no-plot",     action="store_true",
                        help="Skip the 4-panel diagnostic PNG (rho_s, mu_ex, "
                             "mu_total, p_p/P_ext).")
    parser.add_argument("--plot-path",   default=None,
                        help="Override path for the diagnostic PNG.  Default: "
                             "mu_total_diagnostic_{stem}.png in --out-dir.")
    args = parser.parse_args()

    traj_path = Path(args.traj)
    if not traj_path.exists():
        sys.exit(f"Trajectory not found: {traj_path}")

    out_dir = Path(args.out_dir) if args.out_dir else traj_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    stem         = args.out_stem if args.out_stem else traj_path.stem
    frame_path   = out_dir / f"mu_z_cavity_{stem}.dat"
    summary_path = out_dir / f"mu_z_cavity_summary_{stem}.dat"
    plot_path    = (Path(args.plot_path) if args.plot_path
                    else out_dir / f"mu_total_diagnostic_{stem}.png")

    eps_map      = build_eps_map(args.eps_sp, args.eps_ss, args.piston_eps)
    rng          = np.random.default_rng(args.seed)
    cavity_types = set(int(t) for t in args.cavity_types.split(','))

    # Parse --xy-extent: accept either comma- or whitespace-separated.  None
    # means auto-detect from `wall_type` atoms in each frame.
    if args.xy_extent is None:
        xy_extent = None
    else:
        tokens = args.xy_extent.replace(',', ' ').split()
        if len(tokens) != 4:
            parser.error("--xy-extent requires exactly four floats: "
                         "x_lo x_hi y_lo y_hi")
        xy_extent = tuple(float(t) for t in tokens)
        if xy_extent[0] >= xy_extent[1] or xy_extent[2] >= xy_extent[3]:
            parser.error("--xy-extent must satisfy x_lo<x_hi and y_lo<y_hi")

    excl_buf_display = args.exclusion_buffer if args.exclusion_buffer is not None else args.r_cavity
    print(f"Trajectory   : {traj_path}")
    print(f"z-bins       : {args.n_bins}")
    print(f"trials/bin   : {args.n_trial}")
    print(f"r_cavity     : {args.r_cavity} σ")
    print(f"excl_buffer  : {excl_buf_display} σ  "
          f"({'explicit' if args.exclusion_buffer is not None else 'default = r_cavity'})")
    print(f"kT           : {args.temperature}")
    print(f"eps_SP       : {args.eps_sp}   eps_SS : {args.eps_ss}   "
          f"piston_eps : {args.piston_eps} "
          f"({'phantom — piston interior bins active' if args.piston_eps == 0.0 else 'repulsive'})")
    print(f"WCA cutoff   : {args.cutoff} σ")
    print(f"cavity_types : {sorted(cavity_types)}  (piston/support excluded from void detection)")
    if xy_extent is None:
        xy_buf_print = args.xy_buffer if args.xy_buffer is not None else args.cutoff
        print(f"xy_extent    : auto (wall_type={args.wall_type}, buffer={xy_buf_print:.3f} σ; "
              f"falls back to full box if no walls)")
    else:
        print(f"xy_extent    : explicit  x:[{xy_extent[0]:.3f},{xy_extent[1]:.3f}] "
              f"y:[{xy_extent[2]:.3f},{xy_extent[3]:.3f}]")
    print(f"Frame file   : {frame_path}")
    print(f"Summary      : {summary_path}")
    if not args.no_plot:
        print(f"Plot         : {plot_path}")
        print(f"P_ext        : {args.p_ext} eps/sigma^3")
    print()

    all_results = []
    all_rho     = []   # per-frame rho_s(z) on the same bin grid as cavity-Widom
    all_steps   = []   # LAMMPS timesteps, for plot legend

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
                n_bins            = args.n_bins,
                n_trial           = args.n_trial,
                r_cavity          = args.r_cavity,
                temperature       = args.temperature,
                sigma             = args.sigma,
                cutoff            = args.cutoff,
                rng               = rng,
                cavity_types      = cavity_types,
                exclusion_buffer  = args.exclusion_buffer,
                xy_extent         = xy_extent,
                xy_buffer         = args.xy_buffer,
                wall_type         = args.wall_type,
            )

            # ── Solvent density on the SAME z-bins as cavity-Widom ──────────
            # Used to form mu_total(z) = mu_ex(z) + kT*ln(rho_s(z)) in the
            # diagnostic plot.  Cheap relative to the cavity insertions.
            rho_z = bin_solvent_density(atoms_xyz, atom_types, box, args.n_bins)

            write_frame(fh_frame, timestep, bin_results)
            all_results.append(bin_results)
            all_rho.append(rho_z)
            all_steps.append(int(timestep))

            n_nan    = sum(1 for b in bin_results if np.isnan(b['mu_ex']))
            mean_pc  = np.mean([b['p_cav'] for b in bin_results])
            print(f"done  (p_cav_mean={mean_pc:.3f}, {n_nan}/{args.n_bins} bins NaN)")

    write_summary(str(summary_path), all_results, args.n_bins)

    print(f"\nWrote {len(all_results)} frames → {frame_path.name}")
    print(f"Summary → {summary_path.name}")

    # ── Diagnostic plot ────────────────────────────────────────────────────
    if args.no_plot:
        print("[diagnostic plot] --no-plot set; skipping")
    elif len(all_results) == 0:
        print("[diagnostic plot] no frames processed; skipping")
    else:
        try:
            make_diagnostic_plot(
                all_results   = all_results,
                rho_per_frame = np.array(all_rho),
                steps         = all_steps,
                kT            = args.temperature,
                p_ext         = args.p_ext,
                out_path      = plot_path,
                title         = stem,
            )
        except Exception:
            # Don't let a plotting bug throw away the simulation's .dat output.
            # Print the full traceback so the cause is visible in the SLURM log.
            import traceback
            print("[diagnostic plot] FAILED — traceback follows; .dat outputs are unaffected")
            traceback.print_exc()

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
