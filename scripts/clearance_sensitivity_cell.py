# =============================================================================
# Clearance-sensitivity check for ΔV_mix  (standalone notebook cell)
# =============================================================================
# Diagnostic only — does NOT touch the pipeline.
#
# V_mix is currently the geometric bounding box that isolate_gel.py draws around
# the polymer: (polymer extent) + BOX_CLEARANCE on every face. That clearance is
# a *fixed* 0.2σ per face (0.4σ per axis). As the gel compresses with P*, the
# fixed margin becomes a larger FRACTION of a shrinking box, which biases
# ΔV_mix upward with pressure.
#
# This cell recomputes ΔV_mix while shrinking the clearance from the original
# 0.2 down to 0, using the SAME pure-phase references (V_pol, V_sol unchanged).
# If the upward slope mostly flattens as clearance → 0, the trend was the
# clearance artifact. If it persists, it's the real connectivity/packing effect
# and the full NPT-V_mix restructure is what's needed to confirm.
#
# Assumes the notebook has already defined (from the earlier cells):
#   PRESSURES, BASE_DATANAME, INTERACTION, SLAB_STEPS, DATA_DIR,
#   PURE_INTER, PURE_STEPS, load_manifest, avg_box_volume, np, plt
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ORIG_CLEARANCE = 0.2          # BOX_CLEARANCE used in isolate_gel.py (per face)
TEST_CLEARANCES = [0.2, 0.15, 0.10, 0.05, 0.0]   # per face; 0.2 == current pipeline


def _read_box_dims_from_data(path):
    """Return (Lx, Ly, Lz) from a LAMMPS data-file header, or None if not found."""
    lo_hi = {}
    try:
        with open(path) as f:
            for line in f:
                for axis in ("x", "y", "z"):
                    if f"{axis}lo {axis}hi" in line:
                        p = line.split()
                        lo_hi[axis] = (float(p[0]), float(p[1]))
                if len(lo_hi) == 3:
                    break
    except FileNotFoundError:
        return None
    if len(lo_hi) != 3:
        return None
    return tuple(hi - lo for lo, hi in (lo_hi["x"], lo_hi["y"], lo_hi["z"]))


def _find_isolated_dims(p_dir, isolated_dataname):
    """
    Try to get the per-axis box of the isolated mixed gel (exact deflation).
    Falls back to None -> caller uses the cube-root approximation.
    Searches a few plausible local locations for isolated_*.data.
    """
    candidates = [
        p_dir / f"{isolated_dataname}.data",
        p_dir / f"{isolated_dataname}.lammps",
        # add your own local path here if the isolated files live elsewhere:
        # Path.home() / "Documents/lammps_data/slab_with_support" / f"{isolated_dataname}.data",
    ]
    for c in candidates:
        dims = _read_box_dims_from_data(c)
        if dims:
            return dims
    return None


def _deflate_volume(L_box_axes, V_box_scalar, orig_clear, new_clear):
    """
    Recompute V_mix for a new clearance.
    Occupied extent per axis = L_box - 2*orig_clear; new box adds 2*new_clear.
    If per-axis dims are unknown, fall back to an isotropic (cube-root) estimate,
    which UNDER-counts the effect for a slab (thin axis loses a larger fraction).
    """
    delta = 2.0 * (orig_clear - new_clear)   # subtract this length from each axis
    if L_box_axes is not None:
        return float(np.prod([L - delta for L in L_box_axes])), "exact"
    L = V_box_scalar ** (1.0 / 3.0)
    return float((L - delta) ** 3), "approx(cube)"


rows = []
methods_used = set()

for P in PRESSURES:
    pstr = f"{P:.1f}"
    dataname          = f"{BASE_DATANAME}_pstar{pstr}"
    isolated_dataname = f"isolated_{dataname}_{INTERACTION}_{SLAB_STEPS}"
    sol_dataname      = f"{isolated_dataname}_solvent_only"
    pol_dataname      = f"{isolated_dataname}_polymer_only"
    p_dir = DATA_DIR / f"p{pstr}"

    manifest_file = p_dir / f"{isolated_dataname}_volmix_manifest.json"
    sol_vol_file  = p_dir / f"box_dimensions_{sol_dataname}_{PURE_INTER}_{PURE_STEPS}.dat"
    pol_vol_file  = p_dir / f"box_dimensions_{pol_dataname}_{PURE_INTER}_{PURE_STEPS}.dat"
    if not all(f.exists() for f in [manifest_file, sol_vol_file, pol_vol_file]):
        print(f"[SKIP] P*={pstr}: missing manifest/solvent/polymer file")
        continue

    mf        = load_manifest(manifest_file)
    V_box0    = mf["V_mix_isolated"]                       # geometric box at clearance 0.2
    scale_pol = mf["scale_pol"]
    scale_sol = mf["scale_sol"]
    V_pol_ref = avg_box_volume(pol_vol_file, skip_frac=0.5) * scale_pol
    V_sol_ref = avg_box_volume(sol_vol_file, skip_frac=0.5) * scale_sol
    V_ref     = V_pol_ref + V_sol_ref

    L_axes = _find_isolated_dims(p_dir, isolated_dataname)

    row = {"P": P, "V_ref": V_ref}
    for c in TEST_CLEARANCES:
        V_mix_c, method = _deflate_volume(L_axes, V_box0, ORIG_CLEARANCE, c)
        methods_used.add(method)
        row[f"c{c}"] = (V_mix_c - V_ref) / V_ref       # fractional ΔV_mix
    rows.append(row)

if not rows:
    print("No data — check that the staged manifest/box_dimensions files exist.")
else:
    P_arr = np.array([r["P"] for r in rows])

    # --- Plot: ΔV_mix/(V_sol+V_pol) vs P* for each clearance ---
    fig, ax = plt.subplots(figsize=(7, 5))
    print(f"\n{'P*':>5}", *[f"c={c:<5}" for c in TEST_CLEARANCES], sep="  ")
    for c in TEST_CLEARANCES:
        y = np.array([r[f"c{c}"] for r in rows])
        slope = np.polyfit(P_arr, y, 1)[0]
        ax.plot(P_arr, y, "o-", lw=1.8, ms=6,
                label=f"clearance={c}  (slope={slope:+.4f}/P*)")
    ax.axhline(0, color="gray", lw=1, ls="--")
    ax.set_xlabel(r"$P^*$")
    ax.set_ylabel(r"$\Delta V_{\rm mix}\,/\,(V_{\rm sol}+V_{\rm pol})$")
    ax.set_title("Clearance sensitivity of ΔV_mix")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    # --- Table ---
    print()
    for r in rows:
        vals = "  ".join(f"{r[f'c{c}']:+.4f}" for c in TEST_CLEARANCES)
        print(f"{r['P']:>5.1f}  {vals}")

    print(f"\nDeflation method(s) used: {sorted(methods_used)}")
    if "approx(cube)" in methods_used:
        print("NOTE: cube-root fallback under-counts the clearance effect for a "
              "slab geometry.\n      For exact per-axis deflation, make the "
              "isolated_*.data files reachable\n      (edit the candidate paths in "
              "_find_isolated_dims, or stage them from Expanse).")
    print("\nRead-out: if the slope shrinks toward ~0 as clearance → 0, the upward "
          "trend\nwas the fixed-clearance artifact. If it survives, it's the real "
          "packing effect.")
