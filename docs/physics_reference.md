# Physics & Units Reference

[← back to README](../README.md)

## LJ units quick reference

Full details are in [`lj_units_cheat_sheet.md`](../lj_units_cheat_sheet.md). Key conversions for PEG/water:

| Quantity | Multiply LJ value by | To get |
|----------|---------------------|--------|
| Length | 0.76 | nm |
| Time (τ) | 4.6 | ps |
| Pressure (P*) | 9.4 | MPa |
| Force (F*) | 5.4 | pN |
| Temperature (T*) | 300 | K (at T*=1) |
| 10⁶ steps at dt=0.005 | → | ~23 ns |

**Current simulation parameters:**
- Temperature: T* = 1.0 (= 300 K)
- Pressure (NPT target): P* = 1.5
- Timestep: dt = 0.005 τ
- FENE bonds: K=30, R₀=1.5, ε=1, σ=1 (Kremer–Grest standard)
- Pair interactions: WCA (rc = 1.122σ, purely repulsive) — consistent across all active simulation types

---

---

## Barostat choice (equilibration)

Match the barostat to the geometry:

- **Free-swelling gel in a solvent bath** (`slab_with_support` equilibration): use `fix npt … aniso P P pdamp`. Each of x, y, z is barostatted independently to the target pressure, so the box adopts whatever aspect ratio balances σxx = σyy = σzz = P and the gel relaxes to its own equilibrium shape — the analogue of a hydrogel free to swell in all directions. **Do not use `iso` here:** `iso` controls only the mean (hydrostatic) pressure and freezes the box aspect ratio, so any anisotropic stress or z-padding baked into the data file is never relaxed. Switching `slab_with_support` from `iso` → `aniso` (2026-06-24) fixed exactly this: the gel had been stuck artificially swollen along z, and with `aniso` it reaches a noticeably taller, true equilibrium swelling. Use `couple xy` only if you must enforce lateral isotropy (a free gel reaches it anyway); use `tri` only to relax shear stress (lets the box tilt).
- **Piston-driven runs** (`triaxial_compression`, `triaxial_permeation`): production is `fix nvt` with the box fixed and the piston as the sole z-actuator — **no box barostat in production** (a barostat would double-control z and fight the piston). The only barostat is the Phase-0.5 pre-equilibration, run with the piston/support temporarily on `nve`+`setforce` so they scale with the box.

  As of **2026-07-29** this Phase-0.5 barostat is **`fix nph z`** (z-only), replacing the earlier `aniso`/`iso` forms. Rationale: transverse (xx, yy) stresses build up in the polymer network during compression or permeation; a scalar (`iso`) or per-axis (`aniso`) barostat would let those transverse stresses perturb the box and drift the reservoir pressure. Barostatting **z only** targets the zz stress component directly, holding the solvent reservoir at Pzz = P* = 1.5 while x, y box dimensions stay fixed at the periodic slab's equilibrium extent. Because the target is now Pzz (not the full scalar `Press`), the old `+0.41` kinetic offset used by `slab_with_flow`'s `iso` convention is dropped — all three scripts target `npt_P05_target = P_target`.

  These three scripts run on the **periodic** slab geometry (`slab_with_support_periodic.ipynb`): laterally periodic (`p p p`), no side walls, one support+piston sheet per z-period. Because the slab fills the x-y plane, x and y are free to relax and the piston stress is simply `c_piston_fz / (lx*ly)` with no multi-period normalisation.

When adding lateral walls to an equilibrated config, `add_walls_to_slab.ipynb` unwraps the gel via image flags before measuring its extent, so a gel that has drifted across a periodic face (a little "image pollution" in the visualizer) does not corrupt the wall placement. Under `boundary p p p` the wrapped sliver itself is harmless to the run (bonds use the minimum image; `compute com` unwraps). With `RECENTER_LATERAL = True` (default) it also shifts the mobile group (gel + solvent) in x/y so the polymer COM lands at the lateral box center — the support/piston plates and the z-axis are left untouched, which keeps gel↔support/piston contact along the loading axis and relies on the plates being laterally larger than the gel (the notebook checks this and warns if the gel would exceed the support footprint).

---

*Last updated: 2026-07-29. For questions, contact Dylan Pollard (pollard@ucsb.edu).*

## Network stress, M and G from a uniaxial compression

The gel is a two-phase (poroelastic) material, so the measured total stress in each direction splits as (Terzaghi)

    σᵗ_ii(z) = σ′_ii(z) + p_pore,        ii = zz, xx, yy

where p_pore is the solvent (pore) pressure and σ′ is the **network** (effective) stress carried by the polymer. In `triaxial_compression` and `triaxial_permeation` the polymer + solvent partial stress profiles are written for all three diagonal components, and `scripts/lib/triaxial.py` treats every component the same way: p_pore for that component is read from the flat far-reservoir window of the **same** component's total profile (the reservoir fluid is isotropic, so the three baselines agree to noise), and σ′_ii is the total minus that baseline. At ε = 0 mechanical equilibrium with the bath forces σ′_zz ≈ 0; the lateral σ′_xx, σ′_yy need not vanish because the periodic box fixes l_x, l_y, so the moduli below use **increments** relative to the reference state.

Uniaxial strain ε along z with the lateral box fixed (ε_xx = ε_yy = 0), isotropic drained network with Lamé constants λ, G and longitudinal modulus M = λ + 2G:

    σ′_zz = M ε,        σ′_xx = σ′_yy = λ ε = (M − 2G) ε

so

    M = σ′_zz / ε,      σ′_zz / σ′_xx = M / (M − 2G),      G = (σ′_zz − σ′_xx) / (2ε)

(Note the ratio is M/(M − 2G), **not** M/(M − 2G/3): λ = K − 2G/3 and M = K + 4G/3, so M − 2G = λ.) `M_network` uses the plateau network profile averaged over the membrane; `M_piston` = ⟨P⟩/ε from the block-bootstrapped piston force is the independent check. G is formed once from xx and once from yy; by symmetry the two must agree, and their spread is a second error estimate. The cooperative diffusivity D_c from the consolidation fit of u_z(z,t) then gives the hydraulic permeability κ = D_c/M (κ = k/η). All of this is implemented in `scripts/lib/triaxial.py` and drawn by `triaxial_compression_{single,sweep}.ipynb` (see `docs/analysis.md` §7b).
