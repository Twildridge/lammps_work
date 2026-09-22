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
- **Piston-driven runs** (`triaxial_compression`, `triaxial_permeation`): production is `fix nvt` with the box fixed and the two plates (piston down, support up — symmetric drive since 2026-09-18) as the sole z-actuators — **no box barostat in production** (a barostat would double-control z and fight the piston). The only barostat is the Phase-0.5 pre-equilibration, run with the piston/support temporarily on `nve`+`setforce` so they scale with the box.

  As of **2026-07-29** this Phase-0.5 barostat is **`fix nph z`** (z-only), replacing the earlier `aniso`/`iso` forms. Rationale: transverse (xx, yy) stresses build up in the polymer network during compression or permeation; a scalar (`iso`) or per-axis (`aniso`) barostat would let those transverse stresses perturb the box and drift the reservoir pressure. Barostatting **z only** targets the zz stress component directly, holding the solvent reservoir at Pzz = P* = 1.5 while x, y box dimensions stay fixed at the periodic slab's equilibrium extent. Because the target is now Pzz (not the full scalar `Press`), the old `+0.41` kinetic offset used by `slab_with_flow`'s `iso` convention is dropped — all three scripts target `npt_P05_target = P_target`.

  These three scripts run on the **periodic** slab geometry (`slab_with_support_periodic.ipynb`): laterally periodic (`p p p`), no side walls, one support+piston sheet per z-period. Because the slab fills the x-y plane, x and y are free to relax and the piston stress is simply `c_piston_fz / (lx*ly)` with no multi-period normalisation.

When adding lateral walls to an equilibrated config, `add_walls_to_slab.ipynb` unwraps the gel via image flags before measuring its extent, so a gel that has drifted across a periodic face (a little "image pollution" in the visualizer) does not corrupt the wall placement. Under `boundary p p p` the wrapped sliver itself is harmless to the run (bonds use the minimum image; `compute com` unwraps). With `RECENTER_LATERAL = True` (default) it also shifts the mobile group (gel + solvent) in x/y so the polymer COM lands at the lateral box center — the support/piston plates and the z-axis are left untouched, which keeps gel↔support/piston contact along the loading axis and relies on the plates being laterally larger than the gel (the notebook checks this and warns if the gel would exceed the support footprint).

---

*Last updated: 2026-07-29. For questions, contact Dylan Pollard (pollard@ucsb.edu).*

## Network stress, M and G from a uniaxial compression

The gel is a two-phase (poroelastic) material, so the measured total stress in each direction splits as (Terzaghi)

    σᵗ_ii(z) = σ′_ii(z) + p_pore,        ii = zz, xx, yy

where p_pore is the solvent (pore) pressure and σ′ is the **network** (effective) stress carried by the polymer. In `triaxial_compression` and `triaxial_permeation` the polymer + solvent partial stress profiles are written for all three diagonal components, and `scripts/lib/triaxial.py` treats every component the same way: p_pore for that component is read from the flat far-reservoir window of the **same** component's total profile (the reservoir fluid is isotropic, so the three baselines agree to noise), and σ′_ii is the total minus that baseline. In the two-piston geometry that window is the feed reservoir between the gel top and the feed piston, rebuilt per snapshot from the measured piston plane and restricted to bins lying entirely ≥ 3 σ (`Config.res_wall_margin`) from the piston sheet: the solvent next to a piston is depleted/layered and half of the wall pair virial is booked on the piston atoms (not in the solvent profile group), so a bin touching the piston reads σ_zz ≈ 1.43 for a 1.50 bath and would bias every σ′ by that much. At ε = 0 mechanical equilibrium with the bath forces σ′_zz ≈ 0; the lateral σ′_xx, σ′_yy need not vanish because the periodic box fixes l_x, l_y, so the moduli below use **increments** relative to the reference state.

Uniaxial strain ε along z with the lateral box fixed (ε_xx = ε_yy = 0), isotropic drained network with Lamé constants λ, G and longitudinal modulus M = λ + 2G:

    σ′_zz = M ε,        σ′_xx = σ′_yy = λ ε = (M − 2G) ε

so

    M = σ′_zz / ε,      σ′_zz / σ′_xx = M / (M − 2G),      G = (σ′_zz − σ′_xx) / (2ε)

(Note the ratio is M/(M − 2G), **not** M/(M − 2G/3): λ = K − 2G/3 and M = K + 4G/3, so M − 2G = λ.) `M_network` uses the plateau network profile averaged over the membrane; `M_piston` = ⟨P⟩/ε from the block-bootstrapped piston force is the independent check. G is formed once from xx and once from yy; by symmetry the two must agree, and their spread is a second error estimate. The cooperative diffusivity D_c from the consolidation fit of u_z(z,t) then gives the hydraulic permeability κ = D_c/M (κ = k/η). All of this is implemented in `scripts/lib/triaxial.py` and drawn by `triaxial_compression_{single,sweep}.ipynb` (see `docs/analysis.md` §7b).


## NPT-piston reservoirs (two-piston sequence)

*(2026-09-16; `triaxial_permeation_two_pist`, `triaxial_compression_two_pist`, converter `slab_two_pistons.ipynb`.)*

**Piston law.** Marioni et al., *J. Membr. Sci.* 738 (2026) 124837 ("Non-equilibrium simulations of hydraulic
permeation: role of mechanical boundary conditions in dense membranes") hold each solvent reservoir at a prescribed
pressure with a mobile, force-loaded, damped sheet — their Eq. 3, for a sheet of N beads and area A = l_x l_y:

    M_p  dv/dt  =  F_fluid  −  P_i A  −  C v          (per sheet; v = sheet velocity along z)

The LAMMPS port applies it per bead so the sheet translates rigidly (no `fix rigid`: the full-width periodic sheet
is wider than half the box and `rigid` mis-reconstructs it):

    fix <id>_xy  <piston> setforce 0.0 0.0 NULL
    fix <id>_z   <piston> aveforce NULL NULL v_fz_<piston>
    fix <id>_nve <piston> nve
    fz_feed = −P_feed·l_x l_y / N − C_pist · vcm(piston_feed, z)      (top sheet pushes DOWN)
    fz_perm = +P_perm·l_x l_y / N − C_pist · vcm(piston_perm, z)      (bottom sheet pushes UP)

`aveforce` gives every bead the group-average force plus `fz`, so the applied load is −P A on the whole sheet and
the damping on the sheet is N·C_pist·v.  The pistons are not thermostatted (as in the paper); the viscous term
removes their thermal energy.  The bead mass is a deck variable (`mass 5/6/7 ${piston_mass}`, default 1000).

**Why l_x, l_y are fixed and why a converter.** The scheme regulates pressure in z only.  It cannot swell a fresh
lattice laterally, so the two-piston data file is *converted* from the equilibrated aniso-NPH slab
(`slab_with_support`, piston transparent to solvent, σ_p,xx/σ_p,zz = 1.0005): l_x, l_y, the gel dimensions and the
stress-free state are inherited exactly, and the NPT-piston phase in the decks is a z-settle of the reservoirs.
`boundary p p p` is kept: nothing can cross z (wet pistons are WCA walls to solvent and polymer, every wall–wall
pair is off, the margins are vacuum), so the full-box z binning still works.  Thermo `press` is meaningless (vacuum
in V); reservoir pressures come from the pistons and the stress profiles.

**Reservoir pressure readout.** Not `compute reduce sum fz` on a piston — by the time thermo/print evaluates it,
`setforce`/`aveforce` have already modified `f` and the net is ~0 at steady state.  Use the pair force of the
mobile atoms on the sheet: `compute fp_feed piston_feed group/group mobile`, P = ±c_fp[3]/(l_x l_y).

**Critical damping.** The paper's C = 500 is for a few-hundred-atom graphene sheet; with ~23,316 beads and m = 1000
it would be ~50× overdamped (relaxation ~3.6 M steps).  Each deck therefore models the piston on its solvent column
as an oscillator and defaults to critical damping:

    k = K_solv · l_x l_y / L_res,   ω = √(k / (N m)),   period = 2π/ω,   C_crit (per bead) = 2 √(k m / N),
    C_pist = c_pist_frac · C_crit

with `K_solv` the solvent bulk modulus at P* = 1.5 (deck default 10; no measured value is recorded in the repo — a
Carnahan–Starling estimate for the WCA fluid at ρ ≈ 0.43–0.47 gives ≈ 3.5, i.e. C_crit ≈ 1.7× smaller; either way
the settle is critically-to-mildly overdamped).  For the rho04 slab (A ≈ 2077, L_feed ≈ 18, L_perm ≈ 10) the deck
prints ω ≈ 7–9 × 10⁻³/τ, periods ≈ 1.4–1.8 × 10⁵ steps and C_crit ≈ 14–18 per bead — so `NPT_PISTON_STEPS = 1 M`
is ≈ 5–7 periods.

**Modes.** Permeation: P_feed = P_target + dP, P_perm = P_target; flux Q_perm = A·dz_perm/dt from the permeate
piston (bead count crossing the support kept as a cross-check); k = Q_perm L/(A dP).  Compression: both wet pistons
at P_target (drained consolidation at constant bath pressure) and a third, solvent-transparent load piston inside the
feed reservoir loads the network through the usual strain sweep; solvent expelled by the compression raises the
feed piston (and lowers the permeate piston), so the converter's `margin_feed` must cover the deepest strain.
