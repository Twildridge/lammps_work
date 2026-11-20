# Random 10-bead gel: plan, criteria and sizing (2026-11-14)

> AI-generated working note (Claude), like the rest of `docs/`. Numbers marked *proposed* are
> suggestions for the author to accept or change; numbers marked *measured* come from the screen
> of 2026-11-13/14 (`lammps_runs/chain_length_testing/screen_report.md`).

## Goal

An osmotically active gel with no permanent pores, sized so that stress profiles in uniaxial
compression, uniaxial permeation and simple shear are quiet, at a cost per measurement no higher
than the original tetrahedral slab on the old settings (4 nodes, comm cutoff 5.0).

## Phase 0: acceptance criteria (*proposed*)

| property | measure | proposed threshold | random 10 today (*measured*) |
|---|---|---|---|
| no permanent pores | Phi plateau at r_p = 1.0, 1.5 (`screen_voids.py`), flat out to 10 M steps | <= 0.05, the rot45 lattice's value at r_p = 1.0 | 0.08 / 0.09 |
| in-plane isotropy | r = G1/G2 (`screen_isotropy.lmp`) | within 10 % of 1 | 1.04 ± 0.05 (whole box) |
| out-of-plane isotropy | G from xx and from yy in the compression run agree with each other and with G1 | within 10 % | not measured yet |
| osmotically active | chord extension <R/(n b)>; G/M from the compression run | extension <= 0.55; G/M to be set once M is known | extension 0.53; M pending |
| equilibrated | 2 M-step Rg means and the box flat | within 0.15 % over 8 M steps | yes |
| cost | steps/s of the two-piston file on 6 nodes | not below the old tetrahedral two-piston rate (0.55-0.75 M steps/h) | to benchmark |

## Phase 1: sizing (first pass, from existing runs; no new MD)

*Measured* noise of the zero-strain stage of the Part B runs, whole-box stress per gel volume:
the standard deviation of a stress difference averaged over a volume V (sigma^3) and a time t
(tau) is S / sqrt(V t), nearly white from 1k to 100k steps.

| gel | S for sigma_xy | S for sigma_xx - sigma_yy and sigma_zz - sigma_xx |
|---|---|---|
| random 10 | 1.3 - 2.1 | 2.5 - 3.1 |
| random 5 | 1.9 - 2.4 | 2.9 - 6.2 |
| rot45 5 | 1.4 - 1.5 | 1.8 - 4.8 |

Consequences for random 10 (G about 0.04, so the anisotropy signal at 10 % compression is
2 G eps = 0.009):

- Gel-average quantities are not noise-limited: the present slab reaches 3 % of that signal in
  100k steps.
- A z-profile bin (2 sigma thick) on the present 47 x 47 cross-section has 6e-4 of noise per
  1 M steps, 7 % of the signal. Per-bin noise falls as 1 / sqrt(area x time).
- Hold lengths are set by the relaxation time tau_1 = L^2 / (pi^2 D_c), not by noise, and scale
  with thickness squared.

So the gel should be wider and thinner. With N_atoms about A (0.58 L + 82) for the two-piston
file (checked against the rot45 5-bead file to 9 %):

| candidate (x by y by thickness) | two-piston atoms | hold time vs today | per-bin noise at equal time | cost per measurement vs today | thickness in mesh sizes |
|---|---|---|---|---|---|
| today, 47 x 47 x 122 | 0.34 M | 1 | 1 | 1 | 27 |
| 70 x 70 x 80 | 0.63 M | 0.43 | 0.68 | 0.79 | 17 |
| 85 x 85 x 65 | 0.87 M | 0.28 | 0.56 | 0.72 | 14 |
| 95 x 95 x 55 | 1.03 M | 0.20 | 0.50 | 0.61 | 12 |

Cost ceiling: the two-piston decks run at about 1.9e8 atom-steps/s on 6 nodes now, against about
5e7 for the tetrahedral two-piston file on the old layout, so about 1.0 M atoms steps as fast as
the original did. All of these are estimates; phase 3 (one benchmark point) confirms the rate and
the first compression run on the present slab gives D_c and M, which fix the hold lengths.

Open choice: the thickness floor. The surface layers (the trimmed 3 sigma plus one to two mesh
sizes per side) leave about 23 interior bins at 65 sigma and 18 at 55 sigma.

## Phases 2-9

See the table in the conversation of 2026-11-14 (build with `build_beads_sweep.py --beads 10
--box-xy ... --thickness ...`, benchmark, equilibrate in `chain_length_testing`, qualify with the
`screen_*` tools, convert, then compression, permeation and shear on the same `_two_pist_pm40`
file, closure in `cross_deck_comparison.ipynb`).

## Running now on the existing random 10 slab

`simulations/chain_length_testing/r10_measure_bridge.batch` converts equilibration snapshot 12
(step 18.01 M; its box is the closest to the 16-24 M mean) and submits compression, periodic shear
and the frozen-at-reference permeation. Job IDs: `lammps_runs/chain_length_testing/r10_measurements_20261114/jobs.txt`.
The live permeation drive waits for M and D_c from the compression run (dP about 0.24 M, drive
about 4 tau_1).

## Ways to lower the frozen-pore fraction without giving up modulus

*Measured* trend with strand length: G falls fast (0.38, 0.043, 0.018, 0.010, 0.004 for 5, 10,
15, 20, 30 beads) while the Phi plateau drops from 0.32 to about 0.08 between 5 and 10 beads and
then barely moves. Longer strands past 10 beads cost modulus and buy almost nothing in pores.

1. **Intermediate strand length (7 or 8 beads).** The whole drop in Phi happens between 5 and 10
   beads; a strand of 7-8 beads may keep G near 0.1 with Phi near the 10-bead value. Cheapest to
   try: the builder already takes `--beads`.
2. **Interpenetrating networks.** Two or three independent random f = 4 networks of 10- to
   15-bead strands in the same volume. Strand density, and with it G, scales with the number of
   networks, and each network's mesh voids are filled by the others, so the pore scale is set by
   the concentration rather than by one network's mesh. This is how gels made above the overlap
   concentration avoid mesh-sized pores. Needs a builder change (overlay shifted networks).
3. **A more uniform junction placement.** The plateau is the quenched non-uniformity of the
   bond-switched network. Screen builder seeds and switch counts with a static measure before
   any MD (the variance of the crosslinker count in cells of about two mesh sizes) and equilibrate
   only the most uniform one; or place junctions from a relaxed repulsive-particle configuration
   before wiring them.
4. **Less swelling.** A higher bath pressure or a slightly poorer solvent gives a denser gel with
   smaller voids and a higher modulus, at the price of a weaker osmotic response.

Nothing in this list has been tested.

## Status, 2026-11-14 evening: the wide gel pipeline

Submitted with `bash simulations/chain_length_testing/wide10_pipeline.sh`; status with
`bash simulations/chain_length_testing/wide10_status.sh`; log and job IDs in
`lammps_runs/chain_length_testing/wide10_pipeline_20261114/`.

1. `wide10_build.batch`: 10 x 10 x 10 cells (about 78 x 78 sigma, 78 thick as built), the same
   builder recipe as the qualified tall 10-bead slab -> `slab_support_periodic_10beads_wide_rho04_random.data`.
2. `chain_length_testing_wide10.batch`: aniso NPH, Rg halt (1 M-step blocks, 16 M budget), and the
   new mean-box ending of `slab_with_support.lmp` (final configuration at the time-averaged lx, ly).
3. `wide10_bridge.batch`: void persistence on the equilibration tail (label `w10`), a gate
   (Phi <= 0.15 at r_p 1.0 and 1.5, box within 5 % of square), two-piston conversion, then
   in-plane isotropy, compression (hold 4 M + unload 7 M), periodic shear (hold 3 M) and the
   frozen-reference permeation. The live permeation drive waits for M and D_c.

Density-uniformity term in the builder (`www_bond_switch` `w_uniform`, preset `random_u`): first
test negative. On a 6 x 6 x 8 block the polymer-density CoV on 8.7 sigma voxels was 0.088 /
0.084 / 0.079 / 0.078 for weights 0 / 0.3 / 1 / 3 (lattice 0.035) and the largest empty sphere
2.9-3.0 sigma in all four; the greedy clean-up pass after the anneal rewires the long strands
without the term. It is therefore OFF in the wide build. A second comparison with a three times
longer anneal is queued (`builder_uniformity_test.batch`, logs in
`lammps_runs/chain_length_testing/builder_uniformity_20261114/`).

## Status, 2026-11-16: wide gel measurements

All four jobs of the bridge finished on 2026-11-16. Three are valid; the compression is not.

| Measurement | Job | Status | Needs M / D_c from the compression? |
|---|---|---|---|
| In-plane isotropy screen | 55745205 | valid (G_x 0.0346, G_y 0.0338) | no |
| Periodic shear | 55745207 | valid (flat hold, G about 0.040 polymer) | no |
| Frozen-reference permeation (`FROZEN_AT=ref`) | 55745208 | valid | not to run: fixed dP 0.1 through the network frozen undeformed. Only the later comparison of its kappa with D_c/M needs them |
| Compression + unload | 55745206 | NOT usable, rerun | it is the source of M and D_c |
| Live permeation drive | not submitted | waits | yes: dP about 0.24 M, drive about 4 tau_1 (set by the batch from the compression record since 2026-11-20, `scripts/lib/perm_plan.py`) |

- **Why the compression failed.** The gel is equilibrated, but its surface is diffuse: `L_bb - L_rg` is 13-16 sigma (attached
  loops, no dangling ends; rot45 lattices 4-8 sigma). The plates seat on the outermost bead, so the 9 sigma closure at strain
  0.10 went into that fringe: Rg strain 2.5 %, stress 0.0045 +/- 0.0007, no plateau.
- **Rerun.** `SEAT_MODE=rg UNLOAD_MARGIN=8` with the bridge's `NSTEPS_OVERRIDE=4000000 UNLOAD_STEPS=7000000` on
  `uniaxial_compression_single.batch` (plates seated on the Rg face; strain still measured on the bounding box; see
  "SEAT_MODE and the credit guards" in `running_simulations.md`). About 7,600 SU; the guards stop it early if the seating
  pre-strains the network by more than 0.10 or the Rg strain stays below half the applied strain.
  Submitted 2026-11-16 as job 55793779 (`triax_comp2p_w10_rg`, 16 h limit, unload/re-swell on); files `_4000000_c0.10` and `_u0.10`.
- **One permeation deck, two frozen modes.** `FROZEN_AT=end` (the default): the live drive, then the same flow-compacted network
  frozen under the same dP, both in one job; kappa_frozen / kappa_live separates what a rigid network of that shape gives from what
  the network's own motion adds. `FROZEN_AT=ref`: frozen undeformed for the whole production, no live drive (the zero-strain kappa,
  to compare with D_c/M). Only the `ref` run has been done for this gel; the live drive will use the default.
  Since 2026-11-16 the default sequence has a velocity-controlled stage between the drive and the frozen stage (feed piston at
  the drive's mean velocity, dP measured), and the frozen stage keeps that velocity instead of the applied dP; see "Permeation
  sequence since 2026-11-16" in `running_simulations.md`. The converted file needs feed-reservoir room for both stages.
- **Code status.** `seat_mode` and the three guards (compression deck and `isotropic_compression`, commit e9f0e86) have run only on a toy
  system: smoke-test both decks locally before the rerun. The permeation deck gained the velocity-controlled stage on 2026-11-16 (smoke-tested locally, not yet run on Expanse).

