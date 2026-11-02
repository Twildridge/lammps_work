#!/usr/bin/env python3
"""periodic_data_file_generator.py -- the x/y-PERIODIC gel slab builder (2026-11-02; moved here
from scripts/slab_with_support_periodic.ipynb, which now only sets the knobs, the way
data_file_generator.py / data_file_generator.ipynb split on 2026-10-28).

The slab spans the whole x-y plane and its network is crosslinked across the x and y
boundaries (run with `boundary p p p`); in z every period is a finite slab with its own
solvent layer, frozen support sheet (type 4) and parked piston sheet (type 5).  The
reference 14000002 state (slab_support_periodic_5beads_tall_rho04_new) comes from the
'tetrahedral' build; the two-piston converter (slab_two_pistons_periodic.py) turns the
equilibrated final configuration into the triaxial_*_two_pist input.

Two networks, chosen by `network` in build():
    'tetrahedral'  the diamond lattice with the cube axes along the box axes (every periodic
                   slab up to 2026-11-01);
    'random'       the same lattice re-wired into a random f = 4 network by Wooten-Winer-Weaire
                   bond switching (data_file_generator.www_bond_switch with periodic=(box_x,
                   box_y, 0): strands and repulsions across the x/y boundaries are minimum-image,
                   only the two z faces are tethered).  Same recipe and knobs as the isolated
                   'random' gel (3 accepted switches per strand, seed 7); the start is the
                   aligned lattice, not rot45 -- the rotated cell does not tile a 9 x 9 box and
                   the switching leaves no lattice strand either way.  A random f = 4 network
                   spans less volume per strand than the open diamond lattice, so at the lattice
                   box size it starts pre-stretched (~28 % on the isolated block); the aniso NPH
                   equilibration of slab_with_support.lmp sets lx, ly and the thickness.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from data_file_generator import www_bond_switch, verify_network_isotropy  # noqa: E402


def generate_gel_slab_periodic(beads_per_chain=5, units_x=9, units_y=9, units_z=21,
                               solvent_density=0.4,
                               below_thickness=4.0, top_thickness=20.0,
                               support_thickness=2.0,
                               support_spacing=0.2,
                               piston_exclusion=0.1,
                               piston_height=0.8,
                               output_file="gel_slab_periodic.data",
                               pre_swell=1.0,
                               bond_switch=0.0, bond_switch_seed=7, log=print):
    """
    Generate a tetrahedral (diamond-lattice) polymer gel slab that spans the whole
    x-y plane and is crosslinked PERIODICALLY across the x and y boundaries.

    Intended boundary condition for the run: `boundary p p p`.
    The polymer network bonds across x and y PBC only -- it does NOT bond across z.
    Each z-period is therefore one finite-thickness slab with its own solvent layer
    and its own support + piston, repeating in z:

        z-period (bottom -> top, then wraps):
          below_thickness solvent  |  support sheet  |  support_thickness gap
          |  polymer slab (units_z*a tall)  |  top_thickness solvent (piston parked high)

    Parameters
    ----------
    beads_per_chain : beads per chain (crosslink + interior beads + crosslink)
    units_x, units_y : in-plane unit cells; the box is exactly units*a wide (NO side
                       solvent padding) so the diamond lattice tiles seamlessly under PBC.
    units_z          : slab height in unit cells.
    solvent_density  : solvent number density throughout the box (beads/sigma^3).
    below_thickness  : solvent thickness below the support (z, sigma).
    top_thickness    : solvent thickness above the slab to the top of the period (z, sigma).
    support_thickness: support distance from the slab (sigma).
    support_spacing  : target spacing of wall atoms; auto-snapped so an integer number
                       of atoms fits box_x / box_y exactly (seamless tiling, no x/y seam).
    piston_exclusion : solvent exclusion half-width around the piston plane (sigma).
    piston_height    : piston position as a fraction (0..1) from slab top to top of box.
    output_file      : LAMMPS data file name.
    bond_switch      : > 0 (2026-11-02) re-wires the lattice strands into a random f = 4
                       network by Wooten-Winer-Weaire bond switching, this many accepted
                       switches per strand (data_file_generator.www_bond_switch, periodic in
                       x and y); 0 keeps the lattice (byte-identical to the notebook build).
    bond_switch_seed : RNG seed of the switching.
    """

    bead_spacing = 1.2
    chain_length = bead_spacing * (beads_per_chain - 1)
    z_clearance = 0
    # pre_swell (<1) starts the network pre-deswollen near its equilibrium size,
    # shrinking the slow collective-contraction distance. box_x=units*a so tiling
    # stays seamless under PBC. 1.0 = original behaviour.
    a = chain_length * pre_swell  # FCC/diamond lattice constant

    gel_x = units_x * a
    gel_y = units_y * a
    gel_z = units_z * a

    # No side padding: the slab spans the full x-y box.
    box_x = gel_x
    box_y = gel_y

    # z-layering identical to the non-periodic generator.
    offset = np.array([0.0, 0.0,
                       below_thickness + 2 * support_thickness + z_clearance])
    gel_top_z = offset[2] + gel_z
    box_z = gel_top_z + top_thickness
    piston_z = gel_top_z + piston_height * (box_z - gel_top_z)

    particles = []
    bonds = []
    particle_id = 1
    bond_id = 1
    molecule_id = 1

    # ---- Diamond-lattice crosslinks -------------------------------------------------
    # x,y are periodic: include cells i in [0, units_x) only (i = units_x is the
    # periodic image of i = 0). z is NOT periodic for the network: keep both end
    # layers so the slab has real top/bottom faces, exactly like a finite slab.
    crosslinks = {}
    fcc_positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.5, 0.5, 0.0]),
        np.array([0.5, 0.0, 0.5]),
        np.array([0.0, 0.5, 0.5]),
    ]

    for i in range(units_x):            # periodic -> no +1
        for j in range(units_y):        # periodic -> no +1
            for k in range(units_z + 1):  # finite slab -> keep top layer
                for sublattice in (0, 1):
                    for fcc_idx, fcc_frac in enumerate(fcc_positions):
                        frac_pos = fcc_frac + (np.array([0.25, 0.25, 0.25]) if sublattice else 0.0)
                        pos = (np.array([i, j, k]) + frac_pos) * a + offset
                        lx, ly, lz = pos[0] - offset[0], pos[1] - offset[1], pos[2] - offset[2]
                        # x,y: half-open [0, gel) so periodic images aren't duplicated.
                        # z:   closed [0, gel_z] so the top face layer is included.
                        if (-1e-9 <= lx < gel_x - 1e-9 and
                            -1e-9 <= ly < gel_y - 1e-9 and
                            -1e-9 <= lz <= gel_z + 1e-9):
                            key = (i, j, k, sublattice, fcc_idx)
                            crosslinks[key] = particle_id
                            particles.append({'id': particle_id, 'type': 1, 'pos': pos.copy(), 'mol': 0})
                            particle_id += 1

    # ---- Bond search: cell list with minimum-image in x,y only ----------------------
    bond_distance = a * np.sqrt(3) / 4
    tol = 0.1 * bond_distance
    r2_lo = (bond_distance - tol) ** 2
    r2_hi = (bond_distance + tol) ** 2
    r = bond_distance + tol

    cl_ids = list(crosslinks.values())
    cl_pos = np.array([particles[i - 1]['pos'] for i in cl_ids])
    log(f"  Crosslinks placed: {len(cl_ids)}; finding periodic bonds...")

    # periodic cell counts in x,y (>=3 so the 27-neighbour stencil is valid under wrap)
    ncx = max(3, int(np.floor(box_x / r)))
    ncy = max(3, int(np.floor(box_y / r)))
    csx, csy = box_x / ncx, box_y / ncy

    cells = defaultdict(list)
    cix = np.floor((cl_pos[:, 0] % box_x) / csx).astype(np.int64) % ncx
    ciy = np.floor((cl_pos[:, 1] % box_y) / csy).astype(np.int64) % ncy
    ciz = np.floor(cl_pos[:, 2] / r).astype(np.int64)   # z not periodic
    for idx in range(len(cl_ids)):
        cells[(cix[idx], ciy[idx], ciz[idx])].append(idx)

    offsets = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    neighbors = defaultdict(list)
    for idx in range(len(cl_ids)):
        id1 = cl_ids[idx]
        p1 = cl_pos[idx]
        cx, cy, cz = cix[idx], ciy[idx], ciz[idx]
        for dx, dy, dz in offsets:
            for jdx in cells.get(((cx + dx) % ncx, (cy + dy) % ncy, cz + dz), ()):
                id2 = cl_ids[jdx]
                if id2 <= id1:
                    continue
                d = cl_pos[jdx] - p1
                d[0] -= box_x * round(d[0] / box_x)   # min-image x
                d[1] -= box_y * round(d[1] / box_y)   # min-image y
                d2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]
                if r2_lo < d2 < r2_hi:
                    neighbors[id1].append(id2)
    for id1 in neighbors:
        neighbors[id1].sort()

    switch_report = None
    if bond_switch > 0:
        # Random f = 4 network (2026-11-02): re-wire the lattice strands by bond switching
        # (functionality and strand count kept), relax the junctions with minimum-image
        # distances across x/y and the z faces tethered, write the moved junction positions
        # back (x, y wrapped into the box) and rebuild the neighbour lists in the same
        # ascending-id form the chain builder reads.
        id_to_idx = {cid: k for k, cid in enumerate(cl_ids)}
        edges0 = np.array([(id_to_idx[a], id_to_idx[b]) for a in neighbors for b in neighbors[a]], dtype=np.int64)
        lo = offset.copy(); hi = offset + np.array([gel_x, gel_y, gel_z])
        cl_pos, edges_new, switch_report = www_bond_switch(
            cl_pos, edges0, bond_distance, lo, hi, switches_per_strand=bond_switch,
            seed=bond_switch_seed, keep_size=True, periodic=(box_x, box_y, 0.0), log=log)
        cl_pos[:, 0] %= box_x
        cl_pos[:, 1] %= box_y
        for k, cid in enumerate(cl_ids):
            particles[cid - 1]['pos'] = cl_pos[k].copy()
        neighbors = defaultdict(list)
        for i1, i2 in edges_new:
            a_, b_ = sorted((cl_ids[i1], cl_ids[i2]))
            neighbors[a_].append(b_)
        for id1 in neighbors:
            neighbors[id1].sort()

    # ---- Build chains; interior beads along the minimum-image vector, then wrapped --
    for id1 in cl_ids:
        pos1 = particles[id1 - 1]['pos']
        for id2 in neighbors.get(id1, ()):
            pos2 = particles[id2 - 1]['pos']
            delta = pos2 - pos1
            delta[0] -= box_x * round(delta[0] / box_x)
            delta[1] -= box_y * round(delta[1] / box_y)

            chain_ids = [id1]
            for b in range(1, beads_per_chain - 1):
                frac = b / (beads_per_chain - 1)
                pos = pos1 + frac * delta
                pos[0] %= box_x   # wrap into box for periodic dims
                pos[1] %= box_y
                particles.append({'id': particle_id, 'type': 2, 'pos': pos.copy(), 'mol': molecule_id})
                chain_ids.append(particle_id)
                particle_id += 1
            chain_ids.append(id2)

            particles[id1 - 1]['mol'] = molecule_id
            particles[id2 - 1]['mol'] = molecule_id
            for b in range(len(chain_ids) - 1):
                bonds.append({'id': bond_id, 'type': 1, 'atom1': chain_ids[b], 'atom2': chain_ids[b + 1]})
                bond_id += 1
            molecule_id += 1

    # ---- Support + piston sheets, hex packing, auto-snapped to tile seamlessly -------
    nx = max(1, round(box_x / support_spacing))
    snap_dx = box_x / nx
    row_pitch = support_spacing * np.sqrt(3)
    ny = max(2, round(box_y / row_pitch))
    if ny % 2:                      # even row count -> alternating offset wraps cleanly
        ny += 1
    snap_dy = box_y / ny

    def make_sheet(z, ptype):
        nonlocal particle_id, molecule_id
        for j in range(ny):
            y = j * snap_dy
            xshift = snap_dx / 2 if (j % 2) else 0.0
            for i in range(nx):
                x = (i * snap_dx + xshift) % box_x
                particles.append({'id': particle_id, 'type': ptype, 'pos': np.array([x, y, z]), 'mol': molecule_id})
                particle_id += 1
                molecule_id += 1

    make_sheet(below_thickness + support_thickness / 2 + z_clearance, 4)  # support (frozen)
    make_sheet(piston_z, 5)                                               # piston (mobile)

    # ---- Solvent: uniform fill across the whole box (matches original) --------------
    total_volume = box_x * box_y * box_z
    num_solvent_total = int(solvent_density * total_volume)
    log(f"  Placing up to {num_solvent_total} solvent beads...")
    np.random.seed(42)
    solvent_count, attempts, max_attempts = 0, 0, num_solvent_total * 100
    while solvent_count < num_solvent_total and attempts < max_attempts:
        attempts += 1
        pos = np.array([np.random.rand() * box_x, np.random.rand() * box_y, np.random.rand() * box_z])
        if abs(pos[2] - piston_z) < piston_exclusion:
            continue
        too_close = False
        for p in particles[-min(100, len(particles)):]:
            if np.linalg.norm(pos - p['pos']) < 0.8:
                too_close = True
                break
        if not too_close:
            particles.append({'id': particle_id, 'type': 3, 'pos': pos.copy(), 'mol': molecule_id})
            particle_id += 1
            molecule_id += 1
            solvent_count += 1
    if solvent_count < num_solvent_total:
        log(f"Warning: Only placed {solvent_count}/{num_solvent_total} solvent particles")

    # ---- Write LAMMPS data file -----------------------------------------------------
    kind = "random f=4 (bond-switched)" if bond_switch > 0 else "tetrahedral"
    with open(output_file, 'w') as f:
        f.write(f"LAMMPS data file: periodic (x,y) {kind} gel slab with support + piston\n\n")
        f.write(f"{len(particles)} atoms\n{len(bonds)} bonds\n0 angles\n0 dihedrals\n0 impropers\n\n")
        f.write("5 atom types\n1 bond types\n\n")
        f.write(f"0.0 {box_x:.6f} xlo xhi\n")
        f.write(f"0.0 {box_y:.6f} ylo yhi\n")
        f.write(f"0.0 {box_z:.6f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write("1 1.0  # Crosslink\n2 1.0  # Chain bead\n3 1.0  # Solvent\n4 1.0  # Support (frozen)\n5 1.0  # Piston (mobile)\n\n")
        f.write("Atoms\n\n")
        for p in particles:
            f.write(f"{p['id']} {p['mol']} {p['type']} {p['pos'][0]:.6f} {p['pos'][1]:.6f} {p['pos'][2]:.6f}\n")
        f.write("\nBonds\n\n")
        for b in bonds:
            f.write(f"{b['id']} {b['type']} {b['atom1']} {b['atom2']}\n")

    num_solvent_final = sum(1 for p in particles if p['type'] == 3)
    num_polymer = sum(1 for p in particles if p['type'] in (1, 2))
    num_support = sum(1 for p in particles if p['type'] == 4)
    num_piston = sum(1 for p in particles if p['type'] == 5)
    num_chains = len(bonds) // (beads_per_chain - 1) if beads_per_chain > 1 else 0
    log(f"Generated PERIODIC {kind} gel slab (boundary p p p):")
    log(f"  Unit cells: {units_x} x {units_y} x {units_z}   beads/chain: {beads_per_chain}")
    log(f"  Box (= gel in x,y): {box_x:.2f} x {box_y:.2f} x {box_z:.2f}")
    log(f"  Gel top z: {gel_top_z:.3f}   Piston z: {piston_z:.3f}   Bond length: {bond_distance:.4f}")
    log(f"  Wall tiling: {nx} x {ny} atoms/sheet, snapped spacing dx={snap_dx:.4f} dy(row)={snap_dy:.4f}")
    log(f"  Crosslinks: {len(crosslinks)}   Chains: {num_chains}   Polymer beads: {num_polymer}")
    log(f"  Support: {num_support}   Piston: {num_piston}   Solvent: {num_solvent_final}")
    if switch_report is not None:
        log(f"  Bond switching: {switch_report['accepted']} switches ({switch_report['per_strand']:.2f} per strand), "
            f"strands {switch_report['mean_len']:.3f} +/- {switch_report['std_len']:.3f} r0 (max {switch_report['max_len']:.3f}), "
            f"{100 * switch_report['lattice_strands_left']:.1f} % lattice strands left")
    log(f"  Total atoms: {len(particles)}   Total bonds: {len(bonds)}   Output: {output_file}")
    return output_file, dict(box=(box_x, box_y, box_z), bond_distance=bond_distance,
                             n_atoms=len(particles), n_bonds=len(bonds), switch_report=switch_report)


# ---------------------------------------------------------------------------
# THE KNOB: which network
# ---------------------------------------------------------------------------
NETWORKS = {
    'tetrahedral': dict(bond_switch=0.0),
    # the lattice slab randomised by bond switching (3 accepted switches per strand, seed 7,
    # the isolated 'random' gel's recipe)
    'random':      dict(bond_switch=3.0, bond_switch_seed=7),
}


def build(network='tetrahedral', output_file=None, beads_per_chain=5, units_x=9, units_y=9, units_z=23,
          solvent_density=0.4, below_thickness=4.0, top_thickness=25.0, support_thickness=2.0,
          support_spacing=0.2, piston_exclusion=0.1, piston_height=0.85, pre_swell=0.93,
          check_isotropy=None, log=print, **overrides):
    """Build one periodic slab data file.  `network` picks the preset ('tetrahedral' | 'random');
    everything else is generate_gel_slab_periodic's knob of the same name (the defaults are the
    reference 14000002 slab's: 9 x 9 x 23 cells, pre_swell 0.93, rho 0.40, 4 sigma below the
    support, 25 sigma above the gel, piston at 0.85).  `overrides` go straight to the generator
    (e.g. bond_switch=2.0, bond_switch_seed=11).  check_isotropy: run verify_network_isotropy
    on the written file with the x/y box periodic (default: only for 'random').  Returns the path."""
    if network not in NETWORKS:
        raise ValueError(f"network must be one of {list(NETWORKS)} (got {network!r})")
    if output_file is None:
        raise ValueError('output_file is required (e.g. ../../lammps_data_files_local/<name>.data)')
    kw = dict(NETWORKS[network])
    kw.update(overrides)
    log(f"=== periodic_data_file_generator.build: network = {network!r}  ->  {output_file}")
    _, info = generate_gel_slab_periodic(beads_per_chain, units_x, units_y, units_z,
                                         solvent_density, below_thickness, top_thickness,
                                         support_thickness, support_spacing,
                                         piston_exclusion, piston_height, output_file,
                                         pre_swell=pre_swell, log=log, **kw)
    if check_isotropy is None:
        check_isotropy = (network == 'random')
    if check_isotropy:
        log("")
        bx, by, _ = info['box']
        verify_network_isotropy(output_file, f"{network} periodic slab", periodic=(bx, by, 0.0), log=log)
    return output_file


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Build the x/y-periodic gel slab data file (see the module docstring).')
    ap.add_argument('--network', choices=sorted(NETWORKS), default='tetrahedral')
    ap.add_argument('--output', required=True, help='the .data file to write')
    ap.add_argument('--units', type=int, nargs=3, metavar=('X', 'Y', 'Z'), default=(9, 9, 23))
    ap.add_argument('--beads-per-chain', type=int, default=5)
    ap.add_argument('--rho-solv', type=float, default=0.4)
    ap.add_argument('--pre-swell', type=float, default=0.93)
    ap.add_argument('--below-thickness', type=float, default=4.0)
    ap.add_argument('--top-thickness', type=float, default=25.0)
    ap.add_argument('--support-thickness', type=float, default=2.0)
    ap.add_argument('--support-spacing', type=float, default=0.2)
    ap.add_argument('--piston-exclusion', type=float, default=0.1)
    ap.add_argument('--piston-height', type=float, default=0.85)
    ap.add_argument('--bond-switch', type=float, default=None, help='override the preset: accepted switches per strand')
    ap.add_argument('--bond-switch-seed', type=int, default=None)
    ap.add_argument('--no-check', action='store_true', help='skip the isotropy report')
    args = ap.parse_args()
    ov = {}
    if args.bond_switch is not None:
        ov['bond_switch'] = args.bond_switch
    if args.bond_switch_seed is not None:
        ov['bond_switch_seed'] = args.bond_switch_seed
    build(args.network, args.output, args.beads_per_chain, *args.units, solvent_density=args.rho_solv,
          below_thickness=args.below_thickness, top_thickness=args.top_thickness,
          support_thickness=args.support_thickness, support_spacing=args.support_spacing,
          piston_exclusion=args.piston_exclusion, piston_height=args.piston_height, pre_swell=args.pre_swell,
          check_isotropy=False if args.no_check else None, **ov)
