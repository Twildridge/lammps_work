"""data_file_generator.py -- the gel-in-a-bath data-file builder (2026-10-28; moved here from
scripts/slab_with_support.ipynb, now scripts/data_file_generator.ipynb, which only sets the knobs).

Three networks, chosen by `network` in build():
    'tetrahedral'  the diamond-lattice gel with the cube axes along the box axes (the default
                   build of every gel up to 2026-10-24: tall_rho04_nowalls, square_rho10_nowalls, ...);
    'rot45'        the same lattice turned 45 degrees about y (box x = crystal [101], z = [-101]),
                   so an xz shear of the box is a cube-diagonal shear of the lattice -- the
                   (C11 - C12)/2 constant the compression side measures (tall_rho04_nowalls_rot45);
    'random'       the rot45 block with its strands re-wired into a random f = 4 network by
                   Wooten-Winer-Weaire bond switching: no lattice directions, so the shear modulus
                   does not depend on orientation (tall_rho04_nowalls_random).
generate_gel_slab() is the engine (every knob, see its docstring); verify_network_isotropy() reads a
written file back and judges the strand-orientation statistics against random directions.
"""
import numpy as np


# ---------------------------------------------------------------------------
# PERIODIC BOXES (2026-11-02): the bond-switching code below works on an isolated block
# by default.  `periodic` = (Lx, Ly, Lz) with 0 for a non-periodic axis (the x/y-periodic
# slab of periodic_data_file_generator.py passes (box_x, box_y, 0)) makes every strand
# and repulsion distance a minimum-image distance and the neighbour trees toroidal along
# those axes; junctions on a periodic axis have no faces, so none are tethered there.
# periodic=None leaves the isolated path byte-identical to before.
# ---------------------------------------------------------------------------
def _min_image(d, periodic):
    """Minimum-image displacement(s) d (..., 3) for the periodic axes of `periodic`."""
    if periodic is None:
        return d
    d = np.array(d, dtype=float, copy=True)
    for ax in range(3):
        L = periodic[ax]
        if L and L > 0:
            d[..., ax] -= L * np.round(d[..., ax] / L)
    return d


def _wrap(P, periodic):
    """Copy of P with the periodic coordinates wrapped into [0, L) (what cKDTree(boxsize=) needs)."""
    if periodic is None:
        return P
    P = np.array(P, dtype=float, copy=True)
    for ax in range(3):
        L = periodic[ax]
        if L and L > 0:
            P[..., ax] = np.mod(P[..., ax], L)
    return P


def _tree(P, periodic):
    """cKDTree of the junctions, toroidal along the periodic axes (boxsize 0 = open axis)."""
    from scipy.spatial import cKDTree
    if periodic is None:
        return cKDTree(P)
    return cKDTree(_wrap(P, periodic), boxsize=[float(L) if (L and L > 0) else 0.0 for L in periodic])

# ---------------------------------------------------------------------------
# RANDOM f = 4 NETWORK BY BOND SWITCHING (2026-10-25)
#
# Wooten-Winer-Weaire bond switching on the crosslink graph of the diamond gel: pick a
# strand A-B, a neighbour C of A and a neighbour D of B, and replace the strands A-C and
# B-D by A-D and B-C.  Every junction keeps its functionality, the strand count and the
# strand length (beads) are untouched, and the network stays one connected piece (A-B
# is kept, so C and D stay attached through it).  Enough accepted switches erase the
# lattice: no cube axes left, so the shear modulus is the same for every orientation.
#
# Geometry: the junctions are relaxed between switch sweeps on a harmonic strand energy
# ((|r| - r0) / r0)^2, r0 = the lattice strand length a sqrt3 / 4, a soft repulsion
# between non-bonded junctions closer than rep_cut * r0 (an f = 4 spring network is
# under-constrained and would fold without it) and a tether that holds the junctions of
# the block faces at their lattice sites, so the block keeps its size.  Switches are
# proposed on the unrelaxed energy of the four strands involved and accepted with a
# Metropolis rule at a temperature annealed from T_hi to T_lo (in units of the strand
# energy), the Barkema-Mousseau form of the recipe: early sweeps accept almost every
# switch (randomisation), late sweeps only the ones that keep strands near r0.  A
# final stage re-wires every strand longer than max_len * r0 (1.9: 4.3 sigma end-to-end,
# the four beads then start < 1.1 sigma apart); the few the stretched network will not
# give up are allowed to 2.45 r0 (beads 1.4 sigma apart: FENE R0 = 1.5 sigma, and the
# deck's push-off runs on harmonic bonds first).
#
# verify_network_isotropy(data_file) reads a written data file back and reports what
# the shear test needs: the second-rank strand fabric (isotropic for ANY cubic lattice,
# so it cannot tell), the fourth-rank anisotropy (the Born term of the elastic tensor
# is nu k r0^2 <n n n n>), the Steinhardt Q4 / Q6 of the strand directions, the
# Born-affine Zener ratio A_xz = C44 / ((C11 - C12)/2) in the box frame (1 = isotropic;
# the aligned diamond gel has (C11 - C12)/2 = 0 at this level, the rot45 gel C44 = 0),
# and the share of strands that are still lattice strands (memory).  Pass/fail is
# judged against the same numbers for as many random directions.
# ---------------------------------------------------------------------------


def _strand_energy(x, edges, r0, rep_pairs, rep_cut, w_rep, teth_idx, teth_rel, teth_mask, centre, w_teth, all_rel, w_all, all_w, free_scale, periodic=None, grad=True, r0_strand=None):
    r0s = r0 if r0_strand is None else r0_strand      # rest length of the strand term (r0 normalises everything)
    if free_scale:
        lam = x[-1]; P = x[:-1].reshape(-1, 3)
    else:
        lam = 1.0; P = x.reshape(-1, 3)
    d = _min_image(P[edges[:, 1]] - P[edges[:, 0]], periodic)
    L = np.linalg.norm(d, axis=1)
    s = (L - r0s) / r0
    E = np.sum(s * s)
    G = np.zeros_like(P); dlam = 0.0
    if grad:
        f = (2.0 * s / (r0 * L))[:, None] * d          # dE/dP_j
        np.add.at(G, edges[:, 1], f)
        np.add.at(G, edges[:, 0], -f)
    if len(rep_pairs):
        dr = _min_image(P[rep_pairs[:, 1]] - P[rep_pairs[:, 0]], periodic)
        Lr = np.linalg.norm(dr, axis=1)
        m = Lr < rep_cut
        if m.any():
            sr = (rep_cut - Lr[m]) / r0
            E += w_rep * np.sum(sr * sr)
            if grad:
                fr = (-2.0 * w_rep * sr / (r0 * Lr[m]))[:, None] * dr[m]
                np.add.at(G, rep_pairs[m, 1], fr)
                np.add.at(G, rep_pairs[m, 0], -fr)
    if len(teth_idx):
        dt = (P[teth_idx] - centre - lam * teth_rel) * teth_mask / r0      # face-normal components only
        E += w_teth * np.sum(dt * dt)
        if grad:
            G[teth_idx] += 2.0 * w_teth * dt / r0
            dlam = -2.0 * w_teth * np.sum(dt * teth_rel) / r0
    if w_all > 0:
        da = (P - centre - lam * all_rel) / r0
        E += w_all * np.sum(all_w[:, None] * da * da)
        if grad:
            G += 2.0 * w_all * all_w[:, None] * da / r0
            dlam += -2.0 * w_all * np.sum(all_w[:, None] * da * all_rel) / r0
    if not grad:
        return E
    return (E, np.append(G.ravel(), dlam)) if free_scale else (E, G.ravel())


def _relax_junctions(P, edges, r0, teth_idx, teth_rel, teth_mask, centre, rep_cut, w_rep, w_teth, all_rel, w_all, all_w,
                     lam=1.0, free_scale=False, rounds=8, maxiter=80, periodic=None, r0_strand=None):
    """L-BFGS on the junction positions (and the block scale when free_scale), in rounds so
    that the repulsion pair list follows the positions."""
    from scipy.optimize import minimize
    from scipy.spatial import cKDTree
    bonded = set(map(tuple, np.sort(edges, axis=1)))
    Ef = 0.0
    for _ in range(rounds):
        pairs = _tree(P, periodic).query_pairs(1.6 * rep_cut, output_type='ndarray')
        pairs = np.array([p for p in map(tuple, np.sort(pairs, axis=1)) if p not in bonded], dtype=np.int64).reshape(-1, 2)
        x0 = np.append(P.ravel(), lam) if free_scale else P.ravel()
        res = minimize(_strand_energy, x0, jac=True, method='L-BFGS-B',
                       args=(edges, r0, pairs, rep_cut, w_rep, teth_idx, teth_rel, teth_mask, centre, w_teth, all_rel, w_all, all_w, free_scale, periodic, True, r0_strand),
                       options={'maxiter': maxiter})
        Pn = (res.x[:-1] if free_scale else res.x).reshape(-1, 3)
        if free_scale:
            lam = float(res.x[-1])
        moved = np.abs(Pn - P).max()
        P, Ef = Pn, res.fun
        if moved < 0.02 * r0 and res.nit < maxiter:
            break
    return P, Ef, lam


def www_bond_switch(P0, edges0, r0, block_lo, block_hi, switches_per_strand=3.0, seed=7,
                    T_hi=4.0, T_lo=0.01, n_sweeps=16, n_random=6, rep_cut=0.9, w_rep=1.0, w_teth=4.0, w_all=0.1,
                    max_len=1.9, max_len_hard=2.45, face_tol=0.3, keep_size=True, reset_anchors=False,
                    face_full=True, anchor_jitter=0.25, forbid_lattice=True, periodic=None, rest_len=1.0, log=print):
    """Return (P, edges, report).  P0 (n, 3) junction positions, edges0 (m, 2) int junction
    indices (any order), r0 the lattice strand length.  block_lo / block_hi: the block
    faces (junctions within face_tol * r0 of a face are tethered).

    Defaults set on the 13 x 9 x 30 rot45 block (2026-10-25): w_all 0.1 with the lattice
    sites as anchors (reset_anchors False) gave the most isotropic interior (strands one
    strand length or more inside the faces: 4th-rank anisotropy, Q4 and the Born Zener
    ratio all within 1.5 sd of random directions); re-anchoring to the randomised
    positions (reset_anchors True) or a weaker w_all let a few-percent uniaxial texture
    along the long axis develop.  keep_size False (the block scale relaxed as well) is
    experimental: with only the junction repulsion to hold it the block shrinks to about
    half, so leave it True.

    face_full True (2026-10-25, second version) holds the face junctions in all three
    directions, not only along the face normal: held along the normal alone they slid
    sideways and bunched, which left voids and loose loops on the block surface.
    anchor_jitter (in r0) displaces every tether anchor (face and weak) by a fixed random
    offset, so the lattice spacing of neighbouring junctions is not re-formed by the
    anneal while the density stays uniform; 0 anchors on the lattice sites.  Default
    0.25 r0 = 0.6 sigma.  w_teth 4 (was 1): the face junctions that the switching
    connects to partners deep inside were pulled up to 3 sigma into the block, the dents
    and loose loops seen on the surface of the first file; at 4 the outermost junction of
    a 4 x 4 sigma face patch sits 0.25 sigma inside the face on average (1.3 before,
    0.05 for the lattice) and the interior isotropy is unchanged (surface_recess reports it).
    forbid_lattice True: no switch may re-create a strand of the original lattice, and the
    lattice strands that survive the randomising sweeps are switched away afterwards (a
    purge pass).  Without it 7 % of the strands were still lattice pairs, concentrated at
    the tethered faces; the big x-faces hold strands with no x component, which showed as
    a 2 % deficit of x-strands through the whole block (second-rank texture).
    periodic (2026-11-02): (Lx, Ly, Lz) with 0 for an open axis, see _min_image above; the
    x/y-periodic slab passes (box_x, box_y, 0) so that strands and repulsions across the
    x/y boundaries are minimum-image and only the two z faces are tethered.
    rest_len (2026-11-03): the rest length of the strand energy and of the switch acceptance
    in units of r0.  At 1.0 (the default, every random file up to 2026-11-02) the relaxation
    pulled the junctions into clumps -- a random f = 4 network at the lattice's junction
    density wants strands ~1.5 r0 long, and shortening them towards r0 means piling junctions
    together -- which left 3-5 sigma voids inside the gel (nearest other junction 0.84 r0 on
    average, up to 11 within 0.9 r0; polymer density on 5 sigma voxels CoV 0.38 against 0.035
    for the lattice).  The random presets pass RANDOM_RELAX (rest_len 1.4, rep_cut 1.3, w_rep
    3): the junctions are spaced like a soft-sphere liquid (CoV 0.09, largest empty sphere
    2.1 sigma, the junction count within 2 r0 as uniform as the lattice's), strands 1.43 +/-
    0.21 r0, isotropy unchanged, the face junctions recessed a uniform 1.2 sigma.  The
    defaults here stay the old ones so that switch_options={} reproduces the old files."""
    rng = np.random.default_rng(seed)
    P = P0.copy()
    n = len(P)
    rep_cut = rep_cut * r0
    r0s = rest_len * r0                       # rest length of the strand energy and of the switch acceptance
    face_mask = np.zeros((n, 3))
    for ax in range(3):
        if periodic is not None and periodic[ax] and periodic[ax] > 0:
            continue                              # a periodic axis has no faces
        face_mask[(P[:, ax] - block_lo[ax] < face_tol * r0) | (block_hi[ax] - P[:, ax] < face_tol * r0), ax] = 1.0
    teth_idx = np.flatnonzero(face_mask.any(axis=1))
    teth_mask = face_mask[teth_idx]
    centre = 0.5 * (np.asarray(block_lo) + np.asarray(block_hi))
    jit = anchor_jitter * r0 * rng.normal(size=P0.shape) if anchor_jitter > 0 else np.zeros_like(P0)
    jit[teth_idx] *= (1.0 - face_mask[teth_idx])          # never jitter a face junction off its face plane
    all_rel = P0 + jit - centre               # weak tether of every junction to its (jittered) lattice site (w_all)
    teth_rel = all_rel[teth_idx]              # tether target = centre + lam * teth_rel
    if face_full:
        teth_mask = np.ones_like(teth_mask)   # face junctions held in all three directions
    all_w = np.ones(n)                        # per-junction weight on that tether (released for stuck strands)
    free_scale = not keep_size                # keep_size False: the block scale lam is relaxed too
    lam = 1.0

    adj = [set() for _ in range(n)]
    for a, b in edges0:
        adj[a].add(b); adj[b].add(a)
    deg0 = np.array([len(s) for s in adj])
    lattice_edges = set(map(tuple, np.sort(np.asarray(edges0), axis=1)))

    def lkey(a, b):
        return (a, b) if a < b else (b, a)

    def elen(a, b):
        return np.linalg.norm(_min_image(P[a] - P[b], periodic))

    def slen(E):
        return np.linalg.norm(_min_image(P[E[:, 1]] - P[E[:, 0]], periodic), axis=1)

    def sweep(T, n_trials, greedy=False):
        acc = 0
        A_list = rng.integers(0, n, n_trials)
        for A in A_list:
            if len(adj[A]) < 2:
                continue
            B = rng.choice(list(adj[A]))
            if len(adj[B]) < 2:
                continue
            C = rng.choice([c for c in adj[A] if c != B])
            Ds = [d for d in adj[B] if d != A and d != C and d not in adj[A]]
            if not Ds or C in adj[B]:
                continue
            D = Ds[rng.integers(len(Ds))]
            if forbid_lattice and (lkey(A, D) in lattice_edges or lkey(B, C) in lattice_edges):
                continue
            dE = ((elen(A, D) - r0s) ** 2 + (elen(B, C) - r0s) ** 2
                  - (elen(A, C) - r0s) ** 2 - (elen(B, D) - r0s) ** 2) / r0 ** 2
            if dE < 0 or (not greedy and T > 0 and rng.random() < np.exp(-dE / T)):
                adj[A].remove(C); adj[C].remove(A); adj[B].remove(D); adj[D].remove(B)
                adj[A].add(D); adj[D].add(A); adj[B].add(C); adj[C].add(B)
                acc += 1
        return acc

    def edge_array():
        return np.array([(a, b) for a in range(n) for b in adj[a] if b > a], dtype=np.int64)

    m = len(edges0)
    target = int(switches_per_strand * m)
    accepted = 0
    Ts = np.concatenate([np.full(n_random, T_hi), np.geomspace(T_hi, T_lo, n_sweeps)])
    log(f"  bond switching: {n} junctions, {m} strands, r0 = {r0:.3f}, target {target} accepted switches "
        f"({switches_per_strand} per strand), {n_random} sweeps at T {T_hi} then {n_sweeps} sweeps T {T_hi} -> {T_lo}, "
        f"{len(teth_idx)} face junctions tethered")
    for k, T in enumerate(Ts):
        acc = sweep(T, m)
        accepted += acc
        E = edge_array()
        P, Ef, lam = _relax_junctions(P, E, r0, teth_idx, teth_rel, teth_mask, centre, rep_cut, w_rep, w_teth, all_rel, w_all, all_w, lam, free_scale, periodic=periodic, r0_strand=r0s)
        L = slen(E) / r0
        log(f"    sweep {k + 1:2d}  T {T:.3f}  accepted {acc:6d} (total {accepted:7d})  E/strand {Ef / m:.4f}  "
            f"strand length / r0: mean {L.mean():.3f}  max {L.max():.3f}  block scale {lam:.3f}")
        if reset_anchors and k + 1 == n_random:
            # The lattice sites must not survive as tether targets: from here on the weak
            # tether of every junction, and the in-face part of the face tether, anchor to
            # the randomised positions (the face-normal targets keep the block faces).
            all_rel = (P - centre) / lam
            teth_rel = np.where(teth_mask > 0, teth_rel, (P[teth_idx] - centre) / lam)
            log("    tether anchors reset to the randomised positions (lattice sites forgotten)")
    # top up at the final temperature if the target has not been reached
    extra = 0
    while accepted < target and extra < 10:
        acc = sweep(T_lo, m); accepted += acc; extra += 1
        E = edge_array(); P, Ef, lam = _relax_junctions(P, E, r0, teth_idx, teth_rel, teth_mask, centre, rep_cut, w_rep, w_teth, all_rel, w_all, all_w, lam, free_scale, periodic=periodic, r0_strand=r0s)
        log(f"    top-up {extra}  accepted {acc} (total {accepted})")
    # clean-up: every strand longer than max_len * r0 is attacked directly.  To remove a
    # strand A-C the pivot is another strand A-B of A (or C-B of C): A-C and B-D go, A-D
    # and B-C come, and the best of all (B, D) choices by the local energy is taken when
    # it lowers it.  Relax, repeat until no strand exceeds the limit.
    from scipy.spatial import cKDTree

    def shorten(limit, targets=None):
        E = edge_array()
        L = slen(E) / r0
        n_fixed = 0
        for a, c in (E[np.argsort(-L)] if targets is None else targets):
            if c not in adj[a]:          # already re-wired by an earlier move of this pass
                continue
            if targets is None and elen(a, c) / r0 <= limit:
                break
            best = None
            Lac = elen(a, c) if targets is None else limit * r0
            for A, C in ((a, c), (c, a)):
                for B in adj[A]:
                    if B == C:
                        continue
                    for D in adj[B]:
                        if D == A or D == C or D in adj[A] or C in adj[B]:
                            continue
                        if forbid_lattice and (lkey(A, D) in lattice_edges or lkey(B, C) in lattice_edges):
                            continue
                        worst = max(elen(A, D), elen(B, C))      # the pivot A-B stays either way
                        if worst < Lac and (best is None or worst < best[0]):
                            best = (worst, A, B, C, D)
            if best is not None:
                _, A, B, C, D = best
                adj[A].remove(C); adj[C].remove(A); adj[B].remove(D); adj[D].remove(B)
                adj[A].add(D); adj[D].add(A); adj[B].add(C); adj[C].add(B)
                n_fixed += 1
                continue
            # no switch around this strand helps (its ends are low-functionality face or corner
            # junctions as a rule): exchange partners with a strand E-F that has E near a and
            # F near c, (a-c, E-F) -> (a-F, E-c), degrees unchanged; the one move that can split
            # the graph, so it is undone if the network is no longer one piece.
            tree = _tree(P, periodic)
            best = None
            for a_, c_ in ((a, c), (c, a)):
                for E in tree.query_ball_point(_wrap(P[a_], periodic), 2.5 * r0):
                    if E == a_ or E == c_:
                        continue
                    for F in adj[E]:
                        if F in (a_, c_) or F in adj[a_] or c_ in adj[E] or E == F:
                            continue
                        if forbid_lattice and (lkey(a_, F) in lattice_edges or lkey(E, c_) in lattice_edges):
                            continue
                        worst = max(elen(a_, F), elen(E, c_))
                        if worst < Lac and (best is None or worst < best[0]):
                            best = (worst, E, F, a_, c_)
            if best is not None:
                _, E, F, a, c = best
                adj[a].remove(c); adj[c].remove(a); adj[E].remove(F); adj[F].remove(E)
                adj[a].add(F); adj[F].add(a); adj[E].add(c); adj[c].add(E)
                if not connected():
                    adj[a].remove(F); adj[F].remove(a); adj[E].remove(c); adj[c].remove(E)
                    adj[a].add(c); adj[c].add(a); adj[E].add(F); adj[F].add(E)
                else:
                    n_fixed += 1
        return n_fixed

    def connected():
        seen = np.zeros(n, bool); stack = [0]; seen[0] = True
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True; stack.append(v)
        return bool(seen.all())

    if forbid_lattice:
        for k in range(3):
            E = edge_array()
            left = [tuple(e) for e in E if tuple(e) in lattice_edges]
            if not left:
                break
            nf = shorten(max_len, targets=left)
            E = edge_array(); P, Ef, lam = _relax_junctions(P, E, r0, teth_idx, teth_rel, teth_mask, centre, rep_cut, w_rep, w_teth, all_rel, w_all, all_w, lam, free_scale, periodic=periodic, r0_strand=r0s)
            log(f"    lattice purge {k + 1}: {len(left)} lattice strands left, {nf} switched away")
    for k in range(30):
        E = edge_array()
        L = slen(E) / r0
        if L.max() <= max_len:
            break
        nf = shorten(max_len)
        if nf == 0:
            E = edge_array(); L = slen(E) / r0
            for a, c in E[L > max_len]:
                log(f"      stuck strand {a}-{c}: {L[(E[:, 0] == a) & (E[:, 1] == c)][0]:.2f} r0, degrees {len(adj[a])}/{len(adj[c])}, "
                    f"on face {bool(face_mask[a].any())}/{bool(face_mask[c].any())}")
            break
        E = edge_array(); P, Ef, lam = _relax_junctions(P, E, r0, teth_idx, teth_rel, teth_mask, centre, rep_cut, w_rep, w_teth, all_rel, w_all, all_w, lam, free_scale, periodic=periodic, r0_strand=r0s)
        L = slen(E) / r0
        log(f"    clean-up {k + 1}: {nf} long strands re-wired, max strand {L.max():.3f} r0, mean {L.mean():.3f}, "
            f"{int((L > max_len).sum())} still over {max_len}")
    # Strands still over the limit end on tethered low-functionality face junctions whose
    # partner sits deep inside: release those junctions from their tethers (a dent in the
    # face, nothing else) and relax again.
    for k in range(3):
        E = edge_array()
        L = slen(E) / r0
        if L.max() <= max_len:
            break
        rel = np.unique(E[L > max_len].ravel())
        all_w[rel] = 0.0
        pos_in_teth = {j: i for i, j in enumerate(teth_idx)}
        for j in rel:
            if j in pos_in_teth:
                teth_mask[pos_in_teth[j]] = 0.0     # a freed face junction may dent the face
        P, Ef, lam = _relax_junctions(P, E, r0, teth_idx, teth_rel, teth_mask, centre, rep_cut, w_rep, w_teth, all_rel, w_all, all_w, lam, free_scale, periodic=periodic, r0_strand=r0s)
        L = slen(E) / r0
        log(f"    release {k + 1}: {len(rel)} junctions freed from their tethers, max strand {L.max():.3f} r0, {int((L > max_len).sum())} still over {max_len}")
    E = edge_array()
    L = slen(E) / r0
    deg = np.array([len(s) for s in adj])
    assert np.array_equal(deg, deg0), "bond switching changed a functionality"
    assert len(E) == m, "bond switching changed the strand count"
    # connectivity
    seen = np.zeros(n, bool); stack = [0]; seen[0] = True
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if not seen[v]:
                seen[v] = True; stack.append(v)
    n_comp_rest = int((~seen).sum())
    kept = sum(1 for e in map(tuple, E) if e in lattice_edges)
    report = dict(accepted=accepted, per_strand=accepted / m, max_len=float(L.max()), mean_len=float(L.mean()),
                  std_len=float(L.std()), lattice_strands_left=kept / m, unreached_junctions=n_comp_rest,
                  face_shift=float(np.abs((P[teth_idx] - centre - lam * teth_rel) * teth_mask).max()) if len(teth_idx) else 0.0, block_scale=lam)
    log(f"  bond switching done: {accepted} switches ({accepted / m:.2f} per strand), strands {L.mean():.3f} +/- {L.std():.3f} r0 "
        f"(max {L.max():.3f}), {100 * kept / m:.1f} % lattice strands left, {n_comp_rest} junctions unreached, "
        f"largest face-junction shift off its tether {report['face_shift']:.2f} sigma, block scale {lam:.3f}")
    n_over = int((L > max_len).sum())
    report['strands_over_max_len'] = n_over
    if L.max() > max_len_hard:
        raise RuntimeError(f"bond switching: a strand is {L.max():.2f} r0 long after clean-up (hard limit {max_len_hard})")
    if n_over:
        log(f"  WARNING: {n_over} strand(s) of {m} stay between {max_len} and {max_len_hard} r0 (longest {L.max():.2f} r0 = "
            f"{L.max() * r0:.2f} sigma end-to-end, beads {L.max() * r0 / 4:.2f} sigma apart; FENE R0 = 1.5, the push-off is harmonic)")
    if n_comp_rest:
        raise RuntimeError("bond switching: the network is not one connected piece")
    return P, E, report


# ---------------------------------------------------------------------------
# Isotropy report
# ---------------------------------------------------------------------------
def _strand_vectors_from_data(data_file, periodic=None):
    L = open(data_file).read().split('\n')
    ia = L.index('Atoms'); ib = L.index('Bonds')
    A = np.array([l.split()[:6] for l in L[ia + 2:ib] if l.strip()], dtype=float)
    ids = A[:, 0].astype(int); typ = A[:, 2].astype(int); pos = A[:, 3:6]
    B = np.array([l.split()[2:4] for l in L[ib + 2:] if l.strip()], dtype=int)
    idx = {i: k for k, i in enumerate(ids)}
    from collections import defaultdict, Counter
    adj = defaultdict(list)
    for a, b in B:
        adj[a].append(b); adj[b].append(a)
    cl = ids[typ == 1]
    vec = []; pairs = []
    for c in cl:
        for nb in adj[c]:
            p, q = c, nb
            while typ[idx[q]] != 1:
                r = [s for s in adj[q] if s != p][0]; p, q = q, r
            if q > c:
                vec.append(_min_image(pos[idx[q]] - pos[idx[c]], periodic)); pairs.append((c, q))
    deg = Counter(len(adj[c]) for c in cl)
    return np.array(vec), np.array(pairs), dict(sorted(deg.items())), pos[(typ == 1) | (typ == 2)]


def isotropy_metrics(V):
    """V (m, 3) strand vectors.  Returns a dict of orientation metrics."""
    try:
        from scipy.special import sph_harm            # scipy < 1.17: sph_harm(m, l, azimuth, polar)
    except ImportError:                                # scipy >= 1.17 renamed it (2026-11-03)
        from scipy.special import sph_harm_y
        def sph_harm(mm, l, theta, phi):
            return sph_harm_y(l, mm, phi, theta)       # sph_harm_y(l, m, polar, azimuth)
    n = V / np.linalg.norm(V, axis=1)[:, None]
    m = len(n)
    F2 = n.T @ n / m
    T4 = np.einsum('ai,aj,ak,al->ijkl', n, n, n, n) / m
    d = np.eye(3)
    I4 = (np.einsum('ij,kl->ijkl', d, d) + np.einsum('ik,jl->ijkl', d, d) + np.einsum('il,jk->ijkl', d, d)) / 15.0
    aniso4 = np.linalg.norm(T4 - I4) / np.linalg.norm(I4)
    # Steinhardt Q_l of the strand directions (axial: use +n and -n)
    nn = np.vstack([n, -n])
    theta = np.arctan2(nn[:, 1], nn[:, 0]); phi = np.arccos(np.clip(nn[:, 2], -1, 1))
    Q = {}
    for l in (4, 6):
        Ylm = np.array([np.mean(sph_harm(mm, l, theta, phi)) for mm in range(-l, l + 1)])
        Q[l] = float(np.sqrt(4 * np.pi / (2 * l + 1) * np.sum(np.abs(Ylm) ** 2)))
    nx4 = T4[0, 0, 0, 0]; nx2ny2 = T4[0, 0, 1, 1]; nx2nz2 = T4[0, 0, 2, 2]
    den = nx4 - nx2ny2
    A_xz = 2 * nx2nz2 / den if abs(den) > 1e-9 else np.inf
    return dict(n=m, fabric_eig=np.linalg.eigvalsh(F2), aniso4=float(aniso4), Q4=Q[4], Q6=Q[6],
                nx4=float(nx4), nx2ny2=float(nx2ny2), nx2nz2=float(nx2nz2), A_xz=float(A_xz))


def surface_recess(cl, lo, hi, patch=4.0, axes=(0, 1, 2)):
    """Height map of the six faces of a junction cloud: for every patch x patch window of
    a face, how far inside the nominal face (lo / hi of the cloud) the outermost junction
    sits.  Returns mean / sd / max recess (sigma) and the fractions of patches recessed by
    more than 2 and 4 sigma.  A lattice block gives ~0."""
    rec = []
    for ax in axes:
        u, v = [a for a in range(3) if a != ax]
        nu = max(int(round((hi[u] - lo[u]) / patch)), 1); nv = max(int(round((hi[v] - lo[v]) / patch)), 1)
        iu = np.minimum(((cl[:, u] - lo[u]) / (hi[u] - lo[u]) * nu).astype(int), nu - 1)
        iv = np.minimum(((cl[:, v] - lo[v]) / (hi[v] - lo[v]) * nv).astype(int), nv - 1)
        key = iu * nv + iv
        for side, sign, ref in ((0, 1.0, lo[ax]), (1, -1.0, hi[ax])):
            depth = sign * (cl[:, ax] - ref)                   # distance inside the face
            outer = np.full(nu * nv, np.inf)
            np.minimum.at(outer, key, depth)
            rec.append(outer[np.isfinite(outer)])
    rec = np.concatenate(rec)
    return dict(mean=float(rec.mean()), sd=float(rec.std()), max=float(rec.max()),
                frac2=float((rec > 2).mean()), frac4=float((rec > 4).mean()))


def verify_density_uniformity(data_file, label=None, voxel=5.0, margin=6.0, h=0.5, probe=3.0,
                              cov_max=0.15, log=print, periodic=None, scale=1.0):
    """Interior density report of a gel data file (2026-11-03): the polymer beads (types 1, 2)
    are binned on `voxel`-sigma cubes more than `margin` sigma inside the gel (the gel
    extent along each open axis comes from the polymer density profile; a `periodic` axis
    spans the box) and the distance from a grid of points (spacing h) to the nearest polymer
    bead is measured.  Reports the coefficient of variation of the voxel density, the
    emptiest voxel, the largest empty sphere and the volume fraction farther than `probe`
    sigma from any polymer.  Verdict DENSITY OK when CoV <= cov_max and no empty sphere
    reaches `probe` sigma.  The as-built lattice gives CoV 0.035 and 1.8 sigma; the random
    network relaxed with the old knobs 0.35-0.41 and 4-5 sigma (empty 5-sigma voxels), with
    RANDOM_RELAX 0.08-0.10 and 1.9-2.3 sigma; the equilibrated lattice slab 0.13 and 2.2.
    `scale` multiplies voxel, margin and probe: the numbers are calibrated on the 5-bead
    mesh (a = 4.46 sigma pre-swollen); a coarser mesh is judged on the same number of meshes
    per voxel with scale = a / 4.46 (the builders pass it).  Returns the numbers as a dict."""
    from scipy.spatial import cKDTree
    voxel, margin, probe = voxel * scale, margin * scale, probe * scale
    L = open(data_file).read().split('\n')
    box = np.array([[float(v) for v in l.split()[:2]] for l in L[:40] if 'xlo' in l or 'ylo' in l or 'zlo' in l])
    ia = L.index('Atoms'); ib = L.index('Bonds')
    A = np.array([l.split()[:6] for l in L[ia + 2:ib] if l.strip()], dtype=float)
    typ = A[:, 2].astype(int)
    Lb = box[:, 1] - box[:, 0]
    pol = (A[typ <= 2][:, 3:6] - box[:, 0]) % Lb
    per = [bool(periodic is not None and periodic[ax] and periodic[ax] > 0) for ax in range(3)]
    lo = np.zeros(3); hi = Lb.copy()
    for ax in range(3):
        if not per[ax]:
            hist, edges = np.histogram(pol[:, ax], bins=max(10, int(Lb[ax])), range=(0, Lb[ax]))
            on = np.flatnonzero(hist > 0.1 * hist.max())
            lo[ax], hi[ax] = edges[on[0]], edges[on[-1] + 1]
    ilo = np.where(per, 0.0, lo + margin); ihi = np.where(per, Lb, hi - margin)
    tree = cKDTree(pol, boxsize=[Lb[ax] if per[ax] else 0.0 for ax in range(3)])
    axes = [np.arange(h / 2, Lb[ax], h) if per[ax] else np.arange(ilo[ax], ihi[ax], h) for ax in range(3)]
    G = np.array(np.meshgrid(*axes, indexing='ij')).reshape(3, -1).T
    d, _ = tree.query(G, workers=-1)
    inside = np.all((pol > ilo) & (pol < ihi), axis=1)
    nv = [max(1, int(round((ihi[ax] - ilo[ax]) / voxel))) for ax in range(3)]
    H, _ = np.histogramdd(pol[inside], bins=nv, range=[(ilo[ax], ihi[ax]) for ax in range(3)])
    rho = H / np.prod((ihi - ilo) / nv)
    cov = float(rho.std() / rho.mean())
    big = np.flatnonzero(d > probe)
    n_void = 0
    if len(big):
        taken = np.zeros(len(big), bool)
        gt = cKDTree(G[big], boxsize=[Lb[ax] if per[ax] else 0.0 for ax in range(3)])
        for k in np.argsort(-d[big]):
            if taken[k]:
                continue
            taken[gt.query_ball_point(G[big[k]], 4.0)] = True
            n_void += 1
    ok = cov <= cov_max and d.max() < probe
    log(f"=== density report: {label or data_file} ===")
    log(f"  interior = " + ", ".join(f"{'xyz'[ax]} {'box (periodic)' if per[ax] else f'{ilo[ax]:.1f}..{ihi[ax]:.1f}'}" for ax in range(3)) +
        f"; polymer density {inside.sum() / np.prod(ihi - ilo):.3f}")
    log(f"  polymer density on {voxel:g} sigma voxels: CoV {cov:.3f} (lattice 0.035), emptiest / mean {rho.min() / rho.mean():.2f}, fullest / mean {rho.max() / rho.mean():.2f}")
    log(f"  nearest polymer bead from the interior grid: mean {d.mean():.2f}, p99 {np.percentile(d, 99):.2f}, largest empty sphere {d.max():.2f} sigma (lattice 1.8); "
        f"{100 * (d > probe).mean():.3f} % of the volume > {probe:g} sigma from any bead in {n_void} void(s)")
    log(f"  DENSITY {'OK' if ok else 'NOT uniform'} (CoV <= {cov_max} and no empty sphere of {probe:g} sigma)")
    return dict(cov=cov, rho_min_rel=float(rho.min() / rho.mean()), rho_max_rel=float(rho.max() / rho.mean()),
                d_mean=float(d.mean()), d_max=float(d.max()), frac_over_probe=float((d > probe).mean()), n_void=n_void, ok=ok)


def verify_network_isotropy(data_file, label=None, seed=0, n_ref=20, face_depth=None, log=print, periodic=None):
    """Read a gel data file back and print the orientation report of its strands next to
    the same numbers for as many random directions (mean +/- sd over n_ref draws).  The
    strands are also split into INTERIOR (midpoint more than face_depth inside the
    crosslink bounding box; default 7 sigma, three lattice strand lengths) and SKIN: a
    finite block cannot have strands poking out of it, and the tethered, pre-stretched
    faces of the switched network imprint an in-plane bias that reaches two to three
    strand lengths in (measured 2026-10-25: <n_x^2> 0.27-0.31 within 4 sigma of the
    x-faces, <n_y^2> 0.27-0.30 at 2-8 sigma from the y-faces, isotropic beyond 8 sigma).
    The equilibration run relaxes that skin with the pre-stretch; the shear deck's grips
    and core margins sit on it.  The verdict is the interior's: each metric within 3 sd
    of random directions.
    periodic (2026-11-02): (Lx, Ly, 0) for the x/y-periodic slab -- strand vectors are then
    minimum-image and only the z faces count as faces (depth, surface recess).
    Returns (metrics_all, metrics_interior, passed)."""
    V, pairs, deg, pol = _strand_vectors_from_data(data_file, periodic)
    open_axes = tuple(ax for ax in range(3) if periodic is None or not (periodic[ax] and periodic[ax] > 0))
    L = open(data_file).read().split('\n')
    ia = L.index('Atoms'); ib = L.index('Bonds')
    A = np.array([l.split()[:6] for l in L[ia + 2:ib] if l.strip()], dtype=float)
    idx = {int(i): k for k, i in enumerate(A[:, 0])}
    cl = A[A[:, 2] == 1, 3:6]
    lo, hi = cl.min(0), cl.max(0)
    mid = np.array([A[idx[a], 3:6] for a, b in pairs]) + 0.5 * V          # minimum-image midpoint
    depth = np.minimum(mid - lo, hi - mid)[:, list(open_axes)].min(axis=1)
    if face_depth is None:
        face_depth = 7.0       # three lattice strand lengths (2.29 sigma each): the skin the tethered faces imprint
    rng = np.random.default_rng(seed)
    Lr = np.linalg.norm(V, axis=1)
    # The face planes for the recess map are the 0.2 / 99.8 percentiles of the junction cloud
    # (2026-11-03): a handful of junctions pushed out through a face (the random presets'
    # junction repulsion does that to a few of 15 000) used to set lo / hi and read as a
    # 1-2 sigma recess of every patch.  They are counted separately as pop-outs.
    lo_r, hi_r = np.percentile(cl, 0.2, axis=0), np.percentile(cl, 99.8, axis=0)
    rec = surface_recess(cl, lo_r, hi_r, axes=open_axes)
    n_pop = int(sum(((cl[:, ax] < lo_r[ax] - 1.0) | (cl[:, ax] > hi_r[ax] + 1.0)).sum() for ax in open_axes))
    log(f"=== strand orientation report: {label or data_file} ===")
    log(f"  surface: {n_pop} junction(s) poke > 1 sigma out of the faces (the 0.2 / 99.8 percentile planes of the junctions); "
        f"outermost junction per 4 x 4 sigma patch of each face sits {rec['mean']:.2f} sigma inside the face on average "
        f"(sd {rec['sd']:.2f}, worst {rec['max']:.1f}); {100 * rec['frac2']:.1f} % of patches recessed > 2 sigma, {100 * rec['frac4']:.1f} % > 4 sigma  <- voids / dents")
    log(f"  strands {len(V)}, functionality histogram {deg}, strand end-to-end {Lr.mean():.3f} +/- {Lr.std():.3f} sigma (min {Lr.min():.3f}, max {Lr.max():.3f})")
    out = {}
    for name, m in (('ALL strands', depth > -1), (f'INTERIOR (> {face_depth:.1f} sigma inside the faces)', depth > face_depth),
                    (f'SKIN (<= {face_depth:.1f} sigma)', depth <= face_depth)):
        if m.sum() < 50:
            continue
        met = isotropy_metrics(V[m])
        refs = [isotropy_metrics(rng.normal(size=(int(m.sum()), 3))) for _ in range(n_ref)]
        ref = {k: (np.mean([r[k] for r in refs]), np.std([r[k] for r in refs])) for k in ('aniso4', 'Q4', 'Q6', 'A_xz')}
        checks = {k: abs(met[k] - ref[k][0]) <= 3 * ref[k][1] + 1e-3 for k in ('aniso4', 'Q4', 'Q6', 'A_xz')}
        log(f"  -- {name}: N = {int(m.sum())}")
        log(f"     fabric <n n> eigenvalues  {np.round(met['fabric_eig'], 4)}   (isotropic 1/3 each; ANY cubic lattice also gives 1/3)")
        log(f"     4th-rank anisotropy |<nnnn> - iso| / |iso|  {met['aniso4']:.4f}   (random, same N: {ref['aniso4'][0]:.4f} +/- {ref['aniso4'][1]:.4f}) {'ok' if checks['aniso4'] else 'FAIL'}")
        log(f"     Steinhardt Q4 / Q6  {met['Q4']:.4f} / {met['Q6']:.4f}   (random: {ref['Q4'][0]:.4f} +/- {ref['Q4'][1]:.4f} / {ref['Q6'][0]:.4f} +/- {ref['Q6'][1]:.4f}) {'ok' if checks['Q4'] else 'FAIL'} / {'ok' if checks['Q6'] else 'FAIL'}")
        log(f"     <nx^4> {met['nx4']:.4f}  <nx^2 ny^2> {met['nx2ny2']:.4f}  <nx^2 nz^2> {met['nx2nz2']:.4f}   (isotropic 0.2000, 0.0667, 0.0667)")
        log(f"     Born-affine Zener ratio A_xz = C44 / ((C11 - C12)/2), box frame  {met['A_xz']:.3f}   (isotropic 1, random: {ref['A_xz'][0]:.3f} +/- {ref['A_xz'][1]:.3f}; aligned diamond inf, rot45 0) {'ok' if checks['A_xz'] else 'FAIL'}")
        out[name.split()[0]] = (met, all(checks.values()))
    passed = out.get('INTERIOR', out.get('ALL'))[1]
    log(f"  ISOTROPY {'OK' if passed else 'NOT reached'} (verdict = the interior strands, every metric within 3 sd of random directions)")
    return out.get('ALL')[0], out.get('INTERIOR', out.get('ALL'))[0], passed

def generate_gel_slab(beads_per_chain=5, units_x=8, units_y=8, units_z=4,
                                   solvent_density=0.6,
                                   side_thickness=5.0, below_thickness=0.0, top_thickness=None,
                                   support_thickness=2.0,
                                   support_spacing=0.4,
                                   piston_margin=0.0,
                                   piston_exclusion=0.2,
                                   piston_height=0.5,
                                   output_file="gel_slab.data",
                                   pre_swell=1.0, include_walls=True, rot45_y=False,
                                   bond_switch=0.0, bond_switch_seed=7, bond_switch_keep_size=True,
                                   switch_options=None):
    """
    Generate a tetrahedral (diamond lattice) polymer gel slab in solvent bath.

    Parameters:
    - beads_per_chain: number of beads per chain
    - units_x, units_y: number of unit cells in x and y (horizontal)
    - units_z: number of unit cells in z (height, smaller)
    - solvent_density: solvent density throughout entire box (beads/sigma^3)
    - side_thickness: thickness of solvent layer on the x/y sides of the gel (in sigma)
    - below_thickness: thickness of solvent layer below the support (in sigma).
                       Layout bottom->top: below_thickness solvent, support layer,
                       supportThickness gap, gel, top_thickness solvent.
    - top_thickness: thickness of solvent layer above the gel to the top of the box
                     (in sigma); defaults to side_thickness if not set
    - support_thickness: controls the support's distance from the gel (in sigma)
    - support_spacing: spacing between support atoms (smaller = denser)
    - piston_margin: extra extension of piston (and support) beyond the gel+solvent box in x and y (in sigma);
                     the vacuum gap created at the edges fills with solvent during equilibration
    - piston_exclusion: solvent exclusion half-width around piston plane (in sigma)
    - piston_height: piston position as a fraction (0..1) of the way from the top of
                     the gel to the top of the box (0 = at gel top, 1 = at box top)
    - output_file: LAMMPS data file name
    - include_walls: False writes NO support (type 4) and NO piston (type 5) atoms -- an isolated gel
                     free in its bath (2026-09-30).  The file still declares 5 atom types so
                     slab_with_support.lmp's pair_coeff / group lines stay valid (the groups are empty).
                     Without walls there is no wall tension to drive the box shape under `nph aniso`.
                     piston_margin is ignored (see below); the bottom solvent pad is still
                     below_thickness + 2*support_thickness.
    - rot45_y: True (2026-10-25) builds the SAME diamond network rotated 45 degrees about y, so
               the cube axes of the lattice point along the box diagonals of the xz plane
               (box x = crystal [101], box z = crystal [-101], y = [010]).  The rotated lattice
               is periodic in the box frame with the body-centred-tetragonal cell
               a/sqrt(2) x a x a/sqrt(2): units_x and units_z then count THOSE cells
               (3.73 sigma at a = 5.28), units_y still counts cubic cells (a).  Four crosslinks
               per BCT cell (eight per cubic cell), so 13 x 9 x 30 rotated cells hold about as
               many as the 9 x 9 x 21 block.  Built for the shear anisotropy test: an xz shear
               of the box is a cube-face shear (C44) of the aligned gel and a cube-diagonal
               shear ((C11 - C12)/2, the compression family's G) of the rotated one.  The
               default (False) path is unchanged, IDs byte-identical to before.
    - bond_switch: > 0 (2026-10-25) randomises the network by Wooten-Winer-Weaire bond switching
               after the lattice strands are found: that many ACCEPTED switches per strand (3 is
               plenty; the anneal usually accepts more).  Functionality, strand count, strand
               length and the block faces are kept; the junctions are relaxed between sweeps so
               no strand is longer than 1.9 r0 (4.3 sigma end-to-end, bead spacing < 1.1 sigma against FENE R0 = 1.5).  The
               result has no lattice directions left: C44 = (C11 - C12)/2, the shear G is the
               same for every orientation.  Check it with verify_network_isotropy(<file>).
               switch_options (2026-11-03): dict of extra www_bond_switch keywords, the
               junction-relaxation knobs (rest_len, rep_cut, w_rep, ...; see RANDOM_RELAX).
               bond_switch_keep_size False lets the block scale relax as well (the random
               network's own size, ~0.8 of the lattice block: strands then sit at ~1.0 r0
               instead of ~1.17 r0; the default keeps the block the size the lattice had).
    """

    # An isolated gel (include_walls=False) has no sheet to overhang, so a non-zero
    # piston_margin would only write a VACUUM frame of that width around the x/y edges
    # of the box (solvent is placed in 0..box_x, the box is written -margin..box_x+margin).
    # The _03_nowalls file of 2026-09-30 carries such a 5-sigma frame (xlo = -5).
    # Ignore the margin without sheets (2026-10-11).
    if not include_walls and piston_margin:
        print(f"  include_walls=False: ignoring piston_margin={piston_margin} (no sheets to overhang)")
        piston_margin = 0.0

    # Lattice parameters
    bead_spacing = 1.2  # Bond length
    chain_length = bead_spacing * (beads_per_chain - 1)
    z_clearance = 0  # Adjust this value as needed

    # Diamond lattice: use chain_length as the FCC lattice constant
    # pre_swell (<1) starts the network pre-deswollen near its equilibrium size,
    # shrinking the slow collective-contraction distance. 1.0 = original behaviour.
    a = chain_length * pre_swell  # FCC lattice constant

    # Gel dimensions (rot45_y: BCT cell a/sqrt(2) along x and z, see the docstring)
    cell_xz = a / np.sqrt(2) if rot45_y else a
    gel_x = units_x * cell_xz
    gel_y = units_y * a
    gel_z = units_z * cell_xz

    # Box dimensions (gel + solvent thickness + support)
    # top_thickness defaults to side_thickness if not explicitly set
    if top_thickness is None:
        top_thickness = side_thickness

    box_x = gel_x + 2 * side_thickness
    box_y = gel_y + 2 * side_thickness

    # Offset to center gel in x-y, place above support in z.
    # z-layering (bottom -> top): below_thickness of solvent, support layer,
    # supportThickness gap, gel, then top_thickness of solvent. The gel sits
    # 2*support_thickness above the below-solvent base, so its distance from the
    # support plane (at support_thickness/2) stays controlled by support_thickness.
    offset = np.array([side_thickness, side_thickness,
                       below_thickness + 2 * support_thickness + z_clearance])

    # Box height: gel top plus the top solvent thickness
    gel_top_z = offset[2] + gel_z
    box_z = gel_top_z + top_thickness

    # Piston z: halfway between top of gel and top of box
    piston_z  = gel_top_z + piston_height * (box_z - gel_top_z)

    particles = []
    bonds = []
    particle_id = 1
    bond_id = 1
    molecule_id = 1

    # Generate diamond lattice crosslinks
    # Diamond = 2 interpenetrating FCC lattices offset by (a/4, a/4, a/4)
    crosslinks = {}

    fcc1_positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.5, 0.5, 0.0]),
        np.array([0.5, 0.0, 0.5]),
        np.array([0.0, 0.5, 0.5])
    ]

    if rot45_y:
        # Rotated lattice (2026-10-25): place the conventional diamond cells on a cubic
        # index range big enough to cover the block, rotate every site 45 degrees about
        # y (crystal x -> box (1,0,-1)/sqrt2, crystal z -> box (1,0,1)/sqrt2), shift the
        # cloud so that a lattice site sits at the block's lower corner, and keep the
        # sites inside the closed block [0, gel] (the same closed interval as the aligned
        # build, with a tolerance for the rounding of the rotation).  The block faces
        # are then (101)-type lattice planes, the bonds are found by distance below as
        # for the aligned lattice, and the keys carry the cubic indices.
        c, s = np.sqrt(0.5), np.sqrt(0.5)
        R = np.array([[c, 0.0, s],
                      [0.0, 1.0, 0.0],
                      [-s, 0.0, c]])
        tol = 1e-6 * a
        half = int(np.ceil((gel_x + gel_z) / a)) + 2
        for i in range(-half, half + 1):
            for j in range(units_y + 1):
                for k in range(-half, half + 1):
                    for sublattice in [0, 1]:
                        for fcc_idx, fcc_frac in enumerate(fcc1_positions):
                            frac_pos = fcc_frac + (0.25 if sublattice == 1 else 0.0)
                            rel = R @ ((np.array([i, j, k]) + frac_pos) * a)
                            if (-tol <= rel[0] <= gel_x + tol and
                                -tol <= rel[1] <= gel_y + tol and
                                -tol <= rel[2] <= gel_z + tol):
                                key = (i, j, k, sublattice, fcc_idx)
                                crosslinks[key] = particle_id
                                particles.append({
                                    'id': particle_id,
                                    'type': 1,  # Crosslink
                                    'pos': rel + offset,
                                    'mol': 0
                                })
                                particle_id += 1
    else:
      for i in range(units_x + 1):
        for j in range(units_y + 1):
            for k in range(units_z + 1):
                for sublattice in [0, 1]:
                    for fcc_idx, fcc_frac in enumerate(fcc1_positions):
                        if sublattice == 1:
                            frac_pos = fcc_frac + np.array([0.25, 0.25, 0.25])
                        else:
                            frac_pos = fcc_frac

                        pos = (np.array([i, j, k]) + frac_pos) * a + offset

                        if (0 <= pos[0] - offset[0] <= gel_x and
                            0 <= pos[1] - offset[1] <= gel_y and
                            0 <= pos[2] - offset[2] <= gel_z):

                            key = (i, j, k, sublattice, fcc_idx)
                            crosslinks[key] = particle_id
                            particles.append({
                                'id': particle_id,
                                'type': 1,  # Crosslink
                                'pos': pos.copy(),
                                'mol': 0
                            })
                            particle_id += 1

    # Generate tetrahedral bonds.
    # Nearest-neighbor crosslinks sit at distance a*sqrt(3)/4. Instead of the old
    # O(N^2) double loop (~1e8 distance checks for a 9x9x21 slab -> minutes of
    # silent runtime), bin crosslinks into a spatial grid (cell list) and only
    # test the 27 neighbouring cells. Pairs are emitted in the same ascending-id
    # order as the original loop, so atom/bond/molecule IDs are byte-identical.
    from collections import defaultdict

    bond_distance = a * np.sqrt(3) / 4
    tol = 0.1 * bond_distance
    r = bond_distance + tol
    r2_lo = (bond_distance - tol) ** 2
    r2_hi = (bond_distance + tol) ** 2

    cl_ids = list(crosslinks.values())
    cl_pos = np.array([particles[i - 1]['pos'] for i in cl_ids])
    print(f"  Crosslinks placed: {len(cl_ids)}; finding bonds...")

    cell_keys = np.floor(cl_pos / r).astype(np.int64)
    cells = defaultdict(list)
    for idx in range(len(cl_ids)):
        cells[(cell_keys[idx, 0], cell_keys[idx, 1], cell_keys[idx, 2])].append(idx)

    cell_offsets = [(dx, dy, dz)
                    for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]

    neighbors = defaultdict(list)
    for idx in range(len(cl_ids)):
        id1 = cl_ids[idx]
        p1 = cl_pos[idx]
        cx, cy, cz = cell_keys[idx]
        for dx, dy, dz in cell_offsets:
            for jdx in cells.get((cx + dx, cy + dy, cz + dz), ()):
                id2 = cl_ids[jdx]
                if id2 <= id1:
                    continue
                diff = cl_pos[jdx] - p1
                d2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]
                if r2_lo < d2 < r2_hi:
                    neighbors[id1].append(id2)
    for id1 in neighbors:
        neighbors[id1].sort()

    if bond_switch > 0:
        # Random f = 4 network (2026-10-25): re-wire the lattice strands by bond switching
        # (functionality and strand count kept), relax the junctions, write the moved
        # junction positions back and rebuild the neighbour lists in the same ascending-id
        # form the chain builder reads.
        id_to_idx = {cid: k for k, cid in enumerate(cl_ids)}
        edges0 = np.array([(id_to_idx[a], id_to_idx[b]) for a in neighbors for b in neighbors[a]], dtype=np.int64)
        lo = offset; hi = offset + np.array([gel_x, gel_y, gel_z])
        cl_pos, edges_new, switch_report = www_bond_switch(
            cl_pos, edges0, bond_distance, lo, hi, switches_per_strand=bond_switch,
            seed=bond_switch_seed, keep_size=bond_switch_keep_size, **(switch_options or {}))
        for k, cid in enumerate(cl_ids):
            particles[cid - 1]['pos'] = cl_pos[k].copy()
        neighbors = defaultdict(list)
        for i1, i2 in edges_new:
            a, b = sorted((cl_ids[i1], cl_ids[i2]))
            neighbors[a].append(b)
        for id1 in neighbors:
            neighbors[id1].sort()

    # Build chains for each bonded crosslink pair, in ascending-id order.
    for id1 in cl_ids:
        pos1 = particles[id1 - 1]['pos']
        for id2 in neighbors.get(id1, ()):
            pos2 = particles[id2 - 1]['pos']

            chain_ids = [id1]

            for b in range(1, beads_per_chain - 1):
                frac = b / (beads_per_chain - 1)
                pos = pos1 * (1 - frac) + pos2 * frac
                particles.append({
                    'id': particle_id,
                    'type': 2,  # Chain bead
                    'pos': pos.copy(),
                    'mol': molecule_id
                })
                chain_ids.append(particle_id)
                particle_id += 1

            chain_ids.append(id2)

            particles[id1 - 1]['mol'] = molecule_id
            particles[id2 - 1]['mol'] = molecule_id

            for b in range(len(chain_ids) - 1):
                bonds.append({
                    'id': bond_id,
                    'type': 1,
                    'atom1': chain_ids[b],
                    'atom2': chain_ids[b + 1]
                })
                bond_id += 1

            molecule_id += 1

    # Drop crosslink sites no chain reaches (2026-09-29).  The cube keeps every lattice
    # site with 0 <= r <= L, including corner sites whose four neighbours all lie outside
    # it: those beads carry no bond, so they are free monomers that wander through the
    # bath.  They are polymer by TYPE, so every bounding-box measure of the gel (the
    # gel_dimensions output of slab_with_support.lmp, the plate seating of
    # compress_slab.lmp) follows them instead of the network.  4 of them in every cube.
    bonded = set()
    for b in bonds:
        bonded.add(b['atom1'])
        bonded.add(b['atom2'])
    n_before = len(particles)
    particles = [p for p in particles if p['id'] in bonded]
    n_unbonded = n_before - len(particles)
    if n_unbonded:
        remap = {p['id']: i + 1 for i, p in enumerate(particles)}
        for p in particles:
            p['id'] = remap[p['id']]
        for b in bonds:
            b['atom1'], b['atom2'] = remap[b['atom1']], remap[b['atom2']]
        crosslinks = {k: remap[v] for k, v in crosslinks.items() if v in remap}
        particle_id = len(particles) + 1
        print(f"  Removed {n_unbonded} crosslink site(s) with no chain attached (free beads)")

    # Create dense support layers with hexagonal packing
    # Piston and support tile from -piston_margin to box+piston_margin so the
    # wall footprint is larger than the gel+solvent region; the extra vacuum
    # at the edges fills with solvent during equilibration.
    wall_x0, wall_x1 = -piston_margin, box_x + piston_margin
    wall_y0, wall_y1 = -piston_margin, box_y + piston_margin
    nx = int((wall_x1 - wall_x0) / support_spacing) + 1
    ny = int((wall_y1 - wall_y0) / (support_spacing * np.sqrt(3))) + 1

    # Bottom support layer (type 4 - frozen). Sits above the below-solvent layer.
    for i in (range(nx) if include_walls else ()):
        for j in range(ny):
            x = wall_x0 + i * support_spacing
            y = wall_y0 + j * support_spacing * np.sqrt(3) + (support_spacing / 2 if i % 2 else 0)
            z = below_thickness + support_thickness / 2 + z_clearance

            if wall_x0 <= x <= wall_x1 and wall_y0 <= y <= wall_y1:
                particles.append({
                    'id': particle_id,
                    'type': 4,  # Bottom support
                    'pos': np.array([x, y, z]),
                    'mol': molecule_id
                })
                particle_id += 1
                molecule_id += 1

    # Top piston layer (type 5 - mobile, for compression)
    # Placed halfway between the top of the gel and the top of the box so that
    # the feed reservoir (above piston) and permeate region (below piston / gel)
    # are well-defined from the initial data file.
    for i in (range(nx) if include_walls else ()):
        for j in range(ny):
            x = wall_x0 + i * support_spacing
            y = wall_y0 + j * support_spacing * np.sqrt(3) + (support_spacing / 2 if i % 2 else 0)
            z = piston_z

            if wall_x0 <= x <= wall_x1 and wall_y0 <= y <= wall_y1:
                particles.append({
                    'id': particle_id,
                    'type': 5,  # Top piston
                    'pos': np.array([x, y, z]),
                    'mol': molecule_id
                })
                particle_id += 1
                molecule_id += 1

    # Add solvent uniformly throughout entire box
    total_volume = box_x * box_y * box_z
    num_solvent_total = int(solvent_density * total_volume)
    print(f"  Placing up to {num_solvent_total} solvent beads...")

    np.random.seed(42)
    solvent_count = 0
    max_attempts = num_solvent_total * 100
    attempts = 0

    while solvent_count < num_solvent_total and attempts < max_attempts:
        attempts += 1
        pos = np.array([
            np.random.rand() * box_x,
            np.random.rand() * box_y,
            np.random.rand() * box_z
        ])

        # Exclude solvent within piston_exclusion*sigma of the piston plane
        if include_walls and abs(pos[2] - piston_z) < piston_exclusion:
            continue

        min_dist = 0.8
        too_close = False
        for p in particles[-min(100, len(particles)):]:
            if np.linalg.norm(pos - p['pos']) < min_dist:
                too_close = True
                break

        if not too_close:
            particles.append({
                'id': particle_id,
                'type': 3,  # Solvent
                'pos': pos.copy(),
                'mol': molecule_id
            })
            particle_id += 1
            molecule_id += 1
            solvent_count += 1

    if solvent_count < num_solvent_total:
        print(f"Warning: Only placed {solvent_count}/{num_solvent_total} solvent particles")

    # Write LAMMPS data file
    with open(output_file, 'w') as f:
        tag = ""
        if bond_switch > 0:          # the relaxation knobs in the header (see submit_random_periodic_pipeline.sh)
            tag = " (random f=4 bond-switched, relax " + " ".join(f"{k}={v}" for k, v in (switch_options or {}).items()) + ")"
        f.write(f"LAMMPS data file for tetrahedral gel slab with piston{tag}\n\n")
        f.write(f"{len(particles)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write("0 angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n\n")
        f.write("5 atom types\n")
        f.write("1 bond types\n\n")
        lo = -piston_margin if piston_margin else 0.0   # no "-0.000000" when the margin is 0.0
        f.write(f"{lo:.6f} {box_x + piston_margin:.6f} xlo xhi\n")
        f.write(f"{lo:.6f} {box_y + piston_margin:.6f} ylo yhi\n")
        f.write(f"0.0 {box_z} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write("1 1.0  # Crosslink\n")
        f.write("2 1.0  # Chain bead\n")
        f.write("3 1.0  # Solvent\n")
        f.write("4 1.0  # Bottom support (frozen)\n")
        f.write("5 1.0  # Top piston (mobile)\n\n")
        f.write("Atoms\n\n")

        for p in particles:
            f.write(f"{p['id']} {p['mol']} {p['type']} {p['pos'][0]:.6f} {p['pos'][1]:.6f} {p['pos'][2]:.6f}\n")

        f.write("\nBonds\n\n")
        for b in bonds:
            f.write(f"{b['id']} {b['type']} {b['atom1']} {b['atom2']}\n")

    num_solvent_final = sum(1 for p in particles if p['type'] == 3)
    num_polymer = sum(1 for p in particles if p['type'] in [1, 2])
    num_bottom_support = sum(1 for p in particles if p['type'] == 4)
    num_top_piston = sum(1 for p in particles if p['type'] == 5)
    num_crosslinks = len(crosslinks)
    num_chains = len(bonds) // (beads_per_chain - 1) if beads_per_chain > 1 else 0

    print("Generated tetrahedral gel slab with piston:")
    print(f"  Unit cells: {units_x} x {units_y} x {units_z}" +
          (f"  (rot45_y: x, z count BCT cells of a/sqrt2 = {cell_xz:.3f}; lattice rotated 45 deg about y)" if rot45_y else ""))
    if bond_switch > 0:
        print(f"  Bond switching: {switch_report['accepted']} switches ({switch_report['per_strand']:.2f} per strand), "
              f"strands {switch_report['mean_len']:.3f} +/- {switch_report['std_len']:.3f} r0 (max {switch_report['max_len']:.3f}), "
              f"{100 * switch_report['lattice_strands_left']:.1f} % lattice strands left, block scale {switch_report['block_scale']:.3f}")
    print(f"  Beads per chain: {beads_per_chain}")
    print(f"  Gel dimensions: {gel_x:.2f} x {gel_y:.2f} x {gel_z:.2f}")
    print(f"  Box dimensions: {box_x:.2f} x {box_y:.2f} x {box_z:.2f}")
    print(f"  Side thickness (x/y):  {side_thickness:.2f}")
    print(f"  Below thickness (z):   {below_thickness:.2f}")
    print(f"  Top thickness (z):     {top_thickness:.2f}")
    print(f"  Gel top z:          {gel_top_z:.4f}")
    print(f"  Piston z:           {piston_z:.4f}  (piston_height={piston_height:.2f} of gel top -> box top)")
    print(f"  Piston exclusion:   +/-{piston_exclusion:.2f} sigma around piston plane")
    print(f"  Crosslinks: {num_crosslinks}")
    print(f"  Chains: {num_chains}")
    print(f"  Polymer beads: {num_polymer}")
    if include_walls:
        print(f"  Bottom support: {num_bottom_support} (spacing: {support_spacing:.2f}, margin: {piston_margin:.2f})")
        print(f"  Top piston: {num_top_piston} (spacing: {support_spacing:.2f}, margin: {piston_margin:.2f})")
    else:
        print("  NO support / piston atoms (include_walls=False); 5 atom types still declared")
    print(f"  Solvent beads: {num_solvent_final} (density: {solvent_density:.2f})")
    print(f"  Total atoms: {len(particles)}")
    print(f"  Total bonds: {len(bonds)}")
    print(f"  Output: {output_file}")


# ---------------------------------------------------------------------------
# THE KNOB: which network, with the presets the three current gels were built with
# ---------------------------------------------------------------------------
# Junction-relaxation knobs of the random presets (2026-11-03, see the note above
# www_bond_switch): strand rest length 1.4 r0, non-bonded junction repulsion out to 1.3 r0
# at weight 3.  The old defaults (1.0 r0, 0.9 r0, weight 1; switch_options={}) clumped the
# junctions and left 3-5 sigma voids inside the as-built gel.
RANDOM_RELAX = dict(rest_len=1.4, rep_cut=1.3, w_rep=3.0)

NETWORKS = {
    # the aligned diamond lattice: units count CUBIC cells (a = 4.8 pre_swell sigma)
    'tetrahedral': dict(rot45_y=False, bond_switch=0.0),
    # rotated 45 deg about y: units_x / units_z count BCT cells of a / sqrt2 (3.73 sigma at a = 5.28)
    'rot45':       dict(rot45_y=True, bond_switch=0.0),
    # the rot45 block randomised by bond switching (3 accepted switches per strand, seed 7)
    'random':      dict(rot45_y=True, bond_switch=3.0, bond_switch_seed=7, switch_options=RANDOM_RELAX),
}


def build(network='tetrahedral', output_file=None, beads_per_chain=5, units_x=9, units_y=9, units_z=21,
          solvent_density=0.40, side_thickness=15.0, below_thickness=7.0, top_thickness=15.0,
          support_thickness=4.0, support_spacing=0.2, piston_margin=0.0, piston_exclusion=0.1,
          piston_height=0.8, pre_swell=1.10, include_walls=False, check_isotropy=None, **overrides):
    """Build one gel data file.  `network` picks the preset ('tetrahedral' | 'rot45' | 'random');
    everything else is generate_gel_slab's knob of the same name (the defaults are the tall
    isolated rho04 gel's: 9 x 9 x 21 cubic cells, or 13 x 9 x 30 BCT cells for the rotated
    builds -- pass units_x / units_z).  `overrides` go straight to generate_gel_slab (e.g.
    bond_switch=2.0, bond_switch_seed=11, bond_switch_keep_size=False).  check_isotropy: run
    verify_network_isotropy on the written file (default: only for 'random').  Returns the output path."""
    if network not in NETWORKS:
        raise ValueError(f"network must be one of {list(NETWORKS)} (got {network!r})")
    if output_file is None:
        raise ValueError('output_file is required (e.g. ../../lammps_data_files_local/<name>.data)')
    kw = dict(NETWORKS[network])
    kw.update(overrides)
    print(f"=== data_file_generator.build: network = {network!r}  ->  {output_file}")
    generate_gel_slab(beads_per_chain, units_x, units_y, units_z,
                      solvent_density, side_thickness, below_thickness, top_thickness,
                      support_thickness, support_spacing,
                      piston_margin, piston_exclusion, piston_height, output_file,
                      pre_swell=pre_swell, include_walls=include_walls, **kw)
    if check_isotropy is None:
        check_isotropy = (network == 'random')
    if check_isotropy:
        print()
        verify_network_isotropy(output_file, f"{network} network")
        verify_density_uniformity(output_file, f"{network} network", scale=1.2 * (beads_per_chain - 1) * pre_swell / 4.464)
    return output_file
