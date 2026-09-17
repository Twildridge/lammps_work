#!/usr/bin/env python3
"""
make_fixtures.py -- synthetic run directories for testing the plotting scripts
and the analysis notebooks (2026-09-16).  Generates SMALL but format-faithful
output trees:

  <root>/one_piston/   a triaxial_permeation-style run (single piston columns)
  <root>/two_piston_perm/   a triaxial_permeation_two_pist run (feed | perm columns,
                            piston_pressure, permeation_data/permeation)
  <root>/two_piston_comp/   a triaxial_compression_two_pist run, levels c0.05 c0.10
                            (load | feed | perm columns, _ref files, disp_z_polymer, ...)
  <root>/flow_data_local/{compression,permeation}/<RUN_ID>/  +  traj_files.nosync/
                            the same two-piston files laid out as the notebooks expect
                            (flat, plus a 2-frame traj_ref / traj_stress for the box header)

Every file uses the exact header/column layout the decks write, so a reader that
works here works on cluster output.  Usage:  python make_fixtures.py <root>
"""
import os
import sys
import numpy as np
from pathlib import Path

rng = np.random.default_rng(7)
DATANAME = 'fixture_slab'
INTER = '1.0_1.0'
LX, LY, LZ = 40.0, 40.0, 200.0
Z_SUP, Z_GEL_LO, Z_GEL_HI, Z_LOAD, Z_FEED, Z_PERM = 30.0, 45.0, 165.0, 168.0, 185.0, 20.0
BW = 2.0
NB = int(LZ / BW)
ZC = np.arange(NB) * BW + BW / 2


def stem(nsteps, lvl=None):
    return f'{DATANAME}_{INTER}_{nsteps}' + (f'_c{lvl}' if lvl else '')


def mk(root):
    root = Path(root)
    for sub in ('stress_data', 'volume_data', 'piston_data', 'permeation_data', 'displacement_data',
                'pair_data', 'chemical_potential'):
        (root / 'output_files' / sub).mkdir(parents=True, exist_ok=True)
    (root / 'traj_files').mkdir(exist_ok=True)
    (root / 'output_plots').mkdir(exist_ok=True)
    return root


def write_print(path, header, rows, fmt='%.6g'):
    with open(path, 'w') as f:
        if header:
            f.write(header + '\n')
        for r in rows:
            f.write(' '.join(fmt % v for v in r) + '\n')


def write_ave_time_vector(path, fixname, colname, steps, profiles):
    with open(path, 'w') as f:
        f.write(f'# Time-averaged data for fix {fixname}\n# TimeStep Number-of-rows\n# Row {colname}\n')
        for t, prof in zip(steps, profiles):
            f.write(f'{int(t)} {len(prof)}\n')
            for i, v in enumerate(prof, 1):
                f.write(f'{i} {v:.6g}\n')


def write_ave_chunk(path, fixname, steps, cols, ncol_names):
    with open(path, 'w') as f:
        f.write(f'# Chunk-averaged data for fix {fixname} and group all\n# Timestep Number-of-chunks Total-count\n'
                f'# Chunk Coord1 Ncount {" ".join(ncol_names)}\n')
        for t, C in zip(steps, cols):
            f.write(f'{int(t)} {NB} {C[:, 0].sum():.6g}\n')
            for i in range(NB):
                f.write(f'  {i + 1} {ZC[i]:.4f} ' + ' '.join(f'{v:.6g}' for v in C[i]) + '\n')


def gel_mask(zc=ZC, lo=Z_GEL_LO, hi=Z_GEL_HI):
    return (zc > lo) & (zc < hi)


def stress_profiles(steps, eps=0.0, dp=0.0, ref=False):
    """polymer + solvent sigma_zz/xx/yy profiles: reservoir at 1.5 (+dp on the feed side),
    polymer partial ~0.4 in the gel, network stress ~ 0.3*eps."""
    out = {}
    g = gel_mask()
    feed = ZC > Z_GEL_HI
    for comp in ('zz', 'xx', 'yy'):
        P, S = [], []
        for k, t in enumerate(steps):
            relax = 1.0 if ref else (1.0 - 0.5 * np.exp(-k / 3.0))
            p = np.where(g, 0.42 + (0.3 if comp == 'zz' else 0.12) * eps * relax, 0.0) + rng.normal(0, 0.01, NB)
            sv = np.where(g, 1.08, 1.5) + np.where(feed, dp, 0.0) + rng.normal(0, 0.01, NB)
            sv[(ZC < Z_PERM) | (ZC > Z_FEED)] = 0.0
            p[(ZC < Z_PERM) | (ZC > Z_FEED)] = 0.0
            P.append(p)
            S.append(sv)
        out[comp] = (np.array(P), np.array(S))
    return out


def density_cols(steps, eps=0.0):
    g = gel_mask()
    cols = []
    for k, t in enumerate(steps):
        n = np.where(g, 0.19 - 0.02 * eps, 0.43) + rng.normal(0, 0.005, NB)
        n[(ZC < Z_PERM) | (ZC > Z_FEED)] = 0.0
        cols.append(np.column_stack([n * LX * LY * BW, n, n]))
    return cols


def write_traj(path, steps, types_z):
    """2-frame lammpstrj with the wall sheets (4,5,6,7) and a few polymer/solvent atoms."""
    atoms = []
    aid = 1
    for t, zs in types_z.items():
        for _ in range(8):
            atoms.append((aid, t, rng.uniform(0, LX), rng.uniform(0, LY), zs))
            aid += 1
    for _ in range(40):
        atoms.append((aid, 2, rng.uniform(0, LX), rng.uniform(0, LY), rng.uniform(Z_GEL_LO, Z_GEL_HI)))
        aid += 1
    for _ in range(40):
        atoms.append((aid, 3, rng.uniform(0, LX), rng.uniform(0, LY), rng.uniform(Z_PERM + 1, Z_FEED - 1)))
        aid += 1
    with open(path, 'w') as f:
        for t in steps:
            f.write(f'ITEM: TIMESTEP\n{int(t)}\nITEM: NUMBER OF ATOMS\n{len(atoms)}\nITEM: BOX BOUNDS pp pp pp\n'
                    f'0.0 {LX}\n0.0 {LY}\n0.0 {LZ}\nITEM: ATOMS id type mol x y z\n')
            for a in atoms:
                f.write(f'{a[0]} {a[1]} {a[0]} {a[2]:.4f} {a[3]:.4f} {a[4]:.4f}\n')


# ---------------------------------------------------------------------------
def one_piston(root):
    """triaxial_permeation-style (single piston column; piston_force_pressure)."""
    r = mk(root)
    ns = 8000
    st = stem(ns)
    steps = np.arange(0, ns + 1, 500)
    sd, pd_, vd, pm = (r / 'output_files' / s for s in ('stress_data', 'piston_data', 'volume_data', 'permeation_data'))
    z = 180 - 0.002 * steps
    write_print(pd_ / f'piston_position_{st}.dat', '# Fix print output for fix out_piston_pos', np.column_stack([steps, z]))
    write_print(pd_ / f'piston_velocity_{st}.dat', '# Fix print output for fix out_piston_vel', np.column_stack([steps, -0.002 / 0.005 + rng.normal(0, 0.05, len(steps))]))
    F = 160 + rng.normal(0, 40, len(steps))
    write_print(pd_ / f'piston_force_{st}.dat', '# Fix print output for fix out_piston_force', np.column_stack([steps, F]))
    write_print(pd_ / f'piston_force_avg_{st}.dat', '# Time-averaged data for fix out_piston_force_avg\n# TimeStep c_piston_fz', np.column_stack([steps, F]))
    write_print(pd_ / f'piston_force_pressure_{st}.dat', None,
                np.column_stack([steps, F, F / (LX * LY), 1.6 + rng.normal(0, 0.01, len(steps)), 1.5 + rng.normal(0, 0.01, len(steps)), 1.6 * LX * LY + 0 * steps]))
    write_print(sd / f'strain_zz_{st}.dat', None, np.column_stack([steps, 120 + 0 * steps, 120 - 0.0002 * steps]))
    write_print(vd / f'box_dimensions_{st}.dat', None, np.column_stack([steps, LX + 0 * steps, LY + 0 * steps, LZ + 0 * steps]))
    write_print(vd / f'gel_dimensions_rg_{st}.dat', None, np.column_stack([steps, 40 + 0 * steps, 40 + 0 * steps, 120 - 0.0002 * steps]))
    write_print(vd / f'gel_volume_rg_{st}.dat', None, np.column_stack([steps, 40 * 40 * (120 - 0.0002 * steps)]))
    write_print(vd / f'gel_volume_bb_{st}.dat', None, np.column_stack([40 * 40 * (122 - 0.0002 * steps)]))
    write_print(pm / f'permeate_count_{st}.dat', '# Fix print output for fix out_permeate_flux', np.column_stack([steps, 5000 + 0.05 * steps]))
    for f_, v in (('pressure_feed', 1.6), ('pressure_permeate', 1.5)):
        write_print(sd / f'{f_}_{st}.dat', f'# Time-averaged data for fix avg_P\n# TimeStep v_P', np.column_stack([steps, v + rng.normal(0, 0.005, len(steps))]))
    pro = stress_profiles(steps, dp=0.1)
    for comp, (P, S) in pro.items():
        write_ave_time_vector(sd / f'sigma{comp}_polymer_{st}.dat', 'avg', f'c_sigma{comp}_poly', steps, P)
        write_ave_time_vector(sd / f'sigma{comp}_solvent_{st}.dat', 'avg', f'c_sigma{comp}_solv', steps, S)
    for dim in ('x', 'y', 'z'):
        for comp in ('polymer', 'solvent', 'piston', 'support'):
            if comp in ('piston', 'support') and dim != 'z':
                continue
            nb = int((LX if dim == 'x' else LY if dim == 'y' else LZ) / BW)
            prof = [np.abs(rng.normal(0.5, 0.05, nb)) for _ in steps]
            write_ave_time_vector(sd / f'stress_{dim}_{comp}_{st}.dat', 'avg', 'c_press', steps, prof)
    write_ave_chunk(r / 'output_files' / 'chemical_potential' / f'solvent_density_z_{st}.dat', 'avg_solv_density', steps, density_cols(steps), ['density/number'])
    with open(r / f'final_triperm_{st}.data', 'w') as f:
        f.write(f'fixture\n\n10 atoms\n5 atom types\n\n0.0 {LX} xlo xhi\n0.0 {LY} ylo yhi\n0.0 {LZ} zlo zhi\n')
    (r / 'log.lammps').write_text(_log(steps))
    return st


def _log(steps):
    lines = ['LAMMPS (22 Jul 2025)', '   Step     c_mobile_temp      PotEng         KinEng         TotEng         Press      c_mobile_press     Volume           Lx             Ly             Lz      ']
    for t in steps:
        lines.append(f'{int(t):8d} {1 + rng.normal(0, 0.005):.6f} 9.0 1.1 10.1 {1.0 + rng.normal(0, 0.02):.6f} 0.58 {LX * LY * LZ:.2f} {LX} {LY} {LZ}')
    lines.append('Loop time of 1 on 1 procs')
    lines.append('Total wall time: 0:00:01')
    return '\n'.join(lines) + '\n'


# ---------------------------------------------------------------------------
def two_piston_perm(root):
    r = mk(root)
    ns = 8000
    st = stem(ns)
    steps = np.arange(0, ns + 1, 500)
    sd, pd_, vd, pm, cp = (r / 'output_files' / s for s in ('stress_data', 'piston_data', 'volume_data', 'permeation_data', 'chemical_potential'))
    n = len(steps)
    zf = Z_FEED - 0.0015 * steps
    zp = Z_PERM - 0.0015 * steps
    write_print(pd_ / f'piston_position_{st}.dat', '# step z_feed z_perm', np.column_stack([steps, zf, zp]))
    write_print(pd_ / f'piston_velocity_{st}.dat', '# step vz_feed vz_perm', np.column_stack([steps, -0.3 + rng.normal(0, 0.02, n), -0.3 + rng.normal(0, 0.02, n)]))
    Ff = 1.6 * LX * LY + rng.normal(0, 30, n)
    Fp = -1.5 * LX * LY + rng.normal(0, 30, n)
    write_print(pd_ / f'piston_force_{st}.dat', '# step F_fluid_feed F_fluid_perm   (pair force of mobile atoms on each sheet, z)', np.column_stack([steps, Ff, Fp]))
    write_print(pd_ / f'piston_force_avg_{st}.dat', '# Time-averaged data for fix out_piston_force_avg\n# TimeStep F_fluid_feed F_fluid_perm   (block-averaged)', np.column_stack([steps, Ff, Fp]))
    write_print(pd_ / f'piston_force_avg_ref_{st}.dat', '# Time-averaged data for fix ref_piston_force_avg\n# TimeStep F_fluid_feed F_fluid_perm   (zero-flux reference, block-averaged)', np.column_stack([steps[:4], Ff[:4] - 0.1 * LX * LY, Fp[:4]]))
    app = np.where(steps < 1000, 1.5, 1.6)
    write_print(pd_ / f'piston_pressure_{st}.dat', '# step P_feed_meas P_feed_app P_perm_meas P_perm_app',
                np.column_stack([steps, Ff / (LX * LY), app, -Fp / (LX * LY), 1.5 + 0 * steps]))
    Q = LX * LY * 0.0015 / 0.005 + rng.normal(0, 5, n)
    Npm = 8000 + 0.43 * LX * LY * 0.0015 * steps
    write_print(pm / f'permeation_{st}.dat', '# Time-averaged data for fix out_permeation\n# TimeStep z_feed z_perm F_fluid_feed F_fluid_perm P_feed_meas P_perm_meas Q_perm(=A*dz_perm/dt) N_permeate',
                np.column_stack([steps, zf, zp, Ff, Fp, Ff / (LX * LY), -Fp / (LX * LY), Q, Npm]))
    write_print(pm / f'permeate_count_{st}.dat', '# Fix print output for fix out_permeate_flux', np.column_stack([steps, Npm]))
    write_print(sd / f'strain_zz_{st}.dat', None, np.column_stack([steps, 120 + 0 * steps, 120 + rng.normal(0, 0.05, n)]))
    write_print(vd / f'box_dimensions_{st}.dat', None, np.column_stack([steps, LX + 0 * steps, LY + 0 * steps, LZ + 0 * steps]))
    write_print(vd / f'gel_dimensions_rg_{st}.dat', None, np.column_stack([steps, 40 + 0 * steps, 40 + 0 * steps, 120 + rng.normal(0, 0.05, n)]))
    write_print(vd / f'gel_dimensions_bb_{st}.dat', None, np.column_stack([steps, 40 + 0 * steps, 40 + 0 * steps, 122 + rng.normal(0, 0.05, n)]))
    write_print(vd / f'gel_volume_rg_{st}.dat', None, np.column_stack([steps, 40 * 40 * 120 + 0 * steps]))
    write_print(vd / f'gel_volume_bb_{st}.dat', None, np.column_stack([40 * 40 * 122 + 0 * steps]))
    write_print(vd / f'polymer_com_{st}.dat', None, np.column_stack([steps, 20 + 0 * steps, 20 + 0 * steps, 105 + 0 * steps]))
    write_print(sd / f'stress_aniso_{st}.dat', '# Time-averaged data for fix stress_aniso_output\n# TimeStep sig_p_xx sig_p_yy sig_p_zz P_xx P_yy P_zz',
                np.column_stack([steps, 0.4 + rng.normal(0, 0.005, n), 0.4 + rng.normal(0, 0.005, n), 0.4 + rng.normal(0, 0.005, n), 1.0 + 0 * steps, 1.0 + 0 * steps, 1.0 + 0 * steps]))
    for f_, v in (('pressure_feed', 1.6), ('pressure_permeate', 1.5)):
        write_print(sd / f'{f_}_{st}.dat', f'# Time-averaged data for fix avg_P\n# TimeStep P_res rho', np.column_stack([steps, v + rng.normal(0, 0.005, n), 0.43 + 0 * steps]))
    rsteps = np.arange(0, 4) * 500
    pro = stress_profiles(steps, dp=0.1)
    ref = stress_profiles(rsteps, ref=True)
    for comp in ('zz', 'xx', 'yy'):
        write_ave_time_vector(sd / f'sigma{comp}_polymer_{st}.dat', 'avg', f'c_sigma{comp}_poly', steps, pro[comp][0])
        write_ave_time_vector(sd / f'sigma{comp}_solvent_{st}.dat', 'avg', f'c_sigma{comp}_solv', steps, pro[comp][1])
        write_ave_time_vector(sd / f'sigma{comp}_polymer_ref_{st}.dat', 'ref', f'c_sigma{comp}_poly_r', rsteps, ref[comp][0])
        write_ave_time_vector(sd / f'sigma{comp}_solvent_ref_{st}.dat', 'ref', f'c_sigma{comp}_solv_r', rsteps, ref[comp][1])
    for dim in ('x', 'y', 'z'):
        comps = ['polymer', 'solvent'] + (['piston_feed', 'piston_perm', 'support'] if dim == 'z' else [])
        nb = int((LX if dim == 'x' else LY if dim == 'y' else LZ) / BW)
        for comp in comps:
            prof = [np.abs(rng.normal(0.5, 0.05, nb)) for _ in steps]
            write_ave_time_vector(sd / f'stress_{dim}_{comp}_{st}.dat', 'avg', 'c_press', steps, prof)
    write_ave_chunk(cp / f'solvent_density_z_{st}.dat', 'avg_solv_density', steps, density_cols(steps), ['density/number', 'density/mass'])
    write_ave_chunk(cp / f'solvent_density_z_ref_{st}.dat', 'ref_dens', rsteps, density_cols(rsteps), ['density/number', 'density/mass'])
    disp = [np.column_stack([np.where(gel_mask(), 300, 0), np.where(gel_mask(), -0.01 * (ZC - Z_GEL_LO) / 120 * k, 0)]) for k in range(n)]
    write_ave_chunk(r / 'output_files' / 'displacement_data' / f'disp_z_polymer_{st}.dat', 'avg_disp_z_poly', steps, disp, ['v_uz_poly'])
    write_traj(r / 'traj_files' / f'traj_ref_{st}.lammpstrj', rsteps[:2], {4: Z_SUP, 5: Z_FEED, 6: Z_PERM})
    write_traj(r / 'traj_files' / f'traj_stress_{st}.lammpstrj', steps[:2], {4: Z_SUP, 5: Z_FEED, 6: Z_PERM})
    with open(r / f'final_triperm_{st}.data', 'w') as f:
        f.write(f'fixture\n\n10 atoms\n7 atom types\n\n0.0 {LX} xlo xhi\n0.0 {LY} ylo yhi\n0.0 {LZ} zlo zhi\n')
    (r / 'log.lammps').write_text(_log(steps))
    return st


def two_piston_comp(root, levels=('0.05', '0.10')):
    r = mk(root)
    hold = 6000
    sd, pd_, vd, pm, cp, dd = (r / 'output_files' / s for s in ('stress_data', 'piston_data', 'volume_data', 'permeation_data', 'chemical_potential', 'displacement_data'))
    rsteps = np.arange(0, 4) * 500
    ref = stress_profiles(rsteps, ref=True)
    rst = stem(hold)
    for comp in ('zz', 'xx', 'yy'):
        write_ave_time_vector(sd / f'sigma{comp}_polymer_ref_{rst}.dat', 'ref', f'c_sigma{comp}_poly', rsteps, ref[comp][0])
        write_ave_time_vector(sd / f'sigma{comp}_solvent_ref_{rst}.dat', 'ref', f'c_sigma{comp}_solv', rsteps, ref[comp][1])
    write_ave_chunk(cp / f'solvent_density_z_ref_{rst}.dat', 'ref_dens', rsteps, density_cols(rsteps), ['density/number', 'density/mass'])
    Fd0 = 0.002 * LX * LY + rng.normal(0, 20, 4)
    write_print(pd_ / f'piston_force_avg_ref_{rst}.dat', '# Time-averaged data for fix ref_piston_force_avg\n# TimeStep F_load F_fluid_feed F_fluid_perm   (reference preload; block-averaged)',
                np.column_stack([rsteps, Fd0, 1.5 * LX * LY + rng.normal(0, 30, 4), -1.5 * LX * LY + rng.normal(0, 30, 4)]))
    write_traj(r / 'traj_files' / f'traj_ref_{rst}.lammpstrj', rsteps[:2], {4: Z_SUP, 5: Z_FEED, 6: Z_PERM, 7: Z_LOAD})
    write_print(vd / f'box_dimensions_{rst}.dat', None, np.column_stack([rsteps, LX + 0 * rsteps, LY + 0 * rsteps, LZ + 0 * rsteps]))
    write_print(pd_ / f'piston_pressure_{rst}.dat', '# step P_load_meas P_feed_meas P_feed_app P_perm_meas P_perm_app',
                np.column_stack([rsteps, 0.002 + 0 * rsteps, 1.5 + rng.normal(0, 0.01, 4), 1.5 + 0 * rsteps, 1.5 + rng.normal(0, 0.01, 4), 1.5 + 0 * rsteps]))
    write_print(pd_ / f'piston_position_run_{rst}.dat', '# step z_load z_feed z_perm', np.column_stack([rsteps, Z_LOAD + 0 * rsteps, Z_FEED + 0 * rsteps, Z_PERM + 0 * rsteps]))
    t_start = 2000
    L0 = 120.0
    for li, lvl in enumerate(levels):
        eps = float(lvl)
        st = stem(hold, lvl)
        steps = np.arange(t_start, t_start + hold + 1, 300)
        n = len(steps)
        drive = np.minimum(1.0, (steps - t_start) / 800.0)
        zd = Z_LOAD - eps * L0 * drive
        zfeed = Z_FEED + eps * L0 * 0.6 * drive
        zperm = Z_PERM - eps * L0 * 0.4 * drive
        write_print(pd_ / f'piston_position_{st}.dat', '# step z_load z_feed z_perm', np.column_stack([steps, zd, zfeed, zperm]))
        write_print(pd_ / f'piston_velocity_{st}.dat', '# step vz_load vz_feed vz_perm', np.column_stack([steps, -0.04 * (drive < 1), 0 * steps, 0 * steps]))
        Fd = (0.3 * eps * (1 - 0.5 * np.exp(-(steps - t_start) / 1500)) * LX * LY) * (drive >= 1) + rng.normal(0, 15, n)
        Ff = 1.5 * LX * LY + rng.normal(0, 30, n)
        Fp = -1.5 * LX * LY + rng.normal(0, 30, n)
        write_print(pd_ / f'piston_force_{st}.dat', '# step F_load F_fluid_feed F_fluid_perm   (pair force on each sheet, z)', np.column_stack([steps, Fd, Ff, Fp]))
        write_print(pd_ / f'piston_force_avg_{st}.dat', '# Time-averaged data for fix out_piston_force_avg\n# TimeStep F_load F_fluid_feed F_fluid_perm   (block-averaged)', np.column_stack([steps, Fd, Ff, Fp]))
        write_print(pd_ / f'piston_pressure_{st}.dat', '# step P_load_meas P_feed_meas P_feed_app P_perm_meas P_perm_app',
                    np.column_stack([steps, Fd / (LX * LY), Ff / (LX * LY), 1.5 + 0 * steps, -Fp / (LX * LY), 1.5 + 0 * steps]))
        write_print(pm / f'permeation_{st}.dat', '# step z_feed z_perm dV_feed dV_perm dV_total   (solvent expelled since seating, sigma^3; A*dz)',
                    np.column_stack([steps, zfeed, zperm, LX * LY * (zfeed - Z_FEED), -LX * LY * (zperm - Z_PERM), LX * LY * ((zfeed - Z_FEED) - (zperm - Z_PERM))]))
        Lrg = L0 * (1 - eps * drive)
        write_print(sd / f'strain_zz_{st}.dat', None, np.column_stack([steps, L0 + 0 * steps, Lrg]))
        write_print(sd / f'strain_piston_{st}.dat', None, np.column_stack([steps, 140 + 0 * steps, 140 - eps * L0 * drive, eps * L0 * drive / 140]))
        write_print(vd / f'gel_dimensions_bb_{st}.dat', None, np.column_stack([steps, 40 + 0 * steps, 40 + 0 * steps, (L0 + 2) * (1 - eps * drive)]))
        write_print(vd / f'gel_dimensions_rg_{st}.dat', None, np.column_stack([steps, 40 + 0 * steps, 40 + 0 * steps, Lrg]))
        write_print(vd / f'gel_volume_rg_{st}.dat', None, np.column_stack([steps, 1600 * Lrg]))
        write_print(vd / f'gel_volume_bb_{st}.dat', None, np.column_stack([1600 * (Lrg + 2)]))
        write_print(vd / f'box_dimensions_{st}.dat', None, np.column_stack([steps, LX + 0 * steps, LY + 0 * steps, LZ + 0 * steps]))
        write_print(vd / f'polymer_com_{st}.dat', None, np.column_stack([steps, 20 + 0 * steps, 20 + 0 * steps, 105 + 0 * steps]))
        write_print(vd / f'gel_edges_{st}.dat', None, np.column_stack([steps, Z_GEL_LO + 0 * steps, Z_GEL_HI - eps * L0 * drive, zd, Z_SUP + 0 * steps]))
        pro = stress_profiles(steps, eps=eps)
        for comp in ('zz', 'xx', 'yy'):
            write_ave_time_vector(sd / f'sigma{comp}_polymer_{st}.dat', 'avg', f'c_sigma{comp}_poly', steps, pro[comp][0])
            write_ave_time_vector(sd / f'sigma{comp}_solvent_{st}.dat', 'avg', f'c_sigma{comp}_solv', steps, pro[comp][1])
        write_ave_chunk(cp / f'solvent_density_z_{st}.dat', 'avg_solv_density', steps, density_cols(steps, eps), ['density/number', 'density/mass'])
        disp = []
        for k in range(n):
            u = np.where(gel_mask(), -eps * L0 * (ZC - Z_GEL_LO) / L0 * (1 - np.exp(-(k + 1) / 4.0)), 0.0)
            disp.append(np.column_stack([np.where(gel_mask(), 300, 0), u]))
        write_ave_chunk(dd / f'disp_z_polymer_{st}.dat', 'avg_disp_z_poly', steps, disp, ['v_uz_poly'])
        write_traj(r / 'traj_files' / f'traj_stress_{st}.lammpstrj', steps[:2], {4: Z_SUP, 5: Z_FEED, 6: Z_PERM, 7: Z_LOAD - eps * L0})
        for dim in ('x', 'y', 'z'):
            comps = ['polymer', 'solvent'] + (['piston', 'piston_feed', 'piston_perm', 'support'] if dim == 'z' else [])
            nb = int((LX if dim == 'x' else LY if dim == 'y' else LZ) / BW)
            for comp in comps:
                prof = [np.abs(rng.normal(0.5, 0.05, nb)) for _ in steps]
                write_ave_time_vector(sd / f'stress_{dim}_{comp}_{st}.dat', 'avg', 'c_press', steps, prof)
        t_start += hold + 1000
    with open(r / f'final_tricomp_{rst}.data', 'w') as f:
        f.write(f'fixture\n\n10 atoms\n7 atom types\n\n0.0 {LX} xlo xhi\n0.0 {LY} ylo yhi\n0.0 {LZ} zlo zhi\n')
    (r / 'log.lammps').write_text(_log(np.arange(0, 20000, 1000)))
    return rst


def flatten_to_flow_data(run_root, dest, traj_dest):
    """Copy every .dat under output_files/ flat into dest (the notebooks' layout) and
    the trajectories into traj_dest."""
    dest, traj_dest = Path(dest), Path(traj_dest)
    dest.mkdir(parents=True, exist_ok=True)
    traj_dest.mkdir(parents=True, exist_ok=True)
    import shutil
    for p in Path(run_root, 'output_files').rglob('*.dat'):
        shutil.copy(p, dest / p.name)
    for p in Path(run_root, 'traj_files').glob('*.lammpstrj'):
        shutil.copy(p, traj_dest / p.name)


if __name__ == '__main__':
    root = Path(sys.argv[1] if len(sys.argv) > 1 else 'fixtures')
    s1 = one_piston(root / 'one_piston')
    s2 = two_piston_perm(root / 'two_piston_perm')
    s3 = two_piston_comp(root / 'two_piston_comp')
    fdl = root / 'flow_data_local'
    flatten_to_flow_data(root / 'two_piston_perm', fdl / 'permeation' / 'fixture_perm', fdl / 'traj_files.nosync')
    flatten_to_flow_data(root / 'two_piston_comp', fdl / 'compression' / 'fixture_comp', fdl / 'traj_files.nosync')
    print(f'fixtures written under {root}\n  one_piston stem      {s1}\n  two_piston_perm stem {s2}\n  two_piston_comp stem {s3}')
