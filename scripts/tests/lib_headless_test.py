"""Drive lib/uniaxial.py and lib/bulk.py on the synthetic flow_data_local tree exactly as the
two-piston and the bulk-modulus notebooks do (SYNC off).
Usage: python lib_headless_test.py <fixture_root>"""
import sys, importlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.show = lambda *a, **k: None
root = Path(sys.argv[1])
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
import uniaxial as uni
uni.setup_style()
base = str(root / 'flow_data_local')

print('\n################ compression single / sweep (two-piston fixture) ################')
cfg = uni.Config(DATANAME='fixture_slab', INTERACTION='1.0_1.0', RUN_ID='fixture_comp', COMP_LEVELS=['0.05', '0.10'],
                 mode='compression', base_dir=base, VOR_ENABLE=False)
R = uni.load_reference(cfg)
LEVELS = [L for L in (uni.load_level(cfg, R, l) for l in cfg.COMP_LEVELS) if L is not None]
assert len(LEVELS) == 2, 'levels not loaded'
uni.add_volume_fractions(cfg, R, LEVELS)
uni.print_summary(cfg, LEVELS)
L = LEVELS[1]
for f in (uni.fig_strain, uni.fig_volfrac, uni.fig_total_stress, uni.fig_partial_stress, uni.fig_network_stress,
          uni.fig_piston, uni.fig_M, uni.fig_ratio, uni.fig_G, uni.fig_Dc_reswelling, uni.fig_kappa, uni.fig_wet_pistons, uni.fig_solvent_expelled,
          uni.fig_thermo_pressure, uni.fig_osmotic_pressure):
    f(cfg, R, L) if f is not uni.fig_strain else f(cfg, R, [L])
    plt.close('all')
for f in (uni.fig_volfrac_sweep, uni.fig_total_stress_sweep, uni.fig_partial_stress_sweep, uni.fig_network_stress_sweep,
          uni.fig_piston_sweep, uni.fig_M_sweep, uni.fig_stress_strain_sweep, uni.fig_ratio_sweep, uni.fig_G_sweep,
          uni.fig_Dc_sweep, uni.fig_kappa_sweep, uni.fig_wet_pistons_sweep, uni.fig_thermo_pressure_sweep, uni.fig_osmotic_pressure_sweep):
    f(cfg, R, LEVELS)
    plt.close('all')
uni.fig_strain(cfg, R, LEVELS, stem='sweep_strain_diagnostic')
assert L['wet'] is not None and 'P_feed_meas' in L['wet']['plat']

print('\n################ permeation (two-piston fixture) ################')
cfgp = uni.Config(DATANAME='fixture_slab', INTERACTION='1.0_1.0', RUN_ID='fixture_perm', mode='permeation', base_dir=base,
                  P_CAL_MODE='pore')          # the fixture traj has 80 mobile atoms: the Voronoi pass is cheap here
Rp = uni.load_reference(cfgp)
P = uni.load_permeation(cfgp, Rp)
uni.add_perm_volume_fractions(cfgp, Rp, P)
uni.add_perm_psd(cfgp, Rp, P)                # geometric porosity + PSD on the same (tiny) frames
uni.add_perm_displacement(cfgp, Rp, P)
uni.add_perm_vctrl(cfgp, Rp, P)              # deck Phase 2v: feed piston at a prescribed velocity, live network
uni.add_perm_frozen(cfgp, Rp, P)             # deck Phase 2d: the frozen network at the same velocity
uni.print_perm_summary(cfgp, Rp, P)
X = uni.perm_reconcile(cfgp, Rp, P)
# the fixture's stages: one prescribed flux, dP 0.10 live and 0.12 frozen
Sv, Sf = P['vctrl'], P['frozen']
assert Sv is not None and Sv['mode'] == 'velocity' and Sv['primary'] == 'v_measured', 'velocity-controlled stage not read'
assert Sf is not None and Sf['mode'] == 'velocity' and Sf['ratio_ref'].startswith('velocity')
assert abs(Sv['dP_meas'] - 0.10) < 0.03 and abs(Sf['dP_meas'] - 0.12) < 0.03, (Sv['dP_meas'], Sf['dP_meas'])   # 12 samples of 0.019 noise each
assert abs(Sv['Q_ratio'] - 1.0) < 0.05, f"N-slope flux / prescribed flux = {Sv['Q_ratio']}"
assert abs(Sf['ratio']['v_measured']['r'] - Sv['dP_meas'] / Sf['dP_meas']) < 1e-3, Sf['ratio']   # same flux, so the kappa ratio is the inverse dP ratio (L differs by noise only)
assert abs(Sv['k']['v_measured']['k'] - Sv['q_set'] * Sv['L'] / Sv['dP_meas']) < 1e-9
assert {'v', 'f', 'dP'} <= set(X['kappa']) and 'kappa_f over kappa_v' in X['closure'], X
for f in (uni.fig_perm_stages, uni.fig_perm_reconcile, uni.fig_perm_frozen):
    assert f(cfgp, Rp, P) is not None, f.__name__
    plt.close('all')
# a run without the stage files (every run before 2026-11-16) still loads and reconciles on the constant-dP kappa
cfg0 = uni.Config(DATANAME='fixture_slab', INTERACTION='1.0_1.0', RUN_ID='fixture_perm', mode='permeation', base_dir=base,
                  VOR_ENABLE=False, PSD_ENABLE=False, PERM_VCTRL=False, PERM_FROZEN=False)
P0 = uni.load_permeation(cfg0, Rp, verbose=False)
uni.add_perm_displacement(cfg0, Rp, P0, verbose=False)
uni.add_perm_vctrl(cfg0, Rp, P0, verbose=False); uni.add_perm_frozen(cfg0, Rp, P0, verbose=False)
assert P0['vctrl'] is None and P0['frozen'] is None
X0 = uni.perm_reconcile(cfg0, Rp, P0, verbose=False)
assert 'v' not in X0['kappa'] and 'dP' in X0['kappa']
uni.fig_perm_stages(cfg0, Rp, P0); uni.fig_perm_reconcile(cfg0, Rp, P0)
plt.close('all')
for f in (uni.fig_perm_pistons, uni.fig_total_stress, uni.fig_partial_stress, uni.fig_network_stress, uni.fig_perm_density,
          uni.fig_perm_volfrac, uni.fig_perm_volfrac_evolution, uni.fig_perm_psd, uni.fig_perm_flux, uni.fig_perm_permeability,
          uni.fig_thermo_pressure, uni.fig_osmotic_pressure):
    f(cfgp, Rp, P)
    plt.close('all')
assert P['flux'] is not None and any(k.endswith('applied') for k in P['flux']['k'])
assert 'v_applied' not in P['flux']['k'], 'legacy block-mean permeability must not be reported'
assert P['phi_mf'] is not None and P['phi_vor'] is not None
assert P['psd'] is not None and Rp['psd'] is not None, 'PSD not computed on the fixture frames'
por = P['psd']['por'][0]
assert np.nanmax(por) <= 1.0 + 1e-9 and np.nanmin(por) >= 0.0
assert np.isfinite(por[P['interior']]).all(), 'porosity missing inside the membrane'
assert np.all(np.isnan(por[Rp['z'] < P['z_mem_lo'] - 3 * cfgp.binWidth])), 'grid must not cover the far reservoir'
dm = P['psd']['d_mean'][0]
assert np.nanmin(dm[P['interior']]) >= 2 * cfgp.PSD_R_PROBE - 1e-9, 'pore diameter below the probe diameter'
Dc, m, lo, hi = P['psd']['regions']['interior']
assert abs(np.nansum(m) * cfgp.PSD_DBIN - 1.0) < 1e-6, 'PSD not normalised'
if Rp.get('CALIB') is not None:
    assert P['phi_cal'] is not None and P['P_cal'] is not None
    # the pore-pressure ramp: feed baseline above the membrane, permeate baseline below it
    Pl = P['P_cal'][0]
    assert Pl[P['z'] > P['z_mem_hi']].max() >= Pl[P['z'] < P['z_mem_lo']].min() - 1e-9
    # the other two modes must run too (cached tessellation -> instant)
    for mode in ('thermo', 'const'):
        cfgp.P_CAL_MODE = mode
        P2 = uni.load_permeation(cfgp, Rp, verbose=False)
        uni.add_perm_volume_fractions(cfgp, Rp, P2)
        assert P2['phi_cal'] is not None, mode
        plt.close('all')

print('\n################ bulk modulus single / sweep (six-plate fixture) ################')
import bulk
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_fixtures import K_TRUE, B_PI_REF, B_P_REF, B_DC, B_DT, B_DM
cfgb = bulk.Config(DATANAME='fixture_cube', INTERACTION='1.0_1.0', RUN_ID='fixture_bulk', COMP_LEVELS=['0.05', '0.10'],
                   base_dir=base, M_REF=0.30, G_REF=0.04, dt_lj=B_DT, VOR_ENABLE=False)
Rb = bulk.load_reference(cfgb)
LB = [L for L in (bulk.load_level(cfgb, Rb, l) for l in cfgb.COMP_LEVELS) if L is not None]
assert len(LB) == 2, 'bulk levels not loaded'
assert set(Rb['ax']) == {'x', 'y', 'z'}, 'a profile axis is missing'
bulk.add_volume_fractions(cfgb, Rb, LB)        # VOR_ENABLE False -> prints the skip note only
bulk.print_summary(cfgb, LB)
assert abs(Rb['Pi_ref'] - B_PI_REF) < 0.003 and abs(Rb['P_ref'] - B_P_REF) < 0.003, 'zero-strain readings not recovered'
for L in LB:
    assert abs(L['strain']['plate'] - float(L['lvl'])) < 1e-3 * float(L['lvl']) + 1e-6, 'plate strain not recovered'
    assert abs(L['strain']['plate'] - L['strain']['plate_geo']) < 1e-4, "deck and geometric plate strains disagree"
    assert L['eps_src'] == 'disp' and 'disp' in L['strain'], 'displacement-profile strain (disp_*_polymer_cum) not used'
    assert abs(L['strain']['disp'] - float(L['lvl'])) < 2e-3, f"disp-profile eps_vol {L['strain']['disp']:.4f} vs built {L['lvl']}"
    for a in 'xyz':
        assert L['disp_fit'][a]['R2'] > 0.99, f'cumulative displacement profile along {a} not linear'
    for key in ('K_net', 'K_pl'):
        assert abs(L[key] - K_TRUE) < 0.06 * K_TRUE, f"{key} = {L[key]:.4f} at eps_vol {L['lvl']}, built from {K_TRUE}"
        assert L[key + '_lo'] <= L[key] <= L[key + '_hi']
    for a in cfgb.AXES:
        assert abs(L['ax'][a]['K'] - K_TRUE) < 0.10 * K_TRUE, f"K along {a} = {L['ax'][a]['K']:.4f}"
    assert L['dev_rel'] < 0.2, 'isotropic synthetic stress reads as anisotropic'
    assert L['wet'] is not None and abs(L['wet']['plat']['P_feed_meas'] - 1.5) < 0.02
    assert 'Dc' not in L and 'Dc_axis' not in L, 'the compression hold must not produce a D_c any more (2026-10-28)'
    assert L['U'] is not None and L['Dm'] is not None, 'the unload (_u files) was not loaded / fitted'
    F = L['Dm']
    assert abs(F['Dm'] - B_DM) < 0.15 * B_DM and F['edges']['R2'] > 0.95, f"D_m (Rg edges) = {F['Dm']:.4f} (R^2 {F['edges']['R2']:.3f}), built from {B_DM}"
    assert F['Dm_lo'] <= F['Dm'] <= F['Dm_hi'] and not F['edges']['at_bound'], 'D_m interval / bound'
    assert abs(F['volume']['Dm'] - B_DM) < 0.15 * B_DM, f"D_m (Rg volume) = {F['volume']['Dm']:.4f}"
    assert abs(F['stress']['Dm'] - B_DM) < 0.25 * B_DM, f"D_m (interior Pi) = {F['stress']['Dm']:.4f}"
    assert F.get('prof_joint') is not None and not F['prof_joint']['at_bound'], 'profile D_m missing / at a bound'
    assert abs(F['prof_joint']['Dm'] - B_DM) < 0.2 * B_DM, f"D_m (profiles) = {F['prof_joint']['Dm']:.4f}"
    assert np.all(np.abs(F['edges']['recovery'] - 1.0) < 0.1), f"recovery {F['edges']['recovery']}"
    assert abs(L['kappa'] - F['Dm'] / L['K_net']) < 1e-12 and abs(L['Dc_from_Dm'] - F['Dm'] * 0.30 / L['K_net']) < 1e-12
    assert abs(L['G_from_MK'] - 0.75 * (0.30 - L['K_net'])) < 1e-12
Lb = LB[1]
bulk.fig_strain(cfgb, Rb, [Lb])
for f in (bulk.fig_volfrac, bulk.fig_total_stress, bulk.fig_partial_stress, bulk.fig_network_stress, bulk.fig_thermo_pressure,
          bulk.fig_osmotic_pressure):
    f(cfgb, Rb, Lb, axes=cfgb.AXES)
    plt.close('all')
for f in (bulk.fig_plates, bulk.fig_K, bulk.fig_isotropy, bulk.fig_Dm, bulk.fig_closure, bulk.fig_wet_pistons, bulk.fig_ratio):
    assert f(cfgb, Rb, Lb) is not None, f.__name__
    plt.close('all')
bulk.fig_strain(cfgb, Rb, LB, stem='sweep_strain_diagnostic')
assert bulk.fig_anisotropy_run(cfgb, Rb, LB) is not None, 'stress_aniso whole-run figure'
assert bulk.fig_volfrac_evolution(cfgb, Rb, Lb, key='phi_mf') is not None
assert bulk.fig_volfrac_evolution(cfgb, Rb, Lb, key='phi_cal') is None      # VOR_ENABLE False here: skipped, not an error
assert bulk.fig_ratio_sweep(cfgb, Rb, LB) is not None
plt.close('all')
for f in (bulk.fig_volfrac_sweep, bulk.fig_total_stress_sweep, bulk.fig_partial_stress_sweep, bulk.fig_network_stress_sweep,
          bulk.fig_thermo_pressure_sweep, bulk.fig_osmotic_pressure_sweep):
    f(cfgb, Rb, LB, axes=cfgb.AXES)
    plt.close('all')
for f in (bulk.fig_plates_sweep, bulk.fig_K_sweep, bulk.fig_stress_strain_sweep, bulk.fig_isotropy_sweep, bulk.fig_Dm_sweep,
          bulk.fig_closure_sweep, bulk.fig_wet_pistons_sweep):
    assert f(cfgb, Rb, LB) is not None, f.__name__
    plt.close('all')
assert abs(Rb['K_small_net'] - K_TRUE) < 0.06 * K_TRUE and abs(Rb['K_small_pl'] - K_TRUE) < 0.06 * K_TRUE, 'small-strain slopes'
pngs = {q.name for q in cfgb.PLOT_DIR.glob('*.png')}
assert any('_alongx_' in q for q in pngs) and any('_alongy_' in q for q in pngs), 'per-axis figures were not tagged'
assert uni._save.__module__ == 'uniaxial', 'the axis context did not restore uniaxial._save'
# K_STRAIN = 'rg': the same stress over the network's own strain
cfgb.K_STRAIN = 'rg'
L2 = bulk.load_level(cfgb, Rb, '0.10', verbose=False)
assert abs(L2['eps'] - 0.095) < 0.002 and abs(L2['K_net'] - K_TRUE * 0.10 / 0.095) < 0.06 * K_TRUE
print('\nLIB HEADLESS TEST: OK')
