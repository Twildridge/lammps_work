"""Drive lib/triaxial.py and lib/bulk.py on the synthetic flow_data_local tree exactly as the
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
import triaxial as tri
tri.setup_style()
base = str(root / 'flow_data_local')

print('\n################ compression single / sweep (two-piston fixture) ################')
cfg = tri.Config(DATANAME='fixture_slab', INTERACTION='1.0_1.0', RUN_ID='fixture_comp', COMP_LEVELS=['0.05', '0.10'],
                 mode='compression', two_pist=True, base_dir=base, VOR_ENABLE=False)
R = tri.load_reference(cfg)
LEVELS = [L for L in (tri.load_level(cfg, R, l) for l in cfg.COMP_LEVELS) if L is not None]
assert len(LEVELS) == 2, 'levels not loaded'
tri.add_volume_fractions(cfg, R, LEVELS)
tri.print_summary(cfg, LEVELS)
L = LEVELS[1]
for f in (tri.fig_strain, tri.fig_volfrac, tri.fig_total_stress, tri.fig_partial_stress, tri.fig_network_stress,
          tri.fig_piston, tri.fig_M, tri.fig_ratio, tri.fig_G, tri.fig_Dc, tri.fig_kappa, tri.fig_wet_pistons, tri.fig_solvent_expelled,
          tri.fig_thermo_pressure, tri.fig_osmotic_pressure):
    f(cfg, R, L) if f is not tri.fig_strain else f(cfg, R, [L])
    plt.close('all')
for f in (tri.fig_volfrac_sweep, tri.fig_total_stress_sweep, tri.fig_partial_stress_sweep, tri.fig_network_stress_sweep,
          tri.fig_piston_sweep, tri.fig_M_sweep, tri.fig_stress_strain_sweep, tri.fig_ratio_sweep, tri.fig_G_sweep,
          tri.fig_Dc_sweep, tri.fig_kappa_sweep, tri.fig_wet_pistons_sweep, tri.fig_thermo_pressure_sweep, tri.fig_osmotic_pressure_sweep):
    f(cfg, R, LEVELS)
    plt.close('all')
tri.fig_strain(cfg, R, LEVELS, stem='sweep_strain_diagnostic')
assert L['wet'] is not None and 'P_feed_meas' in L['wet']['plat']

print('\n################ permeation (two-piston fixture) ################')
cfgp = tri.Config(DATANAME='fixture_slab', INTERACTION='1.0_1.0', RUN_ID='fixture_perm', mode='permeation', two_pist=True, base_dir=base,
                  P_CAL_MODE='pore')          # the fixture traj has 80 mobile atoms: the Voronoi pass is cheap here
Rp = tri.load_reference(cfgp)
P = tri.load_permeation(cfgp, Rp)
tri.add_perm_volume_fractions(cfgp, Rp, P)
tri.add_perm_psd(cfgp, Rp, P)                # geometric porosity + PSD on the same (tiny) frames
tri.print_perm_summary(cfgp, Rp, P)
for f in (tri.fig_perm_pistons, tri.fig_total_stress, tri.fig_partial_stress, tri.fig_network_stress, tri.fig_perm_density,
          tri.fig_perm_volfrac, tri.fig_perm_volfrac_evolution, tri.fig_perm_psd, tri.fig_perm_flux, tri.fig_perm_permeability,
          tri.fig_thermo_pressure, tri.fig_osmotic_pressure):
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
        P2 = tri.load_permeation(cfgp, Rp, verbose=False)
        tri.add_perm_volume_fractions(cfgp, Rp, P2)
        assert P2['phi_cal'] is not None, mode
        plt.close('all')

print('\n################ bulk modulus single / sweep (six-plate fixture) ################')
import bulk
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_fixtures import K_TRUE, B_PI_REF, B_P_REF, B_DC, B_DT
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
    for key in ('K_net', 'K_pl'):
        assert abs(L[key] - K_TRUE) < 0.06 * K_TRUE, f"{key} = {L[key]:.4f} at eps_vol {L['lvl']}, built from {K_TRUE}"
        assert L[key + '_lo'] <= L[key] <= L[key + '_hi']
    for a in cfgb.AXES:
        assert abs(L['ax'][a]['K'] - K_TRUE) < 0.10 * K_TRUE, f"K along {a} = {L['ax'][a]['K']:.4f}"
    assert L['dev_rel'] < 0.2, 'isotropic synthetic stress reads as anisotropic'
    assert L['wet'] is not None and abs(L['wet']['plat']['P_feed_meas'] - 1.5) < 0.02
    assert abs(L['eps_expelled'] - float(L['lvl'])) < 0.2 * float(L['lvl'])
    assert L['Dc'] is not None and set(L['Dc']['ax']) == {'x', 'y', 'z'}, 'the held-cube D_c fit is missing an axis'
    assert abs(L['Dc']['Dc'] - B_DC) < 0.15 * B_DC and L['Dc']['R2'] > 0.9, f"cube D_c = {L['Dc']['Dc']:.4f} (R^2 {L['Dc']['R2']:.3f}), built from {B_DC}"
    for a, F in L['Dc_axis'].items():
        assert abs(F['Dc'] - B_DC) < 0.15 * B_DC, f"1-D D_c along {a} = {F['Dc']:.4f}"
    assert abs(L['G_from_MK'] - 0.75 * (0.30 - L['K_net'])) < 1e-12
Lb = LB[1]
bulk.fig_strain(cfgb, Rb, [Lb])
for f in (bulk.fig_volfrac, bulk.fig_total_stress, bulk.fig_partial_stress, bulk.fig_network_stress, bulk.fig_thermo_pressure,
          bulk.fig_osmotic_pressure):
    f(cfgb, Rb, Lb, axes=cfgb.AXES)
    plt.close('all')
for f in (bulk.fig_plates, bulk.fig_K, bulk.fig_isotropy, bulk.fig_Dc, bulk.fig_closure, bulk.fig_wet_pistons, bulk.fig_solvent_expelled, bulk.fig_ratio):
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
for f in (bulk.fig_plates_sweep, bulk.fig_K_sweep, bulk.fig_stress_strain_sweep, bulk.fig_isotropy_sweep, bulk.fig_Dc_sweep,
          bulk.fig_closure_sweep, bulk.fig_wet_pistons_sweep):
    assert f(cfgb, Rb, LB) is not None, f.__name__
    plt.close('all')
assert abs(Rb['K_small_net'] - K_TRUE) < 0.06 * K_TRUE and abs(Rb['K_small_pl'] - K_TRUE) < 0.06 * K_TRUE, 'small-strain slopes'
pngs = {q.name for q in cfgb.PLOT_DIR.glob('*.png')}
assert any('_alongx_' in q for q in pngs) and any('_alongy_' in q for q in pngs), 'per-axis figures were not tagged'
assert tri._save.__module__ == 'triaxial', 'the axis context did not restore triaxial._save'
# K_STRAIN = 'rg': the same stress over the network's own strain
cfgb.K_STRAIN = 'rg'
L2 = bulk.load_level(cfgb, Rb, '0.10', verbose=False)
assert abs(L2['eps'] - 0.095) < 0.002 and abs(L2['K_net'] - K_TRUE * 0.10 / 0.095) < 0.06 * K_TRUE
print('\nLIB HEADLESS TEST: OK')
