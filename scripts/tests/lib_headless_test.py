"""Drive lib/triaxial.py on the synthetic flow_data_local tree exactly as the three
two-piston notebooks do (SYNC off).  Usage: python lib_headless_test.py <fixture_root>"""
import sys, importlib
from pathlib import Path
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
tri.print_perm_summary(cfgp, Rp, P)
for f in (tri.fig_perm_pistons, tri.fig_total_stress, tri.fig_partial_stress, tri.fig_network_stress, tri.fig_perm_density,
          tri.fig_perm_volfrac, tri.fig_perm_volfrac_evolution, tri.fig_perm_flux, tri.fig_perm_permeability,
          tri.fig_thermo_pressure, tri.fig_osmotic_pressure):
    f(cfgp, Rp, P)
    plt.close('all')
assert P['flux'] is not None and any(k.endswith('applied') for k in P['flux']['k'])
assert P['phi_mf'] is not None and P['phi_vor'] is not None
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
print('\nLIB HEADLESS TEST: OK')
