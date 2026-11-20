#!/usr/bin/env python3
"""perm_plan.py -- the permeation drive's settings from the compression run of the SAME data file.

    python3 perm_plan.py <DATANAME>            print the plan
    python3 perm_plan.py <DATANAME> --shell    PLAN_* assignments for the batch file to eval

Reads simulations/uniaxial_permeation/compression_records/<DATANAME>.json, which
uniaxial_compression_single.ipynb writes (uni.save_M_record) whenever the level has an unload D_c.
Standard library only: it runs on the cluster's system python.

Rules (edit the constants below, nothing else holds them):
  dP        = 2 EPS_FEED M           feed-face network strain eps_F = dP/(2M); M = network M of the record
  NSTEPS    = N_TAU1 tau_1 / dt      tau_1 = L0^2/(pi^2 D_c), D_c = unload (re-swelling) profile fit;
                                     rounded up to STEP_ROUND
  2v / 2d   = STAGE_STEPS each, cut equally (not below STAGE_MIN) until the feed piston's travel
              q dt (ramp/2 + NSTEPS + 2v + 2d), q = dP D_c/(M L0), times FLUX_SAFETY fits in the
              feed reservoir above the deck's halt thickness.
"""
import json
import math
import os
import sys
from pathlib import Path

EPS_FEED = 0.12          # target network strain at the feed face
N_TAU1 = 4.1             # drive length in tau_1
STEP_ROUND = 100000      # NSTEPS and the stage lengths are multiples of this
STAGE_STEPS = 4000000    # Phase 2v and Phase 2d, each
STAGE_MIN = 1500000      # shortest stage the plan will cut to before it gives up
RAMP_STEPS = 1000000     # the deck's n_ramp x ramp_steps (100 x 10 k)
FEED_HALT = 2.0          # the deck's feed_halt_thick (sigma)
FLUX_SAFETY = 1.10       # the live flux is known to about 10 % before the run
DT = 0.005

RECORD_DIR = Path(__file__).resolve().parents[2] / 'simulations' / 'uniaxial_permeation' / 'compression_records'


def record_path(dataname):
    return RECORD_DIR / f'{dataname}.json'


def plan(rec):
    """dict of the planned settings from a compression record; raises ValueError when a field is missing."""
    need = ('M_small_net', 'Dc_unload', 'L0', 'feed_thickness')
    missing = [k for k in need if not isinstance(rec.get(k), (int, float)) or not rec[k] > 0]
    if missing:
        raise ValueError('record lacks ' + ', '.join(missing))
    M, Dc, L0, feed = (float(rec[k]) for k in need)
    dt = float(rec.get('dt', DT))
    dP = float(f'{2 * EPS_FEED * M:.2g}')
    tau1 = L0 ** 2 / (math.pi ** 2 * Dc)
    nsteps = int(math.ceil(N_TAU1 * tau1 / dt / STEP_ROUND) * STEP_ROUND)
    q = dP * Dc / (M * L0)
    room = feed - FEED_HALT

    def travel(stage):
        return FLUX_SAFETY * q * dt * (0.5 * RAMP_STEPS + nsteps + 2 * stage)

    stage = STAGE_STEPS
    while stage > STAGE_MIN and travel(stage) > room:
        stage -= STEP_ROUND
    return dict(dP=dP, eps_F=dP / (2 * M), M=M, Dc=Dc, L0=L0, tau1=tau1, nsteps=nsteps, n_tau1=nsteps * dt / tau1,
                q=q, stage=stage, travel=travel(stage), room=room, fits=travel(stage) <= room,
                run_id=rec.get('run_id', '?'), saved=rec.get('saved', '?'))


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    dataname, shell = argv[1], '--shell' in argv[2:]
    p = record_path(dataname)
    if not p.exists():
        print(f'perm_plan: no compression record for this data file ({p.name} not in {RECORD_DIR})', file=sys.stderr)
        return 2
    try:
        P = plan(json.loads(p.read_text()))
    except ValueError as e:
        print(f'perm_plan: {p}: {e}', file=sys.stderr)
        return 3
    out = sys.stderr if shell else sys.stdout
    print(f">>> permeation plan from compression run {P['run_id']} (record saved {P['saved']}):\n"
          f"    M = {P['M']:.4f}, unload D_c = {P['Dc']:.4f}, L0 = {P['L0']:.2f}  ->  tau_1 = {P['tau1']:.0f} tau\n"
          f"    DP_PISTON = {P['dP']:g}  (eps_F = dP/(2M) = {P['eps_F']:.3f})\n"
          f"    NSTEPS = {P['nsteps']}  ({P['n_tau1']:.2f} tau_1)\n"
          f"    VCTRL_STEPS = FROZEN_STEPS = {P['stage']}"
          + ('' if P['stage'] == STAGE_STEPS else f"  (cut from {STAGE_STEPS} for the feed reservoir)") + '\n'
          f"    feed piston travel {P['travel']:.1f} sigma (q = {P['q']:.2e} sigma/tau, x{FLUX_SAFETY:g}) of {P['room']:.1f} available",
          file=out)
    if not P['fits']:
        print(f"perm_plan: the feed reservoir is too short even with {STAGE_MIN}-step stages; rebuild the data file "
              f"with a larger --margin-feed or set the lengths by hand", file=sys.stderr)
        return 4
    if shell:
        print(f"PLAN_DP={P['dP']:g} PLAN_NSTEPS={P['nsteps']} PLAN_STAGE={P['stage']}")
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
