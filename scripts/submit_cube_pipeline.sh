#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit_cube_pipeline.sh (2026-11-11; random 10-bead cube 2026-11-22) -- ONE-SHOT cube pipeline
# on Expanse: build, equilibrate, plate, compress and re-swell. Run it once on a login node and log out:
#
#   bash scripts/submit_cube_pipeline.sh [NETWORK] [N] [BEADS]      # defaults: random 10 10
#   bash scripts/submit_cube_pipeline.sh rot45 12 5                 # the 2026-11 lattice cube
#   BUILD_ONLY=1 bash scripts/submit_cube_pipeline.sh               # dry run: stage 1 only
#
#   NETWORK  random | rot45 | tetrahedral   (a preset in data_file_generator.NETWORKS)
#   N        cubic cells per edge of the gel
#   BEADS    beads per strand (5 = the lattice cubes; 10 = the gel of the 2026-11 uniaxial runs)
#
# The default is the cube twin of the wide random 10-bead slab (random_16M_comp_1): 10 cells of
# a_eq = 7.78 sigma, about 78 sigma per edge as built, about 1.2 M atoms.
#
#   1. builds input_data/slab_support_<BEADS>beads_cube<N>_<rho tag>_nowalls_<NETWORK>.data on the login
#      node (scripts/lib/data_file_generator.build, no sheets). rot45 / random count BCT cells of
#      a/sqrt2 along x and z, so units_x = units_z = round(N sqrt2), units_y = N; tetrahedral uses N
#      for all three. 5 beads: pre_swell 1.10, straight strands. Longer strands: the lattice constant
#      and pre_swell of build_beads_sweep.knobs (a_eq = 4.82 ((BEADS-1)/4)^0.59) and strands laid as
#      random walks between their junctions (chain_placement 'walk'), as the periodic 10-bead slabs.
#      Random networks with longer strands (2026-11-22) also take BRIDGE (default 2) and SEED (default 11),
#      which go into the file name (_br<BRIDGE>_s<SEED>): every face junction with BRIDGE or fewer strands
#      gets one extra strand to a neighbour with a free site (data_file_generator.bridge_dangling: fewer
#      dangling ends, so a thinner surface fringe, at the cost of stiffer face junctions; BRIDGE=0 = none),
#      and SEED is the bond-switching seed (11 was the most isotropic of 7, 11, 13, 17 on the 10-cell cube).
#      The build takes about 10 min on a free machine and far longer on a busy login node; a file built
#      elsewhere and copied to input_data/ under the same name is used as it is.
#      An existing file is never overwritten: the build is skipped;
#   2. sbatch simulations/slab_with_support/slab_with_support_single.batch (DATANAME = that stem; iso
#      barostat, PIN_GEL, CHAIN=1). 5 beads: 8 M steps. Longer strands: a 16 M-step budget with the
#      Rg-plateau halt (RG_HALT, 1 M-step blocks), the settings of the wide 10-bead slab;
#   3+4. sbatch simulations/isotropic_compression/chain_cube_six_plates.batch with --dependency=afterany
#      on the slab job: it copies the final config to input_data/, plates it with
#      lib/add_more_plates_to_gel.py --layout six_plates and submits isotropic_compression_single.batch.
#      That batch's defaults are the protocol: plates seated on the Rg faces (SEAT_MODE=rg), eps_vol
#      0.10 hold -> K, then the free re-swelling (UNLOAD=1) with displacement profiles every 25 tau
#      (DISP_FINE_NFREQ=5000) -> D_m from the joint profile fit (isotropic_compression_single.ipynb,
#      figures 9 and 9b). The deck stops before the unload if the network took less than half the
#      applied strain (RG_STRAIN_MIN).
#
# Check on it later with:   squeue -u $USER
#   slab job output:        simulations/slab_with_support/slab_isolated_tall_random.o<job>.*
#   bridge output:          simulations/isotropic_compression/chain_cube_six_plates.o<job>.*
#   compression output:     simulations/isotropic_compression/isotropic_compression.o<job>.*
# ---------------------------------------------------------------------------
set -u
NETWORK="${1:-random}"
N="${2:-10}"
BEADS="${3:-10}"
BUILD_ONLY="${BUILD_ONLY:-0}"
[[ "$N" =~ ^[0-9]+$ ]] && [ "$N" -gt 0 ] || { echo "ERROR: N must be a positive integer (got '$N')"; exit 1; }
[[ "$BEADS" =~ ^[0-9]+$ ]] && [ "$BEADS" -ge 3 ] || { echo "ERROR: BEADS must be an integer >= 3 (got '$BEADS')"; exit 1; }
# Initial bath (2026-11-12). The cube12 rot45 run of 2026-11-12 was built at solvent density 1.0
# with 15-sigma pads: at P* = 1.5 that bath is far too dense (P = 77 at step 0), so the iso
# barostat dilated the box -- and the gel with it -- 2.3x in volume within 20 k steps, and the
# gel then needed ~2 M steps to contract back to its built size (bath 0.47 next to the gel,
# box-average solvent 0.44 at equilibrium). 0.42 in this larger box is the equilibrium solvent
# count of the 5-bead cube (0.47 in the bath + ~0.23 inside the gel), so the box and gel start
# where they end up. The pads are sized for the six_plates crop of add_more_plates_to_gel.py: it
# wants box >= L_bb + 2 (3 + 12) in x, y and L_bb + 2 (3 + 15) in z (a smaller box is converted
# with a thinner bath and a NOTE). L_bb is the lattice edge plus the outer-bead excursions: +7-8
# sigma on the 5-bead cube, and 13-16 sigma on the random 10-bead gels, whose surface is a diffuse
# fringe (the `>>> unload opening` line of the compression decks).
# Longer strands (2026-11-22): the gel holds more solvent (polymer fraction about 0.4 against the
# 5-bead cube's 0.6+), so the bath starts at 0.44; pads 26 / 28 for the wider fringe; the
# equilibration gets the wide 10-bead slab's budget and Rg-plateau halt (that slab's Rg was flat
# after 14.6 M steps at 78 sigma thickness). All three are estimates: the iso barostat sets the
# real box, and the bridge prints the converter's actual margins.
if [ "$BEADS" -le 5 ]; then
    RHO_SOLV=0.42; RHO_TAG=rho042; PAD_XY=21.0; PAD_Z=20.0
    SLAB_NSTEPS=8000000; SLAB_TIME=""; SLAB_ENV=()
else
    RHO_SOLV=0.44; RHO_TAG=rho044; PAD_XY=26.0; PAD_Z=28.0
    SLAB_NSTEPS=16000000; SLAB_TIME="--time=30:00:00"
    SLAB_ENV=(RG_HALT=1 RG_BLOCK=1000000 RG_EVERY=10000 RG_TOL=0.0015 RG_NCONSEC=3 RG_EXTRA=0.10)
fi
# Bridging and seed of the random long-strand builds (see the header); no effect, and no name tag, otherwise.
BRIDGE="${BRIDGE:-2}"
SEED="${SEED:-11}"
NET_TAG=""
if [ "$NETWORK" = random ] && [ "$BEADS" -gt 5 ]; then NET_TAG="_br${BRIDGE}_s${SEED}"; else BRIDGE=0; SEED=7; fi

LAMMPS_WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$LAMMPS_WORK_DIR/scripts/lib/site.sh"
INPUT_DATA="$INPUT_DATA_DIR"
SLAB_DIR="$LAMMPS_WORK_DIR/simulations/slab_with_support"
COMP_DIR="$LAMMPS_WORK_DIR/simulations/isotropic_compression"
SLAB_BATCH="slab_with_support_single.batch"
BRIDGE_BATCH="chain_cube_six_plates.batch"
GEN="$LAMMPS_WORK_DIR/scripts/lib/data_file_generator.py"

case "$NETWORK" in
    rot45|random) UX=$(python3 -c "print(round($N * 2 ** 0.5))") || { echo "ERROR: could not compute units_x"; exit 1; }
                  UY="$N"; UZ="$UX" ;;
    tetrahedral)  UX="$N"; UY="$N"; UZ="$N" ;;
    *) echo "ERROR: NETWORK must be rot45, tetrahedral or random (got '$NETWORK')"; exit 1 ;;
esac

STEM="slab_support_${BEADS}beads_cube${N}_${RHO_TAG}_nowalls_${NETWORK}${NET_TAG}"
DATA="$INPUT_DATA/$STEM.data"
echo ">>> NETWORK=$NETWORK  N=$N  BEADS=$BEADS  ->  units ${UX} x ${UY} x ${UZ}   STEM=$STEM"
echo ">>> bath $RHO_SOLV, pads $PAD_XY (x, y) / $PAD_Z (z);  equilibration budget $SLAB_NSTEPS steps ${SLAB_ENV[*]:-}"

for f in "$GEN" "$LAMMPS_WORK_DIR/scripts/lib/add_more_plates_to_gel.py" \
         "$SLAB_DIR/$SLAB_BATCH" "$COMP_DIR/$BRIDGE_BATCH" "$COMP_DIR/isotropic_compression_single.batch"; do
    [ -f "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

# ---- 1. the data file (login node) ---------------------------------------------------
mkdir -p "$INPUT_DATA"
if [ -e "$DATA" ]; then
    echo ">>> $DATA already exists; skipping the build (never overwrite a file a run may have used)."
else
    echo ">>> building $DATA (this takes a few minutes)..."
    load_python_module || { echo "ERROR: could not load $PYTHON_MODULE"; exit 1; }
    LOG="$INPUT_DATA/$STEM.build.log"
    cd "$LAMMPS_WORK_DIR/scripts/lib" || exit 1
    python -u - "$NETWORK" "$DATA" "$UX" "$UY" "$UZ" "$RHO_SOLV" "$PAD_XY" "$PAD_Z" "$BEADS" "$BRIDGE" "$SEED" <<'EOF' 2>&1 | tee "$LOG"
import sys
import data_file_generator as gen
network, out, ux, uy, uz = sys.argv[1], sys.argv[2], *map(int, sys.argv[3:6])
rho, pad_xy, pad_z = map(float, sys.argv[6:9])
beads, bridge, seed = int(sys.argv[9]), int(sys.argv[10]), int(sys.argv[11])
extra = dict(bond_switch_seed=seed, bridge_max_f=bridge) if (network == 'random' and beads > 5) else {}
if beads <= 5:                      # the lattice cubes: a = 4.8 x 1.10 = 5.28 sigma, straight strands
    pre_swell, placement = 1.10, 'straight'
else:                               # longer strands: the equilibrium mesh estimate of the beads-per-chain sweep
    import build_beads_sweep as bbs
    k = bbs.knobs(beads, network='random' if network == 'random' else 'rot45')
    pre_swell, placement = k['pre_swell'], 'walk'
    print(f"cell estimate: a_eq {k['a_eq']:.3f} sigma, pre_swell {pre_swell:.4f}")
print(f"CELL_A {1.2 * (beads - 1) * pre_swell:.4f}")
gen.build(network=network, output_file=out, beads_per_chain=beads,
          units_x=ux, units_y=uy, units_z=uz, solvent_density=rho,
          side_thickness=pad_xy, below_thickness=pad_z, top_thickness=pad_z,
          pre_swell=pre_swell, include_walls=False, chain_placement=placement, **extra)
EOF
    if [ "${PIPESTATUS[0]}" -ne 0 ] || [ ! -s "$DATA" ]; then
        echo "ERROR: the build of $DATA failed (see $LOG); nothing submitted."
        rm -f "$DATA"
        exit 1
    fi
    # the three gel edges must agree to within one cell (a = 1.2 (BEADS - 1) pre_swell; 5.28 sigma for 5 beads)
    CELL_A=$(grep -m1 '^CELL_A ' "$LOG" | awk '{print $2}')
    DIMS=$(grep -m1 "Gel dimensions:" "$LOG")
    echo ">>> $DIMS"
    grep -m1 "Box dimensions:" "$LOG"
    if ! echo "$DIMS" | python3 -c "
import sys, re
v = [float(x) for x in re.findall(r'[0-9.]+', sys.stdin.read().split(':', 1)[1])]
sys.exit(0 if len(v) == 3 and max(v) - min(v) <= float('${CELL_A:-5.28}') else 1)"; then
        echo "ERROR: the gel edges in $LOG are not cubic to within one cell ($DIMS); nothing submitted."
        echo "       The file is kept at $DATA -- move it aside before re-running."
        exit 1
    fi
fi
ls -la "$DATA"

if [ "$BUILD_ONLY" = "1" ]; then
    echo ">>> BUILD_ONLY=1: stopping before any submission."
    exit 0
fi

# ---- 2. the equilibration ------------------------------------------------------------
cd "$SLAB_DIR" || exit 1
SLAB_JOB=$(env DATANAME="$STEM" NSTEPS="$SLAB_NSTEPS" ${SLAB_ENV[@]+"${SLAB_ENV[@]}"} sbatch --parsable $SLAB_TIME "$SLAB_BATCH") \
    || { echo "ERROR: sbatch of $SLAB_DIR/$SLAB_BATCH failed; nothing submitted."; exit 1; }
SLAB_JOB="${SLAB_JOB%%;*}"
echo ">>> slab job: $SLAB_JOB   ($SLAB_DIR/$SLAB_BATCH, DATANAME $STEM, NSTEPS $SLAB_NSTEPS ${SLAB_ENV[*]:-} $SLAB_TIME)"

# ---- 3+4. the bridge (waits for the slab, plates, submits isotropic_compression) --------------
cd "$COMP_DIR" || exit 1
BRIDGE_JOB=$(sbatch --parsable --dependency=afterany:"$SLAB_JOB" \
        --export=ALL,SLAB_JOB="$SLAB_JOB",STEM="$STEM",SLAB_NSTEPS="$SLAB_NSTEPS",CHAIN_HOPS=0 "$BRIDGE_BATCH") \
    || { echo "ERROR: sbatch of $BRIDGE_BATCH failed; the slab job $SLAB_JOB is running alone -- submit the bridge by hand:"; \
         echo "       cd $COMP_DIR && sbatch --dependency=afterany:$SLAB_JOB --export=ALL,SLAB_JOB=$SLAB_JOB,STEM=$STEM,SLAB_NSTEPS=$SLAB_NSTEPS $BRIDGE_BATCH"; exit 1; }
BRIDGE_JOB="${BRIDGE_JOB%%;*}"
echo ">>> bridge job: $BRIDGE_JOB  (afterany:$SLAB_JOB; it plates the cube and submits isotropic_compression_single.batch: SEAT_MODE=rg, eps_vol 0.10, unload with fine profiles)"
echo
squeue -u "$USER" -o "%.10i %.22j %.4t %.10M %.6D %R"
echo
echo "Done. Nothing else to do: check squeue -u $USER later."
