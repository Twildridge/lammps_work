#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit_cube_pipeline.sh (2026-11-11) -- ONE-SHOT cube pipeline on Expanse:
# build, equilibrate, plate, compress. Run it once on a login node and log out:
#
#   bash scripts/submit_cube_pipeline.sh [NETWORK] [N]          # defaults: rot45 12
#   BUILD_ONLY=1 bash scripts/submit_cube_pipeline.sh rot45 12  # dry run: stage 1 only
#
#   NETWORK  rot45 | tetrahedral | random  (a preset in data_file_generator.NETWORKS)
#   N        cubic cells per edge of the gel (12 = the 2026-10 cube, ~63 sigma)
#
#   1. builds input_data/slab_support_5beads_cube<N>_rho042_nowalls_<NETWORK>.data on the login
#      node (scripts/lib/data_file_generator.build: 5 beads per chain, no sheets, solvent density
#      0.42, bath pads 21 sigma in x/y and 20 sigma below/above). rot45 / random count BCT cells of a/sqrt2 along x and z, so
#      units_x = units_z = round(N sqrt2), units_y = N; tetrahedral uses N for all three.
#      An existing file is never overwritten: the build is skipped;
#   2. sbatch simulations/slab_with_support/slab_with_support_single.batch (DATANAME = that stem,
#      NSTEPS = 8 M; iso barostat, PIN_GEL, CHAIN=1);
#   3+4. sbatch simulations/compress_slab/chain_cube_six_plates.batch with --dependency=afterany
#      on the slab job: it copies the final config to input_data/, plates it with
#      lib/add_more_plates_to_gel.py --layout six_plates and submits compress_slab_single.batch
#      (eps 0.10 hold -> K, then UNLOAD=1 free re-swelling -> D_m).
#
# Check on it later with:   squeue -u $USER
#   slab job output:        simulations/slab_with_support/slab_isolated_tall_random.o<job>.*
#   bridge output:          simulations/compress_slab/chain_cube_six_plates.o<job>.*
#   compression output:     simulations/compress_slab/compress_slab.o<job>.*
# ---------------------------------------------------------------------------
set -u
NETWORK="${1:-rot45}"
N="${2:-12}"
BUILD_ONLY="${BUILD_ONLY:-0}"
SLAB_NSTEPS=8000000
# Initial bath (2026-11-12). The cube12 rot45 run of 2026-11-12 was built at solvent density 1.0
# with 15-sigma pads: at P* = 1.5 that bath is far too dense (P = 77 at step 0), so the iso
# barostat dilated the box -- and the gel with it -- 2.3x in volume within 20 k steps, and the
# gel then needed ~2 M steps to contract back to its built size (bath 0.47 next to the gel,
# box-average solvent 0.44 at equilibrium). 0.42 in this larger box is the equilibrium solvent
# count (0.47 in the bath + ~0.23 inside the gel), so the box and gel start where they end up
# (~40 % fewer atoms too). The pads are sized for the six_plates crop of
# add_more_plates_to_gel.py: it needs box >= L_bb + 2 (3 + 12) in x, y and L_bb + 2 (3 + 15) in z,
# with L_bb ~70-71 sigma for the 12-cell cube (the lattice edge plus outer-bead excursions).
RHO_SOLV=0.42
RHO_TAG=rho042
PAD_XY=21.0
PAD_Z=20.0

LAMMPS_WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$LAMMPS_WORK_DIR/scripts/lib/site.sh"
INPUT_DATA="$INPUT_DATA_DIR"
SLAB_DIR="$LAMMPS_WORK_DIR/simulations/slab_with_support"
COMP_DIR="$LAMMPS_WORK_DIR/simulations/compress_slab"
SLAB_BATCH="slab_with_support_single.batch"
BRIDGE_BATCH="chain_cube_six_plates.batch"
GEN="$LAMMPS_WORK_DIR/scripts/lib/data_file_generator.py"

[[ "$N" =~ ^[0-9]+$ ]] && [ "$N" -gt 0 ] || { echo "ERROR: N must be a positive integer (got '$N')"; exit 1; }
case "$NETWORK" in
    rot45|random) UX=$(python3 -c "print(round($N * 2 ** 0.5))") || { echo "ERROR: could not compute units_x"; exit 1; }
                  UY="$N"; UZ="$UX" ;;
    tetrahedral)  UX="$N"; UY="$N"; UZ="$N" ;;
    *) echo "ERROR: NETWORK must be rot45, tetrahedral or random (got '$NETWORK')"; exit 1 ;;
esac

STEM="slab_support_5beads_cube${N}_${RHO_TAG}_nowalls_${NETWORK}"
DATA="$INPUT_DATA/$STEM.data"
echo ">>> NETWORK=$NETWORK  N=$N  ->  units ${UX} x ${UY} x ${UZ}   STEM=$STEM"

for f in "$GEN" "$LAMMPS_WORK_DIR/scripts/lib/add_more_plates_to_gel.py" \
         "$SLAB_DIR/$SLAB_BATCH" "$COMP_DIR/$BRIDGE_BATCH" "$COMP_DIR/compress_slab_single.batch"; do
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
    python - "$NETWORK" "$DATA" "$UX" "$UY" "$UZ" "$RHO_SOLV" "$PAD_XY" "$PAD_Z" <<'EOF' 2>&1 | tee "$LOG"
import sys
import data_file_generator as gen
network, out, ux, uy, uz = sys.argv[1], sys.argv[2], *map(int, sys.argv[3:6])
rho, pad_xy, pad_z = map(float, sys.argv[6:9])
gen.build(network=network, output_file=out, beads_per_chain=5,
          units_x=ux, units_y=uy, units_z=uz, solvent_density=rho,
          side_thickness=pad_xy, below_thickness=pad_z, top_thickness=pad_z,
          pre_swell=1.10, include_walls=False)
EOF
    if [ "${PIPESTATUS[0]}" -ne 0 ] || [ ! -s "$DATA" ]; then
        echo "ERROR: the build of $DATA failed (see $LOG); nothing submitted."
        rm -f "$DATA"
        exit 1
    fi
    # the three gel edges must agree to within one cell (a = 4.8 x pre_swell = 5.28 sigma)
    DIMS=$(grep -m1 "Gel dimensions:" "$LOG")
    echo ">>> $DIMS"
    grep -m1 "Box dimensions:" "$LOG"
    if ! echo "$DIMS" | python3 -c "
import sys, re
v = [float(x) for x in re.findall(r'[0-9.]+', sys.stdin.read().split(':', 1)[1])]
sys.exit(0 if len(v) == 3 and max(v) - min(v) <= 5.28 else 1)"; then
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
SLAB_JOB=$(DATANAME="$STEM" NSTEPS="$SLAB_NSTEPS" sbatch --parsable "$SLAB_BATCH") \
    || { echo "ERROR: sbatch of $SLAB_DIR/$SLAB_BATCH failed; nothing submitted."; exit 1; }
SLAB_JOB="${SLAB_JOB%%;*}"
echo ">>> slab job: $SLAB_JOB   ($SLAB_DIR/$SLAB_BATCH, DATANAME $STEM, NSTEPS $SLAB_NSTEPS)"

# ---- 3+4. the bridge (waits for the slab, plates, submits compress_slab) --------------
cd "$COMP_DIR" || exit 1
BRIDGE_JOB=$(sbatch --parsable --dependency=afterany:"$SLAB_JOB" \
        --export=ALL,SLAB_JOB="$SLAB_JOB",STEM="$STEM",SLAB_NSTEPS="$SLAB_NSTEPS",CHAIN_HOPS=0 "$BRIDGE_BATCH") \
    || { echo "ERROR: sbatch of $BRIDGE_BATCH failed; the slab job $SLAB_JOB is running alone -- submit the bridge by hand:"; \
         echo "       cd $COMP_DIR && sbatch --dependency=afterany:$SLAB_JOB --export=ALL,SLAB_JOB=$SLAB_JOB,STEM=$STEM,SLAB_NSTEPS=$SLAB_NSTEPS $BRIDGE_BATCH"; exit 1; }
BRIDGE_JOB="${BRIDGE_JOB%%;*}"
echo ">>> bridge job: $BRIDGE_JOB  (afterany:$SLAB_JOB; it plates the cube and submits compress_slab_single.batch)"
echo
squeue -u "$USER" -o "%.10i %.22j %.4t %.10M %.6D %R"
echo
echo "Done. Nothing else to do: check squeue -u $USER later."
