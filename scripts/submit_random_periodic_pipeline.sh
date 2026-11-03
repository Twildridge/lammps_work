#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# submit_random_periodic_pipeline.sh (2026-11-02) -- ONE-SHOT submission of the random
# periodic slab pipeline.  Run it once on the login node and log out:
#
#   bash scripts/submit_random_periodic_pipeline.sh
#
#   1. builds ~/Documents/lammps_data/input_data/slab_support_periodic_5beads_tall_rho04_random.data
#      if it is not there yet (scripts/lib/periodic_data_file_generator.py --network random:
#      the reference 9 x 9 x 23 periodic slab with its lattice re-wired into a random f = 4
#      network; ~10 min on the login node with the anaconda3 module).  A file from before
#      2026-11-03 (old junction relaxation, interior voids; no "relax rest_len" tag in its
#      header) is moved to <name>.old_relax.data and rebuilt;
#   2. sbatch simulations/slab_with_support/slab_with_support_periodic_random.batch
#      (14 M steps, 4 nodes, 24 h, CHAIN=1);
#   3. sbatch simulations/triaxial_compression_two_pist/chain_random_periodic.batch with
#      --dependency=afterany on the slab job: when the slab run ends it copies the final
#      configuration to input_data, converts it to the two-piston geometry and submits
#      triaxial_compression_two_pist_random.batch (eps 0.10, 9 M-step hold).
#
# Check on it later with:   squeue -u $USER
#   slab job output:        simulations/slab_with_support/slab_periodic_random.o<job>.*
#   bridge output:          simulations/triaxial_compression_two_pist/chain_random_periodic.o<job>.*
#   compression output:     simulations/triaxial_compression_two_pist/triax_comp2p_random.o<job>.*
# ---------------------------------------------------------------------------
set -u
LAMMPS_WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_DATA="$HOME/Documents/lammps_data/input_data"
DATANAME="slab_support_periodic_5beads_tall_rho04_random"
SLAB_DIR="$LAMMPS_WORK_DIR/simulations/slab_with_support"
TRIAX_DIR="$LAMMPS_WORK_DIR/simulations/triaxial_compression_two_pist"
SLAB_BATCH="slab_with_support_periodic_random.batch"
BRIDGE_BATCH="chain_random_periodic.batch"

for f in "$SLAB_DIR/$SLAB_BATCH" "$TRIAX_DIR/$BRIDGE_BATCH" "$TRIAX_DIR/triaxial_compression_two_pist_random.batch" \
         "$LAMMPS_WORK_DIR/scripts/lib/periodic_data_file_generator.py"; do
    [ -f "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

# ---- 1. the data file -----------------------------------------------------------------
mkdir -p "$INPUT_DATA"
# A random file built before 2026-11-03 (no "relax rest_len=..." tag in its header line) has
# the junction clumps / interior voids of the old relaxation: set it aside and rebuild.
if [ -s "$INPUT_DATA/$DATANAME.data" ] && ! head -1 "$INPUT_DATA/$DATANAME.data" | grep -q "relax rest_len"; then
    echo ">>> $INPUT_DATA/$DATANAME.data was built with the old junction relaxation (interior voids): moving it to $DATANAME.old_relax.data and rebuilding"
    mv "$INPUT_DATA/$DATANAME.data" "$INPUT_DATA/$DATANAME.old_relax.data"
fi
if [ -s "$INPUT_DATA/$DATANAME.data" ]; then
    echo ">>> data file present: $INPUT_DATA/$DATANAME.data"
else
    echo ">>> building $INPUT_DATA/$DATANAME.data (random periodic slab; this takes a few minutes)..."
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module unload python/3.8.12 2>/dev/null || true
    module load anaconda3/2021.05/q4munrg || { echo "ERROR: anaconda3 module not available"; exit 1; }
    LOG="$INPUT_DATA/$DATANAME.build.log"
    if ! python "$LAMMPS_WORK_DIR/scripts/lib/periodic_data_file_generator.py" --network random \
            --output "$INPUT_DATA/$DATANAME.data" 2>&1 | tee "$LOG"; then
        echo "ERROR: the build failed (see $LOG); nothing submitted."
        rm -f "$INPUT_DATA/$DATANAME.data"
        exit 1
    fi
    grep -q "ISOTROPY OK" "$LOG" || echo "WARNING: the isotropy check did not pass (see $LOG) -- submitting anyway."
    grep -q "DENSITY OK" "$LOG" || echo "WARNING: the density check did not pass (see $LOG) -- submitting anyway."
fi
ls -la "$INPUT_DATA/$DATANAME.data"

# ---- 2. the slab equilibration -------------------------------------------------------
cd "$SLAB_DIR" || exit 1
SLAB_JOB=$(sbatch --parsable "$SLAB_BATCH") || { echo "ERROR: sbatch of $SLAB_BATCH failed"; exit 1; }
SLAB_JOB="${SLAB_JOB%%;*}"
echo ">>> slab job: $SLAB_JOB   ($SLAB_DIR/$SLAB_BATCH)"

# ---- 3. the bridge (waits for the slab, converts, submits the compression) -----------
cd "$TRIAX_DIR" || exit 1
BRIDGE_JOB=$(sbatch --parsable --dependency=afterany:"$SLAB_JOB" --export=ALL,SLAB_JOB="$SLAB_JOB",CHAIN_HOPS=0 "$BRIDGE_BATCH") \
    || { echo "ERROR: sbatch of $BRIDGE_BATCH failed; the slab job $SLAB_JOB is running alone -- submit the bridge by hand:"; \
         echo "       cd $TRIAX_DIR && sbatch --dependency=afterany:$SLAB_JOB --export=ALL,SLAB_JOB=$SLAB_JOB $BRIDGE_BATCH"; exit 1; }
BRIDGE_JOB="${BRIDGE_JOB%%;*}"
echo ">>> bridge job: $BRIDGE_JOB  (afterany:$SLAB_JOB; it submits triaxial_compression_two_pist_random.batch)"
echo
squeue -u "$USER" -o "%.10i %.22j %.4t %.10M %.6D %R"
echo
echo "Done. Nothing else to do: check squeue -u $USER later."
