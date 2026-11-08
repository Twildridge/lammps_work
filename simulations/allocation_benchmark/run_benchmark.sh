#!/usr/bin/env bash
# run_benchmark.sh -- one run of the allocation benchmark (2026-11-05).  Called by
# allocation_benchmark.batch on Expanse, or directly with BENCH_LOCAL=1 on a Mac.
#
#   ROLE=prep   fresh data file -> the deck's full start-up plus PREP_STEPS of production ->
#               prepared/prep_<N>b_s<S>.data, the common starting state of every timed run.
#   ROLE=bench  prepared state  -> cont=1 -> NSTEPS of production -> the "Loop time" is the result.
#
# The deck is simulations/slab_with_support/slab_with_support.lmp itself, copied into the run
# directory.  Only the COPY changes: comm_modify cutoff becomes a variable (both roles), and for
# ROLE=bench thermo reports every THERMO_EVERY steps at 10 significant figures (accuracy check,
# collect_results.py).
# If the deck ever stops matching those patterns the run aborts rather than time the wrong thing.
#
# Input (environment): ROLE NBEADS SCALE MAP(node|core) COMM_CUT NSTEPS REP TAG
#   NETWORK      random (default) | rot45: the slab family (2026-11-07; prepared states and
#                tags of a non-random network carry a _<network> suffix)
#   BENCH_TASKS / BENCH_CPUS  the ranks per node and threads per rank for THIS run, when they
#                differ from the Slurm allocation (debug_sweep.batch cycles layouts inside one
#                full-node allocation); default: the allocation's own values
#   BENCH_FRESH  path of the fresh data file for ROLE=prep when it is not input_data/<stem>.data
{
set -uo pipefail

ROLE="${ROLE:?prep or bench}"
NBEADS="${NBEADS:?}"
SCALE="${SCALE:-1.0}"
MAP="${MAP:-node}"
COMM_CUT="${COMM_CUT:-5.0}"
NSTEPS="${NSTEPS:-20000}"
PREP_STEPS="${PREP_STEPS:-50000}"
THERMO_EVERY="${THERMO_EVERY:-5000}"
REP="${REP:-1}"
TAG="${TAG:?}"
NETWORK="${NETWORK:-random}"
PARTITION="${SLURM_JOB_PARTITION:-local}"
INTERACTION="1.0_1.0"
PRESS_TARGET=1.5
SEED=12345

if [ "${BENCH_LOCAL:-0}" = "1" ]; then
    REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    BENCH_ROOT="${BENCH_ROOT:?set BENCH_ROOT for a local run}"
    INPUT_DIR="$BENCH_ROOT/input_data"
    RUN_ROOT="$BENCH_ROOT/runs"
    LOG_DIR="$BENCH_ROOT/logs"
    LMP="${LMP:-/Applications/LAMMPS-GUI.app/Contents/bin/lmp}"
    NODES=1; TASKS_PER_NODE=1; CPUS="${CPUS:-1}"
    export OMP_NUM_THREADS="$CPUS"
else
    REPO="$HOME/Documents/lammps_work"
    INPUT_DIR="$HOME/Documents/lammps_data/input_data"
    RUN_ROOT="/expanse/lustre/scratch/$USER/temp_project/allocation_benchmark/runs"
    LOG_DIR="$HOME/Documents/lammps_runs/allocation_benchmark/logs"
    # the same binary as run_lammps.sh
    LMP=/home/dpollard/software/lammps/22Jul2025_update3/mpi-omp/gcc/10.2.0/openmpi/4.1.3/lammps-22Jul2025/build/lmp
    NODES="$SLURM_NNODES"; TASKS_PER_NODE="${BENCH_TASKS:-$SLURM_NTASKS_PER_NODE}"; CPUS="${BENCH_CPUS:-${SLURM_CPUS_PER_TASK:-1}}"
    export OMP_NUM_THREADS="$CPUS"
fi
NTASKS=$((NODES * TASKS_PER_NODE))

STEM="slab_support_periodic_${NBEADS}beads_tall_rho04_${NETWORK}"
NET_SUFFIX=""; [ "$NETWORK" = random ] || NET_SUFFIX="_${NETWORK}"
BENCH_DATA="$INPUT_DIR/allocation_benchmark"
PREP_DIR="$BENCH_DATA/prepared"
PREP_NAME="prep_${NBEADS}b_s${SCALE}${NET_SUFFIX}"
if [ -n "${BENCH_FRESH:-}" ]; then FRESH="$BENCH_FRESH"
elif [ "$SCALE" = "1.0" ]; then FRESH="$INPUT_DIR/$STEM.data"; else FRESH="$BENCH_DATA/scale$SCALE/$STEM.data"; fi

JOB="${SLURM_JOB_ID:-local$$}"
WORK="$RUN_ROOT/${TAG}_$JOB"
mkdir -p "$LOG_DIR" "$PREP_DIR" "$WORK"/{data_files,traj_files,output_files/volume_data,output_files/stress_data}
cd "$WORK" || exit 1

DECK_SRC="$REPO/simulations/slab_with_support/slab_with_support.lmp"
DECK="$WORK/deck.lmp"
cp "$DECK_SRC" "$DECK"
perl -i -pe 's/^comm_modify cutoff 5\.0\s*$/comm_modify cutoff \${comm_cutoff}\n/' "$DECK"
grep -q '^comm_modify cutoff ${comm_cutoff}$' "$DECK" \
    || { echo "ERROR: slab_with_support.lmp no longer matches the benchmark's edits; update run_benchmark.sh"; exit 1; }

if [ "$ROLE" = "prep" ]; then
    [ -s "$FRESH" ] || { echo "ERROR: missing $FRESH (run build_benchmark_data.py)"; exit 1; }
    DATANAME="$PREP_NAME"
    ln -s "$FRESH" "data_files/$DATANAME.data"
    CONT=0; RUN_STEPS="$PREP_STEPS"
else
    PREP_FILE="$PREP_DIR/$PREP_NAME.data"
    [ -s "$PREP_FILE" ] || { echo "ERROR: missing $PREP_FILE (run the prep tier first)"; exit 1; }
    DATANAME="$PREP_NAME"
    ln -s "$PREP_FILE" "data_files/$DATANAME.data"
    CONT=1; RUN_STEPS="$NSTEPS"
    perl -i -pe 's/^variable thermo_freq equal 100000\s*$/variable thermo_freq equal '"$THERMO_EVERY"'\n/; s/^thermo \$\{thermo_freq\}\s*$/thermo \${thermo_freq}\nthermo_modify format float %.10g\n/' "$DECK"
    grep -q "^variable thermo_freq equal $THERMO_EVERY$" "$DECK" \
        && grep -q '^thermo_modify format float %.10g$' "$DECK" \
        || { echo "ERROR: slab_with_support.lmp no longer matches the benchmark's edits; update run_benchmark.sh"; exit 1; }
fi

VARS=(-var dataname "$DATANAME" -var interaction "$INTERACTION" -var epsSS 1.0 -var epsSP 1.0
      -var nsteps "$RUN_STEPS" -var oldsteps 0 -var totsteps "$RUN_STEPS"
      -var press_target "$PRESS_TARGET" -var vel_seed "$SEED"
      -var baro_mode aniso -var pin_gel 0 -var piston_transparent 1
      -var wall_timeout 0 -var n_snaps 1 -var cont "$CONT" -var comm_cutoff "$COMM_CUT")

T0=$(date +%s)
if [ "${BENCH_LOCAL:-0}" = "1" ]; then
    "$LMP" -sf omp -pk omp "$CPUS" "${VARS[@]}" -in "$DECK" -log log.lammps -screen none
    RC=$?
else
    # MAP=node is the production layout (run_lammps.sh): ranks dealt round-robin across the nodes.
    # MAP=core fills each node before the next, so neighbouring ranks usually share a node.
    if [ "$MAP" = "core" ]; then MAPSPEC="ppr:${TASKS_PER_NODE}:node:pe=${CPUS}"; else MAPSPEC="node:pe=${CPUS}"; fi
    mpirun -n "$NTASKS" --bind-to core --map-by "$MAPSPEC" --report-bindings \
        "$LMP" -sf omp -pk omp "$CPUS" "${VARS[@]}" -in "$DECK" -log log.lammps -screen none \
        2> bindings.txt &
    PID=$!
    # Same teardown-hang defence as run_lammps.sh, with a short grace: an idle 8-node job costs real SUs.
    while kill -0 "$PID" 2>/dev/null; do
        if grep -q "Total wall time" log.lammps 2>/dev/null; then
            for _ in 1 2 3 4 5 6; do kill -0 "$PID" 2>/dev/null || break; sleep 5; done
            kill -TERM "$PID" 2>/dev/null; sleep 3; kill -KILL "$PID" 2>/dev/null
            break
        fi
        sleep 5
    done
    wait "$PID" 2>/dev/null
    RC=$?
    # the node a rank landed on, for checking the mapping did what was asked
    grep -m 4 "MCW rank" bindings.txt > bindings_head.txt 2>/dev/null || true
fi
T1=$(date +%s)

COMPLETE=0; grep -q "Total wall time" log.lammps 2>/dev/null && COMPLETE=1

cp log.lammps "$LOG_DIR/$TAG.log" 2>/dev/null
cat > "$LOG_DIR/$TAG.meta" <<EOF
tag=$TAG
role=$ROLE
nbeads=$NBEADS
network=$NETWORK
scale=$SCALE
nodes=$NODES
tasks_per_node=$TASKS_PER_NODE
cpus_per_task=$CPUS
partition=$PARTITION
map=$MAP
comm_cut=$COMM_CUT
rep=$REP
job=$JOB
job_wall_s=$((T1 - T0))
complete=$COMPLETE
mpirun_rc=$RC
deck_git=$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo unknown)
EOF

if [ "$ROLE" = "prep" ]; then
    OUT="final_config_${DATANAME}_${INTERACTION}_${RUN_STEPS}.data"
    if [ "$COMPLETE" = "1" ] && [ -s "$OUT" ]; then
        cp "$OUT" "$PREP_DIR/$PREP_NAME.data"
        echo "prepared: $PREP_DIR/$PREP_NAME.data"
    else
        echo "ERROR: prep did not finish; see $WORK/log.lammps"; exit 1
    fi
fi

if [ "$COMPLETE" = "1" ]; then
    cd "$RUN_ROOT" && rm -rf "$WORK"
    echo "done: $TAG"
else
    echo "FAILED: $TAG; run directory kept: $WORK"
    exit 1
fi
exit 0
}
