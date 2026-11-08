#!/bin/bash
# run_lammps.sh -- run one LAMMPS job on Expanse. Every .batch in simulations/ calls it.
#
#   run_lammps.sh <folder> <dataname> <interaction> <nsteps>
#
#   folder       simulations/<folder>/<folder>.lmp is the deck (DECK_FILE=<name> picks <name>.lmp instead)
#   dataname     <dataname>.data in $INPUT_DATA_DIR (lib/site.sh; LAMMPS_INPUT_DATA_OVERRIDE for another folder)
#   interaction  epsSS_epsSP, e.g. 1.0_1.0
#   nsteps       production steps (also the tag of every output file)
#
# Deck knobs come from the environment, named like their LAMMPS variable in upper case
# (PRESS_TARGET -> -var press_target); lib/run_common.sh lists them. Other environment knobs:
#   CHAIN=1       submit a continuation when a deck with -var wall_timeout reports STATUS INCOMPLETE
#   WALL_MARGIN   seconds before the SLURM limit at which that deck stops (default in lib/site.sh)
#   MPI_MAP_BY    core (default) | node
#   DRY_RUN=1     print the resolved run and the mpirun command, create nothing
#   LAMMPS_RUNS_OVERRIDE   put the run dir flat into this folder instead of $RUNS_DIR/<folder>/
# Legacy positional arguments 5-7 (<type> <press_target> <vel_seed>) are still accepted.
#
# What it does: a timestamped run dir under $RUNS_DIR/<folder>/ with the deck copied into it
# (LAMMPS reads its input lazily, so the repo file must not be the one it runs), the data file
# symlinked, trajectories on scratch; mpirun under the guards of lib/run_common.sh; then the
# continuation (CHAIN=1) and postprocess.sh. Exit status: 0 when LAMMPS and the plots completed.
#
# The { } group makes bash parse the whole file before running it: this script lives on NFS and
# a git pull from another node can replace it while a job runs, which would otherwise end the
# script silently at the next read. Keep the trailing exit inside the group.
{
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/site.sh"
source "$SCRIPT_DIR/lib/run_common.sh"
source "$SCRIPT_DIR/node_git_sync.sh"

if [ $# -lt 4 ]; then
    sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 1
fi
FOLDER=$1; DATANAME=$2; INTERACTION=$3; NSTEPS=$4
[ -n "${6:-}" ] && export PRESS_TARGET="$6"
[ -n "${7:-}" ] && export VEL_SEED="$7"
IFS='_' read -r EPSSS EPSSP <<< "$INTERACTION"

SIM_DIR="$LAMMPS_WORK_DIR/simulations/$FOLDER"
DECK_SRC="$SIM_DIR/${DECK_FILE:-$FOLDER}.lmp"
DATA_SRC="$INPUT_DATA_DIR/$DATANAME.data"
[ -f "$DECK_SRC" ] || { echo "Error: deck not found: $DECK_SRC"; exit 1; }
[ -f "$DATA_SRC" ] || { echo "Error: data file not found: $DATA_SRC"; exit 1; }

[ -n "${LAMMPS_BATCH_SYNCED:-}" ] || sync_lammps_work || true
: "${OMP_NUM_THREADS:=${SLURM_CPUS_PER_TASK:-1}}"; export OMP_NUM_THREADS

WALL_TIMEOUT=$(wall_timeout_seconds)
LAMMPS_VARS=(-var dataname "$DATANAME" -var interaction "$INTERACTION" -var epsSS "$EPSSS" -var epsSP "$EPSSP"
             -var nsteps "$NSTEPS" -var oldsteps 0 -var totsteps "$NSTEPS" -var wall_timeout "$WALL_TIMEOUT")
knob_vars

STAMP=$(date +%Y%m%d_%H%M%S)
RUN_NAME="${FOLDER}_${DATANAME}_${INTERACTION}_${STAMP}"
WORK_DIR="${LAMMPS_RUNS_OVERRIDE:-$RUNS_DIR/$FOLDER}/$RUN_NAME"
TRAJ_DIR="$SCRATCH_DIR/lammps_trajectories/$RUN_NAME"

echo "======================================"
echo "Working directory: $WORK_DIR"
echo "Trajectory directory (scratch): $TRAJ_DIR"
echo "Deck: $DECK_SRC"
echo "Data: $DATA_SRC"
echo "Steps: $NSTEPS   wall budget: ${WALL_TIMEOUT}s (WALL_MARGIN $WALL_MARGIN, CHAIN ${CHAIN:-0})"
echo "Allocation: ${SLURM_NNODES:-?} node(s) x ${SLURM_NTASKS_PER_NODE:-?} tasks x $OMP_NUM_THREADS thread(s), map-by $MPI_MAP_BY"
echo "Vars: ${LAMMPS_VARS[*]}"
echo "======================================"
if [ "${DRY_RUN:-0}" = 1 ]; then
    mpirun_cmd; echo "DRY RUN: ${MPIRUN_CMD[*]} ${LAMMPS_VARS[*]} -in $(basename "$DECK_SRC")"; exit 0
fi

make_run_dirs "$WORK_DIR" "$TRAJ_DIR"
cp "$DECK_SRC" "$WORK_DIR/"
ln -s "$DATA_SRC" "$WORK_DIR/data_files/$DATANAME.data"
cd "$WORK_DIR" || exit 1

launch_lammps "$WORK_DIR/$(basename "$DECK_SRC")" "${LAMMPS_VARS[@]}"
case "$LAMMPS_STATUS" in
    failed)  exit "$LAMMPS_RC" ;;
    partial) export_knobs
             bash "$SCRIPT_DIR/postprocess.sh" "$WORK_DIR" "$FOLDER" "$DATANAME" "$INTERACTION" "$NSTEPS" 0 \
                 || echo "WARNING: post-processing of the partial run failed."
             echo ">>> PARTIAL RUN: exiting $LAMMPS_RC (the simulation did not finish)."
             exit "$LAMMPS_RC" ;;
esac

chain_if_incomplete "${SLURM_SUBMIT_DIR:-$SIM_DIR}"
export_knobs
bash "$SCRIPT_DIR/postprocess.sh" "$WORK_DIR" "$FOLDER" "$DATANAME" "$INTERACTION" "$NSTEPS" 0
exit $?
}
