#!/usr/bin/env bash
# submit_benchmark.sh -- submit one tier of the allocation benchmark (2026-11-05), modelled on
# ../chain_length_testing/submit_chain_length_sweep.sh.  Run from the login node:
#
#     bash submit_benchmark.sh <tier> --dry-run     # list the runs and the SU estimate; submits nothing
#     bash submit_benchmark.sh <tier>
#
# Tiers, in the order to run them (README.md explains why):
#   prep        one job per N and scale: fresh data file -> common starting state (needed by the rest).
#               Cutoff 7.0: fresh N >= 25 files hold bonds up to ~6 sigma until the push-off shortens them.
#   pilot       N = 5 on 1 and 4 nodes with both settings: checks the mapping and corrects the estimates
#   production  all N, 1/2/4/8 nodes, 128 x 1, production settings (map node, cutoff 5): the baseline
#   standard    N = 5, 1/2/4/8 nodes, every tasks x cpus layout at cutoff 5 with map node and map core:
#               the fastest Slurm configuration of the standard gel with no change to the deck
#   tuned       all N, 1/2/4/8 nodes, 128x1 64x2 32x4 16x8, map core, cutoff 2.5: the full scan
#   factors     N = 5, 30, 2/4/8 nodes: each setting changed alone, so each one's gain is known
#   subnode     N = 5, 15, 30: 32 and 64 cores on the shared partition (charged per core)
#   scale       N = 5, 15, 30 at 0.5x and 2x beads, 1/2/4/8 nodes
#   repeat      repeat N NODES TASKS CPUS MAP CUT [REPS]: extra copies of one layout, to measure noise
#   one         one N SCALE NODES TASKS CPUS MAP CUT: one more layout, e.g. 16 nodes if 8 is still
#               the fastest, or the tasks x cpus layouts of a 0.5x or 2x slab
#
# A run whose log already exists, or that is already queued, is skipped, so a tier can be re-submitted
# after a failure.  A tier that would take the total past SU_CAP is refused.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p logs_slurm

SU_CAP="${SU_CAP:-4000}"
RATE_MULT="${RATE_MULT:-1.0}"      # scale the assumed speeds below after the pilot
TIER="${1:?usage: submit_benchmark.sh <tier> [--dry-run]; see the top of this file}"
shift
DRY=0
ARGS=()
for a in "$@"; do if [ "$a" = "--dry-run" ]; then DRY=1; else ARGS+=("$a"); fi; done

ALL_N=(5 10 15 20 25 30)
ENDS=(5 30)
THREE=(5 15 30)
NODE_COUNTS=(1 2 4 8)
LOG_DIR="$HOME/Documents/lammps_runs/allocation_benchmark/logs"
INPUT_DIR="$HOME/Documents/lammps_data/input_data"
QUEUED=$(squeue -u "$USER" -h -o %j 2>/dev/null || true)
PLANNED=0
SKIPPED=0
EST_SU=0

# Assumed steps/s at scale 1: anchored on ~310 steps/s for 4 x 128 (chain_length_testing_sweep.batch);
# the other values are guesses until the pilot replaces them through RATE_MULT.
rate() { case "$1" in 32c) echo 80;; 64c) echo 130;; 1) echo 200;; 2) echo 260;; 4) echo 310;; 8|16|32) echo 350;; esac; }

# add_run ROLE N SCALE PARTITION NODES TASKS CPUS MAP CUT REP STEPS
add_run() {
    local role=$1 n=$2 scale=$3 part=$4 nodes=$5 tasks=$6 cpus=$7 map=$8 cut=$9 rep=${10} steps=${11}
    local tag="${role}_N${n}_s${scale}_${part}_n${nodes}_t${tasks}_c${cpus}_${map}_cc${cut}_r${rep}"
    if grep -qs '^complete=1$' "$LOG_DIR/$tag.meta" || grep -qx "$tag" <<<"$QUEUED"; then
        SKIPPED=$((SKIPPED + 1)); return 0
    fi
    # memory: all of the node on compute; 1 GB per core on shared, where the charge follows max(cores, GB/2)
    local mem=(--mem=0) key=$nodes cores=$((nodes * 128)) limit=00:20:00 startup=60
    if [ "$part" = "shared" ]; then mem=(--mem-per-cpu=1G); key="${tasks}c"; cores=$((tasks * cpus)); fi
    if [ "$role" = "prep" ]; then limit=01:30:00; startup=600; fi
    if [ "$scale" = "2.0" ]; then limit=00:40:00; [ "$role" = "prep" ] && limit=03:00:00; fi
    local est
    est=$(awk -v s="$steps" -v r="$(rate "$key")" -v m="$RATE_MULT" -v sc="$scale" -v c="$cores" -v u="$startup" \
          'BEGIN { printf "%.1f", (s * sc / (r * m) + u) * c / 3600 }')
    EST_SU=$(awk -v a="$EST_SU" -v b="$est" 'BEGIN { printf "%.1f", a + b }')
    PLANNED=$((PLANNED + 1))
    if [ "$DRY" = "1" ]; then
        printf '%-68s %-6s %8s SU\n' "$tag" "$limit" "$est"
    else
        PENDING+=("$tag|$role|$n|$scale|$part|$nodes|$tasks|$cpus|$map|$cut|$rep|$steps|$limit|${mem[*]}")
    fi
}

PENDING=()
prod()  { add_run bench "$1" "${3:-1.0}" compute "$2" 128 1 node 5.0 1 20000; }
tuned() { add_run bench "$1" "${5:-1.0}" compute "$2" "$3" "$4" core 2.5 1 20000; }

case "$TIER" in
    prep)
        for n in "${ALL_N[@]}"; do add_run prep "$n" 1.0 compute 1 128 1 core 7.0 1 50000; done
        for s in 0.5 2.0; do for n in "${THREE[@]}"; do add_run prep "$n" "$s" compute 1 128 1 core 7.0 1 50000; done; done
        ;;
    pilot)
        for nodes in 1 4; do prod 5 "$nodes"; tuned 5 "$nodes" 128 1; done
        ;;
    production)
        for n in "${ALL_N[@]}"; do for nodes in "${NODE_COUNTS[@]}"; do prod "$n" "$nodes"; done; done
        ;;
    standard)
        for nodes in "${NODE_COUNTS[@]}"; do for map in node core; do
            add_run bench 5 1.0 compute "$nodes" 128 1 "$map" 5.0 1 20000
            add_run bench 5 1.0 compute "$nodes" 64 2 "$map" 5.0 1 20000
            add_run bench 5 1.0 compute "$nodes" 32 4 "$map" 5.0 1 20000
            add_run bench 5 1.0 compute "$nodes" 16 8 "$map" 5.0 1 20000
        done; done
        ;;
    one)
        [ "${#ARGS[@]}" -eq 7 ] || { echo "usage: submit_benchmark.sh one N SCALE NODES TASKS CPUS MAP CUT"; exit 1; }
        add_run bench "${ARGS[0]}" "${ARGS[1]}" compute "${ARGS[2]}" "${ARGS[3]}" "${ARGS[4]}" "${ARGS[5]}" "${ARGS[6]}" 1 20000
        ;;
    tuned)
        for n in "${ALL_N[@]}"; do for nodes in "${NODE_COUNTS[@]}"; do
            tuned "$n" "$nodes" 128 1; tuned "$n" "$nodes" 64 2; tuned "$n" "$nodes" 32 4; tuned "$n" "$nodes" 16 8
        done; done
        ;;
    factors)
        for n in "${ENDS[@]}"; do for nodes in 2 4 8; do
            add_run bench "$n" 1.0 compute "$nodes" 128 1 core 5.0 1 20000
            add_run bench "$n" 1.0 compute "$nodes" 128 1 node 2.5 1 20000
        done; done
        ;;
    subnode)
        for n in "${THREE[@]}"; do for tasks in 32 64; do
            add_run bench "$n" 1.0 shared 1 "$tasks" 1 core 2.5 1 20000
        done; done
        ;;
    scale)
        for s in 0.5 2.0; do for n in "${THREE[@]}"; do for nodes in "${NODE_COUNTS[@]}"; do
            tuned "$n" "$nodes" 128 1 "$s"
        done; done; done
        ;;
    repeat)
        [ "${#ARGS[@]}" -ge 6 ] || { echo "usage: submit_benchmark.sh repeat N NODES TASKS CPUS MAP CUT [REPS]"; exit 1; }
        for ((r = 2; r <= ${ARGS[6]:-3}; r++)); do
            add_run bench "${ARGS[0]}" 1.0 compute "${ARGS[1]}" "${ARGS[2]}" "${ARGS[3]}" "${ARGS[4]}" "${ARGS[5]}" "$r" 20000
        done
        ;;
    *) echo "unknown tier: $TIER"; exit 1 ;;
esac

SPENT=0
if [ -d "$LOG_DIR" ]; then SPENT=$(python3 collect_results.py --logs "$LOG_DIR" --su-total 2>/dev/null || echo 0); fi
echo "tier $TIER: $PLANNED run(s) to submit, $SKIPPED already done or queued; about $EST_SU SU (assumed speeds x $RATE_MULT); $SPENT SU spent so far; cap $SU_CAP"
OVER=$(awk -v a="$SPENT" -v b="$EST_SU" -v c="$SU_CAP" 'BEGIN { print (a + b > c) ? 1 : 0 }')
if [ "$DRY" = "1" ]; then
    [ "$OVER" = "0" ] || echo "NOTE: submitting this tier would be refused: it passes SU_CAP."
    exit 0
fi
[ "$OVER" = "0" ] || { echo "REFUSED: this tier would pass SU_CAP=$SU_CAP (raise it on purpose with SU_CAP=<n>)."; exit 1; }

# Check every input before submitting anything, as submit_chain_length_sweep.sh does.
for item in "${PENDING[@]:-}"; do
    [ -n "$item" ] || continue
    IFS='|' read -r tag role n scale _ <<<"$item"
    stem="slab_support_periodic_${n}beads_tall_rho04_random"
    if [ "$role" = "prep" ]; then
        if [ "$scale" = "1.0" ]; then f="$INPUT_DIR/$stem.data"; else f="$INPUT_DIR/allocation_benchmark/scale$scale/$stem.data"; fi
        hint="python build_benchmark_data.py"
    else
        f="$INPUT_DIR/allocation_benchmark/prepared/prep_${n}b_s${scale}.data"
        hint="bash submit_benchmark.sh prep, and wait for it"
    fi
    [ -s "$f" ] || { echo "missing input file: $f  ($hint); nothing submitted."; exit 1; }
done

for item in "${PENDING[@]:-}"; do
    [ -n "$item" ] || continue
    IFS='|' read -r tag role n scale part nodes tasks cpus map cut rep steps limit mem <<<"$item"
    # shellcheck disable=SC2086  # $mem is one flag
    jid=$(sbatch --parsable --job-name="$tag" --partition="$part" --nodes="$nodes" --ntasks-per-node="$tasks" \
        --cpus-per-task="$cpus" $mem --time="$limit" \
        --export=ALL,ROLE="$role",NBEADS="$n",SCALE="$scale",MAP="$map",COMM_CUT="$cut",NSTEPS="$steps",PREP_STEPS="$steps",REP="$rep",TAG="$tag" \
        allocation_benchmark.batch)
    echo "$tag: job $jid"
done
