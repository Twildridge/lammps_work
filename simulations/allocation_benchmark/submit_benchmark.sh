#!/usr/bin/env bash
# submit_benchmark.sh -- submit one tier of the allocation benchmark (2026-11-05), modelled on
# ../chain_length_testing/submit_chain_length_sweep.sh.  Run from the login node of Expanse or Pod:
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
#   nodescan    (2026-11-07, after the debug layout scan) the node-count question at 128 x 1 on the
#               compute partition, 10 k steps like the debug scan: --map-by core + cutoff 2.5 on
#               1/2/3/4/6/8 nodes three times each and once on 16; core + cutoff 5.0 (Slurm-only
#               answer) once on 1/3/4/6/8; node + 2.5 on 2/4 and 64 x 2 core + 2.5 on 2/4, once
#               (SCAN_NODES="1 2 3 4 6 8" and BIG_NODES="16" set the node counts)
#   layouts     (2026-11-18) the ranks x threads question without a debug partition: every layout
#               of the cluster on LAYOUT_NODES (default "1 2 4") with core + 2.5, plus the full
#               x 1 layout with node + 5.0 (the settings of the old runners) as the reference
#
# CLUSTER (2026-11-18).  Every "128 x 1" above is "all the cores of a node x 1" and every layout
# list is the cluster's own (scripts/lib/site.sh: CORES_PER_NODE, BENCH_LAYOUTS, BENCH_SBATCH,
# BENCH_TIME_MULT).  Expanse: 128 cores, 128x1 64x2 32x4 16x8.  Pod: 40 cores, 40x1 20x2 10x4 8x5,
# default partition, no account, time limits x 4.  "SU" is core-hours: charged on Expanse, free
# on Pod, where SU_CAP only bounds the queue load.  The assumed Pod speeds are guesses: correct
# them with RATE_MULT after the pilot.  subnode needs Expanse's shared partition.
#
# NETWORK=rot45 (default random) selects the slab family: prepared states, tags and inputs then
# carry a _<network> suffix, as in debug_sweep.batch / run_benchmark.sh.  NS="5 15" limits the N
# of nodescan and layouts; with NS set, prep prepares only those N at scale 1.0.
# BENCH_STEM=<stem> (2026-11-23) benchmarks one data file outside that family, input_data/<stem>.data:
# NS must then name a single N (its beads per strand) and NETWORK labels the gel in the tags, the
# prepared state and results.csv (plot with --network <label>), as in
#   NETWORK=random_wide NS=10 BENCH_STEM=slab_support_periodic_10beads_wide_rho04_random bash submit_benchmark.sh prep
#   NETWORK=random_cube NS=10 BENCH_STEM=slab_support_10beads_cube10_rho044_nowalls_random_br2_s11 bash submit_benchmark.sh prep
#
# A run whose log already exists, or that is already queued, is skipped, so a tier can be re-submitted
# after a failure.  A tier that would take the total past SU_CAP is refused.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ../../scripts/lib/site.sh
mkdir -p logs_slurm
C="$CORES_PER_NODE"                          # cores of one node: the "full x 1" layout is C x 1
PART="$BENCH_PARTITION"
read -r -a LAYOUTS <<<"$BENCH_LAYOUTS"       # ranks x threads per node, full x 1 first
HYBRID="${LAYOUTS[1]}"                       # the first hybrid split (64x2 on Expanse, 20x2 on Pod)

SU_CAP="${SU_CAP:-4000}"
RATE_MULT="${RATE_MULT:-1.0}"      # scale the assumed speeds below after the pilot
TIER="${1:?usage: submit_benchmark.sh <tier> [--dry-run]; see the top of this file}"
shift
DRY=0
ARGS=()
for a in "$@"; do if [ "$a" = "--dry-run" ]; then DRY=1; else ARGS+=("$a"); fi; done

ALL_N=(5 10 15 20 25 30)
NETWORK="${NETWORK:-random}"
NET_SUFFIX=""; [ "$NETWORK" = random ] || NET_SUFFIX="_${NETWORK}"
read -r -a SCAN_N <<<"${NS:-5 15}"
BENCH_STEM="${BENCH_STEM:-}"
if [ -n "$BENCH_STEM" ] && [ "${#SCAN_N[@]}" -ne 1 ]; then
    echo "BENCH_STEM=$BENCH_STEM needs exactly one N in NS (NS='${NS:-5 15}'); nothing submitted."; exit 1
fi
ENDS=(5 30)
THREE=(5 15 30)
NODE_COUNTS=(1 2 4 8)
read -r -a SCAN_NODES <<<"${SCAN_NODES:-1 2 3 4 6 8}"
read -r -a LAYOUT_NODES <<<"${LAYOUT_NODES:-1 2 4}"
if [ "$CLUSTER" = expanse ]; then BIG_NODES="${BIG_NODES-16}"; else BIG_NODES="${BIG_NODES-}"; fi
LOG_DIR="$RUNS_DIR/allocation_benchmark/logs"
INPUT_DIR="$INPUT_DATA_DIR"
QUEUED=$(squeue -u "$USER" -h -o %j 2>/dev/null || true)
PLANNED=0
SKIPPED=0
EST_SU=0

# Assumed steps/s at scale 1.  Expanse: anchored on ~310 steps/s for 4 x 128 (chain_length_testing_sweep.batch);
# the other values are guesses until the pilot replaces them through RATE_MULT.  Pod: guesses throughout
# (one 40-core node at roughly a third of an Expanse node), nothing measured with these settings yet.
if [ "$CLUSTER" = pod ]; then
    rate() { case "$1" in 1) echo 60;; 2) echo 100;; 3) echo 130;; 4) echo 150;; *) echo 180;; esac; }
else
    rate() { case "$1" in 32c) echo 80;; 64c) echo 130;; 1) echo 200;; 2) echo 260;; 3) echo 290;; 4) echo 310;; 6|8|16|32) echo 350;; esac; }
fi
hms() { printf '%02d:%02d:00' $(( $1 / 3600 )) $(( $1 % 3600 / 60 )); }

# add_run ROLE N SCALE PARTITION NODES TASKS CPUS MAP CUT REP STEPS
add_run() {
    local role=$1 n=$2 scale=$3 part=$4 nodes=$5 tasks=$6 cpus=$7 map=$8 cut=$9 rep=${10} steps=${11}
    local tag="${role}_N${n}_s${scale}${NET_SUFFIX}_${part}_n${nodes}_t${tasks}_c${cpus}_${map}_cc${cut}_r${rep}"
    if grep -qs '^complete=1$' "$LOG_DIR/$tag.meta" || grep -qx "$tag" <<<"$QUEUED"; then
        SKIPPED=$((SKIPPED + 1)); return 0
    fi
    # sbatch options of the cluster (Expanse: partition, account, all of the node's memory); on
    # Expanse's shared partition 1 GB per core, where the charge follows max(cores, GB/2)
    local sb="${BENCH_SBATCH[*]-}" key=$nodes cores=$((nodes * C)) limit_s=1200 startup=60 limit
    if [ "$part" = "shared" ]; then sb="--partition=shared --account=csb197 --mem-per-cpu=1G"; key="${tasks}c"; cores=$((tasks * cpus)); fi
    if [ "$role" = "prep" ]; then limit_s=5400; startup=600; fi
    if [ "$scale" = "2.0" ]; then limit_s=2400; [ "$role" = "prep" ] && limit_s=10800; fi
    limit=$(hms $((limit_s * BENCH_TIME_MULT)))
    local est
    est=$(awk -v s="$steps" -v r="$(rate "$key")" -v m="$RATE_MULT" -v sc="$scale" -v c="$cores" -v u="$startup" \
          'BEGIN { printf "%.1f", (s * sc / (r * m) + u) * c / 3600 }')
    EST_SU=$(awk -v a="$EST_SU" -v b="$est" 'BEGIN { printf "%.1f", a + b }')
    PLANNED=$((PLANNED + 1))
    if [ "$DRY" = "1" ]; then
        printf '%-68s %-6s %8s SU\n' "$tag" "$limit" "$est"
    else
        PENDING+=("$tag|$role|$n|$scale|$part|$nodes|$tasks|$cpus|$map|$cut|$rep|$steps|$limit|$sb")
    fi
}

PENDING=()
prod()  { add_run bench "$1" "${3:-1.0}" "$PART" "$2" "$C" 1 node 5.0 1 20000; }
tuned() { add_run bench "$1" "${5:-1.0}" "$PART" "$2" "$3" "$4" core 2.5 1 20000; }

case "$TIER" in
    prep)
        if [ -n "${NS:-}" ]; then
            for n in "${SCAN_N[@]}"; do add_run prep "$n" 1.0 "$PART" 1 "$C" 1 core 7.0 1 50000; done
        else
            for n in "${ALL_N[@]}"; do add_run prep "$n" 1.0 "$PART" 1 "$C" 1 core 7.0 1 50000; done
            for s in 0.5 2.0; do for n in "${THREE[@]}"; do add_run prep "$n" "$s" "$PART" 1 "$C" 1 core 7.0 1 50000; done; done
        fi
        ;;
    pilot)
        for nodes in 1 4; do prod "${SCAN_N[0]}" "$nodes"; tuned "${SCAN_N[0]}" "$nodes" "$C" 1; done
        ;;
    production)
        for n in "${ALL_N[@]}"; do for nodes in "${NODE_COUNTS[@]}"; do prod "$n" "$nodes"; done; done
        ;;
    standard)
        for nodes in "${NODE_COUNTS[@]}"; do for map in node core; do for lay in "${LAYOUTS[@]}"; do
            add_run bench 5 1.0 "$PART" "$nodes" "${lay%x*}" "${lay#*x}" "$map" 5.0 1 20000
        done; done; done
        ;;
    one)
        [ "${#ARGS[@]}" -eq 7 ] || { echo "usage: submit_benchmark.sh one N SCALE NODES TASKS CPUS MAP CUT"; exit 1; }
        add_run bench "${ARGS[0]}" "${ARGS[1]}" "$PART" "${ARGS[2]}" "${ARGS[3]}" "${ARGS[4]}" "${ARGS[5]}" "${ARGS[6]}" 1 20000
        ;;
    tuned)
        for n in "${ALL_N[@]}"; do for nodes in "${NODE_COUNTS[@]}"; do for lay in "${LAYOUTS[@]}"; do
            tuned "$n" "$nodes" "${lay%x*}" "${lay#*x}"
        done; done; done
        ;;
    factors)
        for n in "${ENDS[@]}"; do for nodes in 2 4 8; do
            add_run bench "$n" 1.0 "$PART" "$nodes" "$C" 1 core 5.0 1 20000
            add_run bench "$n" 1.0 "$PART" "$nodes" "$C" 1 node 2.5 1 20000
        done; done
        ;;
    subnode)
        [ "$CLUSTER" = expanse ] || { echo "subnode needs Expanse's shared partition; not available on $CLUSTER."; exit 1; }
        for n in "${THREE[@]}"; do for tasks in 32 64; do
            add_run bench "$n" 1.0 shared 1 "$tasks" 1 core 2.5 1 20000
        done; done
        ;;
    scale)
        for s in 0.5 2.0; do for n in "${THREE[@]}"; do for nodes in "${NODE_COUNTS[@]}"; do
            tuned "$n" "$nodes" "$C" 1 "$s"
        done; done; done
        ;;
    nodescan)
        for n in "${SCAN_N[@]}"; do
            for nodes in "${SCAN_NODES[@]}"; do for r in 1 2 3; do
                add_run bench "$n" 1.0 "$PART" "$nodes" "$C" 1 core 2.5 "$r" 10000
            done; done
            for nodes in $BIG_NODES; do add_run bench "$n" 1.0 "$PART" "$nodes" "$C" 1 core 2.5 1 10000; done
            for nodes in "${SCAN_NODES[@]}"; do
                [ "$nodes" = 2 ] || add_run bench "$n" 1.0 "$PART" "$nodes" "$C" 1 core 5.0 1 10000
            done
            for nodes in 2 4; do
                add_run bench "$n" 1.0 "$PART" "$nodes" "$C" 1 node 2.5 1 10000
                add_run bench "$n" 1.0 "$PART" "$nodes" "${HYBRID%x*}" "${HYBRID#*x}" core 2.5 1 10000
            done
        done
        ;;
    layouts)
        for n in "${SCAN_N[@]}"; do for nodes in "${LAYOUT_NODES[@]}"; do
            for lay in "${LAYOUTS[@]}"; do
                add_run bench "$n" 1.0 "$PART" "$nodes" "${lay%x*}" "${lay#*x}" core 2.5 1 10000
            done
            add_run bench "$n" 1.0 "$PART" "$nodes" "$C" 1 node 5.0 1 10000
        done; done
        ;;
    repeat)
        [ "${#ARGS[@]}" -ge 6 ] || { echo "usage: submit_benchmark.sh repeat N NODES TASKS CPUS MAP CUT [REPS]"; exit 1; }
        for ((r = 2; r <= ${ARGS[6]:-3}; r++)); do
            add_run bench "${ARGS[0]}" 1.0 "$PART" "${ARGS[1]}" "${ARGS[2]}" "${ARGS[3]}" "${ARGS[4]}" "${ARGS[5]}" "$r" 20000
        done
        ;;
    *) echo "unknown tier: $TIER"; exit 1 ;;
esac

SPENT=0
if [ -d "$LOG_DIR" ]; then SPENT=$(python3 collect_results.py --logs "$LOG_DIR" --su-total 2>/dev/null || echo 0); fi
echo "cluster $CLUSTER ($C cores per node; layouts ${LAYOUTS[*]})"
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
    stem="${BENCH_STEM:-slab_support_periodic_${n}beads_tall_rho04_${NETWORK}}"
    if [ "$role" = "prep" ]; then
        if [ "$scale" = "1.0" ]; then f="$INPUT_DIR/$stem.data"; else f="$INPUT_DIR/allocation_benchmark/scale$scale/$stem.data"; fi
        hint="python build_benchmark_data.py"; [ -z "$BENCH_STEM" ] || hint="copy it from the machine that built it"
    else
        f="$INPUT_DIR/allocation_benchmark/prepared/prep_${n}b_s${scale}${NET_SUFFIX}.data"
        hint="bash submit_benchmark.sh prep, and wait for it"
    fi
    [ -s "$f" ] || { echo "missing input file: $f  ($hint); nothing submitted."; exit 1; }
done

for item in "${PENDING[@]:-}"; do
    [ -n "$item" ] || continue
    IFS='|' read -r tag role n scale part nodes tasks cpus map cut rep steps limit sb <<<"$item"
    # shellcheck disable=SC2086  # $sb is the cluster's list of sbatch flags (may be empty)
    jid=$(sbatch --parsable --job-name="$tag" $sb --nodes="$nodes" --ntasks-per-node="$tasks" \
        --cpus-per-task="$cpus" --time="$limit" \
        --export=ALL,CLUSTER="$CLUSTER",ROLE="$role",NETWORK="$NETWORK",NBEADS="$n",SCALE="$scale",MAP="$map",COMM_CUT="$cut",NSTEPS="$steps",PREP_STEPS="$steps",REP="$rep",TAG="$tag",BENCH_STEM="$BENCH_STEM" \
        allocation_benchmark.batch)
    echo "$tag: job $jid"
done
