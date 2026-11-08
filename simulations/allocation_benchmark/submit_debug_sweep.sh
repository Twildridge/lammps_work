#!/usr/bin/env bash
# submit_debug_sweep.sh -- the MPI x OpenMP layout scan of a slab family on the DEBUG partition
# (2026-11-07), for the rot45 standard gel (5 and 15 beads per strand) by default.
#
#     bash submit_debug_sweep.sh --dry-run        # print the plan, submit nothing
#     bash submit_debug_sweep.sh                  # submit
#     NS="5" bash submit_debug_sweep.sh           # one N only
#     COMPUTE=0 bash submit_debug_sweep.sh        # debug jobs only, no 4- and 8-node follow-ups
#
# Debug limits: 2 nodes, 30 min, 2 jobs per user (debug-normal QOS).  So, per N:
#   job 1  1 node   prep (start-up + 20 k steps, if the prepared state is missing), then the
#                   layouts 128x1 64x2 32x4 16x8 at the production cutoff 5.0, plus 128x1 at
#                   cutoff 2.5 -- all in one allocation (debug_sweep.batch)
#   job 2  2 nodes  the same four layouts with --map-by node AND --map-by core, plus the 2.5 cutoff
#                   -- queued in chains/debug_queue.txt and submitted by chains/debug_waiter.sh (a
#                   nohup process on the login node) whenever fewer than two debug jobs are queued.
#                   A job cannot submit its successor itself: MaxSubmitJobsPerAccount counts the
#                   running job too (the first attempt, 2026-11-07, was rejected that way).
#   compute follow-ups (COMPUTE=1, default): 4 nodes 128x1 node cutoff 5.0 (the production
#                   baseline), 4 nodes 64x2 core, 8 nodes 128x1 node -- the node counts the debug
#                   partition cannot give, and 1 node 128x1 node (the cross-check that debug and
#                   compute nodes run alike); one allocation_benchmark.batch job each, 20 min limit,
#                   submitted by job 1 when it ends (the normal QOS allows 64 queued).  They sit in
#                   the normal queue and may start later.
# Every job refreshes results.csv and layouts_<network>.png / .txt in
# ~/Documents/lammps_runs/allocation_benchmark/ when it ends, so the figure is waiting whenever
# you come back; `python3 collect_results.py` on the login node does the same.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p logs_slurm chains
QUEUE="$PWD/chains/debug_queue.txt"

NETWORK="${NETWORK:-rot45}"
NS="${NS:-5 15}"
COMPUTE="${COMPUTE:-1}"
NSTEPS="${NSTEPS:-10000}"
PREP_STEPS="${PREP_STEPS:-10000}"
DRY=0; [ "${1:-}" = "--dry-run" ] && DRY=1
INPUT_DIR="$HOME/Documents/lammps_data/input_data"
RUNS="$HOME/Documents/lammps_runs/allocation_benchmark"
mkdir -p "$RUNS/logs"

ONE_NODE="128x1:node:5.0 64x2:node:5.0 32x4:node:5.0 16x8:node:5.0 128x1:node:2.5"
TWO_NODE="128x1:node:5.0 128x1:core:5.0 64x2:node:5.0 64x2:core:5.0 32x4:node:5.0 32x4:core:5.0 16x8:node:5.0 16x8:core:5.0 128x1:core:2.5"
NET_SUFFIX=""; [ "$NETWORK" = random ] || NET_SUFFIX="_${NETWORK}"

for n in $NS; do
    stem="slab_support_periodic_${n}beads_tall_rho04_${NETWORK}"
    fresh="$INPUT_DIR/$stem.data"
    # the rot45 15-bead file in input_data on 2026-11-07 was a byte-identical copy of the 5-bead
    # one; the benchmark uses the file rebuilt with build_beads_sweep.py --network rot45 --beads 15
    rebuilt="$INPUT_DIR/allocation_benchmark/rot45_rebuilt/$stem.data"
    if [ -s "$rebuilt" ]; then fresh="$rebuilt"; fi
    [ -s "$fresh" ] || { echo "missing input file: $fresh"; exit 1; }
    echo "N=$n ($NETWORK): fresh file $fresh"

    j1="dbg_N${n}${NET_SUFFIX}_1node"
    j2="dbg_N${n}${NET_SUFFIX}_2node"
    chain="$PWD/chains/after_${j1}.sh"
    cat > "$chain" <<EOF
#!/usr/bin/env bash
# written by submit_debug_sweep.sh $(date '+%F %T'); run by $j1 when it ends (compute follow-ups only)
cd "$PWD"
EOF
    TWO_CMD="sbatch --parsable --job-name=$j2 --partition=debug --nodes=2 --ntasks-per-node=128 --cpus-per-task=1 --mem=0 --time=00:30:00 --export=ALL,NETWORK=$NETWORK,NBEADS=$n,NSTEPS=$NSTEPS,PREP_STEPS=$PREP_STEPS,TIME_LIMIT_S=1800,LAYOUTS=\"$TWO_NODE\" debug_sweep.batch"
    if [ "$COMPUTE" = "1" ]; then
        # one job per layout on compute, as submit_benchmark.sh does; the prepared state exists by now
        for spec in "4 128 1 node 5.0" "4 64 2 core 5.0" "8 128 1 node 5.0" "1 128 1 node 5.0"; do
            read -r nodes tasks cpus map cut <<<"$spec"
            tag="bench_N${n}_s1.0${NET_SUFFIX}_compute_n${nodes}_t${tasks}_c${cpus}_${map}_cc${cut}_r1"
            cat >> "$chain" <<EOF
sbatch --job-name="$tag" --partition=compute --nodes=$nodes --ntasks-per-node=$tasks --cpus-per-task=$cpus --mem=0 --time=00:20:00 \\
    --export=ALL,ROLE=bench,NETWORK="$NETWORK",NBEADS="$n",SCALE=1.0,MAP=$map,COMM_CUT=$cut,NSTEPS=$NSTEPS,PREP_STEPS=$NSTEPS,REP=1,TAG="$tag" \\
    allocation_benchmark.batch
EOF
        done
    fi

    cmd=(sbatch --parsable --job-name="$j1" --partition=debug --nodes=1 --ntasks-per-node=128 --cpus-per-task=1 --mem=0 --time=00:30:00
         --export=ALL,NETWORK="$NETWORK",NBEADS="$n",NSTEPS="$NSTEPS",PREP_STEPS="$PREP_STEPS",TIME_LIMIT_S=1800,BENCH_FRESH="$fresh",LAYOUTS="$ONE_NODE",NEXT_SCRIPT="$chain"
         debug_sweep.batch)
    if [ "$DRY" = "1" ]; then
        echo "  would submit: ${cmd[*]}"
        echo "  would queue for the waiter: $TWO_CMD"
        echo "  then (from the job): "; sed 's/^/    /' "$chain" | grep -v '^    #'
    else
        jid=$("${cmd[@]}")
        echo "  $j1: job $jid  (1 node; the compute follow-ups are submitted when it ends)"
        # the 2-node job must wait for job 1 (the prepared state) AND for a free debug slot
        echo "$j2|$TWO_CMD --dependency=afterany:$jid" >> "$QUEUE"
        echo "  $j2: queued for chains/debug_waiter.sh (dependency afterany:$jid)"
    fi
done
if [ "$DRY" = "0" ] && [ -s "$QUEUE" ] && ! pgrep -u "$USER" -f chains/debug_waiter.sh >/dev/null; then
    nohup bash chains/debug_waiter.sh >/dev/null 2>&1 &
    disown
    echo "started chains/debug_waiter.sh (nohup); it submits the 2-node jobs as debug slots free up; log: chains/debug_waiter.log"
fi
echo "results: $RUNS/results.csv, $RUNS/layouts_${NETWORK}.png, $RUNS/layouts_${NETWORK}.txt"
