#!/usr/bin/env bash
# submit_node_screen.sh -- health screen of every node of Pod's batch partition (2026-11-19).
# Submits node_screen.sh once per node and mode (core: 1 core, 15 min; node: whole node, 25 min),
# pinned to its node with --nodelist, plus one summary job that runs when all of them have
# finished and writes summary.txt and exclude.txt with summarize_node_screen.py.  Nodes that
# Slurm reports down or drained are listed in nodes.tsv and skipped.
#
#   bash submit_node_screen.sh [--dry-run] [--modes "core node"] [--nodes "node3 node14"] [--out DIR]
#                              [--chain N] [--begin TIME] [--no-summary]
#   --chain N      SCREEN_CHAIN=N for the jobs: postpone up to N times while the node exporter spins, and
#                  let a core-mode job submit the whole-node run once the exporter is quiet (node_screen.sh)
#   --begin TIME   sbatch --begin (e.g. now+6hours)
#   --no-summary   do not add a summary job (chained resubmissions)
#
# Output folder (default): ~/Documents/lammps_runs/allocation_benchmark/node_screen/<YYYYmmdd>/
# Resubmitting part of a screen (--nodes/--modes with the same --out) appends to nodes.tsv, and the new
# summary job also waits for screen jobs of the earlier submission that are still queued.
# The summary can be rerun at any time, also while jobs are still pending:
#   python3 summarize_node_screen.py --out <folder>
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DRY=0; MODES="core node"; NODES=""; OUT=""; CHAIN=""; BEGIN=""; SUMMARY=1
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) DRY=1 ;;
        --modes) MODES="$2"; shift ;;
        --nodes) NODES="$2"; shift ;;
        --out) OUT="$2"; shift ;;
        --chain) CHAIN="$2"; shift ;;
        --begin) BEGIN="$2"; shift ;;
        --no-summary) SUMMARY=0 ;;
        *) echo "usage: bash submit_node_screen.sh [--dry-run] [--modes \"core node\"] [--nodes \"...\"] [--out DIR]" >&2; exit 2 ;;
    esac
    shift
done
OUT="${OUT:-$HOME/Documents/lammps_runs/allocation_benchmark/node_screen/$(date +%Y%m%d)}"
PARTITION=batch
run() { if [ "$DRY" = 1 ]; then echo "DRY: $*"; else "$@"; fi; }

[ "$DRY" = 1 ] || mkdir -p "$OUT" "$HERE/logs_slurm" || exit 1
NODES_TSV="$OUT/nodes.tsv"
[ "$DRY" = 1 ] || [ -s "$NODES_TSV" ] || printf 'node\tstate_at_submit\tjobs\n' > "$NODES_TSV"

# node list with the state Slurm reports now; one line per node
mapfile -t LINES < <(sinfo -p "$PARTITION" -h -N -o "%N %t" | sort -u -V)
IDS=(); NSUB=0; NSKIP=0
for line in "${LINES[@]}"; do
    read -r node state <<< "$line"
    if [ -n "$NODES" ]; then case " $NODES " in *" $node "*) ;; *) continue ;; esac; fi
    case "$state" in
        *down*|*drain*|*drng*|*fail*|*maint*|*unk*)
            echo "skip $node ($state)"; NSKIP=$((NSKIP + 1))
            [ "$DRY" = 1 ] || printf '%s\t%s\t%s\n' "$node" "$state" "-" >> "$NODES_TSV"
            continue ;;
    esac
    jobs=""
    for mode in $MODES; do
        case "$mode" in
            core) RES=(-n 1 -c 1 --time=00:15:00) ;;
            node) RES=(-n 40 --exclusive --time=00:25:00 --signal=B:USR1@120) ;;
            *) echo "unknown mode $mode" >&2; exit 2 ;;
        esac
        EXTRA=(); [ -n "$BEGIN" ] && EXTRA+=(--begin="$BEGIN")
        CMD=(sbatch --parsable -p "$PARTITION" --nodelist="$node" -N 1 "${RES[@]}" "${EXTRA[@]}" -J "screen_${mode}_$node"
             --chdir="$HERE" --export=ALL,CLUSTER=pod,SCREEN_MODE="$mode",SCREEN_OUT="$OUT",SCREEN_CHAIN="$CHAIN" "$HERE/node_screen.sh")
        if [ "$DRY" = 1 ]; then echo "DRY: ${CMD[*]}"; id="dry"
        else id=$("${CMD[@]}") || { echo "sbatch failed for $node $mode" >&2; continue; }; fi
        IDS+=("$id"); jobs="$jobs${jobs:+,}$mode:$id"; NSUB=$((NSUB + 1))
    done
    [ "$DRY" = 1 ] || printf '%s\t%s\t%s\n' "$node" "$state" "$jobs" >> "$NODES_TSV"
done

if [ "${#IDS[@]}" -gt 0 ] && [ "$SUMMARY" = 1 ]; then
    # the summary waits for the new jobs and for screen jobs of an earlier submission still in the queue
    mapfile -t OLD < <(squeue -u "$USER" -h -o "%i %j" | awk '$2 ~ /^screen_(core|node)_/ {print $1}')
    for id in "${OLD[@]}"; do case " ${IDS[*]} " in *" $id "*) ;; *) IDS+=("$id") ;; esac; done
    DEP=$(IFS=:; echo "${IDS[*]}")
    SUM=(sbatch --parsable -p "$PARTITION" -N 1 -n 1 --time=00:10:00 -J screen_summary --chdir="$HERE"
         --dependency="afterany:$DEP" --export=ALL,SCREEN_OUT="$OUT" --output="logs_slurm/%x.o%j.%N"
         --wrap "source /etc/profile.d/modules.sh 2>/dev/null; command -v python3 >/dev/null || module load miniconda; python3 '$HERE/summarize_node_screen.py' --out '$OUT' > '$OUT/summary.txt' 2>&1; python3 '$HERE/summarize_node_screen.py' --current --quiet")
    if [ "$DRY" = 1 ]; then echo "DRY: sbatch ... -J screen_summary --dependency=afterany:<${#IDS[@]} job ids> --wrap \"python3 summarize_node_screen.py --out $OUT > $OUT/summary.txt\""
    else SID=$("${SUM[@]}") && echo "summary job $SID (runs when all ${#IDS[@]} screen jobs have ended)"; fi
fi
echo "submitted $NSUB jobs on $((NSUB / $(wc -w <<< "$MODES"))) nodes, skipped $NSKIP nodes; output: $OUT"
echo "check:  squeue -u $USER -o '%.9i %.22j %.2t %.10M %R' | sort -k2"
echo "table:  python3 $HERE/summarize_node_screen.py --out $OUT"
