#!/usr/bin/env bash
# debug_waiter.sh -- submit the debug jobs listed in debug_queue.txt (name|sbatch command, one per
# line) one at a time, whenever fewer than 2 debug jobs of this user are queued or running
# (debug-normal MaxSubmitPU = 2, and a job cannot submit its successor while it still counts).
# Started by submit_debug_sweep.sh with nohup; it exits when the queue file is empty.
cd "$(dirname "${BASH_SOURCE[0]}")/.."   # the benchmark folder: the sbatch lines name debug_sweep.batch relative to it
Q=chains/debug_queue.txt; LOG=chains/debug_waiter.log
while [ -s "$Q" ]; do
    n=$(squeue -h -u "$USER" -p debug 2>/dev/null | wc -l)
    if [ "$n" -lt 2 ]; then
        line=$(head -n1 "$Q"); name="${line%%|*}"; cmd="${line#*|}"
        if squeue -h -u "$USER" -n "$name" 2>/dev/null | grep -q .; then
            echo "$(date '+%F %T') $name already queued" >> "$LOG"
        else
            jid=$(eval "$cmd" 2>>"$LOG") && echo "$(date '+%F %T') submitted $name: $jid" >> "$LOG" \
                || { echo "$(date '+%F %T') submit of $name failed; retrying in 60 s" >> "$LOG"; sleep 60; continue; }
        fi
        tail -n +2 "$Q" > "$Q.tmp" && mv "$Q.tmp" "$Q"
    fi
    sleep 60
done
echo "$(date '+%F %T') queue empty, waiter done" >> "$LOG"
