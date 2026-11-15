#!/usr/bin/env bash
# wide10_status.sh -- where the wide random 10-bead pipeline stands. Read-only.
PIPE=$HOME/Documents/lammps_runs/chain_length_testing/wide10_pipeline_20261114
echo "=== pipeline log"; cat "$PIPE/status.log" 2>/dev/null
echo; echo "=== jobs (sacct)"
IDS=$(cut -d'|' -f1 "$PIPE/jobs.txt" 2>/dev/null | paste -sd, -)
[ -n "$IDS" ] && timeout 120 sacct -X -n -j "$IDS" -o jobid,jobname%24,partition,state%16,start,elapsed
echo; echo "=== still queued or running"; timeout 60 squeue -u "$USER" -o "%.10i %.9P %.24j %.2t %.10M %.20S %R"
RUN=$(ls -dt "$HOME"/Documents/lammps_runs/chain_length_testing/chain_length_testing_slab_support_periodic_10beads_wide_rho04_random_1.0_1.0_* 2>/dev/null | head -1)
if [ -n "$RUN" ]; then
    echo; echo "=== equilibration: $RUN"
    grep -E '^>>> (RG_BLOCK|RG_PLATEAU|MEAN_BOX|PRODUCTION_STEPS_DONE)' "$RUN/log.lammps" 2>/dev/null | cut -c1-200 | tail -8
fi
echo; echo "=== report: $HOME/Documents/lammps_runs/chain_length_testing/screen_report.md (row 'random 10 WIDE' once the bridge has run)"
