#!/bin/bash
# submit.sh -- sbatch a .batch file on whichever cluster this is (Expanse or Pod), from any folder.
# `bash install.sh` (repo root, once per machine) installs it as the shell command `lsbatch`:
#
#   lsbatch [sbatch options] <file.batch>
#   lsbatch uniaxial_compression_single               the name alone is enough, with or without .batch
#   DATANAME=<stem> lsbatch --nodes=2 uniaxial_compression_single.batch
#   lsbatch --dry-run <file.batch>                    print what would be submitted, submit nothing
#   CLUSTER=pod lsbatch ...                           force the cluster (default: detected, lib/site.sh)
#   (without the install: bash ~/Documents/lammps_work/scripts/submit.sh ...)
#
# The batch file is looked up by name under simulations/ when it is not a path that exists here,
# and the job is submitted from the batch file's own folder, so its SLURM output file lands
# there (where continue_sim.sh looks for it) wherever the command was typed.
#
# A batch file holds its Expanse allocation as #SBATCH lines and its Pod allocation on one line:
#     #POD --nodes=4 --ntasks-per-node=40 --cpus-per-task=1 --time=96:00:00
# sbatch reads #SBATCH lines before any of the script runs, so the script itself cannot choose
# between the two; this wrapper does. On Expanse it is exactly `sbatch [options] <file.batch>`.
# On Pod it submits the same script with the #SBATCH allocation lines replaced by the #POD ones
# (--job-name and --output are kept). Options given here win over both, as with sbatch:
#     bash scripts/submit.sh --nodes=2 --time=12:00:00 <file.batch>
# Paths, modules, the LAMMPS binary and the MPI settings of each cluster are in lib/site.sh.
{
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/site.sh" || exit 1

DRY=0; OPTS=()
while [ $# -gt 1 ]; do
    case "$1" in --dry-run) DRY=1 ;; *) OPTS+=("$1") ;; esac
    shift
done
BATCH="${1:-}"
if [ -z "$BATCH" ] || [ "$BATCH" = -h ] || [ "$BATCH" = --help ]; then
    sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1
fi
if [ ! -f "$BATCH" ]; then                     # not a path from here: find it by name under simulations/
    NAME="$(basename "$BATCH")"; NAME="${NAME%.batch}.batch"
    FOUND=$(ls "$LAMMPS_WORK_DIR"/simulations/*/"$NAME" 2>/dev/null)
    case "$(printf '%s' "$FOUND" | grep -c .)" in
        1) BATCH="$FOUND" ;;
        0) echo "Error: no $NAME here or in $LAMMPS_WORK_DIR/simulations/*/"; exit 1 ;;
        *) echo "Error: $NAME exists in several folders; give the path:"; printf '%s\n' "$FOUND" | sed 's/^/    /'; exit 1 ;;
    esac
fi
BATCH="$(cd "$(dirname "$BATCH")" && pwd)/$(basename "$BATCH")"
cd "$(dirname "$BATCH")" || exit 1
export CLUSTER

if [ "$CLUSTER" = expanse ]; then
    echo ">>> $CLUSTER: sbatch ${OPTS[*]+"${OPTS[*]}"} $BATCH"
    [ "$DRY" = 1 ] && { grep '^#SBATCH' "$BATCH"; exit 0; }
    exec sbatch ${OPTS[@]+"${OPTS[@]}"} "$BATCH"
fi

TAG="#$(printf '%s' "$CLUSTER" | tr '[:lower:]' '[:upper:]')"
HEADER=$(batch_header "$BATCH") || {
    echo "Error: $(basename "$BATCH") has no '$TAG ...' line, so it has no allocation for $CLUSTER."
    echo "       Add one under its #SBATCH lines, e.g.:  $TAG --nodes=1 --ntasks-per-node=$CORES_PER_NODE --cpus-per-task=1 --time=24:00:00"
    exit 1
}
# The job script goes to sbatch on stdin, so Slurm does not record the file's path; the prologue
# takes it from LAMMPS_BATCH_FILE to re-run the freshly pulled copy of the batch.
export LAMMPS_BATCH_FILE="$BATCH"
JOB_SCRIPT="$(printf '%s\n%s\n' '#!/usr/bin/env bash' "$HEADER"; sed '1{/^#!/d;}' "$BATCH" | grep -v '^#SBATCH')"
echo ">>> $CLUSTER: $(basename "$BATCH") with"
printf '%s\n' "$HEADER" | sed 's/^/      /'
[ ${#OPTS[@]} -gt 0 ] && echo "    and on the command line: ${OPTS[*]}"
[ "$DRY" = 1 ] && exit 0
printf '%s\n' "$JOB_SCRIPT" | sbatch ${OPTS[@]+"${OPTS[@]}"}
exit $?
}
