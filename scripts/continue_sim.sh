#!/bin/bash
# continue_sim.sh -- resume a run from its final .data file, or finish a run that hit its
# wall-clock budget. Run it from the simulation folder that holds the job's SLURM output file.
#
#   continue_sim.sh <job_id> <nsteps>     extend a FINISHED run by <nsteps> production steps
#                                         (outputs tagged <nsteps>, oldsteps 0)
#   continue_sim.sh --chain <job_id>      finish an INCOMPLETE run: resume from the newest
#                                         snapshot, exactly the remaining steps, original tag.
#                                         run_lammps.sh calls this at the end of a CHAIN=1 job;
#                                         by hand it re-launches a chain that broke. Nothing is
#                                         submitted when the log says STATUS COMPLETE.
#   --batch <file>   the .batch whose #SBATCH lines and knobs the new job borrows (default: the
#                    one here whose --job-name matches the output file, else the batch a chained
#                    job inherited, else <folder>_single.batch)
#   --dry-run        print the batch script that would be submitted, submit nothing
#
# Folders with a cont=1 path: slab_with_support, chain_length_testing, polymer_pure, solvent_pure
# (the two-piston decks have none). The continuation runs in <run dir>/continuation_<stamp>/;
# a continuation of a continuation is a sibling there. Knobs (lib/run_common.sh, plus CHAIN,
# WALL_MARGIN, DECK_FILE): the environment wins, then the `export NAME=` lines of the batch.
#
# The { } group makes bash parse the whole file first (see run_lammps.sh); keep the exit inside.
{
set -eu
SELF="$(realpath "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd "$(dirname "$SELF")" && pwd)"
source "$SCRIPT_DIR/lib/site.sh"
source "$SCRIPT_DIR/lib/run_common.sh"

declare -A FOLDER_PREFIX=([slab_with_support]=final_config [chain_length_testing]=final_config
                          [polymer_pure]=purepol [solvent_pure]=puresolv)

usage() { sed -n '2,21p' "$SELF" | sed 's/^# \{0,1\}//'; exit 1; }

MODE=extend; DRY_RUN=0; BATCH_OVERRIDE=""; POSITIONAL=()
while [ $# -gt 0 ]; do
    case "$1" in
        --chain)   MODE=chain; shift ;;
        --run)     MODE=run; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --batch)   BATCH_OVERRIDE="$2"; shift 2 ;;
        -h|--help) usage ;;
        --*)       echo "Error: unknown option $1"; usage ;;
        *)         POSITIONAL+=("$1"); shift ;;
    esac
done
set -- "${POSITIONAL[@]+"${POSITIONAL[@]}"}"

SIM_FOLDER="$PWD"
FOLDER="$(basename "$PWD")"
DECK_SRC="$PWD/${DECK_FILE:-$FOLDER}.lmp"
PREFIX="${FOLDER_PREFIX[$FOLDER]:-}"
[ -n "$PREFIX" ] || { echo "Error: '$FOLDER' has no continuation path. Supported: ${!FOLDER_PREFIX[*]}"; exit 1; }
[ -f "$DECK_SRC" ] || { echo "Error: no deck at $DECK_SRC. Run this from the simulation folder."; exit 1; }

# ${PREFIX}_<dataname>_<epsSS>_<epsSP>_<totsteps>[_snap<N>].data -> DATANAME INTERACTION EPSSS EPSSP
parse_resume_name() {
    local base parts n
    base="$(basename "$1" .data)"; base="${base%_snap[0-9]*}"; base="${base#${PREFIX}_}"
    IFS='_' read -ra parts <<< "$base"; n=${#parts[@]}
    [ "$n" -ge 4 ] || { echo "Error: cannot parse dataname/epsSS/epsSP/totsteps from $base"; exit 1; }
    EPSSP="${parts[$((n-2))]}"; EPSSS="${parts[$((n-3))]}"
    DATANAME=$(IFS='_'; echo "${parts[*]:0:$((n-3))}")
    INTERACTION="${EPSSS}_${EPSSP}"
}
# value NAME takes in the borrowed batch: its plain or `export` assignment lines (several per line
# allowed, last wins) are evaluated in a subshell, with the current environment as fallback.
batch_value() {
    local line v val=""
    while IFS= read -r line; do
        v=$(set +eu; unset "$1"; eval "$line" 2>/dev/null; printf '%s' "${!1-}")
        if [ -n "$v" ]; then val="$v"; fi
    done < <(grep -E "^(export[[:space:]]+)?[A-Z_0-9]+=" "$BATCH_FILE" 2>/dev/null | sed 's/[[:space:]]#.*//')
    printf '%s' "$val"
}
newest_snapshot() { ls "$1"/${PREFIX}_*_snap[0-9]*.data 2>/dev/null | sort -V | tail -1 || true; }

# ===== submit modes (extend / chain) =============================================================
if [ "$MODE" != run ]; then
    if [ "$MODE" = chain ]; then [ $# -eq 1 ] || usage; else [ $# -eq 2 ] || usage; fi
    JOB_ID="$1"; NSTEPS="${2:-}"

    OUTPUT_FILE=$(ls ./*.o"${JOB_ID}".* 2>/dev/null | head -1 || true)
    [ -z "$OUTPUT_FILE" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ] && OUTPUT_FILE=$(ls "$SLURM_SUBMIT_DIR"/*.o"${JOB_ID}".* 2>/dev/null | head -1 || true)
    [ -n "$OUTPUT_FILE" ] || { echo "Error: no SLURM output file *.o${JOB_ID}.* in $PWD"; exit 1; }
    JOBNAME="$(basename "$OUTPUT_FILE")"; JOBNAME="${JOBNAME%%.o${JOB_ID}.*}"

    WORK_DIR=$(grep "^Working directory:" "$OUTPUT_FILE" | awk '{print $NF}' | head -1 || true)
    [ -d "$WORK_DIR" ] || { echo "Error: run directory from $OUTPUT_FILE not found: '$WORK_DIR'"; exit 1; }
    echo "Mode: $MODE   job $JOB_ID ($JOBNAME)   run dir: $WORK_DIR"

    if [ "$MODE" = chain ]; then
        LINE=$(grep -E '^>>> PRODUCTION_STEPS_DONE [0-9]+ OF [0-9]+ OLDSTEPS [0-9]+ TOTSTEPS [0-9]+ STATUS (COMPLETE|INCOMPLETE)' "$WORK_DIR/log.lammps" 2>/dev/null | tail -1 || true)
        [ -n "$LINE" ] || { echo "Error: no PRODUCTION_STEPS_DONE line in $WORK_DIR/log.lammps (deck without a wall budget, or the job died early). To extend by hand: continue_sim.sh $JOB_ID <nsteps>"; exit 1; }
        read -r _ _ DONE _ _ _ PREV_OLD _ TOTSTEPS _ STATUS <<< "$LINE"
        echo "Log: $LINE"
        [ "$STATUS" = INCOMPLETE ] || { echo "Production complete: nothing to chain."; exit 0; }
        OLDSTEPS=$(( PREV_OLD + DONE )); NSTEPS=$(( TOTSTEPS - OLDSTEPS ))
        [ "$NSTEPS" -gt 0 ] || { echo "Nothing left to run (oldsteps $OLDSTEPS >= totsteps $TOTSTEPS)."; exit 0; }
        RESUME_FILE=$(newest_snapshot "$WORK_DIR")
        [ -n "$RESUME_FILE" ] || { echo "Error: incomplete run but no ${PREFIX}_*_snap<N>.data in $WORK_DIR"; exit 1; }
    else
        OLDSTEPS=0; TOTSTEPS="$NSTEPS"
        RESUME_FILE=$(ls "$WORK_DIR"/${PREFIX}_*.data 2>/dev/null | grep -v '_snap[0-9]*\.data$' | head -1 || true)
        if [ -z "$RESUME_FILE" ]; then
            RESUME_FILE=$(newest_snapshot "$WORK_DIR")
            [ -n "$RESUME_FILE" ] || { echo "Error: no ${PREFIX}_*.data in $WORK_DIR. Did the run finish?"; exit 1; }
            echo "NOTE: no final ${PREFIX}_*.data; resuming from the newest snapshot (for a run that hit its wall budget use --chain)."
        fi
    fi
    parse_resume_name "$RESUME_FILE"
    echo "Resume from: $RESUME_FILE"
    echo "This job: nsteps=$NSTEPS oldsteps=$OLDSTEPS totsteps=$TOTSTEPS"

    if [ -n "$BATCH_OVERRIDE" ]; then
        BATCH_FILE="$(realpath "$BATCH_OVERRIDE")"
    else
        BATCH_FILE=$(grep -lE -- "--job-name=${JOBNAME}([[:space:]]|\$)" ./*.batch 2>/dev/null | head -1 || true)
        [ -n "$BATCH_FILE" ] && BATCH_FILE="$(realpath "$BATCH_FILE")"
        [ -z "$BATCH_FILE" ] && [ -f "${CONT_BATCH_FILE:-}" ] && BATCH_FILE="$CONT_BATCH_FILE"
        if [ -z "$BATCH_FILE" ]; then
            if [ -f "$PWD/${FOLDER}_single.batch" ]; then BATCH_FILE="$PWD/${FOLDER}_single.batch"
            elif [ "$(ls ./*.batch 2>/dev/null | wc -l)" -eq 1 ]; then BATCH_FILE="$(realpath ./*.batch)"; fi
            [ -n "$BATCH_FILE" ] && echo "WARNING: no .batch here has --job-name=$JOBNAME; borrowing $(basename "$BATCH_FILE") (--batch <file> picks another)."
        fi
    fi
    [ -f "${BATCH_FILE:-}" ] || { echo "Error: no .batch to borrow SLURM settings from; pass --batch <file>."; exit 1; }
    echo "Borrowing SLURM settings from: $BATCH_FILE"

    # knobs: environment first, then the batch file
    KNOB_EXPORTS=""
    for k in "${KNOBS[@]}" chain wall_margin deck_file; do
        env="${k^^}"; val="${!env:-$(batch_value "$env")}"
        [ -n "$val" ] && KNOB_EXPORTS+="export $env=\"$val\""$'\n'
    done
    echo "Knobs:"; echo "$KNOB_EXPORTS" | sed 's/^/  /'

    SBATCH_LINES=$(grep '^#SBATCH' "$BATCH_FILE" | grep -v -e '--output' -e '--error' -e '--job-name' || true)
    TMPBATCH=$(mktemp /tmp/cont_XXXXXX.sh)
    cat > "$TMPBATCH" <<BATCHEOF
#!/bin/bash
${SBATCH_LINES}
#SBATCH --job-name=cont_${JOB_ID}
#SBATCH --output=%x.o%j.%N
export LAMMPS_BATCH_SYNCED=1
source "${LAMMPS_WORK_DIR}/scripts/batch_prologue.sh"
${KNOB_EXPORTS}export CONT_BATCH_FILE="${BATCH_FILE}"
cd "${SIM_FOLDER}" || exit 1
bash "${SELF}" --run "${RESUME_FILE}" "${NSTEPS}" "${OLDSTEPS}" "${TOTSTEPS}"
BATCHEOF
    if [ "$DRY_RUN" = 1 ]; then
        echo "DRY RUN: would submit from $SIM_FOLDER:"; echo "------"; cat "$TMPBATCH"; echo "------"
    else
        ( cd "$SIM_FOLDER" && sbatch "$TMPBATCH" )
    fi
    rm -f "$TMPBATCH"
    exit 0
fi

# ===== run mode (inside the SLURM job) ===========================================================
[ $# -eq 4 ] || { echo "Error: --run expects <resume_data_file> <nsteps> <oldsteps> <totsteps>"; exit 1; }
RESUME_FILE="$(realpath "$1")"; NSTEPS="$2"; OLDSTEPS="$3"; TOTSTEPS="$4"
[ -f "$RESUME_FILE" ] || { echo "Error: resume file not found: $RESUME_FILE"; exit 1; }
parse_resume_name "$RESUME_FILE"

ORIG_DIR="$(dirname "$RESUME_FILE")"
case "$(basename "$ORIG_DIR")" in continuation_*) ORIG_DIR="$(dirname "$ORIG_DIR")" ;; esac
STAMP=$(date +%Y%m%d_%H%M%S)
CONT_DIR="$ORIG_DIR/continuation_$STAMP"
TRAJ_DIR="$SCRATCH_DIR/lammps_trajectories/$(basename "$ORIG_DIR")/continuation_$STAMP"
make_run_dirs "$CONT_DIR" "$TRAJ_DIR"
ln -s "$RESUME_FILE" "$CONT_DIR/data_files/$DATANAME.data"
cp "$DECK_SRC" "$CONT_DIR/"
cd "$CONT_DIR" || exit 1

: "${OMP_NUM_THREADS:=${SLURM_CPUS_PER_TASK:-1}}"; export OMP_NUM_THREADS
WALL_TIMEOUT=$(wall_timeout_seconds)
LAMMPS_VARS=(-var dataname "$DATANAME" -var interaction "$INTERACTION" -var epsSS "$EPSSS" -var epsSP "$EPSSP"
             -var nsteps "$NSTEPS" -var oldsteps "$OLDSTEPS" -var totsteps "$TOTSTEPS" -var cont 1
             -var wall_timeout "$WALL_TIMEOUT")
knob_vars
echo "======================================"
echo "Working directory: $CONT_DIR"
echo "Trajectory directory (scratch): $TRAJ_DIR"
echo "Resuming from: $RESUME_FILE"
echo "This job: nsteps=$NSTEPS oldsteps=$OLDSTEPS totsteps=$TOTSTEPS   wall budget ${WALL_TIMEOUT}s   CHAIN ${CHAIN:-0}"
echo "Vars: ${LAMMPS_VARS[*]}"
echo "======================================"

set +e
launch_lammps "$CONT_DIR/$(basename "$DECK_SRC")" "${LAMMPS_VARS[@]}"
[ "$LAMMPS_STATUS" = failed ] && exit "$LAMMPS_RC"
chain_if_incomplete "$SIM_FOLDER"
export_knobs
bash "$SCRIPT_DIR/postprocess.sh" "$CONT_DIR" "$FOLDER" "$DATANAME" "$INTERACTION" "$TOTSTEPS" "$OLDSTEPS" \
    || echo "WARNING: postprocess.sh reported a failure; regenerate by hand (see scripts/postprocess.sh)."
[ "$LAMMPS_STATUS" = partial ] && exit "$LAMMPS_RC"
echo "Done: $CONT_DIR"
exit 0
}
