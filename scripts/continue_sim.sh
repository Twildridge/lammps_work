#!/bin/bash
# continue_sim.sh — Continue a LAMMPS simulation from a final .data file, or from
# the newest snapshot of a run that hit its wall-clock budget.
#
# Run this from inside the simulation folder (the one holding the SLURM output
# file of the job you are continuing). The script reads the working directory
# from that output file, picks the .data file to resume from, and launches a
# continuation run — skipping setup phases via cont=1. This is a REAL restart
# (skip setup, keep going), not a fresh resubmission — see
# docs/running_simulations.md §5e for how this differs from editing a .batch
# file's NSTEPS and resubmitting.
#
# Usage (from inside the simulation folder):
#   continue_sim.sh <job_id> <nsteps>      extend a FINISHED run by <nsteps> more
#                                          production steps (outputs tagged <nsteps>,
#                                          oldsteps = 0 — the classic use)
#   continue_sim.sh --chain <job_id>       finish a run that ran out of wall-clock
#                                          (2026-09-30): reads the deck's
#                                          `PRODUCTION_STEPS_DONE ... INCOMPLETE`
#                                          line from the job's log, resumes from
#                                          the newest final_config_*_snap<N>.data
#                                          and asks for exactly the remaining
#                                          steps, keeping the original totsteps
#                                          tag on every output file. run_lammps.sh
#                                          (and this script) call it at the end of
#                                          a CHAIN=1 job; by hand it re-launches a
#                                          chain that broke. Nothing is submitted
#                                          when the log says COMPLETE.
# Options:
#   --batch <file>   the .batch whose #SBATCH lines and exports the continuation
#                    job borrows (default: the .batch in this folder whose
#                    --job-name matches the output file's prefix, else the
#                    CONT_BATCH_FILE a chained job inherited, else <folder>.batch)
#   --dry-run        print the batch script that would be submitted; submit nothing
#
#   job_id  — the SLURM job ID, e.g. 49772594
#             (from the output file <jobname>.o49772594.exp-14-05)
#   nsteps  — additional timesteps to run
#
# Examples:
#   cd ~/Documents/lammps_work/simulations/slab_with_support
#   continue_sim.sh 49800123 500000
#   continue_sim.sh --chain 54398786
#
#   cd ~/Documents/lammps_work/simulations/shear_slab
#   continue_sim.sh 49900456 2000000
#
# Supported folders and what cont=1 means for each (triaxial_compression and
# triaxial_permeation are the ONE-PISTON decks, archived 2026-09-22 under
# archive/simulations/ -- still supported from a folder holding the original
# SLURM output file next to the archived .lmp/.batch; the two-piston decks
# triaxial_*_two_pist have no continuation path at all):
#   slab_with_support               — skip push-off/minimize/gentle-NVT/NPT warm-up; run more NPT production
#   solvent_pure, polymer_pure      — skip pre-relax/gentle-ramp stages; run more NPT production
#   triaxial_compression (sweep)    — auto-detects the last _c<level> reached, skips the
#                                     non-equilibrium drive, extends the equilibration
#                                     hold at THAT level only (never re-sweeps)
#   shear_slab (sweep)              — auto-detects the last _g<strain> reached, skips the
#                                     non-equilibrium shear drive, extends the production
#                                     hold at THAT strain only (never re-sweeps)
#   triaxial_permeation             — NOT a sweep; skips Phase 0/0.5/1.5 AND the piston
#                                     reposition/WCA-relax/force-ramp, resuming the
#                                     constant-pressure forcing drive for more steps
#                                     (continuation here means KEEP FORCING, never a hold)
#
# Not supported: solvent_phase/polymer_phase (internal P-sweeps complete in one
# invocation — "continuing" isn't a meaningful operation for them) or volmix_sweep
# (its own SLURM-chained orchestration, not a fit for this tool). compress_slab
# is a separate project — ask its owner before adding support here.
#
# Output goes into a continuation subfolder inside the original run directory;
# a continuation of a continuation is a SIBLING there, never nested:
#   ~/Documents/lammps_runs/{folder}/{original_run_dir}/continuation_{timestamp}/
#
# Environment the continuation job inherits from the batch it borrows (parsed
# from `export NAME=` / `NAME=` lines unless already set in the environment):
#   BARO_MODE PIN_GEL PISTON_TRANSPARENT PRESS_TARGET CHAIN WALL_MARGIN
# Inside the job the same wall-clock budget as run_lammps.sh applies
# (-var wall_timeout = SLURM end - WALL_MARGIN, default 2400 s), so a chained
# job that runs out of time chains again.

# ── Whole-file parse guard ────────────────────────────────────────────────────
# Same guard as run_lammps.sh (see the note there): this script runs mpirun for
# hours on a compute node while git pulls from other NFS clients can replace
# the file underneath the running bash, which then hits a stale handle, treats
# it as EOF and exits silently before the post-processing lines. The { ... }
# group makes bash parse the whole file up front; the trailing `exit` stops it
# reading past the group. (2026-09-26)
{
set -e

# ── Per-folder output-file prefix (set by each .lmp script's write_data/
#    write_restart at the very end) and, for sweep folders, the per-level
#    output tag used to auto-detect which level to continue.
declare -A FOLDER_PREFIX=(
    [slab_with_support]="final_config"
    [triaxial_compression]="final_tricomp"
    [triaxial_permeation]="final_triperm"
    [shear_slab]="final_shear"
    [solvent_pure]="puresolv"
    [polymer_pure]="purepol"
)
declare -A FOLDER_SWEEP_TAG=(   # empty/unset = not a sweep folder
    [triaxial_compression]="_c"
    [shear_slab]="_g"
)
declare -A FOLDER_SWEEP_VAR=(   # LAMMPS -var name that carries the sweep value
    [triaxial_compression]="compressions"
    [shear_slab]="strains"
)

LMP_BIN=/home/dpollard/software/lammps/22Jul2025_update3/mpi-omp/gcc/10.2.0/openmpi/4.1.3/lammps-22Jul2025/build/lmp

usage() {
    echo "Usage: continue_sim.sh [--batch <file>] [--dry-run] <job_id> <nsteps>"
    echo "       continue_sim.sh [--batch <file>] [--dry-run] --chain <job_id>"
    echo ""
    echo "  job_id  — SLURM job ID (e.g. 49772594)"
    echo "  nsteps  — additional timesteps to run (extend mode)"
    echo "  --chain — finish a run that hit its wall-clock budget: remaining steps,"
    echo "            resume from the newest snapshot, original totsteps tag kept"
    echo ""
    echo "Run from inside the simulation folder:"
    echo "  cd ~/Documents/lammps_work/simulations/slab_with_support"
    echo "  continue_sim.sh 49772594 500000"
    echo "  continue_sim.sh --chain 54398786"
    exit 1
}

# ── Arguments ─────────────────────────────────────────────────────────────────
MODE=extend          # extend | chain | run (run = internal, inside the SLURM job)
DRY_RUN=0
BATCH_OVERRIDE=""
POSITIONAL=()
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
set -- "${POSITIONAL[@]}"

# ── Locate scripts and LAMMPS file ────────────────────────────────────────────
SELF="$(realpath "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$( cd "$( dirname "$SELF" )" && pwd )"
LAMMPS_WORK_DIR="$(dirname "$SCRIPT_DIR")"

# Folder is the name of the directory we're running from
SIM_FOLDER="$PWD"
FOLDER="$(basename "$PWD")"
LAMMPS_FILE="$PWD/${FOLDER}.lmp"

if [ -z "${FOLDER_PREFIX[$FOLDER]:-}" ]; then
    echo "Error: unsupported folder '$FOLDER'."
    echo "  Supported: ${!FOLDER_PREFIX[*]}"
    echo "  Not supported: solvent_phase, polymer_phase (internal P-sweeps, no"
    echo "  meaningful 'continue'), volmix_sweep (its own SLURM-chained pipeline)."
    exit 1
fi

if [ ! -f "$LAMMPS_FILE" ]; then
    echo "Error: no LAMMPS script found at $LAMMPS_FILE"
    echo "  Are you inside the correct simulation folder?"
    exit 1
fi

PREFIX="${FOLDER_PREFIX[$FOLDER]}"
SWEEP_TAG="${FOLDER_SWEEP_TAG[$FOLDER]:-}"

# ── Helpers ───────────────────────────────────────────────────────────────────
# Naming convention set by write_data at the end of each LAMMPS script:
#   ${prefix}_{dataname}_{epsSS}_{epsSP}_{totsteps}[_snap<N>].data
# Split by '_' from the right (after dropping a _snap<N> suffix):
#   index -1 = totsteps tag (integer)
#   index -2 = epsSP       (float, e.g. 1.0)
#   index -3 = epsSS       (float, e.g. 1.0)
#   rest     = dataname parts
# Sets DATANAME, INTERACTION, EPSSS, EPSSP, FILE_TAG.
parse_resume_name() {
    local basename rest n
    basename="$(basename "$1" .data)"
    basename="${basename%_snap[0-9]*}"
    rest="${basename#${PREFIX}_}"
    local parts
    IFS='_' read -ra parts <<< "$rest"
    n=${#parts[@]}
    if [ "$n" -lt 4 ]; then
        echo "Error: cannot parse dataname/epsSS/epsSP/totsteps from: $basename"
        exit 1
    fi
    FILE_TAG="${parts[$((n-1))]}"
    EPSSP="${parts[$((n-2))]}"
    EPSSS="${parts[$((n-3))]}"
    DATANAME=$(IFS='_'; echo "${parts[*]:0:$((n-3))}")
    INTERACTION="${EPSSS}_${EPSSP}"
}

# Value of `export NAME=...` or `NAME=...` in the batch file (last one wins,
# trailing comment stripped, quotes dropped); empty if absent.
batch_value() {
    [ -n "${BATCH_FILE:-}" ] && [ -f "$BATCH_FILE" ] || { echo ""; return 0; }
    grep -E "^(export )?$1=" "$BATCH_FILE" 2>/dev/null | tail -1 \
        | sed -E "s/^(export )?$1=//; s/#.*//" | awk '{print $1}' | tr -d '"'"'" || true
}

# Newest snapshot of a run directory (highest _snap<N>), empty if none.
newest_snapshot() {
    ls "$1"/${PREFIX}_*_snap[0-9]*.data 2>/dev/null | sort -V | tail -1 || true
}

# Same wall-clock budget as run_lammps.sh: seconds LAMMPS may run before the
# deck ends its production loop and reports INCOMPLETE (0 = unlimited).
wall_timeout_seconds() {
    local margin="${WALL_MARGIN:-2400}" budget=0
    if [ -n "${SLURM_JOB_END_TIME:-}" ]; then
        budget=$(( SLURM_JOB_END_TIME - $(date +%s) - margin ))
        [ "$budget" -gt 0 ] || budget=0
    fi
    echo "$budget"
}

# =============================================================================
#   SUBMIT MODES (extend / chain): resolve everything, write a batch, sbatch it
# =============================================================================
if [ "$MODE" != "run" ]; then
    if [ "$MODE" = "chain" ]; then
        [ $# -eq 1 ] || usage
    else
        [ $# -eq 2 ] || usage
    fi
    JOB_ID="$1"
    NSTEPS="${2:-}"

    # ── Find the SLURM output file (here, or where the job was submitted from) ──
    OUTPUT_FILE=$(ls *.o${JOB_ID}.* 2>/dev/null | head -1 || true)
    if [ -z "$OUTPUT_FILE" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
        OUTPUT_FILE=$(ls "${SLURM_SUBMIT_DIR}"/*.o${JOB_ID}.* 2>/dev/null | head -1 || true)
    fi
    if [ -z "$OUTPUT_FILE" ]; then
        echo "Error: no SLURM output file matching *.o${JOB_ID}.* found in $PWD"
        echo "  Make sure you are in the folder where you ran sbatch."
        exit 1
    fi
    JOBNAME="$(basename "$OUTPUT_FILE")"
    JOBNAME="${JOBNAME%%.o${JOB_ID}.*}"

    echo "======================================"
    echo "Mode              : $MODE"
    echo "SLURM output file : $OUTPUT_FILE  (job name: $JOBNAME)"

    # ── Parse working directory from output file ──────────────────────────────
    # run_lammps.sh and this script print: "Working directory: /path/to/work_dir"
    WORK_DIR=$(grep "^Working directory:" "$OUTPUT_FILE" | awk '{print $NF}' | head -1 || true)
    if [ -z "$WORK_DIR" ]; then
        echo "Error: could not find working directory line in output file."
        echo "  Looked for line: 'Working directory: /path/...'"
        echo "  Check that the job started and the output file is intact."
        exit 1
    fi
    # New layout: lammps_runs/{FOLDER}/{FOLDER}_... — fall back if old flat path not found
    if [ ! -d "$WORK_DIR" ]; then
        WORK_DIR_SUB="$(dirname "$WORK_DIR")/${FOLDER}/$(basename "$WORK_DIR")"
        if [ -d "$WORK_DIR_SUB" ]; then
            echo "Note: path in log was flat; found run dir under ${FOLDER}/ subfolder."
            WORK_DIR="$WORK_DIR_SUB"
        else
            echo "Error: could not find working directory."
            echo "  From log  : $WORK_DIR"
            echo "  Also tried: $WORK_DIR_SUB"
            exit 1
        fi
    fi
    echo "Resume from dir   : $WORK_DIR"

    # ── Pick the resume file and the step bookkeeping ─────────────────────────
    if [ "$MODE" = "chain" ]; then
        LOG="$WORK_DIR/log.lammps"
        PROD_LINE=$(grep -E '^>>> PRODUCTION_STEPS_DONE [0-9]+ OF [0-9]+ OLDSTEPS [0-9]+ TOTSTEPS [0-9]+ STATUS (COMPLETE|INCOMPLETE)' "$LOG" 2>/dev/null | tail -1 || true)
        if [ -z "$PROD_LINE" ]; then
            echo "Error: no 'PRODUCTION_STEPS_DONE ... STATUS ...' line in $LOG"
            echo "  Either the deck has no wall-clock budget (pre-2026-09-30 slab_with_support.lmp,"
            echo "  or another folder), or the job died before the production loop ended."
            echo "  To extend by hand: continue_sim.sh $JOB_ID <nsteps>"
            exit 1
        fi
        read -r _ _ DONE _ REQUESTED _ PREV_OLD _ PREV_TOT _ STATUS <<< "$PROD_LINE"
        echo "Log bookkeeping   : done ${DONE} of ${REQUESTED} (oldsteps ${PREV_OLD}, totsteps ${PREV_TOT}) — ${STATUS}"
        if [ "$STATUS" = "COMPLETE" ]; then
            echo "Production complete — nothing to chain."
            exit 0
        fi
        OLDSTEPS=$(( PREV_OLD + DONE ))
        TOTSTEPS="$PREV_TOT"
        NSTEPS=$(( TOTSTEPS - OLDSTEPS ))
        if [ "$NSTEPS" -le 0 ]; then
            echo "Nothing left to run (oldsteps ${OLDSTEPS} >= totsteps ${TOTSTEPS})."
            exit 0
        fi
        RESUME_FILE=$(newest_snapshot "$WORK_DIR")
        if [ -z "$RESUME_FILE" ]; then
            echo "Error: production incomplete but no ${PREFIX}_*_snap<N>.data in $WORK_DIR"
            exit 1
        fi
    else
        OLDSTEPS=0
        TOTSTEPS="$NSTEPS"
        RESUME_FILE=$(ls "$WORK_DIR"/${PREFIX}_*.data 2>/dev/null | grep -v '_snap[0-9]*\.data$' | head -1 || true)
        if [ -z "$RESUME_FILE" ]; then
            RESUME_FILE=$(newest_snapshot "$WORK_DIR")
            if [ -z "$RESUME_FILE" ]; then
                echo "Error: no final .data file found in $WORK_DIR"
                echo "  Expected: ${PREFIX}_*.data"
                echo "  Did the LAMMPS run finish successfully?"
                exit 1
            fi
            echo "NOTE: no final ${PREFIX}_*.data (run cut short?) — resuming from the newest snapshot."
            echo "      If the run hit its wall-clock budget, 'continue_sim.sh --chain ${JOB_ID}' finishes"
            echo "      it with the original step count instead of adding ${NSTEPS} on top."
        fi
    fi
    echo "Resume data file  : $RESUME_FILE"

    parse_resume_name "$RESUME_FILE"
    echo "Data name         : $DATANAME"
    echo "Interaction       : $INTERACTION  (epsSS=$EPSSS  epsSP=$EPSSP)"
    echo "File tag (steps)  : $FILE_TAG"
    echo "This job          : nsteps=${NSTEPS}  oldsteps=${OLDSTEPS}  totsteps=${TOTSTEPS}"

    # ── Sweep folders: auto-detect the last level reached ─────────────────────
    # For triaxial_compression/shear_slab, continuing must extend the equilibration
    # hold at whatever level the original run last completed — never re-drive
    # through the whole strain/compression ladder. Scan the run's
    # output_files/stress_data/ for the highest _c<level>/_g<level> tag present.
    SWEEP_LEVEL=""
    if [ -n "$SWEEP_TAG" ]; then
        SWEEP_LEVEL=$(ls "$WORK_DIR"/output_files/stress_data/*"${SWEEP_TAG}"[0-9]*.dat 2>/dev/null \
            | sed -E "s/.*${SWEEP_TAG}([0-9.]+)\.dat\$/\1/" \
            | sort -g | uniq | tail -1 || true)
        if [ -z "$SWEEP_LEVEL" ]; then
            echo "Error: could not auto-detect the last ${SWEEP_TAG}<level> reached in"
            echo "  $WORK_DIR/output_files/stress_data/"
            exit 1
        fi
        echo "Last level reached : ${SWEEP_TAG}${SWEEP_LEVEL}"
    fi

    # ── The .batch to borrow SLURM settings + exports from ────────────────────
    BATCH_FILE=""
    if [ -n "$BATCH_OVERRIDE" ]; then
        BATCH_FILE="$(realpath "$BATCH_OVERRIDE")"
    else
        # the batch in this folder whose --job-name is the output file's prefix
        BATCH_FILE=$(grep -lE -- "--job-name=${JOBNAME}([[:space:]]|\$)" ./*.batch 2>/dev/null | head -1 || true)
        [ -n "$BATCH_FILE" ] && BATCH_FILE="$(realpath "$BATCH_FILE")"
        if [ -z "$BATCH_FILE" ] && [ -n "${CONT_BATCH_FILE:-}" ] && [ -f "$CONT_BATCH_FILE" ]; then
            BATCH_FILE="$CONT_BATCH_FILE"      # a chained job: the batch its ancestor used
        fi
        if [ -z "$BATCH_FILE" ] && [ -f "$PWD/${FOLDER}.batch" ]; then
            BATCH_FILE="$PWD/${FOLDER}.batch"
            echo "WARNING: no .batch with --job-name=${JOBNAME} here — borrowing ${FOLDER}.batch"
            echo "         (check its --nodes/--time; pass --batch <file> to choose another)."
        fi
    fi
    if [ -z "$BATCH_FILE" ] || [ ! -f "$BATCH_FILE" ]; then
        echo "Error: no .batch file to borrow SLURM settings from (looked for --job-name=${JOBNAME},"
        echo "  \$CONT_BATCH_FILE and ${FOLDER}.batch in $PWD). Pass --batch <file>."
        exit 1
    fi
    echo "Borrowing SLURM   : $BATCH_FILE"

    # Deck knobs: the environment (what the calling job ran with) wins, then the
    # batch file, then the run_lammps.sh defaults.
    BARO_MODE="${BARO_MODE:-$(batch_value BARO_MODE)}";                 BARO_MODE="${BARO_MODE:-aniso}"
    PIN_GEL="${PIN_GEL:-$(batch_value PIN_GEL)}";                       PIN_GEL="${PIN_GEL:-0}"
    PISTON_TRANSPARENT="${PISTON_TRANSPARENT:-$(batch_value PISTON_TRANSPARENT)}"; PISTON_TRANSPARENT="${PISTON_TRANSPARENT:-0}"
    PRESS_TARGET="${PRESS_TARGET:-$(batch_value PRESS_TARGET)}";        PRESS_TARGET="${PRESS_TARGET:-1.5}"
    CHAIN="${CHAIN:-$(batch_value CHAIN)}";                             CHAIN="${CHAIN:-0}"
    WALL_MARGIN="${WALL_MARGIN:-$(batch_value WALL_MARGIN)}";           WALL_MARGIN="${WALL_MARGIN:-2400}"
    echo "Deck knobs        : BARO_MODE=${BARO_MODE} PIN_GEL=${PIN_GEL} PISTON_TRANSPARENT=${PISTON_TRANSPARENT} PRESS_TARGET=${PRESS_TARGET} CHAIN=${CHAIN} WALL_MARGIN=${WALL_MARGIN}"
    echo "======================================"

    # Extract #SBATCH directives (skip --output/--error/--job-name so we set our own)
    SBATCH_LINES=$(grep '^#SBATCH' "$BATCH_FILE" | grep -v -e '--output' -e '--error' -e '--job-name' || true)
    MODULE_LINES=$(grep '^module' "$BATCH_FILE" || true)
    ENV_LINES=$(grep '^declare' "$BATCH_FILE" || true)

    TMPBATCH=$(mktemp /tmp/cont_XXXXXX.sh)
    cat > "$TMPBATCH" << BATCHEOF
#!/bin/bash
${SBATCH_LINES}
#SBATCH --job-name=cont_${JOB_ID}
#SBATCH --output=%x.o%j.%N

${ENV_LINES}

${MODULE_LINES}

# knobs inherited from ${BATCH_FILE}
export BARO_MODE="${BARO_MODE}"
export PIN_GEL="${PIN_GEL}"
export PISTON_TRANSPARENT="${PISTON_TRANSPARENT}"
export PRESS_TARGET="${PRESS_TARGET}"
export CHAIN="${CHAIN}"
export WALL_MARGIN="${WALL_MARGIN}"
export CONT_BATCH_FILE="${BATCH_FILE}"

cd "${SIM_FOLDER}" || exit 1
bash "${SELF}" --run "${RESUME_FILE}" "${NSTEPS}" "${OLDSTEPS}" "${TOTSTEPS}" "${SWEEP_LEVEL}"
BATCHEOF

    if [ "$DRY_RUN" = "1" ]; then
        echo "DRY RUN — would submit this batch script from ${SIM_FOLDER}:"
        echo "--------------------------------------"
        cat "$TMPBATCH"
        echo "--------------------------------------"
        rm -f "$TMPBATCH"
        exit 0
    fi

    echo "Submitting continuation as SLURM batch job (from ${SIM_FOLDER})..."
    ( cd "$SIM_FOLDER" && sbatch "$TMPBATCH" )
    rm -f "$TMPBATCH"
    exit 0
fi

# =============================================================================
#   RUN MODE (inside the SLURM job): set up the continuation dir and run LAMMPS
# =============================================================================
if [ $# -lt 4 ]; then
    echo "Error: --run expects <resume_data_file> <nsteps> <oldsteps> <totsteps> [sweep_level]"
    exit 1
fi
RESUME_FILE="$(realpath "$1")"
NSTEPS="$2"
OLDSTEPS="$3"
TOTSTEPS="$4"
SWEEP_LEVEL="${5:-}"
if [ ! -f "$RESUME_FILE" ]; then
    echo "Error: resume data file not found: $RESUME_FILE"
    exit 1
fi
parse_resume_name "$RESUME_FILE"
PRESS_TARGET="${PRESS_TARGET:-1.5}"
BARO_MODE="${BARO_MODE:-aniso}"
PIN_GEL="${PIN_GEL:-0}"
PISTON_TRANSPARENT="${PISTON_TRANSPARENT:-0}"

# ── Set up continuation directory beside any earlier continuations ───────────
WORK_DIR="$(dirname "$RESUME_FILE")"
ORIG_DIR="$WORK_DIR"
case "$(basename "$WORK_DIR")" in
    continuation_*) ORIG_DIR="$(dirname "$WORK_DIR")" ;;
esac
RUN_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SCRATCH_DIR="/expanse/lustre/scratch/$USER/temp_project"

CONT_DIR="${ORIG_DIR}/continuation_${RUN_TIMESTAMP}"
ORIG_TRAJ_BASE=$(basename "$ORIG_DIR")
TRAJ_DIR="$SCRATCH_DIR/lammps_trajectories/${ORIG_TRAJ_BASE}/continuation_${RUN_TIMESTAMP}"

mkdir -p "$CONT_DIR"/{data_files,output_files/{stress_data,volume_data,piston_data,permeation_data,displacement_data,pair_data,chemical_potential},output_plots}
mkdir -p "$TRAJ_DIR"
ln -s "$TRAJ_DIR" "$CONT_DIR/traj_files"

# Symlink the data file under the name LAMMPS expects
ln -s "$RESUME_FILE" "$CONT_DIR/data_files/${DATANAME}.data"

# Snapshot the deck into the run dir and run THAT copy (see run_lammps.sh: LAMMPS
# reads its input lazily, so a mid-run edit of the repo file would corrupt it).
cp "$LAMMPS_FILE" "$CONT_DIR/$(basename "$LAMMPS_FILE")"
LAMMPS_FILE="$CONT_DIR/$(basename "$LAMMPS_FILE")"

echo "======================================"
echo "Working directory: $CONT_DIR"
echo "Trajectory directory (scratch): $TRAJ_DIR"
echo "Resuming from     : $RESUME_FILE"
echo "Data name         : $DATANAME"
echo "Interaction       : $INTERACTION  (epsSS=$EPSSS  epsSP=$EPSSP)"
echo "This job          : nsteps=${NSTEPS}  oldsteps=${OLDSTEPS}  totsteps=${TOTSTEPS}"
[ -n "$SWEEP_LEVEL" ] && echo "Sweep level       : ${SWEEP_TAG}${SWEEP_LEVEL}"
echo "Deck knobs        : BARO_MODE=${BARO_MODE} PIN_GEL=${PIN_GEL} PISTON_TRANSPARENT=${PISTON_TRANSPARENT} PRESS_TARGET=${PRESS_TARGET} CHAIN=${CHAIN:-0}"
echo "======================================"

cd "$CONT_DIR" || exit 1

# ── Build LAMMPS variable list ────────────────────────────────────────────────
WALL_TIMEOUT=$(wall_timeout_seconds)
echo "LAMMPS wall-clock budget: ${WALL_TIMEOUT}s (-var wall_timeout = SLURM end - ${WALL_MARGIN:-2400}s)"
LAMMPS_VARS=(
    -var dataname    "$DATANAME"
    -var interaction "$INTERACTION"
    -var epsSS       "$EPSSS"
    -var epsSP       "$EPSSP"
    -var nsteps      "$NSTEPS"
    -var oldsteps    "$OLDSTEPS"
    -var totsteps    "$TOTSTEPS"
    -var cont        1
    -var press_target "$PRESS_TARGET"
    -var baro_mode   "$BARO_MODE"
    -var pin_gel     "$PIN_GEL"
    -var piston_transparent "$PISTON_TRANSPARENT"
    -var wall_timeout "$WALL_TIMEOUT"
)
if [ -n "$SWEEP_TAG" ] && [ -n "$SWEEP_LEVEL" ]; then
    LAMMPS_VARS+=(-var "${FOLDER_SWEEP_VAR[$FOLDER]}" "$SWEEP_LEVEL")
fi

# ── Run LAMMPS (same teardown-hang guard as run_lammps.sh) ───────────────────
LAMMPS_LOG="${CONT_DIR}/log.lammps"
MPIRUN_TIMEOUT=""
if [ -n "${SLURM_JOB_END_TIME:-}" ]; then
    REMAIN=$(( SLURM_JOB_END_TIME - $(date +%s) - 300 ))   # 5-min buffer before wall limit
    if [ "$REMAIN" -gt 0 ]; then
        MPIRUN_TIMEOUT="timeout -k 60 ${REMAIN}s"
        echo "mpirun guarded by: ${MPIRUN_TIMEOUT} (force-returns ~5 min before wall limit)"
    fi
fi
echo "Running LAMMPS continuation ($SLURM_NTASKS tasks)..."
set +e
$MPIRUN_TIMEOUT mpirun -n "${SLURM_NTASKS}" \
    --bind-to "${OMPI_UNIT}" \
    --map-by "node:pe=${OMP_NUM_THREADS}" \
    "$LMP_BIN" \
    -sf omp -pk omp "$SLURM_CPUS_PER_TASK" \
    "${LAMMPS_VARS[@]}" \
    -in "$LAMMPS_FILE" &
MPIRUN_PID=$!

TEARDOWN_GRACE=180   # seconds to allow a clean MPI_Finalize after completion
POLL_INTERVAL=15
while kill -0 "$MPIRUN_PID" 2>/dev/null; do
    if grep -q "Total wall time" "$LAMMPS_LOG" 2>/dev/null; then
        echo ">>> Detected 'Total wall time' in log — LAMMPS complete; all output written."
        WAITED=0
        while kill -0 "$MPIRUN_PID" 2>/dev/null && [ "$WAITED" -lt "$TEARDOWN_GRACE" ]; do
            sleep 10; WAITED=$((WAITED + 10))
        done
        if kill -0 "$MPIRUN_PID" 2>/dev/null; then
            echo ">>> mpirun still running ${TEARDOWN_GRACE}s after completion — assuming MPI/UCX teardown hang; killing it."
            kill -TERM "$MPIRUN_PID" 2>/dev/null; sleep 10
            kill -KILL "$MPIRUN_PID" 2>/dev/null; sleep 2
        fi
        break
    fi
    sleep "$POLL_INTERVAL"
done
wait "$MPIRUN_PID" 2>/dev/null
LAMMPS_RC=$?
set -e

if ! grep -q "Total wall time" "$LAMMPS_LOG" 2>/dev/null; then
    echo "LAMMPS did not reach 'Total wall time' in ${LAMMPS_LOG} (mpirun RC ${LAMMPS_RC})."
    echo "Treating as a genuine failure — no chaining, no post-processing."
    if [ -z "${LAMMPS_RC}" ] || [ "${LAMMPS_RC}" -eq 0 ]; then exit 1; fi
    exit "${LAMMPS_RC}"
fi
[ "$LAMMPS_RC" -ne 0 ] && echo "NOTE: mpirun returned ${LAMMPS_RC} after 'Total wall time' — output is complete; continuing."

# ── Chain again if this job, too, ran out of wall-clock ──────────────────────
if [ "${CHAIN:-0}" = "1" ]; then
    if grep -qE '^>>> PRODUCTION_STEPS_DONE [0-9]+ .*STATUS INCOMPLETE' "$LAMMPS_LOG"; then
        echo ">>> Production still incomplete — submitting the next continuation job..."
        ( cd "$SIM_FOLDER" && bash "$SELF" --chain "$SLURM_JOB_ID" ) \
        || echo ">>> WARNING: continuation submission FAILED — submit by hand from $SIM_FOLDER: continue_sim.sh --chain $SLURM_JOB_ID"
    else
        echo ">>> Production complete: chain finished. Final configuration in $CONT_DIR"
    fi
fi

# ── Post-processing (same pipeline as run_lammps.sh) ─────────────────────────
echo "======================================"
echo "Post-processing..."
echo "======================================"
if [ -n "$SWEEP_LEVEL" ]; then
    case "$FOLDER" in
        shear_slab)           export STRAINS="$SWEEP_LEVEL" ;;
        triaxial_compression) export COMPRESSIONS="$SWEEP_LEVEL" ;;
    esac
fi
export SKIP_WIDOM="${SKIP_WIDOM:-1}"
bash "$SCRIPT_DIR/postprocess.sh" \
    "$CONT_DIR" "$FOLDER" "$DATANAME" "$INTERACTION" "$TOTSTEPS" "$OLDSTEPS" "$PRESS_TARGET" \
    || echo "WARNING: postprocess.sh reported a failure — regenerate by hand (see scripts/postprocess.sh)."

echo "======================================"
echo "Done! Results in: $CONT_DIR"
echo "======================================"

# Never fall through the closing brace (see the parse guard above).
exit 0
}
