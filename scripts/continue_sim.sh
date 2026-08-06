#!/bin/bash
# continue_sim.sh — Continue a LAMMPS simulation from a final .data file.
#
# Run this from inside the simulation folder. The script finds the SLURM
# output file for the given job ID, reads the working directory from it,
# locates the final .data file, and launches a continuation run — skipping
# setup phases via cont=1. This is a REAL restart (skip setup, keep going),
# not a fresh resubmission — see README §5e for how this differs from
# editing a .batch file's NSTEPS and resubmitting.
#
# Usage (from inside the simulation folder):
#   continue_sim.sh <job_id> <nsteps>
#
#   job_id  — the SLURM job ID, e.g. 49772594
#             (from the output file <folder>.o49772594.exp-14-05)
#   nsteps  — additional timesteps to run
#
# Examples:
#   cd ~/Documents/lammps_work/simulations/slab_with_support
#   continue_sim.sh 49800123 500000
#
#   cd ~/Documents/lammps_work/simulations/triaxial_compression
#   continue_sim.sh 49900456 2000000
#
# Supported folders and what cont=1 means for each:
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
# Output goes into a continuation subfolder inside the original run directory:
#   ~/Documents/lammps_runs/{folder}/{original_run_dir}/continuation_{timestamp}/

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

# ── Self-submit via SLURM if run from a login node ────────────────────────────
if [ -z "$SLURM_JOB_ID" ]; then
    # Locate the simulation's .batch file to borrow SLURM resource settings
    SELF="$(realpath "${BASH_SOURCE[0]}")"
    BATCH_FILE="$(pwd)/$(basename "$(pwd)").batch"
    if [ ! -f "$BATCH_FILE" ]; then
        echo "Error: not inside a SLURM job and no .batch file found at $BATCH_FILE"
        echo "  Either run from inside the simulation folder, or submit via sbatch manually."
        exit 1
    fi

    # Extract #SBATCH directives (skip --output so we set our own)
    SBATCH_LINES=$(grep '^#SBATCH' "$BATCH_FILE" | grep -v -- '--output' | grep -v -- '--error')

    # Build module-load lines from the .batch file
    MODULE_LINES=$(grep '^module' "$BATCH_FILE" || true)
    ENV_LINES=$(grep '^declare' "$BATCH_FILE" || true)

    TMPBATCH=$(mktemp /tmp/cont_XXXXXX.sh)
    cat > "$TMPBATCH" << BATCHEOF
#!/bin/bash
${SBATCH_LINES}
#SBATCH --job-name=cont_$1
#SBATCH --output=%x.o%j.%N

${ENV_LINES}

${MODULE_LINES}

cd "$(pwd)"
"$SELF" "$1" "$2"
BATCHEOF

    echo "Submitting continuation as SLURM batch job..."
    sbatch "$TMPBATCH"
    rm -f "$TMPBATCH"
    exit 0
fi

# ── Arguments ─────────────────────────────────────────────────────────────────
if [ $# -ne 2 ]; then
    echo "Usage: continue_sim.sh <job_id> <nsteps>"
    echo ""
    echo "  job_id  — SLURM job ID (e.g. 49772594)"
    echo "  nsteps  — additional timesteps to run"
    echo ""
    echo "Run from inside the simulation folder:"
    echo "  cd ~/Documents/lammps_work/simulations/slab_with_support"
    echo "  continue_sim.sh 49772594 500000"
    exit 1
fi

JOB_ID="$1"
NSTEPS="$2"

# ── Locate scripts and LAMMPS file ────────────────────────────────────────────
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
LAMMPS_WORK_DIR="$(dirname "$SCRIPT_DIR")"

# Folder is the name of the directory we're running from
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

# ── Find the SLURM output file ────────────────────────────────────────────────
OUTPUT_FILE=$(ls *.o${JOB_ID}.* 2>/dev/null | head -1)

if [ -z "$OUTPUT_FILE" ]; then
    echo "Error: no SLURM output file matching *.o${JOB_ID}.* found in $PWD"
    echo "  Make sure you are in the folder where you ran sbatch."
    exit 1
fi

echo "======================================"
echo "SLURM output file : $OUTPUT_FILE"

# ── Parse working directory from output file ──────────────────────────────────
# run_lammps.sh prints: "Working directory: /path/to/work_dir"
WORK_DIR=$(grep "^Working directory:" "$OUTPUT_FILE" | awk '{print $NF}' | head -1)

if [ -z "$WORK_DIR" ]; then
    echo "Error: could not find working directory line in output file."
    echo "  Looked for line: 'Working directory: /path/...'"
    echo "  Check that the job completed and the output file is intact."
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
        echo "  Check that the job completed and the directory still exists."
        exit 1
    fi
fi

echo "Original work dir : $WORK_DIR"

# ── Find the final data file ──────────────────────────────────────────────────
PREFIX="${FOLDER_PREFIX[$FOLDER]}"
DATA_FILE=$(ls "$WORK_DIR"/${PREFIX}_*.data 2>/dev/null | head -1)

if [ -z "$DATA_FILE" ] || [ ! -f "$DATA_FILE" ]; then
    echo "Error: no final .data file found in $WORK_DIR"
    echo "  Expected: ${PREFIX}_*.data"
    echo "  Did the LAMMPS run finish successfully?"
    exit 1
fi

echo "Final data file   : $DATA_FILE"

# ── Parse dataname, epsSS, epsSP from data file name ─────────────────────────
# Naming convention set by write_data at the end of each LAMMPS script:
#   ${prefix}_{dataname}_{epsSS}_{epsSP}_{nsteps}.data
# Split by '_' from the right:
#   index -1 = old nsteps  (integer)
#   index -2 = epsSP       (float, e.g. 1.0)
#   index -3 = epsSS       (float, e.g. 1.0)
#   rest     = dataname parts

BASENAME="$(basename "$DATA_FILE" .data)"
REST="${BASENAME#${PREFIX}_}"

IFS='_' read -ra PARTS <<< "$REST"
N=${#PARTS[@]}
if [ "$N" -lt 4 ]; then
    echo "Error: cannot parse dataname/epsSS/epsSP/nsteps from: $BASENAME"
    exit 1
fi

OLD_NSTEPS="${PARTS[$((N-1))]}"
EPSSP="${PARTS[$((N-2))]}"
EPSSS="${PARTS[$((N-3))]}"
DATANAME=$(IFS='_'; echo "${PARTS[*]:0:$((N-3))}")
INTERACTION="${EPSSS}_${EPSSP}"

echo "Data name         : $DATANAME"
echo "Interaction       : $INTERACTION  (epsSS=$EPSSS  epsSP=$EPSSP)"
echo "Prev nsteps       : $OLD_NSTEPS"

# ── Sweep folders: auto-detect the last level reached ─────────────────────────
# For triaxial_compression/shear_slab, continuing must extend the equilibration
# hold at whatever level the original run last completed — never re-drive
# through the whole strain/compression ladder. Scan the original run's
# output_files/stress_data/ for the highest _c<level>/_g<level> tag present.
SWEEP_TAG="${FOLDER_SWEEP_TAG[$FOLDER]:-}"
SWEEP_LEVEL=""
if [ -n "$SWEEP_TAG" ]; then
    SWEEP_LEVEL=$(ls "$WORK_DIR"/output_files/stress_data/*"${SWEEP_TAG}"[0-9]*.dat 2>/dev/null \
        | sed -E "s/.*${SWEEP_TAG}([0-9.]+)\.dat\$/\1/" \
        | sort -g | uniq | tail -1)
    if [ -z "$SWEEP_LEVEL" ]; then
        echo "Error: could not auto-detect the last ${SWEEP_TAG}<level> reached in"
        echo "  $WORK_DIR/output_files/stress_data/"
        exit 1
    fi
    echo "Last level reached : ${SWEEP_TAG}${SWEEP_LEVEL}"
fi
echo "======================================"

# ── Set up continuation directory inside original work dir ───────────────────
OLDSTEPS=0
TOTSTEPS="$NSTEPS"
RUN_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SCRATCH_DIR="/expanse/lustre/scratch/$USER/temp_project"

CONT_DIR="${WORK_DIR}/continuation_${RUN_TIMESTAMP}"
ORIG_TRAJ_BASE=$(basename "$WORK_DIR")
TRAJ_DIR="$SCRATCH_DIR/lammps_trajectories/${ORIG_TRAJ_BASE}/continuation_${RUN_TIMESTAMP}"

mkdir -p "$CONT_DIR"/{data_files,output_files/{stress_data,volume_data,piston_data,permeation_data,displacement_data,pair_data,chemical_potential},output_plots}
mkdir -p "$TRAJ_DIR"
ln -s "$TRAJ_DIR" "$CONT_DIR/traj_files"

# Symlink the data file under the name LAMMPS expects
ln -s "$(realpath "$DATA_FILE")" "$CONT_DIR/data_files/${DATANAME}.data"

echo "Continuation dir  : $CONT_DIR"
echo "Traj dir (scratch): $TRAJ_DIR"
echo "======================================"

cd "$CONT_DIR" || exit 1

# ── Build LAMMPS variable list ────────────────────────────────────────────────
LAMMPS_VARS=(
    -var dataname    "$DATANAME"
    -var interaction "$INTERACTION"
    -var epsSS       "$EPSSS"
    -var epsSP       "$EPSSP"
    -var nsteps      "$NSTEPS"
    -var oldsteps    "$OLDSTEPS"
    -var totsteps    "$TOTSTEPS"
    -var cont        1
)
if [ -n "$SWEEP_TAG" ]; then
    LAMMPS_VARS+=(-var "${FOLDER_SWEEP_VAR[$FOLDER]}" "$SWEEP_LEVEL")
fi

# ── Run LAMMPS ────────────────────────────────────────────────────────────────
echo "Running LAMMPS continuation ($SLURM_NTASKS tasks)..."

mpirun -n "${SLURM_NTASKS}" \
    --bind-to "${OMPI_UNIT}" \
    --map-by "node:pe=${OMP_NUM_THREADS}" \
    /home/dpollard/software/lammps/22Jul2025_update3/mpi-omp/gcc/10.2.0/openmpi/4.1.3/lammps-22Jul2025/build/lmp \
    -sf omp -pk omp "$SLURM_CPUS_PER_TASK" \
    "${LAMMPS_VARS[@]}" \
    -in "$LAMMPS_FILE"

# ── Post-processing ───────────────────────────────────────────────────────────
echo "======================================"
echo "Post-processing..."
echo "======================================"

source /etc/profile.d/modules.sh
module unload python/3.8.12
module load anaconda3/2021.05/q4munrg

STEM="${DATANAME}_${INTERACTION}_${TOTSTEPS}"

python "$SCRIPT_DIR/plot_lammps_log.py" "." "$STEM"

WIDOM_MIN_STEPS=500000
WIDOM_TRAJ="${CONT_DIR}/traj_files/widom_${STEM}.lammpstrj"
if [ ! -f "$WIDOM_TRAJ" ]; then
    echo "WARNING: Widom trajectory not found — skipping cavity_widom.py"
elif [ "$NSTEPS" -lt "$WIDOM_MIN_STEPS" ]; then
    echo "Skipping cavity_widom.py — only ${NSTEPS} steps (need >=${WIDOM_MIN_STEPS} for decorrelated frames)."
else
    echo "Running cavity-biased Widom insertion..."
    python -u "$SCRIPT_DIR/cavity_widom.py" \
        --traj        "$WIDOM_TRAJ" \
        --out-dir     "output_files/chemical_potential" \
        --out-stem    "$STEM" \
        --eps-sp      "$EPSSP" \
        --eps-ss      "$EPSSS" \
        --n-bins      40 \
        --n-trial     50000 \
        --r-cavity    0.5 \
        --temperature 1.0 \
        --p-ext       "1.5" \
        --piston-eps  "1.0"
    python "$SCRIPT_DIR/plot_lammps_log.py" "." "$STEM" --p-ext "1.5"
fi

python "$SCRIPT_DIR/plot_stress_profiles.py"  "." "$STEM" "$OLDSTEPS"
python "$SCRIPT_DIR/plot_piston_data.py"      "." "$STEM" "$OLDSTEPS"

echo "======================================"
echo "Done! Results in: $CONT_DIR"
echo "======================================"
