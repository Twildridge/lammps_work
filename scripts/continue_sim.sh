#!/bin/bash
# continue_sim.sh — Continue a LAMMPS simulation from a final .data file.
#
# Run this from inside the simulation folder (e.g. simulations/slab_with_flow/).
# The script finds the SLURM output file for the given job ID, reads the working
# directory and compression mode from it, locates the final .data file, and
# launches a continuation run — skipping all setup phases (cont=1).
#
# Usage (from inside simulations/slab_with_flow/ or simulations/slab_with_support/):
#   continue_sim.sh <job_id> <nsteps>
#
#   job_id  — the SLURM job ID, e.g. 49772594
#             (from the output file slab_flow.o49772594.exp-14-05)
#   nsteps  — additional timesteps to run
#
# Examples:
#   cd ~/Documents/lammps_work/simulations/slab_with_flow
#   continue_sim.sh 49772594 500000
#
#   cd ~/Documents/lammps_work/simulations/slab_with_support
#   continue_sim.sh 49800123 500000
#
# What cont=1 skips vs. runs:
#   slab_with_flow compression:
#     SKIPPED  — Phase 0 NVT, Phase 0.5 Langevin+NPH, compression piston drive,
#                epsilon=0 reference stress recording
#     RUNS     — all analysis computes/output fixes, then pure stress-relaxation
#                (halt_compress fires immediately → full run = relax_steps)
#   slab_with_flow permeation:
#     SKIPPED  — Phase 0 NVT
#     RUNS     — piston velocity re-applied at v_piston_perm, all observables,
#                halt_perm condition
#   slab_with_support:
#     SKIPPED  — soft push-off, minimize, gentle NVT thermalization, NPT warm-up
#     RUNS     — NPT production with all observables
#
# Output goes into a continuation subfolder inside the original run directory:
#   ~/Documents/lammps_runs/{original_run_dir}/continuation_{timestamp}/

set -e

# ── Arguments ─────────────────────────────────────────────────────────────────
if [ $# -ne 2 ]; then
    echo "Usage: continue_sim.sh <job_id> <nsteps>"
    echo ""
    echo "  job_id  — SLURM job ID (e.g. 49772594)"
    echo "  nsteps  — additional timesteps to run"
    echo ""
    echo "Run from inside the simulation folder:"
    echo "  cd ~/Documents/lammps_work/simulations/slab_with_flow"
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

if [ ! -f "$LAMMPS_FILE" ]; then
    echo "Error: no LAMMPS script found at $LAMMPS_FILE"
    echo "  Are you inside the correct simulation folder?"
    echo "  e.g. cd ~/Documents/lammps_work/simulations/slab_with_flow"
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

if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "Error: could not find (or access) working directory from output file."
    echo "  Looked for line: 'Working directory: /path/...'"
    echo "  Got: '$WORK_DIR'"
    echo "  Check that the job completed and the directory still exists."
    exit 1
fi

echo "Original work dir : $WORK_DIR"

# ── Find the final data file ──────────────────────────────────────────────────
case "$FOLDER" in
    slab_with_flow)
        DATA_FILE=$(ls "$WORK_DIR"/final_flow_*.data 2>/dev/null | head -1) ;;
    slab_with_support)
        DATA_FILE=$(ls "$WORK_DIR"/final_config_*.data 2>/dev/null | head -1) ;;
    *)
        echo "Error: unsupported folder '$FOLDER'. Supported: slab_with_flow, slab_with_support"
        exit 1 ;;
esac

if [ -z "$DATA_FILE" ] || [ ! -f "$DATA_FILE" ]; then
    echo "Error: no final .data file found in $WORK_DIR"
    echo "  Expected: final_flow_*.data or final_config_*.data"
    echo "  Did the LAMMPS run finish successfully?"
    exit 1
fi

echo "Final data file   : $DATA_FILE"

# ── Parse dataname, epsSS, epsSP from data file name ─────────────────────────
# Naming convention set by write_data at end of each LAMMPS script:
#   slab_with_flow:    final_flow_{dataname}_{epsSS}_{epsSP}_{nsteps}.data
#   slab_with_support: final_config_{dataname}_{epsSS}_{epsSP}_{nsteps}.data
#
# Split by '_' from the right:
#   index -1 = old nsteps  (integer)
#   index -2 = epsSP       (float, e.g. 1.0)
#   index -3 = epsSS       (float, e.g. 1.0)
#   rest     = dataname parts

BASENAME="$(basename "$DATA_FILE" .data)"
case "$FOLDER" in
    slab_with_flow)    REST="${BASENAME#final_flow_}" ;;
    slab_with_support) REST="${BASENAME#final_config_}" ;;
esac

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

# ── Auto-detect compression mode (slab_with_flow only) ───────────────────────
# run_lammps.sh runs LAMMPS, which prints:
#   >>> Mode: compression_mode=1  (0=permeation, 1=compression)
# Parse that line from the output file.
COMPRESSION_MODE=""
if [ "$FOLDER" = "slab_with_flow" ]; then
    COMPRESSION_MODE=$(grep ">>> Mode: compression_mode=" "$OUTPUT_FILE" \
                       | grep -o 'compression_mode=[0-9]' | cut -d= -f2 | head -1)
    if [ -z "$COMPRESSION_MODE" ]; then
        echo "Warning: could not detect compression_mode from output file."
        echo "  Falling back to value in LAMMPS script..."
        COMPRESSION_MODE=$(grep -E '^variable[[:space:]]+compression_mode[[:space:]]+equal' \
                           "$LAMMPS_FILE" | awk '{print $NF}')
    fi
    if [ -z "$COMPRESSION_MODE" ]; then
        echo "Error: could not determine compression_mode. Check $LAMMPS_FILE."
        exit 1
    fi
    MODE_LABEL=$([ "$COMPRESSION_MODE" = "0" ] && echo "permeation" || echo "compression")
    echo "Mode              : $COMPRESSION_MODE ($MODE_LABEL)"
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
    -var nsteps_eq   200000
    -var nsteps_prod 100000
    -var cont        1
)

if [ -n "$COMPRESSION_MODE" ]; then
    LAMMPS_VARS+=(-var compression_mode "$COMPRESSION_MODE")
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

if [ "$FOLDER" = "slab_with_flow" ]; then
    python "$SCRIPT_DIR/plot_lammps_log.py" "." "$STEM" --p-ext 1.8
else
    python "$SCRIPT_DIR/plot_lammps_log.py" "." "$STEM"
fi

WIDOM_TRAJ="${CONT_DIR}/traj_files/widom_${STEM}.lammpstrj"
if [ -f "$WIDOM_TRAJ" ]; then
    echo "Running cavity-biased Widom insertion..."
    if [ "$FOLDER" = "slab_with_flow" ]; then
        WIDOM_PEXT="1.8"; WIDOM_EXCL="2.0"; WIDOM_PISTON_EPS="0.0"
    else
        WIDOM_PEXT="1.5"; WIDOM_EXCL="";    WIDOM_PISTON_EPS="1.0"
    fi
    WIDOM_EXCL_ARGS=()
    [ -n "$WIDOM_EXCL" ] && WIDOM_EXCL_ARGS=(--exclusion-buffer "$WIDOM_EXCL")
    python "$SCRIPT_DIR/cavity_widom.py" \
        --traj        "$WIDOM_TRAJ" \
        --out-dir     "output_files/chemical_potential" \
        --out-stem    "$STEM" \
        --eps-sp      "$EPSSP" \
        --eps-ss      "$EPSSS" \
        --n-bins      40 \
        --n-trial     50000 \
        --r-cavity    0.5 \
        --temperature 1.0 \
        --p-ext       "$WIDOM_PEXT" \
        --piston-eps  "$WIDOM_PISTON_EPS" \
        "${WIDOM_EXCL_ARGS[@]}"
    python "$SCRIPT_DIR/plot_lammps_log.py" "." "$STEM" --p-ext "$WIDOM_PEXT"
else
    echo "WARNING: Widom trajectory not found — skipping cavity_widom.py"
fi

python "$SCRIPT_DIR/plot_stress_profiles.py"  "." "$STEM" "$OLDSTEPS"
python "$SCRIPT_DIR/plot_piston_data.py"      "." "$STEM" "$OLDSTEPS"
python "$SCRIPT_DIR/write_tracking.py"        "." "$STEM" ""

echo "======================================"
echo "Done! Results in: $CONT_DIR"
echo "======================================"
