#!/bin/bash
# continue_sim.sh — Continue a LAMMPS simulation from a final .data file.
#
# Parses dataname / epsSS / epsSP from the output .data file's name, finds the
# original run's working directory, creates a continuation subfolder inside it,
# and runs LAMMPS with cont=1 — which skips all setup phases and goes directly
# to production measurement.
#
# What cont=1 skips vs. runs:
#   slab_with_flow compression:
#     SKIPPED  — Phase 0 NVT, Phase 0.5 Langevin+NPH, compression piston drive,
#                ε=0 reference recording
#     RUNS     — all analysis computes/chunks/output fixes, then stress-relaxation
#                run with piston frozen (halt fires immediately → run = relax_steps)
#   slab_with_flow permeation:
#     SKIPPED  — Phase 0 NVT
#     RUNS     — piston velocity re-applied at v_piston_perm, all observables, halt_perm
#   slab_with_support:
#     SKIPPED  — soft push-off, minimize, gentle NVT thermalization, NPT warm-up
#     RUNS     — NPT production with all observables
#
# Output goes into a continuation subfolder inside the original run's working dir:
#   ~/Documents/lammps_runs/{FOLDER}_{DATANAME}_{INTERACTION}_{ORIG_TIMESTAMP}/
#       continuation_{RUN_TIMESTAMP}/
#           output_files/...
#           traj_files -> scratch/...
#
# Usage:
#   ./continue_sim.sh <script_folder> <data_file> <nsteps> [compression_mode]
#
# Arguments:
#   script_folder    : slab_with_flow | slab_with_support
#   data_file        : path to the .data file from the previous run
#                      slab_with_flow    → final_flow_{dataname}_{epsSS}_{epsSP}_{nsteps}.data
#                      slab_with_support → final_config_{dataname}_{epsSS}_{epsSP}_{nsteps}.data
#   nsteps           : number of timesteps for this continuation run
#   compression_mode : (slab_with_flow ONLY) 0=permeation, 1=compression
#                      if omitted, auto-detected from the LAMMPS script file
#
# Examples:
#   ./continue_sim.sh slab_with_flow \
#       final_flow_slab_support_5beads_tall_rho04_1.0_1.0_500000.data \
#       500000 1
#
#   ./continue_sim.sh slab_with_support \
#       final_config_slab_support_5beads_tall_rho04_1.0_1.0_500000.data \
#       500000

set -e

# ── Argument parsing ──────────────────────────────────────────────────────────
if [ $# -lt 3 ]; then
    echo "Usage: $0 <script_folder> <data_file> <nsteps> [compression_mode]"
    echo ""
    echo "  script_folder    : slab_with_flow | slab_with_support"
    echo "  data_file        : path to final .data file from previous run"
    echo "  nsteps           : continuation timesteps"
    echo "  compression_mode : (slab_with_flow only) 0=permeation 1=compression"
    echo "                     auto-detected from script if omitted"
    echo ""
    echo "Examples:"
    echo "  $0 slab_with_flow final_flow_slab_support_5beads_tall_rho04_1.0_1.0_500000.data 500000 1"
    echo "  $0 slab_with_support final_config_slab_support_5beads_tall_rho04_1.0_1.0_500000.data 500000"
    exit 1
fi

FOLDER="$1"
DATA_FILE="$2"
NSTEPS="$3"
COMPRESSION_MODE_ARG="${4:-}"   # optional; slab_with_flow only

# ── Locate scripts and LAMMPS input file ─────────────────────────────────────
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
LAMMPS_WORK_DIR="$(dirname "$SCRIPT_DIR")"
SIM_DIR="$LAMMPS_WORK_DIR/simulations/$FOLDER"
LAMMPS_FILE="$SIM_DIR/${FOLDER}.lmp"

if [ ! -d "$SIM_DIR" ]; then
    echo "Error: simulation folder not found: $SIM_DIR"
    exit 1
fi
if [ ! -f "$LAMMPS_FILE" ]; then
    echo "Error: LAMMPS script not found: $LAMMPS_FILE"
    exit 1
fi
if [ ! -f "$DATA_FILE" ]; then
    echo "Error: data file not found: $DATA_FILE"
    exit 1
fi

DATA_FILE="$(realpath "$DATA_FILE")"

# ── Parse dataname, epsSS, epsSP from the data file name ─────────────────────
# Naming convention (set by write_data at the end of each LAMMPS script):
#   slab_with_flow:    final_flow_{dataname}_{epsSS}_{epsSP}_{nsteps}.data
#   slab_with_support: final_config_{dataname}_{epsSS}_{epsSP}_{nsteps}.data
#
# Strategy: strip the prefix, then split by '_' and read from the right:
#   index -1 = old nsteps  (integer)
#   index -2 = epsSP       (float, e.g. 1.0)
#   index -3 = epsSS       (float, e.g. 1.0)
#   index 0..-4 = dataname parts

BASENAME="$(basename "$DATA_FILE" .data)"

case "$FOLDER" in
    slab_with_flow)    PREFIX="final_flow_" ;;
    slab_with_support) PREFIX="final_config_" ;;
    *)
        echo "Error: unknown folder '$FOLDER'. Supported: slab_with_flow, slab_with_support"
        exit 1
        ;;
esac

REST="${BASENAME#$PREFIX}"
if [ "$REST" = "$BASENAME" ]; then
    echo "Warning: filename does not start with '$PREFIX' — attempting to parse anyway."
    REST="$BASENAME"
fi

IFS='_' read -ra PARTS <<< "$REST"
N=${#PARTS[@]}

if [ "$N" -lt 4 ]; then
    echo "Error: cannot parse dataname/epsSS/epsSP/nsteps from: $BASENAME"
    echo "  Expected: ${PREFIX}{dataname}_{epsSS}_{epsSP}_{oldnsteps}.data"
    exit 1
fi

OLD_NSTEPS="${PARTS[$((N-1))]}"
EPSSP="${PARTS[$((N-2))]}"
EPSSS="${PARTS[$((N-3))]}"
DATANAME=$(IFS='_'; echo "${PARTS[*]:0:$((N-3))}")
INTERACTION="${EPSSS}_${EPSSP}"

echo "======================================"
echo "Parsed from data file name:"
echo "  dataname    = $DATANAME"
echo "  epsSS       = $EPSSS"
echo "  epsSP       = $EPSSP"
echo "  interaction = $INTERACTION"
echo "  prev nsteps = $OLD_NSTEPS"
echo "======================================"

# ── Compression mode (slab_with_flow only) ───────────────────────────────────
COMPRESSION_MODE=""
if [ "$FOLDER" = "slab_with_flow" ]; then
    if [ -n "$COMPRESSION_MODE_ARG" ]; then
        COMPRESSION_MODE="$COMPRESSION_MODE_ARG"
        echo "compression_mode = $COMPRESSION_MODE (from argument)"
    else
        # Auto-detect from the hardcoded value in the LAMMPS script
        COMPRESSION_MODE=$(grep -E '^variable[[:space:]]+compression_mode[[:space:]]+equal' \
                           "$LAMMPS_FILE" | awk '{print $NF}')
        if [ -z "$COMPRESSION_MODE" ]; then
            echo "Error: could not auto-detect compression_mode from $LAMMPS_FILE"
            echo "  Pass it as the 4th argument: 0=permeation, 1=compression"
            exit 1
        fi
        echo "compression_mode = $COMPRESSION_MODE (auto-detected from script)"
    fi
    MODE_LABEL=$([ "$COMPRESSION_MODE" = "0" ] && echo "permeation" || echo "compression")
    echo "  → mode: $MODE_LABEL"
fi

# ── Working directory ─────────────────────────────────────────────────────────
# Fresh timestep counters: this is a new run from the data file, not cumulative.
# totsteps = nsteps so output filenames are self-consistent within this run.
OLDSTEPS=0
TOTSTEPS="$NSTEPS"

RUN_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SCRATCH_DIR="/expanse/lustre/scratch/$USER/temp_project"

# Find the most recent original run directory for this folder/dataname/interaction.
# Continuation outputs go as a subfolder inside it so all data from the same
# physical system stays together.
ORIG_RUNS_DIR="$HOME/Documents/lammps_runs"
ORIG_DIR=$(ls -dt "${ORIG_RUNS_DIR}/${FOLDER}_${DATANAME}_${INTERACTION}_"[0-9]* 2>/dev/null \
           | grep -v "continuation_" | head -1)

if [ -n "$ORIG_DIR" ]; then
    WORK_DIR="${ORIG_DIR}/continuation_${RUN_TIMESTAMP}"
    # Mirror the scratch traj path under the original run's scratch dir
    ORIG_TRAJ_BASE=$(basename "$ORIG_DIR")
    TRAJ_DIR="$SCRATCH_DIR/lammps_trajectories/${ORIG_TRAJ_BASE}/continuation_${RUN_TIMESTAMP}"
    echo "Original run dir  : $ORIG_DIR"
    echo "Continuation dir  : $WORK_DIR"
else
    # Fall back to a standalone directory if no matching original run is found
    WORK_DIR="${ORIG_RUNS_DIR}/${FOLDER}_${DATANAME}_${INTERACTION}_cont_${RUN_TIMESTAMP}"
    TRAJ_DIR="$SCRATCH_DIR/lammps_trajectories/${FOLDER}_${DATANAME}_${INTERACTION}_cont_${RUN_TIMESTAMP}"
    echo "No original run dir found — using standalone directory"
    echo "Working directory : $WORK_DIR"
fi

mkdir -p "$WORK_DIR"/{data_files,output_files/{stress_data,volume_data,piston_data,permeation_data,displacement_data,pair_data,chemical_potential},output_plots}
mkdir -p "$TRAJ_DIR"
ln -s "$TRAJ_DIR" "$WORK_DIR/traj_files"

echo "Traj directory    : $TRAJ_DIR"
echo "======================================"

# Symlink data file into the run directory under the expected name
ln -s "$DATA_FILE" "$WORK_DIR/data_files/${DATANAME}.data"

cd "$WORK_DIR" || exit 1

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
    -var cont        1          # skip setup phases; go directly to production
)

# -var compression_mode overrides the hardcoded script value (slab_with_flow only)
if [ -n "$COMPRESSION_MODE" ]; then
    LAMMPS_VARS+=(-var compression_mode "$COMPRESSION_MODE")
fi

# ── Print run summary ─────────────────────────────────────────────────────────
echo "Running LAMMPS continuation:"
echo "  folder      = $FOLDER"
echo "  dataname    = $DATANAME"
echo "  interaction = $INTERACTION  (epsSS=$EPSSS epsSP=$EPSSP)"
echo "  nsteps      = $NSTEPS  (oldsteps=$OLDSTEPS totsteps=$TOTSTEPS)"
[ -n "$COMPRESSION_MODE" ] && echo "  mode        = $COMPRESSION_MODE ($MODE_LABEL)"
echo "  LAMMPS file = $LAMMPS_FILE"
echo ""

# ── Run LAMMPS ────────────────────────────────────────────────────────────────
echo "SLURM tasks: $SLURM_NTASKS  CPUs/task: $SLURM_CPUS_PER_TASK"

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

# Cavity-biased Widom insertion
WIDOM_TRAJ="${WORK_DIR}/traj_files/widom_${STEM}.lammpstrj"
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
    echo "  Expected: $WIDOM_TRAJ"
fi

python "$SCRIPT_DIR/plot_stress_profiles.py"  "." "$STEM" "$OLDSTEPS"
python "$SCRIPT_DIR/plot_piston_data.py"      "." "$STEM" "$OLDSTEPS"
python "$SCRIPT_DIR/write_tracking.py"        "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" ""

echo "======================================"
echo "Done! Results in: $WORK_DIR"
echo "======================================"
