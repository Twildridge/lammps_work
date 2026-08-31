#!/bin/bash
if [ $# -lt 4 ]; then
    echo "Usage: ./run_lammps.sh <folder_name> <dataname> <interaction> <nsteps> [type]"
    echo "Example (fresh run):  ./run_lammps.sh slab_with_support slab_support_5beads_... 1.5_1.4 20000"
    echo "  interaction format: epsSS_epsSP (e.g., 1.5_0.4)"
    echo "  type: optional, 'stress' (adds 1), 'volume' (adds 2), or 'stressvol' (adds 3) to dataname"
    echo "  to continue a finished run, use continue_sim.sh instead of resubmitting this script"
    exit 1
fi

FOLDER=$1
DATANAME=$2
INTERACTION=$3
NSTEPS=$4
TOTSTEPS=$NSTEPS

# COMPRESS_STAGES: space-separated CUMULATIVE volumetric-strain ladder
# (compress_slab only) — one value per compression stage, passed through to
# LAMMPS as an index variable. Set in compress_slab_pod.batch via
# COMPRESSION_1, COMPRESSION_2, ... -> STAGE_TARGETS -> COMPRESS_STAGES. The
# number of stages is just the length of this list. Default reproduces
# compress_slab.lmp's original fixed 3-stage ladder. Same stale-batch caveat
# as elsewhere: the batch must self-sync BEFORE exporting this.
if [ -z "${COMPRESS_STAGES:-}" ]; then
    echo ">>> NOTE: COMPRESS_STAGES unset — compress_slab falls back to the default 0.015 0.030 0.045 ladder."
fi
COMPRESS_STAGES=${COMPRESS_STAGES:-"0.015 0.030 0.045"}

# Scratch directory for trajectories
SCRATCH_DIR="/scratch/$USER"

# Get the directory where this script lives (should be lammps_work/scripts/)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
LAMMPS_WORK_DIR="$(dirname "$SCRIPT_DIR")"  # Parent directory (lammps_work/)

# Parse interaction into epsSS and epsSP
IFS='_' read -r EPSSS EPSSP <<< "$INTERACTION"

# Check if simulation folder exists
SIM_DIR="$LAMMPS_WORK_DIR/simulations/$FOLDER"
if [ ! -d "$SIM_DIR" ]; then
    echo "Error: Simulation folder $SIM_DIR not found"
    exit 1
fi

# Check if LAMMPS file exists
LAMMPS_FILE="$SIM_DIR/${FOLDER}.lmp"
if [ ! -f "$LAMMPS_FILE" ]; then
    echo "Error: $LAMMPS_FILE not found"
    exit 1
fi

# Create a working directory for this run in home (for small files)
WORK_DIR="$HOME/Documents/lammps_runs/${FOLDER}_${DATANAME}_${INTERACTION}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$WORK_DIR"/{data_files,output_files/{stress_data,volume_data,piston_data,permeation_data,displacement_data,pair_data},output_plots}

# Create trajectory directory in scratch and symlink to it
TRAJ_DIR="$SCRATCH_DIR/lammps_trajectories/${FOLDER}_${DATANAME}_${INTERACTION}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TRAJ_DIR"
ln -s "$TRAJ_DIR" "$WORK_DIR/traj_files"

echo "======================================"
echo "Working directory: $WORK_DIR"
echo "Trajectory directory (scratch): $TRAJ_DIR"
echo "======================================"

# Copy or link the data file
DATA_FILE_SOURCE="$HOME/Documents/lammps_data/input_data/${DATANAME}.data"
if [ ! -f "$DATA_FILE_SOURCE" ]; then
    echo "Error: Data file $DATA_FILE_SOURCE not found"
    echo "Please ensure your .data file is in ~/Documents/lammps_data/input_data/"
    exit 1
fi

# Create symlink to data file
ln -s "$DATA_FILE_SOURCE" "$WORK_DIR/data_files/${DATANAME}.data"

# Change to working directory
cd "$WORK_DIR" || exit 1

echo "Running LAMMPS in $FOLDER with:"
echo "  dataname=$DATANAME"
echo "  epsSS=$EPSSS, epsSP=$EPSSP"
echo "  nsteps=$NSTEPS, totsteps=$TOTSTEPS"
echo "SLURM tasks per node: $SLURM_NTASKS_PER_NODE"
echo "SLURM CPUs per task: $SLURM_CPUS_PER_TASK"
echo "DEBUG: SLURM_NTASKS_PER_NODE: $SLURM_NTASKS_PER_NODE"
echo "       SLURM_NTASKS: $SLURM_NTASKS"
echo "       SLURM_NNODES: $SLURM_NNODES"

# Run LAMMPS-22Jul2025 — has these installed packages:
# ASPHERE COLVARS DIELECTRIC DIPOLE DRUDE EFF EXTRA-FIX EXTRA-PAIR FEP GRANULAR 
# INTERLAYER KOKKOS KSPACE MACHDYN MANYBODY MC MEAM MISC ML-SNAP MOLECULE OPENMP 
# OPT PHONON PYTHON QEQ REAXFF REPLICA RIGID

# check within lammps build directory with: 
# grep "PKG_.*:BOOL=\(yes\|ON\)$" CMakeCache.txt | sed 's/PKG_//; s/:BOOL.*//'

# Check if GPUs are allocated
NGPUS=${SLURM_GPUS_ON_NODE:-0}
echo "SLURM_GPUS_ON_NODE: $SLURM_GPUS_ON_NODE"
echo "GPUs allocated: $SLURM_GPUS"

# ── Pull latest scripts from GitHub before running ────────────────────────────
# Serialised, fast-forward-only, best-effort. See scripts/node_git_sync.sh.
source "$SCRIPT_DIR/node_git_sync.sh"
echo ">>> Syncing lammps_work from GitHub..."
sync_lammps_work || true

# CPU-only mode (NOT USING OMP FOR NOW (not on Expanse 2021 version)
echo "Running CPU-only with $SLURM_NTASKS tasks"
mpirun -n "${SLURM_NTASKS}" --bind-to "${OMPI_UNIT}" --map-by "node:pe=${OMP_NUM_THREADS}" \
    ~/MD/lammps/lammps-22Jul2025/build/lmp \
    -sf omp -pk omp $SLURM_CPUS_PER_TASK \
    -var dataname $DATANAME \
    -var interaction $INTERACTION \
    -var epsSS $EPSSS \
    -var epsSP $EPSSP \
    -var nsteps $NSTEPS \
    -var oldsteps 0 \
    -var totsteps $TOTSTEPS \
    -var stage_targets $COMPRESS_STAGES \
    -var stage_targets_list "$COMPRESS_STAGES" \
    -in $LAMMPS_FILE


# Determine suffix based on 5th argument (type) - moved from 6th position after OLDSTEPS removal
SUFFIX=""
if [ $# -ge 5 ]; then
    case "$5" in
        stress)
            SUFFIX="1"
            ;;
        volume)
            SUFFIX="2"
            ;;
        stressvol)
            SUFFIX="3"
            ;;
    esac
fi

# Run post-processing Python scripts
echo "======================================"
echo "Running post-processing..."
echo "======================================"

cd "$WORK_DIR" || exit 1

module load miniconda
source $(conda info --base)/etc/profile.d/conda.sh # loads Conda's shell functions into bash session
conda activate lammps_analysis

echo "Generating convergence plot..."
python "$SCRIPT_DIR/plot_lammps_log.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}"

echo "Generating stress profiles..."
python "$SCRIPT_DIR/plot_stress_profiles.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0

echo "Generating piston plots..."
python "$SCRIPT_DIR/plot_piston_data.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0


echo "======================================"
echo "Done! Results are in: $WORK_DIR"
echo "======================================"