#!/bin/bash
if [ $# -lt 4 ]; then
    echo "Usage: ./run_lammps.sh <folder_name> <dataname> <interaction> <nsteps> [oldsteps] [type]"
    echo "Example (fresh run):  ./run_lammps.sh slab_with_support slab_support_5beads_... 1.5_1.4 20000"
    echo "Example (continuation): ./run_lammps.sh slab_with_support slab_support_5beads_... 1.5_1.4 20000 20000"
    echo "Example (pure solvent P-sweep): ./run_lammps.sh pure_solvent pure_solvent_1000 1p0 0"
    echo "  interaction format: epsSS_epsSP (e.g., 1.5_0.4), or epsSS only for pure_solvent (e.g., 1p0)"
    echo "  oldsteps: total timesteps from previous run (defaults to 0 for fresh runs)"
    echo "  type: optional, 'stress' (adds 1), 'volume' (adds 2), or 'stressvol' (adds 3) to dataname"
    exit 1
fi

FOLDER=$1
DATANAME=$2
INTERACTION=$3
NSTEPS=$4
OLDSTEPS=${5:-0}  # Default to 0 for fresh runs
TOTSTEPS=$((OLDSTEPS + NSTEPS))

# P-sweep parameters (only used for pure_solvent; ignored by other scripts)
NSTEPS_EQ=200000    # equilibration steps per state point
NSTEPS_PROD=100000  # production/averaging steps per state point (must be divisible by 100)

# Scratch directory for trajectories
SCRATCH_DIR="/expanse/lustre/scratch/$USER/temp_project"

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
# Single timestamp captured once — both WORK_DIR and TRAJ_DIR use the same value
# so the symlink traj_files -> TRAJ_DIR is never stale.
RUN_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
WORK_DIR="$HOME/Documents/lammps_runs/${FOLDER}_${DATANAME}_${INTERACTION}_${RUN_TIMESTAMP}"
mkdir -p "$WORK_DIR"/{data_files,output_files/{stress_data,volume_data,piston_data,permeation_data,displacement_data,pair_data,chemical_potential},output_plots}

# Create trajectory directory in scratch and symlink to it
TRAJ_DIR="$SCRATCH_DIR/lammps_trajectories/${FOLDER}_${DATANAME}_${INTERACTION}_${RUN_TIMESTAMP}"
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
echo "  nsteps=$NSTEPS, oldsteps=$OLDSTEPS, totsteps=$TOTSTEPS"
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

# Pull latest scripts from GitHub before running (safe: skips if no network/conflict)
echo ">>> Syncing lammps_work from GitHub..."
git -C "$LAMMPS_WORK_DIR" pull --rebase --autostash || true

# CPU-only mode (NOT USING OMP FOR NOW (not on Expanse 2021 version)
echo "Running CPU-only with $SLURM_NTASKS tasks"
mpirun -n "${SLURM_NTASKS}" --bind-to "${OMPI_UNIT}" --map-by "node:pe=${OMP_NUM_THREADS}" \
    /home/dpollard/software/lammps/22Jul2025_update3/mpi-omp/gcc/10.2.0/openmpi/4.1.3/lammps-22Jul2025/build/lmp \
    -sf omp -pk omp $SLURM_CPUS_PER_TASK \
    -var dataname $DATANAME \
    -var interaction $INTERACTION \
    -var epsSS $EPSSS \
    -var epsSP $EPSSP \
    -var nsteps $NSTEPS \
    -var oldsteps $OLDSTEPS \
    -var totsteps $TOTSTEPS \
    -var nsteps_eq $NSTEPS_EQ \
    -var nsteps_prod $NSTEPS_PROD \
    -in $LAMMPS_FILE


# Determine suffix based on 6th argument (type) - moved from 7th position
SUFFIX=""
if [ $# -ge 6 ]; then
    case "$6" in
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

source /etc/profile.d/modules.sh
module unload python/3.8.12
module load anaconda3/2021.05/q4munrg
python -c "import numpy; print(numpy.__version__)"

echo "Generating convergence plot..."
python "$SCRIPT_DIR/plot_lammps_log.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}"

# ── Cavity-biased Widom insertion (auto-detected: only runs if dump exists) ──
# slab_with_support dumps an all-atom trajectory at nfreq_widom cadence for
# post-processing by cavity_widom.py.  Standard Widom (fix widom) returns zero
# inside the gel; cavity-biased Widom restricts insertions to void space and
# applies a bias correction, giving valid μ_ex estimates even in dense regions.
#
# Output written to output_files/chemical_potential/ alongside the fix widom
# and solvent density files, using the same dataname_interaction_totsteps stem:
#   mu_z_cavity_${DATANAME}_${INTERACTION}_${TOTSTEPS}.dat       (per-frame)
#   mu_z_cavity_summary_${DATANAME}_${INTERACTION}_${TOTSTEPS}.dat (time-averaged ± stderr)
WIDOM_TRAJ="${WORK_DIR}/traj_files/widom_${DATANAME}_${INTERACTION}_${TOTSTEPS}.lammpstrj"
echo "======================================"
echo "Cavity Widom check:"
echo "  Looking for: $WIDOM_TRAJ"
echo "  traj_files/ contents:"
ls -lh "${WORK_DIR}/traj_files/" 2>&1 | head -20
echo "======================================"
if [ -f "$WIDOM_TRAJ" ]; then
    echo "Running cavity-biased Widom insertion..."
    echo "  Trajectory: $WIDOM_TRAJ"
    echo "  Output dir: output_files/chemical_potential/"
    python "$SCRIPT_DIR/cavity_widom.py" \
        --traj      "$WIDOM_TRAJ" \
        --out-dir   "output_files/chemical_potential" \
        --out-stem  "${DATANAME}_${INTERACTION}_${TOTSTEPS}" \
        --eps-sp    "$EPSSP" \
        --eps-ss    "$EPSSS" \
        --n-bins    40 \
        --n-trial   20000 \
        --r-cavity  0.5 \
        --temperature 1.0

    echo "Re-generating convergence plot with cavity Widom panel..."
    python "$SCRIPT_DIR/plot_lammps_log.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}"
else
    echo "WARNING: Widom trajectory not found — skipping cavity_widom.py"
    echo "  Expected path: $WIDOM_TRAJ"
    echo "  TRAJ_DIR (scratch): $TRAJ_DIR"
    echo "  If the file is missing, check that:"
    echo "    1. The run completed without error"
    echo "    2. dump widom_traj is in slab_with_support.lmp (check git pulled correctly)"
    echo "    3. Scratch dir is accessible: ls $TRAJ_DIR"
fi

# Pure solvent P-sweep: run EOS plot instead of stress/piston/tracking scripts
if [ "$FOLDER" = "solvent_phase" ]; then
    echo "Solvent phase-sweep run detected — generating EOS plot..."
    python "$SCRIPT_DIR/plot_eos.py" "." "$DATANAME" "$INTERACTION"
    echo "======================================"
    echo "Done! Results are in: $WORK_DIR"
    echo "======================================"
    exit 0
fi

# Pure polymer P-sweep: run EOS plot instead of stress/piston/tracking scripts
if [ "$FOLDER" = "polymer_phase" ]; then
    echo "Polymer phase-sweep run detected — generating EOS plot..."
    python "$SCRIPT_DIR/plot_eos.py" "." "$DATANAME" "$INTERACTION"
    echo "======================================"
    echo "Done! Results are in: $WORK_DIR"
    echo "======================================"
    exit 0
fi

echo "Generating stress profiles..."
python "$SCRIPT_DIR/plot_stress_profiles.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" "$OLDSTEPS"

echo "Generating piston plots..."
python "$SCRIPT_DIR/plot_piston_data.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" "$OLDSTEPS"

echo "Generating computational efficiency plot..."
python "$SCRIPT_DIR/write_tracking.py" "." "${DATANAME}${SUFFIX}_${INTERACTION}_${TOTSTEPS}" "$SUFFIX"

echo "======================================"
echo "Done! Results are in: $WORK_DIR"
echo "======================================"