#!/bin/bash
if [ $# -lt 4 ]; then
    echo "Usage: ./run_lammps.sh <folder_name> <dataname> <interaction> <nsteps> [oldsteps] [datasteps] [type]"
    echo "Example: ./run_lammps.sh slab_with_support slab_support_10beads_10x10x5_rho6_extra_padding3 1.5_1.6 40000 40000 40000 stress"
    echo "  interaction format: epsSS_epsSP (e.g., 1.5_0.4)"
    echo "  oldsteps: previous timesteps for restart (defaults to nsteps)"
    echo "  datasteps: number of timesteps in the data file (defaults to nsteps)"
    echo "  type: optional, 'stress' (adds 1), 'volume' (adds 2), or 'stressvol' (adds 3) to dataname"
    exit 1
fi

FOLDER=$1
DATANAME=$2
INTERACTION=$3
NSTEPS=$4
OLDSTEPS=${5:-$NSTEPS}
DATASTEPS=${6:-$NSTEPS}
TOTSTEPS=$((OLDSTEPS + NSTEPS))

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

# Create a working directory for this run
WORK_DIR="$HOME/Documents/lammps_runs/${FOLDER}_${DATANAME}_${INTERACTION}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$WORK_DIR"/{data_files,output_files/{stress_data,volume_data},output_plots,traj_files}

echo "======================================"
echo "Working directory: $WORK_DIR"
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
echo "  datasteps=$DATASTEPS"
echo "SLURM tasks per node: $SLURM_NTASKS_PER_NODE"
echo "SLURM CPUs per task: $SLURM_CPUS_PER_TASK"

# Run LAMMPS — has these installed packages:
# ASPHERE COLVARS DIELECTRIC DIPOLE DRUDE EFF EXTRA-FIX EXTRA-PAIR FEP GRANULAR 
# INTERLAYER KOKKOS KSPACE MACHDYN MANYBODY MC MEAM MISC ML-SNAP MOLECULE OPENMP 
# OPT PHONON PYTHON QEQ REAXFF REPLICA RIGID

mpirun -n $SLURM_NTASKS \
    /opt/packages/LAMMPS/lammps-22Jul2025/build-RM-gcc13.3.1/lmp \
    -sf omp -pk omp $SLURM_CPUS_PER_TASK  \
    -var dataname $DATANAME \
    -var interaction $INTERACTION \
    -var epsSS $EPSSS \
    -var epsSP $EPSSP \
    -var nsteps $NSTEPS \
    -var oldsteps $OLDSTEPS \
    -var totsteps $TOTSTEPS \
    -var datasteps $DATASTEPS \
    -in $LAMMPS_FILE

# Determine suffix based on 7th argument (type)
SUFFIX=""
if [ $# -ge 7 ]; then
    case "$7" in
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

module load anaconda3/2024.10-1

echo "Generating convergence plot..."
python "$SCRIPT_DIR/plot_lammps_log.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}"

echo "Generating stress profiles..."
python "$SCRIPT_DIR/plot_stress_profiles.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" "$DATASTEPS"

echo "Generating computational efficiency plot..."
python "$SCRIPT_DIR/write_tracking.py" "." "${DATANAME}${SUFFIX}_${INTERACTION}_${TOTSTEPS}" "$SUFFIX"

echo "======================================"
echo "Done! Results are in: $WORK_DIR"
echo "======================================"
