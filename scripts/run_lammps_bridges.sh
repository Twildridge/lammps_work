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
STRAINS=${STRAINS:-0.1}  # Space-separated shear-strain list (shear_slab only); LAMMPS index var
# COMPRESSIONS: space-separated cumulative strain targets, a LAMMPS index variable
# (compress_slab: VOLUMETRIC strain).  Default 0.1 = one level.  Same stale-batch caveat
# as elsewhere: the batch must self-sync BEFORE exporting this.
if [ -z "${COMPRESSIONS:-}" ]; then
    echo ">>> NOTE: COMPRESSIONS unset -- the compression decks fall back to the single level 0.1."
fi
COMPRESSIONS=${COMPRESSIONS:-0.1}
PRESS_TARGET=${PRESS_TARGET:-1.5}              # bath pressure on the wet pistons
PISTON_MASS=${PISTON_MASS:-1000}               # mass of every wet-piston bead
C_PIST_FRAC=${C_PIST_FRAC:-1.0}                # piston damping as a fraction of the critical value
NPT_PISTON_STEPS=${NPT_PISTON_STEPS:-1000000}  # Phase-1 NPT-piston settle length
SETTLE_HALT=${SETTLE_HALT:-0}                  # 1 = end the settle early once both pistons are at rest
DRIVE_SPLIT=${DRIVE_SPLIT:-0.5}                # hi-plate share of each gap closure (0.5 = symmetric)
PISTON_TRANSPARENT=${PISTON_TRANSPARENT:-0}    # slab_with_support: 1 = piston sheet invisible to solvent
# compress_slab knobs (2026-09-29; the deck is the bulk-modulus counterpart of
# triaxial_compression_two_pist and reads the SAME variables: COMPRESSIONS is its list of
# cumulative VOLUMETRIC strain targets, and PRESS_TARGET / PISTON_MASS / C_PIST_FRAC /
# NPT_PISTON_STEPS / SETTLE_HALT / DRIVE_SPLIT mean what they mean there).  Two more:
#   SEAT_MODE  rg (default) | bb : the gel face the six plates seat on (compress_slab.lmp)
#   HOLD_AUTO  1 (default) = size every hold from tau_1 = L^2/(pi^2 D_c); 0 = NSTEPS flat
# BARO_MODE (slab_with_support only): aniso (default) | iso | z -- barostat coupling; use
#   iso for an ISOLATED gel (see slab_with_support.lmp).
SEAT_MODE=${SEAT_MODE:-rg}
HOLD_AUTO=${HOLD_AUTO:-1}
BARO_MODE=${BARO_MODE:-aniso}
PIN_GEL=${PIN_GEL:-0}                # slab_with_support: 1 = zero the gel's linear + angular momentum (isolated gel)

# Scratch directory for trajectories
SCRATCH_DIR="/ocean/projects/chm250028p/$USER"

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
mkdir -p "$WORK_DIR"/{data_files,output_files/{stress_data,volume_data,piston_data,permeation_data,displacement_data,pair_data,chemical_potential},output_plots}

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

if [ $NGPUS -gt 0 ]; then
    # GPU mode with Kokkos
    echo "Running with $NGPUS GPU(s) and $SLURM_CPUS_PER_TASK threads per GPU"
    mpirun -n $SLURM_NTASKS \
        /opt/packages/LAMMPS/lammps-22Jul2025/build-V100-gcc13.3.1/lmp \
        -k on g $NGPUS t $SLURM_CPUS_PER_TASK -sf kk -pk kokkos newton on neigh half comm device \
        -var dataname $DATANAME \
        -var interaction $INTERACTION \
        -var epsSS $EPSSS \
        -var epsSP $EPSSP \
        -var nsteps $NSTEPS \
        -var oldsteps 0 \
        -var totsteps $TOTSTEPS \
        -var strains $STRAINS \
        -var strains_list "$STRAINS" \
        -var compressions $COMPRESSIONS \
        -var compressions_list "$COMPRESSIONS" \
        -var press_target "$PRESS_TARGET" \
        -var piston_mass "$PISTON_MASS" \
        -var c_pist_frac "$C_PIST_FRAC" \
        -var npt_piston_steps "$NPT_PISTON_STEPS" \
        -var settle_halt "$SETTLE_HALT" \
        -var drive_split "$DRIVE_SPLIT" \
        -var piston_transparent "$PISTON_TRANSPARENT" \
        -var seat_mode "$SEAT_MODE" \
        -var hold_auto "$HOLD_AUTO" \
        -var baro_mode "$BARO_MODE" \
        -var pin_gel "$PIN_GEL" \
        -in $LAMMPS_FILE
else
    # CPU-only mode
    echo "Running CPU-only with $SLURM_NTASKS tasks"
    mpirun -n $SLURM_NTASKS \
        /opt/packages/LAMMPS/lammps-22Jul2025/build-RM-gcc13.3.1/lmp \
        -sf omp -pk omp $SLURM_CPUS_PER_TASK  \
        -var dataname $DATANAME \
        -var interaction $INTERACTION \
        -var epsSS $EPSSS \
        -var epsSP $EPSSP \
        -var nsteps $NSTEPS \
        -var oldsteps 0 \
        -var totsteps $TOTSTEPS \
        -var strains $STRAINS \
        -var strains_list "$STRAINS" \
        -var compressions $COMPRESSIONS \
        -var compressions_list "$COMPRESSIONS" \
        -var press_target "$PRESS_TARGET" \
        -var piston_mass "$PISTON_MASS" \
        -var c_pist_frac "$C_PIST_FRAC" \
        -var npt_piston_steps "$NPT_PISTON_STEPS" \
        -var settle_halt "$SETTLE_HALT" \
        -var drive_split "$DRIVE_SPLIT" \
        -var piston_transparent "$PISTON_TRANSPARENT" \
        -var seat_mode "$SEAT_MODE" \
        -var hold_auto "$HOLD_AUTO" \
        -var baro_mode "$BARO_MODE" \
        -var pin_gel "$PIN_GEL" \
        -in $LAMMPS_FILE
fi



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

module load anaconda3/2024.10-1

echo "Generating convergence plot..."
python "$SCRIPT_DIR/lib/plot_lammps_log.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}"

if [ "$FOLDER" = "compress_slab" ]; then
    # compress_slab writes one file set per level, tagged _<hold>_c<level>, in the layout of
    # triaxial_compression_two_pist: the consolidated sweep plotter overlays the levels.
    echo "Generating consolidated compression sweep plots (levels: $COMPRESSIONS)..."
    python "$SCRIPT_DIR/lib/plot_compression_strain_sweep.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" "$COMPRESSIONS" 0 \
        || echo "  WARNING: plot_compression_strain_sweep.py failed (skipping)"
    python "$SCRIPT_DIR/lib/plot_piston_data.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0 \
        || echo "  WARNING: plot_piston_data.py failed (skipping)"
elif [ "$FOLDER" = "shear_slab" ]; then
    # shear_slab writes per-strain (_g<strain>) files in its own schema; use the
    # dedicated sweep plotter instead of the compress/flow stress + piston plots.
    echo "Generating shear stress-strain sweep plots (per strain: $STRAINS)..."
    python "$SCRIPT_DIR/lib/plot_shear_strain_sweep.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" "$STRAINS"
else
    echo "Generating stress profiles..."
    python "$SCRIPT_DIR/lib/plot_stress_profiles.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0

    echo "Generating piston plots..."
    python "$SCRIPT_DIR/lib/plot_piston_data.py" "." "${DATANAME}_${INTERACTION}_${TOTSTEPS}" 0
fi


echo "======================================"
echo "Done! Results are in: $WORK_DIR"
echo "======================================"