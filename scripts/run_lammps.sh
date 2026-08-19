#!/bin/bash
if [ $# -lt 4 ]; then
    echo "Usage: ./run_lammps.sh <folder_name> <dataname> <interaction> <nsteps> [type] [press_target]"
    echo "Example (fresh run):  ./run_lammps.sh slab_with_support slab_support_5beads_... 1.5_1.4 20000"
    echo "Example (pure solvent P-sweep): ./run_lammps.sh pure_solvent pure_solvent_1000 1p0 0"
    echo "Example (pressure sweep): ./run_lammps.sh slab_with_support slab_support_pstar0.8 1.0_1.0 600000 \"\" 0.8"
    echo "  interaction format: epsSS_epsSP (e.g., 1.5_0.4), or epsSS only for pure_solvent (e.g., 1p0)"
    echo "  type: optional, 'stress' (adds 1), 'volume' (adds 2), or 'stressvol' (adds 3) to dataname"
    echo "  to continue a finished run, use continue_sim.sh instead of resubmitting this script"
    echo "  press_target: optional, overrides press_target in .lmp file (default: 1.5)"
    echo "  vel_seed: optional, RNG seed for create_velocity and fix langevin (default: 12345)"
    exit 1
fi

FOLDER=$1
DATANAME=$2
INTERACTION=$3
NSTEPS=$4
TOTSTEPS=$NSTEPS
PRESS_TARGET=${6:-1.5}  # Default pressure; overrides press_target in .lmp file
VEL_SEED=${7:-12345}    # RNG seed for create_velocity and fix langevin; vary per replica
SKIP_WIDOM=${SKIP_WIDOM:-0}  # Set to 1 (via env) to minimize Widom output and skip cavity_widom.py
if [ -z "${STRAINS:-}" ]; then
    echo ">>> WARNING: STRAINS is unset — falling back to single strain 0.1."
    echo ">>>          For shear_slab this means NO sweep. If you intended a sweep,"
    echo ">>>          your shear_slab.batch on the cluster was likely stale at sbatch"
    echo ">>>          time. Pull lammps_work and resubmit (the batch now self-syncs)."
fi
STRAINS=${STRAINS:-0.1}      # Space-separated shear-strain list (shear_slab only); passed as a
                             # LAMMPS index variable. Default 0.1 = single operating point.
# COMPRESSIONS: space-separated cumulative compression-strain list (triaxial_compression
# only), passed as a LAMMPS index variable exactly like STRAINS. Set in
# triaxial_compression.batch via COMPRESSIONS=(...). Default 0.1 = single operating
# point (reproduces the original single-run behaviour). Same stale-batch caveat as
# STRAINS applies: the batch must self-sync BEFORE exporting COMPRESSIONS.
if [ -z "${COMPRESSIONS:-}" ]; then
    echo ">>> NOTE: COMPRESSIONS unset — triaxial_compression falls back to single strain 0.1."
fi
COMPRESSIONS=${COMPRESSIONS:-0.1}

# COMPRESS_STAGES: space-separated CUMULATIVE volumetric-strain ladder
# (compress_slab only) — one value per compression stage, passed as a LAMMPS
# index variable exactly like STRAINS/COMPRESSIONS. Set in compress_slab*.batch
# via COMPRESSION_1, COMPRESSION_2, ... -> STAGE_TARGETS -> COMPRESS_STAGES.
# The number of stages is just the length of this list. Default reproduces
# compress_slab.lmp's original fixed 3-stage ladder. Same stale-batch caveat
# as STRAINS/COMPRESSIONS: the batch must self-sync BEFORE exporting this.
if [ -z "${COMPRESS_STAGES:-}" ]; then
    echo ">>> NOTE: COMPRESS_STAGES unset — compress_slab falls back to the default 0.015 0.030 0.045 ladder."
fi
COMPRESS_STAGES=${COMPRESS_STAGES:-"0.015 0.030 0.045"}

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
# LAMMPS_RUNS_OVERRIDE: when set, all run dirs go flat into that folder (no /${FOLDER}/ subdir).
# Used by volmix_sweep.sh to consolidate everything into lammps_runs/volmix_sweep/.
if [ -n "${LAMMPS_RUNS_OVERRIDE:-}" ]; then
    WORK_DIR="${LAMMPS_RUNS_OVERRIDE}/${FOLDER}_${DATANAME}_${INTERACTION}_${RUN_TIMESTAMP}"
else
    WORK_DIR="$HOME/Documents/lammps_runs/${FOLDER}/${FOLDER}_${DATANAME}_${INTERACTION}_${RUN_TIMESTAMP}"
fi
mkdir -p "$WORK_DIR"/{data_files,output_files/{stress_data,volume_data,piston_data,permeation_data,displacement_data,pair_data,chemical_potential},output_plots}

# Create trajectory directory in scratch and symlink to it
TRAJ_DIR="$SCRATCH_DIR/lammps_trajectories/${FOLDER}_${DATANAME}_${INTERACTION}_${RUN_TIMESTAMP}"
mkdir -p "$TRAJ_DIR"
ln -s "$TRAJ_DIR" "$WORK_DIR/traj_files"

echo "======================================"
echo "Working directory: $WORK_DIR"
echo "Trajectory directory (scratch): $TRAJ_DIR"
echo "======================================"

# Copy or link the data file — all input data lives in input_data/
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
echo "  press_target=$PRESS_TARGET"
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
if command -v git &>/dev/null; then
    git -C "$LAMMPS_WORK_DIR" pull --rebase --autostash || true
else
    echo ">>> git not available on this node — skipping sync"
fi

# CPU-only mode (NOT USING OMP FOR NOW (not on Expanse 2021 version)
echo "Running CPU-only with $SLURM_NTASKS tasks"

# ── Teardown-hang guard ───────────────────────────────────────────────────────
# On Expanse, ranks occasionally deadlock in MPI_Finalize / UCX cleanup AFTER
# LAMMPS has already printed "Total wall time" and written every output file.
# When that happens mpirun never returns, SLURM eventually kills the whole step,
# and the post-processing below never runs even though the science is complete.
# (Symptom: log ends at "Total wall time" with neither a post-processing banner
# nor a "LAMMPS failed" line — the script never regained control.)
#
# Defense: run mpirun under `timeout`, sized to fire a few minutes before the
# SLURM wall limit, so a hung teardown is force-killed and control returns here.
# Completion is then judged from the LAMMPS log (below), not the mpirun RC.
MPIRUN_TIMEOUT=""
if [ -n "${SLURM_JOB_END_TIME:-}" ]; then
    REMAIN=$(( SLURM_JOB_END_TIME - $(date +%s) - 300 ))   # 5-min buffer before wall limit
    if [ "$REMAIN" -gt 60 ]; then
        MPIRUN_TIMEOUT="timeout -k 60 ${REMAIN}s"           # SIGTERM, then SIGKILL after 60s
        echo "mpirun guarded by: ${MPIRUN_TIMEOUT} (force-returns ~5 min before wall limit)"
    fi
fi

# ── Launch LAMMPS in the background and watch the log ─────────────────────────
# See the teardown-hang note above. Rather than blocking on mpirun's return (which
# can hang for hours in MPI_Finalize/UCX cleanup AFTER the science is complete),
# we background it and poll log.lammps. LAMMPS prints "Total wall time:" as the
# very last line, AFTER every write_data/write_restart, so once it appears all
# output is on disk. We then allow a short grace period for a clean teardown and,
# if mpirun is still stuck, kill it and proceed to post-processing anyway.
LAMMPS_LOG="${WORK_DIR}/log.lammps"

$MPIRUN_TIMEOUT mpirun -n "${SLURM_NTASKS}" --bind-to "${OMPI_UNIT}" --map-by "node:pe=${OMP_NUM_THREADS}" \
    /home/dpollard/software/lammps/22Jul2025_update3/mpi-omp/gcc/10.2.0/openmpi/4.1.3/lammps-22Jul2025/build/lmp \
    -sf omp -pk omp $SLURM_CPUS_PER_TASK \
    -var dataname $DATANAME \
    -var interaction $INTERACTION \
    -var epsSS $EPSSS \
    -var epsSP $EPSSP \
    -var nsteps $NSTEPS \
    -var oldsteps 0 \
    -var totsteps $TOTSTEPS \
    -var nsteps_eq $NSTEPS_EQ \
    -var nsteps_prod $NSTEPS_PROD \
    -var press_target $PRESS_TARGET \
    -var vel_seed $VEL_SEED \
    -var skip_widom $SKIP_WIDOM \
    -var strains $STRAINS \
    -var strains_list "$STRAINS" \
    -var compressions $COMPRESSIONS \
    -var compressions_list "$COMPRESSIONS" \
    -var stage_targets $COMPRESS_STAGES \
    -var stage_targets_list "$COMPRESS_STAGES" \
    \
    -in $LAMMPS_FILE &
MPIRUN_PID=$!

TEARDOWN_GRACE=180   # seconds to allow a clean MPI_Finalize after completion
POLL_INTERVAL=15     # seconds between log checks
while kill -0 "$MPIRUN_PID" 2>/dev/null; do
    if grep -q "Total wall time" "$LAMMPS_LOG" 2>/dev/null; then
        echo ">>> Detected 'Total wall time' in log — LAMMPS complete; all output written."
        echo ">>> Allowing ${TEARDOWN_GRACE}s for a clean MPI teardown before proceeding..."
        WAITED=0
        while kill -0 "$MPIRUN_PID" 2>/dev/null && [ "$WAITED" -lt "$TEARDOWN_GRACE" ]; do
            sleep 10; WAITED=$((WAITED + 10))
        done
        if kill -0 "$MPIRUN_PID" 2>/dev/null; then
            echo ">>> mpirun still running ${TEARDOWN_GRACE}s after completion — assuming MPI/UCX teardown hang."
            echo ">>> Killing mpirun (PID ${MPIRUN_PID}) and proceeding to post-processing."
            kill -TERM "$MPIRUN_PID" 2>/dev/null; sleep 10
            kill -KILL "$MPIRUN_PID" 2>/dev/null; sleep 2
        fi
        break
    fi
    sleep "$POLL_INTERVAL"
done
wait "$MPIRUN_PID" 2>/dev/null
LAMMPS_RC=$?

# Judge completion from the LAMMPS log, not the mpirun exit code. A nonzero RC (or
# 128+signal from our own kill, or 124 from the timeout guard) can occur during a
# teardown that happens AFTER all output is on disk; those must NOT discard a good
# run. Only the absence of "Total wall time" means a genuine failure.
if [ "$LAMMPS_RC" -eq 124 ]; then
    echo "WARNING: mpirun hit the timeout guard (RC 124) — probable MPI/UCX teardown hang after completion."
fi
if grep -q "Total wall time" "$LAMMPS_LOG" 2>/dev/null; then
    if [ "$LAMMPS_RC" -ne 0 ]; then
        echo "NOTE: mpirun returned ${LAMMPS_RC}, but '${LAMMPS_LOG}' reached 'Total wall time' —"
        echo "      LAMMPS ran to completion; proceeding with post-processing."
    fi
else
    echo "LAMMPS did not reach 'Total wall time' in ${LAMMPS_LOG} (mpirun RC ${LAMMPS_RC})."
    echo "Treating as a genuine failure — skipping post-processing."
    exit "${LAMMPS_RC:-1}"
fi

# ── Post-processing ───────────────────────────────────────────────────────────
# All plot generation lives in postprocess.sh so the automatic pipeline and a
# manual re-run use the EXACT same module + script order. If a teardown hang or
# job kill ever skips this step, regenerate every plot by hand (no MD, no MPI):
#
#   bash scripts/postprocess.sh <run_dir> <folder> <dataname> <interaction> <totsteps> [oldsteps] [press_target]
#
# SKIP_WIDOM, STRAINS and COMPRESSIONS are read from the environment by postprocess.sh.
export SKIP_WIDOM STRAINS COMPRESSIONS
bash "$SCRIPT_DIR/postprocess.sh" \
    "$WORK_DIR" "$FOLDER" "$DATANAME" "$INTERACTION" "$TOTSTEPS" 0 "$PRESS_TARGET"