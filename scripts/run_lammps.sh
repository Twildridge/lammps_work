#!/bin/bash
# ── Whole-file parse guard ────────────────────────────────────────────────────
# Bash reads a script INCREMENTALLY, re-reading from the file after every
# command it forks. This script lives in the git checkout on NFS home and runs
# for days on a compute node, while git pulls from OTHER NFS clients (lsync on
# the login node, the self-sync at the start of every other job) replace the
# file's inode underneath it. Cross-client NFS does not keep an unlinked open
# file alive, so the next read after the replacement fails (stale handle), bash
# treats it as EOF and exits SILENTLY with status 0. Job 54398786
# (triaxial_compression_two_pist, 43 h) lost its whole post-processing step this
# way: commit f9e0772 rewrote this file 7 h into the run, and the shell died the
# moment the mpirun polling loop ended (log ends at "Total wall time", no
# post-processing banner, no error, exit 0). Wrapping the body in ONE { ... }
# group makes bash parse the entire file before executing a single command,
# and the trailing `exit` stops it from ever reading past the group. Nothing
# else changes: a group runs in the current shell, so cd/export/exit behave
# exactly as before. Same guard in continue_sim.sh. (2026-09-26)
{
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
    echo "  env knobs: PISTON_TRANSPARENT BARO_MODE (slab_with_support); STRAINS / COMPRESSIONS (sweeps);"
    echo "             DP_PISTON PISTON_MASS C_PIST_FRAC NPT_PISTON_STEPS SETTLE_HALT (triaxial_*_two_pist)"
    echo "             DRIVE_SPLIT (triaxial_compression*: piston share of the gap closure, 0.5 = symmetric)"
    echo "             UNLOAD / UNLOAD_STEPS / UNLOAD_MARGIN (triaxial_compression_two_pist: re-swelling hold after each level)"
    echo "             compress_slab reads the two-piston knobs too, plus SEAT_MODE (rg|bb) and HOLD_AUTO (1|0)"
    exit 1
fi

FOLDER=$1
DATANAME=$2
INTERACTION=$3
NSTEPS=$4
TOTSTEPS=$NSTEPS
PRESS_TARGET=${6:-1.5}  # Default pressure; overrides press_target in .lmp file
VEL_SEED=${7:-12345}    # RNG seed for create_velocity and fix langevin; vary per replica
SKIP_WIDOM=${SKIP_WIDOM:-1}  # Cavity-Widom output is OFF by default (feature archived 2026-08;
                             # see lammps_work/archive/). Set to 0 (via env) to re-enable the
                             # widom_traj dumps and cavity_widom.py post-processing.
# PISTON_TRANSPARENT=1 (env) zeroes the solvent-piston pair (slab_with_support only):
# trial to test whether the piston sheet's wall tension is what makes the polymer
# partial stress anisotropic under the aniso barostat (2026-09-07). Default 0.
PISTON_TRANSPARENT=${PISTON_TRANSPARENT:-0}
# Two-piston (feed/permeate, NPT-piston) knobs (2026-09-16), read only by
# triaxial_permeation_two_pist / triaxial_compression_two_pist (other engines
# ignore them).  Set in the *_two_pist.batch files; same env-var pattern as
# PISTON_TRANSPARENT.  DP_PISTON is the ONE permeation pressure step (no sweep).
DP_PISTON=${DP_PISTON:-0.1}                 # P_feed = P_target + DP_PISTON (permeation only)
PISTON_MASS=${PISTON_MASS:-1000}            # mass of every piston bead (types 5/6/7)
C_PIST_FRAC=${C_PIST_FRAC:-1.0}             # piston damping as a fraction of the critical value
NPT_PISTON_STEPS=${NPT_PISTON_STEPS:-1000000}  # Phase-1 NPT-piston settle length
SETTLE_HALT=${SETTLE_HALT:-0}               # 1 = halt the settle early once both pistons are at rest
DRIVE_SPLIT=${DRIVE_SPLIT:-0.5}             # compression decks: piston share of each level's gap closure
                                            # (0.5 = symmetric drive, piston down + support up; 1.0 = old top-only)
# triaxial_permeation_two_pist (2026-10-06): dP ramp (N_RAMP stages x RAMP_STEPS steps), profile snapshots over
# the production run (PERM_STRESS_CURVES) and the displacement-profile block length (DISP_NFREQ steps)
N_RAMP=${N_RAMP:-100}
RAMP_STEPS=${RAMP_STEPS:-10000}
PERM_STRESS_CURVES=${PERM_STRESS_CURVES:-10}
DISP_NFREQ=${DISP_NFREQ:-50000}
DISP_FINE_NFREQ=${DISP_FINE_NFREQ:-0}       # triaxial_compression_two_pist: block length (steps) of the fine-cadence
DISP_FINE_BIN=${DISP_FINE_BIN:-1.0}         # hold displacement profile disp_z_polymer_fine and its bin width; 0 = off
UNLOAD=${UNLOAD:-0}                         # triaxial_compression_two_pist (2026-10-09): 1 = after every level's hold drive the
UNLOAD_STEPS=${UNLOAD_STEPS:-0}             # plates back to the seated gap and record a free re-swelling hold of UNLOAD_STEPS
UNLOAD_MARGIN=${UNLOAD_MARGIN:-0.0}         # (0 = the level's hold length), each plate UNLOAD_MARGIN sigma beyond its seat
CALIB_FRAMES=${CALIB_FRAMES:-5}          # calibration-dump frames near run end (polymer_pure /
CALIB_DUMP_EVERY=${CALIB_DUMP_EVERY:-2000}  # solvent_pure only; other engines ignore these vars)
# PRERELAXED=1 tells polymer_pure to skip its Stage 0 harmonic pre-relaxation:
# the input data file already holds an equilibrated FENE configuration (e.g. the
# isolated_* files from isolate_gel.py) rather than a fresh diamond lattice.
# Velocities are still seeded from VEL_SEED and Stage 1 still runs, so replicas
# stay independent. Not the same as continue_sim.sh's cont=1, which also assumes
# the file carries its own Velocities section. Engines that don't define the
# variable simply ignore it.
PRERELAXED=${PRERELAXED:-0}
if [ -z "${STRAINS:-}" ]; then
    echo ">>> WARNING: STRAINS is unset — falling back to single strain 0.1."
    echo ">>>          For shear_slab this means NO sweep. If you intended a sweep,"
    echo ">>>          your shear_slab.batch on the cluster was likely stale at sbatch"
    echo ">>>          time. Pull lammps_work and resubmit (the batch now self-syncs)."
fi
STRAINS=${STRAINS:-0.1}      # Space-separated shear-strain list (shear_slab only); passed as a
                             # LAMMPS index variable. Default 0.1 = single operating point.
# COMPRESSIONS: space-separated cumulative compression-strain list (triaxial_compression_two_pist,
# and the archived one-piston triaxial_compression), passed as a LAMMPS index variable
# exactly like STRAINS. Set in the compression .batch via STRAIN_TARGETS=(...) ->
# COMPRESSIONS. Default 0.1 = single operating point (reproduces the original single-run behaviour). Same stale-batch caveat as
# STRAINS applies: the batch must self-sync BEFORE exporting COMPRESSIONS.
if [ -z "${COMPRESSIONS:-}" ]; then
    echo ">>> NOTE: COMPRESSIONS unset — the compression decks fall back to single strain 0.1."
fi
COMPRESSIONS=${COMPRESSIONS:-0.1}

# compress_slab knobs (2026-09-29; the deck is the bulk-modulus counterpart of
# triaxial_compression_two_pist and reads the SAME variables: COMPRESSIONS is its list of
# cumulative VOLUMETRIC strain targets, and PRESS_TARGET / PISTON_MASS / C_PIST_FRAC /
# NPT_PISTON_STEPS / SETTLE_HALT / DRIVE_SPLIT mean what they mean there).  Two more:
#   SEAT_MODE  rg (default) | bb : the gel face the six plates seat on (compress_slab.lmp)
#   HOLD_AUTO  1 (default) = size every hold from tau_1 = L^2/(pi^2 D_c); 0 = NSTEPS flat
# BARO_MODE (slab_with_support only): aniso (default) | iso | z -- barostat coupling (see
#   slab_with_support.lmp).  PIN_GEL=1 (slab_with_support, isolated gel): zero the polymer's
#   linear + angular momentum so the block neither drifts nor turns away from the box axes.
SEAT_MODE=${SEAT_MODE:-rg}
HOLD_AUTO=${HOLD_AUTO:-1}
BARO_MODE=${BARO_MODE:-aniso}
PIN_GEL=${PIN_GEL:-0}                # slab_with_support: 1 = zero the gel's linear + angular momentum (isolated gel)

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

# Snapshot the input script into the run dir and run THAT copy, not the file in
# the git working tree. LAMMPS reads its input lazily -- it keeps the file open
# and only reads the next command when the current `run` returns, so a multi-day
# job holds the repo file open for its whole life. Editing that file mid-run
# shifts the byte offsets under LAMMPS: job 54220814 (slab_with_support, 14M
# steps) had `compute mobile_temp` inserted 16 lines above the read point while
# its last `run` block was in flight, and when that block finished LAMMPS resumed
# mid-line and spliced two unrelated lines into a garbage command --
# "ERROR: Unknown command: jumpniso_${dataname}..." -- after 12 h of good MD.
# Copying first makes mid-run edits harmless, and archives the exact script that
# produced each run's results alongside them. (2026-09-08)
cp "$LAMMPS_FILE" "$WORK_DIR/$(basename "$LAMMPS_FILE")"
LAMMPS_FILE="$WORK_DIR/$(basename "$LAMMPS_FILE")"

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

# ── Pull latest scripts from GitHub before running ────────────────────────────
# Serialised, fast-forward-only, best-effort. See scripts/node_git_sync.sh.
source "$SCRIPT_DIR/node_git_sync.sh"
echo ">>> Syncing lammps_work from GitHub..."
sync_lammps_work || true

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

# ── Wall-clock budget for the deck (2026-09-30) ──────────────────────────────
# A deck that honours `-var wall_timeout` (slab_with_support) ends its production
# loop on its own this many seconds after LAMMPS starts, writes the snapshot it
# was in and reports `PRODUCTION_STEPS_DONE ... STATUS INCOMPLETE` in its log.
# Budget = time left in the SLURM job minus WALL_MARGIN: the final write_data of
# a ~1 M-atom system, the MPI teardown grace above, the continuation submission
# and post-processing all fit comfortably in 40 min.  0 = unlimited (no SLURM).
# Decks without the variable ignore it.  With CHAIN=1 (env, set by the .batch)
# the INCOMPLETE line makes this script submit the continuation (see below).
WALL_MARGIN=${WALL_MARGIN:-2400}
WALL_TIMEOUT=0
if [ -n "${SLURM_JOB_END_TIME:-}" ]; then
    WALL_TIMEOUT=$(( SLURM_JOB_END_TIME - $(date +%s) - WALL_MARGIN ))
    [ "$WALL_TIMEOUT" -gt 0 ] || WALL_TIMEOUT=0
    echo "LAMMPS wall-clock budget: ${WALL_TIMEOUT}s (-var wall_timeout = SLURM end - ${WALL_MARGIN}s; CHAIN=${CHAIN:-0})"
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
    -var prerelaxed $PRERELAXED \
    -var skip_widom $SKIP_WIDOM \
    -var calib_frames $CALIB_FRAMES \
    -var calib_dump_every $CALIB_DUMP_EVERY \
    -var strains $STRAINS \
    -var strains_list "$STRAINS" \
    -var compressions $COMPRESSIONS \
    -var compressions_list "$COMPRESSIONS" \
    -var seat_mode "$SEAT_MODE" \
    -var hold_auto "$HOLD_AUTO" \
    -var baro_mode "$BARO_MODE" \
    -var pin_gel "$PIN_GEL" \
    -var piston_transparent "${PISTON_TRANSPARENT:-0}" \
    -var dp_piston "$DP_PISTON" \
    -var piston_mass "$PISTON_MASS" \
    -var c_pist_frac "$C_PIST_FRAC" \
    -var npt_piston_steps "$NPT_PISTON_STEPS" \
    -var settle_halt "$SETTLE_HALT" \
    -var drive_split "$DRIVE_SPLIT" \
    -var disp_fine_nfreq "$DISP_FINE_NFREQ" \
    -var disp_fine_bin "$DISP_FINE_BIN" \
    -var unload "$UNLOAD" \
    -var unload_steps "$UNLOAD_STEPS" \
    -var unload_margin "$UNLOAD_MARGIN" \
    -var n_ramp "$N_RAMP" \
    -var ramp_steps "$RAMP_STEPS" \
    -var perm_stress_curves "$PERM_STRESS_CURVES" \
    -var disp_nfreq "$DISP_NFREQ" \
    -var wall_timeout "$WALL_TIMEOUT" \
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
    # Never propagate a zero here: ${LAMMPS_RC:-1} substitutes only when RC is
    # unset or empty, so an RC of 0 on an incomplete log would exit 0 and report
    # a crashed run as a success. The log check is the authority, so a failure
    # it detects must always exit nonzero.
    if [ -z "${LAMMPS_RC}" ] || [ "${LAMMPS_RC}" -eq 0 ]; then
        exit 1
    fi
    exit "${LAMMPS_RC}"
fi

# ── Chain a continuation when the deck ran out of wall-clock (2026-09-30) ─────
# CHAIN=1 in the environment (exported by the .batch) plus the deck's
# `PRODUCTION_STEPS_DONE ... STATUS INCOMPLETE` line: continue_sim.sh --chain reads
# that line, resumes from the newest snapshot and submits the job for the
# remaining steps NOW, before post-processing, so the submission never competes
# with the wall limit.  The new job chains again if it, too, runs out of time.
# Run from the submit directory: --chain keys on this job's SLURM output file.
if [ "${CHAIN:-0}" = "1" ]; then
    if grep -qE '^>>> PRODUCTION_STEPS_DONE [0-9]+ .*STATUS INCOMPLETE' "$LAMMPS_LOG"; then
        echo ">>> Production incomplete (wall-clock budget reached) — submitting the continuation job..."
        ( cd "${SLURM_SUBMIT_DIR:-$SIM_DIR}" && \
          PRESS_TARGET="$PRESS_TARGET" bash "$SCRIPT_DIR/continue_sim.sh" --chain "$SLURM_JOB_ID" ) \
        || echo ">>> WARNING: continuation submission FAILED — submit it by hand from $SIM_DIR: continue_sim.sh --chain $SLURM_JOB_ID"
    else
        echo ">>> CHAIN=1: production complete (or the deck has no wall-clock budget) — nothing to chain."
    fi
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

# Never fall through the closing brace: bash would try to read the file again
# (see the parse guard at the top). Propagate postprocess.sh's status.
exit $?
}
