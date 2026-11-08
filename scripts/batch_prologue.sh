# batch_prologue.sh -- the first command of every .batch, right after its #SBATCH header:
#     source "$HOME/Documents/lammps_work/scripts/batch_prologue.sh"
# It pulls lammps_work (serialised, fast-forward only: node_git_sync.sh) and, on the first pass,
# re-execs the freshly pulled copy of the batch, so the knobs below the source line come from the
# current file and not from the copy sbatch took at submission. Then it sets the MPI environment,
# loads the LAMMPS modules and exports LAMMPS_WORK_DIR. Outside SLURM it only sets the paths.

LAMMPS_WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export LAMMPS_WORK_DIR
source "$LAMMPS_WORK_DIR/scripts/lib/site.sh"
source "$LAMMPS_WORK_DIR/scripts/node_git_sync.sh"

if [ -n "${SLURM_JOB_ID:-}" ]; then
    sync_lammps_work || true
    if [ -z "${LAMMPS_BATCH_SYNCED:-}" ]; then
        export LAMMPS_BATCH_SYNCED=1
        # the submitted script's path, as sbatch recorded it
        _batch=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/.*Command=\([^ ]*\).*/\1/p' | head -1)
        case "$_batch" in /*) ;; ?*) _batch="${SLURM_SUBMIT_DIR:-$PWD}/$_batch" ;; esac
        case "$_batch" in
            "$LAMMPS_WORK_DIR"/*) [ -f "$_batch" ] && { echo ">>> re-running the synced copy of $_batch"; exec bash "$_batch"; } ;;
        esac
        unset _batch
    fi
    mpi_env
    load_lammps_modules
    echo ">>> job ${SLURM_JOB_ID} (${SLURM_JOB_NAME:-?}): ${SLURM_NNODES:-?} node(s) x ${SLURM_NTASKS_PER_NODE:-?} tasks x ${OMP_NUM_THREADS:-1} thread(s); lammps_work @ $(git -C "$LAMMPS_WORK_DIR" rev-parse --short HEAD 2>/dev/null || echo '?')"
fi
