# Site settings shared by run_lammps.sh, continue_sim.sh, postprocess.sh and batch_prologue.sh.
# Every path and MPI setting lives here; the scripts read these names and nothing else.
# Each value can be overridden from the environment.

LAMMPS_WORK_DIR="${LAMMPS_WORK_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
LMP_BIN="${LMP_BIN:-/home/dpollard/software/lammps/22Jul2025_update3/mpi-omp/gcc/10.2.0/openmpi/4.1.3/lammps-22Jul2025/build/lmp}"
INPUT_DATA_DIR="${LAMMPS_INPUT_DATA_OVERRIDE:-$HOME/Documents/lammps_data/input_data}"
RUNS_DIR="$HOME/Documents/lammps_runs"          # run dirs: $RUNS_DIR/<folder>/<folder>_<dataname>_<interaction>_<stamp>
SCRATCH_DIR="${SCRATCH_DIR:-/expanse/lustre/scratch/$USER/temp_project}"   # trajectories

# MPI layout (allocation benchmark, 2026-11-07): ranks fill a node before the next one
# (MPI_MAP_BY=core) beat round-robin placement (node) by 1.4-1.6x; 128 x 1 beat every hybrid split.
MPI_MAP_BY="${MPI_MAP_BY:-core}"               # core | node
MPI_BIND_TO="${MPI_BIND_TO:-core}"

# Wall-clock margin (s) before the SLURM limit at which a deck with -var wall_timeout stops and,
# with CHAIN=1, run_lammps.sh submits the continuation. Covers the final write_data, the MPI
# teardown grace, the sbatch and the plots.
WALL_MARGIN="${WALL_MARGIN:-2400}"

LAMMPS_MODULES=(gcc/10.2.0 openmpi/4.1.3 python/3.8.12)
PYTHON_MODULE=anaconda3/2021.05/q4munrg          # numpy + matplotlib for the plotters

# Open MPI / UCX settings for Expanse compute nodes.
mpi_env() {
    export OMPI_MCA_btl='self,vader'
    export UCX_LOG_LEVEL='ERROR'
    export UCX_TLS='self,cma,shm,rc,ud,dc'
    export UCX_NET_DEVICES='mlx5_2:1'
    export UCX_MAX_RNDV_RAILS=1
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
}

load_lammps_modules() {
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module reset 2>/dev/null || true
    module load "${LAMMPS_MODULES[@]}"
}

load_python_module() {
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module unload python/3.8.12 2>/dev/null || true
    module load "$PYTHON_MODULE"
}
