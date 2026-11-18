# Site settings shared by run_lammps.sh, continue_sim.sh, postprocess.sh, batch_prologue.sh,
# submit.sh and the allocation benchmark. Every path, module and MPI setting lives here, one
# block per cluster; the scripts read these names and nothing else. Each value can be overridden
# from the environment.
#
# CLUSTER = expanse | pod. Detected from the Slurm cluster name and the host name; set it in the
# environment to force one (CLUSTER=pod bash scripts/submit.sh ...). submit.sh exports it into
# the job, so compute nodes never have to guess.

detect_cluster() {
    local name host
    name="${SLURM_CLUSTER_NAME:-}"
    [ -n "$name" ] || name=$(scontrol show config 2>/dev/null | awk '/^ClusterName/ {print $3}')
    host=$(hostname -f 2>/dev/null || hostname)
    case "$(printf '%s:%s' "$name" "$host" | tr '[:upper:]' '[:lower:]')" in
        expanse:*|*:*expanse*|*:exp-*) echo expanse ;;
        pod*:*|*:pod*|*:*cnsi*)        echo pod ;;
        *) if [ -d /expanse/lustre ]; then echo expanse; else echo unknown; fi ;;
    esac
}
CLUSTER="$(printf '%s' "${CLUSTER:-$(detect_cluster)}" | tr '[:upper:]' '[:lower:]')"
if [ "$CLUSTER" = unknown ]; then
    # a Mac (dry runs, tests) lands here; so would a cluster this file does not know
    command -v sbatch >/dev/null 2>&1 \
        && echo "WARNING: lib/site.sh does not recognise this cluster; using the Expanse settings. Set CLUSTER=expanse|pod." >&2
    CLUSTER=expanse
fi

LAMMPS_WORK_DIR="${LAMMPS_WORK_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
INPUT_DATA_DIR="${LAMMPS_INPUT_DATA_OVERRIDE:-$HOME/Documents/lammps_data/input_data}"
RUNS_DIR="$HOME/Documents/lammps_runs"          # run dirs: $RUNS_DIR/<folder>/<folder>_<dataname>_<interaction>_<stamp>

# MPI layout: MPI_MAP_BY=core fills a node with ranks before the next one, node deals them
# round-robin. Expanse allocation benchmark (2026-11-07): core beat node by 1.4-1.6x on several
# nodes, and 128 ranks x 1 thread beat every hybrid split. Pod: not measured yet with these
# settings (the archived Pod runner used node); simulations/allocation_benchmark measures it.
MPI_MAP_BY="${MPI_MAP_BY:-core}"               # core | node
MPI_BIND_TO="${MPI_BIND_TO:-core}"

# Wall-clock margin (s) before the SLURM limit at which a deck with -var wall_timeout stops and,
# with CHAIN=1, run_lammps.sh submits the continuation. Covers the final write_data, the MPI
# teardown grace, the sbatch and the plots.
WALL_MARGIN="${WALL_MARGIN:-2400}"

case "$CLUSTER" in
# ---- Expanse (SDSC): 128-core AMD EPYC 7742 nodes, HDR InfiniBand, charged in SUs ---------------
expanse)
    CORES_PER_NODE=128
    LMP_BIN="${LMP_BIN:-/home/dpollard/software/lammps/22Jul2025_update3/mpi-omp/gcc/10.2.0/openmpi/4.1.3/lammps-22Jul2025/build/lmp}"
    SCRATCH_DIR="${SCRATCH_DIR:-/expanse/lustre/scratch/$USER/temp_project}"   # trajectories
    LAMMPS_MODULES=(gcc/10.2.0 openmpi/4.1.3 python/3.8.12)
    PYTHON_MODULE=anaconda3/2021.05/q4munrg          # numpy + matplotlib for the plotters

    # allocation benchmark (simulations/allocation_benchmark/submit_benchmark.sh)
    BENCH_PARTITION=compute                                        # label in the run tags
    BENCH_SBATCH=(--partition=compute --account=csb197 --mem=0)    # on every benchmark sbatch
    BENCH_LAYOUTS="128x1 64x2 32x4 16x8"                           # ranks x threads per node
    BENCH_TIME_MULT=1                                              # scales the job time limits

    # Open MPI 4 / UCX settings for Expanse compute nodes.
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
    ;;

# ---- Pod (CNSI, UCSB): 40-core nodes (2 x 20-core Xeon Gold 6148), Omni-Path, no SU charge ------
# Partitions (sinfo, 2026-11-18): batch = 62 nodes, 58 d limit, 3 d default; short = 1 node, 2 h limit;
# largemem (1 TB) and gpu (V100) are not used here. /scratch is a 19 TB SSD array, /bigscratch
# 70 TB NVMe; neither is backed up.
# The binary is the Pod build of lammps_build/build_lammps_pod.sh (22Jul2025 update 3, the Expanse
# package list, gnu14 + openmpi5, built 2026-11-18). The fabric settings are from
# archive/pod_bridges/run_lammps_pod.sh (last used 2026-10). `bash scripts/cluster_check.sh`
# compares this block with the machine.
pod)
    CORES_PER_NODE=40
    LMP_BIN="${LMP_BIN:-$HOME/software/lammps/22Jul2025_update3/mpi-omp/gnu14/14.2.0/openmpi5/5.0.7/lammps-22Jul2025/build/lmp}"
    SCRATCH_DIR="${SCRATCH_DIR:-/scratch/$USER}"
    LAMMPS_MODULES=(gnu14 openmpi5)
    PYTHON_MODULE=miniconda                          # plus the conda env below
    PYTHON_CONDA_ENV="${PYTHON_CONDA_ENV:-lammps_analysis}"

    BENCH_PARTITION=batch                            # label in the run tags; batch is Pod's regular partition
    BENCH_SBATCH=()
    BENCH_LAYOUTS="40x1 20x2 10x4 8x5"
    BENCH_TIME_MULT=4                                # one Pod node has 40 cores against Expanse's 128

    # Open MPI 5 over Intel Omni-Path (PSM2).
    mpi_env() {
        export OMPI_MCA_pml='cm'
        export OMPI_MCA_mtl='psm2'
        export OMPI_MCA_btl='self,vader'
        export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
    }
    load_lammps_modules() {
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load "${LAMMPS_MODULES[@]}"
    }
    load_python_module() {
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load "$PYTHON_MODULE"
        # shellcheck disable=SC1091
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$PYTHON_CONDA_ENV"
    }
    ;;
*)  echo "ERROR: lib/site.sh has no settings for CLUSTER='$CLUSTER' (expanse | pod)." >&2
    return 1 2>/dev/null || exit 1 ;;
esac

# batch_header <file.batch>: the #SBATCH lines of a batch file for THIS cluster. A batch file
# carries its Expanse allocation as ordinary #SBATCH lines (so plain `sbatch` works there) and the
# allocation of every other cluster on one line tagged with the cluster's name:
#     #POD --nodes=4 --ntasks-per-node=40 --cpus-per-task=1 --time=96:00:00
# For Expanse this prints the #SBATCH lines as written. For another cluster it keeps
# --job-name / --output / --error and replaces the rest with the options of the tagged line.
# Returns 1 when the file has no line for this cluster.
batch_header() {
    local tag line opt
    if [ "$CLUSTER" = expanse ]; then grep '^#SBATCH' "$1"; return 0; fi
    tag="#$(printf '%s' "$CLUSTER" | tr '[:lower:]' '[:upper:]')"
    line=$(grep -m1 "^${tag}[[:space:]]" "$1") || return 1
    grep -E '^#SBATCH[[:space:]]+--(job-name|output|error)=' "$1"
    for opt in ${line#"$tag"}; do
        case "$opt" in '#'*) break ;; esac          # trailing comment on the line
        echo "#SBATCH $opt"
    done
}

# Older Slurm versions do not export SLURM_JOB_END_TIME, which the wall-clock budget
# (wall_timeout_seconds, launch_lammps in lib/run_common.sh) is computed from: ask squeue.
set_job_end_time() {
    local end
    [ -n "${SLURM_JOB_ID:-}" ] && [ -z "${SLURM_JOB_END_TIME:-}" ] || return 0
    end=$(squeue -h -j "$SLURM_JOB_ID" -o %e 2>/dev/null | head -1)
    end=$(date -d "$end" +%s 2>/dev/null) || return 0
    [ -n "$end" ] && export SLURM_JOB_END_TIME="$end"
    return 0
}
