#!/usr/bin/env bash
# =============================================================================
# volmix_sweep.sh
# =============================================================================
# Submits a chained SLURM pipeline for slab_with_support at 11 pressures
# P* = 1.0, 1.1, 1.2, ..., 2.0 (0.1 increments).
#
# Widom/chempot is suppressed for all sweep runs (SKIP_WIDOM=1).
#
# Pipeline per pressure:
#   [1] slab_with_support   (600k steps, NPT at P*)
#         → copies final_config to lammps_data/slab_with_support/
#   [2] split_gel.py        (splits polymer-only and solvent-only .data files)
#         → polymer_only → lammps_data/polymer_pure/
#         → solvent_only → lammps_data/solvent_pure/
#   [3a] solvent_pure       (100k steps, NPT at same P*)
#   [3b] polymer_pure       (100k steps, NPT at same P*)
#         [3a] and [3b] run concurrently after [2]
#
# Concurrency: at most 3 pressure pipelines run simultaneously.
# slab_p{n} is gated on afterok of sol_p{n-3} and pol_p{n-3}.
#
# Run from scripts/ on the Expanse login node:
#   bash volmix_sweep.sh
#
# Job logs  → ~/Documents/lammps_runs/volmix_sweep_logs/
# Manifests → ~/Documents/lammps_runs/sweep_manifest/
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$(cd "${SCRIPT_DIR}/../scripts" 2>/dev/null && pwd || echo "${SCRIPT_DIR}")"
# SCRIPTS_DIR resolves to lammps_work/scripts/ whether this file lives in
# simulations/ (one level up) or scripts/ itself (same dir, fallback).

BASE_DATANAME="slab_support_5beads_tall_rho04"
INTERACTION="1.0_1.0"
PURE_INTERACTION="1.0_0.0"
SLAB_STEPS=600000
PURE_STEPS=100000

LAMMPS_DATA="$HOME/Documents/lammps_data"
LAMMPS_RUNS="$HOME/Documents/lammps_runs"
SLAB_DATA_DIR="${LAMMPS_DATA}/slab_with_support"
SOL_DATA_DIR="${LAMMPS_DATA}/solvent_pure"
POL_DATA_DIR="${LAMMPS_DATA}/polymer_pure"
MANIFEST_DIR="${LAMMPS_RUNS}/sweep_manifest"
LOG_DIR="${LAMMPS_RUNS}/volmix_sweep_logs"

mkdir -p "$MANIFEST_DIR" "$LOG_DIR" "$SLAB_DATA_DIR" "$SOL_DATA_DIR" "$POL_DATA_DIR"

if ! command -v sbatch &>/dev/null; then
    echo "ERROR: sbatch not found. Are you on the Expanse login node?"
    exit 1
fi

PRESSURES=(1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)
WINDOW=3   # max concurrent pressure pipelines

# Track sol/pol JIDs in order so we can gate later pipelines on them
SOL_JIDS=()
POL_JIDS=()

echo "======================================"
echo "Volume of mixing sweep submission"
echo "Pressures: ${PRESSURES[*]}"
echo "Window: ${WINDOW} concurrent pipelines"
echo "======================================"

for i in "${!PRESSURES[@]}"; do
    P="${PRESSURES[$i]}"

    DATANAME="${BASE_DATANAME}_pstar${P}"
    TOTSTEPS=$SLAB_STEPS
    SOL_DATANAME="final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_solvent_only"
    POL_DATANAME="final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_polymer_only"

    # Symlink so run_lammps.sh finds the right input file by DATANAME
    SRC_DATA="${SLAB_DATA_DIR}/${BASE_DATANAME}.data"
    LNK_DATA="${SLAB_DATA_DIR}/${DATANAME}.data"
    if [ ! -e "$LNK_DATA" ] && [ ! -L "$LNK_DATA" ]; then
        ln -s "$SRC_DATA" "$LNK_DATA"
    fi

    # ------------------------------------------------------------------
    # Sliding-window gate: slab[i] waits for sol[i-WINDOW] and pol[i-WINDOW]
    # ------------------------------------------------------------------
    SLAB_DEP=""
    if [ "$i" -ge "$WINDOW" ]; then
        GATE_SOL="${SOL_JIDS[$((i - WINDOW))]}"
        GATE_POL="${POL_JIDS[$((i - WINDOW))]}"
        SLAB_DEP="--dependency=afterok:${GATE_SOL}:${GATE_POL}"
    fi

    # ------------------------------------------------------------------
    # Job 1: slab_with_support  (Widom suppressed via SKIP_WIDOM=1)
    # ------------------------------------------------------------------
    SLAB_BATCH=$(mktemp /tmp/slab_volmix_p${P}_XXXX.batch)
    cat > "$SLAB_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=slab_p${P}
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=5:00:00
#SBATCH --output=${LOG_DIR}/slab_p${P}.o%j.%N

declare -xr OMPI_UNIT='core'
declare -xr OMPI_MCA_btl='self,vader'
declare -xr UCX_LOG_LEVEL='ERROR'
declare -xr UCX_TLS='self,cma,shm,rc,ud,dc'
declare -xr UCX_NET_DEVICES='mlx5_2:1'
declare -xir UCX_MAX_RNDV_RAILS=1
declare -xir OMP_NUM_THREADS="\${SLURM_CPUS_PER_TASK}"

module reset
module load gcc/10.2.0
module load openmpi/4.1.3
module load python/3.8.12

export SKIP_WIDOM=1

cd "${SCRIPTS_DIR}" || { echo "ERROR: cd to SCRIPTS_DIR failed"; exit 1; }
./run_lammps.sh "slab_with_support" "${DATANAME}" "${INTERACTION}" \
    "${SLAB_STEPS}" "0" "" "${P}"

WORK_DIR=\$(ls -dt "\$HOME/Documents/lammps_runs/slab_with_support_${DATANAME}_${INTERACTION}_"* 2>/dev/null | head -1)
if [ -z "\$WORK_DIR" ]; then
    echo "ERROR: Could not find slab work directory for ${DATANAME}"
    exit 1
fi

FINAL_CONFIG="\${WORK_DIR}/final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}.data"
if [ ! -f "\$FINAL_CONFIG" ]; then
    echo "ERROR: Final config not found: \$FINAL_CONFIG"
    exit 1
fi

cp "\$FINAL_CONFIG" "${SLAB_DATA_DIR}/"
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_slab.workdir"
HEREDOC

    # shellcheck disable=SC2086
    SLAB_JID=$(sbatch --parsable $SLAB_DEP "$SLAB_BATCH")
    echo "P=${P}: submitted slab_with_support  JID=${SLAB_JID}${SLAB_DEP:+  (gated on ${SLAB_DEP#--dependency=afterok:})}"

    # ------------------------------------------------------------------
    # Job 2: split_gel.py  (depends on slab)
    # ------------------------------------------------------------------
    SPLIT_BATCH=$(mktemp /tmp/split_volmix_p${P}_XXXX.batch)
    cat > "$SPLIT_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=split_p${P}
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --output=${LOG_DIR}/split_p${P}.o%j.%N

module reset
module load gcc/10.2.0
module load python/3.8.12

INPUT_FILE="${SLAB_DATA_DIR}/final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}.data"
if [ ! -f "\$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: \$INPUT_FILE"
    exit 1
fi

python3 "${SCRIPTS_DIR}/split_gel.py" "\$INPUT_FILE" \
    --polymer-dir "${POL_DATA_DIR}" \
    --solvent-dir "${SOL_DATA_DIR}"

echo "Split complete:"
ls -lh "${POL_DATA_DIR}/final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_polymer_only.data"
ls -lh "${SOL_DATA_DIR}/final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_solvent_only.data"
HEREDOC

    SPLIT_JID=$(sbatch --parsable --dependency=afterok:${SLAB_JID} "$SPLIT_BATCH")
    echo "P=${P}: submitted split_gel           JID=${SPLIT_JID}  (after ${SLAB_JID})"

    # ------------------------------------------------------------------
    # Job 3a: solvent_pure  (depends on split)
    # ------------------------------------------------------------------
    SOL_BATCH=$(mktemp /tmp/solvent_volmix_p${P}_XXXX.batch)
    cat > "$SOL_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=sol_p${P}
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=1:00:00
#SBATCH --output=${LOG_DIR}/solvent_p${P}.o%j.%N

declare -xr OMPI_UNIT='core'
declare -xr OMPI_MCA_btl='self,vader'
declare -xr UCX_LOG_LEVEL='ERROR'
declare -xr UCX_TLS='self,cma,shm,rc,ud,dc'
declare -xr UCX_NET_DEVICES='mlx5_2:1'
declare -xir UCX_MAX_RNDV_RAILS=1
declare -xir OMP_NUM_THREADS="\${SLURM_CPUS_PER_TASK}"

module reset
module load gcc/10.2.0
module load openmpi/4.1.3
module load python/3.8.12

cd "${SCRIPTS_DIR}" || { echo "ERROR: cd to SCRIPTS_DIR failed"; exit 1; }
./run_lammps.sh "solvent_pure" "${SOL_DATANAME}" "${PURE_INTERACTION}" \
    "${PURE_STEPS}" "0" "" "${P}"

WORK_DIR=\$(ls -dt "\$HOME/Documents/lammps_runs/solvent_pure_${SOL_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_solvent.workdir"
HEREDOC

    SOL_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} "$SOL_BATCH")
    echo "P=${P}: submitted solvent_pure         JID=${SOL_JID}   (after ${SPLIT_JID})"

    # ------------------------------------------------------------------
    # Job 3b: polymer_pure  (depends on split, concurrent with solvent_pure)
    # ------------------------------------------------------------------
    POL_BATCH=$(mktemp /tmp/polymer_volmix_p${P}_XXXX.batch)
    cat > "$POL_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=pol_p${P}
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=5:00:00
#SBATCH --output=${LOG_DIR}/polymer_p${P}.o%j.%N

declare -xr OMPI_UNIT='core'
declare -xr OMPI_MCA_btl='self,vader'
declare -xr UCX_LOG_LEVEL='ERROR'
declare -xr UCX_TLS='self,cma,shm,rc,ud,dc'
declare -xr UCX_NET_DEVICES='mlx5_2:1'
declare -xir UCX_MAX_RNDV_RAILS=1
declare -xir OMP_NUM_THREADS="\${SLURM_CPUS_PER_TASK}"

module reset
module load gcc/10.2.0
module load openmpi/4.1.3
module load python/3.8.12

cd "${SCRIPTS_DIR}" || { echo "ERROR: cd to SCRIPTS_DIR failed"; exit 1; }
./run_lammps.sh "polymer_pure" "${POL_DATANAME}" "${PURE_INTERACTION}" \
    "${PURE_STEPS}" "0" "" "${P}"

WORK_DIR=\$(ls -dt "\$HOME/Documents/lammps_runs/polymer_pure_${POL_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_polymer.workdir"
HEREDOC

    POL_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} "$POL_BATCH")
    echo "P=${P}: submitted polymer_pure         JID=${POL_JID}   (after ${SPLIT_JID})"

    # Record JIDs for sliding-window gating of future pressures
    SOL_JIDS+=("$SOL_JID")
    POL_JIDS+=("$POL_JID")

    echo "--------------------------------------"
done

echo "======================================"
echo "All jobs submitted."
echo "Monitor with: squeue -u \$USER"
echo "Results manifest: ${MANIFEST_DIR}/"
echo "Volume of mixing analysis: volume_of_mixing.ipynb"
echo "======================================"
