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
# Batch submission: submits WINDOW=1 pipelines (12 jobs) at a time.
# The last polymer job in each batch re-invokes this script with --from N
# to submit the next batch automatically.
#
# Usage:
#   bash volmix_sweep.sh              # start from P=1.0
#   bash volmix_sweep.sh --from 6    # resume from index 6 (P=1.6)
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

# Parse --from N argument (defaults to 0)
FROM=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

PRESSURES=(1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)
WINDOW=1
TOTAL=${#PRESSURES[@]}

# Clamp END to array bounds
END=$(( FROM + WINDOW ))
if [ "$END" -gt "$TOTAL" ]; then END=$TOTAL; fi

if [ "$FROM" -ge "$TOTAL" ]; then
    echo "All pressures already submitted."
    exit 0
fi

echo "======================================"
echo "Volume of mixing sweep — batch submission"
echo "Submitting indices ${FROM}–$((END-1)): ${PRESSURES[*]:$FROM:$((END-FROM))}"
echo "Remaining after this batch: $((TOTAL - END)) pressure(s)"
echo "======================================"

for (( i=FROM; i<END; i++ )); do
    P="${PRESSURES[$i]}"
    IS_LAST_IN_BATCH=$([ "$i" -eq "$((END-1))" ] && echo "yes" || echo "no")

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
    # Job 1: slab_with_support  (Widom suppressed via SKIP_WIDOM=1)
    # ------------------------------------------------------------------
    SLAB_BATCH=$(mktemp /tmp/slab_volmix_p${P}_XXXX.batch)
    cat > "$SLAB_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=slab_p${P}
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=5
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

WORK_DIR=\$(ls -dt "\$HOME/Documents/lammps_runs/slab_with_support/slab_with_support_${DATANAME}_${INTERACTION}_"* 2>/dev/null | head -1)
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

    SLAB_JID=$(sbatch --parsable "$SLAB_BATCH")
    echo "P=${P}: submitted slab_with_support  JID=${SLAB_JID}"

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

WORK_DIR=\$(ls -dt "\$HOME/Documents/lammps_runs/solvent_pure/solvent_pure_${SOL_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_solvent.workdir"
HEREDOC

    SOL_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} "$SOL_BATCH")
    echo "P=${P}: submitted solvent_pure         JID=${SOL_JID}   (after ${SPLIT_JID})"

    # ------------------------------------------------------------------
    # Job 3b: polymer_pure  (depends on split, concurrent with solvent_pure)
    # If this is the last job in the batch, trigger next batch on completion.
    # ------------------------------------------------------------------
    POL_BATCH=$(mktemp /tmp/polymer_volmix_p${P}_XXXX.batch)
    NEXT_FROM=$END   # value at script-write time

    # Build the optional next-batch trigger lines (expanded now, embedded literally)
    if [ "$IS_LAST_IN_BATCH" = "yes" ] && [ "$NEXT_FROM" -lt "$TOTAL" ]; then
        NEXT_BATCH_TRIGGER="echo \"Submitting next batch (--from ${NEXT_FROM})...\"
bash \"${SCRIPT_DIR}/volmix_sweep.sh\" --from ${NEXT_FROM}"
    else
        NEXT_BATCH_TRIGGER=""
    fi

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

WORK_DIR=\$(ls -dt "\$HOME/Documents/lammps_runs/polymer_pure/polymer_pure_${POL_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_polymer.workdir"

${NEXT_BATCH_TRIGGER}
HEREDOC

    POL_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} "$POL_BATCH")
    echo "P=${P}: submitted polymer_pure         JID=${POL_JID}   (after ${SPLIT_JID})${IS_LAST_IN_BATCH:+ [triggers next batch]}"

    echo "--------------------------------------"
done

echo "======================================"
echo "Batch submitted: ${PRESSURES[*]:$FROM:$((END-FROM))}"
if [ "$END" -lt "$TOTAL" ]; then
    echo "Next batch (indices ${END}+) will auto-submit when this batch completes."
fi
echo "Monitor with: squeue -u \$USER"
echo "Results manifest: ${MANIFEST_DIR}/"
echo "======================================"
