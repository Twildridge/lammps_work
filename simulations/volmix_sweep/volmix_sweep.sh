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
#         → copies final_config to lammps_data/input_data/
#   [2] split_gel.py        (writes polymer_only and solvent_only to lammps_data/input_data/)
#   [3a] solvent_pure       (100k steps, NPT at same P*)
#   [3b] polymer_pure       (100k steps, NPT at same P*)
#         [3a] and [3b] run concurrently after [2]
#
# All run directories land in ~/Documents/lammps_runs/volmix_sweep/
# (via LAMMPS_RUNS_OVERRIDE exported in each job script).
#
# Batch submission: WINDOW=1 pipeline at a time.  The last polymer job triggers
# the next batch via a lightweight shared-partition launcher job.
#
# Usage:
#   bash volmix_sweep.sh              # start from P=1.0
#   bash volmix_sweep.sh --from 6    # resume from index 6 (P=1.6)
#
# SLURM log  → simulations/volmix_sweep/volmix_sweep_YYYYMMDD_HHMMSS.log  (new file per run)
# Run data   → ~/Documents/lammps_runs/volmix_sweep/
# Manifests  → ~/Documents/lammps_runs/volmix_sweep/sweep_manifest/
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$(cd "${SCRIPT_DIR}/../../scripts" 2>/dev/null && pwd || echo "${SCRIPT_DIR}")"
# SCRIPTS_DIR resolves to lammps_work/scripts/ (two levels up from simulations/volmix_sweep/)

BASE_DATANAME="slab_support_5beads_tall_rho04"
INTERACTION="1.0_1.0"
PURE_INTERACTION="1.0_0.0"
SLAB_STEPS=600000
PURE_STEPS=100000

LAMMPS_DATA="$HOME/Documents/lammps_data"
LAMMPS_RUNS="$HOME/Documents/lammps_runs"
INPUT_DATA_DIR="${LAMMPS_DATA}/input_data"
SLAB_DATA_DIR="${LAMMPS_DATA}/slab_with_support"

# All volmix run dirs, manifests, and the single SLURM log go here
VOLMIX_RUNS="${LAMMPS_RUNS}/volmix_sweep"
MANIFEST_DIR="${VOLMIX_RUNS}/sweep_manifest"
LOG_DIR="${SCRIPT_DIR}"   # SLURM log lives alongside this script in simulations/volmix_sweep/

mkdir -p "$VOLMIX_RUNS" "$MANIFEST_DIR" "$SLAB_DATA_DIR" "$INPUT_DATA_DIR"

# Unique log file per sweep run so reruns never intermix
SWEEP_TS="$(date +%Y%m%d_%H%M%S)"
SWEEP_LOG="${LOG_DIR}/volmix_sweep_${SWEEP_TS}.log"

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
echo "Run data  → ${VOLMIX_RUNS}/"
echo "SLURM log → ${SWEEP_LOG}"
echo "======================================"

for (( i=FROM; i<END; i++ )); do
    P="${PRESSURES[$i]}"
    IS_LAST_IN_BATCH=$([ "$i" -eq "$((END-1))" ] && echo "yes" || echo "no")

    DATANAME="${BASE_DATANAME}_pstar${P}"
    TOTSTEPS=$SLAB_STEPS

    # --- Adaptive rank packing -------------------------------------------
    # High P* (iso-NPT) compresses the box: more atoms + fatter ghost shells
    # per subdomain, which OOM-kills the densest rank at 128/node. Halve to
    # 64/node above P*=1.5 (4x64=256 ranks for slab) so each rank gets ~2x
    # RAM. Node count is unchanged, so inter-node communication stays put.
    TPN=$(awk -v p="$P" 'BEGIN{print (p>=1.6)?64:128}')

    ISOLATED_DATANAME="isolated_${DATANAME}_${INTERACTION}_${TOTSTEPS}"
    SOL_DATANAME="${ISOLATED_DATANAME}_solvent_only"
    POL_DATANAME="${ISOLATED_DATANAME}_polymer_only"

    # Symlink pstar-specific DATANAME → base file so run_lammps.sh finds it in input_data/
    SRC_DATA="${INPUT_DATA_DIR}/${BASE_DATANAME}.data"
    LNK_DATA="${INPUT_DATA_DIR}/${DATANAME}.data"
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
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=${TPN}
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=5:00:00
#SBATCH --output=${SWEEP_LOG}

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

echo ""; echo "====== slab_p${P} | \$(date) | \$(hostname) ======"
export SKIP_WIDOM=1
export LAMMPS_RUNS_OVERRIDE="${VOLMIX_RUNS}"

# Retry loop — transient node OOM failures are retried up to 3 times.
# Each attempt creates a fresh timestamped work dir; the glob below picks the latest.
MAX_ATTEMPTS=6
for attempt in \$(seq 1 \$MAX_ATTEMPTS); do
    echo "slab_p${P}: attempt \$attempt of \$MAX_ATTEMPTS"
    if bash "${SCRIPTS_DIR}/run_lammps.sh" "slab_with_support" "${DATANAME}" "${INTERACTION}" \
            "${SLAB_STEPS}" "0" "" "${P}"; then
        echo "slab_p${P}: succeeded on attempt \$attempt"
        break
    fi
    if [ "\$attempt" -eq "\$MAX_ATTEMPTS" ]; then
        echo "ERROR: slab_p${P} failed after \$MAX_ATTEMPTS attempts"
        exit 1
    fi
    echo "slab_p${P}: attempt \$attempt failed, retrying in 120s..."
    sleep 120
done

WORK_DIR=\$(ls -dt "${VOLMIX_RUNS}/slab_with_support_${DATANAME}_${INTERACTION}_"* 2>/dev/null | head -1)
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
    # Job 2: isolate_gel.py  (depends on slab)
    #   Strips bath solvent, support, piston → isolated_*.data
    #   Output does NOT overwrite final_config_*.data
    # ------------------------------------------------------------------
    ISOLATE_BATCH=$(mktemp /tmp/isolate_volmix_p${P}_XXXX.batch)
    cat > "$ISOLATE_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=isolate_p${P}
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --output=${SWEEP_LOG}

set -euo pipefail
module reset
module load gcc/10.2.0
module load python/3.8.12

echo ""; echo "====== isolate_p${P} | \$(date) | \$(hostname) ======"
FINAL_CONFIG="${SLAB_DATA_DIR}/final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}.data"
ISOLATED_OUT="${SLAB_DATA_DIR}/${ISOLATED_DATANAME}.data"

if [ ! -f "\$FINAL_CONFIG" ]; then
    echo "ERROR: final_config not found: \$FINAL_CONFIG"
    exit 1
fi
if [ -f "\$ISOLATED_OUT" ]; then
    echo "WARNING: isolated file already exists, skipping: \$ISOLATED_OUT"
    exit 0
fi

python3 "${SCRIPTS_DIR}/isolate_gel.py" \
    --input  "\$FINAL_CONFIG" \
    --output "\$ISOLATED_OUT"

echo "Isolated gel written: \$ISOLATED_OUT"
echo "\$ISOLATED_OUT" > "${MANIFEST_DIR}/p${P}_isolated.path"
HEREDOC

    ISOLATE_JID=$(sbatch --parsable --dependency=afterok:${SLAB_JID} "$ISOLATE_BATCH")
    echo "P=${P}: submitted isolate_gel          JID=${ISOLATE_JID} (after ${SLAB_JID})"

    # ------------------------------------------------------------------
    # Job 3: split_gel.py  (depends on isolate; uses isolated_*.data)
    #   Writes isolated_*_polymer_only.data and isolated_*_solvent_only.data
    #   Does NOT touch final_config_*.data
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
#SBATCH --output=${SWEEP_LOG}

module reset
module load gcc/10.2.0
module load python/3.8.12

echo ""; echo "====== split_p${P} | \$(date) | \$(hostname) ======"
ISOLATED_FILE="${SLAB_DATA_DIR}/${ISOLATED_DATANAME}.data"
if [ ! -f "\$ISOLATED_FILE" ]; then
    echo "ERROR: Isolated file not found: \$ISOLATED_FILE"
    exit 1
fi

python3 "${SCRIPTS_DIR}/split_gel.py" "\$ISOLATED_FILE" \
    --polymer-dir "${INPUT_DATA_DIR}" \
    --solvent-dir "${INPUT_DATA_DIR}"

echo "Split complete:"
ls -lh "${INPUT_DATA_DIR}/${ISOLATED_DATANAME}_polymer_only.data"
ls -lh "${INPUT_DATA_DIR}/${ISOLATED_DATANAME}_solvent_only.data"
HEREDOC

    SPLIT_JID=$(sbatch --parsable --dependency=afterok:${ISOLATE_JID} "$SPLIT_BATCH")
    echo "P=${P}: submitted split_gel           JID=${SPLIT_JID}  (after ${ISOLATE_JID})"

    # ------------------------------------------------------------------
    # Job 4a: solvent_pure  (depends on split)
    # ------------------------------------------------------------------
    SOL_BATCH=$(mktemp /tmp/solvent_volmix_p${P}_XXXX.batch)
    cat > "$SOL_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=sol_p${P}
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${TPN}
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=1:00:00
#SBATCH --output=${SWEEP_LOG}

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

echo ""; echo "====== sol_p${P} | \$(date) | \$(hostname) ======"
export LAMMPS_RUNS_OVERRIDE="${VOLMIX_RUNS}"

bash "${SCRIPTS_DIR}/run_lammps.sh" "solvent_pure" "${SOL_DATANAME}" "${PURE_INTERACTION}" \
    "${PURE_STEPS}" "0" "" "${P}"

WORK_DIR=\$(ls -dt "${VOLMIX_RUNS}/solvent_pure_${SOL_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_solvent.workdir"
HEREDOC

    SOL_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} "$SOL_BATCH")
    echo "P=${P}: submitted solvent_pure         JID=${SOL_JID}   (after ${SPLIT_JID})"

    # ------------------------------------------------------------------
    # Job 4b: polymer_pure  (depends on split, concurrent with solvent_pure)
    # ------------------------------------------------------------------
    POL_BATCH=$(mktemp /tmp/polymer_volmix_p${P}_XXXX.batch)
    NEXT_FROM=$END

    cat > "$POL_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=pol_p${P}
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${TPN}
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=5:00:00
#SBATCH --output=${SWEEP_LOG}

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

echo ""; echo "====== pol_p${P} | \$(date) | \$(hostname) ======"
export LAMMPS_RUNS_OVERRIDE="${VOLMIX_RUNS}"

bash "${SCRIPTS_DIR}/run_lammps.sh" "polymer_pure" "${POL_DATANAME}" "${PURE_INTERACTION}" \
    "${PURE_STEPS}" "0" "" "${P}"

WORK_DIR=\$(ls -dt "${VOLMIX_RUNS}/polymer_pure_${POL_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_polymer.workdir"
HEREDOC

    POL_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} "$POL_BATCH")
    echo "P=${P}: submitted polymer_pure         JID=${POL_JID}   (after ${SPLIT_JID})"

    # Submit launcher for next batch (from login node, avoids sbatch-from-compute issues)
    if [ "$IS_LAST_IN_BATCH" = "yes" ] && [ "$NEXT_FROM" -lt "$TOTAL" ]; then
        LAUNCH_BATCH=$(mktemp /tmp/volmix_launch_XXXX.batch)
        cat > "$LAUNCH_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=volmix_launch
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=0:05:00
#SBATCH --output=${SWEEP_LOG}

echo ""; echo "====== launcher --from ${NEXT_FROM} | \$(date) | \$(hostname) ======"
module reset
bash "${SCRIPT_DIR}/volmix_sweep.sh" --from ${NEXT_FROM}
HEREDOC
        LAUNCH_JID=$(sbatch --parsable --dependency=afterok:${POL_JID} "$LAUNCH_BATCH")
        echo "P=${P}: submitted next-batch launcher JID=${LAUNCH_JID}  (after ${POL_JID}, --from ${NEXT_FROM})"
    fi

    echo "--------------------------------------"
done

echo "======================================"
echo "Batch submitted: ${PRESSURES[*]:$FROM:$((END-FROM))}"
if [ "$END" -lt "$TOTAL" ]; then
    echo "Next batch (indices ${END}+) will auto-submit when this batch completes."
fi
echo "Monitor with: squeue -u \$USER"
echo "Run data: ${VOLMIX_RUNS}/"
echo "======================================"
