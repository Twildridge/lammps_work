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
# SLURM log  → simulations/volmix_sweep/volmix_sweep.log  (this directory)
# Run data   → ~/Documents/lammps_runs/volmix_sweep/
# Manifests  → ~/Documents/lammps_runs/volmix_sweep/sweep_manifest/
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$(cd "${SCRIPT_DIR}/../../scripts" && pwd)"
# SCRIPTS_DIR = lammps_work/scripts/ (two levels up from simulations/volmix_sweep/)

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
echo "SLURM log → ${LOG_DIR}/volmix_sweep.log"
echo "======================================"

for (( i=FROM; i<END; i++ )); do
    P="${PRESSURES[$i]}"
    IS_LAST_IN_BATCH=$([ "$i" -eq "$((END-1))" ] && echo "yes" || echo "no")

    DATANAME="${BASE_DATANAME}_pstar${P}"
    TOTSTEPS=$SLAB_STEPS
    SOL_DATANAME="final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_solvent_only"
    POL_DATANAME="final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_polymer_only"

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
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=5:00:00
#SBATCH --output=${LOG_DIR}/volmix_sweep.log
#SBATCH --open-mode=append

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

cd "${SCRIPTS_DIR}" || { echo "ERROR: cd to SCRIPTS_DIR failed"; exit 1; }

# Retry loop — transient node OOM failures are retried up to 3 times.
# Each attempt creates a fresh timestamped work dir; the glob below picks the latest.
MAX_ATTEMPTS=3
for attempt in \$(seq 1 \$MAX_ATTEMPTS); do
    echo "slab_p${P}: attempt \$attempt of \$MAX_ATTEMPTS"
    if ./run_lammps.sh "slab_with_support" "${DATANAME}" "${INTERACTION}" \
            "${SLAB_STEPS}" "0" "" "${P}"; then
        echo "slab_p${P}: succeeded on attempt \$attempt"
        break
    fi
    if [ "\$attempt" -eq "\$MAX_ATTEMPTS" ]; then
        echo "ERROR: slab_p${P} failed after \$MAX_ATTEMPTS attempts"
        exit 1
    fi
    echo "slab_p${P}: attempt \$attempt failed, retrying in 60s..."
    sleep 60
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
#SBATCH --output=${LOG_DIR}/volmix_sweep.log
#SBATCH --open-mode=append

module reset
module load gcc/10.2.0
module load python/3.8.12

echo ""; echo "====== split_p${P} | \$(date) | \$(hostname) ======"
INPUT_FILE="${SLAB_DATA_DIR}/final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}.data"
if [ ! -f "\$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: \$INPUT_FILE"
    exit 1
fi

python3 "${SCRIPTS_DIR}/split_gel.py" "\$INPUT_FILE" \
    --polymer-dir "${INPUT_DATA_DIR}" \
    --solvent-dir "${INPUT_DATA_DIR}"

echo "Split complete:"
ls -lh "${INPUT_DATA_DIR}/final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_polymer_only.data"
ls -lh "${INPUT_DATA_DIR}/final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_solvent_only.data"
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
#SBATCH --output=${LOG_DIR}/volmix_sweep.log
#SBATCH --open-mode=append

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

cd "${SCRIPTS_DIR}" || { echo "ERROR: cd to SCRIPTS_DIR failed"; exit 1; }
./run_lammps.sh "solvent_pure" "${SOL_DATANAME}" "${PURE_INTERACTION}" \
    "${PURE_STEPS}" "0" "" "${P}"

WORK_DIR=\$(ls -dt "${VOLMIX_RUNS}/solvent_pure_${SOL_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_solvent.workdir"
HEREDOC

    SOL_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} "$SOL_BATCH")
    echo "P=${P}: submitted solvent_pure         JID=${SOL_JID}   (after ${SPLIT_JID})"

    # ------------------------------------------------------------------
    # Job 3b: polymer_pure  (depends on split, concurrent with solvent_pure)
    # ------------------------------------------------------------------
    POL_BATCH=$(mktemp /tmp/polymer_volmix_p${P}_XXXX.batch)
    NEXT_FROM=$END

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
#SBATCH --output=${LOG_DIR}/volmix_sweep.log
#SBATCH --open-mode=append

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

cd "${SCRIPTS_DIR}" || { echo "ERROR: cd to SCRIPTS_DIR failed"; exit 1; }
./run_lammps.sh "polymer_pure" "${POL_DATANAME}" "${PURE_INTERACTION}" \
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
#SBATCH --partition=shared
#SBATCH --account=csb197
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=0:05:00
#SBATCH --output=${LOG_DIR}/volmix_sweep.log
#SBATCH --open-mode=append

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
