#!/usr/bin/env bash
# =============================================================================
# volmix_sweep.sh
# =============================================================================
# Submits a chained SLURM pipeline for slab_with_support at 11 pressures
# P* = 1.0, 1.1, 1.2, ..., 2.0 (0.1 increments).
#
# Pipeline per pressure:
#   [1] slab_with_support   (600k steps, NPT at P*)
#         → writes N_SNAPS snapshot configs at equal intervals during production
#         → copies all snap files to lammps_data/slab_with_support/
#   For each snapshot snap=1..N_SNAPS:
#   [2] isolate_gel.py      (strips bath solvent, support, piston → isolated_*_snap<N>.data)
#   [3] split_gel.py        → polymer_only and solvent_only for this snapshot
#   For each replica rep=1..NREPS (concurrent per snapshot):
#   [3a] solvent_pure       (100k steps, NPT at same P*)
#   [3b] polymer_pure       (100k steps, NPT at same P*)  concurrent with 3a
#   [3c] mixed NPT          (100k steps, NPT at same P*)  concurrent with 3a/b
#
# This nested design gives N_SNAPS × NREPS ΔV_mix measurements per pressure,
# capturing both snapshot-to-snapshot composition variance (from different gel
# swelling states in the slab) AND thermal (NPT-replica) variance.
#
# Naming convention:
#   Slab DATANAME:     ${BASE_DATANAME}_pstar${P}           (no rep suffix)
#   Snapshot config:   final_config_..._snap${snap}.data
#   Isolated stem:     isolated_..._snap${snap}
#   NPT run dataname:  ${isolated_stem}_rep${rep}            (symlink in INPUT_DATA_DIR)
#   Manifest keys:     p${P}_slab.workdir
#                      p${P}_snap${snap}_isolated.path
#                      p${P}_snap${snap}_rep${rep}_mixed.workdir
#                      p${P}_snap${snap}_rep${rep}_solvent.workdir
#                      p${P}_snap${snap}_rep${rep}_polymer.workdir
#
# Usage:
#   bash volmix_sweep.sh                    # full pipeline from P=1.0
#   bash volmix_sweep.sh --skip-slab        # skip slab runs (snap files must exist)
#   bash volmix_sweep.sh --from 6           # resume from index 6 (P=1.6)
#   bash volmix_sweep.sh --from 6 --skip-slab
#   bash volmix_sweep.sh --only 1.3         # smoke-test single pressure
#
# SLURM log  → simulations/volmix_sweep/volmix_sweep_YYYYMMDD_HHMMSS.log
# Run data   → ~/Documents/lammps_runs/volmix_sweep/
# Manifests  → ~/Documents/lammps_runs/volmix_sweep/sweep_manifest/
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$(cd "${SCRIPT_DIR}/../../scripts" 2>/dev/null && pwd || echo "${SCRIPT_DIR}")"

BASE_DATANAME="slab_support_5beads_tall_rho04"
INTERACTION="1.0_1.0"
PURE_INTERACTION="1.0_0.0"
SLAB_STEPS=600000
PURE_STEPS=100000

LAMMPS_DATA="$HOME/Documents/lammps_data"
LAMMPS_RUNS="$HOME/Documents/lammps_runs"
INPUT_DATA_DIR="${LAMMPS_DATA}/input_data"
SLAB_DATA_DIR="${LAMMPS_DATA}/slab_with_support"

VOLMIX_RUNS="${LAMMPS_RUNS}/volmix_sweep"
MANIFEST_DIR="${VOLMIX_RUNS}/sweep_manifest"
LOG_DIR="${SCRIPT_DIR}"

mkdir -p "$VOLMIX_RUNS" "$MANIFEST_DIR" "$SLAB_DATA_DIR" "$INPUT_DATA_DIR"

SWEEP_TS="$(date +%Y%m%d_%H%M%S)"
SWEEP_LOG="${LOG_DIR}/volmix_sweep_${SWEEP_TS}.log"

if ! command -v sbatch &>/dev/null; then
    echo "ERROR: sbatch not found. Are you on the Expanse login node?"
    exit 1
fi

# Parse arguments
FROM=0
SKIP_SLAB=0
ONLY=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)      FROM="$2"; shift 2 ;;
        --skip-slab) SKIP_SLAB=1; shift ;;
        --only)      ONLY="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

PRESSURES=(1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)

if [ -n "$ONLY" ]; then
    PRESSURES=("$ONLY")
    FROM=0
    echo ">>> --only ${ONLY}: single-pressure smoke test (auto-chaining disabled)"
fi

WINDOW=1
N_SNAPS=3   # snapshots extracted from a single slab production run
NREPS=3     # NPT thermal replicas per snapshot  →  9 ΔV_mix samples per pressure
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
echo "Snapshots per slab:  ${N_SNAPS}"
echo "NPT reps per snap:   ${NREPS}"
echo "ΔV_mix samples/P*:  $((N_SNAPS * NREPS))"
echo "Run data  → ${VOLMIX_RUNS}/"
echo "SLURM log → ${SWEEP_LOG}"
echo "======================================"

for (( i=FROM; i<END; i++ )); do
    P="${PRESSURES[$i]}"
    IS_LAST_IN_BATCH=$([ "$i" -eq "$((END-1))" ] && echo "yes" || echo "no")
    NEXT_FROM=$END

    TPN=$(awk -v p="$P" 'BEGIN{print (p>=1.6)?64:128}')

    # One DATANAME per pressure — slab runs once (no rep suffix)
    DATANAME="${BASE_DATANAME}_pstar${P}"
    TOTSTEPS=$SLAB_STEPS

    # Symlink once per pressure so run_lammps.sh finds the base data file
    SRC_DATA="${INPUT_DATA_DIR}/${BASE_DATANAME}.data"
    LNK_DATA="${INPUT_DATA_DIR}/${DATANAME}.data"
    if [ ! -e "$LNK_DATA" ] && [ ! -L "$LNK_DATA" ]; then
        ln -s "$SRC_DATA" "$LNK_DATA"
    fi

    # ------------------------------------------------------------------
    # Job 1: slab_with_support  (runs ONCE per pressure)
    #   Writes N_SNAPS snapshot configs during production:
    #     final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_snap1.data
    #     final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_snap2.data
    #     final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_snap3.data
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
export LAMMPS_RUNS_OVERRIDE="${VOLMIX_RUNS}"

MAX_ATTEMPTS=6
for attempt in \$(seq 1 \$MAX_ATTEMPTS); do
    echo "slab_p${P}: attempt \$attempt of \$MAX_ATTEMPTS"
    if bash "${SCRIPTS_DIR}/run_lammps.sh" "slab_with_support" "${DATANAME}" "${INTERACTION}" \
            "${SLAB_STEPS}" "0" "" "${P}" "12345"; then
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

# Copy all ${N_SNAPS} snapshot files to SLAB_DATA_DIR
for (( snap=1; snap<=${N_SNAPS}; snap++ )); do
    SNAP_FILE="\${WORK_DIR}/final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_snap\${snap}.data"
    if [ ! -f "\$SNAP_FILE" ]; then
        echo "ERROR: Snapshot \$snap not found: \$SNAP_FILE"
        exit 1
    fi
    cp "\$SNAP_FILE" "${SLAB_DATA_DIR}/"
    echo "Copied snap\${snap}: \$SNAP_FILE"
done

echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_slab.workdir"
HEREDOC

    if [ "$SKIP_SLAB" = "1" ]; then
        # Verify all N_SNAPS snapshot files exist
        all_snaps_ok=1
        for snap in $(seq 1 $N_SNAPS); do
            SNAP_CHECK="${SLAB_DATA_DIR}/final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_snap${snap}.data"
            if [ ! -f "$SNAP_CHECK" ]; then
                echo "ERROR: --skip-slab set but snapshot not found: $SNAP_CHECK"
                all_snaps_ok=0
            fi
        done
        if [ "$all_snaps_ok" = "0" ]; then exit 1; fi
        echo "P=${P}: --skip-slab: all ${N_SNAPS} snapshots found"
        SLAB_DEP=""
        rm -f "$SLAB_BATCH"
    else
        SLAB_JID=$(sbatch --parsable "$SLAB_BATCH")
        echo "P=${P}: submitted slab_with_support    JID=${SLAB_JID}"
        SLAB_DEP="--dependency=afterok:${SLAB_JID}"
    fi

    LAST_POL_JID=""

    # ------------------------------------------------------------------
    # Snapshot loop — one isolate+split per snapshot, then NREPS NPT sets
    # ------------------------------------------------------------------
    for (( snap=1; snap<=N_SNAPS; snap++ )); do

    SNAP_CONFIG="final_config_${DATANAME}_${INTERACTION}_${TOTSTEPS}_snap${snap}.data"
    ISOLATED_DATANAME="isolated_${DATANAME}_${INTERACTION}_${TOTSTEPS}_snap${snap}"
    SOL_DATANAME="${ISOLATED_DATANAME}_solvent_only"
    POL_DATANAME="${ISOLATED_DATANAME}_polymer_only"

    # ------------------------------------------------------------------
    # Job 2: isolate_gel.py  (once per snapshot)
    # ------------------------------------------------------------------
    ISOLATE_BATCH=$(mktemp /tmp/isolate_volmix_p${P}_s${snap}_XXXX.batch)
    cat > "$ISOLATE_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=isolate_p${P}_s${snap}
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

echo ""; echo "====== isolate_p${P}_s${snap} | \$(date) | \$(hostname) ======"
SNAP_FILE="${SLAB_DATA_DIR}/${SNAP_CONFIG}"
ISOLATED_OUT="${SLAB_DATA_DIR}/${ISOLATED_DATANAME}.data"

if [ ! -f "\$SNAP_FILE" ]; then
    echo "ERROR: snapshot not found: \$SNAP_FILE"
    exit 1
fi
if [ -f "\$ISOLATED_OUT" ]; then
    echo "WARNING: isolated file already exists, skipping: \$ISOLATED_OUT"
    echo "\$ISOLATED_OUT" > "${MANIFEST_DIR}/p${P}_snap${snap}_isolated.path"
    exit 0
fi

python3 "${SCRIPTS_DIR}/isolate_gel.py" \
    --input  "\$SNAP_FILE" \
    --output "\$ISOLATED_OUT"

echo "Isolated gel written: \$ISOLATED_OUT"
echo "\$ISOLATED_OUT" > "${MANIFEST_DIR}/p${P}_snap${snap}_isolated.path"
HEREDOC

    ISOLATE_JID=$(sbatch --parsable ${SLAB_DEP} "$ISOLATE_BATCH")
    echo "P=${P} snap=${snap}: submitted isolate_gel     JID=${ISOLATE_JID} (dep: ${SLAB_DEP:-none})"

    # ------------------------------------------------------------------
    # Job 3: split_gel.py  (once per snapshot)
    # ------------------------------------------------------------------
    SPLIT_BATCH=$(mktemp /tmp/split_volmix_p${P}_s${snap}_XXXX.batch)
    cat > "$SPLIT_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=split_p${P}_s${snap}
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

echo ""; echo "====== split_p${P}_s${snap} | \$(date) | \$(hostname) ======"
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
    echo "P=${P} snap=${snap}: submitted split_gel       JID=${SPLIT_JID}  (after ${ISOLATE_JID})"

    # ------------------------------------------------------------------
    # Replica loop — NREPS independent NPT runs from the same snapshot
    # VEL_SEED unique across all (snap, rep) combos:
    #   snap1 rep1=11111  rep2=22222  rep3=33333
    #   snap2 rep1=44444  rep2=55555  rep3=66666
    #   snap3 rep1=77777  rep2=88888  rep3=99999
    # ------------------------------------------------------------------
    for (( rep=1; rep<=NREPS; rep++ )); do

    VEL_SEED=$(( ((snap-1)*NREPS + rep) * 11111 ))

    # Rep-tagged datanames for unique work directories per (snap, rep)
    MIX_DATANAME="${ISOLATED_DATANAME}_rep${rep}"
    SOL_RUN_DATANAME="${SOL_DATANAME}_rep${rep}"
    POL_RUN_DATANAME="${POL_DATANAME}_rep${rep}"

    # ------------------------------------------------------------------
    # Job 3a: mixed-gel NPT  (depends on isolate; concurrent with sol/pol)
    # ------------------------------------------------------------------
    MIX_BATCH=$(mktemp /tmp/mix_volmix_p${P}_s${snap}_r${rep}_XXXX.batch)
    cat > "$MIX_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=mix_p${P}_s${snap}_r${rep}
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

echo ""; echo "====== mix_p${P}_s${snap}_r${rep} | \$(date) | \$(hostname) ======"
export LAMMPS_RUNS_OVERRIDE="${VOLMIX_RUNS}"

# Symlink isolated gel → rep-tagged name so work dirs are unique per (snap, rep)
ISO_SRC="${SLAB_DATA_DIR}/${ISOLATED_DATANAME}.data"
ISO_LNK="${INPUT_DATA_DIR}/${MIX_DATANAME}.data"
if [ ! -f "\$ISO_SRC" ]; then echo "ERROR: isolated gel not found: \$ISO_SRC"; exit 1; fi
[ -e "\$ISO_LNK" ] || ln -s "\$ISO_SRC" "\$ISO_LNK"

bash "${SCRIPTS_DIR}/run_lammps.sh" "polymer_pure" "${MIX_DATANAME}" "${INTERACTION}" \
    "${PURE_STEPS}" "0" "" "${P}" "${VEL_SEED}"

WORK_DIR=\$(ls -dt "${VOLMIX_RUNS}/polymer_pure_${MIX_DATANAME}_${INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_snap${snap}_rep${rep}_mixed.workdir"
HEREDOC

    MIX_JID=$(sbatch --parsable --dependency=afterok:${ISOLATE_JID} "$MIX_BATCH")
    echo "P=${P} snap=${snap} rep=${rep}: submitted mixed NPT   JID=${MIX_JID}  (after ${ISOLATE_JID})"

    # ------------------------------------------------------------------
    # Job 3b: solvent_pure  (depends on split; concurrent with mix/pol)
    # ------------------------------------------------------------------
    SOL_BATCH=$(mktemp /tmp/solvent_volmix_p${P}_s${snap}_r${rep}_XXXX.batch)
    cat > "$SOL_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=sol_p${P}_s${snap}_r${rep}
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

echo ""; echo "====== sol_p${P}_s${snap}_r${rep} | \$(date) | \$(hostname) ======"
export LAMMPS_RUNS_OVERRIDE="${VOLMIX_RUNS}"

# Symlink solvent_only → rep-tagged name
SOL_SRC="${INPUT_DATA_DIR}/${SOL_DATANAME}.data"
SOL_LNK="${INPUT_DATA_DIR}/${SOL_RUN_DATANAME}.data"
if [ ! -f "\$SOL_SRC" ]; then echo "ERROR: solvent_only not found: \$SOL_SRC"; exit 1; fi
[ -e "\$SOL_LNK" ] || ln -s "\$SOL_SRC" "\$SOL_LNK"

bash "${SCRIPTS_DIR}/run_lammps.sh" "solvent_pure" "${SOL_RUN_DATANAME}" "${PURE_INTERACTION}" \
    "${PURE_STEPS}" "0" "" "${P}" "${VEL_SEED}"

WORK_DIR=\$(ls -dt "${VOLMIX_RUNS}/solvent_pure_${SOL_RUN_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_snap${snap}_rep${rep}_solvent.workdir"
HEREDOC

    SOL_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} "$SOL_BATCH")
    echo "P=${P} snap=${snap} rep=${rep}: submitted solvent NPT  JID=${SOL_JID}   (after ${SPLIT_JID})"

    # ------------------------------------------------------------------
    # Job 3c: polymer_pure  (depends on split; concurrent with mix/sol)
    # ------------------------------------------------------------------
    POL_BATCH=$(mktemp /tmp/polymer_volmix_p${P}_s${snap}_r${rep}_XXXX.batch)
    cat > "$POL_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=pol_p${P}_s${snap}_r${rep}
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

echo ""; echo "====== pol_p${P}_s${snap}_r${rep} | \$(date) | \$(hostname) ======"
export LAMMPS_RUNS_OVERRIDE="${VOLMIX_RUNS}"

# Symlink polymer_only → rep-tagged name
POL_SRC="${INPUT_DATA_DIR}/${POL_DATANAME}.data"
POL_LNK="${INPUT_DATA_DIR}/${POL_RUN_DATANAME}.data"
if [ ! -f "\$POL_SRC" ]; then echo "ERROR: polymer_only not found: \$POL_SRC"; exit 1; fi
[ -e "\$POL_LNK" ] || ln -s "\$POL_SRC" "\$POL_LNK"

bash "${SCRIPTS_DIR}/run_lammps.sh" "polymer_pure" "${POL_RUN_DATANAME}" "${PURE_INTERACTION}" \
    "${PURE_STEPS}" "0" "" "${P}" "${VEL_SEED}"

WORK_DIR=\$(ls -dt "${VOLMIX_RUNS}/polymer_pure_${POL_RUN_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_snap${snap}_rep${rep}_polymer.workdir"
HEREDOC

    POL_JID=$(sbatch --parsable --dependency=afterok:${SPLIT_JID} "$POL_BATCH")
    echo "P=${P} snap=${snap} rep=${rep}: submitted polymer NPT  JID=${POL_JID}   (after ${SPLIT_JID})"

    LAST_POL_JID="${POL_JID}"

    done  # end replica loop

    done  # end snapshot loop

    # Submit launcher for next batch — depends on last polymer job so all
    # snap×rep combos finish before the next pressure starts.
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
SKIP_FLAG=$([ "${SKIP_SLAB}" = "1" ] && echo "--skip-slab" || echo "")
bash "${SCRIPT_DIR}/volmix_sweep.sh" --from ${NEXT_FROM} \${SKIP_FLAG}
HEREDOC
        LAUNCH_JID=$(sbatch --parsable --dependency=afterok:${LAST_POL_JID} "$LAUNCH_BATCH")
        echo "P=${P}: submitted next-batch launcher JID=${LAUNCH_JID}  (after ${LAST_POL_JID}, --from ${NEXT_FROM})"
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
