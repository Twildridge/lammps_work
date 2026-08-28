#!/usr/bin/env bash
# =============================================================================
# calibration_sweep.sh
# =============================================================================
# Submits the chained SLURM pipeline for the volume-fraction CALIBRATION sweep
# (2026-08 plan). Supersedes volmix_sweep.sh (archived in lammps_work/archive/;
# its sbatch-chaining / manifest / resume machinery is inherited verbatim).
#
# Physics: homogeneous periodic mixed-gel NPT runs over a fixed composition
# grid (N_f set DIRECTLY by adjust_solvent.py — identical composition coverage
# at every P) × a pressure grid. Analysis fits ln<V> vs ln N_f per P →
# φ_f^th = ∂ln<V>/∂ln N_f, partial molar volumes, and the Voronoi calibration
# λ(φ_p, P). ΔV_mix(P) falls out per grid point via the scale=1 convention
# (per-loading pure-solvent companions hold EXACTLY the mixed box's atoms).
#
# Pipeline:
#   [prep] (ONCE, not per pressure; skipped with --skip-prep)
#       isolate_gel.py    BASE_SNAPSHOT → ISOLATED_STEM.data   (skip if exists)
#       for NF in NF_GRID:
#           adjust_solvent.py → ISOLATED_STEM_nf${NF}.data      (skip if exists)
#           split_gel.py      → ..._nf${NF}_solvent_only.data   (+ polymer_only)
#   per pressure P (chained, WINDOW pressures per batch):
#       polymer_pure companion (ONCE per P — N_p is loading-independent)
#       for NF in NF_GRID: for rep in 1..NREPS:
#           mixed NPT      (polymer_pure engine, aniso)   concurrent
#           solvent_pure companion (same N_f)             concurrent
#
# Naming convention:
#   Isolated stem:     isolated_<BASE_SNAPSHOT stem without final_config_>
#   Loading file:      ${ISOLATED_STEM}_nf${NF}.data
#   NPT run dataname:  ${ISOLATED_STEM}_nf${NF}_pstar${P}_rep${rep}   (symlink
#                      in INPUT_DATA_DIR; pure companions analogous)
#   Manifest keys:     prep_isolated.path
#                      prep_nf${NF}.path
#                      p${P}_polymer.workdir
#                      p${P}_nf${NF}_rep${rep}_mixed.workdir
#                      p${P}_nf${NF}_rep${rep}_solvent.workdir
#
# Usage:
#   bash calibration_sweep.sh                  # full sweep from the first pressure
#   bash calibration_sweep.sh --only 1.5       # single-pressure smoke test
#   bash calibration_sweep.sh --from 3         # resume from pressure index 3
#   bash calibration_sweep.sh --skip-prep      # prep artifacts already exist
#
# SLURM log  → simulations/calibration_sweep/calibration_sweep_YYYYMMDD_HHMMSS.log
# Run data   → ~/Documents/lammps_runs/calibration_sweep/
# Manifests  → ~/Documents/lammps_runs/calibration_sweep/sweep_manifest/
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$(cd "${SCRIPT_DIR}/../../scripts" 2>/dev/null && pwd || echo "${SCRIPT_DIR}")"

# =============================================================================
# GRIDS AND COUNTS — everything a sweep is made of lives HERE, nothing is
# buried below and nothing is hard-coded in any .lmp (all -var driven).
# =============================================================================

# Base snapshot: the 14M-step P=1.5 slab_with_support_periodic equilibration
# (must exist in SLAB_DATA_DIR on the cluster). NO per-pressure slab runs —
# NPT re-equilibration at each grid P erases the snapshot's origin.
BASE_SNAPSHOT="final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000000.data"

# Pressure grid (plan: 0.5 … 2.0, step 0.25 or finer).
PRESSURES=(0.50 0.75 1.00 1.25 1.50 1.75 2.00)

# Composition grid: exact solvent counts N_f, identical at every pressure.
# The base snapshot isolates to N_p = 104283, N_f = 53407 (equilibrium-swollen
# at P = 1.5). Defaults: 10 geometric steps (ratio ≈ 0.90) from N_f_eq down to
# ~0.39×N_f_eq. EDIT the lower end to match the max φ_p reached in the c0.40
# triaxial_compression profiles before production (validation needs the grid
# to span it).
NF_GRID=(53400 48100 43300 39000 35100 31600 28400 25600 23000 20800)

# NPT thermal replicas per (P, N_f) grid point (velocity-seed variation only —
# same input file via symlink).
NREPS=2

# Run lengths (steps). Closed boxes from an equilibrated snapshot — no
# swelling/diffusion mode. Budget ~100–300k (≤500k at low P); VERIFY from the
# <V> trace in calibration_analysis.ipynb, don't assume. Must be divisible by
# n_blocks*100 = 1000 (engine block averaging).
CALIB_STEPS=300000        # P above LOWP_MAX
CALIB_STEPS_LOWP=500000   # P at or below LOWP_MAX (slower V relaxation)
LOWP_MAX=0.75

# Calibration dump: engines write CALIB_FRAMES+1 all-atom frames near the run
# end, CALIB_DUMP_EVERY steps apart (Voronoi input; exported to run_lammps.sh).
CALIB_FRAMES=5
CALIB_DUMP_EVERY=2000

# Interaction labels (epsSS_epsSP): athermal WCA mix / pure conventions.
INTERACTION="1.0_1.0"
PURE_INTERACTION="1.0_0.0"

# Pressures per submission batch (auto-chained launcher submits the next batch).
WINDOW=1

# =============================================================================

LAMMPS_DATA="$HOME/Documents/lammps_data"
LAMMPS_RUNS="$HOME/Documents/lammps_runs"
INPUT_DATA_DIR="${LAMMPS_DATA}/input_data"
SLAB_DATA_DIR="${LAMMPS_DATA}/slab_with_support"

CALIB_RUNS="${LAMMPS_RUNS}/calibration_sweep"
MANIFEST_DIR="${CALIB_RUNS}/sweep_manifest"
LOG_DIR="${SCRIPT_DIR}"

mkdir -p "$CALIB_RUNS" "$MANIFEST_DIR" "$INPUT_DATA_DIR"

SWEEP_TS="$(date +%Y%m%d_%H%M%S)"
SWEEP_LOG="${LOG_DIR}/calibration_sweep_${SWEEP_TS}.log"

# Isolated stem derived from the snapshot name (final_config_X.data → isolated_X)
ISOLATED_STEM="isolated_${BASE_SNAPSHOT#final_config_}"
ISOLATED_STEM="${ISOLATED_STEM%.data}"

if ! command -v sbatch &>/dev/null; then
    echo "ERROR: sbatch not found. Are you on the Expanse login node?"
    exit 1
fi

# Parse arguments
FROM=0
SKIP_PREP=0
ONLY=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)      FROM="$2"; shift 2 ;;
        --skip-prep) SKIP_PREP=1; shift ;;
        --only)      ONLY="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [ -n "$ONLY" ]; then
    PRESSURES=("$ONLY")
    FROM=0
    echo ">>> --only ${ONLY}: single-pressure smoke test (auto-chaining disabled)"
fi

TOTAL=${#PRESSURES[@]}
END=$(( FROM + WINDOW ))
if [ "$END" -gt "$TOTAL" ]; then END=$TOTAL; fi

if [ "$FROM" -ge "$TOTAL" ]; then
    echo "All pressures already submitted."
    exit 0
fi

echo "======================================"
echo "Calibration sweep — batch submission"
echo "Submitting indices ${FROM}–$((END-1)): ${PRESSURES[*]:$FROM:$((END-FROM))}"
echo "Composition grid:   ${#NF_GRID[@]} loadings (${NF_GRID[0]} … ${NF_GRID[$((${#NF_GRID[@]}-1))]})"
echo "NPT reps per point: ${NREPS}"
echo "Runs per pressure:  $(( ${#NF_GRID[@]} * NREPS * 2 + 1 )) (mixed + solvent per point, polymer once)"
echo "Run data  → ${CALIB_RUNS}/"
echo "SLURM log → ${SWEEP_LOG}"
echo "======================================"

# Record the sweep configuration for the analysis notebook
{
    echo "date        ${SWEEP_TS}"
    echo "snapshot    ${BASE_SNAPSHOT}"
    echo "isolated    ${ISOLATED_STEM}"
    echo "pressures   ${PRESSURES[*]}"
    echo "nf_grid     ${NF_GRID[*]}"
    echo "nreps       ${NREPS}"
    echo "steps       ${CALIB_STEPS} (lowP ${CALIB_STEPS_LOWP} at P<=${LOWP_MAX})"
    echo "interaction ${INTERACTION} / ${PURE_INTERACTION}"
    echo "calib_dump  ${CALIB_FRAMES}+1 frames every ${CALIB_DUMP_EVERY}"
} > "${MANIFEST_DIR}/sweep_config_${SWEEP_TS}.txt"

# ------------------------------------------------------------------
# [prep] isolate + adjust + split (ONCE for the whole sweep)
# ------------------------------------------------------------------
PREP_DEP=""
if [ "$SKIP_PREP" = "1" ]; then
    prep_ok=1
    for NF in "${NF_GRID[@]}"; do
        for f in "${INPUT_DATA_DIR}/${ISOLATED_STEM}_nf${NF}.data" \
                 "${INPUT_DATA_DIR}/${ISOLATED_STEM}_nf${NF}_solvent_only.data"; do
            if [ ! -f "$f" ]; then
                echo "ERROR: --skip-prep set but prep artifact not found: $f"
                prep_ok=0
            fi
        done
    done
    if [ ! -f "${INPUT_DATA_DIR}/${ISOLATED_STEM}_polymer_only.data" ]; then
        echo "ERROR: --skip-prep set but polymer_only not found"
        prep_ok=0
    fi
    if [ "$prep_ok" = "0" ]; then exit 1; fi
    echo "prep: --skip-prep: all $(( ${#NF_GRID[@]} * 2 + 1 )) artifacts found"
else
    PREP_BATCH=$(mktemp "${TMPDIR:-/tmp}/calib_prep_XXXXXX.batch")
    cat > "$PREP_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=calib_prep
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=${SWEEP_LOG}

set -euo pipefail
module reset
module load gcc/10.2.0
module load python/3.8.12

echo ""; echo "====== calib_prep | \$(date) | \$(hostname) ======"

SNAP_FILE="${SLAB_DATA_DIR}/${BASE_SNAPSHOT}"
ISOLATED_OUT="${SLAB_DATA_DIR}/${ISOLATED_STEM}.data"

if [ ! -f "\$SNAP_FILE" ]; then
    echo "ERROR: base snapshot not found: \$SNAP_FILE"
    exit 1
fi

# isolate (skip if present)
if [ -f "\$ISOLATED_OUT" ]; then
    echo "isolated file exists, skipping: \$ISOLATED_OUT"
else
    python3 "${SCRIPTS_DIR}/isolate_gel.py" --input "\$SNAP_FILE" --output "\$ISOLATED_OUT"
fi
echo "\$ISOLATED_OUT" > "${MANIFEST_DIR}/prep_isolated.path"

# adjust + split per loading (skip-if-present per artifact)
for NF in ${NF_GRID[*]}; do
    NF_OUT="${INPUT_DATA_DIR}/${ISOLATED_STEM}_nf\${NF}.data"
    if [ -f "\$NF_OUT" ]; then
        echo "nf\${NF}: exists, skipping adjust"
    else
        python3 "${SCRIPTS_DIR}/adjust_solvent.py" \
            --input "\$ISOLATED_OUT" --target-nf "\$NF" --output "\$NF_OUT" --seed 12345
    fi
    if [ -f "${INPUT_DATA_DIR}/${ISOLATED_STEM}_nf\${NF}_solvent_only.data" ]; then
        echo "nf\${NF}: solvent_only exists, skipping split"
    else
        python3 "${SCRIPTS_DIR}/split_gel.py" "\$NF_OUT" --output-dir "${INPUT_DATA_DIR}"
    fi
    echo "\$NF_OUT" > "${MANIFEST_DIR}/prep_nf\${NF}.path"
done

# Canonical polymer_only (identical for every loading — N_p fixed): keep the
# first loading's copy under the loading-independent name.
POL_CANON="${INPUT_DATA_DIR}/${ISOLATED_STEM}_polymer_only.data"
if [ ! -f "\$POL_CANON" ]; then
    cp "${INPUT_DATA_DIR}/${ISOLATED_STEM}_nf${NF_GRID[0]}_polymer_only.data" "\$POL_CANON"
fi
echo "prep complete"
HEREDOC

    PREP_JID=$(sbatch --parsable "$PREP_BATCH")
    echo "prep: submitted isolate+adjust+split    JID=${PREP_JID}"
    PREP_DEP="--dependency=afterok:${PREP_JID}"
fi

# ------------------------------------------------------------------
# Pressure loop
# ------------------------------------------------------------------
for (( i=FROM; i<END; i++ )); do
    P="${PRESSURES[$i]}"
    IS_LAST_IN_BATCH=$([ "$i" -eq "$((END-1))" ] && echo "yes" || echo "no")
    NEXT_FROM=$END

    TPN=$(awk -v p="$P" 'BEGIN{print (p>=1.6)?64:128}')
    NSTEPS=$(awk -v p="$P" -v lo="$LOWP_MAX" -v a="$CALIB_STEPS_LOWP" -v b="$CALIB_STEPS" \
             'BEGIN{print (p<=lo)?a:b}')

    LAST_JID=""

    # ------------------------------------------------------------------
    # polymer_pure companion — ONCE per pressure (N_p loading-independent)
    # ------------------------------------------------------------------
    POL_RUN_DATANAME="${ISOLATED_STEM}_polymer_only_pstar${P}"
    POL_BATCH=$(mktemp "${TMPDIR:-/tmp}/calib_pol_p${P}_XXXXXX.batch")
    cat > "$POL_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=cpol_p${P}
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

echo ""; echo "====== cpol_p${P} | \$(date) | \$(hostname) ======"
export LAMMPS_RUNS_OVERRIDE="${CALIB_RUNS}"
export CALIB_FRAMES="${CALIB_FRAMES}" CALIB_DUMP_EVERY="${CALIB_DUMP_EVERY}"

POL_SRC="${INPUT_DATA_DIR}/${ISOLATED_STEM}_polymer_only.data"
POL_LNK="${INPUT_DATA_DIR}/${POL_RUN_DATANAME}.data"
if [ ! -f "\$POL_SRC" ]; then echo "ERROR: polymer_only not found: \$POL_SRC"; exit 1; fi
[ -e "\$POL_LNK" ] || ln -s "\$POL_SRC" "\$POL_LNK"

bash "${SCRIPTS_DIR}/run_lammps.sh" "polymer_pure" "${POL_RUN_DATANAME}" "${PURE_INTERACTION}" \
    "${NSTEPS}" "" "${P}" "12345"

WORK_DIR=\$(ls -dt "${CALIB_RUNS}/polymer_pure_${POL_RUN_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_polymer.workdir"
HEREDOC

    POL_JID=$(sbatch --parsable ${PREP_DEP} "$POL_BATCH")
    echo "P=${P}: submitted polymer companion   JID=${POL_JID}  (dep: ${PREP_DEP:-none})"
    LAST_JID="$POL_JID"

    # ------------------------------------------------------------------
    # Loading × replica loop — mixed + solvent companion, all concurrent
    # ------------------------------------------------------------------
    NF_IDX=0
    for NF in "${NF_GRID[@]}"; do
    for (( rep=1; rep<=NREPS; rep++ )); do

    VEL_SEED=$(( (NF_IDX*NREPS + rep) * 11111 ))

    MIX_RUN_DATANAME="${ISOLATED_STEM}_nf${NF}_pstar${P}_rep${rep}"
    SOL_RUN_DATANAME="${ISOLATED_STEM}_nf${NF}_solvent_only_pstar${P}_rep${rep}"

    # ---- mixed NPT (polymer_pure engine, aniso) ----
    MIX_BATCH=$(mktemp "${TMPDIR:-/tmp}/calib_mix_p${P}_nf${NF}_r${rep}_XXXXXX.batch")
    cat > "$MIX_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=cmix_p${P}_nf${NF}_r${rep}
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

echo ""; echo "====== cmix_p${P}_nf${NF}_r${rep} | \$(date) | \$(hostname) ======"
export LAMMPS_RUNS_OVERRIDE="${CALIB_RUNS}"
export CALIB_FRAMES="${CALIB_FRAMES}" CALIB_DUMP_EVERY="${CALIB_DUMP_EVERY}"

MIX_SRC="${INPUT_DATA_DIR}/${ISOLATED_STEM}_nf${NF}.data"
MIX_LNK="${INPUT_DATA_DIR}/${MIX_RUN_DATANAME}.data"
if [ ! -f "\$MIX_SRC" ]; then echo "ERROR: loading file not found: \$MIX_SRC"; exit 1; fi
[ -e "\$MIX_LNK" ] || ln -s "\$MIX_SRC" "\$MIX_LNK"

bash "${SCRIPTS_DIR}/run_lammps.sh" "polymer_pure" "${MIX_RUN_DATANAME}" "${INTERACTION}" \
    "${NSTEPS}" "" "${P}" "${VEL_SEED}"

WORK_DIR=\$(ls -dt "${CALIB_RUNS}/polymer_pure_${MIX_RUN_DATANAME}_${INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_nf${NF}_rep${rep}_mixed.workdir"
HEREDOC

    MIX_JID=$(sbatch --parsable ${PREP_DEP} "$MIX_BATCH")
    echo "P=${P} nf=${NF} rep=${rep}: submitted mixed NPT    JID=${MIX_JID}"

    # ---- solvent_pure companion (same N_f — scale=1 ΔV_mix reference) ----
    SOL_BATCH=$(mktemp "${TMPDIR:-/tmp}/calib_sol_p${P}_nf${NF}_r${rep}_XXXXXX.batch")
    cat > "$SOL_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=csol_p${P}_nf${NF}_r${rep}
#SBATCH --partition=compute
#SBATCH --account=csb197
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${TPN}
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=2:00:00
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

echo ""; echo "====== csol_p${P}_nf${NF}_r${rep} | \$(date) | \$(hostname) ======"
export LAMMPS_RUNS_OVERRIDE="${CALIB_RUNS}"
export CALIB_FRAMES="${CALIB_FRAMES}" CALIB_DUMP_EVERY="${CALIB_DUMP_EVERY}"

SOL_SRC="${INPUT_DATA_DIR}/${ISOLATED_STEM}_nf${NF}_solvent_only.data"
SOL_LNK="${INPUT_DATA_DIR}/${SOL_RUN_DATANAME}.data"
if [ ! -f "\$SOL_SRC" ]; then echo "ERROR: solvent_only not found: \$SOL_SRC"; exit 1; fi
[ -e "\$SOL_LNK" ] || ln -s "\$SOL_SRC" "\$SOL_LNK"

bash "${SCRIPTS_DIR}/run_lammps.sh" "solvent_pure" "${SOL_RUN_DATANAME}" "${PURE_INTERACTION}" \
    "${NSTEPS}" "" "${P}" "${VEL_SEED}"

WORK_DIR=\$(ls -dt "${CALIB_RUNS}/solvent_pure_${SOL_RUN_DATANAME}_${PURE_INTERACTION}_"* 2>/dev/null | head -1)
echo "\$WORK_DIR" > "${MANIFEST_DIR}/p${P}_nf${NF}_rep${rep}_solvent.workdir"
HEREDOC

    SOL_JID=$(sbatch --parsable ${PREP_DEP} "$SOL_BATCH")
    echo "P=${P} nf=${NF} rep=${rep}: submitted solvent NPT  JID=${SOL_JID}"

    LAST_JID="$SOL_JID"

    done  # end replica loop
    NF_IDX=$((NF_IDX+1))
    done  # end loading loop

    # Submit launcher for the next batch — depends on the last job submitted
    # for this pressure (inherited volmix pattern; jobs run concurrently, so
    # this is an approximation of after-all that has worked in practice).
    if [ "$IS_LAST_IN_BATCH" = "yes" ] && [ "$NEXT_FROM" -lt "$TOTAL" ]; then
        LAUNCH_BATCH=$(mktemp "${TMPDIR:-/tmp}/calib_launch_XXXXXX.batch")
        cat > "$LAUNCH_BATCH" << HEREDOC
#!/usr/bin/env bash
#SBATCH --job-name=calib_launch
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
bash "${SCRIPT_DIR}/calibration_sweep.sh" --from ${NEXT_FROM} --skip-prep
HEREDOC
        LAUNCH_JID=$(sbatch --parsable --dependency=afterok:${LAST_JID} "$LAUNCH_BATCH")
        echo "P=${P}: submitted next-batch launcher JID=${LAUNCH_JID}  (after ${LAST_JID}, --from ${NEXT_FROM})"
    fi

    echo "--------------------------------------"
done

echo "======================================"
echo "Batch submitted: ${PRESSURES[*]:$FROM:$((END-FROM))}"
if [ "$END" -lt "$TOTAL" ]; then
    echo "Next batch (indices ${END}+) will auto-submit when this batch completes."
fi
echo "Monitor with: squeue -u \$USER"
echo "Run data: ${CALIB_RUNS}/"
echo "======================================"
