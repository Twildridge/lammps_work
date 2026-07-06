#!/bin/bash
# ==============================================================================
# postprocess.sh — regenerate ALL analysis plots for a completed LAMMPS run.
#
# This is the single source of truth for post-processing. run_lammps.sh calls it
# automatically at the end of every run, and you can call it BY HAND on any
# finished run directory — no MD, no MPI, just the plotting scripts in the right
# order with the right module loaded. Use the manual path whenever a teardown
# hang or job kill skipped the automatic step (log ends at "Total wall time"
# with no post-processing banner) — all the output files are already on disk.
#
# Usage:
#   bash postprocess.sh <run_dir> <folder> <dataname> <interaction> <totsteps> [oldsteps] [press_target]
#
# Example (regenerate the non-periodic slab_with_support run):
#   bash scripts/postprocess.sh \
#     ~/Documents/lammps_runs/slab_with_support/slab_with_support_slab_support_5beads_tall_rho04_1.0_1.0_20260705_124556 \
#     slab_with_support slab_support_5beads_tall_rho04 1.0_1.0 3000000
#
# Optional env vars (default to the run_lammps.sh defaults):
#   SKIP_WIDOM=1   minimize/skip cavity_widom.py (slab_with_flow only)
#   STRAINS="..."  space-separated shear-strain list (shear_slab only)
# ==============================================================================
set -u

if [ $# -lt 5 ]; then
    echo "Usage: bash postprocess.sh <run_dir> <folder> <dataname> <interaction> <totsteps> [oldsteps] [press_target]"
    echo "Example:"
    echo "  bash postprocess.sh ~/Documents/lammps_runs/slab_with_support/<run_dir> \\"
    echo "       slab_with_support slab_support_5beads_tall_rho04 1.0_1.0 3000000"
    exit 1
fi

RUN_DIR=$1
FOLDER=$2
DATANAME=$3
INTERACTION=$4
TOTSTEPS=$5
OLDSTEPS=${6:-0}
PRESS_TARGET=${7:-1.5}
SKIP_WIDOM=${SKIP_WIDOM:-0}
STRAINS=${STRAINS:-0.1}

# Directory this script lives in (lammps_work/scripts/) — holds the plotters.
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Derive per-species epsilons from the interaction string (epsSS_epsSP).
IFS='_' read -r EPSSS EPSSP <<< "$INTERACTION"

# Plot-file stem shared by every analysis script.
STEM="${DATANAME}_${INTERACTION}_${TOTSTEPS}"

if [ ! -d "$RUN_DIR" ]; then
    echo "Error: run directory not found: $RUN_DIR"
    exit 1
fi
cd "$RUN_DIR" || { echo "Error: cannot cd into $RUN_DIR"; exit 1; }

echo "======================================"
echo "Running post-processing..."
echo "  run_dir     = $RUN_DIR"
echo "  folder      = $FOLDER"
echo "  stem        = $STEM"
echo "======================================"

# Load the Python stack used for plotting (same order as before).
source /etc/profile.d/modules.sh
module unload python/3.8.12 2>/dev/null
module load anaconda3/2021.05/q4munrg
python -c "import numpy; print(numpy.__version__)"

echo "Generating convergence plot..."
# Pass --p-ext for slab_with_flow so the pore-pressure panel uses the barostat target
if [ "$FOLDER" = "slab_with_flow" ]; then
    python "$SCRIPT_DIR/plot_lammps_log.py" "." "$STEM" --p-ext 1.8
else
    python "$SCRIPT_DIR/plot_lammps_log.py" "." "$STEM"
fi

# ── Cavity-biased Widom insertion (slab_with_flow only) ───────────────────────
if [ "$FOLDER" != "slab_with_flow" ] || [ "$SKIP_WIDOM" = "1" ]; then
    echo "Cavity Widom post-processing applies only to slab_with_flow (and not when SKIP_WIDOM=1) — skipping for ${FOLDER}."
else
    WIDOM_TRAJ="${RUN_DIR}/traj_files/widom_${STEM}.lammpstrj"
    echo "======================================"
    echo "Cavity Widom check:"
    echo "  Looking for: $WIDOM_TRAJ"
    echo "  traj_files/ contents:"
    ls -lh "${RUN_DIR}/traj_files/" 2>&1 | head -20
    echo "======================================"
    if [ -f "$WIDOM_TRAJ" ]; then
        echo "Running cavity-biased Widom insertion..."
        WIDOM_PEXT="1.8"
        WIDOM_EXCL="2.0"
        WIDOM_PISTON_EPS="0.0"

        WIDOM_EXCL_ARGS=()
        if [ -n "$WIDOM_EXCL" ]; then
            WIDOM_EXCL_ARGS=(--exclusion-buffer "$WIDOM_EXCL")
        fi

        python "$SCRIPT_DIR/cavity_widom.py" \
            --traj      "$WIDOM_TRAJ" \
            --out-dir   "output_files/chemical_potential" \
            --out-stem  "$STEM" \
            --eps-sp    "$EPSSP" \
            --eps-ss    "$EPSSS" \
            --n-bins    40 \
            --n-trial   50000 \
            --r-cavity  0.5 \
            --temperature 1.0 \
            --p-ext     "$WIDOM_PEXT" \
            --piston-eps "$WIDOM_PISTON_EPS" \
            "${WIDOM_EXCL_ARGS[@]}"

        echo "Re-generating convergence plot with cavity Widom panel..."
        python "$SCRIPT_DIR/plot_lammps_log.py" "." "$STEM" --p-ext "$WIDOM_PEXT"
    else
        echo "WARNING: Widom trajectory not found — skipping cavity_widom.py"
        echo "  Expected path: $WIDOM_TRAJ"
    fi
fi

# Pure solvent / polymer P-sweep: EOS plot instead of stress/piston scripts.
if [ "$FOLDER" = "solvent_phase" ] || [ "$FOLDER" = "polymer_phase" ]; then
    echo "${FOLDER} phase-sweep run detected — generating EOS plot..."
    python "$SCRIPT_DIR/plot_eos.py" "." "$DATANAME" "$INTERACTION"
    echo "======================================"
    echo "Done! Results are in: $RUN_DIR"
    echo "======================================"
    exit 0
fi

if [ "$FOLDER" = "shear_slab" ]; then
    echo "Generating shear stress-strain sweep plots (per strain: $STRAINS)..."
    python "$SCRIPT_DIR/plot_shear_strain_sweep.py" "." "$STEM" "$STRAINS"
else
    echo "Generating stress profiles..."
    python "$SCRIPT_DIR/plot_stress_profiles.py" "." "$STEM" "$OLDSTEPS"

    echo "Generating piston plots..."
    python "$SCRIPT_DIR/plot_piston_data.py" "." "$STEM" "$OLDSTEPS"
fi

echo "======================================"
echo "Done! Results are in: $RUN_DIR"
echo "======================================"
