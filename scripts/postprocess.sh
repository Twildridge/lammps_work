#!/bin/bash
# postprocess.sh -- the plots of a finished run directory. run_lammps.sh and continue_sim.sh call
# it; run it by hand on any run dir whose job skipped the step (no MD, no MPI):
#
#   bash postprocess.sh <run_dir> <folder> <dataname> <interaction> <totsteps> [oldsteps]
#
# STRAINS / COMPRESSIONS (environment, space-separated) are the sweep levels whose _g<level> /
# _c<level> files the sweep plotters overlay.
set -u
if [ $# -lt 5 ]; then
    sed -n '2,9p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 1
fi
RUN_DIR=$1; FOLDER=$2; DATANAME=$3; INTERACTION=$4; TOTSTEPS=$5; OLDSTEPS=${6:-0}
STEM="${DATANAME}_${INTERACTION}_${TOTSTEPS}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/site.sh"
LIB="$SCRIPT_DIR/lib"

cd "$RUN_DIR" || { echo "Error: run directory not found: $RUN_DIR"; exit 1; }
echo "======================================"
echo "Post-processing $FOLDER  stem $STEM  in $RUN_DIR"
echo "======================================"
load_python_module
python -c "import numpy; print('numpy', numpy.__version__)"

FAILED=0
try() { "$@" || { echo "  WARNING: $(basename "$2") failed (skipping)"; FAILED=1; }; }

try python "$LIB/plot_lammps_log.py" . "$STEM"
case "$FOLDER" in
    solvent_phase|polymer_phase)
        try python "$LIB/plot_eos.py" . "$DATANAME" "$INTERACTION" ;;
    shear_slab_two_pist)
        try python "$LIB/plot_shear_strain_sweep.py" . "$STEM" "${STRAINS:-0.1}"
        try python "$LIB/plot_piston_data.py" . "$STEM" "$OLDSTEPS" ;;
    triaxial_compression_two_pist|compress_slab)
        # every level's _c<level> curves on one figure per quantity, plus the whole-run piston file
        try python "$LIB/plot_compression_strain_sweep.py" . "$STEM" "${COMPRESSIONS:-0.1}" "$OLDSTEPS"
        try python "$LIB/plot_piston_data.py" . "$STEM" "$OLDSTEPS" ;;
    *)  # slab_with_support, chain_length_testing, triaxial_permeation_two_pist, pure phases
        try python "$LIB/plot_stress_profiles.py" . "$STEM" "$OLDSTEPS"
        try python "$LIB/plot_piston_data.py" . "$STEM" "$OLDSTEPS" ;;
esac
echo "Done: $RUN_DIR"
exit $FAILED
