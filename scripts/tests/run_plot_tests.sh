#!/bin/bash
# run_plot_tests.sh -- regenerate the synthetic fixtures and run the three
# plotting scripts on the one-piston AND two-piston formats (2026-09-16).
# Usage: bash scripts/tests/run_plot_tests.sh [python] [fixture_root]
set -u
PY=${1:-python3}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$(dirname "$HERE")"
ROOT=${2:-${TMPDIR:-/tmp}/lammps_plot_fixtures}
rm -rf "$ROOT"; mkdir -p "$ROOT"
$PY "$HERE/make_fixtures.py" "$ROOT" || exit 1
fail=0
for case in one_piston two_piston_perm two_piston_comp; do
    dir="$ROOT/$case"
    stem=$(ls "$dir"/final_*.data | head -1 | sed -E 's|.*/final_(triperm\|tricomp)_||; s|\.data$||')
    echo "=== $case  (stem $stem)"
    ( cd "$dir" && $PY "$SCRIPTS/plot_piston_data.py" . "$stem" 0 && $PY "$SCRIPTS/plot_stress_profiles.py" . "$stem" 0 \
        && $PY "$SCRIPTS/plot_lammps_log.py" . "$stem" ) 2>&1 | grep -v -E "UserWarning|from pandas|findfont|tight_layout|set_ylim" | sed 's/^/    /'
    if [ "$case" = two_piston_comp ]; then
        # the sweep plotter + per-level piston plots, as postprocess.sh does
        ( cd "$dir" && $PY "$SCRIPTS/plot_compression_strain_sweep.py" . "$stem" "0.05 0.10" 0 \
            && for l in 0.05 0.10; do $PY "$SCRIPTS/plot_piston_data.py" . "${stem}_c$l" 0; done ) 2>&1 | grep -v -E "UserWarning|from pandas|findfont|tight_layout|set_ylim" | sed 's/^/    /'
    fi
    if [ "$case" = two_piston_comp ]; then
        want=("$dir"/output_plots/*_piston.png "$dir"/output_plots/*_piston_sweep.png "$dir"/output_plots/*_stress_profiles_sweep.png "$dir"/output_plots/convergence_plots/*_convergence.png)
    else
        want=("$dir"/output_plots/*_piston.png "$dir"/output_plots/*_stress.png "$dir"/output_plots/convergence_plots/*_convergence.png "$dir"/output_plots/convergence_plots/*_flow_diagnostics.png)
    fi
    for png in "${want[@]}"; do
        [ -s "$png" ] && echo "    OK  $(basename "$png")" || { echo "    MISSING $(basename "$png")"; fail=1; }
    done
done
echo "fixtures + plots under $ROOT"
[ $fail = 0 ] && echo "PLOT TESTS: OK" || { echo "PLOT TESTS: FAILED"; exit 1; }
