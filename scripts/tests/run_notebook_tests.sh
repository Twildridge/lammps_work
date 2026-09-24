#!/bin/bash
# run_notebook_tests.sh -- execute the two-piston notebooks headlessly (nbconvert)
# against the synthetic flow_data_local tree from make_fixtures.py, with the
# sync cell disabled and the Config pointed at the fixtures (2026-09-16).
# Usage: bash scripts/tests/run_notebook_tests.sh [python] [fixture_root]
set -u
PY=${1:-python3}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$(dirname "$HERE")"
ROOT=${2:-${TMPDIR:-/tmp}/lammps_plot_fixtures}
[ -d "$ROOT/flow_data_local" ] || $PY "$HERE/make_fixtures.py" "$ROOT" || exit 1
OUT="$ROOT/notebooks"; mkdir -p "$OUT"
fail=0
run_nb () {   # name  python-rewrite-expression
    local nb=$1 rewrite=$2
    $PY - "$SCRIPTS/$nb" "$OUT/$nb" "$ROOT" <<PYEOF
import sys, json
src, dst, root = sys.argv[1], sys.argv[2], sys.argv[3]
nb = json.load(open(src))
for c in nb['cells']:
    if c['cell_type'] != 'code':
        continue
    s = ''.join(c['source'])
    s = s.replace('SYNC, FORCE_SYNC = True, False', 'SYNC, FORCE_SYNC = False, False')
    # nbconvert starts the kernel in the copy's directory: make the scripts/ paths absolute
    scripts = str(__import__('pathlib').Path(src).resolve().parent)
    s = s.replace("LIB = Path('lib').resolve()", 'LIB = Path(r"' + scripts + '/lib")')
    s = s.replace('import slab_two_pistons as s2p', 'sys.path.insert(0, r"' + scripts + '"); import slab_two_pistons as s2p')
    s = s.replace("input_file = s2p.DEFAULT_INPUT", 'input_file = r"' + scripts + '/" + s2p.DEFAULT_INPUT')
    s = s.replace("input_file = s2p.SMOKE_INPUT", 'input_file = r"' + scripts + '/" + s2p.SMOKE_INPUT')
    $rewrite
    c['source'] = s
json.dump(nb, open(dst, 'w'))
PYEOF
    ( cd "$SCRIPTS" && $PY -m jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=600 \
        --output "$OUT/executed_$nb" "$OUT/$nb" ) > "$OUT/$nb.log" 2>&1
    if [ $? -eq 0 ]; then echo "  OK   $nb"; else echo "  FAIL $nb  (see $OUT/$nb.log)"; grep -m3 -E "Error|Exception" "$OUT/$nb.log"; fail=1; fi
}
echo "=== notebooks (fixtures under $ROOT)"
run_nb triaxial_compression_single_two_pist.ipynb "s = __import__('re').sub(r'RUN_ID\\s*=\\s*\"[^\"]*\",', 'RUN_ID      = \"fixture_comp\", base_dir=\"' + root + '/flow_data_local\",', s).replace('DATANAME    = \"final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000002_two_pist\"', 'DATANAME    = \"fixture_slab\"').replace('VOR_ENABLE = True,', 'VOR_ENABLE = False,')"
run_nb triaxial_compression_sweep_two_pist.ipynb "s = __import__('re').sub(r'RUN_ID\\s*=\\s*\"[^\"]*\",', 'RUN_ID      = \"fixture_comp\", base_dir=\"' + root + '/flow_data_local\",', s).replace('DATANAME    = \"final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000002_two_pist\"', 'DATANAME    = \"fixture_slab\"').replace('COMP_LEVELS = [\"0.05\", \"0.10\", \"0.15\", \"0.20\"]', 'COMP_LEVELS = [\"0.05\", \"0.10\"]').replace('VOR_ENABLE = True,', 'VOR_ENABLE = False,')"
run_nb triaxial_permeation_single_two_pist.ipynb "s = __import__('re').sub(r'RUN_ID\\s*=\\s*\"[^\"]*\",', 'RUN_ID      = \"fixture_perm\", base_dir=\"' + root + '/flow_data_local\",', s).replace('DATANAME    = \"final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000002_two_pist\"', 'DATANAME    = \"fixture_slab\"')"
# the converter notebook runs the real converter (needs a slab file); skipped when absent
if [ -f "$SCRIPTS/../../lammps_data_files_local/final_config_slab_support_periodic_5beads_tall_rho04_new_1.0_1.0_14000000.data" ]; then
    run_nb slab_two_pistons.ipynb "s = s.replace('log_info           = True,', 'log_info           = False,').replace('output_file        = None,', 'output_file        = \"' + root + '/converter_test.data\",')"
else
    echo "  SKIP slab_two_pistons.ipynb (no slab data file on this machine)"
fi
[ $fail = 0 ] && echo "NOTEBOOK TESTS: OK" || { echo "NOTEBOOK TESTS: FAILED"; exit 1; }
