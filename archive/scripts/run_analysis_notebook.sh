#!/bin/bash
# ==============================================================================
# run_analysis_notebook.sh — headlessly rerun an analysis notebook and report
# back its text output + which plots it (re)wrote, without bloating the
# tracked .ipynb with embedded outputs.
#
# Executes the notebook into a throwaway temp copy (jupyter nbconvert
# --execute), so the tracked notebook file on disk is never touched. Prints
# every cell's text/error output, then lists any .png under flow_data_local/
# whose mtime changed during the run (works regardless of each notebook's own
# PLOT_DIR convention). The temp copy is discarded on exit.
#
# Usage:
#   bash run_analysis_notebook.sh <notebook.ipynb> [timeout_seconds]
#
# Example:
#   bash scripts/run_analysis_notebook.sh triaxial_compression.ipynb
#
# Notes:
#   - Run from anywhere; the notebook path is resolved relative to this
#     script's own directory (scripts/) unless given as an absolute path.
#   - Some notebooks (e.g. triaxial_compression.ipynb) opportunistically
#     rsync missing files from a cluster over SSH; if the cluster/VPN is
#     unreachable those cells print a warning and continue rather than
#     hanging indefinitely, but budget extra wall time for that case.
#   - default per-cell timeout is 1800s (30 min); pass a second argument to
#     override for notebooks with heavier trajectory parsing.
# ==============================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_ROOT="$(cd "$REPO_ROOT/../flow_data_local" 2>/dev/null && pwd || true)"

if [ $# -lt 1 ]; then
    echo "Usage: bash run_analysis_notebook.sh <notebook.ipynb> [timeout_seconds]"
    exit 1
fi

NOTEBOOK="$1"
CELL_TIMEOUT="${2:-1800}"

if [[ "$NOTEBOOK" != /* ]]; then
    NOTEBOOK="$SCRIPT_DIR/$NOTEBOOK"
fi
if [ ! -f "$NOTEBOOK" ]; then
    echo "Notebook not found: $NOTEBOOK" >&2
    exit 1
fi

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT
MARKER="$TMPDIR/marker"
touch "$MARKER"
sleep 1   # ensure mtime resolution reliably marks new/rewritten plots as newer

EXECUTED="$TMPDIR/executed.ipynb"

echo ">>> Executing $(basename "$NOTEBOOK") (timeout ${CELL_TIMEOUT}s/cell)..."
if ! jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout="$CELL_TIMEOUT" \
        --output "$EXECUTED" \
        "$NOTEBOOK" 2> "$TMPDIR/stderr.log"; then
    echo ">>> EXECUTION FAILED"
    echo "--- stderr tail ---"
    tail -n 80 "$TMPDIR/stderr.log"
    exit 1
fi

echo ""
echo "--- text output ---"
python3 - "$EXECUTED" <<'PY'
import json, sys
nb = json.load(open(sys.argv[1]))
for i, c in enumerate(nb.get("cells", [])):
    for out in c.get("outputs", []):
        if out.get("output_type") == "stream":
            sys.stdout.write("".join(out.get("text", [])))
        elif out.get("output_type") == "error":
            print(f"[cell {i} error] {out.get('ename')}: {out.get('evalue')}")
PY

echo ""
echo "--- plots written/updated during this run ---"
if [ -n "$DATA_ROOT" ] && [ -d "$DATA_ROOT/plots" ]; then
    found=$(find "$DATA_ROOT/plots" -name '*.png' -newer "$MARKER")
    if [ -n "$found" ]; then
        echo "$found"
    else
        echo "(none newer than run start — check the notebook actually reached its savefig cells)"
    fi
else
    echo "(no flow_data_local/plots directory found at $DATA_ROOT/plots)"
fi
