#!/usr/bin/env bash
# pod_setup.sh -- the Pod-only part of the one-time setup. install.sh runs it when lib/site.sh
# detects Pod; it can also be run on its own at any time:
#     bash scripts/pod_setup.sh
# Idempotent: it creates what is missing and reports what is there.
#   1. module load miniconda and create the conda environment lammps_analysis (python 3.11, numpy,
#      scipy, matplotlib, pandas) if it does not exist; verify the four imports either way, adding a
#      package that is missing. The jobs activate this environment for the plotters
#      (load_python_module in lib/site.sh), so an incomplete one breaks post-processing, not just notebooks.
#   2. Report whether the LAMMPS binary of the Pod block of lib/site.sh (LMP_BIN) exists. Pod has no
#      site LAMMPS; every account builds its own once with lammps_build/build_lammps_pod.sh.
set -o pipefail      # not -u: Lmod's init file, the module function and conda.sh read unset variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
export CLUSTER=pod
source "$SCRIPT_DIR/lib/site.sh" || exit 1

# `module` is a shell function in a login shell; a clean shell (cron, ssh host cmd) needs Lmod's init file
type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh 2>/dev/null || source /etc/profile.d/modules.sh 2>/dev/null \
    || { echo "!! no 'module' command: is this a Pod login node?"; exit 1; }
module load "$PYTHON_MODULE" || { echo "!! module load $PYTHON_MODULE failed"; exit 1; }
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

ENV_NAME="$PYTHON_CONDA_ENV"
PKGS=(python=3.11 numpy scipy matplotlib pandas)
if conda env list 2>/dev/null | awk 'NF && $1 !~ /^#/ {print $1}' | grep -qx "$ENV_NAME"; then
    echo ">>> conda env $ENV_NAME exists"
else
    echo ">>> creating conda env $ENV_NAME (${PKGS[*]}), a few minutes..."
    conda create -n "$ENV_NAME" "${PKGS[@]}" -y || { echo "!! conda create failed"; exit 1; }
fi
if conda run -n "$ENV_NAME" python -c "import numpy, scipy, matplotlib, pandas" >/dev/null 2>&1; then
    echo ">>> $ENV_NAME: numpy, scipy, matplotlib, pandas import"
else
    echo ">>> $ENV_NAME is missing a package; installing ${PKGS[*]:1}"
    conda install -n "$ENV_NAME" "${PKGS[@]:1}" -y || { echo "!! conda install failed"; exit 1; }
    conda run -n "$ENV_NAME" python -c "import numpy, scipy, matplotlib, pandas" \
        || { echo "!! $ENV_NAME still cannot import numpy, scipy, matplotlib, pandas"; exit 1; }
    echo ">>> $ENV_NAME: numpy, scipy, matplotlib, pandas import"
fi

if [ -x "$LMP_BIN" ]; then
    echo ">>> LAMMPS binary: $LMP_BIN"
else
    echo "!! no LAMMPS binary at $LMP_BIN"
    echo "   Build it once, on this login node (20-40 min, every thread):"
    echo "       bash $REPO_ROOT/lammps_build/build_lammps_pod.sh"
    echo "   then check the whole setup with:  bash $REPO_ROOT/scripts/cluster_check.sh"
    exit 2
fi
