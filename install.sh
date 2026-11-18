#!/usr/bin/env bash
# One-time setup of this clone on this machine (your computer, Expanse or Pod). After `git clone`:
#
#     bash install.sh
#
# then open a new terminal. On every machine it gives you `lsync` (pull, commit and push in one
# step) and `lsbatch` (submit a .batch file from any folder), and stops notebook plots from being
# committed. On a cluster it also sets up what the jobs need (on Pod: the conda environment for the
# plotters and, if missing, the LAMMPS binary, built on the login node in about 10 min) and ends
# with scripts/cluster_check.sh, whose last line says whether the setup is complete. Safe to run
# again at any time. Options and details: scripts/setup.sh (--no-claude, --no-build).
exec bash "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts/setup.sh" "$@"
