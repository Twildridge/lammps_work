#!/usr/bin/env bash
# One-time setup of this clone on this machine (your computer or a cluster). After `git clone`:
#
#     bash install.sh
#
# then open a new terminal. It gives you the commands `lsync` (pull, commit and push in one
# step) and `lsbatch` (submit a .batch file from any folder), and stops notebook plots from
# being committed. Safe to run again at any time. Details: scripts/install_lsync.sh.
exec bash "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts/install_lsync.sh" "$@"
