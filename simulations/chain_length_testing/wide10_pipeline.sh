#!/usr/bin/env bash
# wide10_pipeline.sh (2026-11-14) -- ONE-SHOT submission of the wide random 10-bead gel, from the
# data file to the first measurements. Run once on the login node and log out:
#     bash simulations/chain_length_testing/wide10_pipeline.sh
#   1. wide10_build.batch            build the slab (debug, plus a shared fallback copy)
#   2. chain_length_testing_wide10.batch   equilibrate (Rg halt, ends at the mean box), after 1
#   3. wide10_bridge.batch           after 2: void check + gate, convert, submit isotropy,
#                                    compression, periodic shear, frozen-reference permeation
# Status at any time:  bash simulations/chain_length_testing/wide10_status.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
PIPE=$HOME/Documents/lammps_runs/chain_length_testing/wide10_pipeline_20261114
mkdir -p "$PIPE"
B1=$(sbatch --parsable wide10_build.batch)
B2=$(sbatch --parsable --partition=shared --time=02:00:00 --job-name=wide10_build_sh wide10_build.batch)
EQ=$(sbatch --parsable --dependency="afterok:${B1}?afterok:${B2}" chain_length_testing_wide10.batch)
BR=$(sbatch --parsable --dependency=afterany:"$EQ" --export=ALL,SLAB_JOB="$EQ" wide10_bridge.batch)
{ echo "$B1|build (debug)"; echo "$B2|build (shared fallback; exits at once if the file exists)"
  echo "$EQ|equilibration, 6 nodes, 16 M-step budget, Rg halt, mean-box ending"; echo "$BR|bridge: void gate, convert, submit measurements"; } >> "$PIPE/jobs.txt"
echo "$(date +%Y-%m-%dT%H:%M) pipeline submitted: build $B1 / $B2, equilibration $EQ, bridge $BR" >> "$PIPE/status.log"
echo "build $B1 (debug) / $B2 (shared), equilibration $EQ, bridge $BR"
