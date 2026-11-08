#!/usr/bin/env bash
# Submit the chain-length sweep: one SLURM job per beads-per-strand value (default 10 15 20 30
# for the random network, 5 15 for rot45), each running chain_length_testing_sweep.batch with
# NBEADS and NETWORK exported and the job named chain_len_<N>b[_<network>] (so the SLURM output
# file, the continuation chain and the run directory tell the runs apart).  Run from this folder
# on the login node:
#     bash submit_chain_length_sweep.sh                    # the random sweep, all four
#     bash submit_chain_length_sweep.sh 30                 # one N
#     NETWORK=rot45 bash submit_chain_length_sweep.sh      # the rot45 lattice, 5 and 15 beads (2026-11-07)
#     NETWORK=rot45 bash submit_chain_length_sweep.sh 15   # one of them
# Checks that every input file exists in ~/Documents/lammps_data/input_data/ before
# submitting anything.  2026-11-03.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
NETWORK="${NETWORK:-random}"
NS=("$@")
if [ ${#NS[@]} -eq 0 ]; then
    if [ "$NETWORK" = rot45 ]; then NS=(5 15); else NS=(10 15 20 30); fi
fi
INPUT_DIR="$HOME/Documents/lammps_data/input_data"
for n in "${NS[@]}"; do
    f="$INPUT_DIR/slab_support_periodic_${n}beads_tall_rho04_${NETWORK}.data"
    [ -f "$f" ] || { echo "missing input file: $f  (scripts/lib/build_beads_sweep.py --network $NETWORK --beads $n)"; exit 1; }
done
SUFFIX=""; [ "$NETWORK" = random ] || SUFFIX="_${NETWORK}"
for n in "${NS[@]}"; do
    jid=$(sbatch --parsable --job-name="chain_len_${n}b${SUFFIX}" --export=ALL,NBEADS="$n",NETWORK="$NETWORK" chain_length_testing_sweep.batch)
    echo "N = ${n} (${NETWORK}): job ${jid} (chain_len_${n}b${SUFFIX})"
done
