#!/usr/bin/env bash
# Submit the chain-length sweep: one SLURM job per beads-per-strand value (default 10 15 20 30),
# each running chain_length_testing_sweep.batch with NBEADS exported and the job named
# chain_len_<N>b (so the SLURM output file, the continuation chain and the run directory
# tell the N apart).  Run from this folder on the login node:
#     bash submit_chain_length_sweep.sh            # all four
#     bash submit_chain_length_sweep.sh 30         # one N
# Checks that every input file exists in ~/Documents/lammps_data/input_data/ before
# submitting anything.  2026-11-03.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
NS=("$@"); [ ${#NS[@]} -eq 0 ] && NS=(10 15 20 30)
INPUT_DIR="$HOME/Documents/lammps_data/input_data"
for n in "${NS[@]}"; do
    f="$INPUT_DIR/slab_support_periodic_${n}beads_tall_rho04_random.data"
    [ -f "$f" ] || { echo "missing input file: $f  (scripts/lib/build_beads_sweep.py --beads $n)"; exit 1; }
done
for n in "${NS[@]}"; do
    jid=$(sbatch --parsable --job-name="chain_len_${n}b" --export=ALL,NBEADS="$n" chain_length_testing_sweep.batch)
    echo "N = ${n}: job ${jid} (chain_len_${n}b)"
done
