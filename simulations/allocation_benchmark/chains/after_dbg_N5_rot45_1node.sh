#!/usr/bin/env bash
# written by submit_debug_sweep.sh 2026-11-07 20:55:50; run by dbg_N5_rot45_1node when it ends (compute follow-ups only)
cd "/home/dpollard/Documents/lammps_work/simulations/allocation_benchmark"
sbatch --job-name="bench_N5_s1.0_rot45_compute_n4_t128_c1_node_cc5.0_r1" --partition=compute --nodes=4 --ntasks-per-node=128 --cpus-per-task=1 --mem=0 --time=00:20:00 \
    --export=ALL,ROLE=bench,NETWORK="rot45",NBEADS="5",SCALE=1.0,MAP=node,COMM_CUT=5.0,NSTEPS=10000,PREP_STEPS=10000,REP=1,TAG="bench_N5_s1.0_rot45_compute_n4_t128_c1_node_cc5.0_r1" \
    allocation_benchmark.batch
sbatch --job-name="bench_N5_s1.0_rot45_compute_n4_t64_c2_core_cc5.0_r1" --partition=compute --nodes=4 --ntasks-per-node=64 --cpus-per-task=2 --mem=0 --time=00:20:00 \
    --export=ALL,ROLE=bench,NETWORK="rot45",NBEADS="5",SCALE=1.0,MAP=core,COMM_CUT=5.0,NSTEPS=10000,PREP_STEPS=10000,REP=1,TAG="bench_N5_s1.0_rot45_compute_n4_t64_c2_core_cc5.0_r1" \
    allocation_benchmark.batch
sbatch --job-name="bench_N5_s1.0_rot45_compute_n8_t128_c1_node_cc5.0_r1" --partition=compute --nodes=8 --ntasks-per-node=128 --cpus-per-task=1 --mem=0 --time=00:20:00 \
    --export=ALL,ROLE=bench,NETWORK="rot45",NBEADS="5",SCALE=1.0,MAP=node,COMM_CUT=5.0,NSTEPS=10000,PREP_STEPS=10000,REP=1,TAG="bench_N5_s1.0_rot45_compute_n8_t128_c1_node_cc5.0_r1" \
    allocation_benchmark.batch
sbatch --job-name="bench_N5_s1.0_rot45_compute_n1_t128_c1_node_cc5.0_r1" --partition=compute --nodes=1 --ntasks-per-node=128 --cpus-per-task=1 --mem=0 --time=00:20:00 \
    --export=ALL,ROLE=bench,NETWORK="rot45",NBEADS="5",SCALE=1.0,MAP=node,COMM_CUT=5.0,NSTEPS=10000,PREP_STEPS=10000,REP=1,TAG="bench_N5_s1.0_rot45_compute_n1_t128_c1_node_cc5.0_r1" \
    allocation_benchmark.batch
