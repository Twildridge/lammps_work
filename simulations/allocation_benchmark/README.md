# Allocation benchmark

Finds the nodes, tasks per node and CPUs per task that run the x/y-periodic random slab fastest
at 5 to 30 beads per chain, while checking the results stay correct. Production files are not edited:
each run works on a temporary copy of `simulations/slab_with_support/slab_with_support.lmp`.

## Handoff: if you are Claude running this on Expanse

This work was designed in a Mac session whose memory is not available here. Status on 2026-11-05: everything
below was tested on a Mac, but nothing has run on Expanse yet, and `ppr:<tasks>:node:pe=<cpus>` has never run
under real MPI.

**The goal, from the project mentor:**
1. Find the exact Slurm configuration that is fastest for the standard gel.
2. Find the optimal settings for each variation: N = 10, 15, 20, 25, 30, and the 0.5x and 2x slabs.
3. Compare everything with the baseline, which is the current production setup.

**Definitions:**
- Standard gel: N = 5, scale 1.0, `slab_support_periodic_5beads_tall_rho04_random` (9 x 9 x 23 cells, 238,717 atoms).
- Baseline: `--nodes=4 --ntasks-per-node=128 --cpus-per-task=1`, `--map-by node`, `comm_modify cutoff 5.0`. Its run is
  `bench_N5_s1.0_compute_n4_t128_c1_node_cc5.0_r1`.

**How to work:**
- Stop after the pilot and check it, as described below, before spending the rest.
- Stay under 4,000 SU in total, the limit the user set. Only the user can raise it.
- Do not edit `run_lammps.sh`, `slab_with_support.lmp` or any production batch file without asking.
- If 8 nodes is the fastest for a case, speed may still be rising: try 16 with the `one` tier if the budget allows.
- The `scale` tier runs only 128 x 1. To call a 0.5x or 2x setting optimal, use `one` to try the other tasks x cpus
  layouts at its best node count.
- Repeat each winner with `repeat`, at scale 1.0, and treat differences within the spread of the repeats as ties.

**Report for each case:** the exact `#SBATCH` lines, steps/s, hours and SUs per 1 M steps, and the gain over the
baseline. Keep two answers apart, because `--map-by` lives in `run_lammps.sh` and the cutoff in the deck: the
fastest configuration from Slurm settings alone, and the fastest with those changes too. Mac results in
`local_test_results/` are local tests only and say nothing about MPI layouts.

## Running it on Expanse

Run from this folder on the login node, after `git pull`. Add `--dry-run` to any `submit` line to see its
runs and SU estimate without submitting. Total estimate: about 3,500 SU before repeats, against a 4,000 SU
cap that the script enforces (`SU_CAP`).

```bash
# 1. Build the missing data files. Existing files are kept. The anaconda3 module provides numpy and scipy.
# Run only one copy: two at once on a Mac were still unfinished after 35 min on a build that takes 2 min alone.
# The largest file (N = 5, 2x) took about 10 min on the Mac.
module load anaconda3/2021.05/q4munrg
python build_benchmark_data.py

# 2. Prepare the common starting states (about 380 SU). Wait until squeue is empty.
bash submit_benchmark.sh prep

# 3. Pilot (about 40 SU). Check the mapping and the speeds before spending the rest.
bash submit_benchmark.sh pilot
python3 collect_results.py

# 4. The scan. If the pilot speeds are, say, twice the assumed ones, add RATE_MULT=2.
bash submit_benchmark.sh production     # about 390 SU
bash submit_benchmark.sh standard       # about 450 SU: the standard gel, every layout, unchanged deck
bash submit_benchmark.sh tuned          # about 1,600 SU
bash submit_benchmark.sh factors        # about 240 SU
bash submit_benchmark.sh subnode        # about 20 SU
bash submit_benchmark.sh scale          # about 450 SU

# 5. Results, then repeat the winner of each N to measure the noise.
python3 collect_results.py
bash submit_benchmark.sh repeat 30 2 64 2 core 2.5     # N NODES TASKS CPUS MAP CUT [copies, default 3]
```

Bring `~/Documents/lammps_runs/allocation_benchmark/results.csv` back to the Mac, then:

```bash
python3 plot_results.py results.csv --out allocation_benchmark.png
```

The logs stay in `~/Documents/lammps_runs/allocation_benchmark/logs/`. Run directories on scratch are deleted
after a successful run and kept after a failure.

## What the pilot should confirm

- `bench_N5_..._core_cc2.5_r1` finished. `ppr:<tasks>:node:pe=<cpus>` could not be tested on the Mac, which has
  no MPI, so a mapping failure would first appear here. It fails within seconds, at little cost.
- `python3 collect_results.py` lists no excluded runs.
- The measured speeds are close to the assumed ones. If not, use `RATE_MULT` so the SU check stays honest.

## What varies

| Tier | N | Scale | Nodes | Tasks x CPUs per node | Mapping, comm cutoff |
|------|---|-------|-------|-----------------------|----------------------|
| production | 5, 10, 15, 20, 25, 30 | 1x | 1, 2, 4, 8 | 128 x 1 | node, 5.0 (today's settings) |
| standard | 5 | 1x | 1, 2, 4, 8 | 128 x 1, 64 x 2, 32 x 4, 16 x 8 | node + 5.0, core + 5.0 |
| tuned | 5, 10, 15, 20, 25, 30 | 1x | 1, 2, 4, 8 | 128 x 1, 64 x 2, 32 x 4, 16 x 8 | core, 2.5 |
| factors | 5, 30 | 1x | 2, 4, 8 | 128 x 1 | core + 5.0, then node + 2.5 |
| subnode | 5, 15, 30 | 1x | shared partition | 32 x 1, 64 x 1 | core, 2.5 |
| scale | 5, 15, 30 | 0.5x, 2x | 1, 2, 4, 8 | 128 x 1 | core, 2.5 |

**Scale** is the number of lattice cells, so roughly the number of beads. The 0.5x and 2x files are for timing
only; a network two cells wide is too narrow for physics.

**Starting state.** The `prep` tier takes each fresh file through the deck's full start-up plus 50 k production
steps. Every timed run of that N and scale then restarts from that one file with `cont=1` and runs 20 k steps, so
only the allocation differs. The result is LAMMPS's `Loop time` of that block; start-up, file reading and final
file writing are excluded.

**Output cost.** The 10 k-step cadence of the deck's `fix print` and stress sampling is the same as production.
The `dump` every 100 k steps never falls inside a 20 k-step loop, so its cost is not in these timings. Check the
`Output` share of a production log if dumps become more frequent.

## How correctness is checked

`collect_results.py` excludes a run from the results if any of these fail:

1. **Complete:** LAMMPS reached `Total wall time`.
2. **Safe neighbor lists:** `Dangerous builds = 0`.
3. **Same forces:** PE and pressure at step 0 agree to 1 part in 10^6 with every other run of the same starting
   state. Layouts and cutoffs start from identical coordinates, so a lost interaction would show up here. After
   step 0, Langevin noise makes trajectories differ between layouts even when all of them are correct.

These checks cover a short run. Before adopting a new layout for a long campaign, compare a longer run's mean
pressure, box size and Rg against the existing result of the same N.

## The comm cutoff applies after start-up only

`comm_modify cutoff` has to reach every bonded partner. Fresh files built with `chain_placement='walk'` hold a few
long bonds before the push-off shortens them. Longest bond in each file as built on 2026-11-05:

| N | 5 | 10 | 15 | 20 | 25 | 30 |
|---|---|----|----|----|----|----|
| longest bond (sigma) | 0.92 | 2.39 | 1.86 | 2.66 | 3.72 | 5.95 |

So `prep` uses a cutoff of 7.0, and only the timed runs, which start after the push-off, use 2.5. If 2.5 is adopted
in production, set it just before the production loop, not at the top of the deck. N = 30 has 4 bonds longer than
the deck's current 5.0, so a fresh N = 30 start at 5.0 can stop with `Bond atoms missing`, depending on where the
subdomain boundaries fall.
