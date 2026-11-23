# Allocation benchmark

Finds the nodes, tasks per node and CPUs per task that run the x/y-periodic slab fastest
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

## Running it on Pod (set up 2026-11-18, not yet run)

Pod is the backup for when the Expanse allocation runs out. The question is the same as on Expanse: which
nodes x ranks x threads runs the rot45 standard gel (5 and 15 beads per strand) fastest, with `--map-by core`
and comm cutoff 2.5, the two settings that gave 1.4-1.6x and 1.4x on Expanse. The archived Pod runner used
`--map-by node` and cutoff 5.0, and its notes (`archive/pod_bridges/documenting_pod_runs.md`) show 2 nodes no
faster than 1 with those settings.

`submit_benchmark.sh` reads the cluster from `scripts/lib/site.sh`: on Pod a node is 40 cores, the layouts are
40x1, 20x2, 10x4 and 8x5, jobs go to the `batch` partition (Pod's default) with no account, and time limits are four times the
Expanse ones. "SU" in the output is core-hours; nothing is charged on Pod. `submit_debug_sweep.sh` and the
`subnode` tier are Expanse-only.

```bash
# 0. On the Mac, VPN on: the two fresh files (354k and 660k atoms), after pushing the repo
scp lammps_data_files_local/slab_support_periodic_{5,15}beads_tall_rho04_rot45.data \
    pod-login1.cnsi.ucsb.edu:Documents/lammps_data/input_data/

# 1. On Pod: does lib/site.sh match the machine? Fix its pod block until this is clean.
cd ~/Documents/lammps_work && git pull && bash scripts/cluster_check.sh

# 2. The common starting states (2 jobs, 1 node each), then wait for them.  Each tier below takes
#    NETWORK and NS from the environment; the node health screen (section further down) decides
#    the --exclude list on its own through lsbatch's batch_header, so nothing is pinned by hand.
cd simulations/allocation_benchmark
export NETWORK=rot45 NS="5 15"
bash submit_benchmark.sh prep --dry-run && bash submit_benchmark.sh prep

# 3. Pilot (4 jobs): the first run of ppr:40:node:pe=1 under Open MPI 5. Check it before the rest.
bash submit_benchmark.sh pilot
python3 collect_results.py                      # no excluded runs; note the steps/s against the assumed ones

# 4. Node count at 40 x 1 (54 jobs) and ranks x threads on 1, 2 and 4 nodes (30 jobs)
bash submit_benchmark.sh nodescan               # SCAN_NODES="1 2 3 4" if 6 and 8 nodes do not start
bash submit_benchmark.sh layouts

# 5. Results
python3 collect_results.py
R=~/Documents/lammps_runs/allocation_benchmark
python plot_nodescan.py $R/results.csv --network rot45 --out $R/nodescan_rot45_pod.png > $R/nodescan_rot45_pod.txt
python plot_layouts.py  $R/results.csv --network rot45 --out $R/layouts_rot45_pod.png  > $R/layouts_rot45_pod.txt
```

The plots need matplotlib (`module load miniconda; conda activate lammps_analysis`). What to check in the pilot:
`logs/<tag>.bindings` shows 40 ranks on each node for `core`; a mapping error appears within seconds, and
`MAP=node` runs are the fallback. Threads above 1 are safe for timing, but `bond_style hybrid` with `-sf omp`
and more than one thread under-tallies the per-atom bond virial (found 2026-09-22), so a hybrid layout that
wins here needs the stress output checked before it is used for production.

What to do with the result: set the node count and `--time` of the `#POD` line in each batch file (they are
4 nodes and four times the Expanse time until then), set `MPI_MAP_BY` in the pod block of `lib/site.sh` if
`node` wins, and add a Pod table next to the Expanse one in `docs/running_simulations.md` §5f.

### The random 10-bead gels on Pod (2026-11-23)

The production gel since 2026-11 is the random 10-bead network, in two forms: the periodic wide slab
`slab_support_periodic_10beads_wide_rho04_random.data` (the two-piston runs; about 750 k atoms) and its
isolated cube twin `slab_support_10beads_cube10_rho044_nowalls_random_br2_s11.data` (the cube pipeline,
`scripts/submit_cube_pipeline.sh`). Neither stem fits the `periodic_<N>beads_tall_rho04_<network>`
family the tiers assume, so `BENCH_STEM=<stem>` names the fresh file directly and `NETWORK` becomes a
free label that keeps the two gels apart in the tags, the prepared states (`prep_10b_s1.0_<label>.data`)
and `results.csv`; `NS` must name the single N (10). Both files are built on the laptop, so copy them
first:

```bash
# on the laptop
scp ~/Documents/lammps_data/input_data/slab_support_periodic_10beads_wide_rho04_random.data \
    ~/Documents/lammps_data/input_data/slab_support_10beads_cube10_rho044_nowalls_random_br2_s11.data \
    pod-login1.cnsi.ucsb.edu:Documents/lammps_data/input_data/

# on Pod, one gel at a time; the same tier order as above
cd ~/Documents/lammps_work/simulations/allocation_benchmark
export NS=10 NETWORK=random_wide BENCH_STEM=slab_support_periodic_10beads_wide_rho04_random
bash submit_benchmark.sh prep --dry-run && bash submit_benchmark.sh prep
export NETWORK=random_cube BENCH_STEM=slab_support_10beads_cube10_rho044_nowalls_random_br2_s11
bash submit_benchmark.sh prep --dry-run && bash submit_benchmark.sh prep
# ... pilot, nodescan, layouts per gel; plots with --network random_wide / --network random_cube
```

Dry runs on 2026-11-23 (per gel): prep 1 job, pilot 4, nodescan 27, layouts 15, about 250 core-hours.

## Debug-partition layout scan of the rot45 standard gel (2026-11-07)

`submit_debug_sweep.sh` answers the MPI x OpenMP question for one slab family quickly and almost
for free, on the debug partition (2 nodes, 30 min, 2 jobs per user): per N it submits one 1-node
job that preps the starting state (deck start-up + 10 k steps) and then times 128x1, 64x2, 32x4
and 16x8 ranks x threads per node, each as its own `mpirun` inside the same allocation
(`debug_sweep.batch`, layout via `BENCH_TASKS` / `BENCH_CPUS` of `run_benchmark.sh`), plus 128x1
at comm cutoff 2.5. The 2-node job (the same layouts with `--map-by node` and `--map-by core`) is
submitted by `chains/debug_waiter.sh`, a `nohup` process the submit script leaves on the login
node, as soon as a debug slot is free (a job cannot submit its successor: the 2-job submit limit
counts the running job too). When the 1-node job ends it submits, unless `COMPUTE=0`, the
compute-partition follow-ups for the node counts debug cannot give: 4 nodes 128x1 node (the
production baseline), 4 nodes 64x2 core, 8 nodes 128x1, and 1 node 128x1 as the cross-check that
debug and compute nodes (both zen2, 128 cores, 256 GB; debug is `exp-9-[55-56]` beside the compute
rack 9) time alike.
Every job refreshes `~/Documents/lammps_runs/allocation_benchmark/results.csv`,
`layouts_rot45.png` and `layouts_rot45.txt` (ranked table, gain over the baseline, hours and SUs
of the full production budget, the `#SBATCH` lines to copy) when it finishes.

```bash
bash submit_debug_sweep.sh --dry-run
bash submit_debug_sweep.sh                       # NS="5 15" NETWORK=rot45 COMPUTE=1 by default
python3 collect_results.py                       # later, on the login node: csv + verdict
module load anaconda3/2021.05/q4munrg && python plot_layouts.py ~/Documents/lammps_runs/allocation_benchmark/results.csv
```

Layouts a job could not fit before its 30 min are listed in `logs/<jobname>.skipped`; re-running
`submit_debug_sweep.sh` skips finished runs and times only the missing ones. `NETWORK=random`
runs the same scan on the random slab (`prep_<N>b_s1.0.data` states are shared with the tiers
below). The 15-bead rot45 file in `input_data/` on 2026-11-07 was a byte-identical copy of the
5-bead one; the scan uses the file rebuilt by `scripts/lib/build_beads_sweep.py --network rot45
--beads 15` in `input_data/allocation_benchmark/rot45_rebuilt/` when that exists.

## Node-count scan on the compute partition (2026-11-07, evening)

The debug scan settled the ranks x threads question (128 x 1 wins everywhere) and the two
settings (`--map-by core`, comm cutoff 2.5) but could only reach 1 and 2 nodes. `nodescan` fills
in the node axis on the compute partition, with the same 10 k-step blocks as the debug jobs:

```bash
NETWORK=rot45 bash submit_benchmark.sh nodescan --dry-run   # 56 runs, about 800 SU at the assumed speeds
NETWORK=rot45 bash submit_benchmark.sh nodescan             # submitted 2026-11-07 23:35, jobs 55604446-55604501
```

Per N (5 and 15): 128 x 1 core cc2.5 on 1/2/3/4/6/8 nodes three times each and once on 16;
128 x 1 core cc5.0 (the Slurm-only answer, deck unchanged) once on 1/3/4/6/8; 128 x 1 node
cc2.5 and 64 x 2 core cc2.5 once on 2 and 4 nodes. `submit_benchmark.sh` now honours
`NETWORK` and `NS` for every tier. `allocation_benchmark.batch` does not refresh the figure, so
after the last job run on the login node:

```bash
python3 collect_results.py
module load cpu/0.17.3b anaconda3/2021.05/q4munrg
python plot_layouts.py ~/Documents/lammps_runs/allocation_benchmark/results.csv --network rot45 \
    --out ~/Documents/lammps_runs/allocation_benchmark/layouts_rot45.png > ~/Documents/lammps_runs/allocation_benchmark/layouts_rot45.txtpython plot_nodescan.py ~/Documents/lammps_runs/allocation_benchmark/results.csv --network rot45 \
    --out ~/Documents/lammps_runs/allocation_benchmark/nodescan_rot45.png > ~/Documents/lammps_runs/allocation_benchmark/nodescan_rot45.txt
```

`plot_nodescan.py` is the figure for the node-count decision: steps/s and SU per 1 M steps against
nodes at 128 x 1, one line per setting (production today, map core only, map core + cutoff 2.5),
the bar through each point spanning the slowest and fastest repeat.  Completed 2026-11-07 23:42,
56 runs, all valid, about 420 SU.  Repeat spread on the compute partition reached 20-30 % at 3-6
nodes (node placement differs between jobs), so neighbouring node counts inside the bars are ties.

**Deck change made the same evening.** Every active deck under `simulations/` now sets
`variable comm_cutoff index 2.5` and issues `comm_modify cutoff ${comm_cutoff}` once the bonds
are relaxed: `slab_with_support.lmp` (and its symlink `chain_length_testing.lmp`) and
`pure_polymer.lmp` keep 5.0 through the push-off of a fresh file and switch right after it; the
two-piston, phase and solvent decks, which read relaxed files, switch at start-up. `-var
comm_cutoff 5.0` restores the old behaviour. **Decision 2026-11-08:** `MPI_MAP_BY=core` in `scripts/lib/site.sh` (the runner maps `ppr:128:node:pe=1`) and 6 nodes x 128 x 1 in
every production batch; the summary table is in `docs/running_simulations.md` §5f.

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
