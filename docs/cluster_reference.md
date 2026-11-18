# Cluster Reference

[← back to README](../README.md)

Production runs on Expanse (SDSC). Pod (CNSI, UCSB) is the backup since 2026-11-18. Both use the same runner (`scripts/run_lammps.sh`), decks and batch files. `scripts/lib/site.sh` holds one block of settings per cluster and picks the block from the Slurm cluster name and the host name (`CLUSTER=pod` or `CLUSTER=expanse` in the environment forces it). `bash scripts/cluster_check.sh` on a login node compares that block with the machine and changes nothing.

## Settings side by side

| Setting | Expanse | Pod |
|---------|---------|-----|
| Login node | `login.expanse.sdsc.edu` | `pod-login1.cnsi.ucsb.edu` (campus VPN) |
| Cores/node (`CORES_PER_NODE`) | 128 (2 x AMD EPYC 7742) | 40 (2 x Intel Xeon Gold 6148), 64 nodes, Omni-Path |
| Partition, account | `--partition=compute --account=csb197 --mem=0` | `batch` (node[1-55,57-63]), no account; `short` = nodes 56 and 64, 2 h limit; `largemem`, `gpu` (V100) unused |
| Modules (`LAMMPS_MODULES`) | `gcc/10.2.0 openmpi/4.1.3 python/3.8.12` | `gnu14 openmpi5` |
| LAMMPS binary (`LMP_BIN`) | `/home/dpollard/software/lammps/22Jul2025_update3/.../build/lmp` | `~/MD/lammps/lammps-22Jul2025/build/lmp` |
| MPI fabric (`mpi_env`) | UCX over InfiniBand (`UCX_TLS`, `UCX_NET_DEVICES=mlx5_2:1`) | PSM2 over Omni-Path (`OMPI_MCA_pml=cm`, `OMPI_MCA_mtl=psm2`) |
| Rank placement (`MPI_MAP_BY`) | `core` = `ppr:128:node:pe=1`, measured 1.4-1.6x over `node` | `core` = `ppr:40:node:pe=1`, not measured; the archived runner used `node` |
| Comm cutoff | 2.5 (deck default since 2026-11-07) | the same decks, so 2.5 |
| Scratch (`SCRATCH_DIR`) | `/expanse/lustre/scratch/$USER/temp_project` | `/scratch/$USER` (19 TB SSD; `/bigscratch` is 70 TB NVMe; neither backed up) |
| Python for the plotters | `module load anaconda3/2021.05/q4munrg` | `module load miniconda` + `conda activate lammps_analysis` |
| Default production allocation | 6 nodes x 128 x 1 | 4 nodes x 40 x 1, provisional |
| Charge | SUs (core-hours) | none |

The Pod column was copied from the runner archived on 2026-10-11 and from the CNSI hardware page. It had not been re-run on Pod when this was written.

## Submitting a job

`sbatch` reads the `#SBATCH` lines before the script runs, so a batch file cannot choose its own allocation. Each ported batch file keeps its Expanse allocation as `#SBATCH` lines and has one more line for Pod:

```
#POD --nodes=4 --ntasks-per-node=40 --cpus-per-task=1 --time=112:00:00
```

The command `lsbatch` (`scripts/submit.sh`, installed by `bash install.sh`) reads the line of the cluster it is on. It works from any folder: the batch file is found by name under `simulations/`, and the job is submitted from that file's folder, so the SLURM output file lands there as it does with `sbatch`.

```bash
lsbatch uniaxial_compression_single.batch              # Pod: the #POD line; Expanse: plain sbatch
DATANAME=<stem> lsbatch --nodes=2 uniaxial_compression_single.batch
lsbatch --dry-run uniaxial_compression_single          # print the allocation, submit nothing (.batch is optional)
```

On Expanse `sbatch <file>.batch` from the batch file's folder also works. On Pod a plain `sbatch` sends the Expanse lines and is rejected.

Ported to Pod: the `_single` and `_sweep` batch of every simulation folder, `uniaxial_compression_random` and `slab_with_support_periodic_random`. Not ported: the pipelines (`submit_*_pipeline.sh`, the `chain_*` bridges), `chain_length_testing/` and `calibration_sweep/`, which name Expanse partitions and paths.

## Expanse

| | |
|---|---|
| Data transfer node | `data.expanse.sdsc.edu` |
| Partitions | `compute` (48 h limit); `debug` (2 nodes, 30 min, 2 jobs per user) for smoke tests and layout scans; `shared` with one core for bridge jobs that only copy and convert files |
| Account | `--account=csb197`, required |
| Build | `lammps_build/` (22Jul2025 update 3, MPI + OpenMP) |

### Allocation size (benchmark of 2026-11-07)

`simulations/allocation_benchmark/` timed the rot45 slab over every ranks × threads split, both MPI placements, two comm cutoffs and 1 to 16 nodes. Three settings matter, in order of gain. The runner and the decks apply all three; nothing is set by hand.

1. One MPI rank per core, 128 × 1 per node. Every hybrid split loses: 64 × 2 by 7-30 %, 32 × 4 and 16 × 8 by 2-3x.
2. Ranks fill a node before the next one (`--map-by ppr:128:node:pe=1`, `MPI_MAP_BY=core` in `site.sh`). Round-robin placement (`--map-by node`) costs 1.4-1.6x at 2 nodes and made 2 to 8 nodes slower than 1.
3. `comm_modify cutoff 2.5` in place of 5.0 once the bonds are relaxed (`variable comm_cutoff index 2.5` in every deck): another 1.4x with identical step-0 forces. `-var comm_cutoff 5.0` restores the old value.

With all three in place (steps/s / SU per 1 M steps):

| Beads per strand | 1 node | 4 nodes | 6 nodes | 8 nodes | 16 nodes |
|---|---|---|---|---|---|
| 5 (354 k atoms) | 362 / 98 | 710 / 200 | **893 / 239** | 753 / 378 | 849 / 670 |
| 15 (660 k atoms) | 189 / 188 | 368 / 387 | 460 / 463 | **501 / 567** | 537 / 1059 |

The setup before that (4 nodes, map-by node, cutoff 5.0) ran the 5-bead gel at 236 steps/s for 604 SU per 1 M steps. Repeat spread at 3 to 6 nodes is 20-30 %, because node placement differs between jobs, so neighbouring node counts are ties.

- 6 nodes is the production default (`--nodes=6 --ntasks-per-node=128 --cpus-per-task=1 --mem=0`).
- 1 node is the cheapest per step (2.4x fewer SUs than 6 nodes for the 5-bead gel) and suits anything not on the critical path. 4 nodes is the compromise: `lsbatch --nodes=4 <file>.batch`.
- Bigger gels scale further. The 15-bead gel gains up to 16 nodes at a steep SU price; 8 to 10 nodes is the sensible ceiling (`--nodes=10` for the 15-bead two-piston runs, 831 k atoms).
- SUs are charged per allocated core-hour: nodes × 128 × hours. `sbatch` reserves the full `--time` from the allocation at submit time, so keep `--time` close to the expected length when several jobs are queued. The `Performance:` lines in `log.lammps` give the rate to check against.

The raw tables are in `~/Documents/lammps_runs/allocation_benchmark/{nodescan,layouts}_rot45.txt` on Expanse; the method is in the benchmark's `README.md`.

## Pod

Pod requires the UCSB campus VPN. Install Ivanti Secure Access from [it.ucsb.edu/ivanti-secure-access-campus-vpn/get-connected-campus-vpn](https://it.ucsb.edu/ivanti-secure-access-campus-vpn/get-connected-campus-vpn) and connect before `ssh pod-login1.cnsi.ucsb.edu`.

The plotters need a conda environment that is created once:

```bash
module load miniconda
conda create -n lammps_analysis python=3.11 numpy scipy matplotlib pandas -y
```

The jobs activate it themselves (`load_python_module` in `site.sh`).

Open points until the benchmark has run on Pod (`simulations/allocation_benchmark/README.md`, "Running it on Pod"):

- The `#POD` lines are placeholders: 4 nodes (1 for the pure-phase runs) and four times the Expanse time.
- `--map-by ppr:40:node:pe=1` has not run under Pod's Open MPI 5. The archived runner used `--map-by node` with cutoff 5.0, and its notes show 2 nodes no faster than 1.
- The time limit of the `batch` partition is not known. The two-piston decks cannot be continued, so a level has to fit inside one job.

The runner and batches of 2026-10 are kept in `archive/pod_bridges/` with the benchmark notes [`documenting_pod_runs.md`](../archive/pod_bridges/documenting_pod_runs.md).

## Bridges-2 (PSC), not in use since 2026-10-11

Its runner and batches are in `archive/pod_bridges/`; `archive/README.md` says how to bring them back.

| Setting | Value |
|---------|-------|
| Login / data transfer node | `bridges2.psc.edu` / `data.bridges2.psc.edu` |
| CPU partition | `RM`, 128 cores/node (120 ranks ran faster than 128), 5 day limit |
| Scratch | `/ocean/projects/chm250028p/dpollard/` |
| Home quota | 10 GB, so trajectories must go to scratch (`traj_files/` is a symlink there) |
| LAMMPS binary | `/opt/packages/LAMMPS/lammps-22Jul2025/build-RM-gcc13.3.1/lmp` |
| Modules | `python/3.8.6 openmpi/5.0.8-gcc13.3.1 cuda/12.6.1 intel-mkl/2023.2.0 LAMMPS/22Jul25-gcc` after `module purge` |
