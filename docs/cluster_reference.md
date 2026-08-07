# Cluster Reference

[← back to README](../README.md)

### Expanse (SDSC) — `run_lammps.sh`

| Setting | Value |
|---------|-------|
| Login node | `login.expanse.sdsc.edu` |
| Data transfer node | `data.expanse.sdsc.edu` |
| CPU partition | `compute` |
| GPU partition | `gpu` |
| Account flag | Required: `--account=csb197` |
| Cores/node | 128 |
| Max walltime (compute) | 48 hours |
| Scratch | `/expanse/lustre/scratch/$USER/temp_project` |
| LAMMPS binary | loaded via `module load LAMMPS/...` then `lmp` |

Load modules (in this order):
```bash
module reset
module load gcc/10.2.0
module load openmpi/4.1.3
module load python/3.8.12
```

---

### Bridges-2 (PSC) — `run_lammps_bridges.sh`

| Setting | Value |
|---------|-------|
| Login node | `bridges2.psc.edu` |
| Data transfer node | `data.bridges2.psc.edu` |
| CPU partition | `RM` |
| GPU partition | `GPU` |
| Cores/node | 128 (use 120 for optimal MPI performance) |
| Max walltime | 5 days |
| Scratch | `/ocean/projects/chm250028p/dpollard/` |
| Home quota | 10 GB — trajectory files must go to scratch |
| LAMMPS binary | `/opt/packages/LAMMPS/lammps-22Jul2025/build-RM-gcc13.3.1/lmp` |

Load modules (in this order):
```bash
module purge
module load python/3.8.6
module load openmpi/5.0.8-gcc13.3.1
module load cuda/12.6.1
module load intel-mkl/2023.2.0
module load LAMMPS/22Jul25-gcc
```

> **Bridges-2 quota:** If you see "Disk quota exceeded", trajectory files have filled your 10 GB home directory. Check that `traj_files/` is a symlink pointing to scratch: `ls -la ~/Documents/lammps_runs/*/traj_files`

---

### Pod (CNSI, UCSB) — `run_lammps_pod.sh`

| Setting | Value |
|---------|-------|
| Login node | `pod-login1.cnsi.ucsb.edu` |
| Cores/node | 40 (optimal) |
| Internode bandwidth | 100 Gb/s (half of Bridges-2's 200 Gb/s — avoid multi-node for large systems) |
| Note | GPUs available (L40S); check [`documenting_pod_runs.md`](../documenting_pod_runs.md) for benchmarks |

> **Pod requires the UCSB campus VPN.** You must be connected before you can SSH in. Download and install Ivanti Secure Access from [it.ucsb.edu/ivanti-secure-access-campus-vpn/get-connected-campus-vpn](https://it.ucsb.edu/ivanti-secure-access-campus-vpn/get-connected-campus-vpn), then connect to the UCSB VPN before running `ssh pod-login1.cnsi.ucsb.edu`. Expanse and Bridges-2 do not require a VPN.

---
