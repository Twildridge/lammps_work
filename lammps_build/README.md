# Building LAMMPS on Expanse (the binary `run_lammps.sh` uses)

`run_lammps.sh` hard-codes the binary this build produces:

```
~/software/lammps/22Jul2025_update3/mpi-omp/gcc/10.2.0/openmpi/4.1.3/lammps-22Jul2025/build/lmp
```

To reproduce it on SDSC Expanse (one-time, ~1 h on the `shared` partition):

```bash
# from a login node, in this folder or anywhere
sbatch build_lammps.sh          # edit --account first if you are not on csb197
```

What the job does (`build_lammps.sh`):

1. loads `gcc/10.2.0`, `openmpi/4.1.3`, `cmake/3.21.4`, `fftw/3.3.10`, `hdf5/1.10.7` — the same compiler + MPI modules the run scripts load, so the binary and the job environment match;
2. downloads `lammps-22Jul2025_update3.tar.gz` from download.lammps.org into `~/software/lammps/...`;
3. configures with CMake: MPI + OpenMP (`-sf omp` in the run scripts), KOKKOS (OpenMP backend, Zen2), FFTW3, and the packages the decks need — MOLECULE, EXTRA-FIX, EXTRA-PAIR, EXTRA-MOLECULE, RIGID, MISC, MC, KSPACE, VORONOI — plus a broad set of others (see the `-DPKG_*` list);
4. builds with `make -j 32`.

`build-lammps.o46984163.exp-1-06.log` is the SLURM output of the build that produced the binary in use (2026-03-04), kept as a record of the module versions and CMake summary.

## Pod (CNSI, UCSB)

`build_lammps_pod.sh` is the same build for Pod (binary `~/software/lammps/22Jul2025_update3/mpi-omp/gnu14/14.2.0/openmpi5/5.0.7/lammps-22Jul2025/build/lmp`, the `LMP_BIN` of the Pod block in `scripts/lib/site.sh`). The CMake options and the package list are identical; what differs is listed in its header: modules `gnu14/14.2.0 openmpi5/5.0.7 hdf5/1.14.6`, `Kokkos_ARCH_SKX` for Pod's Skylake Xeons, `FFT=KISS` because Pod has no FFTW (the decks use no kspace), and the source from the GitHub tag `stable_22Jul2025_update3` because download.lammps.org no longer serves that tarball.

It runs directly on the login node, not as a Slurm job: Pod's queues are often full for hours, and the login node has 80 threads, all of which the build uses (`-j $(nproc)` at nice 10), 20-40 min. The script logs itself to `build-lammps.<date>.<host>.log` in this folder, the counterpart of the Expanse job log.

```bash
# on a Pod login node, from anywhere
bash ~/Documents/lammps_work/lammps_build/build_lammps_pod.sh
mpirun -n 1 ~/software/lammps/22Jul2025_update3/mpi-omp/gnu14/14.2.0/openmpi5/5.0.7/lammps-22Jul2025/build/lmp -h | head -3
```

A bare `lmp -h` fails on Pod's login node (Open MPI 5 asks for a PMIx transport key); run it through `mpirun -n 1` as above. The build of 2026-11-18 replaced the old `~/MD/lammps/lammps-22Jul2025/build/lmp` (22Jul2025 update 1 with only MOLECULE, RIGID, KSPACE, MC, OPENMP, EXTRA-FIX, EXTRA-DUMP).

Other machines: the decks only need the MOLECULE, EXTRA-FIX, EXTRA-PAIR, EXTRA-MOLECULE and RIGID packages plus OpenMP. On a Mac the LAMMPS GUI app's bundled `lmp` (LAMMPS_GUI.app/Contents/bin/lmp, with `-sf omp -pk omp N`) runs every deck for local smoke tests. Bridges-2 used its system LAMMPS module (see `docs/cluster_reference.md`).
