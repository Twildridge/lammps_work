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

Other machines: the decks only need the MOLECULE, EXTRA-FIX, EXTRA-PAIR, EXTRA-MOLECULE and RIGID packages plus OpenMP. On a Mac the LAMMPS GUI app's bundled `lmp` (LAMMPS_GUI.app/Contents/bin/lmp, with `-sf omp -pk omp N`) runs every deck for local smoke tests. Bridges-2 and Pod use their system LAMMPS modules (see `docs/cluster_reference.md`).
