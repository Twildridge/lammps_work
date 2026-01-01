## Slurm Command Examples


#### CPU only
```bash
#!/bin/bash -l
#SBATCH --job-name=slab_support
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=120
#SBATCH --cpus-per-task=1
#SBATCH --partition=RM
#SBATCH --time=2:00:00
```

```bash
mpirun -n $SLURM_NTASKS \
    /opt/packages/LAMMPS/lammps-22Jul2025/build-RM-gcc13.3.1/lmp \
    -sf omp -pk omp $SLURM_CPUS_PER_TASK  \
    -var dataname $DATANAME \
    -var interaction $INTERACTION \
    -var epsSS $EPSSS \
    -var epsSP $EPSSP \
    -var nsteps $NSTEPS \
    -var oldsteps $OLDSTEPS \
    -var totsteps $TOTSTEPS \
    -in $LAMMPS_FILE
```


#### GPU
```bash
#!/bin/bash -l
#SBATCH --job-name=slab_gpu_test1
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=GPU
#SBATCH --gpus=v100-16:1
#SBATCH --time=1:00:00
```

Other options:
Format: --gpus=type:n
"n" is the number of GPUs

#SBATCH --gpus=v100-16:8
#SBATCH --gpus=v100-32:16

## Compiling

#### Pod

```bash
cmake ../cmake \
    -D PKG_KOKKOS=on \
    -D Kokkos_ARCH_NATIVE=yes \
    -D Kokkos_ARCH_VOLTA70=yes \
    -D Kokkos_ENABLE_CUDA=yes \
    -D Kokkos_ENABLE_OPENMP=yes \
    -D CMAKE_CXX_COMPILER=~/MD/lammps/lammps-22Jul2025/lib/kokkos/bin/nvcc_wrapper \
    -D PKG_MC=yes \
    -D PKG_MOLECULE=yes \
    -D PKG_KSPACE=yes \
    -D PKG_RIGID=yes \
    -D PKG_EXTRA-DUMP=yes \
    -D PKG_EXTRA-FIX=yes \
    -D BUILD_LAMMPS_GUI=no
```



