#!/usr/bin/env bash
#SBATCH --job-name=build-lammps
#SBATCH --account=csb197
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%x.o%j.%N
declare -xr COMPILER_MODULE='gcc/10.2.0'
declare -xr MPI_MODULE='openmpi/4.1.3'
declare -xr CMAKE_MODULE='cmake/3.21.4'
declare -xr FFTW_MODULE='fftw/3.3.10'
declare -xr HDF5_MODULE='hdf5/1.10.7'
declare -xr LAMMPS_VERSION='22Jul2025_update3'
declare -xr LAMMPS_BUILD='mpi-omp'
declare -xr LAMMPS_ROOT_DIR="${HOME}/software/lammps"
declare -xr LAMMPS_ROOT_URL='https://download.lammps.org/tars'
declare -xr LAMMPS_INSTALL_DIR="${LAMMPS_ROOT_DIR}/${LAMMPS_VERSION}/${LAMMPS_BUILD}/${SOFTWARE_MODULE}/${COMPILER_MODULE}/${MPI_MODULE}"
module reset
module load "${COMPILER_MODULE}"
module load "${MPI_MODULE}"
module load "${CMAKE_MODULE}"
module load "${FFTW_MODULE}"
module load "${HDF5_MODULE}"
module list
printenv
mkdir -p "${LAMMPS_INSTALL_DIR}"
cd "${LAMMPS_INSTALL_DIR}"
wget "${LAMMPS_ROOT_URL}/lammps-${LAMMPS_VERSION}.tar.gz"
tar -xf "lammps-${LAMMPS_VERSION}.tar.gz"
cd "lammps-22Jul2025"
find . -exec touch {} +
mkdir -p build
cd build
cmake \
  -DBUILD_SHARED_LIBS=ON \
  -DLAMMPS_EXCEPTIONS=ON \
  -DBUILD_MPI=ON \
  -DBUILD_OMP=ON \
  -DBUILD_LIB=ON \
  -DWITH_JPEG=OFF \
  -DWITH_PNG=OFF \
  -DWITH_FFMPEG=OFF \
  -DPKG_ASPHERE=ON \
  -DPKG_BODY=ON \
  -DPKG_CLASS2=ON \
  -DPKG_COLLOID=ON \
  -DPKG_COMPRESS=ON \
  -DPKG_CORESHELL=ON \
  -DPKG_DIPOLE=ON \
  -DPKG_GRANULAR=ON \
  -DPKG_KSPACE=ON \
  -DPKG_KOKKOS=ON \
  -DKokkos_ARCH_ZEN2=yes \
  -DKokkos_ENABLE_OPENMP=yes \
  -DPKG_MANYBODY=ON \
  -DPKG_MC=ON \
  -DPKG_MISC=ON \
  -DPKG_MOLECULE=ON \
  -DPKG_MPIIO=ON \
  -DPKG_OPENMP=ON \
  -DPKG_PERI=ON \
  -DPKG_POEMS=ON \
  -DPKG_PYTHON=OFF \
  -DPKG_QEQ=ON \
  -DPKG_REPLICA=ON \
  -DPKG_RIGID=ON \
  -DPKG_SHOCK=ON \
  -DPKG_SNAP=ON \
  -DPKG_SPIN=ON \
  -DPKG_SRD=ON \
  -DPKG_VORONOI=ON \
  -DPKG_REAXFF=ON \
  -DPKG_COLVARS=ON \
  -DPKG_EXTRA-FIX=ON \
  -DPKG_EXTRA-PAIR=ON \
  -DPKG_EXTRA-MOLECULE=ON \
  -DFFT=FFTW3 \
  -DCMAKE_CXX_FLAGS="-std=c++11" \
  -DEXTERNAL_KOKKOS=OFF \
  ../cmake
cmake --build . -- -j ${SLURM_CPUS_PER_TASK}



