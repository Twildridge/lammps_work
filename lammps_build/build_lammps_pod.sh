#!/usr/bin/env bash
# Pod (CNSI, UCSB) copy of build_lammps.sh. Same LAMMPS version, same CMake options and the same
# -DPKG_* list as the Expanse build; only what differs between the clusters changes:
#   where         directly on the login node, not as a Slurm job: Pod's queues are often full for
#                 hours, the login node has 80 threads, and the build uses them all (-j $(nproc), at
#                 nice 10 so interactive users are not starved); about 20-40 min
#   modules       gnu14/14.2.0 + openmpi5/5.0.7 (the modules lib/site.sh loads for the jobs), hdf5/1.14.6
#   Kokkos arch   SKX (Xeon Gold 6148, Skylake) in place of ZEN2
#   FFT           KISS: Pod has no FFTW module or headers. The decks use no kspace_style, so it is unused.
#   source        download.lammps.org no longer serves the update3 tarball (only 22Jul2025 and update6),
#                 so the GitHub tag stable_22Jul2025_update3 is used; the download is skipped when the
#                 tarball is already there. The GitHub archive unpacks as lammps-stable_22Jul2025_update3/,
#                 renamed to lammps-22Jul2025/ to keep the Expanse layout.
# Binary:  ~/software/lammps/22Jul2025_update3/mpi-omp/gnu14/14.2.0/openmpi5/5.0.7/lammps-22Jul2025/build/lmp
# Run on a Pod login node, from anywhere:
#     bash ~/Documents/lammps_work/lammps_build/build_lammps_pod.sh
# It logs itself to lammps_build/build-lammps.<date>.<host>.log (the Pod counterpart of the Expanse
# build-lammps.o<jobid>.<node>.log) and ends by printing the version and the installed packages.
set -o pipefail
declare -xr COMPILER_MODULE='gnu14/14.2.0'
declare -xr MPI_MODULE='openmpi5/5.0.7'
declare -xr HDF5_MODULE='hdf5/1.14.6'
declare -xr LAMMPS_VERSION='22Jul2025_update3'
declare -xr LAMMPS_BUILD='mpi-omp'
declare -xr LAMMPS_ROOT_DIR="${HOME}/software/lammps"
declare -xr LAMMPS_TARBALL_URL="https://github.com/lammps/lammps/archive/refs/tags/stable_${LAMMPS_VERSION}.tar.gz"
declare -xr LAMMPS_INSTALL_DIR="${LAMMPS_ROOT_DIR}/${LAMMPS_VERSION}/${LAMMPS_BUILD}/${COMPILER_MODULE}/${MPI_MODULE}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${SCRIPT_DIR}/build-lammps.$(date +%Y%m%d).$(hostname -s).log"
exec > >(tee -a "$LOG") 2>&1
echo "=== build started: $(date)   log: $LOG"
type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh 2>/dev/null || source /etc/profile.d/modules.sh
module purge
module load "${COMPILER_MODULE}"
module load "${MPI_MODULE}"
module load "${HDF5_MODULE}"
module list
cmake --version | head -1      # system cmake 3.26.5; LAMMPS 22Jul2025 needs >= 3.20
gcc --version | head -1
mpirun --version | head -1
echo "threads: $(nproc)"
mkdir -p "${LAMMPS_INSTALL_DIR}" || exit 1
cd "${LAMMPS_INSTALL_DIR}" || exit 1
if [ ! -f "lammps-${LAMMPS_VERSION}.tar.gz" ]; then
    curl -sSL -o "lammps-${LAMMPS_VERSION}.tar.gz" "${LAMMPS_TARBALL_URL}" || exit 1
fi
[ -d "lammps-22Jul2025" ] || tar -xf "lammps-${LAMMPS_VERSION}.tar.gz" || exit 1
[ -d "lammps-stable_${LAMMPS_VERSION}" ] && mv "lammps-stable_${LAMMPS_VERSION}" "lammps-22Jul2025"
cd "lammps-22Jul2025" || exit 1
find . -exec touch {} +
mkdir -p build
cd build || exit 1
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
  -DKokkos_ARCH_SKX=yes \
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
  -DFFT=KISS \
  -DCMAKE_CXX_FLAGS="-std=c++11" \
  -DEXTERNAL_KOKKOS=OFF \
  ../cmake || { echo "=== CMAKE FAILED"; exit 1; }
nice -n 10 cmake --build . -- -j "$(nproc)" || { echo "=== BUILD FAILED: $(date)"; exit 1; }
echo "=== build finished: $(date)"
mpirun -n 1 ./lmp -h | sed -n '1,3p;/Installed packages/,/^List of/p' | grep -v '^List of'
echo "=== binary: $(pwd)/lmp"
