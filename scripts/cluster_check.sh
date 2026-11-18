#!/bin/bash
# cluster_check.sh -- does this cluster match what lib/site.sh says about it? Run it on a login
# node before the first job on a cluster, or after the cluster's software changed:
#
#   bash scripts/cluster_check.sh
#
# It changes nothing. Each line is OK or a PROBLEM with the setting of lib/site.sh to correct.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/site.sh" || exit 1
BAD=0
ok()  { echo "  OK       $*"; }
bad() { echo "  PROBLEM  $*"; BAD=$((BAD + 1)); }

echo "cluster: $CLUSTER   (host $(hostname -f 2>/dev/null || hostname), Slurm cluster name '$(scontrol show config 2>/dev/null | awk '/^ClusterName/ {print $3}')')"
echo "bash $BASH_VERSION   $(sbatch --version 2>/dev/null || echo 'sbatch: not found')"
if [ "${BASH_VERSINFO[0]}" -ge 4 ]; then ok "bash >= 4 (the runner uses associative arrays)"; else bad "bash ${BASH_VERSINFO[0]} is too old for run_lammps.sh (needs 4+)"; fi

echo "paths"
if [ "$LAMMPS_WORK_DIR" = "$HOME/Documents/lammps_work" ]; then ok "repo at $LAMMPS_WORK_DIR"; else bad "repo is at $LAMMPS_WORK_DIR; every .batch sources \$HOME/Documents/lammps_work/scripts/batch_prologue.sh"; fi
if [ -x "$LMP_BIN" ]; then ok "LMP_BIN $LMP_BIN"; else bad "LMP_BIN not found or not executable: $LMP_BIN$([ "$CLUSTER" = pod ] && echo " (build it: bash lammps_build/build_lammps_pod.sh)")"; fi
if [ -d "$INPUT_DATA_DIR" ]; then ok "INPUT_DATA_DIR $INPUT_DATA_DIR ($(ls "$INPUT_DATA_DIR" 2>/dev/null | grep -c '\.data$') .data files)"; else bad "INPUT_DATA_DIR missing: $INPUT_DATA_DIR (mkdir -p it, then copy the .data files in)"; fi
if mkdir -p "$SCRATCH_DIR/lammps_trajectories" 2>/dev/null && [ -w "$SCRATCH_DIR/lammps_trajectories" ]; then ok "SCRATCH_DIR $SCRATCH_DIR is writable"; else bad "SCRATCH_DIR not writable: $SCRATCH_DIR"; fi
mkdir -p "$RUNS_DIR" 2>/dev/null && ok "RUNS_DIR $RUNS_DIR" || bad "cannot create RUNS_DIR $RUNS_DIR"

echo "modules and the binary"
if ( load_lammps_modules ) >/dev/null 2>&1; then
    ok "module load ${LAMMPS_MODULES[*]}"
    ( load_lammps_modules >/dev/null 2>&1
      echo "           mpirun: $(command -v mpirun || echo 'not found')  ($(mpirun --version 2>/dev/null | head -1))"
      if [ -x "$LMP_BIN" ]; then
          # -h lists the build: version, MPI, installed packages (the decks need MOLECULE, EXTRA-FIX, OPENMP, RIGID)
          HELP=$("$LMP_BIN" -h 2>/dev/null)
          # on Pod a bare `lmp -h` dies on the login node (Open MPI 5 PMIx transport key); a launcher fixes it
          [ -n "$HELP" ] || HELP=$(mpirun -n 1 "$LMP_BIN" -h 2>/dev/null)
          echo "$HELP" | grep -m1 'Large-scale Atomic' | sed 's/^/           /'
          for pkg in MOLECULE OPENMP EXTRA-FIX; do
              echo "$HELP" | grep -qw "$pkg" && echo "  OK       LAMMPS package $pkg" || echo "  PROBLEM  LAMMPS package $pkg not in $LMP_BIN -h"
          done
      fi )
else
    bad "module load ${LAMMPS_MODULES[*]} failed (LAMMPS_MODULES in lib/site.sh); module avail shows what exists"
fi
if ( load_python_module && python -c "import numpy, scipy, matplotlib" ) >/dev/null 2>&1; then
    ok "python for the plotters: numpy, scipy, matplotlib import"
else
    bad "python for the plotters: load_python_module of lib/site.sh failed or numpy / scipy / matplotlib is missing"
fi

echo "Slurm"
if command -v sinfo >/dev/null 2>&1; then
    echo "  partitions (name, availability, time limit, nodes, cores per node, memory MB):"
    sinfo -h -o "    %P %a %l %D %c %m" | sort -u
    # one row per node (-N): per-partition rows would weigh a 14-node gpu partition like the 62-node batch one
    CORES=$(sinfo -h -N -o "%c" | sort -n | uniq -c | sort -rn | awk 'NR==1 {print $2}')
    if [ "$CORES" = "$CORES_PER_NODE" ]; then ok "most nodes have $CORES cores = CORES_PER_NODE"; else bad "most nodes have $CORES cores, lib/site.sh says CORES_PER_NODE=$CORES_PER_NODE"; fi
    echo "  my associations (account, partition, QOS):"
    sacctmgr -nP show assoc user="$USER" format=account,partition,qos 2>/dev/null | sed 's/^/    /' | head -8
    echo "  limits worth knowing before a sweep: scontrol show partition | grep -E 'PartitionName|MaxTime|MaxNodes'; sacctmgr show qos format=name,maxwall,maxtrespu,maxsubmitpu"
else
    bad "sinfo not found: is this a login node?"
fi

echo "batch files with an allocation for $CLUSTER:"
for f in "$LAMMPS_WORK_DIR"/simulations/*/*.batch; do
    [ -L "$f" ] && continue
    batch_header "$f" >/dev/null 2>&1 && echo "    ${f#"$LAMMPS_WORK_DIR"/simulations/}"
done
echo
if [ "$BAD" = 0 ]; then echo "No problems found. A DRY_RUN=1 run of run_lammps.sh inside a short job is the next check."; else echo "$BAD problem(s): fix the $CLUSTER block of scripts/lib/site.sh (or the cluster setup) before submitting."; fi
exit $((BAD > 0))
