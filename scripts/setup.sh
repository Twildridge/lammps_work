#!/usr/bin/env bash
# =============================================================================
# setup.sh — the one-time setup of a clone on this machine (laptop, Expanse or Pod). Run it as
#
#     bash install.sh            (from the repo root; install.sh just runs this file)
#
# then open a new terminal. It needs no options and is idempotent: run it again after moving the
# clone, after a git pull, or whenever something seems off; it refreshes everything to this
# clone's path and skips what is already there.
#
# On every machine
#   1. nbstripout, registered for THIS clone: notebook outputs (plots, PNGs) are stripped before a
#      commit, so re-running a notebook never bloats the repo. .gitattributes (committed) says which
#      files to filter; each clone still has to register the filter once in its own .git/config,
#      which git clone does not carry over. That is this step.
#   2. two shell aliases in ~/.bashrc or ~/.zshrc (whichever matches $SHELL): lsync (scripts/git_sync.sh,
#      pull + commit + push in one step) and lsbatch (scripts/submit.sh, submit a .batch by name from
#      any folder with the allocation of the cluster it is on)
#   3. the Claude Code /lsync command, only on a machine that has Claude Code (--no-claude skips it)
#   4. ~/Documents/lammps_data/input_data, the folder the runner reads .data files from
# On a cluster (scripts/lib/site.sh recognises Expanse and Pod; a laptop skips this part)
#   5. what the jobs need that a clone does not bring:
#      Pod      the conda environment lammps_analysis (python 3.11, numpy, scipy, matplotlib, pandas),
#               which the jobs activate for the plotters: created if missing, verified either way;
#               the LAMMPS binary of lib/site.sh: Pod has no site LAMMPS, so it is built here, on the
#               login node, if missing (lammps_build/build_lammps_pod.sh, about 10 min; --no-build skips it)
#      Expanse  nothing to install: the jobs load anaconda3 themselves and lib/site.sh points at a
#               shared binary, which is checked below
#   6. scripts/cluster_check.sh: does this machine match lib/site.sh (binary and its packages, modules,
#      python, scratch, cores per node, Slurm partitions)? It changes nothing and can be re-run alone
#      at any time. Its last line says whether the setup is complete.
#
# Options:  --no-claude   skip step 3 on a machine that has Claude Code
#           --no-build    Pod: do not build LAMMPS when it is missing (just say how)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
SYNC="$REPO_ROOT/scripts/git_sync.sh"
[ -f "$SYNC" ] || { echo "ERROR: $SYNC not found — is this a full lammps_work clone?"; exit 1; }

DO_CLAUDE=0; DO_BUILD=1
if [ -d "$HOME/.claude" ] || command -v claude >/dev/null 2>&1; then DO_CLAUDE=1; fi
for arg in "$@"; do
    case "$arg" in
        --no-claude) DO_CLAUDE=0 ;;
        --no-build)  DO_BUILD=0 ;;
        *) echo "ERROR: unknown option $arg (options: --no-claude, --no-build)"; exit 1 ;;
    esac
done

# Which machine: lib/site.sh names the cluster (a laptop comes out as expanse; sbatch tells them apart).
# site.sh is sourced in a subshell with set -e/-u off: it and Lmod read variables that may be unset.
CLUSTER_NAME="$( set +eu; source "$REPO_ROOT/scripts/lib/site.sh" >/dev/null 2>&1; echo "${CLUSTER:-}" )"
ON_CLUSTER=0; command -v sbatch >/dev/null 2>&1 && ON_CLUSTER=1
if [ "$ON_CLUSTER" = 1 ]; then WHERE="cluster $CLUSTER_NAME"; else WHERE="this computer (not a cluster)"; CLUSTER_NAME=""; fi

echo "======================================"
echo "lammps_work one-time setup"
echo "clone:   $REPO_ROOT"
echo "machine: $WHERE"
echo "======================================"

# ── 1. nbstripout: strip notebook outputs before every commit ──────────────
export PATH="$HOME/.local/bin:$PATH"            # where pip --user puts nbstripout
if [ "$CLUSTER_NAME" = pod ]; then
    # Pod: pip is miniconda's only in a shell whose ~/.bashrc ran conda init; load the module explicitly.
    # Lmod's init file and the module function read unset variables, which set -u would make fatal.
    set +u
    type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh 2>/dev/null || true
    module load miniconda 2>/dev/null || echo "!! module load miniconda failed; pip is $(command -v pip || echo missing)"
    set -u
fi
if ! command -v nbstripout >/dev/null 2>&1; then
    echo ">>> nbstripout not found — installing via pip..."
    pip install --quiet --user nbstripout 2>/dev/null || pip3 install --quiet --user nbstripout 2>/dev/null || true
fi
if command -v nbstripout >/dev/null 2>&1; then
    ( cd "$REPO_ROOT" && nbstripout --install )
    echo ">>> nbstripout registered for this clone (notebook outputs won't be committed)"
else
    echo "!! could not install nbstripout automatically (no pip on PATH)."
    echo "   Run by hand:  pip install nbstripout && cd \"$REPO_ROOT\" && nbstripout --install"
    echo "   Until then, notebooks will commit WITH their embedded plot images."
fi

# ── 2. shell aliases ───────────────────────────────────────────────────────
case "${SHELL:-}" in
    *zsh) RC="$HOME/.zshrc" ;;
    *)    RC="$HOME/.bashrc" ;;
esac
touch "$RC"
for pair in "lsync:git_sync.sh:one-command git sync" "lsbatch:submit.sh:submit a .batch by name from any folder (Expanse or Pod)"; do
    IFS=: read -r name script what <<<"$pair"
    ALIAS="alias $name='bash \"$REPO_ROOT/scripts/$script\"'"
    if grep -q "^alias $name=" "$RC"; then
        # rewrite in place so a moved/re-cloned repo updates cleanly
        tmp="$(mktemp)"
        grep -v "^alias $name=" "$RC" > "$tmp"
        printf '%s\n' "$ALIAS" >> "$tmp"
        mv "$tmp" "$RC"
        echo ">>> updated existing $name alias in $RC"
    else
        printf '\n# lammps_work: %s (scripts/setup.sh)\n%s\n' "$what" "$ALIAS" >> "$RC"
        echo ">>> added $name alias to $RC"
    fi
done

# ── 3. Claude Code /lsync command (optional) ───────────────────────────────
# Claude Code is an optional AI coding tool; nothing in this repo needs it. The repo ships
# .claude/commands/lsync.md, so inside the lammps_work folder /lsync works with no install; this
# copy makes it work from a parent folder too. A real copy, not a symlink: some command loaders
# skip symlinks. A git pull that changes lsync.md needs this script re-run (cheap, idempotent).
if [ "$DO_CLAUDE" -eq 1 ]; then
    CMD_SRC="$REPO_ROOT/.claude/commands/lsync.md"
    CMD_DIR="$HOME/.claude/commands"
    if [ -f "$CMD_SRC" ]; then
        mkdir -p "$CMD_DIR"
        cp "$CMD_SRC" "$CMD_DIR/lsync.md"
        echo ">>> copied $CMD_SRC -> $CMD_DIR/lsync.md"
    else
        echo ">>> skipped Claude Code command ($CMD_SRC not in this clone)"
    fi
else
    echo ">>> Claude Code not found on this machine (or --no-claude): /lsync command not installed; nothing else needs it"
fi

# ── 4. the folder for simulation input files (not tracked by git) ──────────
mkdir -p "$HOME/Documents/lammps_data/input_data"
echo ">>> input files go in $HOME/Documents/lammps_data/input_data"

# ── 5. what the jobs need on this cluster ──────────────────────────────────
SETUP_OK=1
if [ "$CLUSTER_NAME" = pod ]; then
    echo ">>> Pod: conda environment for the plotters, and the LAMMPS binary"
    # a subshell without -e/-u: Lmod, the module function and conda.sh read unset variables
    ( set +eu
      source "$REPO_ROOT/scripts/lib/site.sh" || exit 1
      type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh 2>/dev/null || source /etc/profile.d/modules.sh 2>/dev/null \
          || { echo "!! no 'module' command: is this a Pod login node?"; exit 1; }
      module load "$PYTHON_MODULE" || { echo "!! module load $PYTHON_MODULE failed"; exit 1; }
      # shellcheck disable=SC1091
      source "$(conda info --base)/etc/profile.d/conda.sh"
      ENV_NAME="$PYTHON_CONDA_ENV"
      PKGS=(python=3.11 numpy scipy matplotlib pandas)
      if conda env list 2>/dev/null | awk 'NF && $1 !~ /^#/ {print $1}' | grep -qx "$ENV_NAME"; then
          echo "    conda env $ENV_NAME exists"
      else
          echo "    creating conda env $ENV_NAME (${PKGS[*]}), a few minutes..."
          conda create -n "$ENV_NAME" "${PKGS[@]}" -y || { echo "!! conda create failed"; exit 1; }
      fi
      if conda run -n "$ENV_NAME" python -c "import numpy, scipy, matplotlib, pandas" >/dev/null 2>&1; then
          echo "    $ENV_NAME: numpy, scipy, matplotlib, pandas import"
      else
          echo "    $ENV_NAME is missing a package; installing ${PKGS[*]:1}"
          conda install -n "$ENV_NAME" "${PKGS[@]:1}" -y || { echo "!! conda install failed"; exit 1; }
          conda run -n "$ENV_NAME" python -c "import numpy, scipy, matplotlib, pandas" \
              || { echo "!! $ENV_NAME still cannot import numpy, scipy, matplotlib, pandas"; exit 1; }
      fi
      if [ -x "$LMP_BIN" ]; then
          echo "    LAMMPS binary: $LMP_BIN"
      elif [ "$DO_BUILD" = 1 ]; then
          echo "    no LAMMPS binary at $LMP_BIN"
          echo "    building it now on this login node with every thread (about 10 min; Ctrl-C and"
          echo "    re-run with --no-build to postpone). Log: lammps_build/build-lammps.<date>.<host>.log"
          bash "$REPO_ROOT/lammps_build/build_lammps_pod.sh" >/dev/null 2>&1 \
              || { echo "!! the build failed: see the newest lammps_build/build-lammps.*.log"; exit 1; }
          [ -x "$LMP_BIN" ] && echo "    LAMMPS binary: $LMP_BIN" \
              || { echo "!! the build finished but $LMP_BIN is not there: see lammps_build/build-lammps.*.log"; exit 1; }
      else
          echo "!! no LAMMPS binary at $LMP_BIN (--no-build). Build it on this login node with:"
          echo "       bash $REPO_ROOT/lammps_build/build_lammps_pod.sh"
          exit 1
      fi
    ) || SETUP_OK=0
elif [ "$CLUSTER_NAME" = expanse ]; then
    echo ">>> Expanse: nothing to install (the jobs load anaconda3 and a shared LAMMPS binary)"
    LMP_SHARED="$( set +eu; source "$REPO_ROOT/scripts/lib/site.sh" >/dev/null 2>&1; echo "${LMP_BIN:-}" )"
    if [ -x "$LMP_SHARED" ]; then
        echo "    LAMMPS binary: $LMP_SHARED"
    else
        echo "!! the LAMMPS binary of lib/site.sh is not readable from this account: $LMP_SHARED"
        echo "   Either ask its owner for read access, or build your own with"
        echo "       cd $REPO_ROOT/lammps_build && sbatch build_lammps.sh      (~1 h; edit --account first)"
        echo "   and put its path in the expanse block of scripts/lib/site.sh (LMP_BIN)."
        SETUP_OK=0
    fi
fi

# ── 6. on a cluster: does the machine match lib/site.sh? ───────────────────
if [ "$ON_CLUSTER" = 1 ]; then
    echo ">>> scripts/cluster_check.sh (re-run it alone at any time)"
    bash "$REPO_ROOT/scripts/cluster_check.sh" | sed 's/^/    /' || SETUP_OK=0
fi

echo "======================================"
if [ "$SETUP_OK" = 1 ]; then
    echo "Done. Open a new terminal (or run:  source $RC ), then:"
else
    echo "Done, with the problem(s) marked !! or PROBLEM above. Fix them and run  bash install.sh  again."
    echo "Meanwhile, open a new terminal (or run:  source $RC ) for:"
fi
echo "    lsync                    pull, commit and push in one step  (lsync \"my message\")"
echo "    lsbatch <name>.batch     submit a batch file from any folder (on a cluster)"
[ "$DO_CLAUDE" -eq 1 ] && echo "    /lsync                   the same sync, inside Claude Code"
echo "======================================"
[ "$SETUP_OK" = 1 ]
