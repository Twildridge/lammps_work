#!/usr/bin/env bash
# =============================================================================
# install_lsync.sh — one-time repo setup for this machine (laptop, Expanse,
# Bridges-2, Pod). Run ONCE per machine, from anywhere, right after cloning:
#     bash /path/to/lammps_work/install.sh          (the same as running this file)
#
# No options are needed on any machine. It is idempotent: re-running it just
# refreshes everything to this clone's path.
#
# What it sets up:
#   0. nbstripout, registered FOR THIS CLONE — strips notebook output images
#      (plots, embedded PNGs) before they're committed, so re-running a notebook
#      never bloats the repo.  This is a per-clone git setting (lives in
#      .git/config, not in anything trackable) — .gitattributes tells git WHICH
#      files to filter, but each clone still has to register the filter once,
#      which is what this step does.  Skipped with a warning if `pip` isn't on
#      PATH; notebooks will still commit, just with their outputs embedded.
#   1. two shell aliases in ~/.zshrc or ~/.bashrc (whichever matches $SHELL):
#      `lsync` (scripts/git_sync.sh) and `lsbatch` (scripts/submit.sh: submit a batch
#      file by name from any folder, with the allocation of the cluster it is on)
#   2. a user-level Claude Code `/lsync` command that points at this clone, ONLY
#      on a machine that has Claude Code (a ~/.claude folder or a `claude` command).
#      Everywhere else, a cluster included, the step is skipped by itself.
#      --no-claude skips it on a machine that does have Claude Code.
#   3. ~/Documents/lammps_data/input_data, the folder the runner reads .data files from
#
# Both `lsync` entry points shell out to scripts/git_sync.sh in THIS clone, so
# there is a single implementation and nothing to keep in sync by hand.
#
# NOTE: the repo already ships .claude/commands/lsync.md, so if you open the
# lammps_work folder itself in Claude Code, /lsync works with no install at all.
# Step 2 is only needed when your Claude Code working directory is a PARENT of
# the clone (or somewhere else entirely).
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
SYNC="$REPO_ROOT/scripts/git_sync.sh"
# Claude Code is an optional AI coding tool; nothing in this repo needs it. Its /lsync command
# is installed only where Claude Code already is.
DO_CLAUDE=0
if [ -d "$HOME/.claude" ] || command -v claude >/dev/null 2>&1; then DO_CLAUDE=1; fi
[ "${1:-}" = "--no-claude" ] && DO_CLAUDE=0

[ -f "$SYNC" ] || { echo "ERROR: $SYNC not found — is this a full lammps_work clone?"; exit 1; }

echo "======================================"
echo "lammps_work one-time setup"
echo "clone: $REPO_ROOT"
echo "======================================"

# ── 0. nbstripout: strip notebook outputs before every commit ──────────────
# nbstripout --install writes filter.nbstripout.* / diff.ipynb.textconv into
# THIS repo's .git/config, pointing at whatever `python` resolves to right now
# on THIS machine. That local registration is exactly what a fresh `git clone`
# does not carry over — .gitattributes (committed, shared) says *.ipynb should
# be filtered; this step is what actually wires the filter up here.
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

# ── 1. shell alias ─────────────────────────────────────────────────────────
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
        printf '\n# lammps_work: %s (scripts/install_lsync.sh)\n%s\n' "$what" "$ALIAS" >> "$RC"
        echo ">>> added $name alias to $RC"
    fi
done

# ── 2. user-level Claude Code command (optional) ───────────────────────────
if [ "$DO_CLAUDE" -eq 1 ]; then
    CMD_SRC="$REPO_ROOT/.claude/commands/lsync.md"
    CMD_DIR="$HOME/.claude/commands"
    if [ -f "$CMD_SRC" ]; then
        mkdir -p "$CMD_DIR"
        # A real copy, not a symlink: some command loaders skip symlinks when
        # scanning for commands (isFile() is false for a symlink without an
        # extra resolve step), which can make /lsync silently undiscoverable.
        # Tradeoff: a `git pull` that changes lsync.md needs this script re-run
        # to pick it up — cheap and idempotent, so just do that when in doubt.
        cp "$CMD_SRC" "$CMD_DIR/lsync.md"
        echo ">>> copied $CMD_SRC -> $CMD_DIR/lsync.md"
    else
        echo ">>> skipped Claude Code command ($CMD_SRC not in this clone)"
    fi
else
    echo ">>> Claude Code not found on this machine (or --no-claude): /lsync command not installed; nothing else needs it"
fi

# ── 3. the folder for simulation input files (not tracked by git) ──────────
mkdir -p "$HOME/Documents/lammps_data/input_data"
echo ">>> input files go in $HOME/Documents/lammps_data/input_data"

echo "======================================"
echo "Done. Open a new terminal (or run:  source $RC ), then:"
echo "    lsync                    pull, commit and push in one step  (lsync \"my message\")"
echo "    lsbatch <name>.batch     submit a batch file from any folder (on a cluster)"
[ "$DO_CLAUDE" -eq 1 ] && echo "    /lsync                   the same sync, inside Claude Code"
echo "======================================"
