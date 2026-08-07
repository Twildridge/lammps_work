#!/usr/bin/env bash
# =============================================================================
# install_lsync.sh — one-time repo setup for this machine (laptop, Expanse,
# Bridges-2, Pod). Run ONCE per machine, from anywhere, right after cloning:
#     bash /path/to/lammps_work/scripts/install_lsync.sh
#
# It is idempotent: re-running it just refreshes everything to this clone's path.
#
# What it sets up:
#   0. nbstripout, registered FOR THIS CLONE — strips notebook output images
#      (plots, embedded PNGs) before they're committed, so re-running a notebook
#      never bloats the repo.  This is a per-clone git setting (lives in
#      .git/config, not in anything trackable) — .gitattributes tells git WHICH
#      files to filter, but each clone still has to register the filter once,
#      which is what this step does.  Skipped with a warning if `pip` isn't on
#      PATH; notebooks will still commit, just with their outputs embedded.
#   1. a shell alias  `lsync`  in ~/.zshrc or ~/.bashrc (whichever matches $SHELL)
#   2. optionally, a user-level Claude Code `/lsync` command that points at this
#      clone (skip with --no-claude; not useful on a cluster)
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
DO_CLAUDE=1
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
ALIAS="alias lsync='bash \"$REPO_ROOT/scripts/git_sync.sh\"'"

touch "$RC"
if grep -q "^alias lsync=" "$RC"; then
    # rewrite in place so a moved/re-cloned repo updates cleanly
    tmp="$(mktemp)"
    grep -v "^alias lsync=" "$RC" > "$tmp"
    printf '%s\n' "$ALIAS" >> "$tmp"
    mv "$tmp" "$RC"
    echo ">>> updated existing lsync alias in $RC"
else
    printf '\n# lammps_work one-command git sync (scripts/install_lsync.sh)\n%s\n' "$ALIAS" >> "$RC"
    echo ">>> added lsync alias to $RC"
fi

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
    echo ">>> skipped Claude Code command (--no-claude)"
fi

echo "======================================"
echo "Done.  Run:  source $RC"
echo "Then:        lsync              (or  lsync \"my message\")"
echo "In Claude Code:  /lsync"
echo "======================================"
