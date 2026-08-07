#!/usr/bin/env bash
# =============================================================================
# git_sync.sh — one-command commit + push for lammps_work
#
# Usage (from anywhere):
#   lsync              (shell alias, or /lsync inside Claude Code)
#   lsync "my message" (custom commit message)
#
# What it does:
#   1. cd to the lammps_work repo root
#   2. git pull --rebase  (pick up any cluster pushes first)
#   3. Stage all changes (git add -A), respecting .gitignore
#   4. Auto-generate a commit message from changed file names (or use your own)
#   5. git push
#
# -----------------------------------------------------------------------------
# SETUP — once per machine (laptop, Expanse, Bridges-2, Pod), from any clone:
#
#     bash scripts/install_lsync.sh
#
# That adds the `lsync` alias to ~/.zshrc or ~/.bashrc and links the Claude Code
# /lsync command.  It is idempotent and derives every path from the clone it
# lives in, so it works for any collaborator without editing anything.
# On a cluster (no Claude Code) use:  bash scripts/install_lsync.sh --no-claude
#
# Entry points, all calling THIS script so they can never drift apart:
#     lsync            shell alias
#     /lsync           Claude Code  (.claude/commands/lsync.md, shipped in-repo)
# -----------------------------------------------------------------------------
#
# WHY RUN IT ON THE CLUSTER TOO: edits made directly on Expanse (e.g. tweaking
# NSTEPS in a .batch before an sbatch) are otherwise uncommitted, and
# triaxial_compression.batch self-syncs with `git pull --rebase --autostash` at
# job start.  Running lsync there commits those edits instead of leaving them to
# be stashed and silently reapplied — and makes them visible on the Mac side.
# =============================================================================

set -euo pipefail

# ── Locate repo root (works wherever this script is called from) ───────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"   # lammps_work/

cd "$REPO_ROOT"

echo "======================================"
echo "lammps_work git sync"
echo "Repo: $REPO_ROOT"
echo "======================================"

# ── 1. Pull (rebase) to incorporate any cluster commits ────────────────────
echo ">>> Pulling latest from remote..."
git pull --rebase --autostash || {
    echo "ERROR: git pull failed. Resolve conflicts manually, then re-run."
    exit 1
}

# ── 2. Check for changes ───────────────────────────────────────────────────
if git diff --quiet && git diff --staged --quiet && [ -z "$(git ls-files --others --exclude-standard)" ]; then
    echo ">>> Nothing to commit. Already up to date."
    exit 0
fi

# ── 3. Stage everything (respects .gitignore) ─────────────────────────────
git add -A

# Warn (don't block) if a notebook is about to be committed WITH its plot
# outputs embedded, i.e. nbstripout isn't registered on this machine yet.
# See scripts/install_lsync.sh step 0 for why this is a per-clone setting.
if git diff --staged --name-only | grep -q '\.ipynb$' \
   && [ -z "$(git config --get filter.nbstripout.clean 2>/dev/null)" ]; then
    echo "!! WARNING: nbstripout is not set up on this machine — the notebook(s)"
    echo "   below will be committed WITH their embedded plot images (large diffs)."
    git diff --staged --name-only | grep '\.ipynb$' | sed 's/^/     /'
    echo "   Fix once, then re-run lsync:  bash \"$REPO_ROOT/scripts/install_lsync.sh\""
fi

# ── 4. Build commit message ────────────────────────────────────────────────
if [ $# -ge 1 ] && [ -n "$1" ]; then
    # Custom message passed as argument
    COMMIT_MSG="$1"
else
    # Auto-generate from changed files
    CHANGED=$(git diff --staged --name-only | head -8 | tr '\n' ', ' | sed 's/,$//')
    NFILES=$(git diff --staged --name-only | wc -l | tr -d ' ')
    TIMESTAMP=$(date "+%Y-%m-%d %H:%M")
    if [ "$NFILES" -le 3 ]; then
        COMMIT_MSG="updated: ${CHANGED} [${TIMESTAMP}]"
    else
        COMMIT_MSG="updated ${NFILES} files: ${CHANGED}... [${TIMESTAMP}]"
    fi
fi

echo ">>> Committing: $COMMIT_MSG"
git commit -m "$COMMIT_MSG"

# ── 5. Push ────────────────────────────────────────────────────────────────
echo ">>> Pushing to origin..."
git push

echo "======================================"
echo "Done! lammps_work is up to date on GitHub."
echo "Run 'git pull' on Expanse / Bridges-2 / Pod to sync clusters."
echo "======================================"
