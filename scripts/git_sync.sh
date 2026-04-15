#!/usr/bin/env bash
# =============================================================================
# git_sync.sh — one-command commit + push for lammps_work
#
# Usage (from anywhere):
#   lsync              (if you added the alias below to ~/.zshrc)
#   lsync "my message" (custom commit message)
#
# What it does:
#   1. cd to the lammps_work repo root
#   2. git pull --rebase  (pick up any cluster pushes first)
#   3. Stage all changes (git add -A), respecting .gitignore
#   4. Auto-generate a commit message from changed file names (or use your own)
#   5. git push
#
# For clusters (Expanse / Bridges-2 / Pod):
#   Just run:  cd ~/Documents/lammps_work && git pull
#   Optionally add this to the top of run_lammps.sh or run_lammps_bridges.sh
#   so every job automatically gets the latest scripts before running.
#
# Add this alias to ~/.zshrc so you can call it from any terminal:
#   alias lsync='bash ~/docs/grad_research/lammps/lammps_work/scripts/git_sync.sh'
# Then run:  source ~/.zshrc
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
