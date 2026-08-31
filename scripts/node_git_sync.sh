#!/bin/bash
# Shared best-effort "pull lammps_work on the compute node" helper.
# Sourced by run_lammps.sh / run_lammps_pod.sh / run_lammps_bridges.sh.
# Requires LAMMPS_WORK_DIR to be set by the caller.

# ── Pull latest scripts from GitHub before running ────────────────────────────
# Best-effort: never fails the job, never leaves the repo in a broken state.
#
# NOTE (2026-08-30): the previous `git pull --rebase --autostash` raced badly.
# A sweep batch starts ~41 jobs at once, all pulling the SAME working copy; they
# fought over .git locks and 22/287 jobs logged "no such ref was fetched" (the
# whole P=1.00 batch). A losing racer could also have been mid-rebase or holding
# someone else's autostash when the next one started. Fixes:
#   - flock serialises the repo access, so only one job touches .git at a time;
#   - a freshness stamp means the rest skip instead of re-fetching 41 times;
#   - fetch + `merge --ff-only` replaces rebase/autostash, so a dirty or diverged
#     tree is left exactly as-is rather than half-rewritten.
sync_lammps_work() {
    local stamp="${LAMMPS_WORK_DIR}/.git/.last_node_sync"
    local lock="${LAMMPS_WORK_DIR}/.git/.node_sync.lock"
    local max_age=300   # s; one sync per batch is plenty

    if ! command -v git &>/dev/null; then
        echo ">>> git not available on this node — skipping sync"; return 0
    fi
    if [ ! -d "${LAMMPS_WORK_DIR}/.git" ]; then
        echo ">>> $LAMMPS_WORK_DIR is not a git repo — skipping sync"; return 0
    fi

    # Fast path: someone already synced recently.
    if [ -f "$stamp" ] && [ $(( $(date +%s) - $(stat -c %Y "$stamp") )) -lt "$max_age" ]; then
        echo ">>> synced $(( $(date +%s) - $(stat -c %Y "$stamp") ))s ago by another job — skipping"
        return 0
    fi

    if ! command -v flock &>/dev/null; then
        echo ">>> flock not available — skipping sync (cannot serialise safely)"; return 0
    fi

    (
        # -w: don't queue 41 deep; a job that can't get the lock just runs.
        if ! flock -w 60 9; then
            echo ">>> another job holds the sync lock — skipping"; exit 0
        fi
        # Re-check under the lock: the holder we waited on probably just synced.
        if [ -f "$stamp" ] && [ $(( $(date +%s) - $(stat -c %Y "$stamp") )) -lt "$max_age" ]; then
            echo ">>> synced by the job ahead of us — skipping"; exit 0
        fi
        # Fetch the checked-out branch by name. A bare `git fetch origin` leaves
        # FETCH_HEAD ambiguous and was part of the old "no such ref" failure.
        local_branch=$(git -C "$LAMMPS_WORK_DIR" symbolic-ref --quiet --short HEAD 2>/dev/null || true)
        if [ -z "$local_branch" ]; then
            echo ">>> detached HEAD — skipping sync"; exit 0
        fi
        if ! timeout 120 git -C "$LAMMPS_WORK_DIR" fetch --quiet origin "$local_branch"; then
            echo ">>> fetch of '$local_branch' failed (no network / auth?) — running with the on-disk checkout"
            exit 0
        fi
        # Fast-forward only: refuses on local commits or a dirty tree instead of
        # rewriting them. Either way the checkout stays coherent for the run.
        if git -C "$LAMMPS_WORK_DIR" merge --ff-only --quiet FETCH_HEAD 2>/dev/null; then
            echo ">>> synced to $(git -C "$LAMMPS_WORK_DIR" rev-parse --short HEAD)"
        else
            echo ">>> cannot fast-forward (local commits or dirty tree) — running with the on-disk checkout"
        fi
        touch "$stamp"
    ) 9>"$lock"
}

