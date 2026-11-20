# run_record.sh -- the run log: one JSON record per job, published to GitHub (lib/run_record.py has
# the format and builds the README tables). Sourced by run_lammps.sh, continue_sim.sh and git_sync.sh;
# needs LAMMPS_WORK_DIR and RUNS_DIR (lib/site.sh).
#
#   run_record_start <work_dir>            when the job starts (after write_run_info)
#   run_record_end   <work_dir> <status>   when LAMMPS has ended: complete | partial | failed
#   publish_run_records                    push every record that GitHub does not have yet
#
# Everything here is best-effort: a failure prints one line and never changes the job's outcome.
# RUN_RECORDS=0 turns it off.
#
# The records live in $RUNS_DIR/.run_records/ on the cluster, not in the working tree. Publishing
# builds a commit on top of origin's branch with git plumbing (a private index, commit-tree, push of
# that commit), so the checkout that running jobs read their scripts from is never touched and a
# dirty or diverged tree does not matter. The records arrive in the working tree with the next pull.
# A push that fails (no network or credentials on a compute node) leaves the records where they are:
# the next job, `lsbatch` or `lsync` on that cluster publishes them.

RUN_RECORD_STORE="${RUN_RECORD_STORE:-$RUNS_DIR/.run_records}"

_run_record_py() {
    command -v python3 >/dev/null 2>&1 || { echo ">>> run record: no python3, skipped"; return 1; }
    python3 "$LAMMPS_WORK_DIR/scripts/lib/run_record.py" "$@"
}

run_record_start() {
    [ "${RUN_RECORDS:-1}" = 0 ] && return 0
    ( _run_record_py start "$RUN_RECORD_STORE" "$1" && publish_run_records ) 2>&1 | sed 's/^/>>> /' || true
    return 0
}

run_record_end() {
    [ "${RUN_RECORDS:-1}" = 0 ] && return 0
    ( _run_record_py end "$RUN_RECORD_STORE" "$1" "$2" && publish_run_records ) 2>&1 | sed 's/^/>>> /' || true
    return 0
}

publish_run_records() {
    [ "${RUN_RECORDS:-1}" = 0 ] && return 0
    local repo="$LAMMPS_WORK_DIR" store="$RUN_RECORD_STORE"
    ls "$store"/*.json >/dev/null 2>&1 || return 0
    command -v git >/dev/null 2>&1 && [ -d "$repo/.git" ] || return 0
    local to=()
    command -v timeout >/dev/null 2>&1 && to=(timeout 45)
    (
        if command -v flock >/dev/null 2>&1; then
            flock -w 120 9 || { echo "run records: the repository is busy, not published now"; exit 0; }
        fi
        branch=$(git -C "$repo" symbolic-ref --quiet --short HEAD 2>/dev/null) || branch=main
        export GIT_AUTHOR_NAME="${GIT_AUTHOR_NAME:-$(git -C "$repo" config user.name || echo run-log)}"
        export GIT_AUTHOR_EMAIL="${GIT_AUTHOR_EMAIL:-$(git -C "$repo" config user.email || echo run-log@localhost)}"
        export GIT_COMMITTER_NAME="$GIT_AUTHOR_NAME" GIT_COMMITTER_EMAIL="$GIT_AUTHOR_EMAIL"
        export GIT_TERMINAL_PROMPT=0
        idx=$(mktemp "${TMPDIR:-/tmp}/run_record_index.XXXXXX") || exit 0
        trap 'rm -f "$idx"' EXIT
        for attempt in 1 2; do
            ${to[@]+"${to[@]}"} git -C "$repo" fetch --quiet origin "$branch" 2>/dev/null \
                || { echo "run records: fetch failed, not published now (they stay in $store)"; exit 0; }
            base=$(git -C "$repo" rev-parse FETCH_HEAD) || exit 0
            GIT_INDEX_FILE="$idx" git -C "$repo" read-tree "$base" || exit 0
            n=0
            for f in "$store"/*.json; do
                name="runs/$(basename "$f")"
                blob=$(git -C "$repo" hash-object -w "$f") || continue
                [ "$blob" = "$(git -C "$repo" rev-parse --quiet --verify "${base}:${name}" 2>/dev/null)" ] && continue
                GIT_INDEX_FILE="$idx" git -C "$repo" update-index --add --cacheinfo "100644,$blob,$name" || continue
                n=$((n + 1))
            done
            [ "$n" -gt 0 ] || exit 0
            tree=$(GIT_INDEX_FILE="$idx" git -C "$repo" write-tree) || exit 0
            commit=$(git -C "$repo" commit-tree "$tree" -p "$base" -m "run log: $n record(s) from ${CLUSTER:-?} [$(date '+%Y-%m-%d %H:%M')]") || exit 0
            if ${to[@]+"${to[@]}"} git -C "$repo" push --quiet origin "${commit}:refs/heads/${branch}" 2>/dev/null; then
                echo "run records: published $n to GitHub"
                exit 0
            fi
            sleep $((attempt * 3))          # someone pushed in between (or no credentials): fetch again and retry
        done
        echo "run records: push failed, not published now (they stay in $store; lsync on this cluster publishes them)"
    ) 9>"$repo/.git/.node_sync.lock"
    return 0
}
