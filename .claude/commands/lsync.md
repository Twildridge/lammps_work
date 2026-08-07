---
description: Sync lammps_work with GitHub — pull --rebase, stage all, commit, push
argument-hint: [commit message]
allowed-tools: Bash(sh -c:*)
---

!`sh -c 'for d in "$LAMMPS_WORK_DIR" "$(git rev-parse --show-toplevel 2>/dev/null)" ./lammps_work ~/Documents/lammps_work ~/docs/grad_research/lammps/lammps_work; do if [ -n "$d" ] && [ -f "$d/scripts/git_sync.sh" ]; then exec bash "$d/scripts/git_sync.sh" "$ARGUMENTS"; fi; done; echo "ERROR: could not find lammps_work/scripts/git_sync.sh from $(pwd)."; echo "Fix: export LAMMPS_WORK_DIR=/path/to/lammps_work"; exit 1'`

Above is the output of `git_sync.sh` — the same script behind the `lsync` shell
alias, so this command and the terminal do exactly the same thing.

Summarise in one or two lines: what was committed and pushed, or that there was
nothing to commit.

If the pull hit a rebase conflict, name the conflicting files and stop — report
it, do not resolve it. This repo is shared and edits are often made directly on
Expanse without being committed, so a conflict normally means the same file was
changed in two places and a human has to pick the winner.

If the script was not found, tell the user to `export LAMMPS_WORK_DIR` to their
clone path.
