# Syncing With GitHub

[← back to README](../README.md)

The repo lives at `https://github.com/Twildridge/lammps_work`. Changes flow like this:

```
MacBook (edit scripts/notebooks)
    ↓  lsync
GitHub
    ↓  git pull (automatic at start of each cluster job)
Expanse / Bridges-2 / Pod
```

### One-time setup (each machine: MacBook, Expanse, Bridges-2, Pod)

```bash
bash scripts/install_lsync.sh
source ~/.zshrc     # or ~/.bashrc on a cluster
```

This single command sets up everything `lsync` needs, on whichever machine you run it: the `lsync` shell command itself, `/lsync` in Claude Code if you use it, and `nbstripout` (see the note below) so notebooks never commit with megabytes of embedded plot images. It's idempotent — safe to re-run any time, e.g. after moving the clone. See [getting_started.md](getting_started.md) for the full first-clone walkthrough.

> **`nbstripout` (notebook output stripping):** the `.gitattributes` file at the repo root marks every `*.ipynb` for filtering, but each machine still has to *register* the filter locally (this is a `.git/config` setting — it doesn't travel with `git clone`). `install_lsync.sh` does this for you. If `lsync` ever warns that a notebook is about to commit with its outputs embedded, that's this step not having run on that machine — just run the command above.

### Using it

```bash
lsync                        # auto-generates commit message from changed files
lsync "my custom message"    # use a specific commit message
```

`lsync` does: `git pull --rebase` → `git add -A` → `git commit` → `git push`. If there are no changes, it exits cleanly.

### Clusters → GitHub

If you edit scripts or `.batch` files directly on a cluster, run `lsync` there too (after the one-time setup above) — otherwise those edits stay local to the cluster and never reach GitHub. This matters more than it sounds: `triaxial_compression.batch` self-syncs with `git pull --rebase --autostash` at job start, so an uncommitted cluster-side edit is only ever *stashed and silently reapplied*, never actually saved anywhere else. Running `lsync` on the cluster commits it for real.

Clusters also pull automatically at the start of every job (built into `run_lammps.sh` and `run_lammps_bridges.sh`). To pull manually at any time:
```bash
cd ~/Documents/lammps_work && git pull
```

### Handling merge conflicts

If you edited the same file on two machines without syncing between them, `git pull --rebase` will stop with a conflict error. Fix it like this:

```bash
git status                  # identify the conflicting file(s)
# Open each conflicting file — look for markers like:
#   <<<<<<< HEAD
#   (your local version)
#   =======
#   (incoming version from GitHub)
#   >>>>>>> origin/main
# Edit the file to keep what you want, removing all the marker lines.
git add <filename>
git rebase --continue
git push
```

The simplest way to avoid conflicts: always run `lsync` (or `git pull`) before you start editing on any machine.

---

### SSH key setup (if `git push` gives a login error)

This is needed once per cluster the first time you try to push. The error usually looks like `Permission denied (publickey)` or `fatal: Authentication failed`.

**1. Generate an SSH key on the cluster:**
```bash
ssh-keygen -t ed25519 -C "pollard@ucsb.edu"
# Press Enter for all prompts (default location, no passphrase needed)
```

**2. Print your public key:**
```bash
cat ~/.ssh/id_ed25519.pub
# Copy the entire line of output
```

**3. Add the key to GitHub:**
- Go to [github.com/settings/ssh/new](https://github.com/settings/ssh/new)
- Title: something descriptive like `"Bridges-2 cluster"` or `"Expanse cluster"`
- Key: paste the full line from step 2
- Click **Add SSH key**

**4. Test the connection:**
```bash
ssh -T git@github.com
# Expected: "Hi Twildridge! You've successfully authenticated..."
```

**5. Switch your repo remote to SSH:**
```bash
git remote set-url origin git@github.com:Twildridge/lammps_work.git
```

**6. Try pushing again:**
```bash
git push
```

Repeat steps 1–6 for each cluster (Expanse, Bridges-2, Pod) — each needs its own key added to GitHub. On your MacBook, this is usually handled automatically by the macOS keychain, so you shouldn't need to do it there.

---

### What is and isn't tracked by git

The `.gitignore` excludes large files so the repo stays fast:

| Excluded | Why |
|----------|-----|
| `*.lammpstrj`, `*.restart*` | Trajectory/restart files can be hundreds of GB |
| `output_files/`, `traj_files/` | Per-run outputs — too large and numerous |
| `log.lammps`, `*.log` | Redundant with run directories |
| `final_*.data` | LAMMPS-written final configs |
| `slurm_*.out`, `*.o*.*`, `*.e*.*` | SLURM stdout/stderr logs |
| `.DS_Store`, `__pycache__/`, `.ipynb_checkpoints/` | OS and IDE clutter |

Everything in `simulations/`, `scripts/`, and the docs is tracked. `.data` input files are **not** tracked (too large) — manage them manually with `scp`/`rsync`.

---
