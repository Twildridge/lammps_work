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

This installs, on whichever machine it runs on, the `lsync` shell command, `/lsync` for Claude Code, and the `nbstripout` filter so notebooks are committed without their embedded plot images. Running it again does no harm, for instance after moving the clone. [getting_started.md](getting_started.md) is the first-clone walkthrough.

`nbstripout`: `.gitattributes` at the repo root marks every `*.ipynb` for filtering, but each machine has to register the filter in its own `.git/config`, which `git clone` does not carry over. `install_lsync.sh` does that. If `lsync` warns that a notebook is about to be committed with its outputs, the install step has not run on that machine.

### Using it

```bash
lsync                        # auto-generates commit message from changed files
lsync "my custom message"    # use a specific commit message
```

`lsync` does: `git pull --rebase` → `git add -A` → `git commit` → `git push`. If there are no changes, it exits cleanly.

### Clusters → GitHub

Scripts or `.batch` files edited directly on a cluster must be committed there with `lsync`, or they never reach GitHub. The batch files self-sync with `git pull --rebase --autostash` at job start, which only stashes an uncommitted cluster-side edit and reapplies it; nothing saves it anywhere else.

Clusters also pull automatically at the start of every job (built into `run_lammps.sh`; the archived Bridges-2 / Pod runners had it too). To pull manually at any time:
```bash
cd ~/Documents/lammps_work && git pull
```

### Handling merge conflicts

When the same file was edited on two machines without a sync in between, `git pull --rebase` stops with a conflict. To resolve it:

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

Conflicts are avoided by running `lsync` (or `git pull`) before editing on any machine.

---

### SSH key setup (if `git push` gives a login error)

Needed once per cluster, the first time a push fails with `Permission denied (publickey)` or `fatal: Authentication failed`.

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

Each cluster needs its own key added to GitHub. On a Mac the keychain handles this, so the steps are not needed there.

---

### What is and isn't tracked by git

`.gitignore` excludes the large and machine-specific files:

| Excluded | Why |
|----------|-----|
| `*.lammpstrj`, `*.restart*` | Trajectory/restart files can be hundreds of GB |
| `output_files/`, `traj_files/` | Per-run outputs — too large and numerous |
| `log.lammps`, `*.log` | Redundant with run directories |
| `final_*.data` | LAMMPS-written final configs |
| `slurm_*.out`, `*.o*.*`, `*.e*.*` | SLURM stdout/stderr logs |
| `.DS_Store`, `__pycache__/`, `.ipynb_checkpoints/` | OS and IDE clutter |

Everything in `simulations/`, `scripts/` and the docs is tracked. `.data` input files are not (too large); they are copied by hand with `scp` or `rsync`.

---
