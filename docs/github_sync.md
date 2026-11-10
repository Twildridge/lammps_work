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

Scripts or `.batch` files edited directly on a cluster must be committed there with `lsync`, or they never reach GitHub. At job start `batch_prologue.sh` pulls fast-forward only (`node_git_sync.sh`) and re-runs the fresh batch; an uncommitted cluster-side edit blocks that pull, so the job then runs with the on-disk checkout. Nothing saves such an edit anywhere else.

Clusters also pull automatically at the start of every job (`batch_prologue.sh`). To pull manually at any time:
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

### Step 0: get collaborator access (do this before anything else)

The repo belongs to Dylan's GitHub account (`Twildridge`). Anyone can *download* it, but only people he has added as **collaborators** can *upload* (push) changes. Without that access, every push fails, however the keys and passwords below are set up.

**If you don't have a GitHub account yet:** sign up for free at [github.com/signup](https://github.com/signup), then send Dylan your GitHub **username** (the name shown at the top right of github.com when you're logged in; it is not your email).

**Requesting access:** GitHub has no "request access" button for a personal repo, so ask Dylan directly (email, Slack or in person) and include your GitHub username.

**What Dylan does (repo owner only):**
1. Go to [github.com/Twildridge/lammps_work/settings/access](https://github.com/Twildridge/lammps_work/settings/access) (repo page → **Settings** → **Collaborators**).
2. Click **Add people**, type the student's GitHub username, and pick them from the list.
3. Click **Add to repository**.

**What you do next:** GitHub sends you an invitation email. Click **Accept invitation** in it, or open [github.com/Twildridge/lammps_work/invitations](https://github.com/Twildridge/lammps_work/invitations) while logged in. The invite expires after 7 days; if it does, ask Dylan to send a new one.

**Checking that it worked:** once you've accepted, the invitation link above says there are no pending invitations, and pushes no longer fail with `403` or `Permission to Twildridge/lammps_work.git denied`.

---

### Setting up pushing from a cluster (Expanse, Bridges-2, Pod), one time per cluster

**Why this is needed.** GitHub has to confirm who you are before it accepts a push. On a Mac this happens automatically. On a cluster nothing is set up, so the first `lsync` usually stops at the push step with one of these errors:

| Error you see | What it means |
|---|---|
| `fatal: could not read Username for 'https://github.com'` | Git wanted to ask for your GitHub login but couldn't. This always happens when `lsync` runs from inside Claude Code, since Claude Code can't show password prompts. |
| `Permission denied (publickey)` | The cluster has no SSH key, or GitHub hasn't been given it yet. |
| `fatal: Authentication failed` / `403` | GitHub doesn't recognise the login, or you aren't a collaborator yet (see Step 0). |

The fix is an **SSH key**. Think of it as a password file that stays on the cluster: you give GitHub the matching public half once, and from then on every push is accepted with no prompts. Unlike a token, it doesn't expire.

**Common questions**
- *Does the cluster's two-factor login (Duo) get in the way?* No. Two-factor only applies to logging **in to** the cluster. Pushing is the cluster connecting **out** to GitHub, and the key handles that.
- *Does it matter that the repo belongs to someone else?* No. The key goes on **your own** GitHub account and proves the push comes from you. What you're allowed to push depends on your collaborator access (Step 0). If you're ever removed as a collaborator, the key stops working for this repo by itself.
- *Where does the key NOT go?* Not in the repo's **Settings → Deploy keys** page. It goes in your personal account settings (step 3 below).

**1. Make the key on the cluster.** Log in to the cluster and run:
```bash
ssh-keygen -t ed25519 -N "" -C "<your-name>@<cluster-name>"
```
Press **Enter** if it asks where to save. `-N ""` means "no passphrase", so `lsync` never stops to ask for one. The text after `-C` is just a label, for example `jkapner@expanse`.

**2. Show the public half of the key:**
```bash
cat ~/.ssh/id_ed25519.pub
```
This prints one long line starting with `ssh-ed25519`. Copy the whole line. This half is safe to share; the other file, `~/.ssh/id_ed25519` (no `.pub`), is the private half, and you should never copy or share it.

**3. Give the public key to GitHub:**
- Log in to github.com **as yourself** and open [github.com/settings/ssh/new](https://github.com/settings/ssh/new) (profile picture → **Settings** → **SSH and GPG keys** → **New SSH key**).
- **Title:** something you'll recognise, such as `Expanse`.
- **Key:** paste the line from step 2.
- Click **Add SSH key**. GitHub may ask for your password or 2FA code first.

**4. Check that GitHub accepts the key:**
```bash
ssh -T git@github.com
```
If it asks "Are you sure you want to continue connecting?", type `yes`. Success looks like:
`Hi <your-github-username>! You've successfully authenticated, but GitHub does not provide shell access.`
The "no shell access" part is normal and means everything is working.

**5. Tell git to use the key for pushing:**
```bash
cd ~/Documents/lammps_work
git remote set-url --push origin git@github.com:Twildridge/lammps_work.git
git remote -v
```
`--push` switches only uploads to the key. Downloads stay on the plain web address, so the automatic `git pull` at the start of each job keeps working without a login. `git remote -v` should now show:
```
origin  https://github.com/Twildridge/lammps_work.git (fetch)
origin  git@github.com:Twildridge/lammps_work.git (push)
```

**6. Run `lsync`.** If a commit is waiting to go up from an earlier failed push, `lsync` (or `git push`) sends it now.

**If step 4 hangs or times out:** that cluster blocks the usual SSH port. GitHub also listens on port 443. Add this to `~/.ssh/config` on the cluster, create the file if it doesn't exist, then repeat step 4:
```
Host github.com
    Hostname ssh.github.com
    Port 443
    User git
```

Each cluster needs its own key, so repeat these steps on each one. On a Mac they are not needed.

**About the "Committer: … Your name and email address were configured automatically" message:** this is only a warning. The commit still works; it is just labelled with your cluster username (for example `jkapner@login01.expanse.sdsc.edu`). To use your real name instead, run once per cluster:
```bash
git config --global user.name "Your Name"
git config --global user.email "the-email-on-your-github-account"
```

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
