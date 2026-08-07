# Getting Started

[← back to README](../README.md)

This walks you through everything from "I have never used this repo" to "I ran a simulation and looked at its output." It assumes no prior git experience. Every command below is copy-pasteable.

If you get stuck anywhere, the goal is to unblock you, not to make you debug alone — ask.

---

## Before you start: two things someone needs to give you

1. **A GitHub account**, added as a collaborator on `Twildridge/lammps_work` (ask Dylan — this has to be done from his GitHub account, on the repo's Settings → Collaborators page).
2. **A cluster account** (Expanse, Bridges-2, and/or Pod, depending on what you'll run). This is a separate request Dylan or your PI has to make on your behalf — it can take a few days, so ask early.

Everything after this point assumes you already have both.

---

## Two sentences on git, so the rest of this makes sense

**Git** tracks every change to every file in this project, and **GitHub** is the website that hosts the shared copy everyone syncs with. You'll only ever need one command — `lsync` — which pulls the latest changes, saves yours, and uploads them, all in one step. You do not need to learn `git commit`, `git push`, or any of the usual git vocabulary to use this repo.

---

## Step 1 — One-time setup on your own computer

Open a terminal (macOS: the **Terminal** app; Windows: use WSL or Git Bash — ask if you're not sure which).

**Clone the repository** (this downloads it to your computer):
```bash
cd ~/Documents
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
chmod +x scripts/*.sh
```

The first time you do this, it will ask for your GitHub username and password — GitHub no longer accepts your real password here. Use a **Personal Access Token** instead: go to [github.com/settings/tokens](https://github.com/settings/tokens) → "Tokens (classic)" → generate one with the `repo` scope checked, and paste that in when prompted. Save it somewhere (a password manager), you'll need it again.

**Install the `lsync` command** (this is the one command you'll actually use — see [github_sync.md](github_sync.md) for what it does under the hood):
```bash
bash scripts/install_lsync.sh
source ~/.zshrc      # or: source ~/.bashrc, if that's what you use
```

Test it:
```bash
lsync
```
It should print `Nothing to commit. Already up to date.` — if you see that, setup worked.

**Make a place for simulation input files** (not tracked by git — see why in [github_sync.md](github_sync.md)):
```bash
mkdir -p ~/Documents/lammps_data/input_data
```

---

## Step 2 — One-time setup on the cluster

SSH into the cluster (example below is Expanse — see [cluster_reference.md](cluster_reference.md) for Bridges-2 / Pod login addresses):
```bash
ssh <your-username>@login.expanse.sdsc.edu
```

Once you're in, repeat the same clone + `lsync` setup, this time with `--no-claude` (there's no Claude Code on a cluster):
```bash
cd ~/Documents
mkdir -p lammps_data/input_data
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
chmod +x scripts/*.sh
bash scripts/install_lsync.sh --no-claude
source ~/.bashrc
```

You'll likely also need to set up git credentials and possibly an SSH key on the cluster the first time you push from there — full instructions (with the exact error messages to watch for) are in [github_sync.md](github_sync.md).

That's every one-time step. **Everything below, you'll do again for every new simulation.**

---

## Step 3 — Run your first simulation, start to finish

We're going to run the *simplest* simulation in this repo — a small box of pure solvent, no polymer gel at all. It is not scientifically interesting on its own, but it exercises every stage of the pipeline (build a data file → copy it to the cluster → submit a job → download results → look at them) in about 15 minutes of your time and 1 hour of cluster time. Once this works, every other simulation type ([running_simulations.md](running_simulations.md)) follows the exact same shape.

### 3a. Build the input data file (on your own computer)

Every simulation starts from a `.data` file describing where every particle is. These are built with Jupyter notebooks in `scripts/`.

Start Jupyter:
```bash
cd ~/Documents/lammps_work/scripts
jupyter lab
```
This opens a browser tab. Open **`pure_solvent_1.ipynb`**, and run every cell in order (Shift+Enter on each, top to bottom — or use the menu: Run → Run All Cells).

The last cell writes a file to `../../lammps_data/solvent_pure/solvent_box_rho04_02.data` — that's your input file. The name encodes the parameters used (density 0.4, box variant "02"); you don't need to change anything for this first run.

### 3b. Get the file onto the cluster

Copy it into the one folder `run_lammps.sh` always looks in, `~/Documents/lammps_data/input_data/`, both locally and on the cluster:

```bash
# On your own computer:
cp ~/Documents/lammps_work/lammps_data/solvent_pure/solvent_box_rho04_02.data \
   ~/Documents/lammps_data/input_data/

# Then upload it to the cluster:
scp ~/Documents/lammps_data/input_data/solvent_box_rho04_02.data \
    <your-username>@login.expanse.sdsc.edu:~/Documents/lammps_data/input_data/
```

### 3c. Submit the job

SSH back into the cluster if you're not still connected, then:
```bash
cd ~/Documents/lammps_work
git pull                      # make sure you have the latest scripts
cd simulations/solvent_pure
```

Open `solvent_pure.batch` in a text editor (`nano solvent_pure.batch` works fine over SSH) and find this block near the bottom:
```bash
FOLDER="solvent_pure"
DATANAME="final_config_slab_support_5beads_tall_rho04_p1.5_1.0_1.0_600000_solvent_only"
INTERACTION="1.0_0.0"
NSTEPS=100000
```
Change `DATANAME` to the name of the file you just built, **without** the `.data` extension:
```bash
DATANAME="solvent_box_rho04_02"
```
Leave everything else as-is (`INTERACTION="1.0_0.0"` is correct for a solvent-only run — the second number, polymer–solvent coupling, has nothing to act on here). Save and exit.

Submit it:
```bash
sbatch solvent_pure.batch
```
This prints a job ID, e.g. `Submitted batch job 12345678`. Your job is now queued — you do not need to keep your terminal open or stay connected.

### 3d. Check on it

```bash
squeue -u $USER
```
Shows `PD` (pending, waiting for resources) or `R` (running). This is a small 1-node, 1-hour job, so it typically starts quickly and finishes well inside that hour.

Once it's running, you can watch the LAMMPS log live:
```bash
tail -f ~/Documents/lammps_runs/solvent_pure/solvent_pure_solvent_box_rho04_02_1.0_0.0_<timestamp>/log.lammps
```
(The exact folder name is printed near the top of the job's `.o<jobid>.*` SLURM output file if you're not sure of it — `Working directory: ...`.) Press Ctrl-C to stop watching (this does not stop the job).

When the job finishes, `squeue -u $USER` will no longer list it, and the log will end with a line like `Total wall time: ...`.

### 3e. Get your results back

From your own computer:
```bash
rsync -avP <your-username>@login.expanse.sdsc.edu:\
  ~/Documents/lammps_runs/solvent_pure/solvent_pure_solvent_box_rho04_02_1.0_0.0_<timestamp>/output_files/ \
  ~/Documents/lammps_work/flow_data_local/solvent_pure/
```
(Full detail on what's downloadable and why trajectories are handled separately: [downloading_output.md](downloading_output.md).)

### 3f. Look at what you ran

Plots are actually generated automatically at the end of every job (no extra step needed) — but here's how to regenerate or inspect them yourself:
```bash
cd ~/Documents/lammps_runs/solvent_pure/solvent_pure_solvent_box_rho04_02_1.0_0.0_<timestamp>
python ~/Documents/lammps_work/scripts/plot_lammps_log.py . solvent_box_rho04_02_1.0_0.0_100000
```
This writes convergence plots (temperature, pressure, volume vs. time) to `./output_plots/`. Open them and confirm temperature and pressure look flat/stable by the end — that's what "the simulation worked" looks like. Full details on every analysis script and notebook: [analysis.md](analysis.md).

**You just ran a full simulation, start to finish.** Every other simulation type in this repo — including the actual hydrogel physics this project is about — uses this same pattern: build a data file, copy it over, edit a `.batch` file, submit, download, analyze.

---

## What to read next

| I want to... | Read |
|---|---|
| Understand what each simulation type measures, and which one to use | [running_simulations.md](running_simulations.md) |
| See the full repo file tree and what lives where | [repository_layout.md](repository_layout.md) |
| Set up git credentials / SSH keys on a cluster, or fix a merge conflict | [github_sync.md](github_sync.md) |
| Understand `.data` files and which notebook builds which geometry | [building_data_files.md](building_data_files.md) |
| Run an analysis notebook or script on real output | [analysis.md](analysis.md) |
| Look up cluster login nodes, partitions, walltimes | [cluster_reference.md](cluster_reference.md) |
| Convert LJ units to real physical units (nm, MPa, ps) | [physics_reference.md](physics_reference.md) |
| Fix a specific error message | [common_issues.md](common_issues.md) |

If something in this repo doesn't work the way this guide says, that's a bug in the guide — say so.
