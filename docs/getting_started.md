# Getting Started

[← back to README](../README.md)

This page goes from a fresh clone of the repo to a first finished simulation and its plots. It assumes no git experience. The commands can be pasted as they are.

If something here does not work, ask Dylan rather than debugging alone. If the page itself is wrong, say so, and it gets fixed.

---

## Two things you need first

1. A GitHub account, added as a collaborator on `Twildridge/lammps_work`. Dylan does this from the repo's Settings → Collaborators page.
2. A cluster account on Expanse (Bridges-2 and Pod are no longer used). Dylan or your PI requests it on your behalf. It can take a few days.

The rest of the page assumes both are in place.

---

## Git in two sentences

Git tracks every change to every file in this project, and GitHub hosts the shared copy that every machine syncs with. The only command you need is `lsync`, which pulls the latest changes, commits yours and pushes them in one step; `git commit` and `git push` are never typed by hand here.

---

## Step 1 — One-time setup on your own computer

Open a terminal (macOS: the Terminal app; Windows: WSL or Git Bash).

Clone the repository, which downloads it to your computer:
```bash
cd ~/Documents
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
chmod +x scripts/*.sh
```

The first clone asks for a GitHub username and password. GitHub does not accept the account password here; it wants a Personal Access Token. Make one at [github.com/settings/tokens](https://github.com/settings/tokens) under "Tokens (classic)" with the `repo` scope checked, paste it at the prompt, and keep it in a password manager because it is needed again.

Install the `lsync` command ([github_sync.md](github_sync.md) describes what it runs):
```bash
bash scripts/install_lsync.sh
source ~/.zshrc      # or: source ~/.bashrc, if that's what you use
```

Test it:
```bash
lsync
```
It should print `Nothing to commit. Already up to date.`

Make a folder for simulation input files. These are not tracked by git ([github_sync.md](github_sync.md) says why):
```bash
mkdir -p ~/Documents/lammps_data/input_data
```

---

## Step 2 — One-time setup on the cluster

SSH into Expanse ([cluster_reference.md](cluster_reference.md) has the login addresses):
```bash
ssh <your-username>@login.expanse.sdsc.edu
```

Repeat the clone and `lsync` setup there, with `--no-claude` because the cluster has no Claude Code:
```bash
cd ~/Documents
mkdir -p lammps_data/input_data
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
chmod +x scripts/*.sh
bash scripts/install_lsync.sh --no-claude
source ~/.bashrc
```

The first push from the cluster usually needs git credentials and possibly an SSH key. [github_sync.md](github_sync.md) has the steps and the error messages that mean each one is needed.

That is all the one-time setup. The steps below are repeated for every simulation.

---

## Step 3 — Run your first simulation, start to finish

The example is the smallest simulation in the repo, a box of pure solvent with no gel. It is not interesting on its own, but it goes through every stage (build a data file, copy it to the cluster, submit a job, download the results, plot them) in about 15 minutes of your time and an hour of cluster time. Every other simulation type ([running_simulations.md](running_simulations.md)) follows the same steps.

### 3a. Build the input data file (on your own computer)

Every simulation starts from a `.data` file that lists every particle's position. These are built with the Jupyter notebooks in `scripts/`.

Start Jupyter:
```bash
cd ~/Documents/lammps_work/scripts
jupyter lab
```
This opens a browser tab. Open `pure_solvent_1.ipynb` and run every cell in order (Shift+Enter on each, or Run → Run All Cells).

The last cell writes the input file to `lammps_data_files_local/solvent_box_rho04_02.data`, in the data folder next to `lammps_work`. The name encodes the parameters (density 0.4, box variant 02). Nothing needs changing for this first run.

### 3b. Get the file onto the cluster

`run_lammps.sh` reads input files from one folder on the cluster, `~/Documents/lammps_data/input_data/`. Upload the file there:

```bash
# On your own computer:
scp ~/Documents/lammps_data_files_local/solvent_box_rho04_02.data \
    <your-username>@login.expanse.sdsc.edu:~/Documents/lammps_data/input_data/
```

### 3c. Submit the job

SSH back into the cluster if you're not still connected, then:
```bash
cd ~/Documents/lammps_work
git pull                      # make sure you have the latest scripts
cd simulations/solvent_pure
```

Open `solvent_pure_single.batch` in a text editor (`nano solvent_pure_single.batch` over SSH) and find this block near the bottom:
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
Leave the rest as it is. `INTERACTION="1.0_0.0"` is right for a solvent-only run: the second number is the polymer–solvent coupling, and there is no polymer. Save and exit.

Submit it:
```bash
sbatch solvent_pure_single.batch
```
This prints a job ID, e.g. `Submitted batch job 12345678`. The job is queued; the terminal can be closed.

### 3d. Check on it

```bash
squeue -u $USER
```
The state column shows `PD` (pending) or `R` (running). This is a 1-node, 1-hour job, so it usually starts within minutes and finishes well inside the hour.

While it runs, the LAMMPS log can be watched live:
```bash
tail -f ~/Documents/lammps_runs/solvent_pure/solvent_pure_solvent_box_rho04_02_1.0_0.0_<timestamp>/log.lammps
```
The exact folder name is on the `Working directory:` line near the top of the job's `.o<jobid>.*` SLURM output file. Ctrl-C stops the watching, not the job.

When the job finishes, `squeue -u $USER` will no longer list it, and the log will end with a line like `Total wall time: ...`.

### 3e. Get your results back

From your own computer:
```bash
rsync -avP <your-username>@login.expanse.sdsc.edu:\
  ~/Documents/lammps_runs/solvent_pure/solvent_pure_solvent_box_rho04_02_1.0_0.0_<timestamp>/output_files/ \
  ~/Documents/flow_data_local/solvent_pure/      # the local output folder next to lammps_work (the notebooks read ../../flow_data_local)
```
[downloading_output.md](downloading_output.md) covers what is worth downloading and why trajectories are handled separately.

### 3f. Look at the results

The plots are generated at the end of every job. To regenerate them by hand:
```bash
cd ~/Documents/lammps_runs/solvent_pure/solvent_pure_solvent_box_rho04_02_1.0_0.0_<timestamp>
python ~/Documents/lammps_work/scripts/lib/plot_lammps_log.py . solvent_box_rho04_02_1.0_0.0_100000
```
This writes convergence plots (temperature, pressure and volume against time) to `./output_plots/`. Temperature and pressure should be flat by the end of the run. [analysis.md](analysis.md) describes every analysis script and notebook.

Every other simulation type in the repo, including the hydrogel runs the project is about, uses the same sequence: build a data file, copy it over, edit a `.batch` file, submit, download, analyse.

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

