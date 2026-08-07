# Cluster Setup (first time on Expanse / Bridges-2 / Pod)

[← back to README](../README.md)

> **What is a cluster?** A cluster is a collection of powerful computers ("nodes") you access remotely. You write a script describing your job (how many CPUs, how long), submit it to a queue manager called SLURM, and SLURM runs it when resources are available. You don't interact with the job while it runs.

SSH into the cluster, then clone the repo:

```bash
# Example for Expanse:
ssh <username>@login.expanse.sdsc.edu

cd ~/Documents
mkdir -p lammps_data/input_data    # for .data files you'll scp over
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
chmod +x scripts/*.sh
```

**Activate the Python environment** (required for post-processing scripts):

The run scripts load Python differently depending on the cluster. No extra steps are needed on Expanse or Bridges-2 — `module load anaconda3/...` is called automatically inside `run_lammps.sh`. On **Pod**, you need to activate a conda environment once before submitting jobs:

```bash
# Pod only — one-time setup (and any time you open a new terminal):
module load miniconda
conda activate lammps_analysis
```

If `lammps_analysis` doesn't exist yet on Pod, create it:
```bash
module load miniconda
conda create -n lammps_analysis python=3.11 numpy scipy matplotlib pandas -y
conda activate lammps_analysis
```

**Configure git credentials** (so `git pull` works without a password prompt):
```bash
git config --global user.name "Twildridge"
git config --global user.email "semiinfiniteslab@icloud.com"
git config --global credential.helper store
# Then do: git pull  — it'll ask for your GitHub PAT once, then store it
```

> A **Personal Access Token (PAT)** replaces your GitHub password for command-line use. Generate one at github.com → Settings → Developer settings → Personal access tokens → Tokens (classic). Check the `repo` scope, copy it, and paste it when prompted by `git pull`.

**Copy your .data files to the cluster:**
```bash
# From MacBook terminal:
scp ~/Documents/lammps_data/input_data/your_file.data \
    <username>@data.expanse.sdsc.edu:~/Documents/lammps_data/input_data/
```

Full details for each cluster are in [`expanse_lammps_guide.md`](../expanse_lammps_guide.md) and [`slurm_commands_and_compiling.md`](../slurm_commands_and_compiling.md).
