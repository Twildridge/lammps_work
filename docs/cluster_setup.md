# Cluster Setup (first time on Expanse / Bridges-2 / Pod)

[← back to README](../README.md)

A cluster is a set of computers (nodes) used remotely. A job is a script that says how many CPUs it needs and for how long; the queue manager, SLURM, runs it when the nodes are free. Nobody interacts with a job while it runs.

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

**Python for the post-processing scripts.** On Expanse and Bridges-2 nothing needs doing: `run_lammps.sh` runs `module load anaconda3/...` itself. On Pod a conda environment has to be activated before submitting jobs:

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

**Git credentials**, so that `git pull` does not prompt for a password every time:
```bash
git config --global user.name "Twildridge"
git config --global user.email "semiinfiniteslab@icloud.com"
git config --global credential.helper store
# Then do: git pull  — it'll ask for your GitHub PAT once, then store it
```

A Personal Access Token (PAT) stands in for the GitHub password on the command line. Generate one at github.com → Settings → Developer settings → Personal access tokens → Tokens (classic) with the `repo` scope checked, and paste it when `git pull` asks.

**Copy your .data files to the cluster:**
```bash
# From MacBook terminal:
scp ~/Documents/lammps_data/input_data/your_file.data \
    <username>@data.expanse.sdsc.edu:~/Documents/lammps_data/input_data/
```

Cluster-specific detail is in [`expanse_lammps_guide.md`](../expanse_lammps_guide.md) and [`slurm_commands_and_compiling.md`](../slurm_commands_and_compiling.md).
