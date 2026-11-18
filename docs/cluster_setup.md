# Cluster Setup (first time on Expanse or Pod)

[← back to README](../README.md)

A cluster is a set of computers (nodes) used remotely. A job is a script that says how many CPUs it needs and for how long; the queue manager, SLURM, runs it when the nodes are free. Nobody interacts with a job while it runs.

SSH into the cluster, then clone the repo:

```bash
# Example for Expanse:
ssh <username>@login.expanse.sdsc.edu

mkdir -p ~/Documents && cd ~/Documents    # the batch files expect the repo at ~/Documents/lammps_work
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
bash install.sh        # lsync, lsbatch, the notebook filter and ~/Documents/lammps_data/input_data; then log out and back in
```

**Python for the post-processing scripts.** On Expanse nothing needs doing: the job loads `anaconda3` itself. On Pod the job activates the conda environment `lammps_analysis`; `install.sh` creates it (python 3.11, numpy, scipy, matplotlib, pandas) and verifies the imports, through `scripts/pod_setup.sh`, which can be re-run alone at any time.

**LAMMPS on Pod.** Pod has no site LAMMPS, so each account builds its own once, directly on the login node with every thread it has, since Pod's queues are often full for hours (`lammps_build/README.md`). `install.sh` says so when the binary is missing. The whole Pod setup is four commands:

```bash
git clone https://github.com/Twildridge/lammps_work.git ~/Documents/lammps_work && cd ~/Documents/lammps_work
bash install.sh                                   # lsync, lsbatch, nbstripout, input_data, the conda env
bash lammps_build/build_lammps_pod.sh             # on the login node, 20-40 min: ~/software/lammps/.../lmp, the LMP_BIN of site.sh
bash scripts/cluster_check.sh
```

**Check the cluster settings.** `bash scripts/cluster_check.sh` compares `scripts/lib/site.sh` with the machine (LAMMPS binary, modules, scratch, cores per node) and lists what does not match.

**Git credentials**, so that `git pull` does not prompt for a password every time:
```bash
git config --global user.name "<your GitHub username>"
git config --global user.email "<the email of your GitHub account>"
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

Partitions, allocation sizes and the settings of each cluster are in [cluster_reference.md](cluster_reference.md). Older cluster notes are in [`expanse_lammps_guide.md`](../expanse_lammps_guide.md) and [`slurm_commands_and_compiling.md`](../slurm_commands_and_compiling.md).
