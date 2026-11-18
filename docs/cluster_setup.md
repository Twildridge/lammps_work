# Cluster Setup (first time on Expanse or Pod)

[← back to README](../README.md)

A cluster is a set of computers (nodes) used remotely. A job is a script that says how many CPUs it needs and for how long; the queue manager, SLURM, runs it when the nodes are free. Nobody interacts with a job while it runs.

SSH into the cluster, clone the repo into `~/Documents` (the batch files expect it at `~/Documents/lammps_work`) and run the one setup command. It is the same two commands on Expanse and on Pod:

```bash
ssh <username>@login.expanse.sdsc.edu        # Pod: ssh <username>@pod-login1.cnsi.ucsb.edu, on the campus VPN

mkdir -p ~/Documents && cd ~/Documents
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
bash install.sh        # everything below; then log out and back in
```

`install.sh` (`scripts/setup.sh`) works out which cluster it is on from `scripts/lib/site.sh` and does, in order:

1. the parts every machine gets: the `lsync` and `lsbatch` commands, the `nbstripout` filter so notebook plots are not committed, and the input folder `~/Documents/lammps_data/input_data`;
2. what the jobs need on this cluster. **Expanse:** nothing to install; the jobs load `anaconda3` themselves and `site.sh` points at a shared LAMMPS binary, which is checked. **Pod:** the conda environment `lammps_analysis` (python 3.11, numpy, scipy, matplotlib, pandas) that the jobs activate for the plotters, created if missing; and the LAMMPS binary, which Pod does not provide, built on the login node with every thread if missing (`lammps_build/build_lammps_pod.sh`, about 10 min; `bash install.sh --no-build` postpones it);
3. `scripts/cluster_check.sh`, which compares `site.sh` with the machine (binary and its packages, modules, python, scratch, cores per node, Slurm partitions) and prints `OK` or `PROBLEM` per item. Its last line says whether the setup is complete.

Running `bash install.sh` again at any time is safe: it skips what exists and repeats the check. `bash scripts/cluster_check.sh` alone repeats only the check, for instance after the cluster's software changed.

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
