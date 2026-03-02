# Finding LAMMPS Build Information on Expanse (SDSC)

This guide walks you through finding what LAMMPS version is installed on Expanse, what packages it was compiled with, and what dependencies it requires — without needing to run a simulation.

---

## Step 1: Log into Expanse

```bash
ssh <your_username>@login.expanse.sdsc.edu
```

---

## Step 2: Find Available LAMMPS Modules

```bash
module spider LAMMPS
# or if that returns nothing:
module avail lammps
```

Note the full module name of the version you want (e.g. `LAMMPS/20230802`).

---

## Step 3: Check Required Dependencies

```bash
module spider LAMMPS/<version>
```

**This is the most important step.** Expanse uses a hierarchical module system — certain prerequisite modules (compiler, MPI) must be loaded *before* LAMMPS becomes available. This command tells you exactly what to load first and in what order.

---

## Step 4: Load the Modules

Load prerequisites first, then LAMMPS:

```bash
module purge
module load gcc/10.2.0
module load openmpi/4.1.3
module load LAMMPS/20230802
```

The exact names come from Step 3.

---

## Step 5: Find the LAMMPS Binary Name

The binary may not be called `lmp_mpi` as on other clusters:

```bash
which lmp
ls $(dirname $(which lmp))
```

It will likely be just `lmp` or `lmp_kokkos`.

---

## Step 6: Find What Packages Were Compiled In

**Do NOT run `lmp -h` on the login node.** If it is a GPU/Kokkos build, it will hang indefinitely trying to initialize GPUs — login nodes have none.

Instead, read the CMakeCache from the build directory directly.

### 6a. Find the CMakeCache

```bash
find /cm/shared/apps/spack/ -name "CMakeCache.txt" -path "*lammps*" 2>/dev/null
```

This returns a path like:

```
/cm/shared/apps/spack/0.17.3/gpu/b/opt/spack/.../lammps-XXXXXXX/.spack/archived-files/spack-build-XXXXXXX/CMakeCache.txt
```

### 6b. List Installed Packages

```bash
grep "PKG_.*:BOOL=ON" /path/to/CMakeCache.txt | sed 's/PKG_//; s/:BOOL.*//'
```

This prints a clean list of every LAMMPS package compiled into the build.

---

## Step 7: Check GPU Hardware

```bash
sinfo -p gpu --format="%N %G %m"
```

Expanse GPU nodes have V100 (32GB) GPUs. The current LAMMPS build uses CUDA 11.2 and Kokkos 3.7, which target V100 architecture.

---

## Notes for Updating SLURM Scripts

Key differences when adapting batch scripts from Bridges-2:

| Item | Bridges-2 | Expanse |
|---|---|---|
| Partition (CPU) | `RM` | `compute` |
| Partition (GPU) | `GPU` | `gpu` |
| Cores per node | 128 | 128 |
| GPU type | V100 | V100 (32GB) |
| Account flag | optional | **required**: `--account=<allocation>` |
| LAMMPS binary | explicit path | just `lmp` after module load |

A minimal GPU batch script header for Expanse:

```bash
#!/bin/bash
#SBATCH --job-name=lammps_run
#SBATCH --account=<your_allocation>
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --gpus=4
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

module purge
module load gcc/10.2.0
module load openmpi/4.1.3
module load cuda/11.2.2
module load LAMMPS/20230802

mpirun -n $SLURM_NTASKS lmp \
    -k on g $SLURM_GPUS_ON_NODE \
    -sf kk -pk kokkos newton on neigh half \
    -var dataname $DATANAME \
    -in your_script.lmp
```

> The module versions above reflect what is currently installed on Expanse. Always re-run Steps 2–3 to confirm, as modules can be updated.

---

## Quick Reference

```bash
module spider LAMMPS                                # 1. find versions
module spider LAMMPS/<version>                      # 2. find dependencies
module purge && module load <deps> LAMMPS/<version> # 3. load modules
which lmp && ls $(dirname $(which lmp))             # 4. find binary
find /cm/shared/apps/spack/ -name "CMakeCache.txt" -path "*lammps*" 2>/dev/null  # 5a. find cache
grep "PKG_.*:BOOL=ON" /path/to/CMakeCache.txt | sed 's/PKG_//; s/:BOOL.*//'      # 5b. list packages
sinfo -p gpu --format="%N %G %m"                   # 6. check GPU hardware
```

---

## GitHub Setup and Basic Workflow

GitHub no longer accepts passwords for Git operations — you need a **Personal Access Token (PAT)** instead.

### Generating a Personal Access Token

1. Go to [github.com](https://github.com) and log in
2. Click your profile picture → **Settings**
3. Scroll down to **Developer settings** (bottom of left sidebar)
4. Click **Personal access tokens** → **Tokens (classic)**
5. Click **Generate new token** → **Generate new token (classic)**
6. Give it a descriptive name (e.g. `expanse-token`)
7. Set an expiration (90 days is reasonable)
8. Check the **repo** scope (this gives full access to your repositories)
9. Click **Generate token** at the bottom
10. **Copy the token immediately** — GitHub will never show it again

### Configuring Git on Expanse

Set your identity to match your GitHub account:

```bash
git config --global user.name "Twildridge"
git config --global user.email "semiinfiniteslab@icloud.com"
```

Save your credentials so you don't have to re-enter the token every time:

```bash
git config --global credential.helper store
```

### Cloning the Repository

```bash
cd ~/Documents   # or wherever you want the repo to live
git clone https://github.com/Twildridge/lammps_work.git
cd lammps_work
```

When prompted for a password, **paste your Personal Access Token** (not your GitHub password). With `credential.helper store` set, this is only needed once.

---

### Basic Git Workflow

The typical day-to-day cycle is: pull → make changes → add → commit → push.

**1. Before starting work, pull any changes from GitHub:**
```bash
git pull
```

**2. After making changes, check what has been modified:**
```bash
git status
```

**3. Stage the files you want to commit:**
```bash
git add filename.txt          # stage a specific file
git add .                     # stage all changed files
```

**4. Commit with a descriptive message:**
```bash
git commit -m "Brief description of what you changed"
```

**5. Push to GitHub:**
```bash
git push
```

---

### Other Useful Commands

```bash
git log --oneline             # see recent commit history
git diff filename.txt         # see what changed in a file before staging
git restore filename.txt      # discard uncommitted changes to a file
git branch                    # see which branch you're on
```

> **Important:** Always `git pull` before starting new work, especially if multiple people are using the same repository. This avoids merge conflicts.