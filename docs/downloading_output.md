# Downloading Output From Clusters

[← back to README](../README.md)

After a job finishes, the output `.dat` files are in the working directory on the cluster. Trajectory files (`.lammpstrj`) are stored separately in scratch storage (the `traj_files/` symlink points there).

**Download output data files to MacBook** (into `flow_data_local/`):
```bash
# From MacBook terminal — adjust path to match your run
rsync -avP <username>@login.expanse.sdsc.edu:\
  ~/Documents/lammps_runs/triaxial_compression_<dataname>_<interaction>_<timestamp>/output_files/ \
  ~/Documents/lammps/flow_data_local/
```

**Download a trajectory file** (large — use the data transfer node):
```bash
rsync -avPz <username>@data.bridges2.psc.edu:\
  /ocean/projects/chm250028p/dpollard/lammps_trajectories/<run_folder>/<file>.lammpstrj.gz \
  ~/Downloads/
# Drop the 'z' flag if the file is already gzip-compressed
```

> **Trajectory files are large** (often 1–50 GB). Only download if you need to analyze atom positions directly (e.g. volume fraction profiles). The `.dat` output files are much smaller and contain pre-computed quantities.

---
