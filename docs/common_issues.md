# Common Issues

[← back to README](../README.md)

**"Data file not found"**
The `.data` file must be in `~/Documents/lammps_data/input_data/` with the exact name matching `DATANAME` in the batch file (no extra suffixes).

**"Disk quota exceeded" on Bridges-2**
Trajectory files ended up in your home directory. Verify the symlink: `ls -la ~/Documents/lammps_runs/*/traj_files` should show `→ /ocean/projects/...`. If not, move the `.lammpstrj` files to scratch manually.

**Job stuck in `CG` state on Bridges-2**
This is a cluster-side issue. Contact PSC support.

**Simulation diverged ("gel evaporated")**
Check that `epsSP` and `epsSS` are physically reasonable. Verify you're reading the correct data file. Ensure `timestep = 0.005` (larger timesteps can cause FENE bond divergence).

**`Invalid syntax in variable formula` on first run of `shear_slab.lmp`**
This was caused by `xztilt`/`yztilt` not being valid LAMMPS thermo keywords — the correct keywords are `xz` and `yz`. Already fixed in the current script.

**`fix ave/chunk norm none` error**
`norm none` requires LAMMPS ≥ March 2020. On older builds, replace the `prof_z_polymer/solvent` fixes with an equivalent `compute reduce` + `fix ave/time` approach (see git history for `slab_with_flow.lmp` pre-2026-08-06 removal for a worked example).

**`git pull` fails on cluster ("merge conflict")**
Run `git stash` to set aside local changes, then `git pull`, then `git stash pop` to restore them. If conflicts persist, resolve them manually or ask for help.

---
