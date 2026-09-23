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

---

### Per-atom stresses too high (or G too low) in a run with `bond_style hybrid` and OpenMP threads

**Symptom:** the polymer partial stresses from `compute stress/atom` come out wrong (the bulk gel pressure in a `shear_slab` run read 1.70–1.76 instead of 1.50) while the thermo pressure, temperature and dynamics look normal.

**Cause (found 2026-09-22):** with `-sf omp` and **more than one OpenMP thread per MPI task**, `bond_style hybrid` tallies the *per-atom* bond virial at 1/N_threads of its true value (only one thread's share survives the reduction). Forces and the global virial are correct, so only `stress/atom`-based outputs are affected — but in shear the network shear stress is mostly FENE bond virial, so G would be badly wrong. Plain `bond_style fene` (the triaxial and `compress_slab` decks) is fine with `/omp`. Verified on the plated 14000002 file: bulk P 1.70 with 4 threads, 1.50 with 1 thread or without the suffix.

**Where it bites:** local smoke tests with `lmp -sf omp -pk omp N` (N > 1), or any batch with `--cpus-per-task > 1`. The Expanse `shear_slab.batch` uses 1 thread per task, so the 2026-04 shear runs were unaffected.

**Fix:** `shear_slab.lmp` now wraps the bond style in `suffix off` / `bond_style hybrid fene harmonic` / `suffix on`, so the bond styles run un-suffixed whatever the thread count (the pair style still runs as `lj/cut/omp`). Do the same in any new deck that needs `bond_style hybrid`.

