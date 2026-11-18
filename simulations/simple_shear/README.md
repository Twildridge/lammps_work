# simple_shear: file names (since 2026-11-14)

> AI-generated working note (Claude), like the rest of the documentation in this repository.

Two decks live here, and every file name says which gel it is for. Both systems have the two
external NPT pistons that hold the bath at `PRESS_TARGET`.

| gel | deck | batches | input file | notebook `DECK` |
|---|---|---|---|---|
| **periodic** in x and y (the default) | `simple_shear_periodic.lmp` | `simple_shear_periodic_single.batch`, `simple_shear_periodic_sweep.batch` | `<stem>_two_pist_pm40.data`, the file the compression and permeation decks read | `"simple_shear"` |
| **isolated** block with grip plates | `simple_shear_isolated.lmp` | `simple_shear_isolated_single.batch`, `simple_shear_isolated_sweep.batch` | `<stem>_x_plates.data` from `lib/add_more_plates_to_gel.py` | `"simple_shear_iso"` |

Defaults, as symlinks to the periodic files: `simple_shear.lmp`,
`simple_shear_single.batch`, `simple_shear_sweep.batch`. So `run_lammps.sh
simple_shear ...` without `DECK_FILE` runs the periodic deck, and every batch names its
deck with `DECK_FILE` anyway.

Old names. Before 2026-11-14 the unsuffixed deck and batches were the isolated ones (now
`..._isolated.lmp`, `..._isolated_single.batch`, `..._isolated_sweep.batch`) and the one
periodic batch was `..._periodic.batch`; `docs/running_simulations.md`, `README.md` and the
shear notebooks' text still describe that layout in places. Before 2026-11-18 this folder and
every file in it were named `shear_slab_two_pist...`, and runs made until then are in
`lammps_runs/shear_slab_two_pist/` on the cluster (the notebooks' sync looks there too).

Differences that matter when choosing: the periodic deck has no unload phase and no wall-clock
guard (`wall_timeout`); the isolated deck has both.
