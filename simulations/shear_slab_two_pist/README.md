# shear_slab_two_pist: file names (since 2026-11-14)

> AI-generated working note (Claude), like the rest of the documentation in this repository.

Two decks live here, and every file name says which gel it is for. Both systems have the two
external NPT pistons that hold the bath at `PRESS_TARGET`.

| gel | deck | batches | input file | notebook `DECK` |
|---|---|---|---|---|
| **periodic** in x and y (the default) | `shear_slab_two_pist_periodic.lmp` | `shear_slab_two_pist_periodic_single.batch`, `shear_slab_two_pist_periodic_sweep.batch` | `<stem>_two_pist_pm40.data`, the file the compression and permeation decks read | `"shear_slab_two_pist"` |
| **isolated** block with grip plates | `shear_slab_two_pist_isolated.lmp` | `shear_slab_two_pist_isolated_single.batch`, `shear_slab_two_pist_isolated_sweep.batch` | `<stem>_x_plates.data` from `lib/add_more_plates_to_gel.py` | `"shear_slab_two_pist_iso"` |

Defaults, as symlinks to the periodic files: `shear_slab_two_pist.lmp`,
`shear_slab_two_pist_single.batch`, `shear_slab_two_pist_sweep.batch`. So `run_lammps.sh
shear_slab_two_pist ...` without `DECK_FILE` runs the periodic deck, and every batch names its
deck with `DECK_FILE` anyway.

Old names: `shear_slab_two_pist.lmp` was the isolated deck (now `..._isolated.lmp`);
`shear_slab_two_pist_single.batch` / `_sweep.batch` were the isolated batches (now
`..._isolated_single.batch` / `..._isolated_sweep.batch`); `shear_slab_two_pist_periodic.batch`
is now `shear_slab_two_pist_periodic_single.batch`. `docs/running_simulations.md`, `README.md`
and the shear notebooks' text still use the old names in places.

Differences that matter when choosing: the periodic deck has no unload phase and no wall-clock
guard (`wall_timeout`); the isolated deck has both.
