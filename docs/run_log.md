# Run log

The "Recent runs" table at the top of the README and the full table in [`runs/README.md`](../runs/README.md) list the jobs started on any cluster, newest first. Nothing has to be entered by hand; job IDs of old runs and notes on what went wrong can be added in `runs/annotations.json` (below).

## What is recorded

One JSON file per job, `runs/<cluster>_<jobid>.json`:

| Column in the table | Source |
|---|---|
| Started, cluster, job | `run_info.txt` of the run folder |
| Deck | the simulation folder (and the deck file when `DECK_FILE` picked another one); "continuation" for a `continue_sim.sh` job |
| Lattice | read from the data file name: `random`, `rot45`, else diamond; beads per strand from `<N>beads` |
| Initial data file | the file the job read (`-var dataname`, or the resume file of a continuation) |
| Final data file(s) | `final_*.data` in the run's `data_files/` when the job ended. The job's data name (the initial file without `.data`) is written `{dataname}` |
| Simulated | steps summed over every `run` in `log.lammps`, and steps × timestep in τ. While the job runs: the planned production steps |
| Wall time | "Total wall time" of the log, else end minus start |
| Status | `running`, `complete`, `partial (wall limit)`, `failed`, `stopped at wall budget, continued`, or `ended without finishing` (the job disappeared without writing its end, usually killed at the SLURM limit). A reason follows the status when there is one: the note in `runs/annotations.json`, else the last `ERROR` line of `log.lammps` for a job that did not complete |
| Laptop RUN_ID | the folder under `flow_data_local/` whose files were pulled from this run; see below |

A job appears when it starts, not when it is submitted: a pending job is not in the table.

## How a record reaches GitHub

`run_lammps.sh` and `continue_sim.sh` call `run_record_start` after `run_info.txt` is written and `run_record_end` when LAMMPS has ended (`scripts/lib/run_record.sh`; the record itself is built by `scripts/lib/run_record.py`). On the cluster the records are kept in `$RUNS_DIR/.run_records/`, outside the repository.

Publishing does not use the working tree. It fetches the branch, adds the changed records to that commit's tree in a private index, makes a commit with `git commit-tree` and pushes that commit. The checkout that running jobs read their scripts from is never modified, and uncommitted edits on the cluster do not get in the way. The record files appear in the working tree with the next pull.

If the push fails (a compute node without network or credentials), the records stay in `.run_records/` and the next job or the next `lsync` on that cluster publishes them. `lsync` on a cluster also adds records for run folders of the last 30 days that have none, and closes records of jobs that are no longer in the queue. A folder without a `run_info.txt` (started before 2026-11-17) is recorded from its name and its `log.lammps`: no job ID, and the data file name is read from the folder name.

Nothing in this path can fail a job: every step is best-effort and prints one line to the job's output. `RUN_RECORDS=0` in the environment turns it off.

## Job IDs and notes added by hand

`runs/annotations.json` is the only hand-written file of the run log. Each key is a run folder name (the `run_name` of a record) or a job ID:

```json
{
 "<run folder name or job ID>": {"job_id": "55745206", "note": "NOT usable. Plates seated in the fringe ..."}
}
```

`job_id` fills the Job column of a record that has none (the backfilled ones). `note` is printed after the status, for a run that failed or finished but is not usable. Both are read when the tables are built. They are not written into the records, because a cluster pushes its own copy of each record over the one in the repository.

## How the tables are rebuilt

`lsync` (`scripts/git_sync.sh`) runs `python3 scripts/lib/run_record.py table` after its pull. That rewrites the block between the `RECENT_RUNS` markers in `README.md` (five rows) and `runs/README.md` (all rows). The output depends only on the records, so a sync with no new record changes nothing.

Run `lsync` on the laptop to see what the clusters have recorded since the last one.

## The laptop RUN_ID

The RUN_ID is the name of a folder under `flow_data_local/raw_data/<deck>/`, the data folder next to `lammps_work` (`install.sh` creates it, with `lammps_data_files_local/`, on a personal computer; neither is in git). The link between a RUN_ID and the cluster run folders is a small file in the repository, `runs/local/<deck>__<RUN_ID>.json`, so the table on GitHub shows it for everyone. It is written in two ways:

- When a notebook's sync cell pulls files from a cluster, the staging step reports which run folders they came from and `sync_pull` (`scripts/lib/uniaxial.py`) writes the link. This one is exact.
- `lsync` on a computer that has `flow_data_local/` cross-references the folders already there with the recorded jobs (`cross_reference` in `scripts/lib/run_record.py`), without a cluster login and without rerunning a notebook. A folder named like a cluster run folder is that run. Otherwise a job matches when it is of the same deck and a file in the folder carries its tag `_<dataname>_<interaction>_<steps>`. When the same data file was run more than once, each file is given the newest of those jobs that had started when the file was downloaded, which is the job the sync would have staged. A link written by a pull is never replaced by this.

The second way is an inference from file names and download times. It can be wrong when two jobs of one deck used the same data file and the older one's files were downloaded after the newer one had started; rerun that notebook's sync cell with `FORCE_SYNC=True` to replace it with the exact link. It also assumes the laptop and the cluster keep the same clock (both Pacific time here).

A job with no RUN_ID has no folder on any laptop that ran `lsync`: its data was never downloaded, or the folder was superseded by a later run of the same data file.

## Limits

- Jobs started before 2026-11-20 are in the log only if `lsync` on the cluster backfilled them (the last 30 days at the time of that `lsync`). Backfilled rows from before 2026-11-17 have no job ID unless `runs/annotations.json` gives one.
- `lsync` on a cluster prints a `>>> run log (...)` line with what it found and whether it published; when a pull changes the sync scripts it restarts itself with the new version.
- Jobs that do not go through `run_lammps.sh` or `continue_sim.sh` (the benchmark, the data-file builders, the screens) are not recorded.
- A job killed at the SLURM limit keeps the status `running` until the next `lsync` on that cluster.
