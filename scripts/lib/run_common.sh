# Shared by run_lammps.sh and continue_sim.sh (sourced; expects lib/site.sh sourced first).
# A sourced file is read whole before it runs, so replacing it mid-job is harmless.

# ---- deck knobs --------------------------------------------------------------------------------
# LAMMPS variable names; the environment variable is the same name in upper case (PRESS_TARGET ->
# -var press_target). A knob is forwarded only when its environment variable is non-empty, so a
# deck keeps its own `variable X index` default otherwise. Decks that do not define a knob ignore it.
KNOBS=(
    press_target vel_seed comm_cutoff prerelaxed skip_widom calib_frames calib_dump_every
    baro_mode pin_gel piston_transparent                                     # slab_with_support
    rg_halt rg_block rg_every rg_tol rg_nconsec rg_extra                    # Rg-plateau halt
    mean_box mean_box_window mean_box_ramp mean_box_hold mean_lx mean_ly    # mean-box ending (slab_with_support, 2026-11-14)
    piston_mass c_pist_frac npt_piston_steps settle_halt                    # NPT pistons
    dp_piston n_ramp ramp_steps perm_stress_curves disp_nfreq               # permeation
    frozen frozen_steps frozen_at                                           # permeation, frozen network
    vctrl vctrl_steps vctrl_avg_steps vfeed                                 # permeation, velocity-controlled flux (Phase 2v)
    drive_split hold_auto seat_mode disp_fine_nfreq disp_fine_bin           # compression decks
    unload unload_steps unload_margin                                       # compression / shear unload
    save_strained                                                           # strained-state restart at the end of each hold (default 1; 0 = off)
    l0_avg_steps                                                            # window of the L0 / L_ss means (compression, compress_slab, permeation)
    fringe_max rg_strain_min prestrain_max                                  # credit guards of triaxial_compression_two_pist (0 = off)
    ref_avg_steps nsteps_shear_base disp_fine_span stress_freq halt_every   # shear
    strains compressions                                                    # sweep lists (see below)
    nsteps_eq nsteps_prod                                                   # phase sweeps
)
# Knobs the decks cannot default themselves: forwarded with these values when the environment
# does not set them. strains / compressions also go out as the quoted string <name>_list.
declare -A KNOB_DEFAULTS=([strains]=0.1 [compressions]=0.1 [nsteps_eq]=200000 [nsteps_prod]=100000)

# Append the -var pairs of every set knob to the array LAMMPS_VARS.
knob_vars() {
    local k env val
    for k in "${KNOBS[@]}"; do
        env="${k^^}"
        val="${!env:-${KNOB_DEFAULTS[$k]:-}}"
        [ -n "$val" ] || continue
        case "$k" in
            strains|compressions)
                # shellcheck disable=SC2206  # word-split on purpose: one -var with several values
                LAMMPS_VARS+=(-var "$k" $val -var "${k}_list" "$val") ;;
            *)  LAMMPS_VARS+=(-var "$k" "$val") ;;
        esac
    done
}

# Export every knob that is set, for a child job (continuation) or postprocess.sh.
export_knobs() {
    local k env
    for k in "${KNOBS[@]}" chain wall_margin deck_file; do
        env="${k^^}"
        [ -n "${!env:-}" ] && export "$env"
    done
    return 0
}

# ---- bookkeeping -------------------------------------------------------------------------------
# Seconds a deck with -var wall_timeout may run: SLURM end minus WALL_MARGIN (0 = unlimited).
wall_timeout_seconds() {
    local budget=0
    if [ -n "${SLURM_JOB_END_TIME:-}" ]; then
        budget=$(( SLURM_JOB_END_TIME - $(date +%s) - WALL_MARGIN ))
        [ "$budget" -gt 0 ] || budget=0
    fi
    echo "$budget"
}

# make_run_dirs <work_dir> <traj_dir>: the output tree, trajectories on scratch behind a symlink,
# and the README.txt that explains the tree to whoever opens the folder later.
make_run_dirs() {
    mkdir -p "$1"/{data_files,output_files/{stress_data,volume_data,piston_data,permeation_data,displacement_data,pair_data,chemical_potential},output_plots/convergence_plots}
    mkdir -p "$2"
    ln -s "$2" "$1/traj_files"
    write_run_readme "$1"
}

# write_run_readme <work_dir>: the layout of a run folder, in the folder. Keep it in step with
# docs/repository_layout.md ("Inside a run folder") when a deck or plotter changes where it writes.
write_run_readme() {
    cat > "$1/README.txt" <<'READMEEOF'
WHAT IS IN THIS RUN FOLDER
(written when the run started; the same guide with examples: lammps_work/docs/repository_layout.md)

  README.txt      this file
  run_info.txt    what was run: SLURM job, deck, input data file, every -var given to LAMMPS, git commit
  log.lammps      the LAMMPS log: thermo output and the deck's own ">>>" progress lines
  <deck>.lmp      the copy of the input deck that this run executed
  data_files/     EVERY .data and .restart file of the run (configurations: read_data / read_restart)
  output_files/   EVERY measured series and profile (.dat text tables; .dump for per-pair data)
  output_plots/   EVERY figure made from this run (.png)
  traj_files/     trajectories (.lammpstrj); a symlink to scratch, which is purged: copy what must last
  continuation_<stamp>/   a continuation of this run (continue_sim.sh); same layout inside

data_files/
  <dataname>.data                            the INPUT configuration (a symlink into lammps_data/input_data)
  final_<deck>_<stem>.data / .restart        the state at the END of the run
  final_config_<stem>_snap<N>.data           snapshots during an equilibration (resume points)
  strained_<deck>_<stem>_<stage>.restart     the strained state at the end of a level's hold, before any
                                             unload: the start of a repeated unload / D_c measurement
  <deck>: config = equilibration (slab_with_support, chain_length_testing), tricomp = triaxial
  compression, triperm = triaxial permeation, shear, compress = six-plate bulk compression;
  the pure-phase runs write purepol_<stem>.data / puresolv_<stem>.data.

output_files/
  stress_data/         stress profiles and series, strain logs, reservoir pressures
  volume_data/         gel and box dimensions and volumes, polymer centre of mass
  piston_data/         piston / plate positions, velocities, forces, pressures
  permeation_data/     solvent crossing the gel (flux, expelled volume)
  displacement_data/   polymer displacement profiles
  pair_data/           per-pair and per-bond dumps (pore and bond analysis)
  chemical_potential/  solvent density profiles

output_plots/
  <stem>_<figure>.png                        the results figures of postprocess.sh
  convergence_plots/                         run-health figures: T, P and volume convergence, diagnostics

FILE NAMES:  <quantity>_<stem>[_<stage>].<ext>
  <stem>  = <dataname>_<epsSS>_<epsSP>_<steps>   (<steps> = the length of one hold / of the production)
  <stage> = _ref        zero-load reference window, before any loading
            _c<strain>  hold at that compression level        _g<strain>  hold at that shear level
            _u<strain>  unload after that level (re-swelling / relaxation; permeation: _u)
            _v, _f      permeation: velocity-controlled stage, frozen-network stage
  no stage = the whole run.

Runs started before 2026-11-17 have their final_* / snapshot files at the top level of the run
folder instead of in data_files/, and no README.txt, run_info.txt or strained_* restart.
READMEEOF
}

# write_run_info <work_dir> <deck_src> <data_src>: run_info.txt, the provenance of the run (reads
# LAMMPS_VARS). The SLURM output file (the job's stdout: guards, verdict, post-processing) stays
# in the folder the job was submitted from; its name is recorded here.
write_run_info() {
    {
        echo "started:      $(date '+%Y-%m-%d %H:%M:%S') on $(hostname)"
        echo "slurm job:    ${SLURM_JOB_ID:-none} (${SLURM_JOB_NAME:-not a batch job})"
        echo "slurm output: ${SLURM_SUBMIT_DIR:-$PWD}/${SLURM_JOB_NAME:-?}.o${SLURM_JOB_ID:-?}.<node>"
        echo "allocation:   ${SLURM_NNODES:-?} node(s) x ${SLURM_NTASKS_PER_NODE:-?} tasks x ${OMP_NUM_THREADS:-1} thread(s)"
        echo "deck:         $2"
        echo "lammps_work:  commit $(git -C "$LAMMPS_WORK_DIR" rev-parse --short HEAD 2>/dev/null || echo '?')$(git -C "$LAMMPS_WORK_DIR" diff --quiet HEAD -- "$2" 2>/dev/null || echo ' + uncommitted changes to the deck (the copy in this folder is what ran)')"
        echo "input data:   $3"
        echo "lammps:       $LMP_BIN"
        echo "lammps vars:  ${LAMMPS_VARS[*]}"
    } > "$1/run_info.txt"
}

# ---- launching ---------------------------------------------------------------------------------
# mpirun_cmd: fills the array MPIRUN_CMD with the mpirun prefix for this allocation.
mpirun_cmd() {
    local ntasks="${SLURM_NTASKS:-1}" threads="${OMP_NUM_THREADS:-1}" mapspec
    local per_node="${SLURM_NTASKS_PER_NODE:-$(( ntasks / ${SLURM_NNODES:-1} ))}"
    case "$MPI_MAP_BY" in
        core) mapspec="ppr:${per_node}:node:pe=${threads}" ;;   # fill each node before the next
        *)    mapspec="node:pe=${threads}" ;;                   # round-robin over the nodes
    esac
    MPIRUN_CMD=(mpirun -n "$ntasks" --bind-to "$MPI_BIND_TO" --map-by "$mapspec"
                "$LMP_BIN" -sf omp -pk omp "$threads")
}

# launch_lammps <deck> <-var ...>: run LAMMPS in the current directory (log.lammps), then judge it.
# Sets LAMMPS_RC and LAMMPS_STATUS = complete | partial | failed.
#   * mpirun runs under `timeout`, ending 5 min before the SLURM limit: ranks sometimes hang in
#     MPI_Finalize / UCX teardown after LAMMPS has written everything, and a hung mpirun would
#     otherwise swallow the post-processing.
#   * Completion is judged from the log, not the exit code: LAMMPS prints "Total wall time" last,
#     after every write_data, so once it is there all output is on disk. A kill during teardown
#     must not discard a good run.
#   * RC 124 without that line means the guard fired while LAMMPS was still running (the job
#     outran its limit); the output written so far is still worth plotting ("partial").
launch_lammps() {
    local deck="$1"; shift
    local log=log.lammps guard=() pid waited remain
    if [ -n "${SLURM_JOB_END_TIME:-}" ]; then
        remain=$(( SLURM_JOB_END_TIME - $(date +%s) - 300 ))
        [ "$remain" -gt 60 ] && guard=(timeout -k 60 "${remain}s")
    fi
    mpirun_cmd
    echo "${guard[*]} ${MPIRUN_CMD[*]} ... -in $(basename "$deck")"
    "${guard[@]}" "${MPIRUN_CMD[@]}" "$@" -in "$deck" &
    pid=$!
    while kill -0 "$pid" 2>/dev/null; do
        if grep -q "Total wall time" "$log" 2>/dev/null; then
            waited=0
            while kill -0 "$pid" 2>/dev/null && [ "$waited" -lt 180 ]; do sleep 10; waited=$((waited + 10)); done
            if kill -0 "$pid" 2>/dev/null; then
                echo ">>> LAMMPS finished but mpirun is still up after 180 s: killing it (MPI teardown hang)."
                kill -TERM "$pid" 2>/dev/null; sleep 10; kill -KILL "$pid" 2>/dev/null; sleep 2
            fi
            break
        fi
        sleep 15
    done
    wait "$pid" 2>/dev/null
    LAMMPS_RC=$?
    if grep -q "Total wall time" "$log" 2>/dev/null; then
        LAMMPS_STATUS=complete
        [ "$LAMMPS_RC" -ne 0 ] && echo ">>> mpirun returned $LAMMPS_RC after 'Total wall time': output complete, continuing."
    elif [ "$LAMMPS_RC" -eq 124 ] && [ -d output_files ]; then
        LAMMPS_STATUS=partial
        echo ">>> PARTIAL RUN: the wall-clock guard killed LAMMPS mid-run. Plotting what was written; the last phase is truncated."
    else
        LAMMPS_STATUS=failed
        echo ">>> LAMMPS did not reach 'Total wall time' (mpirun RC $LAMMPS_RC): failure, no post-processing."
        [ "$LAMMPS_RC" -eq 0 ] && LAMMPS_RC=1
    fi
    return 0
}

# chain_if_incomplete <sim_folder>: with CHAIN=1 and a deck that reported STATUS INCOMPLETE,
# submit the continuation now (before post-processing, so it never competes with the wall limit).
chain_if_incomplete() {
    [ "${CHAIN:-0}" = 1 ] || return 0
    if grep -qE '^>>> PRODUCTION_STEPS_DONE [0-9]+ .*STATUS INCOMPLETE' log.lammps; then
        echo ">>> Production incomplete (wall-clock budget): submitting the continuation..."
        export_knobs
        ( cd "$1" && bash "$SCRIPT_DIR/continue_sim.sh" --chain "$SLURM_JOB_ID" ) \
            || echo ">>> WARNING: continuation not submitted. By hand: cd $1 && continue_sim.sh --chain $SLURM_JOB_ID"
    else
        echo ">>> CHAIN=1: production complete (or the deck has no wall-clock budget); nothing to chain."
    fi
}
