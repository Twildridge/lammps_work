#!/bin/bash
# compress_slab_sequence.sh -- the isolated-cube bulk-modulus pipeline, start to finish:
#
#   1. slab_with_support   free-swelling equilibration of an ISOLATED 5-bead rot45 cube in its bath
#                          (BARO_MODE=iso, PIN_GEL=1, PISTON_TRANSPARENT=1, P* = 1.5)
#   2. add_more_plates_to_gel.py --layout six_plates  on that run's final_config
#   3. compress_slab       drained isotropic compression of the plated file
#
#   compress_slab_sequence.sh            Mac (no sbatch): smoke test, all three steps in a row with
#                                        short step counts, LAMMPS GUI binary + OpenMP threads
#   compress_slab_sequence.sh            Expanse (sbatch found): production. Job A = step 1; job B
#                                        (--dependency=afterok:A) = steps 2 + 3
#   DRY_RUN=1 compress_slab_sequence.sh  print what would run / be submitted, change nothing
#
# Input data (step 1): slab_support_5beads_square_rho04_nowalls_rot45.data, the 13 x 9 x 13 BCT-cell
# rot45 cube (48.5 x 47.5 x 48.5 sigma gel, 15 sigma of rho 0.4 bath on every side, 236 k atoms),
# built with lib/data_file_generator.build('rot45', units_x=13, units_y=9, units_z=13).
#   Mac:     built into lammps_data_files_local/ if it is not there yet.
#   Expanse: copy it by hand into ~/Documents/lammps_data/input_data/ first; the script stops if it
#            is missing (it does not build on the cluster).
#
# Mac outputs: run dirs under lammps_runs_local/compress_slab_sequence_<stamp>/ (next to lammps_work)
# and the plated smoke file <...>_smoke_six_plates.data in lammps_data_files_local/. The smoke file
# comes from a few thousand steps of equilibration, so it tests the plumbing only; it is NOT a
# compress_slab production input.
#
# Written for the macOS bash 3.2 as well as Expanse's bash (no ${x^^}, no associative arrays).
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LAMMPS_WORK_DIR="${LAMMPS_WORK_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
export LAMMPS_WORK_DIR

# ---- what is run (both machines) ---------------------------------------------------------------
DATANAME="${DATANAME:-slab_support_5beads_square_rho04_nowalls_rot45}"
INTERACTION="${INTERACTION:-1.0_1.0}"
PRESS_TARGET="${PRESS_TARGET:-1.5}"
CUBE_UNITS="${CUBE_UNITS:-13 9 13}"            # units_x units_y units_z for the generator (rot45: x, z in BCT cells)

# ---- Expanse production lengths and allocations ------------------------------------------------
SLAB_NSTEPS="${SLAB_NSTEPS:-8000000}"          # as slab_with_support_single.batch
SLAB_SBATCH="${SLAB_SBATCH:---partition=compute --account=csb197 --nodes=6 --ntasks-per-node=128 --cpus-per-task=1 --mem=0 --time=16:00:00}"
COMP_NSTEPS="${COMP_NSTEPS:-3000000}"          # as compress_slab_single.batch
COMP_SBATCH="${COMP_SBATCH:---partition=compute --account=csb197 --nodes=6 --ntasks-per-node=128 --cpus-per-task=1 --mem=0 --time=27:00:00}"

# ---- Mac smoke test ----------------------------------------------------------------------------
LMP_LOCAL="${LMP_LOCAL:-/Applications/LAMMPS-GUI.app/Contents/bin/lmp}"
NTHREADS="${NTHREADS:-8}"
PYTHON="${PYTHON:-python3}"
LOCAL_DATA_DIR="${LOCAL_DATA_DIR:-$(cd "$LAMMPS_WORK_DIR/.." && pwd)/lammps_data_files_local}"
LOCAL_RUNS_DIR="${LOCAL_RUNS_DIR:-$(cd "$LAMMPS_WORK_DIR/.." && pwd)/lammps_runs_local}"
SMOKE_SLAB_NSTEPS="${SMOKE_SLAB_NSTEPS:-3000}"  # after the fixed ~20 k start-up steps of the deck
SMOKE_COMP_NSTEPS="${SMOKE_COMP_NSTEPS:-2000}"
# compress_slab smoke line from docs/running_simulations.md (verified there 2026-10-28)
SMOKE_COMP_VARS=(-var phase0_steps 200 -var npt_piston_steps 400 -var min_iter 50 -var t_seat 1000
    -var nsteps_settle 400 -var ref_avg_steps 400 -var ref_nfreq 200 -var compressions 0.02
    -var compressions_list "0.02" -var hold_auto 0 -var unload 1 -var unload_steps 2000
    -var unload_margin 1.0 -var v_piston_prod 0.4 -var thermo_freq 500 -var volume_freq 200
    -var strain_freq 100 -var aniso_nevery 100 -var aniso_nrepeat 2 -var aniso_nfreq 200)

EPSSS="${INTERACTION%_*}"; EPSSP="${INTERACTION#*_}"
DRY="${DRY_RUN:-0}"
say() { echo ">>> $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

# run_local <step label> <work dir> <deck> <data file> <-var ...>: one LAMMPS run in its own dir,
# laid out like run_lammps.sh makes it (data_files/, output_files/*, traj_files/).
run_local() {
    local label="$1" work="$2" deck="$3" data="$4"; shift 4
    say "$label: $(basename "$deck") on $(basename "$data") in $work"
    if [ "$DRY" = 1 ]; then echo "    DRY RUN: $LMP_LOCAL -sf omp -pk omp $NTHREADS $* -in $(basename "$deck")"; return 0; fi
    mkdir -p "$work"/data_files "$work"/traj_files "$work"/output_plots
    mkdir -p "$work"/output_files/stress_data "$work"/output_files/volume_data "$work"/output_files/piston_data \
             "$work"/output_files/permeation_data "$work"/output_files/displacement_data \
             "$work"/output_files/pair_data "$work"/output_files/chemical_potential
    cp "$deck" "$work/"
    ln -sf "$data" "$work/data_files/$(basename "$data")"
    local t0; t0=$(date +%s)
    ( cd "$work" && "$LMP_LOCAL" -sf omp -pk omp "$NTHREADS" "$@" -in "$(basename "$deck")" \
          -log log.lammps -screen none )
    local rc=$?
    grep -q "Total wall time" "$work/log.lammps" 2>/dev/null \
        || die "$label: LAMMPS stopped early (rc $rc). Last lines of $work/log.lammps:
$(tail -n 25 "$work/log.lammps" 2>/dev/null)"
    say "$label: done in $(( $(date +%s) - t0 )) s"
}

# ================================================================================================
if ! command -v sbatch >/dev/null 2>&1; then
    # ------------------------------------------------------------------ Mac: smoke test, inline
    [ -x "$LMP_LOCAL" ] || die "LAMMPS binary not found: $LMP_LOCAL (set LMP_LOCAL)"
    STAMP=$(date +%Y%m%d_%H%M%S)
    SEQ_DIR="$LOCAL_RUNS_DIR/compress_slab_sequence_$STAMP"
    INPUT="$LOCAL_DATA_DIR/$DATANAME.data"
    say "Mac smoke test: slab $SMOKE_SLAB_NSTEPS steps, compress_slab $SMOKE_COMP_NSTEPS steps, $NTHREADS OpenMP threads"
    say "run dirs: $SEQ_DIR"

    # step 0: the input cube
    if [ -f "$INPUT" ]; then
        say "input: $INPUT (exists, not rebuilt)"
    else
        say "building $INPUT (rot45, units $CUBE_UNITS)"
        set -- $CUBE_UNITS
        [ "$DRY" = 1 ] || ( cd "$LAMMPS_WORK_DIR/scripts/lib" && "$PYTHON" -c "
import data_file_generator as g
g.build('rot45', output_file='$INPUT', units_x=$1, units_y=$2, units_z=$3)" ) || die "data file build failed"
    fi

    # step 1: slab_with_support
    SLAB_DIR="$SEQ_DIR/slab_with_support"
    run_local "step 1/3 slab_with_support" "$SLAB_DIR" \
        "$LAMMPS_WORK_DIR/simulations/slab_with_support/slab_with_support.lmp" "$INPUT" \
        -var dataname "$DATANAME" -var interaction "$INTERACTION" -var epsSS "$EPSSS" -var epsSP "$EPSSP" \
        -var nsteps "$SMOKE_SLAB_NSTEPS" -var oldsteps 0 -var totsteps "$SMOKE_SLAB_NSTEPS" -var wall_timeout 0 \
        -var press_target "$PRESS_TARGET" -var baro_mode iso -var pin_gel 1 -var piston_transparent 1
    FINAL="$SLAB_DIR/final_config_${DATANAME}_${INTERACTION}_${SMOKE_SLAB_NSTEPS}.data"
    [ "$DRY" = 1 ] || [ -f "$FINAL" ] || die "step 1 finished without $FINAL"

    # step 2: six plates
    PLATED_STEM="final_config_${DATANAME}_${INTERACTION}_${SMOKE_SLAB_NSTEPS}_smoke_six_plates"
    PLATED="$LOCAL_DATA_DIR/$PLATED_STEM.data"
    say "step 2/3 add_more_plates_to_gel.py --layout six_plates -> $PLATED"
    if [ "$DRY" != 1 ]; then
        ( cd "$LAMMPS_WORK_DIR/scripts/lib" && "$PYTHON" add_more_plates_to_gel.py --layout six_plates \
              --input "$FINAL" --output "$PLATED" --no-log ) | tee "$SEQ_DIR/add_more_plates.log"
        [ "${PIPESTATUS[0]}" -eq 0 ] || die "converter failed (see $SEQ_DIR/add_more_plates.log)"
        grep -q "RESULT: OK" "$SEQ_DIR/add_more_plates.log" || die "converter validation did not report RESULT: OK"
    fi

    # step 3: compress_slab
    run_local "step 3/3 compress_slab" "$SEQ_DIR/compress_slab" \
        "$LAMMPS_WORK_DIR/simulations/compress_slab/compress_slab.lmp" "$PLATED" \
        -var dataname "$PLATED_STEM" -var interaction "$INTERACTION" -var epsSS "$EPSSS" -var epsSP "$EPSSP" \
        -var nsteps "$SMOKE_COMP_NSTEPS" -var oldsteps 0 -var totsteps "$SMOKE_COMP_NSTEPS" -var wall_timeout 0 \
        -var press_target "$PRESS_TARGET" "${SMOKE_COMP_VARS[@]}"
    if [ "$DRY" != 1 ]; then
        # match printed lines only: LAMMPS also echoes the deck's `if` / `print` commands into the log
        CLOG="$SEQ_DIR/compress_slab/log.lammps"
        grep '^>>> \[unload .*retraction done' "$CLOG" | sed 's/^/    /'
        grep '^>>> *WARNING' "$CLOG" | sed 's/^/    deck /'
        grep -q '^Bulk compression .*sweep complete' "$CLOG" || die "compress_slab log has no 'sweep complete' line"
        say "files tagged _u0.02: $(ls "$SEQ_DIR"/compress_slab/output_files/*/*_u0.02* 2>/dev/null | wc -l | tr -d ' ')"
    fi
    [ "$DRY" = 1 ] && { say "DRY RUN: nothing was run."; exit 0; }
    say "SMOKE TEST PASSED. For Expanse copy $INPUT to ~/Documents/lammps_data/input_data/ and run this script there."
    exit 0
fi

# ================================================================================================
# ---------------------------------------------------------------------- Expanse: two chained jobs
source "$LAMMPS_WORK_DIR/scripts/lib/site.sh"
[ -f "$INPUT_DATA_DIR/$DATANAME.data" ] \
    || die "$INPUT_DATA_DIR/$DATANAME.data not found: copy it from lammps_data_files_local first"
STAMP=$(date +%Y%m%d_%H%M%S)
SEQ_DIR="$RUNS_DIR/compress_slab_sequence/$STAMP"     # job B finds job A's final_config here
FINAL_NAME="final_config_${DATANAME}_${INTERACTION}_${SLAB_NSTEPS}.data"
PLATED_STEM="final_config_${DATANAME}_${INTERACTION}_${SLAB_NSTEPS}_six_plates"
PROLOGUE="source \"$LAMMPS_WORK_DIR/scripts/batch_prologue.sh\""

# Job A: as slab_with_support_single.batch, without CHAIN (a continuation would break afterok:
# the dependency fires when the first job ends). The deck writes the un-suffixed final_config only
# when production completed, so job B checks for it and stops if A ran out of time.
JOB_A="$PROLOGUE
export PRESS_TARGET=$PRESS_TARGET BARO_MODE=iso PIN_GEL=1 PISTON_TRANSPARENT=1
export LAMMPS_RUNS_OVERRIDE=\"$SEQ_DIR\"
bash \"\$LAMMPS_WORK_DIR/scripts/run_lammps.sh\" slab_with_support $DATANAME $INTERACTION $SLAB_NSTEPS"

# Job B: convert, put the plated file into input_data, then compress_slab as compress_slab_single.batch.
JOB_B="$PROLOGUE
set -o pipefail
FINAL=\$(find \"$SEQ_DIR\" -maxdepth 2 -name \"$FINAL_NAME\" | head -1)
[ -n \"\$FINAL\" ] || { echo \"ERROR: no $FINAL_NAME under $SEQ_DIR (step 1 incomplete?)\"; exit 1; }
load_python_module
( cd \"\$LAMMPS_WORK_DIR/scripts/lib\" && python add_more_plates_to_gel.py --layout six_plates \\
      --input \"\$FINAL\" --output \"$INPUT_DATA_DIR/$PLATED_STEM.data\" ) | tee \"$SEQ_DIR/add_more_plates.log\" || exit 1
grep -q 'RESULT: OK' \"$SEQ_DIR/add_more_plates.log\" || { echo 'ERROR: converter validation failed'; exit 1; }
load_lammps_modules
export PRESS_TARGET=$PRESS_TARGET COMPRESSIONS=0.10 UNLOAD=1 UNLOAD_STEPS=0 UNLOAD_MARGIN=3.0
export LAMMPS_RUNS_OVERRIDE=\"$SEQ_DIR\"
bash \"\$LAMMPS_WORK_DIR/scripts/run_lammps.sh\" compress_slab $PLATED_STEM $INTERACTION $COMP_NSTEPS"

say "Expanse: slab_with_support ($SLAB_NSTEPS steps) -> six plates -> compress_slab ($COMP_NSTEPS steps)"
say "run dir: $SEQ_DIR"
if [ "$DRY" = 1 ]; then
    echo "DRY RUN job A: sbatch $SLAB_SBATCH --job-name=seq_slab_rot45_cube job_a_slab_with_support.sh"
    printf '#!/bin/bash\n%s\n' "$JOB_A" | sed 's/^/    /'
    echo "DRY RUN job B: sbatch $COMP_SBATCH --job-name=seq_compress_rot45_cube --dependency=afterok:<A> job_b_plates_compress_slab.sh"
    printf '#!/bin/bash\n%s\n' "$JOB_B" | sed 's/^/    /'
    exit 0
fi
mkdir -p "$SEQ_DIR"
printf '#!/bin/bash\n%s\n' "$JOB_A" > "$SEQ_DIR/job_a_slab_with_support.sh"
printf '#!/bin/bash\n%s\n' "$JOB_B" > "$SEQ_DIR/job_b_plates_compress_slab.sh"
cd "$SEQ_DIR" || exit 1
# shellcheck disable=SC2086  # the allocation strings are word-split on purpose
A=$(sbatch --parsable $SLAB_SBATCH --job-name=seq_slab_rot45_cube --output="%x.o%j.%N" \
        job_a_slab_with_support.sh) || die "sbatch of job A failed"
A=${A%%;*}
# shellcheck disable=SC2086
B=$(sbatch --parsable $COMP_SBATCH --job-name=seq_compress_rot45_cube --output="%x.o%j.%N" \
        --dependency=afterok:"$A" --kill-on-invalid-dep=yes job_b_plates_compress_slab.sh) \
    || die "sbatch of job B failed (job A $A is queued: scancel $A if you do not want it alone)"
B=${B%%;*}
say "submitted: A = $A (slab_with_support), B = $B (plates + compress_slab, after A succeeds)"
say "watch: squeue -u \$USER ; logs in $SEQ_DIR"
