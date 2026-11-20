#!/usr/bin/env bash
# node_screen.sh -- health screen of one Pod node (2026-11-19).  Submitted by submit_node_screen.sh,
# one job per node and mode; do not sbatch it by hand.
#   SCREEN_MODE=core   one core, shares the node: kernel, uptime, clocks, temperatures, throttle
#                      counters, stray processes, scratch write speed, a 1-rank LJ run
#   SCREEN_MODE=node   the whole node (--exclusive): the same, plus LJ runs with 40 and 39 ranks and
#                      1500 steps of the N15 rot45 benchmark deck from its prepared state, with clocks
#                      and temperatures sampled every 5 s while LAMMPS runs.  39 ranks leave one core
#                      free: on a node where a system process spins at 100% (the Prometheus node
#                      exporter does, 2026-11-19) 40 ranks share a core with it and wait in Comm.
# SCREEN_CHAIN=N (submit_node_screen.sh --chain N): if the Prometheus node exporter is spinning on the
# node (CSC is restarting it, 2026-11-20) the run is postponed 6 h, up to N times, instead of measuring a
# node that is about to change; a core-mode job that finds the exporter quiet submits the whole-node run.
# Every LAMMPS run has its own timeout, the deck gets what is left of the job's time, and the job
# writes its row even when Slurm's time limit hits (sbatch --signal=B:USR1@120 from the submitter).
# Output: $SCREEN_OUT/<node>/<mode>/report.txt and the LAMMPS logs, plus one line appended to
# $SCREEN_OUT/rows.tsv (summarize_node_screen.py turns those into the ranking).
# Speed numbers are core-microseconds per atom-step: Loop time x ranks / (steps x atoms), so one
# number compares 1-rank and 40-rank runs and different nodes.  Lower is better.
# SCREEN_TEST=1 shrinks every run for a check on the login node (no Slurm needed).
#SBATCH --output=logs_slurm/%x.o%j.%N

set -o pipefail            # no set -u: Lmod's init reads unset Slurm variables (Pod)

MODE="${SCREEN_MODE:-core}"
OUT="${SCREEN_OUT:?SCREEN_OUT = output folder of the screen}"
TEST="${SCREEN_TEST:-0}"
NODE=$(hostname -s)
JOB="${SLURM_JOB_ID:-local$$}"
D="$OUT/$NODE/$MODE"
mkdir -p "$D" || exit 1
R="$D/report.txt"; : > "$R"
log() { printf '%s\n' "$*" >> "$R"; }
sec() { log ""; log "== $* =="; }
T_START=$(date +%s)
SCHEMA=v2                 # rows_v2.tsv; summarize_node_screen.py reads every rows*.tsv
STATUS=ok
LMP_PID=
on_signal() {             # Slurm's warning (USR1, 120 s before the limit) or TERM: stop LAMMPS, keep the row
    log ""; log "signal received at $(date -Is): stopping LAMMPS and writing the row"
    [ -n "$LMP_PID" ] && { kill -TERM "$LMP_PID" 2>/dev/null; sleep 5; kill -KILL "$LMP_PID" 2>/dev/null; }
    STATUS="${STATUS}+timeout"
    write_row; exit 1
}
trap on_signal USR1 TERM

source "$HOME/Documents/lammps_work/scripts/lib/site.sh" || exit 1
mpi_env
load_lammps_modules
LMP="$LMP_BIN"
BENCH_DIR="$(dirname "$(dirname "$LMP")")/bench"          # .../lammps-22Jul2025/bench
PREP="$INPUT_DATA_DIR/allocation_benchmark/prepared/prep_15b_s1.0_rot45.data"
DECK_SRC="$LAMMPS_WORK_DIR/simulations/slab_with_support/slab_with_support.lmp"

# ---- sensors -------------------------------------------------------------------------------------
mhz() {   # mean min max of the current core clocks (MHz)
    awk '/^cpu MHz/ {s+=$4; n++; if ($4>mx) mx=$4; if (mn=="" || $4<mn) mn=$4}
         END {if (n) printf "%.0f %.0f %.0f", s/n, mn, mx; else print "NA NA NA"}' /proc/cpuinfo
}
temps() { # max package temperature, max core temperature (C) from coretemp
    local h f l v pk=-1 co=-1
    for h in /sys/class/hwmon/hwmon*; do
        [ "$(cat "$h/name" 2>/dev/null)" = coretemp ] || continue
        for f in "$h"/temp*_input; do
            v=$(cat "$f" 2>/dev/null) || continue; v=$((v / 1000))
            l=$(cat "${f%_input}_label" 2>/dev/null)
            case "$l" in Package*) [ "$v" -gt "$pk" ] && pk=$v ;; *) [ "$v" -gt "$co" ] && co=$v ;; esac
        done
    done
    echo "$pk $co"
}
throttle() { # sum of core throttle events, max package throttle events, max package throttle time (ms), since boot
    local f v c=0 p=0 t=0
    for f in /sys/devices/system/cpu/cpu[0-9]*/thermal_throttle; do
        v=$(cat "$f/core_throttle_count" 2>/dev/null) && c=$((c + v))
        v=$(cat "$f/package_throttle_count" 2>/dev/null) && [ "$v" -gt "$p" ] && p=$v
        v=$(cat "$f/package_throttle_total_time_ms" 2>/dev/null) && [ "$v" -gt "$t" ] && t=$v
    done
    echo "$c $p $t"
}
busy() {  # system-wide CPU busy % over $1 seconds
    local a b; a=$(head -1 /proc/stat); sleep "$1"; b=$(head -1 /proc/stat)
    awk -v a="$a" -v b="$b" 'BEGIN {split(a, x, " "); split(b, y, " "); ta=0; tb=0
         for (i=2; i<=9; i++) {ta+=x[i]; tb+=y[i]}; d=tb-ta; di=(y[5]+y[6])-(x[5]+x[6])
         if (d>0) printf "%.1f", 100*(d-di)/d; else print "NA"}'
}
SAMPLER=
sampler_start() {   # $1 = name; one line per 5 s: epoch, mhz mean min max, pkg C, core C, throttle c p ms
    ( while :; do echo "$(date +%s) $(mhz) $(temps) $(throttle)"; sleep 5; done ) > "$D/samples_$1.tsv" 2>/dev/null &
    SAMPLER=$!
}
sampler_stop() { kill "$SAMPLER" 2>/dev/null; wait "$SAMPLER" 2>/dev/null; SAMPLER=; }
sample_summary() {  # $1 = name -> "mhz_mean mhz_max pkg_max" over the run (NA if nothing sampled)
    awk '{s+=$2; n++; if ($4>mx) mx=$4; if ($5>pk) pk=$5} END {if (n) printf "%.0f %.0f %d", s/n, mx, pk; else print "NA NA NA"}' "$D/samples_$1.tsv" 2>/dev/null
}

# ---- static facts --------------------------------------------------------------------------------
sec "node"
log "node=$NODE job=$JOB mode=$MODE date=$(date -Is)"
log "kernel=$(uname -r)"
UPTIME_D=$(awk '{printf "%.1f", $1/86400}' /proc/uptime); log "uptime_days=$UPTIME_D"
log "cmdline=$(cat /proc/cmdline)"
log "cpus_allowed=$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')"
ONLINE=$(grep -c ^processor /proc/cpuinfo); log "cpus_online=$ONLINE"
SMT=$(cat /sys/devices/system/cpu/smt/active 2>/dev/null || echo NA); log "smt_active=$SMT"
GOV=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo NA)
log "cpufreq: driver=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_driver 2>/dev/null || echo NA) governor=$GOV no_turbo=$(cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null || echo NA) max_perf_pct=$(cat /sys/devices/system/cpu/intel_pstate/max_perf_pct 2>/dev/null || echo NA)"
log "model=$(awk -F: '/model name/ {print $2; exit}' /proc/cpuinfo)"
log "microcode=$(awk -F: '/microcode/ {print $2; exit}' /proc/cpuinfo)"
log "thp=$(cat /sys/kernel/mm/transparent_hugepage/enabled 2>/dev/null) numa_balancing=$(cat /proc/sys/kernel/numa_balancing 2>/dev/null)"
sec "memory"
grep -E "MemTotal|MemAvailable|MemFree|Cached|HugePages_Total|Dirty" /proc/meminfo >> "$R"
for n in /sys/devices/system/node/node[0-9]*; do log "$(basename "$n"): $(grep -E 'MemTotal|MemFree' "$n/meminfo" | awk '{printf "%s %d MB  ", $3, $4/1024}')"; done
sec "omni-path hfi1_0 port 1"
P=/sys/class/infiniband/hfi1_0/ports/1
if [ -d "$P" ]; then
    OPA="$(cat $P/state 2>/dev/null | awk '{print $2}')"
    log "state=$(cat $P/state) rate=$(cat $P/rate) phys=$(cat $P/phys_state)"
    for c in link_downed link_error_recovery symbol_error port_rcv_errors port_xmit_discards; do
        [ -f "$P/counters/$c" ] && log "$c=$(cat $P/counters/$c)"
    done
else OPA=none; log "no hfi1_0"; fi

# ---- state before any load -----------------------------------------------------------------------
sec "before load"
read -r MHZ0 MHZ0_MIN MHZ0_MAX <<< "$(mhz)"
read -r PKG0 CORE0 <<< "$(temps)"
read -r THR_C THR_P THR_MS <<< "$(throttle)"
LOAD1=$(awk '{print $1}' /proc/loadavg)
BUSY0=$(busy 5)
log "clock_MHz mean=$MHZ0 min=$MHZ0_MIN max=$MHZ0_MAX"
log "temp_C package_max=$PKG0 core_max=$CORE0"
log "throttle core_events=$THR_C package_events=$THR_P package_time_ms=$THR_MS"
log "loadavg=$(cat /proc/loadavg)  busy_pct_5s=$BUSY0"
sec "processes by CPU (top, second pass)"
top -b -n 2 -d 3 -w 160 2>/dev/null | awk '/^top -/ {n++} n==2' | head -24 > "$D/top.txt"
cat "$D/top.txt" >> "$R"
# the busiest process that is not ours: a spinning daemon steals a core from a 40-rank run
read -r FOREIGN_PCT FOREIGN_CMD <<< "$(awk -v me="$USER" '$1 ~ /^[0-9]+$/ && $2 != me && $2 != substr(me,1,7)"+" {if ($9+0 > m) {m=$9+0; c=$12" ("$2")"}} END {printf "%.0f %s", m, (c==""?"none":c)}' "$D/top.txt")"
log "foreign_top_process=$FOREIGN_CMD cpu_pct=$FOREIGN_PCT"
EXPORTER_PCT=$(awk '$1 ~ /^[0-9]+$/ && ($12 ~ /prometh|node_exporter/ || $2 ~ /^prometh/) {if ($9+0 > m) m=$9+0} END {printf "%.0f", m}' "$D/top.txt")
log "exporter_cpu_pct=$EXPORTER_PCT"
CHAIN="${SCREEN_CHAIN:-}"
if [ -n "$CHAIN" ]; then
    SUBMIT="$LAMMPS_WORK_DIR/simulations/allocation_benchmark/submit_node_screen.sh"
    if [ "$EXPORTER_PCT" -ge 50 ] && [ "$CHAIN" -gt 0 ]; then
        log "node exporter still spinning: this $MODE run is postponed 6 h ($CHAIN tries left); no row written"
        mv "$R" "$D/report_postponed_$(date +%m%dT%H%M).txt"
        if [ "$TEST" = 1 ]; then echo "TEST: would postpone $MODE for $NODE (chain $((CHAIN - 1)))"
        else bash "$SUBMIT" --modes "$MODE" --nodes "$NODE" --out "$OUT" --chain $((CHAIN - 1)) --begin now+6hours --no-summary; fi
        exit 0
    fi
    if [ "$MODE" = core ] && [ "$EXPORTER_PCT" -lt 50 ]; then
        log "node exporter quiet: submitting the whole-node run for $NODE"
        if [ "$TEST" = 1 ]; then echo "TEST: would submit node mode for $NODE"
        else bash "$SUBMIT" --modes node --nodes "$NODE" --out "$OUT" --chain 0 >> "$R" 2>&1; fi
    fi
fi
sec "processes not mine, by accumulated CPU"
ps -eo user:14,pid,pcpu,pmem,etimes,comm --sort=-pcpu | awk -v me="$USER" 'NR==1 || $1!=me' | head -14 >> "$R"

sec "scratch write (512 MB, fsync)"
SCR="$SCRATCH_DIR/allocation_benchmark/node_screen"; mkdir -p "$SCR"
F="$SCR/dd_${NODE}_$JOB"; T0=$(date +%s.%N)
if dd if=/dev/zero of="$F" bs=1M count=512 conv=fsync status=none 2>>"$R"; then
    SCR_MBPS=$(awk -v t0="$T0" -v t1="$(date +%s.%N)" 'BEGIN {printf "%.0f", 512/(t1-t0)}')
else SCR_MBPS=NA; fi
rm -f "$F"; log "write_MB_per_s=$SCR_MBPS"

# ---- LAMMPS --------------------------------------------------------------------------------------
# run_lmp NAME NRANKS LIMIT_S ARGS...  -> globals LOOP PROCS STEPS ATOMS US (core-us/atom-step), the
# avg-time columns PAIR NEIGH COMM MODIFY OUTPUT and COMMFRAC (Comm / Loop); log in $D/NAME.log,
# stderr/bindings in $D/NAME.err.  LIMIT_S is the run's own timeout (exit 124).
run_lmp() {
    local name=$1 n=$2 limit=$3; shift 3
    local t0 t1 bind
    if [ "$n" -gt 1 ]; then bind=(--bind-to core --map-by "ppr:${n}:node:pe=1" --report-bindings); else bind=(--bind-to none); fi
    sampler_start "$name"; t0=$(date +%s)
    timeout -k 15 "$limit" mpirun -n "$n" "${bind[@]}" "$LMP" -sf omp -pk omp 1 "$@" -log "$D/$name.log" -screen none 2> "$D/$name.err" &
    LMP_PID=$!
    wait "$LMP_PID"; RC=$?; LMP_PID=
    t1=$(date +%s); sampler_stop
    [ "$RC" = 124 ] && log "$name: hit its own timeout of $limit s"
    pkill -u "$USER" -x lmp 2>/dev/null     # nothing of ours may linger into the next run
    read -r LOOP PROCS STEPS ATOMS <<< "$(awk '/^Loop time of/ {l=$4; p=$6; s=$9; a=$12} END {if (l!="") print l, p, s, a; else print "NA NA NA NA"}' "$D/$name.log" 2>/dev/null)"
    if [ "$LOOP" != NA ]; then US=$(awk -v l="$LOOP" -v p="$PROCS" -v s="$STEPS" -v a="$ATOMS" 'BEGIN {printf "%.3f", l*p*1e6/(s*a)}'); else US=NA; fi
    PAIR=$(tsec Pair "$name"); NEIGH=$(tsec Neigh "$name"); COMM=$(tsec Comm "$name"); MODIFY=$(tsec Modify "$name"); OUTPUT=$(tsec Output "$name")
    if [ "$LOOP" != NA ] && [ "$COMM" != NA ]; then COMMFRAC=$(awk -v c="$COMM" -v l="$LOOP" 'BEGIN {printf "%.2f", c/l}'); else COMMFRAC=NA; fi
    read -r MHZ_L MHZ_LMAX PKG_L <<< "$(sample_summary "$name")"
    log "$name: ranks=$n rc=$RC wall_s=$((t1 - t0)) loop_s=$LOOP steps=$STEPS atoms=$ATOMS core_us_per_atom_step=$US"
    log "$name: avg_s pair=$PAIR neigh=$NEIGH comm=$COMM modify=$MODIFY output=$OUTPUT comm_fraction=$COMMFRAC"
    log "$name: during run clock_MHz mean=$MHZ_L max=$MHZ_LMAX package_C_max=$PKG_L"
    grep -m 3 "bound to" "$D/$name.err" >> "$R" 2>/dev/null
}
tsec() { awk -F'|' -v k="$1" '$1 ~ "^"k" " {gsub(/ /, "", $3); v=$3} END {print (v=="" ? "NA" : v)}' "$D/$2.log" 2>/dev/null; }

LJ_STEPS=1500; LJ40_STEPS=300; DECK_STEPS=1500; NR=40; NR1=39
if [ "$TEST" = 1 ]; then LJ_STEPS=100; LJ40_STEPS=50; DECK_STEPS=50; NR=2; NR1=1; fi
sed "s/^run .*/run $LJ_STEPS/" "$BENCH_DIR/in.lj" > "$D/in.lj1"
sed "s/^run .*/run $LJ40_STEPS/" "$BENCH_DIR/in.lj" > "$D/in.lj40"
remaining() {   # seconds left before Slurm's limit (a large number when not under Slurm)
    if [ -n "$SLURM_JOB_END_TIME" ]; then echo $((SLURM_JOB_END_TIME - $(date +%s))); else echo 3600; fi
}

LJ1_US=NA; LJ1_MHZ=NA; LJ1_PKG=NA; LJ40_US=NA; LJ40_MHZ=NA; LJ40_PKG=NA; LJ40_CF=NA; LJ39_US=NA; LJ39_CF=NA
DECK_SPS=NA; DECK_US=NA; DECK_MHZ=NA; DECK_PKG=NA; DECK_CF=NA; D_PAIR=NA; D_NEIGH=NA; D_COMM=NA; D_MODIFY=NA; D_OUTPUT=NA

# ---- one row for the summary ---------------------------------------------------------------------
write_row() {
    read -r THR_C2 THR_P2 THR_MS2 <<< "$(throttle)"
    local hdr row
    hdr="node\tmode\tjob\tdate\tkernel\tuptime_d\tcpus_online\tsmt\tgovernor\topa\tmhz_idle_mean\tmhz_idle_max\tpkg_c_idle\tcore_c_idle\tthr_core\tthr_pkg\tthr_pkg_ms\tthr_pkg_ms_new\tbusy_idle_pct\tload1\tforeign_pct\tforeign_cmd\tscratch_mbps\tlj1_us\tlj1_mhz_max\tlj1_pkg_c\tlj40_us\tlj40_comm_frac\tlj40_mhz_mean\tlj40_pkg_c\tlj39_us\tlj39_comm_frac\tdeck_steps_s\tdeck_us\tdeck_comm_frac\tdeck_mhz_mean\tdeck_pkg_c\tdeck_pair_s\tdeck_neigh_s\tdeck_comm_s\tdeck_modify_s\tdeck_output_s\tjob_wall_s\tstatus"
    row="$NODE\t$MODE\t$JOB\t$(date +%FT%H:%M)\t$(uname -r)\t$UPTIME_D\t$ONLINE\t$SMT\t$GOV\t$OPA\t$MHZ0\t$MHZ0_MAX\t$PKG0\t$CORE0\t$THR_C\t$THR_P\t$THR_MS\t$((THR_MS2 - THR_MS))\t$BUSY0\t$LOAD1\t$FOREIGN_PCT\t${FOREIGN_CMD// /_}\t$SCR_MBPS\t$LJ1_US\t$LJ1_MHZ\t$LJ1_PKG\t$LJ40_US\t$LJ40_CF\t$LJ40_MHZ\t$LJ40_PKG\t$LJ39_US\t$LJ39_CF\t$DECK_SPS\t$DECK_US\t$DECK_CF\t$DECK_MHZ\t$DECK_PKG\t$D_PAIR\t$D_NEIGH\t$D_COMM\t$D_MODIFY\t$D_OUTPUT\t$(( $(date +%s) - T_START ))\t$STATUS"
    (
        flock 9
        [ -s "$OUT/rows_$SCHEMA.tsv" ] || printf '%b\n' "$hdr" > "$OUT/rows_$SCHEMA.tsv"
        printf '%b\n' "$row" >> "$OUT/rows_$SCHEMA.tsv"
    ) 9> "$OUT/rows.lock"
    log "row written to rows_$SCHEMA.tsv status=$STATUS"
}

sec "LAMMPS: LJ melt, 1 rank, 32000 atoms, $LJ_STEPS steps"
cd "$D" && run_lmp lj1 1 300 -in "$D/in.lj1"
LJ1_US=$US; LJ1_MHZ=$MHZ_LMAX; LJ1_PKG=$PKG_L; [ "$US" = NA ] && STATUS=lj1_failed

if [ "$MODE" = node ]; then
    sec "LAMMPS: LJ melt, $NR ranks, 32000 atoms per rank, $LJ40_STEPS steps"
    LJV=(-var x 2 -var y 2 -var z 10)        # 1,280,000 atoms, the same box for 40 and 39 ranks
    [ "$TEST" = 1 ] && LJV=(-var x 1 -var y 1 -var z 2)
    cd "$D" && run_lmp lj40 "$NR" 600 "${LJV[@]}" -in "$D/in.lj40"
    LJ40_US=$US; LJ40_MHZ=$MHZ_L; LJ40_PKG=$PKG_L; LJ40_CF=$COMMFRAC; [ "$US" = NA ] && STATUS=lj40_failed

    sec "LAMMPS: LJ melt, $NR1 ranks (one core left free), same box, $LJ40_STEPS steps"
    cd "$D" && run_lmp lj39 "$NR1" 600 "${LJV[@]}" -in "$D/in.lj40"
    LJ39_US=$US; LJ39_CF=$COMMFRAC; [ "$US" = NA ] && STATUS=lj39_failed

    sec "LAMMPS: benchmark deck N15 rot45 from prepared state, $NR ranks, $DECK_STEPS steps, cutoff 7.0"
    if [ -s "$PREP" ] && [ -s "$DECK_SRC" ]; then
        WORK="$SCR/${NODE}_${JOB}_deck"
        mkdir -p "$WORK"/{data_files,traj_files,output_files/volume_data,output_files/stress_data}
        cp "$DECK_SRC" "$WORK/deck.lmp"
        perl -i -pe 's/^comm_modify cutoff 5\.0\s*$/comm_modify cutoff \${comm_cutoff}\n/' "$WORK/deck.lmp"
        if grep -q '^comm_modify cutoff ${comm_cutoff}$' "$WORK/deck.lmp"; then
            ln -s "$PREP" "$WORK/data_files/prep_15b_s1.0_rot45.data"
            # the benchmark's own settings (run_benchmark.sh, ROLE=bench)
            VARS=(-var dataname prep_15b_s1.0_rot45 -var interaction 1.0_1.0 -var epsSS 1.0 -var epsSP 1.0
                  -var nsteps "$DECK_STEPS" -var oldsteps 0 -var totsteps "$DECK_STEPS"
                  -var press_target 1.5 -var vel_seed 12345
                  -var baro_mode aniso -var pin_gel 0 -var piston_transparent 1
                  -var wall_timeout 0 -var n_snaps 1 -var cont 1 -var comm_cutoff 7.0 -var mean_box 0)
            DECK_LIMIT=$(( $(remaining) - 150 ))
            if [ "$DECK_LIMIT" -lt 90 ]; then
                US=NA; STATUS=deck_skipped_no_time; log "only $(remaining) s left in the job: deck skipped"
            else
                cd "$WORK" && run_lmp deck "$NR" "$DECK_LIMIT" "${VARS[@]}" -in "$WORK/deck.lmp"
                cd "$D"
            fi
            if [ "$US" != NA ]; then
                DECK_US=$US; DECK_MHZ=$MHZ_L; DECK_PKG=$PKG_L; DECK_CF=$COMMFRAC
                DECK_SPS=$(awk -v l="$LOOP" -v s="$STEPS" 'BEGIN {printf "%.2f", s/l}')
                D_PAIR=$PAIR; D_NEIGH=$NEIGH; D_COMM=$COMM; D_MODIFY=$MODIFY; D_OUTPUT=$OUTPUT
                rm -rf "$WORK"
            else
                [ "$STATUS" = ok ] && STATUS=deck_failed
                log "deck gave no Loop time (status $STATUS); run directory kept: $WORK"; tail -5 "$D/deck.log" >> "$R" 2>/dev/null
            fi
        else STATUS=deck_pattern; log "ERROR: slab_with_support.lmp no longer matches the benchmark's comm_modify edit"; fi
    else STATUS=deck_missing; log "missing $PREP or $DECK_SRC"; fi
fi

# ---- after load ----------------------------------------------------------------------------------
sec "after load"
read -r THR_C2 THR_P2 THR_MS2 <<< "$(throttle)"
read -r PKG2 CORE2 <<< "$(temps)"
log "temp_C package_max=$PKG2 core_max=$CORE2"
log "throttle core_events=$THR_C2 package_events=$THR_P2 package_time_ms=$THR_MS2 (new during job: core=$((THR_C2 - THR_C)) package=$((THR_P2 - THR_P)) ms=$((THR_MS2 - THR_MS)))"
log "job_wall_s=$(( $(date +%s) - T_START )) status=$STATUS"

# ---- the row for the summary (write_row is defined above, before the LAMMPS runs) -------------
trap - USR1 TERM
write_row
echo "node_screen $NODE $MODE: $STATUS lj1=$LJ1_US lj40=$LJ40_US lj39=$LJ39_US deck=$DECK_US core-us/atom-step"
