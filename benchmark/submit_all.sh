#!/usr/bin/env bash
# =============================================================================
# submit_all.sh  --  submit the whole benchmark matrix to SLURM.
#
# For every preFilt level it submits one prefilter job, then submits
# (scripts x preFilt x callThresh) TSScall jobs, each depending on the matching
# prefilter job. Then a collate job (which WAITS for every condition's metrics
# row before building tables), then an overlap-analysis job.
#
#   conditions = SCRIPTS x PREFILT_LEVELS x CALLTHRESH_LEVELS
#
# Comparisons this enables:
#   * orig vs refac  : same (preFilt, callThresh), different script -> collate
#   * time & memory across callThresh : the CALLTHRESH_LEVELS axis
#   * with vs without prefilter, both scripts : the PREFILT_LEVELS axis (1 = off)
#   * peak-universe overlap across levels : analyze_overlap.py
#
# Usage:  ./submit_all.sh          # submit everything
#         ./submit_all.sh --dry    # print the plan and the sbatch commands only
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$HERE/config.sh"
source "$CONFIG"

DRY=0
[[ "${1:-}" == "--dry" ]] && DRY=1

# sanity checks (skipped in dry mode so you can preview before wiring paths)
if [[ $DRY -eq 0 ]]; then
    for f in "$ORIG_SCRIPT" "$REFAC_SCRIPT" "$FWD" "$REV" "$CHROMSIZES" "$GTF"; do
        [[ -e "$f" ]] || { echo "ERROR: missing input: $f" >&2; exit 1; }
    done
fi

mkdir -p "$RESULTS_DIR"/{inputs,logs,conditions}
MANIFEST="$RESULTS_DIR/manifest.tsv"
METRICS="$RESULTS_DIR/metrics.tsv"
: > "$MANIFEST"
printf 'script\tpreFilt\tcallThresh\trc\twall_s\tmaxrss_kb\tsacct_maxrss\tsacct_elapsed\tbed_lines\tbed_md5\tdetail_md5\n' > "$METRICS"

script_path() { case "$1" in orig) echo "$ORIG_SCRIPT";; refac) echo "$REFAC_SCRIPT";; esac; }

# Submit and return a clean job id. sbatch --parsable prints "jobid" or
# "jobid;cluster" on one line; take the first token, strip whitespace.
submit() {
    if [[ $DRY -eq 1 ]]; then echo "DRYRUN: $*" >&2; echo "DRYJOB"; return; fi
    local out
    out="$("$@")"
    out="${out//$'\n'/}"
    echo "${out%%;*}"
}

# ---- 1. prefilter jobs, one per preFilt level -------------------------------
declare -A PF_JOB PF_FWD PF_REV
for PF in "${PREFILT_LEVELS[@]}"; do
    OUT_FWD="$RESULTS_DIR/inputs/fwd.preFilt${PF}.bedGraph"
    OUT_REV="$RESULTS_DIR/inputs/rev.preFilt${PF}.bedGraph"
    PF_FWD[$PF]="$OUT_FWD"; PF_REV[$PF]="$OUT_REV"
    jid=$(submit sbatch --parsable \
        --job-name="preFilt${PF}" --partition="$PARTITION" \
        --mem="$PREFILT_MEM" --time="$PREFILT_TIME" --cpus-per-task=1 \
        --output="$RESULTS_DIR/logs/prefilter_preFilt${PF}.log" \
        --export=ALL,CONFIG="$CONFIG",PREFILT="$PF",OUT_FWD="$OUT_FWD",OUT_REV="$OUT_REV" \
        "$HERE/prefilter.sh")
    PF_JOB[$PF]="$jid"
    echo "prefilter preFilt=$PF -> job $jid"
done

# ---- 2. TSScall jobs, full matrix -------------------------------------------
ALL_RUN_JOBS=()
for TAG in "${SCRIPTS[@]}"; do
  SP="$(script_path "$TAG")"
  for PF in "${PREFILT_LEVELS[@]}"; do
    for CT in "${CALLTHRESH_LEVELS[@]}"; do
      COND="${TAG}_preFilt${PF}_callThresh${CT}"
      COND_DIR="$RESULTS_DIR/conditions/$COND"
      jid=$(submit sbatch --parsable \
        --job-name="$COND" --partition="$PARTITION" \
        --mem="$JOB_MEM" --time="$JOB_TIME" --cpus-per-task="$CPUS_PER_TASK" \
        --dependency="afterok:${PF_JOB[$PF]}" \
        --output="$RESULTS_DIR/logs/${COND}.log" \
        --export=ALL,CONFIG="$CONFIG",SCRIPT_TAG="$TAG",SCRIPT_PATH="$SP",PREFILT="$PF",CALLTHRESH="$CT",IN_FWD="${PF_FWD[$PF]}",IN_REV="${PF_REV[$PF]}",COND_DIR="$COND_DIR",METRICS_FILE="$METRICS" \
        "$HERE/run_one.sh")
      ALL_RUN_JOBS+=("$jid")
      printf '%s\t%s\t%s\t%s\t%s\n' "$COND" "$TAG" "$PF" "$CT" "$jid" >> "$MANIFEST"
      echo "run $COND -> job $jid (after ${PF_JOB[$PF]})"
    done
  done
done

# ---- 3. collate job ---------------------------------------------------------
# afterany gets it running as the matrix winds down; collate.sh then WAITS for
# the metrics table to be complete before building anything (this is the real
# guard against the early-collate race).
DEP=$(IFS=:; echo "${ALL_RUN_JOBS[*]}")
cjid=$(submit sbatch --parsable \
    --job-name="collate" --partition="$PARTITION" \
    --mem=4G --time=02:00:00 --cpus-per-task=1 \
    --dependency="afterany:${DEP}" \
    --output="$RESULTS_DIR/logs/collate.log" \
    --export=ALL,CONFIG="$CONFIG" \
    "$HERE/collate.sh")
echo "collate -> job $cjid (waits for complete metrics)"

# ---- 4. overlap analysis job ------------------------------------------------
ajid=$(submit sbatch --parsable \
    --job-name="overlap" --partition="$PARTITION" \
    --mem=8G --time=02:00:00 --cpus-per-task=1 \
    --dependency="afterok:${cjid}" \
    --output="$RESULTS_DIR/logs/overlap.log" \
    --export=ALL,CONFIG="$CONFIG" \
    "$HERE/run_overlap.sh")
echo "overlap analysis -> job $ajid (after collate)"

echo
echo "Submitted. Results dir: $RESULTS_DIR"
echo "When finished, see:"
echo "  $RESULTS_DIR/metrics.tsv       (time + memory per condition)"
echo "  $RESULTS_DIR/comparison.tsv    (orig vs refac identity per preFilt,callThresh)"
echo "  $RESULTS_DIR/benchmark.tsv     (tidy time/memory table)"
echo "  $RESULTS_DIR/overlap_summary.tsv / overlap_by_type.tsv / overlap_shift_histogram.tsv"
echo "  $RESULTS_DIR/peak_dumps/       (BEDs of lost/gained/shifted peaks)"
