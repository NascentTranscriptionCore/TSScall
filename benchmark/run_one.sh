#!/usr/bin/env bash
# =============================================================================
# run_one.sh  --  run ONE TSScall condition and record wall time + peak memory.
# Submitted by submit_all.sh as a SLURM job. Env vars expected:
#   CONFIG       : path to config.sh
#   SCRIPT_TAG   : "orig" | "refac"
#   SCRIPT_PATH  : path to the TSScall.py to run
#   PREFILT      : pre-filter cutoff used to build the inputs
#   CALLTHRESH   : value for --set_read_threshold
#   IN_FWD/IN_REV: (pre-filtered) input bedGraphs
#   COND_DIR     : output directory for this condition
#   METRICS_FILE : shared TSV to append one metrics row to
# =============================================================================
set -euo pipefail
source "$CONFIG"
eval "$ENV_SETUP"

mkdir -p "$COND_DIR"
BED_OUT="$COND_DIR/${EXPNAME}_output.bed"
DETAIL_OUT="$COND_DIR/${EXPNAME}_detail.txt"
TIME_LOG="$COND_DIR/time.txt"

echo "[run_one] tag=$SCRIPT_TAG preFilt=$PREFILT callThresh=$CALLTHRESH"
echo "[run_one] script=$SCRIPT_PATH"
echo "[run_one] inputs: $IN_FWD  $IN_REV"

# Prefer GNU time -v for a portable MaxRSS; fall back to sacct afterward.
TIME_BIN=""
for cand in /usr/bin/time /bin/time; do
    if [[ -x "$cand" ]]; then TIME_BIN="$cand"; break; fi
done

CMD=( "$PYTHON" "$SCRIPT_PATH"
      --annotation_file "$GTF"
      --detail_file "$DETAIL_OUT"
      --set_read_threshold "$CALLTHRESH"
      "${COMMON_FLAGS[@]}"
      "$IN_FWD" "$IN_REV" "$CHROMSIZES" "$BED_OUT" )

start_epoch=$(date +%s)
set +e
if [[ -n "$TIME_BIN" ]]; then
    "$TIME_BIN" -v "${CMD[@]}" 2> "$TIME_LOG"
    rc=$?
else
    echo "[run_one] WARNING: /usr/bin/time not found; MaxRSS will rely on sacct" >&2
    "${CMD[@]}"
    rc=$?
fi
set -e
end_epoch=$(date +%s)
wall_s=$(( end_epoch - start_epoch ))

# --- parse peak memory (KB) from GNU time, if present ---
maxrss_kb=""
if [[ -f "$TIME_LOG" ]]; then
    maxrss_kb=$(awk -F': ' '/Maximum resident set size/ {print $2}' "$TIME_LOG" | tr -d ' ')
fi
# --- also record SLURM's own accounting as a cross-check ---
sacct_maxrss=""
sacct_elapsed=""
if command -v sacct >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
    # .batch step carries the MaxRSS for the payload
    read -r sacct_maxrss sacct_elapsed < <(
        sacct -n -P -j "${SLURM_JOB_ID}.batch" \
              --format=MaxRSS,Elapsed 2>/dev/null | head -1 | tr '|' ' ' )
fi

# checksum outputs. Missing or body-less output -> NA (never the md5 of an
# empty stream), so a failed run can't be mistaken for a real result. Skip the
# BED's line 1 (track header echoes the output filename).
if [[ -s "$BED_OUT" ]] && [[ "$(wc -l < "$BED_OUT")" -ge 2 ]]; then
    bed_md5=$(tail -n +2 "$BED_OUT" | md5sum | awk '{print $1}')
    bed_lines=$(wc -l < "$BED_OUT")
else
    bed_md5="NA"
    bed_lines=0
fi
if [[ -s "$DETAIL_OUT" ]]; then
    detail_md5=$(md5sum "$DETAIL_OUT" | awk '{print $1}')
else
    detail_md5="NA"
fi

# Append one row. A flock keeps concurrent jobs from interleaving writes.
{
    flock 9
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$SCRIPT_TAG" "$PREFILT" "$CALLTHRESH" "$rc" "$wall_s" \
        "${maxrss_kb:-NA}" "${sacct_maxrss:-NA}" "${sacct_elapsed:-NA}" \
        "$bed_lines" "${bed_md5:-NA}" "${detail_md5:-NA}" >> "$METRICS_FILE"
} 9>>"${METRICS_FILE}.lock"

echo "[run_one] rc=$rc wall=${wall_s}s maxrss_kb=${maxrss_kb:-NA} bed_lines=${bed_lines}"
exit "$rc"
