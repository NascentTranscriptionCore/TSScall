#!/usr/bin/env bash
# =============================================================================
# collate.sh  --  build the comparison + benchmark tables.
#
# THE TIMING FIX: this job WAITS until metrics.tsv has a row for every expected
# condition (or COLLATE_MAX_WAIT elapses) before building anything. That makes a
# slow straggler unable to produce a partial/empty comparison, independent of
# whatever scheduler dependency launched this job. Safe to re-run by hand at any
# time:  CONFIG=./config.sh bash collate.sh
#
# Produces, under RESULTS_DIR:
#   comparison.tsv  : orig-vs-refac output identity per (preFilt, callThresh)
#   benchmark.tsv   : tidy time + peak-memory table for every condition
#   diffs/          : full diffs for any (preFilt, callThresh) that did NOT match
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CONFIG:-$HERE/config.sh}"
eval "$ENV_SETUP"

METRICS="$RESULTS_DIR/metrics.tsv"
EXPECTED=$(( ${#SCRIPTS[@]} * ${#PREFILT_LEVELS[@]} * ${#CALLTHRESH_LEVELS[@]} ))

# ---- wait for a complete metrics table --------------------------------------
data_rows() { [[ -f "$METRICS" ]] && awk 'NR>1' "$METRICS" | grep -c . || echo 0; }
waited=0
while : ; do
    have=$(data_rows)
    if [[ "$have" -ge "$EXPECTED" ]]; then
        echo "[collate] metrics complete: $have/$EXPECTED rows"
        break
    fi
    if [[ "${COLLATE_MAX_WAIT:-0}" -le 0 || "$waited" -ge "${COLLATE_MAX_WAIT}" ]]; then
        echo "[collate] WARNING: metrics incomplete ($have/$EXPECTED) after ${waited}s;" \
             "building tables from what is present. Re-run collate.sh once the" \
             "remaining jobs finish." >&2
        break
    fi
    echo "[collate] waiting for metrics: $have/$EXPECTED (${waited}s elapsed)"
    sleep 30
    waited=$(( waited + 30 ))
done

mkdir -p "$RESULTS_DIR/diffs"

"$PYTHON" - "$RESULTS_DIR" "$EXPNAME" <<'PY'
import csv, os, sys, subprocess

results = sys.argv[1]
EXP = sys.argv[2]
metrics = os.path.join(results, "metrics.tsv")

# read metrics, de-duplicating by (script, preFilt, callThresh) keeping the
# last occurrence (so re-running a condition doesn't double-count).
seen = {}
with open(metrics) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        seen[(row["script"], row["preFilt"], row["callThresh"])] = row
rows = list(seen.values())

def to_gb(kb):
    try:
        return "%.2f" % (float(kb) / (1024.0 * 1024.0))
    except (ValueError, TypeError):
        return "NA"

# ---- benchmark.tsv : tidy time + memory --------------------------------------
bench = os.path.join(results, "benchmark.tsv")
with open(bench, "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["script", "preFilt", "callThresh", "rc",
                "wall_s", "maxrss_kb", "maxrss_gb", "bed_lines"])
    for row in sorted(rows, key=lambda d: (d["script"], int(d["preFilt"]),
                                           int(d["callThresh"]))):
        w.writerow([row["script"], row["preFilt"], row["callThresh"],
                    row["rc"], row["wall_s"], row["maxrss_kb"],
                    to_gb(row["maxrss_kb"]), row["bed_lines"]])
print("wrote", bench)

# ---- comparison.tsv : orig vs refac identity per (preFilt, callThresh) --------
by_key = {}
for row in rows:
    by_key.setdefault((row["preFilt"], row["callThresh"]), {})[row["script"]] = row

comp = os.path.join(results, "comparison.tsv")
n_mismatch = 0
with open(comp, "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["preFilt", "callThresh", "bed_identical",
                "detail_identical", "orig_rc", "refac_rc",
                "orig_bed_lines", "refac_bed_lines",
                "orig_wall_s", "refac_wall_s", "orig_maxrss_gb", "refac_maxrss_gb"])
    for (pf, ct), pair in sorted(by_key.items(),
                                 key=lambda kv: (int(kv[0][0]), int(kv[0][1]))):
        o = pair.get("orig"); r = pair.get("refac")
        if not o or not r:
            w.writerow([pf, ct, "MISSING", "MISSING",
                        o["rc"] if o else "NA", r["rc"] if r else "NA",
                        o["bed_lines"] if o else "NA",
                        r["bed_lines"] if r else "NA",
                        o["wall_s"] if o else "NA", r["wall_s"] if r else "NA",
                        to_gb(o["maxrss_kb"]) if o else "NA",
                        to_gb(r["maxrss_kb"]) if r else "NA"])
            continue
        # A YES certifies BOTH ran cleanly (rc==0) AND produced identical,
        # real output. Reject "NA" and the md5-of-empty-input so a failed or
        # body-less run can never be reported as a match.
        EMPTY_MD5 = "d41d8cd98f00b204e9800998ecf8427e"
        def _real(h):
            return h not in ("NA", "", EMPTY_MD5)
        rc_ok = (o.get("rc") == "0" and r.get("rc") == "0")
        bed_ok = rc_ok and _real(o["bed_md5"]) and o["bed_md5"] == r["bed_md5"]
        det_ok = rc_ok and _real(o["detail_md5"]) and o["detail_md5"] == r["detail_md5"]
        w.writerow([pf, ct, "YES" if bed_ok else "NO", "YES" if det_ok else "NO",
                    o["rc"], r["rc"], o["bed_lines"], r["bed_lines"],
                    o["wall_s"], r["wall_s"],
                    to_gb(o["maxrss_kb"]), to_gb(r["maxrss_kb"])])
        if not (bed_ok and det_ok):
            n_mismatch += 1
            cdir = os.path.join(results, "conditions")
            od = os.path.join(cdir, "orig_preFilt%s_callThresh%s" % (pf, ct))
            rd = os.path.join(cdir, "refac_preFilt%s_callThresh%s" % (pf, ct))
            dd = os.path.join(results, "diffs", "preFilt%s_callThresh%s" % (pf, ct))
            os.makedirs(dd, exist_ok=True)
            ob = os.path.join(od, "%s_output.bed" % EXP)
            rb = os.path.join(rd, "%s_output.bed" % EXP)
            with open(os.path.join(dd, "bed.diff"), "w") as out:
                subprocess.call(
                    "diff <(tail -n +2 %s) <(tail -n +2 %s)" % (ob, rb),
                    shell=True, executable="/bin/bash", stdout=out)
            subprocess.call(
                ["bash", "-c",
                 "diff %s/%s_detail.txt %s/%s_detail.txt > %s/detail.diff" %
                 (od, EXP, rd, EXP, dd)])
print("wrote", comp)
print("mismatched (preFilt,callThresh) cells:", n_mismatch)
PY

echo
echo "===================== comparison.tsv ====================="
column -t -s$'\t' "$RESULTS_DIR/comparison.tsv" 2>/dev/null || cat "$RESULTS_DIR/comparison.tsv"
echo
echo "====================== benchmark.tsv ====================="
column -t -s$'\t' "$RESULTS_DIR/benchmark.tsv" 2>/dev/null || cat "$RESULTS_DIR/benchmark.tsv"
