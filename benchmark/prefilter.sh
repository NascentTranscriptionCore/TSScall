#!/usr/bin/env bash
# =============================================================================
# prefilter.sh  --  produce a pre-filtered copy of the fwd/rev bedGraphs.
# Keeps a bedGraph interval only if its read count (column 4) is >= PREFILT.
# PREFILT=1 is a pass-through (no filtering) but still materializes the files
# so the downstream jobs have a uniform input path.
#
# Submitted by submit_all.sh as a SLURM job. Env vars expected:
#   CONFIG   : path to config.sh
#   PREFILT  : pre-filter cutoff
#   OUT_FWD  : destination forward bedGraph
#   OUT_REV  : destination reverse bedGraph
# =============================================================================
set -euo pipefail
source "$CONFIG"
eval "$ENV_SETUP"

filter_one() {
    local in="$1" out="$2"
    if [[ "$PREFILT" -le 1 ]]; then
        # No filtering: copy through unchanged (preserves any track lines).
        cp -f "$in" "$out"
    else
        # Keep header/track lines; keep data rows with column 4 >= PREFILT.
        awk -v k="$PREFILT" '
            /^track/ || /^browser/ || NF<4 { print; next }
            ($4+0) >= k              { print }
        ' "$in" > "$out"
    fi
}

echo "[prefilter] preFilt=$PREFILT"
echo "[prefilter] $FWD -> $OUT_FWD"
filter_one "$FWD" "$OUT_FWD"
echo "[prefilter] $REV -> $OUT_REV"
filter_one "$REV" "$OUT_REV"

echo "[prefilter] done. line counts:"
wc -l "$OUT_FWD" "$OUT_REV"
