#!/usr/bin/env bash
# =============================================================================
# run_overlap.sh  --  SLURM wrapper that runs analyze_overlap.py on the results
# tree. Submitted by submit_all.sh after collate. Env: CONFIG (path to config.sh).
# Re-runnable by hand:  CONFIG=./config.sh bash run_overlap.sh
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CONFIG:-$HERE/config.sh}"
eval "$ENV_SETUP"

"$PYTHON" "$HERE/analyze_overlap.py" \
    --results-dir "$RESULTS_DIR" \
    --expname "$EXPNAME" \
    --script "$OVERLAP_SCRIPT" \
    --max-shift "$OVERLAP_MAX_SHIFT" \
    --dump-dir "$RESULTS_DIR/peak_dumps"
