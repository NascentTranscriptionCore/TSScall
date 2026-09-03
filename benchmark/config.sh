#!/usr/bin/env bash
# =============================================================================
# TSScall benchmark configuration.  EDIT THIS FILE, then run ./submit_all.sh
# =============================================================================
# Everything the harness needs is set here. All other scripts source this file.
#
# Two independent read cutoffs are swept, named consistently everywhere:
#   preFilt    -- pre-filter cutoff applied to the bedGraphs BEFORE TSScall.
#                 A position is kept if reads >= preFilt.  preFilt=1 means
#                 "no pre-filtering" (keep everything).
#   callThresh -- value passed to TSScall's --set_read_threshold (the internal
#                 calling threshold), applied to whatever input survived.

# ---- Scripts under test ------------------------------------------------------
# ORIG_SCRIPT  : the current/original TSScall.py
# REFAC_SCRIPT : the refactored TSScall.py
ORIG_SCRIPT="/path/to/original/TSScall.py"
REFAC_SCRIPT="/path/to/TSScall_refactored.py"

# ---- Input data (a single, lower-depth full PRO-seq dataset) -----------------
# FWD/REV are the UNFILTERED bedGraphs (all reads). The harness derives the
# pre-filtered versions itself so that "with vs without pre-filter" is a
# controlled axis rather than something baked into the inputs.
FWD="/path/to/proseq.forward.bedGraph"
REV="/path/to/proseq.reverse.bedGraph"
CHROMSIZES="/path/to/chrom.sizes"
GTF="/path/to/annotation.gtf"

# A short label used in output directory names.
EXPNAME="proseq_test"

# ---- Condition matrix --------------------------------------------------------
# SCRIPTS          : which scripts to run (keys map to the two paths above).
# PREFILT_LEVELS   : pre-filter cutoffs (see preFilt above). 1 = no pre-filter.
# CALLTHRESH_LEVELS: values passed to --set_read_threshold (see callThresh).
# The axes are independent: preFilt trims the input bedGraph on disk;
# callThresh is the cutoff TSScall applies internally to whatever survived.
SCRIPTS=(orig refac)
PREFILT_LEVELS=(1 4 8)
CALLTHRESH_LEVELS=(4 6 8 10)

# ---- Fixed TSScall flags (mirror the production invocation) ------------------
# --set_read_threshold and the input/output files are added per-job.
COMMON_FLAGS=(--annotation_join_distance 500 --annotation_search_window 1000)

# ---- Overlap analysis (analyze_overlap.py, run after collate) ---------------
# OVERLAP_SCRIPT    : whose BEDs to analyze (orig and refac are identical).
# OVERLAP_MAX_SHIFT : bp within which two peaks count as "same, shifted".
OVERLAP_SCRIPT="refac"
OVERLAP_MAX_SHIFT=100

# ---- SLURM resources ---------------------------------------------------------
PARTITION="norm"                 # sbatch --partition
JOB_MEM="120G"                   # sbatch --mem  (matches production ceiling)
JOB_TIME="24:00:00"              # sbatch --time (jobs should finish well under)
PREFILT_MEM="8G"
PREFILT_TIME="04:00:00"
CPUS_PER_TASK=1                  # TSScall is single-threaded

# collate waits for every condition's metrics row before building tables, so a
# slow straggler can never produce a partial comparison. This caps that wait.
COLLATE_MAX_WAIT=3600            # seconds (0 disables waiting)

# ---- Environment activation --------------------------------------------------
# These lines are executed at the top of every job, before python runs.
# The refactored TSScall and analyze_overlap.py need numpy on the path.
# Adapt to your cluster (module names, conda env, etc.).
ENV_SETUP='
module load miniconda/3 2>/dev/null || true
source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate your_env_with_numpy 2>/dev/null || true
'
# Python interpreter to invoke (after ENV_SETUP has run).
PYTHON="python"

# ---- Output layout -----------------------------------------------------------
# All artifacts land under RESULTS_DIR: prefiltered inputs, per-condition
# TSScall outputs, per-job metrics, logs, and the final collated tables.
RESULTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/results_${EXPNAME}"
