# TSScall benchmark harness

Submits a matrix of SLURM jobs to compare the **original** and **refactored**
`TSScall.py` on real data, to characterize how time and peak memory scale with
the calling threshold and with input pre-filtering, and to analyze how the
called-peak universe changes across those levels.

Two independent read cutoffs are swept, named consistently everywhere:

- **preFilt** — pre-filter cutoff applied to the bedGraphs *before* TSScall runs.
  A position is kept if reads ≥ preFilt. `preFilt=1` means no pre-filtering.
- **callThresh** — value passed to TSScall's `--set_read_threshold` (the internal
  calling threshold), applied to whatever input survived pre-filtering.

It answers:

1. **Do the two scripts agree?** For each `(preFilt, callThresh)` cell, the
   original and refactored BED body and detail file are compared by md5.
2. **Time and memory across calling thresholds** — the `callThresh` axis.
3. **With vs. without pre-filtering, for both scripts** — the `preFilt` axis
   (`preFilt=1` is the "without" baseline).
4. **How the peak universe changes** — `analyze_overlap.py` reports exact /
   shifted / unique peaks across levels, split by TSS type.

## Files

| file | role |
|------|------|
| `config.sh`         | all paths, sweeps, SLURM resources, env setup. **Edit this.** |
| `submit_all.sh`     | builds pre-filtered inputs and submits the whole pipeline |
| `prefilter.sh`      | one SLURM job: strips bedGraph rows below `preFilt` |
| `run_one.sh`        | one SLURM job: runs a single condition, records time + MaxRSS |
| `collate.sh`        | one SLURM job: waits for complete metrics, then builds tables |
| `run_overlap.sh`    | one SLURM job: runs the overlap analysis after collate |
| `analyze_overlap.py`| peak-universe analytics (invoked by `run_overlap.sh`) |
| `equivalence/`      | standalone fuzz tests proving the rewritten hot paths match upstream |
| `RESULTS.md`        | distilled headline numbers from the validation run |

## Quick start

1. Edit `config.sh`:
   - `ORIG_SCRIPT`, `REFAC_SCRIPT`
   - `FWD`, `REV`, `CHROMSIZES`, `GTF`, `EXPNAME`
   - `ENV_SETUP` — your cluster's `module load` / `conda activate` lines.
     The refactored TSScall **and** `analyze_overlap.py` need numpy in that env.
   - `PARTITION`, `JOB_MEM`, `JOB_TIME`
   - optionally the sweeps `PREFILT_LEVELS`, `CALLTHRESH_LEVELS`, and the overlap
     match window `OVERLAP_MAX_SHIFT`
2. Preview without submitting: `./submit_all.sh --dry`
3. Submit: `./submit_all.sh`

Everything lands under `results_<EXPNAME>/`. Jobs are chained: each pre-filter
job runs first; the TSScall jobs that use its output wait on it (`afterok`);
`collate` runs as the matrix winds down and **waits until every condition has
written its metrics row** before building tables; `run_overlap` runs after
collate succeeds.

## The matrix

```
conditions = SCRIPTS x PREFILT_LEVELS x CALLTHRESH_LEVELS
```

Defaults: `{orig, refac} x {1, 4, 8} x {4, 6, 8, 10}` = 24 jobs, plus 3
pre-filter jobs, 1 collate job, and 1 overlap job. The fixed flags mirror the
production invocation: `--annotation_join_distance 500
--annotation_search_window 1000`, with `--annotation_file`, `--detail_file`,
`--set_read_threshold`, and the input/output paths supplied per job.

Note: when `callThresh < preFilt`, the pre-filter floor masks the calling
threshold (e.g. every `preFilt=8, callThresh<8` cell is identical), so only
`callThresh ≥ preFilt` cells are a clean threshold-matched comparison.

## Outputs

Under `results_<EXPNAME>/`:

- `metrics.tsv` — raw per-job row: script, preFilt, callThresh, return code,
  wall seconds, MaxRSS (`/usr/bin/time -v`), sacct MaxRSS/Elapsed cross-check,
  BED line count, BED/detail md5s.
- `comparison.tsv` — one row per `(preFilt, callThresh)`: `bed_identical`,
  `detail_identical` (YES/NO), plus orig/refac wall time and peak memory.
- `benchmark.tsv` — tidy time + peak-memory table (MaxRSS also in GB); the one
  to plot.
- `overlap_summary.tsv` / `overlap_by_type.tsv` / `overlap_shift_histogram.tsv`
  — peak-universe analytics (see below).
- `peak_dumps/<comparison>/` — `lost_from_reference.bed`, `gained_in_other.bed`,
  `shifted_reference_side.bed` for each preFilt comparison; load in IGV.
- `conditions/<script>_preFilt<PF>_callThresh<CT>/` — the BED, detail, and
  `time.txt` for each run.
- `diffs/preFilt<PF>_callThresh<CT>/` — full `bed.diff` / `detail.diff`, written
  **only** for cells that did not match (empty `diffs/` = everything agreed).
- `logs/` — per-job stdout/stderr.

## Peak-overlap analytics (`analyze_overlap.py`)

Runs automatically as the final pipeline step, and can be re-run by hand:

```
CONFIG=./config.sh bash run_overlap.sh
# or directly:
python analyze_overlap.py --results-dir results_<EXPNAME> --expname <EXPNAME> \
    --script refac --max-shift 100 --dump-dir results_<EXPNAME>/peak_dumps
```

For any pair of BED outputs it reports how many called TSSs are at identical
coordinates, how many are the same peak shifted by a few bp (with the shift
distribution), and how many are unique to one side — split by TSS type
(`obsTSS` = annotated, `uTSS` = unannotated / enhancer-class). Three comparison
families: **preFilt effect** (reference preFilt vs each other level at matched
callThresh — the read for "is no pre-filtering a reasonable standard?");
**callThresh effect** (adjacent thresholds at fixed preFilt); and **orig vs
refac** coordinate identity at every cell.

Key options: `--max-shift` sets how far apart two peaks can be and still count as
"the same, shifted" (the shift histogram always shows the true structure
regardless); `--ref-prefilt` picks the reference level (default: the lowest
present); `--dump-dir` writes the lost/gained/shifted BEDs. Needs numpy; no
bedtools required.

## Equivalence tests (`equivalence/`)

Two standalone fuzz tests prove the refactor's two rewritten hot paths return
exactly what the original code returns, on randomized inputs. They need only
Python and numpy — no SLURM, no cluster, no data files — and are the
machine-checkable core of the "byte-identical" claim.

```bash
cd equivalence
python binwinner_equiv.py     # bin_winner two-pointer vs. original nested loop
python filter_equiv.py        # grouped filter vs. original filterBedGraphListByWindows
```

Each prints a trial count and a mismatch count and ends in `PASS` (0
mismatches) or `FAIL`. `binwinner_equiv.py` also runs a fixed-tiling model — the
plausible-but-wrong refactor — and reports how often it diverges from the
original, as a control that the test can actually detect a bad implementation.
Expected: `bin_winner` 60,000 trials / 0 mismatches (tiling diverges on
~19,800); filter 40,000 trials / 0 mismatches. See
[`RESULTS.md`](RESULTS.md) for the recorded numbers.

## Notes

- **The collate timing guard.** `collate.sh` counts expected conditions from the
  config and waits (up to `COLLATE_MAX_WAIT` seconds) for `metrics.tsv` to be
  complete before building tables, so a slow straggler can never yield a partial
  comparison. If it ever times out it warns loudly and you can just re-run
  `CONFIG=./config.sh bash collate.sh` once the last job finishes.
- **Memory capture.** `run_one.sh` prefers `/usr/bin/time -v` for MaxRSS; if it
  is not on `PATH`, it falls back to SLURM `sacct`. Both are recorded.
- **Pre-filtering is lossy under `bin_winner`.** Expect `preFilt > 1` cells to
  differ from the `preFilt = 1` baseline even for the same script — that is the
  effect being measured. The orig-vs-refac comparison holds preFilt fixed.
- **Threads.** TSScall is single-threaded; jobs request `--cpus-per-task=1`.
- Re-running is safe: inputs are regenerated and each condition writes to its own
  directory. Delete `results_<EXPNAME>/` for a clean slate.
