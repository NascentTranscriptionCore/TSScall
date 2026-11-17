# Performance refactor

This document describes the performance refactor of `TSScall.py` in this fork:
what changed, why the output is identical to upstream, the measured speed and
memory gains, and two upstream behaviors that were deliberately preserved.

## Summary

The refactor targets two independent scaling walls in the original: a runtime
bottleneck in per-window TSS calling, and a memory bottleneck in how coverage
is represented. Neither the calling logic nor the statistics were changed. On
every configuration tested — `bin_winner`, `global`, FDR-selected threshold,
and no-annotation — and on both synthetic and full-depth real PRO-seq data, the
refactored script produces **byte-identical BED and detail output** to upstream
(verified by md5).

The practical outcome: full-depth calling with no pre-filtering, previously
impossible on our nodes (OOM at 90 GB), now runs in ~30 GB and well under an
hour; and where the original completes at all, the refactor is tens to
hundreds of times faster with identical results.

## What changed

All changes are result-preserving. Each is a local edit; every method not
listed here is byte-identical to upstream.

### 1. Per-window calling: O(H²) → O(H)

The default `bin_winner` caller scores, for a window of `H` covered positions,
one overlapping bin anchored at each position (covering `[pos_i, pos_i +
bin_size]`), then returns the top position in the highest-scoring bin. The
original computes every bin's total with a nested scan that never terminates
early, so it is O(H²) per window and re-sums overlapping ranges repeatedly —
the dominant cost on dense windows.

The refactor computes the same bins with a **prefix sum** (each bin total
becomes one subtraction, O(1)) and a **two-pointer** sweep whose right edge
only advances (never resets), making the whole pass O(H). It reproduces the
original's overlapping, hit-anchored windows exactly — it does **not** tile the
axis — and preserves the original tie-breaking (first strict maximum in sorted
order, for both bin selection and the within-bin winner). `global` mode was
already O(H) and is unchanged.

### 2. Coverage representation: per-nucleotide dicts → columnar arrays

The original expands the input into one Python dict per single-nucleotide
position (`~340` bytes each), so a few hundred million covered positions
require 100+ GB. The refactor stores coverage as compact `numpy` `int64` arrays
of positions and reads, grouped by `(strand, chromosome)` — dropping the
redundant per-row `chromosome`, `strand`, and `end` fields — at roughly
`16–24` bytes per position. Grouping by `(strand, chromosome)` is equivalent to
the original global `(strand, chromosome, start)` sort, so every downstream
merge-walk (threshold scan, window intersection, filtering, uTSS seeding)
operates on the grouped arrays and yields identical results.

### 3. Vectorized threshold computation

`countLoci` and the total read count (used only on the FDR path) are computed
with vectorized `numpy` reductions instead of per-entry Python loops. Same
values, same selected threshold.

### 4. Linear window merges

Two window-merge routines drained a list with `list.pop(0)`, which is O(n) per
call and therefore O(n²) overall — a hidden quadratic that could dominate on
deep data. Both now use an index cursor (O(n)). The merge logic and output
order are unchanged.

## Correctness: how "identical" was verified

- **Unit fuzz tests.** The two rewritten hot paths were checked against the
  original implementations on randomized inputs: 60,000 trials for the
  `bin_winner` two-pointer and 40,000 trials for the grouped filter, **zero
  mismatches** in both. A fixed-tiling model (the tempting but wrong refactor)
  diverged on ~19,800 of the `bin_winner` trials, confirming the test detects
  the difference. See [`benchmark/equivalence/`](benchmark/equivalence/).
- **End-to-end, synthetic.** BED body and detail file are md5-identical to
  upstream across the production configuration, `global`, FDR, and
  no-annotation modes.
- **End-to-end, full-depth real data.** On a high-depth PRO-seq dataset, every
  condition where the original completed is md5-identical (BED and detail),
  including the largest (~1.5 M called TSSs). See
  [`benchmark/RESULTS.md`](benchmark/RESULTS.md).

## Benchmarks (full-depth real data)

Forward + reverse bedGraphs ~3.6 GB each. `preFilt` = pre-filter cutoff applied
to the input (1 = none); `callThresh` = `--set_read_threshold`. Original jobs
were given a 90 GB memory ceiling.

| preFilt | callThresh | original time | refactored time | speedup | original mem | refactored mem | identical |
|:--:|:--:|--:|--:|--:|--:|--:|:--:|
| 1 | 8 | OOM (>90 GB) | 21.1 min | — | — | 32.7 GB | orig cannot run |
| 1 | 12 | OOM (>90 GB) | 14.5 min | — | — | 27.1 GB | orig cannot run |
| 8 | 8 | 13.2 h | 2.6 min | 306× | 10.1 GB | 5.6 GB | yes |
| 8 | 12 | 4.1 h | 3.2 min | 77× | 8.5 GB | 4.0 GB | yes |
| 12 | 8 | 3.9 h | 2.2 min | 105× | 5.7 GB | 3.5 GB | yes |
| 12 | 12 | 2.7 h | 2.2 min | 74× | 5.4 GB | 3.2 GB | yes |

Two things to read from this. First, the original **cannot** perform
full-depth, no-pre-filter calling on these nodes — it is OOM-killed — while the
refactor does it in ~30 GB. Second, where the original does complete, output is
identical and the refactor is 74–306× faster, with the largest gains on the
densest windows (where the old O(H²) cost was worst).

Memory scales the way the mechanism predicts: peak usage is set by loading the
coverage, so it is flat across thresholds and dominated by the representation.
The dict-per-position expansion is exactly what the columnar arrays replace,
which is why the gap between original and refactored widens with depth — modest
on small pre-filtered inputs, decisive at full depth.

## Pre-filtering is unnecessary and lossy

Because the original's memory forced many users to strip low-count positions
from the bedGraphs before running, this fork measured what that pre-filtering
actually does to the calls. Two findings:

- **It perturbs the peak set.** Under `bin_winner`, removing sub-threshold
  positions changes the bin sums that pick winners, so pre-filtering shifts and
  fragments a small fraction of calls. At a matched threshold it produced
  ~1–2 % *more* peaks than full data on our dataset — a net change, not a clean
  subset. No pre-filtered result matched the full-data result at any threshold.
- **It silently overrides your threshold.** Pre-filtering at level `k` removes
  everything below `k`, so any `--set_read_threshold` below `k` is a no-op —
  e.g. pre-filtering at 8 makes `--set_read_threshold` of 4 and 6 produce
  identical output.

Since the refactor makes full-depth calling feasible, the recommendation is to
**run on full bedGraphs with no pre-filtering** and control sensitivity with
`--set_read_threshold`. That gives the undistorted `bin_winner` result and
makes the threshold mean exactly what it says. The peak-overlap analyzer used
for this analysis is in [`benchmark/`](benchmark/) (`analyze_overlap.py`).

## Two terminal-region bugs: preserved through the refactor, then fixed

The performance changes above were kept byte-identical to upstream, which meant
**reproducing two pre-existing bugs of the original that silently drop uTSS
signal.** They were then **fixed in a separate, clearly-labeled commit** (search
the code for `BEHAVIOR CHANGE`). This is the only intended difference from the
original's output.

Scope of both: **uTSS only** (annotated obsTSS calls are unaffected — the
annotated pass uses a different window builder and the unfiltered coverage), and
both are **end-of-sort-order effects**, confined to the terminal region of the
last lexically-sorted chromosome (minus strand, plus any scaffolds that sort
after it). They are disjoint, sitting on opposite sides of the window masked
around the last known TSS.

1. **Filter tail-drop** (`filterBedGraphListByWindows`). The original merge-walk
   terminated as soon as the filter windows were exhausted, discarding every
   coverage position ordered *after* the last filter window — so uTSS candidates
   past the last annotated/called TSS on the terminal chromosome were never
   seeded. **Fixed:** the vectorized filter masks and keeps every group (no
   truncation).

2. **Dropped final uTSS window**
   (`createUnannotatedSearchWindowsFromBedgraph`). The window-merge loop did not
   append its final accumulated window, dropping the single last (highest-sorted)
   uTSS search window. **Fixed:** the final window is now appended.

### What the fix changes (and what it does not)

The fix is a **strict, confined improvement**: it recovers the terminal uTSS the
original silently dropped, appended with the highest uTSS IDs, and wires their
bidirectional partners. **Every other call — ID, position, reads, type — is
byte-identical.** Verified by diffing against the byte-identical build: all
changes fall in the last chromosome's terminal cluster (minus strand), the only
existing-call changes are divergent/convergent partner fields pointing at the
recovered neighbors, and the per-cluster member count updates accordingly. On
the synthetic sets here the effect was +1 to +2 uTSS per run. The differential
fuzz test (`benchmark/equivalence/filter_equiv.py`) shows the removed tail-drop
would have changed the filter output in ~74% of random cases, so the footprint
is real, just confined to the sort tail.

**For the original-vs-improved comparison:** check out the byte-identical commit
to compare performance against the pristine original with no output confound
(`comparison.tsv` all YES/YES); check out this commit to adopt the corrected
calls and quantify the recovered uTSS with `benchmark/analyze_overlap.py`.

## Reproducing the validation

The [`benchmark/`](benchmark/) harness submits the full comparison matrix
(original vs. refactored × pre-filter levels × thresholds) to SLURM, records
wall time and peak memory, verifies output identity by md5, and runs the
peak-overlap analysis. Start from [`benchmark/README.md`](benchmark/README.md).
The equivalence fuzz tests in
[`benchmark/equivalence/`](benchmark/equivalence/) run standalone (Python +
numpy) and reproduce the zero-mismatch results above.
