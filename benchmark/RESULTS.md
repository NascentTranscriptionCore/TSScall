# Benchmark results

Distilled headline numbers from validating the refactored `TSScall.py` against
upstream. This file is the permanent, small record; the raw outputs
(prefiltered bedGraphs, per-condition BEDs, `results_*/` trees) are not
committed. To regenerate, see [`README.md`](README.md).

Naming: `preFilt` = pre-filter cutoff applied to the input bedGraphs (1 = no
pre-filter); `callThresh` = value of `--set_read_threshold`.

## Correctness

Output is **byte-identical** to upstream (md5 of the BED body and of the detail
file) in every comparable condition.

- **Synthetic data:** identical across the production configuration,
  `global`, FDR-selected threshold, and no-annotation modes.
- **Full-depth real PRO-seq:** identical in all four conditions where the
  original completed (the largest with ~1.5 M called TSSs). The remaining two
  conditions have no comparison because the original was OOM-killed — see below.

### Unit equivalence (fuzz)

The two rewritten hot paths were checked against the original implementations on
randomized inputs (see [`equivalence/`](equivalence/)):

| Test | Trials | Mismatches vs. original |
|---|--:|--:|
| `bin_winner` two-pointer (`binwinner_equiv.py`) | 60,000 | **0** |
| grouped filter (`filter_equiv.py`) | 40,000 | **0** |

As a control, `binwinner_equiv.py` also runs a fixed-tiling model (the plausible
but incorrect refactor): it diverges from the original on **19,812** of the
60,000 trials, confirming the test distinguishes a wrong implementation from a
correct one. The filter test exercised an actual removal/tail-drop in 32,477 of
its 40,000 trials.

## Performance (full-depth real data)

Forward + reverse bedGraphs ~3.6 GB each. Original jobs were given a 90 GB
memory ceiling.

| preFilt | callThresh | original time | refactored time | speedup | original mem | refactored mem | identical |
|:--:|:--:|--:|--:|--:|--:|--:|:--:|
| 1 | 8 | OOM (>90 GB) | 21.1 min | — | — | 32.7 GB | orig cannot run |
| 1 | 12 | OOM (>90 GB) | 14.5 min | — | — | 27.1 GB | orig cannot run |
| 8 | 8 | 13.2 h | 2.6 min | 306× | 10.1 GB | 5.6 GB | yes |
| 8 | 12 | 4.1 h | 3.2 min | 77× | 8.5 GB | 4.0 GB | yes |
| 12 | 8 | 3.9 h | 2.2 min | 105× | 5.7 GB | 3.5 GB | yes |
| 12 | 12 | 2.7 h | 2.2 min | 74× | 5.4 GB | 3.2 GB | yes |

- **Full-depth, no pre-filter is only possible with the refactor.** The
  original is OOM-killed at 90 GB on the `preFilt 1` inputs; the refactor runs
  them in ~27–33 GB in 15–21 minutes.
- **Where the original completes, output is identical and the refactor is
  74–306× faster**, with the largest speedups on the densest windows.
- Peak memory is flat across thresholds (set by loading coverage), and the
  original-vs-refactored gap widens with depth, as expected from replacing the
  per-position dict expansion with columnar arrays.

## Pre-filtering effect

At a matched threshold, pre-filtering did not reproduce the full-data peak set —
it produced slightly **more** calls (a net change, not a subset), consistent
with `bin_winner` bin sums shifting when sub-threshold positions are removed:

| comparison | full-data peaks | pre-filtered peaks | change |
|---|--:|--:|--:|
| `callThresh 8`: preFilt 1 vs preFilt 8 | 1,494,847 | 1,513,829 | +1.3 % |
| `callThresh 12`: preFilt 1 vs preFilt 8 | 856,680 | 869,076 | +1.4 % |
| `callThresh 12`: preFilt 1 vs preFilt 12 | 856,680 | 873,020 | +1.9 % |

Pre-filtering also masks the calling threshold: at `preFilt 8`, a `callThresh`
below 8 is a no-op (`callThresh 4` and `callThresh 6` produced identical
output). Recommendation: run full-depth with no pre-filter and set sensitivity
with `--set_read_threshold`. See
[`../PERFORMANCE.md`](../PERFORMANCE.md#pre-filtering-is-unnecessary-and-lossy).
