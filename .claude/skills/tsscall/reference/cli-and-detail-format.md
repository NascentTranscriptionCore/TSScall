# TSScall — CLI options, detail-file format, benchmark harness

`TSScall.py` is a single standard-library + numpy Python program (argparse `main`
at the bottom). No scheduler, no flag folder — wrap it in `sbatch`/rAP yourself.

## Positional arguments (order matters)
| # | Arg | Meaning |
|--|--|--|
| 1 | `forward_bedgraph` | forward-strand 5' bedGraph (1-nt counts, **0-based**, 4-col) |
| 2 | `reverse_bedgraph` | reverse-strand 5' bedGraph (same format) |
| 3 | `chrom_sizes` | tab-delimited `chrom <TAB> size` |
| 4 | `output_bed` | output TSS BED6 (written with a `track` line) |

## Optional flags (with code defaults)
| Flag | Default | Effect |
|--|--|--|
| `--annotation_file`, `-a` | none | GTF for obsTSS calling. **Omit → uTSS-only (no obsTSS).** |
| `--detail_file` | none | write the tab-delimited detail TSV (see columns below) |
| `--cluster_bed` | none | write clusters to a separate BED |
| `--set_read_threshold` FLOAT | none | fixed calling threshold; **skips** data-driven estimation |
| `--fdr` FLOAT | none | set threshold by FDR (default method targets < 0.001) |
| `--false_positives` INT | none | set threshold by an allowed false-positive count |
| `--call_method` {`bin_winner`,`global`} | `bin_winner` | per-bin winner vs. global max |
| `--bin_winner_size` INT | `200` | bin size for `bin_winner` |
| `--annotation_join_distance` INT | `200` | join adjacent annotation search windows (**prod = 500**) |
| `--annotation_search_window` INT | `1000` | obsTSS search window around each annotated TSS |
| `--utss_filter_size` INT | `750` | drop reads within this of an obsTSS/annoTSS before uTSS calling |
| `--utss_search_window` INT | `250` | uTSS search window |
| `--bidirectional_threshold` INT | `1000` | max distance to pair divergent/convergent bidirectional TSSs |
| `--cluster_threshold` INT | `1000` | max distance to group TSSs into a cluster |

Threshold precedence: an explicit `--set_read_threshold` wins; otherwise the model
estimates one (`--fdr` / `--false_positives` tune it). `requirements.txt`: `numpy>=1.20`
(Python 2.7 or 3). Refactored run ≈ 30 GB / 15–20 min full depth; original ≈ 90 GB (OOMs).

## Output BED6
One line per called TSS, each a **1-nt feature** spanning `start-1 .. start`
(0-based half-open), the TSS ID in the name column, score, and strand. A `track`
header line is emitted first.

## Detail file (TSV) column layout
Header is 16 fixed columns, then **one extra column per GTF attribute field** found
in the annotation (appended after column 16). Order is load-bearing — `DOBBE_F`
slices it with `cut -f1,5-8`.

| # | Column | Notes |
|--|--|--|
| 1 | `TSS ID` | unique id (used in the BED name too) |
| 2 | `Type` | `obsTSS` or `uTSS` |
| 3 | `Transcripts` | `;`-joined transcript_id(s), else `NA` |
| 4 | `Gene ID` | `;`-joined gene_id(s), else `NA` |
| 5 | `Strand` | |
| 6 | `Chromosome` | |
| 7 | `Position` | start coordinate |
| 8 | `Reads` | read count at the call |
| 9 | `Divergent?` | `True`/`False` |
| 10 | `Divergent partner` | partner TSS ID or `NA` |
| 11 | `Divergent distance` | or `NA` |
| 12 | `Convergent?` | `True`/`False` |
| 13 | `Convergent partner` | or `NA` |
| 14 | `Convergent distance` | or `NA` |
| 15 | `TSS cluster` | cluster id |
| 16 | `TSSs in associated cluster` | member count |
| 17+ | GTF attribute fields | one column each; values `;`-joined, `None`→`None`, missing→`NA` |

`cut -f1,5-8` (DOBBE_F) therefore yields: TSS ID, Strand, Chromosome, Position,
Reads — enough to build a uTSS BED. Unobserved annotated TSSs (windows with no hit)
are also emitted when a GTF is supplied.

## utils/ (adjuncts, not run by core)
Need `TSScall.py` importable (set `PYTHONPATH`); `TSSclassify.pl` additionally needs Perl.
- `bedgraph_oneBasedToZeroBased.py` — convert 1-based → required 0-based bedGraphs.
- `CLUSTERclassify.py`, `TSSclassify.pl` — classify TSSs/clusters.
- `TSScheckForORIO.py`, `TSScompare.py`, `TSSfilterToClusterRepresentatives.py` — ORIO
  check, compare two TSS call sets, reduce clusters to representatives.

## benchmark/ harness
Raw-`sbatch` matrix (no rAP) comparing **original vs. refactored** `TSScall.py` on
real data. Two independent cutoffs are swept:
- **preFilt** — bedGraph pre-filter applied *before* TSScall (`prefilter.sh`; keep
  rows with col-4 reads ≥ preFilt; **`preFilt=1` = pass-through, no filtering**).
- **callThresh** — passed to `--set_read_threshold` after pre-filtering.

Files: `config.sh` (paths, sweeps, SLURM resources, `ENV_SETUP` — the one file you
edit), `submit_all.sh` (builds inputs, submits the pipeline), `prefilter.sh`,
`run_one.sh` (one condition; records wall time + MaxRSS), `collate.sh` (tables),
`run_overlap.sh` + `analyze_overlap.py` (peak-universe: exact/shifted/unique by TSS
type), `equivalence/` (standalone fuzz tests proving rewritten hot paths match
upstream), `RESULTS.md` (headline numbers). It answers: do the two scripts agree
(BED body + detail compared by md5 per cell), how time/memory scale, with-vs-without
pre-filter, and how the peak universe shifts.

## Why you should NOT pre-filter for production calls
The benchmark exposes pre-filter as an *analysis knob*, not a recommended step.
Stripping low-count bedGraph positions is **lossy under `bin_winner`**: it changes
which position wins each bin and therefore the called-peak universe. This contradicts
older GGA README advice to pre-filter — for faithful calls, feed unfiltered 0-based
bedGraphs and control sensitivity with `--set_read_threshold`/`--fdr` instead.

## Evolution, rationale & gotchas
- **Lineage:** NTC fork of `lavenderca/TSScall` (Andy Lavender / Burkholder, NIEHS);
  upstream calling logic from ~2016–2018 is stable and unchanged here.
- **Sep 2026 performance refactor** (Seth Goldman, PR #1): AI-assisted, NTC-supervised
  rewrite of the hot paths. Output is **byte-identical to upstream** (md5-verified via
  the `benchmark/equivalence/` tests and per-cell BED+detail md5). It cut peak memory
  from ~90 GB (OOM-prone) to ~30 GB and runs in ~15–20 min at full depth.
- **Preserved legacy bugs (deliberate, for byte-identity):** (1) uTSS candidates after
  the last filter window on the **lexically-last chromosome** are dropped; (2) the
  **final uTSS search window is dropped** (uTSS only). "Fixing" either changes output
  and breaks equivalence — re-run the harness if you ever touch calling logic.
- **Lexical chromosome sort:** internal ordering is lexical, so chrom-name conventions
  must match across bedGraph, chrom.sizes, and GTF, or calls shift/drop silently.
- **0-based bedGraph requirement:** 1-based input gives wrong coordinates; convert first.
- **Vendored drift:** GeneAnnotationScripts carries its own `TSScall-master/` copy and
  also runs TSScall; it predates this refactor and may differ. **TODO: diff the
  vendored copy against this repo** before assuming parity.
- Latest repo activity: 2026-09-11.
