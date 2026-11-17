---
name: tsscall
description: >
  Use when calling transcription start sites from single-nucleotide 5' coverage
  bedGraphs (PRO-seq / Start-seq) — annotated obsTSS plus novel/unannotated uTSS
  (enhancer-class) — or when a task in a related NTC repo touches TSScall, its
  bedGraph inputs, or its detail file. Covers running TSScall.py, its CLI/options,
  the detail-TSV column layout, utils, and the benchmark harness. Upstream are the
  1-nt strand-split bedGraphs from NTC PRO-seq mapping (pipeline2/RNAseqMapping) and
  the GTF + chrom.sizes from GeneAnnotationScripts/Genome_Onboarding; downstream are
  GeneAnnotationScripts (vendors AND runs it) and DOBBE_F (reads the detail file).
  NOTE: TSScall is NOT proTSScall (different tool, different mechanism).
  Aliases: TSScall, obsTSS, uTSS (ex-nuTSS).
---

# TSScall

## Purpose
TSScall calls **transcription start sites** from single-nucleotide 5' coverage
bedGraphs (PRO-seq, Start-seq). It reports **obsTSS** (observed starts at annotated
gene TSSs, requires a GTF) and **uTSS** (novel/unannotated starts, enhancer-class)
and pairs divergent/convergent bidirectional sites into clusters. Repo: `TSScall`
(branch `master`); an **NTC fork of lavenderca/TSScall** (Lavender/Burkholder,
NIEHS). Primary use here: understanding/integrating its I/O, not routinely running it.

## Architecture
**NOT rAP, NOT smartSlurm** — a plain **single-threaded Python program**, no
scheduler logic, no `ssbatch`, no flag folder. You wrap it in `sbatch` (or an rAP
step) yourself. Core is `TSScall.py` (~1180 lines, argparse). Pipeline inside one
process: read fwd/rev bedGraphs into numpy arrays keyed by (strand, chrom) → set
read threshold (fixed or FDR/Poisson) → call obsTSS in windows around GTF TSSs →
filter reads near obsTSS, then call uTSS → pair bidirectional + assign clusters →
write BED + optional detail TSV + optional cluster BED.

Supporting `utils/` (not run by core; need `TSScall.py` importable via PYTHONPATH):
`CLUSTERclassify.py`, `TSScheckForORIO.py`, `TSSclassify.pl` (Perl),
`TSScompare.py`, `TSSfilterToClusterRepresentatives.py`,
`bedgraph_oneBasedToZeroBased.py`. A `benchmark/` harness (raw `sbatch` bash +
Python) sweeps thresholds/pre-filter and proves refactor equivalence. Full option
table, detail columns, and harness: `reference/cli-and-detail-format.md`.

## I/O contract
**Inputs:** forward + reverse 5' bedGraphs (1-nt counts, **0-based**, 4-col);
tab-delimited `chrom.sizes`; optional GTF `--annotation_file` (**no obsTSS without
it** — uTSS-only otherwise). **Outputs:** BED6 with a track line (each TSS as a
1-nt feature `start-1 .. start`, TSS ID in name, strand); optional detail TSV
(16 columns + appended GTF attribute fields — see reference); optional
`--cluster_bed`.

Invocation (positional order: `forward_bg reverse_bg chrom_sizes output_bed`):
```
python TSScall.py \
  --annotation_file annotation.gtf --detail_file EXP_detail.txt \
  --set_read_threshold 8 --annotation_join_distance 500 \
  --annotation_search_window 1000 \
  forward.bedGraph reverse.bedGraph chrom.sizes EXP_output.bed
```
Key options (defaults in parens): `--set_read_threshold` (fixed threshold; else
FDR/Poisson), `--fdr`, `--call_method bin_winner|global` (bin_winner),
`--annotation_join_distance` (200; **prod uses 500**), `--annotation_search_window`
(1000), `--utss_filter_size` (750), `--bidirectional_threshold` (1000),
`--cluster_threshold` (1000). Full table: `reference/cli-and-detail-format.md`.

## Dependencies & environment
- Python 2.7 or 3 + **numpy>=1.20** (`requirements.txt`). No conda env pinned here.
- `utils/` need `TSScall.py` on PYTHONPATH; `TSSclassify.pl` needs Perl.
- Resources (refactored): ~**30 GB, 15–20 min** at full depth; the original OOMs
  near 90 GB. No module/partition baked in — you choose them in your sbatch wrapper.

## Modifying it safely (gotchas & assumptions)
- **Preserved legacy bugs are intentional** (kept for byte-identity with upstream):
  uTSS candidates after the last filter window on the lexically-last chromosome are
  dropped; the final uTSS search window is dropped (uTSS only). Do not "fix" these
  without re-running the equivalence harness.
- **Lexical chromosome sort** — chrom names must be consistent across bedGraph,
  chrom.sizes, and GTF or calls silently shift/drop.
- **bedGraphs must be 0-based** — convert 1-based input with
  `utils/bedgraph_oneBasedToZeroBased.py`.
- **Do NOT pre-filter low-count positions** before calling: it is lossy under
  `bin_winner` and changes the peak universe (contradicts old GGA README advice).
- No scheduler integration — single-threaded; you own the sbatch/resource wrapper.
- The copy **vendored inside GeneAnnotationScripts** (`TSScall-master/`) **may
  drift** from this repo (TODO: diff before assuming parity).

## Integration points
- **Upstream:** 1-nt strand-split 5' bedGraphs from NTC PRO-seq mapping
  (`RNAseqMapping`/pipeline2); GTF + chrom.sizes from `GeneAnnotationScripts` /
  `Genome_Onboarding`.
- **Downstream:** `GeneAnnotationScripts` (GGA) — vendors an **old copy** at
  `TSScall-master/` AND runs TSScall via sbatch (`--tsscallMin/Window/Time/Partition`),
  then builds dominant TSS/TES from the detail file + kallisto; `DOBBE_F`
  (`filter_and_annotate.sh` reads the detail file via `cut -f1,5-8`, intersects uTSS
  with dREG peaks for eRNA — **depends on detail-file column order**).
- **Shared formats:** uTSS/obsTSS naming (uTSS was formerly `nuTSS`; GGA docs still
  say `nuTSS`), BED6 + track header, the detail TSV.
- **NOT related:** `proTSScall` is a conceptual sibling only (standalone R,
  different mechanism) — not a fork or variant of TSScall.

## Reference files
- `reference/cli-and-detail-format.md` — full CLI option table, positional args,
  detail-TSV columns, benchmark harness, no-prefilter advice, evolution & gotchas.
- `.claude/reference/runAsPipeline-grammar.md` — smartSlurm/rAP grammar (TSScall is
  *not* rAP, but you often wrap it in an rAP step; useful ecosystem context).
