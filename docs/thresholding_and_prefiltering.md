# TSScall thresholding and pre-filtering — policy and evidence

**Audience:** anyone running TSScall on PRO-seq data (directly or via GGA), and
anyone deciding whether to pre-filter bedGraphs or how to set the calling
threshold.

**Scope:** this is about *noise handling and threshold choice* for PRO-seq TSS
calling. It does not change TSScall's output for a fixed input — it is a
recommendation about which inputs and parameters to use, backed by the analysis
summarized below.

---

## Recommendation (TL;DR)

For PRO-seq TSS calling with TSScall:

1. **Run on full-depth (raw) bedGraphs. Do not pre-filter** the bedGraphs by a
   minimum read count.
2. **Set the calling threshold explicitly** with `--set_read_threshold`
   (via GGA: `--tsscallMin`), chosen for library depth. Describe it honestly as
   a *transparent sensitivity cutoff*, not as an FDR-controlled statistic.
3. **Keep `--call_method bin_winner`** (the default).
4. Treat the residual edge-case differences (~2% of genes) as an inherent limit
   of TSS calling from this data type, not a bug to tune away.

One-line rationale: pre-filtering is a non-monotonic perturbation that changes
calls unpredictably; the global-Poisson FDR null is not valid for PRO-seq; and
no available TSScall knob reliably improves the hard cases. Full-depth input
plus an explicit threshold is the simplest defensible choice — and the evidence
below shows the alternatives do not help.

---

## Background: why the data type matters

TSScall was developed for **Start-seq**, where signal is sharp and confined to
TSSs and background is sparse and roughly spatially uniform. Under those
conditions a single genome-wide Poisson rate — the basis of TSScall's `--fdr`
threshold (after Jothi et al. 2008, *Nucleic Acids Res* 36:5221) — is a
reasonable null model for "what counts as signal above background."

**PRO-seq violates that assumption:**

- Signal spans gene bodies, pausing, termination, DoG readthrough, antisense and
  enhancers — most of the transcribed genome carries real, non-TSS signal.
- Coverage is spatially non-uniform and overdispersed; a single global rate is
  not a meaningful null.
- 5′ ends are not TSS-enriched the way Start-seq's are.

**Consequence:** a global Poisson FDR on PRO-seq is calibrated against the wrong
null. Because signal is pervasive, the whole transcribed genome looks "enriched"
relative to a uniform-random null, so FDR tends to select a *permissive*
threshold and admits diffuse gene-body / DoG positions. This is why `--fdr` is
not recommended for PRO-seq, and is consistent with the lab's existing practice
of setting `--set_read_threshold` explicitly.

---

## Evidence

All evidence is from one PRO-seq dataset (Armstrong005), run through GGA two
ways:

- **master** = original TSScall + **min-4 pre-filtered** bedGraphs
- **refac**  = refactored TSScall + **raw** bedGraphs

The two TSScall builds are byte-identical (md5-verified; see `PERFORMANCE.md`),
so **every difference below is caused by pre-filtering, not by code.**

### 1. Pre-filtering is a non-monotonic perturbation, not a noise knob

Pre-filtering removes low-count positions *before* bin selection, so it changes
which bin `bin_winner` picks — discontinuously. Sweeping the pre-filter level
(min0/raw → min4) at a fixed threshold (`--set_read_threshold 8`) on 8 hard
loci, scored by the distance (bp) of the best obsTSS call to the nearest
annotated TSS (lower = closer to annotation):

| gene        | min0 (raw) | min2 | min3 | min4 | direction |
|-------------|-----------:|-----:|-----:|-----:|-----------|
| EFCAB8      | **58**     | 641  | 641  | 641  | raw best; pre-filter destroys the correct call |
| CHD3        | **0**      | **0**| 550  | 550  | raw/min2 exact; pre-filter destroys it |
| NPIPB5      | 385        | 385  | **108** | 688 | non-monotonic (min3 best, min4 worst) |
| AC020929.1  | 202        | 143  | 202  | 143  | flip-flops min2↔min3 |
| AL035587.1  | 331        | **64** | 64 | 64   | pre-filter helps |
| SPON2       | 57         | 57   | 57   | 57→worse | neutral; min4 degrades one call |
| ARHGAP27    | 22         | 22   | 22   | 22   | neutral; pre-filter only adds spurious calls |
| ELF3        | 190        | 190  | 229  | 203  | mild hurt |

Across the 8: pre-filter **helps 2, hurts 3, neutral/erratic 3**. Two takeaways:

- **No pre-filter level is reliably best** — adjacent levels disagree
  (NPIPB5: 385 → 385 → 108 → 688), and raw is frequently the best.
- Pre-filtering **destroys correct raw calls** (EFCAB8 58→641, CHD3 0→550) about
  as often as it rescues wrong ones. That is not noise *reduction*.

Identical results at ±20 kb and ±50 kb windows — obsTSS calling is
window-size-invariant (the search window is ±1 kb around the annotated TSS), so
this is not a window-clipping artifact.

### 2. Call method: `global` does not help

`bin_winner` (highest read-*sum* bin) vs `global` (highest single position), on
raw data at a fixed threshold, 8 loci, scored by distance-to-annotation of the
best call: **6 ties, `bin_winner` closer in 2, `global` closer in 0.** Within
the ±1 kb obsTSS search window the competing positions all sit near the
annotation, and `bin_winner`'s cluster-support tracks the real start at least as
well as plain-max. Switching call methods buys nothing here.

### 3. FDR is the wrong null for PRO-seq

See **Background**. Use an explicit fixed threshold instead.

---

## Where the GGA "shifts" actually come from

In the Armstrong005 comparison, ~2.2% of genes differ between min4 and raw:
166 `tss_shift`, 97 `lost_from_baseline`, 21 `new_in_test`, 11
`tss_and_tes_shift`, 3 `tes_shift` (of 13,769 genes).

Decomposing **all 166 TSS shifts** using GGA's per-obsTSS read tables and
dominance flags (`*_Annotated_Dominant_and_Nondominant_obsTSS_fordREG.txt`):

- **97% (161/166) are reads-driven.** GGA's dominant-TSS rule picks, for most
  genes, the obsTSS with the most reads. Pre-filtering perturbs `bin_winner`'s
  obsTSS calls (lever 1 above), changing *which* obsTSS has the most reads, so
  GGA's dominant pick flips. GGA is a faithful messenger; the shift is the
  pre-filter artifact propagating through it.
- **3% (5/166) are kallisto-tiebreak.** When no obsTSS has a ≥60% read majority,
  GGA breaks the near-tie by proximity to the kallisto (RNA-seq TPM) dominant
  TSS, overriding reads. These are genuine dual-promoter near-ties (e.g.
  AC020929.1, ARHGAP27, SERPINF1) — inherently ambiguous biology, not a
  thresholding effect. (See the GGA companion note.)

So the shifts are overwhelmingly the pre-filter→`bin_winner` perturbation,
traced end-to-end. **Standardize on raw and 97% of these differences do not
arise** — there is no pre-filter to reshuffle the read ranking — and raw uses
full-depth reads, so the surviving call is the better-supported one.

---

## Scope and caveats

- Evidence is from one PRO-seq dataset. The 8-locus coverage analysis is a small
  illustrative sample; the 166-shift decomposition is the **full** shifted set.
- The annotation yardstick is the Ensembl GTF GGA used. A RefSeq cross-check is
  pending, but the non-monotonicity (§1) and the branch decomposition are
  yardstick-independent.
- Distance-to-annotation is a *proxy* for correctness; the closest-to-annotation
  call is not always biologically correct (alternative promoters, DoG regions).
- The affected genes are ~2% of the total, concentrated in low-expression and
  ambiguous loci. The recommendation deliberately does **not** over-optimize for
  them.

---

## Reproducing the analysis

The analysis scripts (prefilter sweep, `bin_winner`-vs-`global`, the 166-shift
branch decomposition) and the per-locus coverage windows used here are kept with
the investigation notes; each runs on Python + numpy and calls `TSScall.py`
directly on gene-window bedGraphs with the GGA parameters
(`--set_read_threshold 8 --annotation_search_window 1000
--annotation_join_distance 500 --call_method bin_winner --bin_winner_size 200`).
See `benchmark/` for the end-to-end identity harness.
