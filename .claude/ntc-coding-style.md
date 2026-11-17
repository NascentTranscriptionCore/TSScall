# NTC coding style — specification

**Maintained by the ecosystem overseer (EO); pushed to agents via the per-repo overlay.**
Source of truth = this file in `ntc-ecosystem-context`; it is stamped into each fork as
`.claude/ntc-coding-style.md` and referenced from the repo's CLAUDE.md. Edit it here and
re-propagate — do not hand-edit the copies. Companion: `docs/ntc-review-checklist.md` (the
concise, citable rules; stamped as each fork's `REVIEW.md`).

## Who this code is for (the governing principle)
NTC pipelines are **maintained and used primarily by biologists without formal programming or
statistical training.** Therefore:

> **Prefer clear, simple, well-commented code over clever solutions.** Cleverness is justified
> *only* when it yields a **significant** gain in efficiency or clarity — and when it does, a
> comment must explain the trade-off. A dense one-liner that saves three readable lines is not
> justified; one that turns an O(n²) scan into an O(n) pass, or replaces 40 tangled lines with
> one well-named idiom, may be.

When in doubt, write the version a careful biologist could read in six months and modify safely.

## Hard requirements (MUST — these are review-blocking)
1. **Validate inputs.** Required arguments, input files, and reference assets are checked before
   work begins; a missing/empty/malformed input fails **early** with a clear message, never a
   cryptic downstream crash or a silent wrong result.
2. **Informative errors.** On failure, say *what* failed, *why*, and *what to do* — with a
   non-zero exit. No bare `exit 1`, no swallowed errors, no "command not found" as the only clue.
3. **Helpful comments.** Comment the *why* and anything non-obvious — especially the biology/
   statistics (what a threshold means, why a strand is flipped, what a column is), and any
   deviation from the obvious. Don't narrate the obvious (`i = i + 1  # add one`).
4. **No silent failures.** A step that can't do its job must not report success. Guard the cases
   that would otherwise "succeed" with empty/zero output (the exit-0-empty trap).
5. **Match the surrounding code.** New code reads like the file it lives in — same naming,
   comment density, idioms, error style. Consistency beats personal preference.

## Run the linter (strongly recommended, pre-handoff)
Agents should run the standard linter for the language as part of the static self-check before
handoff (see `docs/update-workflow.md`), and address or justify its findings:
- **bash →** `shellcheck` (and `bash -n` for a parse check).
- **R →** `lintr` (and `Rscript -e 'parse(...)'`).
- **Python →** `ruff` (or `flake8`) + run with `python -I` when reading untrusted/downloaded data.
Linters are advisory, not a CI gate (for now) — but "shellcheck is clean" is the expected bar for
new/edited bash. Note genuine false positives inline (`# shellcheck disable=SCxxxx — reason`)
rather than leaving them unexplained.

## Bash (the dominant language — drivers, workers, utilities)
House idioms (already standard across the repos; keep them):
- **Strict mode** at the top: `set -euo pipefail` (confirmed across DOBBE_F, Genome_Onboarding,
  RNAseqMapping, MHT, …). rAP `#@`-driver bodies are the one exception — follow the
  `smartslurm-pipelines` skill there (blank-line/`#@end` closure, no heredocs in blocks).
- **A `die()`/`fail()` helper** for errors: `die(){ echo "ERROR: $*" >&2; exit 1; }` — and use it
  with context (`die "genome FASTA not found: $fa"`), never a bare `exit 1`.
- **A `usage()` function** and real argument parsing (long flags preferred; see NTC-016's move
  from positional to flagged args).
- **Validate required args** with the house idiom `: "${manifest:?-m|--manifest is required}"`
  (confirmed widespread), and **check files/assets exist** before use
  (`[ -s "$fa" ] || die "empty/missing FASTA: $fa"`).
- **Quote all expansions** (`"$var"`, `"${arr[@]}"`); shellcheck enforces this.
- **Avoid clever one-liners** that obscure intent — a readable `for`/`if` beats a nested
  `awk`/`sed`/parameter-expansion puzzle unless the puzzle is a real, commented win.
- Paths that must survive into jobs/other nodes: **absolute**, and remember files are the
  inter-step currency (see the smartSlurm skill).

## R (analysis + plotting: DESeq2, stats, figures)
- **Parse args explicitly** (`commandArgs(trailingOnly=TRUE)` is the house pattern) and
  **validate them**: check `length(args)`, `file.exists()` each input, and `stop()` with a clear
  message if wrong — don't index `args[[3]]` and let it error cryptically. Consider `stopifnot()`
  or a small `need(cond, msg)` helper for preconditions.
- **`stop()` for errors, `message()`/`warning()` for progress/caveats** (not bare `print`), so
  failures are loud and logs are informative.
- **Comment the statistics.** State the model, the contrast, the filter thresholds, and what each
  output column means — the reader may not know DESeq2.
- Prefer **readable, named steps** over deeply nested tidyverse/`data.table` chains; if a dense
  chain is genuinely clearer or materially faster, comment what it does.
- Set a seed where results depend on RNG, and say so.

## Python (parsers, converters: gtf_parse.py, TSScall.py)
- **`argparse`** with helpful `help=` per argument (the house pattern) and sensible defaults;
  **validate** (files exist, values in range) and `sys.exit` with a message (or `raise` a clear
  exception) on bad input — don't `assert` for user-facing validation (asserts vanish under `-O`).
- A module/function **docstring** saying what it does, its inputs, and its outputs.
- Favor clear, straightforward code over clever comprehensions/generators when the clever form is
  harder to read; a comment earns a dense expression.
- Read untrusted/downloaded data with `python -I`; don't run interpreters from inside a data dir.

## Comments & errors — the standard, concretely
- **Good comment:** `# minus-strand coverage is stored positive; flip sign for the browser track`
- **Weak comment:** `# loop over samples` (the loop already says that).
- **Good error:** `die "no bowtie1 index for '$g' at $indexdir/$g — build it or fix --genome"`
- **Weak error:** `exit 1` / `echo "failed"`.

## How this is enforced
1. **Agents** load `.claude/ntc-coding-style.md` (referenced from CLAUDE.md) and self-check +
   lint before handoff (`docs/update-workflow.md`).
2. **Review** (`/code-review` / a human) cites `REVIEW.md` (the checklist) — a MUST violation is a
   finding; a cleverness call is weighed against the "significant gain" bar.
3. **EO** keeps this file current as conventions evolve and re-propagates on refresh.

*Grounded in a 2026-11 survey of the NTC repos; update the examples as the code evolves.*
