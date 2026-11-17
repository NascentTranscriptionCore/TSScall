# NTC review checklist (REVIEW.md)

Concise, checkable rules for reviewing NTC pipeline code (human or `/code-review`). Full rationale
+ per-language detail: the NTC coding-style spec (`.claude/ntc-coding-style.md`). Maintained by
the ecosystem overseer; stamped into each fork as `REVIEW.md`. Audience: **code read and
maintained by biologists** — clarity beats cleverness.

## Blocking (a violation is a finding)
- [ ] **Inputs validated** — required args checked; input files/assets verified to exist and be
      non-empty before use; bad input fails early with a clear message.
- [ ] **Errors are informative** — what failed + why + what to do, non-zero exit; no bare
      `exit 1`, no swallowed/`|| true`-hidden failures, no silent exit-0-empty.
- [ ] **Comments explain the non-obvious** — especially biology/statistics (thresholds, strand
      handling, column meanings) and any deviation from the obvious. No obvious-narration noise.
- [ ] **No silent wrong result** — a step that can't do its job doesn't report success.
- [ ] **Matches surrounding code** — naming, idioms, error style, comment density consistent with
      the file/repo.
- [ ] **Clever only when justified** — dense/non-obvious code carries a comment and a *significant*
      efficiency/clarity win; otherwise rewrite it plainly.

## Language specifics
- [ ] **bash:** `set -euo pipefail` (non-rAP-body scripts); `die()`/`fail()` with context; a
      `usage()` + flagged args; `: "${x:?msg}"` for required args; all expansions quoted;
      **`shellcheck` clean** (or disables justified inline).
- [ ] **R:** args parsed + validated (`length`, `file.exists`, `stop()` with message); `stop`/
      `message`/`warning` (not bare `print`); statistics commented; seed set where RNG matters;
      **`lintr`-checked**.
- [ ] **Python:** `argparse` with `help=`; validate + `sys.exit`/`raise` (not `assert`) on bad
      input; module/function docstring; `python -I` for untrusted data; **`ruff`/`flake8`-checked**.

## Checked before handoff (per `docs/update-workflow.md`)
- [ ] Parse/syntax: `bash -n` / `Rscript -e 'parse()'` / `python -m py_compile`.
- [ ] Linter run and addressed/justified.
- [ ] Dry-run invocation; for rAP drivers, the relevant `slurmless` block.

> A reviewer flags MUST violations by rule; the "clever only when justified" line is a judgment
> call weighed against the significant-gain bar, not an automatic reject.
