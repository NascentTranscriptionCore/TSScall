# NTC context stub (generated — do not hand-edit)

These files are stamped into this repo by `make-repo-stub.sh` in the
`ntc-ecosystem-context` meta-repo, so a Claude Code session opened here acts as this
repo's deep-analysis agent *and* is aware of the wider NTC ecosystem:

- `.claude/skills/tsscall/` — this repo's own deep-analysis skill (architecture, I/O
  contract, gotchas; deeper tables under its `reference/`).
- `.claude/ntc-ecosystem.md` — the always-loaded ecosystem catalog (the other NTC tools,
  the data-flow DAG, the shared file formats). Imported by the repo's `CLAUDE.md`.
- `.claude/reference/runAsPipeline-grammar.md` — the smartSlurm grammar cheat sheet.
- `.claude/ntc-coding-style.md` — the NTC coding-style spec (clarity over cleverness; validate
  inputs; informative errors; comments). EO-maintained; follow it for any code you write here.
- `REVIEW.md` (repo root) — the concise review checklist derived from the style spec; cited by
  `/code-review` and human reviewers. EO-maintained.
- `.claude/ntc-todo.md` — this repo's slice of the ecosystem BACKLOG (a generated projection;
  advisory — surface to the user, act only on confirmation).
- `.claude/ntc-changelog.md` — **repo-owned** (NOT overwritten): append a line here when you
  make a notable revision, so the ecosystem manager can track it.
- `.claude/ntc-brief-*.md` — **ad-hoc task briefs from EO** (a specific job spec: what to change,
  why, how, acceptance test). Posted on this repo's **default branch**; a brief names which branch
  the actual *work* targets. `ntc-todo.md` indexes any open briefs.

## Comms channel — how notes pass between EO and this repo's agent
All EO ↔ agent comms live in `.claude/` on this repo's **default branch** (the one branch every
fresh session sees — so notes are findable without knowing a feature-branch name):
- **EO → you:** `ntc-ecosystem.md` (standing context) · `ntc-todo.md` (backlog projection + an
  index of open briefs) · `ntc-brief-<topic>.md` (ad-hoc job specs).
- **you → EO:** append to `ntc-changelog.md` (repo-owned; EO harvests it).
- **At session start, read** `ntc-ecosystem.md` + `ntc-todo.md` and `ls .claude/ntc-brief-*.md`
  on the default branch — even if you then switch to a feature branch to do the work.

Edit the originals in the meta-repo and regenerate (`make-repo-stub.sh` /
`refresh-ntc-context`); changes here will be overwritten. Generated from
ntc-ecosystem-context @ 006b425 on 2026-11-12.
