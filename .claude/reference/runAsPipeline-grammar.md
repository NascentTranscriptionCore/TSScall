# smartSlurm `runAsPipeline` / `ssbatch` grammar cheat sheet

Self-contained reference for analyzing NTC bash driver scripts. Source of truth:
`/opt/ntc/smartSlurm` (`bin/runAsPipeline`, `bin/ssbatch`, `bin/jobRecordSchema.sh`,
`README.md`). Verified against code, not just docs.

> **Authoritative branch is now `master`.** The `ux-update` work was merged to `master`
> (`e1432b9`, 2026-10-22) and the `ux-update` branch deleted (NTC decision 2026-10 is
> satisfied — current `master`, verified at `52c849e`, *is* the former ux-update line).
> Author and validate against `master`; the O2 deployment at
> `/n/data1/cores/ntc/scripts/SmartSlurm` should track `master`. Anything in the NTC
> pipeline READMEs or older notes that says "use branch `ux-update`" is now stale → `master`
> (backlog NTC-008).

> **How to classify a script.** A script is a **runAsPipeline (rAP) pipeline** only
> if it contains `#@` decorator comment lines. `grep -n '#@' script.sh`. No `#@`
> lines ⇒ it is a **standard script** — read it by ordinary shell semantics, do NOT
> impose a DAG. Classification is **per script**, not per repo; one repo may hold both.
> The true rAP pipelines are **RNAseqMapping, pipeline2, DOBBE_F, and now ChIPseq**
> (its new `ChIPseq_rAP_pipeline_submission.sh`); GGA/MHT/peakMorph are hand-rolled
> `sbatch --dependency` DAGs, not rAP.

---

## 1. The `#@` job block

A *job block* = one `#@` marker line + the command line(s) directly beneath it, up to
its terminator (blank line or `#@end` — see §1.1).

```
#@ stepID , dependIDs , name , reference , inputs , sbatchOptions
<command(s) to run for this step>
```

Fields are comma-split (`IFS=','`). Order and meaning:

| # | Field | Meaning | Required | Example |
|--|-------|---------|:--:|---------|
| 0 | `stepID` | Unique integer id for the step. Duplicate ⇒ run aborts. | **Yes** | `1` |
| 1 | `dependIDs` | `0` = no dependency, else upstream stepIDs **dot-joined**. | **Yes** | `0`, `1`, `1.3` |
| 2 | `name` | Program label; groups `ssbatch` resource stats + names the job. | **Yes** | `bowtie2` |
| 3 | `reference` | Ref file(s)/dir(s), **dot-joined**; rsynced to node `/tmp` under `useTmp`. | No | `genome.fa`, `db1.db2` |
| 4 | `inputs` | Input file(s), **dot-joined**; size drives resource estimation. | No | `reads.fq`, `in1.in2` |
| 5 | `sbatchOptions` | Per-step sbatch string; must start with `sbatch`. Omit ⇒ CLI default. | No | `sbatch -c 4 -t 2:0:0 --mem 24G` |

Rules verified in `runAsPipeline`:
- **Only the first three fields are required**; keep empty commas as placeholders
  (`#@2,1,merge,,,`). Empty step/dep/name ⇒ hard error.
- `sbatchOptions`, if present, **must** begin with `sbatch` or the run aborts.
- **Dot** (`.`) is the separator *inside* the dependID, reference and inputs fields
  (e.g. depends on steps 1 and 3 = `1.3`; two refs = `db1.db2`).
- `name` may carry suffixes: `name:...` (`%%:*` strips trailing `:`), `.NoCheckpoint`
  disables checkpointing for that step; otherwise the global `checkpoint` tag applies.
- **Job arrays are NOT supported** — `--array=`/`-a` in a step's sbatch options ⇒ error.
- `ssbatch` must NOT appear in a command line inside the script ⇒ error (rAP adds it).

### 1.1 Where a block ends — two modes

`runAsPipeline` picks the mode **per script**:

**Explicit-closure mode (`#@end`).** If the script contains **any** `#@end` line, the WHOLE
script is parsed in explicit-closure mode:
- blank lines are treated as *body* (no longer terminators),
- only `#@end` closes a block,
- an unclosed `#@` at end-of-script is an error.
- In this mode each block body is written out as a **heredoc** to an on-disk step script
  `smartSlurmLog/jobs/<flag>/cmd.sh` (so multi-line bodies, `done`, and trailing commands
  are safe); literal `$` that must reach the node is escaped `\$` in the generated heredoc.
- `--lint` warns where a block relies on legacy blank-line closure.
- **Prefer `#@ … #@end`** for new/edited pipelines — it removes the blank-line foot-guns.
  ChIPseq's rAP driver uses this mode throughout.

**Legacy blank-line mode (no `#@end` anywhere in the script):**
- A `#@` block greedily collects **every following line** (joined with `;`) until a
  **blank/whitespace-only line**. Always leave a blank line before `done`, the next
  `#@`, or any following plain command, or they get swallowed into the job command.
- Full-line `# comment`: does **not** end the block; passed through to the converted
  script but not part of the job command.
- Trailing `code # comment`: the ` #…` tail is stripped **textually** (first space-`#`
  to EOL). Avoid a literal ` #` even inside quotes (`sed 's/ #/x/'` gets truncated). A
  `#` with no preceding space (`grep '#'`) is safe.
- **Heredocs break blocks** — the parser is line-based and flattens to one line; `<<EOF`
  then swallows the rest. Don't use heredocs (or `: <<EOF` comments) inside a `#@` block.

---

## 2. The two phases (when does a line run?)

```
PHASE 1 CONVERT (once, submit host):  runAsPipeline reads script top→bottom and writes
   smartSlurmLog/slurmPipeLine.<md5>.sh   (<md5> = md5sum of your script contents)
     • for / while / #loopStart  → remember loop var (becomes part of each job flag)
     • #@ marker + command        → emit an  ssbatch --wrap "command"  call
                                     (or, in #@end mode, a jobs/<flag>/cmd.sh heredoc)
     • any other line             → copied through unchanged
PHASE 2 SUBMIT (run the converted script):
     • plain lines  → run NOW, on the submit host
     • ssbatch call → SUBMIT a job (capture id); command runs LATER on a node once
                      -d afterok:<deps> is satisfied. Independent jobs beyond
                      firstBatchCount (5) are submitted HELD (-H), released after early
                      jobs produce records so later jobs can be resource-estimated.
```

Consequences to watch for when reading/modifying a pipeline:
- **Plain lines that inspect a step's output run before that output exists.** Any
  `[ -f out ]`, `$(samtools view -c out)` etc. in submit scope evaluates at submission,
  before the job runs. Move such logic into a *later* `#@` step (runs on a node).
- Unchanged script ⇒ reconversion skipped (md5 match); already-successful steps skipped
  unless rerun chosen (detected by `<flag>.success` in `smartSlurmLogDir`).
- Functions defined in plain scope are **not** visible inside `#@` blocks — `export -f`
  them first.
- **`/bin/sh` must be bash.** The generated script and step heredocs rely on bash features;
  a dash `/bin/sh` breaks the run.

### Two variable scopes
| Scope | What lives here | Evaluated | Visible to |
|-------|-----------------|-----------|-----------|
| **Submit** | plain lines + assignments outside any `#@` block | once, submit host (Phase 2) | later plain lines; the ref/inputs/sbatch fields + command text of every `#@` block (substituted at submit) |
| **Node** | code inside a `#@` block | later, on the compute node | only that one job |

- `$var` assigned in submit scope, used in a `#@` block ⇒ **baked in at submission**.
- `$var` assigned *and* used inside the same `#@` block ⇒ deferred, expands **on the node**
  (in `#@end` mode this is the `\$`-escaped form in `cmd.sh`).
- A var assigned inside a `#@` block exists nowhere else.
- **Steps pass data through files on shared storage, never through shell variables.**
  Convention: name each output path as a submit-scope var *immediately above* the `#@`
  that first uses it; producer writes it, consumer reads it, consumer declares the dep.
- **Submit-scope vars persist across steps, branches, loop iterations.** A missing
  assignment silently reuses a stale value — assign on *every* branch, and *inside*
  per-sample loops (not above them).

---

## 3. Loops

rAP preserves loop structure but extracts the `#@` blocks inside; the **loop variable
becomes part of each job's flag** (`1.0.findNumber.1`, `1.0.findNumber.2`, …).

Loop detection is **first-word** on the line (hardened on master): a line whose first token
is `for` opens a loop, `done` closes it, and `#loopStart:<var>` declares the var for a
`while`/`read` loop. A `for`/`done`/`#loopStart` appearing mid-line or mid-command is not
treated as loop structure — so a driver line must not *start* with `for` unless it is a real
loop (ChIPseq works around this by keeping manifest `for`-logic out of the driver, in
`prepare_manifest.sh`).

- **`for` loops**: the var right after `for` is auto-detected.
  ```bash
  for file in `ls dir`; do
      #@1,0,process,,file
      process.sh $file
  done
  ```
- **`while` loops** need a hint — declare the loop var on the line above:
  ```bash
  #loopStart:f1
  while read -r f1 f2 f3 f4; do
      #@1,0,process,,f1
      process.sh "$f1"
  done < samples.txt
  ```

---

## 4. `runAsPipeline` invocation

```
runAsPipeline --script "SCRIPT [ARGS]" --tmp {useTmp|noTmp} \
   [--sbatch-options "sbatch -p short -c 1 --mem 2G -t 50:0"] \
   [--mode dryrun] [--lint] [--pipeline NAME] [--email noEmail|noSuccEmail] \
   [--special checkpoint|excludeFailedNodes] [-v|--verbose | -q|--quiet] [-A|--account=ACCT]
```

| Option | Meaning | Default |
|--------|---------|---------|
| `--script "…"` | annotated script + args, one quoted string | required |
| `--tmp useTmp\|noTmp` | copy each step's `reference` to node-local `/tmp` or not | required |
| `--sbatch-options "…"` | default sbatch for steps with no own options; must start `sbatch` | `sbatch -p short -c 1 --mem 2G -t 50:0` |
| `--mode dryrun` | build + log plan, submit nothing | submit |
| `--lint` | parse + warn (e.g. legacy blank-line closure) without submitting | off |
| `--pipeline NAME` | label this run/pipeline (groups logs/records under a pipeline name) | derived |
| `--email noEmail\|noSuccEmail` | silence all / success emails | all |
| `--special checkpoint\|excludeFailedNodes` | checkpoint long jobs / avoid previously-failing nodes | none |
| `-v\|--verbose`, `-q\|--quiet` | output verbosity | normal |

- **Named flags only.** `runAsPipeline` parses named flags and aborts on any unknown token.
  The **legacy positional form** `runAsPipeline "SCRIPT ARGS" useTmp|noTmp run` in the
  pipeline READMEs is **not accepted** and is stale (backlog NTC-008).
- `--pipeline NAME` is new on master and not yet in the README — note it when documenting.

Outputs land in `smartSlurmLog/` (relative to cwd): per-step `.sh`/`jobs/<flag>/cmd.sh`
scripts, `.out` logs, `.success`/`.failed` flag files, `allJobs.txt`, and `.smartSlurm.log`.

### 4.1 Companion tooling (master)
- **`checkRun`** — browse step status/logs from the launch dir.
- **`checkRunGrid`** — grid/matrix status view across steps/samples (new).
- **`diagnose.sh`** — post-run diagnostics helper (new).
- **`cancelAllJobs`** — cancel this dir's running/pending jobs before a rerun.
- **`smartSlurmSetup.sh`** + `smartslurm activate|deactivate` — set up / toggle smartSlurm on
  `PATH` for a shell (the old `unExport` path is **deprecated**).
- **`jrRepair.sh`** — upgrade an old v1 `jobRecord.txt` in place to the v2 schema (§5.1).
- Retired: the plotting PNGs and `workflowPlot.py` (no longer part of the toolchain).

---

## 5. `ssbatch` (smart `sbatch`) — what each step becomes

```
ssbatch [-P PROGRAM] [-I INPUTS] [-F FLAG] [SBATCH_OPTIONS] --wrap="CMD"   [dryrun]
ssbatch [-P PROGRAM] [-I INPUTS] [-F FLAG] [SBATCH_OPTIONS] SCRIPT.sh [ARGS] [dryrun]
```

- `-P` program label (stats grouping); `-I` input file(s)/dir(s) or `jobSize:N`; `-F`
  unique flag. All three optional but `-P`/`-I` make estimation good.
- **Partition** auto-chosen from requested run-time via `adjustPartition` in
  `config.txt` (`partition{1,2,3}TimeLimit`); a bad `-p` is corrected for you.
- **OOM/OOT auto-resubmit**: `cleanUp.sh` records usage; a job that died out-of-memory/time
  is resubmitted using the **banded next-step** `nextMemAfterFail` / `nextTimeAfterFail`
  (step up to the next band / roughly double) with its formula cleared, up to a cap.
- **Checkpointing** (opt-in): snapshot + resume a long job before it hits its limit.
- **Emails**: `cleanUp.sh` emails the job script, exact submit command, and logs. Dampen with
  `noSuccEmail` / `noEmail`.
- Key `config.txt` knobs: `defaultMem=4096`M, `defaultTime=120`min, `defaultExtraMem`,
  `defaultExtraTime`, `firstBatchCount=5`.

> **`firstBatchCount` (5) ≠ the records-minimum rule.** A few records are the minimum before
> the size-fit model engages. `firstBatchCount` = how many independent rAP jobs run
> immediately before the rest are held pending early records.

### 5.1 jobRecords v2 (the resource history)
The schema moved from the old positional 18/21-column CSV to a **named v2 schema** defined in
`bin/jobRecordSchema.sh` — **17 fields**, written by column *name*, not position:

```
pipeline  step  ref  sample  jobID  inputSize  memUsed  memAllocated
timeUsed  timeAllocated  cores  status  user  node  epoch  isoDate  flag
```

- Default location `~/.smartSlurm/jobRecord.txt` (+ `stats/`). One row per completed job.
- **Reference-class pooling:** the `ref` key is normalized into a class (e.g. a mammalian
  genome → `mammal`) so records from like-sized references pool for estimation.
- **`jrRepair.sh`** upgrades a legacy v1 file to v2 in place; run it once on an old history.
- `sbatchBin/sbatch` (the bundled mock/passthrough) is **stale** — it still writes a v1
  header, so don't rely on it to produce v2 records; the real cluster path writes v2.

### 5.2 Resource estimation model (`estimateResource.sh` / `calculateMemTime.sh`)
Estimation keys on **program + reference-class**, reading `memUsed`/`timeUsed`, with
`inputSize` driving the size fit:
- **With input size:** a linear fit of size → mem and size → time, taken as `line × (1+frac)`
  evaluated at the **93rd percentile** (a headroom margin above the fitted line).
- **Slope ≤ 0:** falls back to a **constant model** (percentile of past usage, no size term).
- **Thin data** (too few records to fit): a **trend-scaled fallback** — scale a baseline by
  the limited trend available rather than fitting.
- **Failure floor:** estimates never drop below a floor derived from prior failures for that
  key, so a job that once OOM'd isn't re-sized too small.

---

## 6. Using `ssbatch` under other managers
`ssbatch` takes standard sbatch syntax, so Snakemake / Nextflow / Cromwell can use it as
their cluster submit command and inherit smart sizing (see smartSlurm `bin/Snakefile`,
`bin/nextflow.nf`, and README "Using ssbatch with other pipeline managers").

---

## 7. slurmless — the Slurm-free test harness (now in-tree)
slurmless is a Slurm-free regression harness for smartSlurm *itself* (mock `sbatch`/`sacct`,
fixtures, resource models; **~64 test blocks**), not wired to the NTC pipelines. It is **merged
in-tree at `smartSlurm/slurmless/`** (PR #144, 2026-11-11; the standalone repo
`ntc-agnt/slurmless_smartSlurm_testing` is now **deprecated**) — layout under `slurmless/`:
`bootstrap.sh`, `test/run_tests.sh`, `test/checkrun_ui_tests.sh`, `test/lib/harness.sh`,
`test/mockSlurm/`, fixtures + models, `skills/` (the harness, `smartslurm-adversary`, and
`migration-assistant` skills), and `agents/`. One copy, versioned with the engine, unable to drift.
Use it to prototype an rAP conversion or a resource-estimation change without a cluster.
**Co-versioning rule:** any smartSlurm engine change must ship with a paired slurmless harness
update and be run through `smartslurm-adversary` before on-cluster testing or an upstream push.
→ skills `slurmless-harness`, `smartslurm-adversary`.
