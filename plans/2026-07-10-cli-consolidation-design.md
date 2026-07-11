# CLI Consolidation Design — `__main__.py` + `cli.py`

Date: 2026-07-10
Author: Nathan Mathieu (with Claude Code)
Status: design — awaiting review (not yet committed or implemented)

## Purpose

RelocaTE3 ships **two** argparse front-ends that have diverged:

- `src/RelocaTE3/__main__.py` — the registered console-script entry point
  (`relocaTE3 = "RelocaTE3.__main__:main"`). This is the path the real-rice
  validation harness runs and the target of recent TSD/genotype/excision parity
  work. **Behavioral source of truth.**
- `src/RelocaTE3/cli.py` — documented in `README.md`/`docs/source/`, imported by
  `tests/pipeline_test.py`, but never wired as the entry point and never
  validated. Drives a different (function-based) backend with different
  semantics.

This has already caused two null-result SLURM runs during parity work (wrong
path edited). Goal: **one CLI, one code path, no change to the installed,
validated `__main__.py` command behavior.** (Replacing the old `cli.py` parser
intentionally removes its alternate semantics — that is not a regression, since
that path was never wired or validated.)

## Consolidation direction (approved)

`__main__.py` remains the canonical CLI implementation. `cli.py` becomes a thin
compatibility import that re-exports `main` from `__main__`. The validated
handler bodies are **not moved** during this change.

```python
# src/RelocaTE3/cli.py  (after)
"""Compatibility shim: the canonical CLI now lives in RelocaTE3.__main__."""
from RelocaTE3.__main__ import main

__all__ = ["main"]
```

Rationale: moving the validated class-based handlers into `cli.py` is the single
change most likely to perturb real-rice output. Since "any change to validation
output is a regression," we keep the validated code exactly where it is and let
`cli.py` delegate. The `RelocaTE3.cli` import path keeps working for any external
consumer and for `tests/pipeline_test.py`.

**Compatibility scope.** `RelocaTE3.cli.main` remains importable and reaches the
canonical parser. This is *not* a promise that callers of the old `cli.py`
subcommands and flags (`--te`, `--sample`, `-i/--insertions`, the function-based
`find-reference`, the full-pipeline `run`, etc.) keep working — those alternate
semantics are intentionally removed. It is labeled a "compatibility shim" (not a
"deprecated shim") because no formal deprecation policy/timeline is being
introduced here; it simply preserves the import entry point.

## Scope boundaries (approved)

- **No new public pipeline semantics.** The consolidation does not add a
  full-pipeline CLI `run` or a `find-reference` command.
- `cli.py`'s full-pipeline `run` (via `run_sample`) and `find-reference`
  (steps 0/6 reference/shared calls) are **deferred**. Follow-up work items are
  filed (see "Follow-ups").
- **No compatibility aliases by default.** The installed `__main__.py` flags are
  preserved exactly. Tests and docs are updated to match those flags — not the
  other way around. An alias is added only if an existing repository consumer or
  supported workflow requires it and the semantics are strictly identical.

## What the validated interface is (must be preserved verbatim)

From `validation/real_rice/run_relocate3.sh`, the harness invokes:

| Command | Flags used by the harness | Backend | Meaning |
|---|---|---|---|
| `index-genome` | `-g` | `Aligner.index_genome` | samtools faidx + minimap2 index |
| `run` | `--left --right --te-library --name --outdir --threads --aligner --min-match --min-trimmed --mismatch` | `RelocaTE.identify_TE_reads` | **map + trim only** (TE-read ID + flank generation) — NOT the full pipeline |
| `align-genome` | `-g -f -n -o --threads` | `Aligner.map_genome_minimap` | trimmed flanks -> genome BAM |
| `find-insertions` | `-b --read-repeat --tsd --target --name --outdir --te-name --reference-ins --mismatch --min-mapq` | `InsertionFinder` class | cluster junctions -> non-ref calls |
| `characterize` | `-s -b -g -o --samtools --bcftools` | `Characterizer` class | zygosity (+ `-x` excision available) |

`map` and `annotate-ref` also exist in `__main__.py` (not used by the harness but
part of the installed surface). The harness performs reference handling via
`find-insertions --reference-ins`, not via any `find-reference` command.

## Changes

### 1. `src/RelocaTE3/cli.py` -> thin re-export shim
Replace the ~442-line parallel parser with the shim above. Removes the divergent
function-based subcommands from the installed/importable surface while keeping
`RelocaTE3.cli.main` importable and functional (it reaches the canonical parser).

### 2. `src/RelocaTE3/__main__.py` -> behavior-preserving only
No changes to handlers, flags, defaults, or dispatch behavior. **Help strings and
docstrings may be corrected** where they are inaccurate:
- The installed `run` subcommand currently advertises itself as the "complete
  RelocaTE3 pipeline" (`help="Run full pipeline: map + trim"`, description "Run
  the complete RelocaTE3 pipeline...", and `cmd_run.__doc__`). This is
  misleading. Reword all three to describe it explicitly as **TE-read
  identification and flank generation** (map reads to the TE library, then trim
  the TE-matching portion and emit flanking reads) — **not** the complete
  pipeline. This is a text-only change; parsing and behavior are untouched.

No aliases. (If, during doc/test updates, a specific existing consumer is found
to require an alias with strictly identical semantics, it will be added narrowly
and called out in review.)

### 3. Tests
- `tests/pipeline_test.py`
  - Migrate `test_trim_cli` (currently drives `cli.py`'s FASTQ `trim`) to the
    installed `run` command, which is the actual semantic equivalent (both call
    `identify_TE_reads`, both emit `flanking/<name>.{left,right}.flankingReads.fq`).
    Rename to reflect it (e.g. `test_run_cli_generates_flanking_reads`). Import
    `main` from `RelocaTE3.cli` (which now reaches the canonical parser) to also
    cover the shim path, or from `RelocaTE3.__main__` directly — decided in the
    implementation plan.
  - Keep `test_run_sample_produces_outputs` (library-level `run_sample`,
    unaffected by the CLI change).
  - Optional: add a test for the installed BAM-consuming `trim` command **only
    if a suitable BAM fixture already exists** under `tests/data/`. Not required
    to preserve current coverage.
- Shim coverage (imported-main): add a test that importing and invoking
  `RelocaTE3.cli.main` reaches the canonical parser successfully (e.g.
  `cli.main(["--version"])` returns/prints the version, or a known subcommand
  parses). Do **not** assert `cli.main is __main__.main` — that is an
  implementation detail.
- Entry-point smoke (subprocess): add lightweight subprocess tests that exercise
  the actually-installed entry points end to end. These **complement** (do not
  replace) the imported-main tests:
  - `relocaTE3 --version` (the console script)
  - `python -m RelocaTE3 --version` (the module entry point)
  Each should exit 0 and print the version. Guard the console-script test with a
  `shutil.which("relocaTE3")` skip so it does not fail in a bare checkout where
  the entry point is not installed.
- `tests/main_test.py` — unchanged (already drives `__main__`).

### 4. Documentation
Update `README.md`, `docs/source/usage.rst`, `docs/source/migration.rst`, and the
"two CLI files" note in `AGENTS.md` to describe the **installed staged
interface** accurately:
- Correct the subcommand table/flags to the `__main__.py` surface (add `map`,
  `index-genome`, `annotate-ref`; fix `trim` `-b` BAM input, `characterize`
  `-s/-b/-g/-x`, `find-insertions` `--tsd/--reference-ins`).
- Replace the misleading one-command full-pipeline `run` example with the real
  staged workflow (index-genome -> run -> align-genome -> find-insertions ->
  characterize).
- **State explicitly** that the present `run` performs TE-read identification and
  flank generation (map + trim), **not** the complete pipeline.
- Remove/relabel `find-reference` from user docs (it is not in the installed
  CLI); reference it only as planned/deferred.

## Validation protocol (the gate)

The harness skips already-completed outputs, so both runs use `--force` to
recompute from scratch.

1. **Baseline (before any edit):**
   - `pixi run test` — record pass/fail summary.
   - `pixi run validate-rice --local --force B_10` — capture the normalized
     non-reference call records, the comparison report (recall, precision,
     Jaccard), the genotype/status fields, and the list of expected output files.
   - **Copy the baseline normalized records and reports to a scratch directory
     OUTSIDE the validation output tree** (e.g. the session scratchpad), so the
     second `--force` run cannot overwrite them before comparison.
2. Implement changes 1–4.
3. **After:**
   - `pixi run test` — must be green.
   - `pixi run validate-rice --local --force B_10` — recompute.
   - **Compare sorted, normalized call records** (not raw file checksums) against
     the copied baseline: non-reference calls, recall, precision, Jaccard,
     genotype/status fields, and the set of expected output files. Any difference
     is a regression to fix or explicitly justify in the plan/PR.

## Follow-ups (filed, out of scope here)

Deferred capability:

- `todo/`: expose a full-pipeline CLI `run` (wraps `run_sample`) behind its own
  tests — new public semantics, deliberately not added in this consolidation.
- `todo/`: expose `find-reference` (steps 0/6 reference/shared insertion calls)
  as an installed, tested subcommand.

Behavioral defects discovered during review (do **not** fix in this
consolidation — each requires its own tested, reviewed change, and both are
relevant to the planned simulated-aligner benchmark):

- `todo/`: `cmd_run` parses `--min-match`, `--min-trimmed`, and `--mismatch` but
  does not forward them to `identify_TE_reads` — it calls
  `identify_TE_reads(reads, out, search_tool=aligner)`, silently dropping the
  three trim thresholds. Fixing changes trim output, so it must be validated
  separately.
- `todo/`: `run --aligner bwa` is accepted by argparse (`choices=["minimap2",
  "bwa"]`), but `identify_TE_reads` currently rejects non-minimap search tools.
  The advertised `bwa` path is therefore non-functional.

## Files touched

- `src/RelocaTE3/cli.py` (-> compatibility shim)
- `src/RelocaTE3/__main__.py` (**text only**: correct the `run` help/description
  and `cmd_run` docstring; no parsing/behavior change)
- `tests/pipeline_test.py` (migrate + rename `trim` test; add shim-reachability
  and subprocess `--version` smoke tests)
- `README.md`, `docs/source/usage.rst`, `docs/source/migration.rst`, `AGENTS.md`
- `todo/` (four follow-up items: two deferred capabilities, two behavioral defects)

**Not touched:** every validated `__main__.py` handler body, its flags, defaults,
and dispatch; the validation harness scripts; `tests/main_test.py`.

## Risks

- Removing `cli.py`'s subcommands from the importable surface could break an
  external caller that imported one of its private `_run_*`/`_add_*` helpers.
  Mitigation: those were never public API; the public `RelocaTE3.cli.main` is
  preserved.
- Doc rewrite could drift from actual flags. Mitigation: cross-check every
  documented flag against `__main__.py` and the harness before finishing.
