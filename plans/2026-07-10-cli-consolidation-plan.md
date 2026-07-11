# CLI Consolidation Implementation Plan

> **SUPERSEDED (2026-07-11):** This plan made `__main__.py` canonical and `cli.py` a
> shim. The CLI has since been relocated — **`cli.py` is now canonical, `__main__.py`
> is a thin launcher, entry point `RelocaTE3.cli:main`** (see
> `plans/2026-07-10-cli-canonical-relocation-plan.md` and `AGENTS.md`). History only.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Collapse the two divergent argparse front-ends into one — `__main__.py` stays the canonical, validated CLI and `cli.py` becomes a compatibility shim — without changing any installed `__main__.py` command behavior.

**Architecture:** `src/RelocaTE3/cli.py` is replaced by a thin module that re-exports `main` from `RelocaTE3.__main__`. The validated handlers, flags, defaults, and dispatch in `__main__.py` are untouched except for inaccurate help/docstring text. Tests and docs are corrected to describe the installed staged interface. A before/after B_10 real-rice validation gate proves no behavioral regression.

**Tech Stack:** Python 3.10+, argparse, pytest, pixi (pins minimap2/samtools/bedtools), the `validate-rice` harness.

**Design doc:** `plans/2026-07-10-cli-consolidation-design.md` (approved).

**Note on location:** Per repo convention, plans live in `plans/` (not the skill's default `docs/plans/`), matching `plans/2026-06-19-characterize-validation-design.md`.

**Reference skills:** @superpowers:test-driven-development for each code task; @superpowers:executing-plans to drive the plan.

---

## Task 0: Capture the B_10 validation baseline (BEFORE any edit)

This must run on the current, unmodified tree so the baseline reflects validated
behavior. Nothing here modifies the repo.

**Files:** none modified. Baseline copied to the session scratchpad.

**Step 1: Record the baseline test result**

Run: `pixi run test 2>&1 | tee /scratch/nmath020/claude-4628/-bigdata-stajichlab-nmath020-github-github-tools-RelocaTE-RelocaTE3-jason-RelocaTE3/0b8a6ff7-eae4-4efb-8e91-8cade8bcacff/scratchpad/baseline_pytest.txt`
Expected: suite passes (record the summary line, e.g. "N passed").

**Step 2: Resolve the harness output + report directories**

Run (from repo root):
```bash
cd validation/real_rice
CONFIG="$(ls config*.toml | grep -v example | head -1)"; [ -z "$CONFIG" ] && CONFIG="config.example.toml"
python3 -c "from _config import load_config; c=load_config('$CONFIG'); print(c['paths']['relocate3_outdir']); print(c['paths']['report_dir'])"
cd -
```
Expected: two absolute paths — `relocate3_outdir` (per-sample outputs) and `report_dir` (normalized `.tsv` + `summary.txt`). Note both.

**Step 3: Run the forced B_10 validation**

Run: `pixi run validate-rice --local --force B_10`
Expected: completes; writes normalized non-ref + characterized records and
`<report_dir>/nonref/summary.txt` and `<report_dir>/characterized/summary.txt`.
(`--force` is required because the harness skips already-completed outputs.)

**Step 4: Copy the baseline OUT of the validation tree**

The second (post-change) `--force` run will overwrite these, so copy first.
```bash
BASE=/scratch/nmath020/claude-4628/-bigdata-stajichlab-nmath020-github-github-tools-RelocaTE-RelocaTE3-jason-RelocaTE3/0b8a6ff7-eae4-4efb-8e91-8cade8bcacff/scratchpad/cli-consolidation-baseline
mkdir -p "$BASE"
cp -r "<report_dir>" "$BASE/report"                       # normalized records + summaries
cp -r "<relocate3_outdir>/B_10/results" "$BASE/B_10_results"   # per-sample call outputs
ls -R "$BASE" | head -40
```
Expected: baseline `report/` and `B_10_results/` present under `$BASE`.

**Step 5: Commit (no repo change — baseline is external; skip commit).**
This task produces no tracked changes. Proceed to Task 1.

---

## Task 1: Replace `cli.py` with a compatibility shim (TDD)

**Files:**
- Test: `tests/pipeline_test.py` (add one test)
- Modify: `src/RelocaTE3/cli.py` (replace entire file)

**Step 1: Write the failing shim-reachability test**

Add to `tests/pipeline_test.py` (keep the existing `from RelocaTE3.cli import main`):
```python
def test_cli_shim_reaches_canonical_parser():
    """RelocaTE3.cli.main must reach the canonical __main__ parser.

    `index-genome` exists only in the canonical parser, so this passes only
    once cli.py delegates to __main__. We assert reachability, not identity.
    """
    assert main(["index-genome", "--help"]) == 0
```

**Step 2: Run test to verify it fails**

Run: `pixi run pytest tests/pipeline_test.py::test_cli_shim_reaches_canonical_parser -v`
Expected: FAIL/ERROR — the old `cli.py` has no `index-genome` subcommand, so
`argparse` raises `SystemExit(2)` (its `parse_args` is outside the try/except).

**Step 3: Replace `cli.py` with the shim**

Overwrite `src/RelocaTE3/cli.py` with exactly:
```python
"""Compatibility shim for the RelocaTE3 command-line interface.

The canonical CLI now lives in :mod:`RelocaTE3.__main__`. This module is kept
so that ``from RelocaTE3.cli import main`` continues to work and reaches the
canonical parser.

Compatibility scope: the ``RelocaTE3.cli.main`` entry point is preserved. The
old ``cli.py`` subcommands and flags (e.g. ``--te``, ``--sample``,
``-i/--insertions``, the function-based ``find-reference``, and the
full-pipeline ``run``) are intentionally removed in favor of the single,
validated ``__main__`` interface.
"""

from __future__ import annotations

from RelocaTE3.__main__ import main

__all__ = ["main"]
```

**Step 4: Run test to verify it passes**

Run: `pixi run pytest tests/pipeline_test.py::test_cli_shim_reaches_canonical_parser -v`
Expected: PASS (returns 0 — `--help` triggers `SystemExit(0)`, which `__main__.main` swallows and returns 0).

**Step 5: Confirm no other importer broke**

Run: `grep -rn "from RelocaTE3.cli import\|RelocaTE3.cli\." --include="*.py" . | grep -v .pixi`
Expected: only `tests/pipeline_test.py`. If anything imports a removed symbol
(`build_parser`, `_run_*`, `_add_*`), STOP and report — the design assumed none exist.

**Step 6: Commit**

```bash
git add src/RelocaTE3/cli.py tests/pipeline_test.py
git commit -m "refactor(cli): make cli.py a compatibility shim over __main__"
```

---

## Task 2: Migrate the FASTQ `trim` test to the installed `run` (TDD)

The old `cli.py trim` (FASTQ in → `identify_TE_reads`) is semantically the
installed `run` (FASTQ in → `identify_TE_reads`). Migrate + rename.

**Files:**
- Modify: `tests/pipeline_test.py` (replace `test_trim_cli`)

**Step 1: Rewrite the test**

Replace the whole `test_trim_cli` function with:
```python
def test_run_cli_generates_flanking_reads(tmp_path: Path):
    """The installed `run` subcommand (map + trim) emits flanking FASTQs.

    This is the semantic equivalent of the old cli.py FASTQ `trim` command;
    both drive RelocaTE.identify_TE_reads. Exercised through RelocaTE3.cli.main
    to also cover the shim path.
    """
    rc = main(
        [
            "run",
            "-l",
            str(R1),
            "-r",
            str(R2),
            "-T",
            str(TELIB),
            "-n",
            "HEG4",
            "-o",
            str(tmp_path),
            "--threads",
            "2",
        ]
    )
    assert rc == 0
    assert (tmp_path / "flanking" / "HEG4.right.flankingReads.fq").exists()
```

**Step 2: Run the test**

Run: `pixi run pytest tests/pipeline_test.py::test_run_cli_generates_flanking_reads -v`
Expected: PASS. (`run` uses default `--aligner minimap2`; `identify_TE_reads`
writes `flanking/HEG4.{left,right}.flankingReads.fq`.)

**Step 3: Run the whole file to confirm nothing else regressed**

Run: `pixi run pytest tests/pipeline_test.py -v`
Expected: `test_run_sample_produces_outputs`, `test_run_cli_generates_flanking_reads`,
`test_cli_shim_reaches_canonical_parser` all PASS.

**Step 4: Commit**

```bash
git add tests/pipeline_test.py
git commit -m "test(cli): migrate FASTQ trim test to installed run command"
```

---

## Task 3: Add subprocess entry-point smoke tests

These complement (do not replace) the imported-`main` tests by exercising the
real installed entry points end to end.

**Files:**
- Modify: `tests/main_test.py` (add imports + two tests)

**Step 1: Add the tests**

Add at the top of `tests/main_test.py` (imports):
```python
import shutil
import subprocess
import sys
```
Add at the end of `tests/main_test.py`:
```python
def test_module_entry_point_version():
    """`python -m RelocaTE3 --version` runs and prints a version."""
    proc = subprocess.run(
        [sys.executable, "-m", "RelocaTE3", "--version"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert proc.stdout.strip()


@pytest.mark.skipif(
    shutil.which("relocaTE3") is None,
    reason="console script not installed on PATH",
)
def test_console_script_version():
    """The installed `relocaTE3 --version` console script runs."""
    proc = subprocess.run(
        ["relocaTE3", "--version"], capture_output=True, text=True
    )
    assert proc.returncode == 0
    assert proc.stdout.strip()
```

**Step 2: Run the new tests**

Run: `pixi run pytest tests/main_test.py -v`
Expected: existing tests PASS; `test_module_entry_point_version` PASS;
`test_console_script_version` PASS (or SKIP if the script isn't on PATH).

**Step 3: Commit**

```bash
git add tests/main_test.py
git commit -m "test(cli): add subprocess smoke tests for both entry points"
```

---

## Task 4: Correct inaccurate `run` help/description/docstring (text only)

No parsing or behavior change — only strings that call `run` the "complete
pipeline."

**Files:**
- Modify: `src/RelocaTE3/__main__.py` (three text sites)

**Step 1: Fix the subparser help + description**

In `main()` where the `run` subparser is added (~lines 651-658), change:
```python
    _menu_run(
        subparsers.add_parser(
            "run",
            formatter_class=CustomHelpFormatter,
            help="Run full pipeline: map + trim",
            description="Run the complete RelocaTE3 pipeline: align reads to the TE library then trim TE sequences.",
        )
    )
```
to:
```python
    _menu_run(
        subparsers.add_parser(
            "run",
            formatter_class=CustomHelpFormatter,
            help="Identify TE-containing reads and generate flanking reads (map + trim)",
            description="Map reads to the TE library then trim the TE-matching "
            "portion, emitting flanking reads and the read_repeat_name table. "
            "This is TE-read identification and flank generation (steps 2-3), "
            "NOT the complete insertion-calling pipeline.",
        )
    )
```

**Step 2: Fix the `cmd_run` docstring**

Change `cmd_run`'s docstring (~line 482) from:
```python
    """Run full pipeline: align reads to TE library then trim TE sequences."""
```
to:
```python
    """Identify TE-containing reads and generate flanking reads (map + trim).

    Maps reads to the TE library then trims the TE-matching portion. This is
    steps 2-3 (TE-read identification and flank generation), not the complete
    insertion-calling pipeline.
    """
```

**Step 3: Verify help text + behavior unchanged**

Run: `pixi run relocaTE3 run --help`
Expected: help now describes map + trim / flank generation; the flag list
(`-l/-r/-T/-n/-o/--threads/--aligner/--min-match/--min-trimmed/--mismatch`) is unchanged.

Run: `pixi run pytest tests/ -q`
Expected: all PASS (no behavior change).

**Step 4: Commit**

```bash
git add src/RelocaTE3/__main__.py
git commit -m "docs(cli): correct run help/docstring (map+trim, not full pipeline)"
```

---

## Task 5: Rewrite user docs to the installed staged interface

**Files:**
- Modify: `README.md`
- Modify: `docs/source/usage.rst`
- Modify: `docs/source/migration.rst`

**Step 1: README — subcommand table + flags**

- Replace the flagship one-command `relocaTE3 run ... --genome --repeatmasker --genotype`
  example with the real staged workflow:
  ```
  relocaTE3 index-genome -g reference.fa
  relocaTE3 run -l reads_1.fq.gz -r reads_2.fq.gz -T RiceTE.fa -n HEG4 -o HEG4_out --threads 8 --mismatch 2
  relocaTE3 align-genome -g reference.fa -f HEG4_out/flanking/HEG4.5.flankingReads.fq HEG4_out/flanking/HEG4.3.flankingReads.fq -n HEG4 -o HEG4_out --threads 8
  relocaTE3 find-insertions -b HEG4_out/genome_aln/HEG4.genome.bam --read-repeat HEG4_out/te_containing/HEG4.read_repeat_name.txt --tsd TTA -n HEG4 -o HEG4_out --reference-ins reference.fa.RepeatMasker.out
  relocaTE3 characterize -s HEG4_out/results/HEG4.all_nonref_insert.txt -b reads_to_genome.bam -g reference.fa -o HEG4_out/results
  ```
  (Cross-check every flag against `src/RelocaTE3/__main__.py` and
  `validation/real_rice/run_relocate3.sh` before finalizing exact filenames.)
- Update the "Individual steps" table to the installed subcommands: `map`,
  `index-genome`, `annotate-ref`, `trim` (BAM in), `run` (map+trim), `align-genome`,
  `find-insertions`, `characterize`. **State that `run` = TE-read identification +
  flank generation, not the full pipeline.**
- Remove `find-reference` from the installed table; if mentioned, mark "planned".
- Fix the Migration table rows whose flags don't match the installed CLI.

**Step 2: usage.rst — same corrections in reStructuredText**

- Replace the `find-reference` example (`docs/source/usage.rst:50`) — not an
  installed command.
- Correct `trim` (`-b` BAM), `characterize` (`-s/-b/-g`, `-x` excision),
  `find-insertions` (`--tsd`, `--reference-ins`) examples.
- Add the explicit note that `run` performs map + trim only.

**Step 3: migration.rst — reconcile the RelocaTE2→3 mapping**

- Ensure `characterizer.pl` maps to `relocaTE3 characterize` with the installed
  flags; drop any `run --genotype` claim that implies a one-command full pipeline.

**Step 4: Verify docs build**

Run: `pixi run docs`
Expected: Sphinx build succeeds with no new errors/warnings about the edited files.

**Step 5: Commit**

```bash
git add README.md docs/source/usage.rst docs/source/migration.rst
git commit -m "docs: describe the installed staged CLI accurately"
```

---

## Task 6: Update the AGENTS.md "two CLI files" warning

**Files:**
- Modify: `AGENTS.md` (the "Two CLI files — read this before editing" section, ~lines 49-56)

**Step 1: Replace the section**

Rewrite it to state there is now ONE canonical CLI (`__main__.py`), that
`cli.py` is a compatibility shim re-exporting `main`, and that new/renamed
subcommands only need to be added in `__main__.py` + the README table.

**Step 2: Commit**

```bash
git add AGENTS.md
git commit -m "docs(agents): single canonical CLI; cli.py is now a shim"
```

---

## Task 7: File the four follow-up items in `todo/`

**Files:**
- Create: `todo/cli-full-pipeline-run.md`
- Create: `todo/cli-find-reference.md`
- Create: `todo/cmd-run-drops-trim-thresholds.md`
- Create: `todo/run-aligner-bwa-unsupported.md`
- Modify: `todo/TODO_REGISTRY.md` (add four rows)

**Step 1: Create each item** using `todo/TODO_ITEM_TEMPLATE.md` as the format.
Content per item:
- `cli-full-pipeline-run` (priority medium, category feature): expose a
  full-pipeline CLI `run` wrapping `run_sample`, behind its own tests. New public
  semantics — deliberately excluded from the consolidation.
- `cli-find-reference` (priority medium, category feature): expose
  `find-reference` (steps 0/6 reference/shared calls) as an installed, tested
  subcommand.
- `cmd-run-drops-trim-thresholds` (priority high, category bug): `cmd_run` parses
  `--min-match/--min-trimmed/--mismatch` but calls
  `identify_TE_reads(reads, out, search_tool=aligner)`, silently dropping all
  three. Fixing changes trim output → validate separately. Relevant to the
  simulated-aligner benchmark.
- `run-aligner-bwa-unsupported` (priority high, category bug): `run --aligner bwa`
  is accepted by argparse (`choices=["minimap2","bwa"]`) but `identify_TE_reads`
  rejects non-minimap search tools; the advertised bwa path is non-functional.

**Step 2: Add four rows to `todo/TODO_REGISTRY.md`** (above the "Add new entries"
marker), each linking its file, dated 2026-07-10.

**Step 3: Commit**

```bash
git add todo/
git commit -m "docs(todo): file follow-ups (deferred CLI features + 2 cmd_run defects)"
```

---

## Task 8: Post-change validation gate

**Files:** none modified. Proves no behavioral regression.

**Step 1: Full test suite + lint**

Run: `pixi run test`
Expected: all PASS.
Run: `pre-commit run --all-files`
Expected: ruff/black/codespell/pydocstyle clean (fix any formatting the hooks flag, then re-commit).

**Step 2: Re-run the forced B_10 validation**

Run: `pixi run validate-rice --local --force B_10`
Expected: completes and regenerates `<report_dir>` + `<relocate3_outdir>/B_10/results`.

**Step 3: Compare sorted normalized records against the baseline**

Do NOT compare raw checksums. Compare sorted record content + metrics:
```bash
BASE=/scratch/nmath020/claude-4628/-bigdata-stajichlab-nmath020-github-github-tools-RelocaTE-RelocaTE3-jason-RelocaTE3/0b8a6ff7-eae4-4efb-8e91-8cade8bcacff/scratchpad/cli-consolidation-baseline
# Non-ref + characterized normalized records (sort to ignore ordering):
for f in $(cd "$BASE/report" && find . -name '*.tsv'); do
  echo "== $f =="
  diff <(sort "$BASE/report/$f") <(sort "<report_dir>/$f") && echo "  identical" || echo "  DIFF (investigate)"
done
# Summaries (recall / precision / Jaccard / counts):
diff "$BASE/report/nonref/summary.txt" "<report_dir>/nonref/summary.txt"
diff "$BASE/report/characterized/summary.txt" "<report_dir>/characterized/summary.txt"
# Genotype/status fields + expected output files:
diff <(sort "$BASE/B_10_results"/*.txt) <(sort "<relocate3_outdir>/B_10/results"/*.txt)
ls "$BASE/B_10_results" && ls "<relocate3_outdir>/B_10/results"
```
Expected: normalized records identical; recall/precision/Jaccard unchanged;
genotype/status fields unchanged; same set of expected output files. **Any diff
is a regression** — stop and fix, or explicitly justify in the PR.

**Step 4: Record the result** in the design doc's validation section (append a
short "Result: <date>, B_10 identical, N passed" note) and commit.

```bash
git add plans/2026-07-10-cli-consolidation-design.md
git commit -m "docs(plans): record B_10 validation parity for CLI consolidation"
```

---

## Done criteria

- One canonical CLI; `cli.py` is a shim; `RelocaTE3.cli.main` still reaches it.
- No installed `__main__.py` command behavior changed (only `run` help/docstring text).
- `pixi run test` green; both entry-point subprocess smoke tests present.
- Docs describe the installed staged interface; `run` explicitly documented as map+trim.
- Four follow-ups filed.
- B_10 normalized records + metrics identical before/after.
