# CLI Dispatch Cleanup Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Modernize `src/RelocaTE3/cli.py` internals — extract `build_parser()`, dispatch by passing the `Namespace`, and drop the 8 `**kwargs` sinks — with zero command-behavior change.

**Architecture:** Three edits inside `cli.py`: (1) move `main()`'s inline subparser assembly into `build_parser()`; (2) change `parsed.func(**vars(parsed))` → `return parsed.func(parsed)`; (3) convert each `cmd_x(a, b, …, **kwargs)` to `cmd_x(args)` reading `args.foo`. The `cmd_run` threshold bug is preserved on purpose (documented) so B_10 stays byte-identical. Behavior-preserving refactor only; the `cmd_run` and `bwa` fixes are separate branches.

**Tech Stack:** Python 3.10+, argparse, pytest, pixi, ruff (F821 undefined-name is the safety net for missed `args.` prefixes).

**Design doc:** `plans/2026-07-11-cli-dispatch-cleanup-design.md` (approved).

**Reference skills:** @superpowers:executing-plans to drive; @superpowers:test-driven-development discipline (keep the 48-test suite green across each task).

---

## Task 0: Capture the B_10 baseline (before any edit)

**Files:** none modified.

**Step 1:** `pixi run test` → record (expect 48 passed).
**Step 2:** Resolve harness dirs (from repo root):
```bash
cd validation/real_rice && CONFIG="config.toml"
python3 -c "from _config import load_config; c=load_config('$CONFIG'); print(c['paths']['relocate3_outdir']); print(c['paths']['report_dir'])"; cd -
```
**Step 3:** `pixi run validate-rice --local --force B_10`
**Step 4:** Copy baseline OUT of the validation tree:
```bash
BASE=<scratchpad>/cli-dispatch-baseline; mkdir -p "$BASE"
cp -r validation/real_rice/report "$BASE/report"
cp -r validation/real_rice/results/B_10/results "$BASE/B_10_results"
```
(Use the session scratchpad dir for `<scratchpad>`.)

---

## Task 1: Extract `build_parser()` (behavior-identical; dispatch unchanged)

**Files:** Modify `src/RelocaTE3/cli.py`; Test `tests/main_test.py`.

**Step 1: Write the failing test** — add to `tests/main_test.py`:
```python
def test_build_parser_registers_subcommand():
    """build_parser() returns a parser whose subcommands dispatch to cmd_*."""
    from RelocaTE3.cli import build_parser, cmd_index_genome

    parsed = build_parser().parse_args(["index-genome", "-g", "ref.fa"])
    assert parsed.func is cmd_index_genome
```

**Step 2: Run — expect FAIL** (`ImportError: cannot import name 'build_parser'`):
`pixi run pytest tests/main_test.py::test_build_parser_registers_subcommand -v`

**Step 3: Implement** — in `src/RelocaTE3/cli.py`, add a `build_parser()` function (place it directly above `main()`) containing the parser + subparser assembly currently inside `main()` (the block from `parser = argparse.ArgumentParser(...)` through the last `_menu_characterize(...)` call, i.e. current lines ~629–710):
```python
def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level parser with all subcommands."""
    prog = __entry_points__.get(__name__, "relocaTE3")
    parser = argparse.ArgumentParser(
        prog=prog,
        formatter_class=CustomHelpFormatter,
        description=main.__doc__,
        epilog=f"Written by {__author__}",
    )
    parser.add_argument("-V", "--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.required = True
    # ... move the 8 _menu_*(subparsers.add_parser(...)) blocks here VERBATIM ...
    return parser
```
Then rewrite `main()` to call it, keeping the try/except and dispatch EXACTLY as they are for now:
```python
def main(args: list[str] | None = None) -> int:
    """Tool for identifying Transposable transposition from WGS data by comparison to a reference genome."""
    parser = build_parser()

    try:
        cli_args = args or sys.argv[1:]
        if not cli_args:
            parser.print_help(sys.stderr)
            raise SystemExit(0)
        parsed = parser.parse_args(cli_args)
        if getattr(parsed, "verbose", False):
            logger.setLevel("DEBUG")
            for handler in logger.handlers:
                handler.setLevel("DEBUG")
            logger.debug("Debug mode enabled.")
        parsed.func(**vars(parsed))          # UNCHANGED in this task
    except KeyboardInterrupt:
        logger.warning("Terminated by user.")
        return 1
    except SystemExit as err:
        if err.code != 0:
            logger.error(err)
            return 1
    except Exception as err:
        logger.error(err)
        return 1
    return 0
```
(`build_parser` references `main.__doc__` at call time, so `main` being defined below it is fine.)

**Step 4: Run — expect PASS**, and the whole suite green:
`pixi run pytest tests/main_test.py::test_build_parser_registers_subcommand -v` then `pixi run test` (48 + 1 = 49 passed).

**Step 5: Verify wiring unchanged:** `pixi run relocaTE3 --help` (8 subcommands, prog `relocaTE3`); `pixi run relocaTE3 run --help` (same flags).

**Step 6: Ruff** — `~/.local/bin/ruff check src/RelocaTE3/cli.py` → clean.

**Step 7: Commit**
```bash
git add src/RelocaTE3/cli.py tests/main_test.py
git commit -m "refactor(cli): extract build_parser() from main()"
```

---

## Task 2: Namespace dispatch + convert all 8 handlers (drop **kwargs)

**Files:** Modify `src/RelocaTE3/cli.py`, `tests/main_test.py`.

Dispatch and handler signatures MUST change together (dispatch passes the Namespace positionally, so handlers must accept one positional).

**Step 1: Change dispatch** in `main()`: `parsed.func(**vars(parsed))` → `return parsed.func(parsed)`. Keep the trailing `return 0` (it still handles the no-args/`--help` `SystemExit(0)` fall-through). Behavior identical because every handler returns 0.

**Step 2: Convert every `cmd_*` handler.** For each of the 8 handlers (`cmd_map`, `cmd_trim`, `cmd_run`, `cmd_characterize`, `cmd_annotate_ref`, `cmd_index_genome`, `cmd_align_genome`, `cmd_find_insertions`): change the signature to `def cmd_x(args) -> int:` and prefix every former-parameter reference in the body with `args.`. Delete the `**kwargs`. Bodies are otherwise unchanged.

Worked example — `cmd_index_genome` (smallest):
```python
# BEFORE
def cmd_index_genome(genome_fasta, force, verbose, **kwargs) -> int:
    """..."""
    from RelocaTE3.align import Aligner
    aln = Aligner()
    aln.verbose = bool(verbose)
    aln.index_genome(genome_fasta, force=force)
    logger.info("Indexed genome %s", genome_fasta)
    return 0
# AFTER
def cmd_index_genome(args) -> int:
    """..."""
    from RelocaTE3.align import Aligner
    aln = Aligner()
    aln.verbose = bool(args.verbose)
    aln.index_genome(args.genome_fasta, force=args.force)
    logger.info("Indexed genome %s", args.genome_fasta)
    return 0
```

Worked example — `cmd_run` (ADD the bug-preservation comment):
```python
def cmd_run(args) -> int:
    """Identify TE-containing reads and generate flanking reads (map + trim).

    Maps reads to the TE library then trims the TE-matching portion. This is
    steps 2-3 (TE-read identification and flank generation), not the complete
    insertion-calling pipeline.
    """
    from RelocaTE3.librelocate import RelocaTE
    from RelocaTE3.ReadLibrary import ReadLibrary

    fileset = [args.left] + ([args.right] if args.right else [])
    reads = ReadLibrary(fileset, args.name)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    relocate = RelocaTE(TElib=args.te_library, threads=args.threads, verbose=int(args.verbose))
    # BUG preserved intentionally (see todo/cmd-run-drops-trim-thresholds.md):
    # args.minimum_match_length / minimum_trimmed_length / mismatch_allowance are
    # parsed but NOT forwarded to identify_TE_reads. Fixing changes trim output;
    # do it in its own validated branch, not this refactor.
    n = relocate.identify_TE_reads(reads, out, search_tool=args.aligner)
    logger.info("%d read(s) written", n)
    return 0
```

Apply the same mechanical transform to the remaining 6. The `dest` names (= `args.` attribute names) come from each `_menu_*`; e.g. `cmd_characterize` reads `args.sites_file`, `args.bam`, `args.genome_fasta`, `args.outdir`, `args.excision`, `args.samtools`, `args.bcftools`, `args.verbose`.

**Step 3: Update `tests/main_test.py` mock** to receive a `Namespace`:
```python
def mockreturn(args):
    """Stand-in command handler that records its Namespace and succeeds."""
    mockreturn.args = args
    return 0
```
and update the assertions in the two tests that use it:
```python
# test_main_map
assert mockreturn.args.name == "HEG4"
assert mockreturn.args.left == "r1.fq"
assert mockreturn.args.te_library == "te.fa"
# test_main_characterize
assert mockreturn.args.sites_file == "sites.txt"
assert mockreturn.args.bam == ["a.bam", "b.bam"]
```
(The `monkeypatch.setattr(cli, "cmd_map"/"cmd_characterize", mockreturn)` targets are already correct.)

**Step 4: Run the suite** — `pixi run test` → all pass (49). The monkeypatched tests confirm `map`/`characterize` dispatch a Namespace; `test_run_cli_generates_flanking_reads` and the acceptance test exercise real handler bodies.

**Step 5: Ruff F821 gate (critical).** Several CLI handler bodies (`align-genome`, `find-insertions`, `annotate-ref`, `map`, `trim`) are NOT run by pytest — only by B_10. A missed `args.` prefix becomes an undefined name. Run:
`~/.local/bin/ruff check src/RelocaTE3/cli.py`
Expected: **All checks passed** — in particular ZERO `F821 undefined name`. If ruff reports F821 for a bare param name, you missed an `args.` prefix — fix it. Also `~/.local/bin/ruff format --check src/RelocaTE3/cli.py`.

**Step 6: Manual entry smoke** for a handler pytest doesn't cover, to catch a NameError ruff might miss:
`pixi run relocaTE3 index-genome --help` and `pixi run relocaTE3 align-genome --help` (parse OK), and confirm `pixi run relocaTE3 --version` works.

**Step 7: Commit**
```bash
git add src/RelocaTE3/cli.py tests/main_test.py
git commit -m "refactor(cli): dispatch Namespace to cmd_* handlers; drop **kwargs sinks"
```

---

## Task 3: B_10 byte-identical validation gate

**Files:** none modified.

**Step 1:** `pixi run test` green (49); `~/.local/bin/ruff check src/RelocaTE3/cli.py` clean.
**Step 2:** `pixi run validate-rice --local --force B_10`
**Step 3:** Compare sorted normalized records against the Task-0 baseline (NOT checksums):
```bash
BASE=<scratchpad>/cli-dispatch-baseline
for f in $(cd "$BASE/report" && find . -name '*.tsv'); do
  diff <(sort "$BASE/report/$f") <(sort "validation/real_rice/report/$f") >/dev/null && echo "OK $f" || echo "DIFF $f"
done
diff "$BASE/report/nonref/summary.txt" validation/real_rice/report/nonref/summary.txt
diff "$BASE/report/characterized/summary.txt" validation/real_rice/report/characterized/summary.txt
for f in "$BASE"/B_10_results/*; do b=$(basename "$f"); diff <(sort "$f") <(sort "validation/real_rice/results/B_10/results/$b") >/dev/null && echo "OK $b" || echo "DIFF $b"; done
```
Expected: every file OK / identical; recall/precision/Jaccard/TSD/status unchanged. **Any diff is a regression** — stop and fix.
**Step 4:** Append the result to the design doc and commit.

---

## Done criteria

- `cli.py` has `build_parser()`; `main()` dispatches `return parsed.func(parsed)`; all 8 handlers are `cmd_x(args)` with no `**kwargs`.
- `cmd_run` still drops the trim thresholds, with the guard comment.
- `pixi run test` green (49); ruff clean (no F821); B_10 byte-identical.
- `cmd_run` + `bwa` fixes remain unstarted (separate branches).
