# User-Selectable Outputs Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a `relocaTE3 pipeline` command that runs the whole pipeline in one process and lets the user choose which outputs to retain via `--keep`, defaulting to final-only, while always writing a small run-summary.

**Architecture:** A new `outputs.py` module owns the output-category registry, `--keep` parsing, post-run cleanup, and the summary writer. `run_sample()` gains `te_aligner`/`genome_aligner` passthrough. A new `cmd_pipeline`/`_menu_pipeline` in `cli.py` wires `run_sample()` + summary + cleanup and registers the `pipeline` subcommand. Existing subcommands are untouched.

**Tech Stack:** Python 3.10+, argparse, pysam, pytest. Env: `pixi run --manifest-path pixi.toml test` (runs `pytest -ra -q`). Design: `plans/2026-07-22-selectable-outputs-design.md`.

**Conventions:** Tests live in `tests/<module>_test.py` (suffix `_test.py`, not `test_` prefix). Test data: `tests/data/{sim_reads,sim_genome,TE_lib}` (see `tests/pipeline_test.py`). Run a single test with `pixi run --manifest-path pixi.toml pytest tests/outputs_test.py::test_name -v`. Commit after each task. Use `-l` login shell for any SLURM/module commands (not needed here).

---

## Task 1: Category registry + `--keep` parsing

**Files:**
- Create: `src/RelocaTE3/outputs.py`
- Test: `tests/outputs_test.py`

**Step 1: Write the failing test**

```python
"""Tests for output-category selection, cleanup, and run summary."""
from __future__ import annotations

import pytest

from RelocaTE3.outputs import OUTPUT_CATEGORIES, KEEP_ALL, resolve_keep


def test_default_is_final_only():
    assert resolve_keep(None) == {"final"}
    assert resolve_keep("") == {"final"}


def test_comma_list_parsed():
    assert resolve_keep("final,te_containing") == {"final", "te_containing"}


def test_whitespace_tolerated():
    assert resolve_keep(" final , flanking ") == {"final", "flanking"}


def test_all_token_returns_keep_all_sentinel():
    assert resolve_keep("all") is KEEP_ALL


def test_unknown_category_errors():
    with pytest.raises(ValueError, match="unknown output category"):
        resolve_keep("final,bogus")


def test_registry_has_expected_categories():
    assert set(OUTPUT_CATEGORIES) == {
        "final", "flanking", "te_containing", "read_repeat",
        "te_portions", "alignments",
    }
```

**Step 2: Run test to verify it fails**

Run: `pixi run --manifest-path pixi.toml pytest tests/outputs_test.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'RelocaTE3.outputs'`.

**Step 3: Write minimal implementation**

```python
"""Output-category selection, post-run cleanup, and run summary.

The pipeline command writes every artifact, then retains only the categories the
user asked to keep (default: ``final``). Categories map to glob patterns relative
to the sample output directory. ``te_containing`` and ``read_repeat`` share the
``te_containing/`` directory, so cleanup deletes matched paths only (never rmdir a
shared parent).
"""
from __future__ import annotations

# Glob patterns are relative to the sample outdir.
OUTPUT_CATEGORIES: dict[str, tuple[str, ...]] = {
    "final": ("results",),  # directory: nonref/ref insert + characterized
    "flanking": ("flanking/*.flankingReads.fq",),
    "te_containing": ("te_containing/*.ContainingReads.fq",),
    "read_repeat": ("te_containing/*.read_repeat_name.txt",),
    "te_portions": ("te_portions/*.fa",),
    "alignments": ("genome_aln/*",),
}

# Sentinel: keep everything (no cleanup). Distinct from any category set.
KEEP_ALL = object()


def resolve_keep(spec: str | None) -> set[str] | object:
    """Parse a ``--keep`` spec into a set of category names.

    ``None``/empty -> ``{"final"}``. ``"all"`` -> :data:`KEEP_ALL`. Unknown names
    raise ``ValueError``.
    """
    if not spec or not spec.strip():
        return {"final"}
    names = {part.strip() for part in spec.split(",") if part.strip()}
    if "all" in names:
        return KEEP_ALL
    unknown = names - set(OUTPUT_CATEGORIES)
    if unknown:
        raise ValueError(
            f"unknown output category: {', '.join(sorted(unknown))}; "
            f"valid: {', '.join(sorted(OUTPUT_CATEGORIES))}, all"
        )
    return names
```

**Step 4: Run test to verify it passes**

Run: `pixi run --manifest-path pixi.toml pytest tests/outputs_test.py -v`
Expected: PASS (6 tests).

**Step 5: Commit**

```bash
git add src/RelocaTE3/outputs.py tests/outputs_test.py
git commit -m "feat(outputs): output-category registry + --keep parsing"
```

---

## Task 2: Post-run cleanup

**Files:**
- Modify: `src/RelocaTE3/outputs.py`
- Test: `tests/outputs_test.py`

**Step 1: Write the failing test** (append to `tests/outputs_test.py`)

```python
from pathlib import Path

from RelocaTE3.outputs import cleanup_outputs


def _populate(outdir: Path, name: str = "S") -> None:
    """Create one file per category (mirrors real layout)."""
    (outdir / "results").mkdir(parents=True)
    (outdir / "results" / f"{name}.all_nonref_insert.txt").write_text("x\n")
    (outdir / "flanking").mkdir()
    (outdir / "flanking" / f"{name}.left.flankingReads.fq").write_text("@r\nA\n+\nI\n")
    (outdir / "te_containing").mkdir()
    (outdir / "te_containing" / f"{name}.left.ContainingReads.fq").write_text("@r\nA\n+\nI\n")
    (outdir / "te_containing" / f"{name}.read_repeat_name.txt").write_text("r\tTE\n")
    (outdir / "te_portions").mkdir()
    (outdir / "te_portions" / f"{name}.five_prime.fa").write_text(">r\nAC\n")
    (outdir / "genome_aln").mkdir()
    (outdir / "genome_aln" / f"{name}.genome.bam").write_text("bam")
    (outdir / f"{name}.summary.tsv").write_text("metric\tvalue\n")


def test_cleanup_keeps_only_requested(tmp_path: Path):
    _populate(tmp_path)
    cleanup_outputs(tmp_path, {"final"}, summary_name="S.summary.tsv")
    assert (tmp_path / "results" / "S.all_nonref_insert.txt").exists()
    assert not (tmp_path / "flanking").exists() or not list((tmp_path / "flanking").glob("*.fq"))
    assert not list((tmp_path / "genome_aln").glob("*"))
    # summary always survives
    assert (tmp_path / "S.summary.tsv").exists()


def test_cleanup_shared_dir_partial(tmp_path: Path):
    """Keeping read_repeat but not te_containing removes only ContainingReads."""
    _populate(tmp_path)
    cleanup_outputs(tmp_path, {"read_repeat"}, summary_name="S.summary.tsv")
    assert (tmp_path / "te_containing" / "S.read_repeat_name.txt").exists()
    assert not (tmp_path / "te_containing" / "S.left.ContainingReads.fq").exists()


def test_cleanup_all_sentinel_noop(tmp_path: Path):
    from RelocaTE3.outputs import KEEP_ALL
    _populate(tmp_path)
    cleanup_outputs(tmp_path, KEEP_ALL, summary_name="S.summary.tsv")
    assert (tmp_path / "flanking" / "S.left.flankingReads.fq").exists()
    assert (tmp_path / "genome_aln" / "S.genome.bam").exists()
```

**Step 2: Run test to verify it fails**

Run: `pixi run --manifest-path pixi.toml pytest tests/outputs_test.py -k cleanup -v`
Expected: FAIL — `ImportError: cannot import name 'cleanup_outputs'`.

**Step 3: Write minimal implementation** (append to `src/RelocaTE3/outputs.py`)

```python
import logging
import shutil
from pathlib import Path

logger = logging.getLogger("RelocaTE3")


def cleanup_outputs(
    outdir: str | Path,
    keep: set[str] | object,
    summary_name: str | None = None,
) -> None:
    """Delete outputs for categories not in ``keep``.

    ``keep`` may be :data:`KEEP_ALL` (no-op). ``summary_name`` (a filename under
    ``outdir``) is never deleted. Matched files and directories are removed;
    shared parent directories are left in place.
    """
    if keep is KEEP_ALL:
        return
    outdir = Path(outdir)
    drop = set(OUTPUT_CATEGORIES) - set(keep)  # type: ignore[arg-type]
    for cat in sorted(drop):
        for pattern in OUTPUT_CATEGORIES[cat]:
            for path in outdir.glob(pattern):
                if summary_name and path.name == summary_name:
                    continue
                if path.is_dir():
                    shutil.rmtree(path, ignore_errors=True)
                else:
                    path.unlink(missing_ok=True)
        logger.info("Cleanup: dropped output category %s", cat)
```

**Step 4: Run test to verify it passes**

Run: `pixi run --manifest-path pixi.toml pytest tests/outputs_test.py -v`
Expected: PASS (all tests).

**Step 5: Commit**

```bash
git add src/RelocaTE3/outputs.py tests/outputs_test.py
git commit -m "feat(outputs): post-run cleanup of non-kept categories"
```

---

## Task 3: Run summary

**Files:**
- Modify: `src/RelocaTE3/outputs.py`
- Test: `tests/outputs_test.py`

**Step 1: Write the failing test** (append)

```python
from RelocaTE3.outputs import gather_counts, write_summary


def test_gather_counts_from_files(tmp_path: Path):
    _populate(tmp_path)
    counts = gather_counts(tmp_path, "S")
    assert counts["reads_te_containing"] == 1
    assert counts["flanking_reads_written"] == 1
    assert counts["nonref_insertions"] == 1  # one data line


def test_write_summary_creates_tsv(tmp_path: Path):
    path = write_summary(tmp_path, "S", {"reads_te_containing": 42})
    assert path == tmp_path / "S.summary.tsv"
    text = path.read_text()
    assert text.splitlines()[0] == "metric\tvalue"
    assert "reads_te_containing\t42" in text
```

**Step 2: Run test to verify it fails**

Run: `pixi run --manifest-path pixi.toml pytest tests/outputs_test.py -k summary -v`
Expected: FAIL — `ImportError: cannot import name 'gather_counts'`.

**Step 3: Write minimal implementation** (append to `outputs.py`)

```python
def _count_fastq_records(paths) -> int:
    total = 0
    for p in paths:
        with open(p) as fh:
            total += sum(1 for _ in fh)
    return total // 4


def _count_data_lines(path: Path) -> int:
    if not path.exists():
        return 0
    with open(path) as fh:
        lines = [ln for ln in fh if ln.strip()]
    return max(0, len(lines) - 1)  # minus header


def gather_counts(outdir: str | Path, name: str) -> dict[str, int]:
    """Collect summary counts from on-disk outputs (call BEFORE cleanup)."""
    outdir = Path(outdir)
    counts: dict[str, int] = {}
    counts["reads_te_containing"] = _count_fastq_records(
        (outdir / "te_containing").glob(f"{name}.*.ContainingReads.fq")
    )
    counts["flanking_reads_written"] = _count_fastq_records(
        (outdir / "flanking").glob(f"{name}.*.flankingReads.fq")
    )
    counts["nonref_insertions"] = _count_data_lines(
        outdir / "results" / f"{name}.all_nonref_insert.txt"
    )
    return counts


def write_summary(outdir: str | Path, name: str, counts: dict[str, int]) -> Path:
    """Write ``<outdir>/<name>.summary.tsv`` (always kept)."""
    outdir = Path(outdir)
    path = outdir / f"{name}.summary.tsv"
    with open(path, "w") as fh:
        fh.write("metric\tvalue\n")
        for metric, value in counts.items():
            fh.write(f"{metric}\t{value}\n")
    return path
```

**Step 4: Run test to verify it passes**

Run: `pixi run --manifest-path pixi.toml pytest tests/outputs_test.py -v`
Expected: PASS.

**Step 5: Commit**

```bash
git add src/RelocaTE3/outputs.py tests/outputs_test.py
git commit -m "feat(outputs): run-summary counts + writer"
```

---

## Task 4: Thread aligners through `run_sample()`

**Files:**
- Modify: `src/RelocaTE3/pipeline.py` (`run_sample` signature + the two call sites)
- Test: `tests/pipeline_test.py`

**Context:** `run_sample` currently calls `identify_TE_reads(...)` and
`align_to_genome(...)` with no aligner arg (both default to minimap2).
`identify_TE_reads` already accepts `te_aligner=`; `align_to_genome` already
accepts `genome_aligner=`. Add two params to `run_sample` and forward them.

**Step 1: Write the failing test** (append to `tests/pipeline_test.py`)

```python
import inspect
from RelocaTE3.pipeline import run_sample as _run_sample_fn


def test_run_sample_accepts_aligner_params():
    params = inspect.signature(_run_sample_fn).parameters
    assert "te_aligner" in params and params["te_aligner"].default == "minimap2"
    assert "genome_aligner" in params and params["genome_aligner"].default == "minimap2"


def test_run_sample_bwa_te_aligner_runs(tmp_path):
    """run_sample honors an explicit te_aligner end-to-end on the test data."""
    reads = ReadLibrary([str(R1), str(R2)], "HEG4")
    gff = run_sample(
        reads, str(TELIB), str(GENOME), tmp_path, threads=2, te_aligner="bwa"
    )
    assert gff.exists()
```

**Step 2: Run test to verify it fails**

Run: `pixi run --manifest-path pixi.toml pytest tests/pipeline_test.py -k aligner -v`
Expected: FAIL — `assert 'te_aligner' in params`.

**Step 3: Write minimal implementation**

In `src/RelocaTE3/pipeline.py`, add to the `run_sample` signature (after
`mismatch_allowance` / near `required_junction_reads`):

```python
    te_aligner: str = "minimap2",
    genome_aligner: str = "minimap2",
```

Update the `identify_TE_reads` call to pass `te_aligner=te_aligner`:

```python
    relocate.identify_TE_reads(
        reads,
        outdir,
        TE_library=te_library,
        te_aligner=te_aligner,
        len_cut_match=len_cut_match,
        len_cut_trim=len_cut_trim,
        mismatch_allowance=mismatch_allowance,
    )
```

Update the `align_to_genome` call to pass `genome_aligner=genome_aligner`:

```python
    genome_bam, fullreads_bam = align_to_genome(
        reads, genome, outdir, threads=threads, genome_aligner=genome_aligner
    )
```

**Step 4: Run test to verify it passes**

Run: `pixi run --manifest-path pixi.toml pytest tests/pipeline_test.py -v`
Expected: PASS (existing + new tests).

**Step 5: Commit**

```bash
git add src/RelocaTE3/pipeline.py tests/pipeline_test.py
git commit -m "feat(pipeline): thread te_aligner/genome_aligner through run_sample"
```

---

## Task 5: `pipeline` CLI command (wiring)

**Files:**
- Modify: `src/RelocaTE3/cli.py` (add `_menu_pipeline`, `cmd_pipeline`, register subparser)
- Test: `tests/main_test.py`

**Step 1: Write the failing test** (append to `tests/main_test.py`)

```python
def test_pipeline_help_registered():
    from RelocaTE3.cli import main
    assert main(["pipeline", "--help"]) == 0
```

(If `tests/main_test.py` lacks imports, add `from RelocaTE3.cli import main` at top.)

**Step 2: Run test to verify it fails**

Run: `pixi run --manifest-path pixi.toml pytest tests/main_test.py -k pipeline -v`
Expected: FAIL — argparse errors "invalid choice: 'pipeline'" (SystemExit, non-zero).

**Step 3: Write minimal implementation** in `src/RelocaTE3/cli.py`

Add `cmd_pipeline` (near the other `cmd_*` handlers). It imports the outputs
helpers, resolves `--keep` first (so a bad value fails fast before compute), runs
the pipeline, writes the summary, then cleans up:

```python
def cmd_pipeline(args: argparse.Namespace) -> int:
    """Run the full pipeline in one process and retain selected outputs."""
    from pathlib import Path

    from RelocaTE3.outputs import (
        cleanup_outputs,
        gather_counts,
        resolve_keep,
        write_summary,
    )
    from RelocaTE3.pipeline import run_sample
    from RelocaTE3.ReadLibrary import ReadLibrary

    keep = resolve_keep(args.keep)  # raises ValueError on unknown category

    fileset = [args.left] + ([args.right] if args.right else [])
    reads = ReadLibrary(fileset, args.name)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    run_sample(
        reads,
        args.te_library,
        args.genome,
        outdir,
        repeatmasker=args.repeatmasker,
        genotype=args.genotype,
        threads=args.threads,
        len_cut_match=args.minimum_match_length,
        len_cut_trim=args.minimum_trimmed_length,
        mismatch_allowance=args.mismatch_allowance,
        te_aligner=args.te_aligner,
        genome_aligner=args.genome_aligner,
        verbose=int(args.verbose),
    )

    counts = gather_counts(outdir, args.name)
    summary = write_summary(outdir, args.name, counts)
    cleanup_outputs(outdir, keep, summary_name=summary.name)
    logger.info("Pipeline outputs retained: %s; summary at %s", args.keep or "final", summary)
    return 0
```

Add `_menu_pipeline` (mirror `_menu_run`, add `--genome`, `--repeatmasker`,
`--genotype`, `--genome-aligner`, `--keep`):

```python
def _menu_pipeline(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Arguments for the end-to-end 'pipeline' command."""
    parser.add_argument("-l", "--left", "--r1", required=True, metavar="R1",
                        help="Left/R1 read file (FASTQ)")
    parser.add_argument("-r", "--right", "--r2", metavar="R2",
                        help="Right/R2 read file for paired-end (FASTQ)")
    parser.add_argument("-T", "--te-library", required=True, dest="te_library",
                        metavar="FASTA", help="Transposon library FASTA")
    parser.add_argument("-g", "--genome", required=True, metavar="FASTA",
                        help="Reference genome FASTA")
    parser.add_argument("-n", "--name", required=True,
                        help="Sample name (output prefix)")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory")
    parser.add_argument("--threads", type=int, default=1, help="CPU threads")
    parser.add_argument("--repeatmasker", default=None,
                        help="RepeatMasker .out for reference-TE filtering (optional)")
    parser.add_argument("--genotype", action="store_true",
                        help="Characterize insertions (zygosity/excision)")
    parser.add_argument("--te-aligner", "--aligner", dest="te_aligner",
                        default="minimap2", choices=list(TE_ALIGNERS),
                        help="Aligner for TE-library search")
    parser.add_argument("--genome-aligner", dest="genome_aligner",
                        default="minimap2", choices=list(GENOME_ALIGNERS),
                        help="Aligner for genome re-alignment")
    parser.add_argument("--min-match", type=int, default=10,
                        dest="minimum_match_length", help="Min TE match length")
    parser.add_argument("--min-trimmed", type=int, default=10,
                        dest="minimum_trimmed_length", help="Min trimmed flank length")
    parser.add_argument("--mismatch", type=int, default=0,
                        dest="mismatch_allowance", help="Allowed TE mismatches")
    parser.add_argument("--keep", default=None,
                        help="Comma-separated output categories to retain "
                             "(final,flanking,te_containing,read_repeat,"
                             "te_portions,alignments,all). Default: final")
    _add_common_args(parser)
    parser.set_defaults(func=cmd_pipeline)
    return parser
```

Register the subparser in `build_parser()` (alongside the others, after `_menu_run(...)`):

```python
    _menu_pipeline(
        subparsers.add_parser(
            "pipeline",
            formatter_class=CustomHelpFormatter,
            help="Run the full pipeline end-to-end and retain selected outputs",
            description="Run steps 3-5 (plus optional reference annotation and "
            "genotyping) in one process, write a run summary, and retain only the "
            "output categories requested via --keep (default: final).",
        )
    )
```

**Step 4: Run test to verify it passes**

Run: `pixi run --manifest-path pixi.toml pytest tests/main_test.py -k pipeline -v`
Expected: PASS.

**Step 5: Commit**

```bash
git add src/RelocaTE3/cli.py tests/main_test.py
git commit -m "feat(cli): end-to-end 'pipeline' command with --keep"
```

---

## Task 6: End-to-end integration test

**Files:**
- Test: `tests/pipeline_test.py`

**Step 1: Write the failing test** (append)

```python
def test_pipeline_command_keep_final_only(tmp_path):
    """`pipeline --keep final` keeps results + summary, drops intermediates."""
    rc = main([
        "pipeline", "-l", str(R1), "-r", str(R2), "-T", str(TELIB),
        "-g", str(GENOME), "-n", "HEG4", "-o", str(tmp_path), "--threads", "2",
    ])
    assert rc == 0
    assert (tmp_path / "results" / "HEG4.all_nonref_insert.txt").exists()
    assert (tmp_path / "HEG4.summary.tsv").exists()
    # intermediates dropped
    assert not list((tmp_path / "flanking").glob("*.fq")) if (tmp_path / "flanking").exists() else True
    assert not list((tmp_path / "genome_aln").glob("*")) if (tmp_path / "genome_aln").exists() else True


def test_pipeline_command_keep_te_containing(tmp_path):
    """`--keep te_containing` retains TE-containing reads and drops final."""
    rc = main([
        "pipeline", "-l", str(R1), "-r", str(R2), "-T", str(TELIB),
        "-g", str(GENOME), "-n", "HEG4", "-o", str(tmp_path), "--threads", "2",
        "--keep", "te_containing",
    ])
    assert rc == 0
    assert list((tmp_path / "te_containing").glob("HEG4.*.ContainingReads.fq"))
    assert (tmp_path / "HEG4.summary.tsv").exists()  # summary survives
    assert not (tmp_path / "results").exists() or not list((tmp_path / "results").glob("*"))
```

**Step 2: Run test to verify it fails then passes**

Run: `pixi run --manifest-path pixi.toml pytest tests/pipeline_test.py -k pipeline_command -v`
Expected: PASS (the feature is implemented; this test guards it). If it FAILS on
the `alignments`/`te_containing` glob, verify the real on-disk filenames produced
by `run_sample` on the test data and adjust `OUTPUT_CATEGORIES` in `outputs.py`
to match, then re-run. **Verification command:**
`ls -R <a tmp outdir from test_run_sample_produces_outputs>` — or add a throwaway
`print` — to confirm `genome_aln/*.bam`, `te_containing/*.ContainingReads.fq`,
and any TE-library map BAM location. If map-stage BAMs land at the outdir root,
add their glob to the `alignments` category.

**Step 3: Commit**

```bash
git add tests/pipeline_test.py
git commit -m "test(pipeline): end-to-end --keep integration coverage"
```

---

## Task 7: Docs

**Files:**
- Modify: `README.md` (add a `pipeline` usage section)
- Modify: `plans/FEATURES.md` (note the new command, if that file tracks features)

**Step 1:** Add a short README section showing:

```
relocaTE3 pipeline -l R1.fq -r R2.fq -T te.fa -g ref.fa -n SAMPLE -o OUT \
    [--genotype] [--te-aligner bwa] [--genome-aligner bwa] [--keep final]
```

with the category list and the default (`final`), and a one-line note that the
`<name>.summary.tsv` is always written.

**Step 2: Commit**

```bash
git add README.md plans/FEATURES.md
git commit -m "docs: document the pipeline command and --keep"
```

---

## Task 8: Full test suite + wrap-up

**Step 1:** Run the whole suite:

Run: `pixi run --manifest-path pixi.toml test`
Expected: all pass (no regressions in existing subcommand/pipeline tests).

**Step 2:** If green, the branch `feat-selectable-outputs` is ready. Use
superpowers:finishing-a-development-branch to open the PR (user pushes/merges).

**Step 3 (follow-up, separate):** After merge, re-pin the benchmark
(`callers/relocate3/pixi.toml`) to the new RelocaTE3 rev if/when the benchmark
should adopt the `pipeline` command. Not required for this feature.

---

## Notes / risks

- **`OUTPUT_CATEGORIES` globs must match real filenames.** Task 6's verification
  step is where this is confirmed on the test data. The genome BAM from
  `run_sample` is `genome_aln/<name>.genome.bam` (see existing
  `test_run_sample_produces_outputs`); adjust if map-stage TE BAMs also need
  sweeping under `alignments`.
- **Cleanup only on success:** `cmd_pipeline` runs cleanup after `run_sample`
  returns; if `run_sample` raises, the exception propagates and cleanup is
  skipped (outputs kept for debugging) — no try/except needed.
- **Backward compatibility:** existing subcommands and `run_sample`'s existing
  callers (`run_samples`) are unchanged (new params default to minimap2).
```
