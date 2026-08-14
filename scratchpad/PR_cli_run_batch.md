# feat(cli): `run-batch` — run a cohort from a sample sheet or a FASTQ directory

Closes the last RelocaTE2 usability gap: RelocaTE2 took `--fq_dir` and processed
every sample in it, while RelocaTE3 required one invocation per sample with
explicit `--left`/`--right`. Running a cohort meant hand-writing a loop.

Stacks on `feat/release-readiness` (branch base: `b033d57`).

```bash
# a directory of paired FASTQs (the --fq_dir equivalent)
relocaTE3 run-batch --fq-dir reads/ -T RiceTE.fa -g reference.fa \
  -o results --threads 8 --jobs 4 --repeatmasker reference.fa.RepeatMasker.out

# or an explicit sample sheet
relocaTE3 run-batch --samples samples.csv -T RiceTE.fa -g reference.fa -o results
```

Each sample runs the ordinary single-sample pipeline into `<outdir>/<sample>/`.

## What's here

**`src/RelocaTE3/samples.py` (new)** — resolves either input to the same explicit
`Sample` list.

- `discover_fastq_dir()` pairs `<sample>_R1`/`_R2` and `<sample>_1`/`_2` (also with
  `.` or `-` before the marker).
- `read_sample_sheet()` reads CSV or TSV, skips `#` comments and blank lines, and
  accepts optional per-row `te_library` / `reference_genome` / `repeatmasker`
  columns that override the global flags.

**`cmd_run_batch` in `cli.py`** — `--jobs` for sample-level concurrency (each
sample still uses `--threads`), resume, fail-fast/`--keep-going`, and a non-zero
exit with a per-sample failure summary.

**Docs** — README "Many samples" section, `run-batch` row in the subcommand
table, and the RelocaTE2 migration table now points `--fq_dir` at it.

## Design decisions worth reviewing

**Not built on `pipeline.run_samples()`.** For the same reason `run-all` isn't
built on `run_sample`: those route through the module-level
`insertions.find_insertions`, not the `InsertionFinder` path the rice benchmark
validated. Building on them would have shipped a second, unvalidated caller under
a friendlier name. A batch is exactly N `cmd_run_all` invocations, so cohort
results are identical to what users get today one sample at a time.

**Sheet schema matches `plans/FEATURES.md:81`** — `sample_id,r1_fq,r2_fq`, with
nf-core spellings (`sample`, `fastq_1`, `fastq_2`) accepted as aliases. A sheet
written now stays valid when the planned Nextflow entry point lands.

**Fails before starting, not twenty samples in.** Every referenced read file is
checked up front; a typo in row 20 is reported before any compute is spent.

**An orphaned `_R1` raises instead of silently becoming single-end.** A
half-present pair is far more often a staging mistake than a deliberate SE run,
and the silent version costs a full cohort re-run to notice. The error names the
sample and points at `run-all --left` for a genuine SE case.

**Resumable, per the repo's idempotence convention.** A sample whose results
table already exists is skipped unless `--force`.

**Per-stage tuning flags are defined once** (`_add_pipeline_tuning_args`) and
shared by `run-all` and `run-batch`, so the two entry points cannot drift apart.

## Verification

- Full suite: **126 passed, 1 skipped** (`blat` not on PATH).
- End-to-end on real data via `--fq-dir --jobs 2`: two samples, **158 calls
  each** — identical to the single-sample baseline — both with `all_ref_insert`
  output.
- Re-running the same batch skipped both samples in **0.26 s**.
- Sample-sheet path with a third sample: 158 calls.
- `relocaTE3 --help` lists `run-all`, `run-batch`, `find-reference`.

New tests: `tests/samples_test.py` (13, sheet parsing + `--fq-dir` discovery) and
8 `run-batch` cases in `tests/main_test.py` covering resume, `--force`,
fail-fast, `--keep-going`, and per-row overrides beating global flags.
