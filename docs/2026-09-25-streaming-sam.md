# Incremental PSL-to-SAM writing

Recorded: 2026-09-25, 12:36 PDT (America/Los_Angeles).
Project: RelocaTE3, branch `main`, following merged PR #49.
Base Git revision: `942a5ea266bd000a47e27df7ef9e2536c877acbd` plus uncommitted
memory improvements. No commit, release pin, or tag created.
Status: implemented; focused tests passed; single-sample HPC validation prepared,
submission unconfirmed.

## Change and compatibility

The preceding validated patch reduced full-workflow peak RSS from 34.45 to
20.36 GiB without changing calls. The remaining peak occurs immediately before
SAM-to-BAM conversion, where `psl_to_sam` builds a complete list of SAM strings.

`src/RelocaTE3/aligners.py` now provides a private `_iter_psl_to_sam` generator.
It uses the same conversion/admission logic but yields each accepted alignment
instead of appending to a list. `_blat_side` writes directly from this iterator.
The public `psl_to_sam` function still returns a list by materializing that same
iterator, preserving existing callers' indexing, length, and eager-error behavior.
There are not two independent conversion implementations to keep synchronized.

No alignment thresholds, CIGAR construction, coordinate conventions, strand
handling, NM tags, read names, ordering, or duplicate/multimapping behavior were
changed. File and CLI contracts are unchanged. The previous sequence-dictionary
release before sorting remains in place; the exhausted generator does not keep
that dictionary alive. Conversion errors propagate before BAM sorting starts.

The production converter no longer retains a whole SAM-output list. It still
needs the selected-read sequence dictionary, and trimming still uses about
19 GiB in the most recent profile. Do not infer a particular whole-workflow
memory reduction until the new replay finishes. No Rust or Nextflow code was
added; the single-record conversion path is compatible with future native
acceleration behind unchanged file-level pipeline boundaries.

## Validation completed

All **28 focused tests passed**, with no skips:

```bash
PATH="$PWD/.pixi/envs/default/bin:$PATH" python -m pytest -q \
  tests/psl_streaming_test.py tests/memory_optimization_test.py \
  tests/aligners_test.py::TestPslToSam tests/aligners_test.py::TestBlatCommand
ruff check tests/psl_streaming_test.py tests/memory_optimization_test.py
git diff --check
bash -n scripts/validate_memory_changes.slurm
```

New tests cover exact plus/minus-strand SAM strings, reverse complementation,
soft clipping, two-block insertion/deletion CIGARs, NM tags, missing sequences,
all five alignment-admission limits, headers/empty input, duplicate preservation,
input order, lazy input consumption, late parse failures, and list compatibility.
The production writer test checks that the first yielded record is actually
written before the next record is requested, and that the sequence dictionary
is released before sorting. It would reject materializing the iterator back
into a list in the writing path.

A separate lightweight differential check imported the converter from the
previously validated source snapshot and used 2,000 seeded synthetic PSL lines
(seed `20260925`) plus two header/empty lines. The previous converter, public
wrapper, and new iterator produced identical 569-record outputs, both with and
without a query-sequence mapping. This check supplements the independent exact
SAM expectations in the tests; it is not a large-data performance measurement.

## Prepared replay and next action

New candidate:
`results/memory-profile/2026-09-25-streaming-sam-cov30x-rep1`.
Baseline for this incremental change:
`results/memory-profile/2026-09-23-memory-optimized-cov30x-rep1`.

All 63 checksummed candidate snapshot files were verified after preparation;
the saved aligner source matches the current working tree. Earlier benchmark
and profiling outputs remain untouched.

Reuse the established validation job: full test suite, baseline/candidate
trim-output checks, then a complete single-sample BLAT/bwa-aln replay with memory
sampling and normalized-call/truth-score comparisons. The unchanged trimming
checks are retained by this existing driver; no new benchmark framework was
introduced. Compare final resources with the latest 20.36-GiB profile, not just
the original 34.45-GiB run, to isolate this additional change.

Resources: eight CPUs, 64 GB, 24 hours on `epyc`. The readable local SLURM
configuration checked today gives `epyc MaxTime=30-00:00`; the request fits.
Submission with a 20-second timeout returned exit 124 without a numeric job ID.
No validation-start claim existed at inspection. The job is **not confirmed
queued or running**. Check the queue in a normal cluster terminal before retrying:

```bash
squeue -u "$USER" -o '%.18i %.40j %.10T %.10M %.30R'
```

If no matching job exists, submit from the RelocaTE3 repository root:

```bash
sbatch --parsable scripts/validate_memory_changes.slurm \
  results/memory-profile/2026-09-25-streaming-sam-cov30x-rep1 \
  results/memory-profile/2026-09-23-memory-optimized-cov30x-rep1
```

From `relocate-benchmark`, first `cd ../../RelocaTE3_jason/RelocaTE3`.
No cleanup or repeat preparation is needed. Output guards reject duplicate
or incomplete reruns instead of overwriting results. The next action is to
obtain a confirmed job ID and review memory and call-equality results when it
finishes. The new full suite and large-sample replay have not yet been verified.
