# Memory-lifetime and BLAT buffering changes

Recorded: 2026-09-23, 22:04 PDT (America/Los_Angeles).
Project: RelocaTE3, branch `main`, following merged PR #49.
Base commit: `942a5ea266bd000a47e27df7ef9e2536c877acbd`.
Status: implementation and focused tests complete; HPC validation prepared,
but submission unconfirmed. Changes are uncommitted. No release tag or pin.

## Why these changes

Completed baseline profile `29005979` reproduced all benchmark calls and
localized the largest peak (34.45 GiB) to Python trimming, with an earlier
20.80-GiB peak in Python BLAT-result sequence recovery/conversion. See
`2026-09-23-memory-profile-results.md` for measured baseline evidence.

This patch changes memory ownership and buffering, not alignment settings,
read admission, family assignment, junction evidence, TSD inference, or
insertion filtering. New large-sample memory values have not yet been measured.

## Implementation

1. `src/RelocaTE3/librelocate.py`, `write_trimmed_reads`: explicitly release
   each completed mate's `coord` dictionary after writing its outputs. Python
   otherwise builds the next dictionary before replacing the previous one,
   keeping both live during parsing. Function signatures and ordering are
   unchanged. No forced garbage collection or allocator-specific hooks.

2. `src/RelocaTE3/aligners.py`, `_query_sequences`: write seqtk output to an
   automatically cleaned temporary file and iterate its FASTA lines. This
   removes the extra full stdout string and `splitlines()` list. Write the
   selected names incrementally too. The no-seqtk fallback also iterates its
   input instead of using `read().splitlines()`. Subprocess failures still
   propagate; the temporary file is closed on success and failure.

3. `BlatBackend._blat_side`: release the sequence dictionary after writing
   SAM and before BAM conversion/sorting. That dictionary is no longer needed
   once SAM contains the sequences.

Tradeoff: the seqtk path now uses temporary disk space and additional I/O
instead of redundant Python buffers. The required selected-sequence dictionary
and one mate's selected-record dictionary still scale with matched reads.
This is not yet a fully bounded-memory or chunk-streaming implementation.

## What the RelocaTE2 source showed

Inspected the user-specified checkout at `../references/RelocaTE2`, revision
`b7746986af5b662f8537432ed1b1ed819a046406`:

- `scripts/relocaTE2.py:91` splits input into 200,000-read chunks when
  `--split` is enabled. Some nearby legacy comments say one million; the
  executable command specifies 200,000.
- Its step-3 loop, around line 533, writes a separate BLAT-and-trim command
  sequence per input chunk. The benchmark RelocaTE2 adapter passes `--split`.
- `scripts/relocaTE_trim.py:40` parses PSL incrementally and retains selected
  alignment metadata by read. Around line 340, it streams the original FASTQ
  record-by-record, using current sequence/quality to write trimmed outputs.
- The BLAT path trims directly from PSL plus original FASTQ, avoiding
  RelocaTE3's PSL-to-SAM sequence reconstruction. Separate chunk processes
  also provide clear memory-lifetime boundaries.

RelocaTE2 is not constant-memory in every stage: its selected-hit dictionaries
still grow with each chunk, and several processes may run concurrently.
The useful lessons are bounded work units, streamed sequences, and explicit
ownership, not copying its shell orchestration or changing read filters.

## Rust and Nextflow considerations

No Rust or Nextflow implementation was added. Existing CLI subcommands,
function signatures, and FASTQ/BAM/read-table contracts remain stable. This
keeps the validated Python reference usable while a future native
implementation replaces internal parsing/scoring at those boundaries.

For future work, prefer file/record-batch interfaces with bounded ownership
over transferring multi-million-entry Python dictionaries across a Rust
boundary. Rust alone does not make unbounded retention bounded. A later
chunk-based map/trim implementation can build on RelocaTE2's design, but must
preserve deterministic selection, read names, original qualities, ordering,
and the all-TE-hit artifact required for downstream mate-state decisions.
That is a larger, separately benchmarked change, not part of this patch.

Nextflow should orchestrate the existing explicit file-producing steps;
biological decisions remain in the caller, not the workflow. Keeping a
working single-node CLI also avoids coupling correctness to a scheduler.

## Tests

Added six focused cases in `tests/memory_optimization_test.py`:

- completed mate container is not live when the next mate is parsed;
- fallback FASTA parsing is iterative and preserves multiline selected reads;
- seqtk output goes to a real file instead of captured Python stdout;
- seqtk failure propagates and its temporary handle closes;
- empty PSL needs neither query reads nor a subprocess;
- query-sequence dictionary is released before BAM sorting.

Before the code patch, four of the first five cases failed on the old behavior.
Afterward, 31 selected tests pass, including real seqtk conversion/loading,
reverse-strand trimming, original-quality restoration, deterministic family
selection, and profiler/validation helpers. No selected tests were skipped.

```bash
PATH="$PWD/.pixi/envs/default/bin:$PATH" python -m pytest -q \
  tests/memory_optimization_test.py tests/memory_profile_test.py \
  tests/trim_reverse_strand_test.py tests/te_family_determinism_test.py \
  tests/aligners_test.py::TestBlatCommand
ruff check scripts/profile_benchmark_memory.py scripts/validate_trim_memory.py \
  tests/memory_profile_test.py tests/memory_optimization_test.py
bash -n scripts/validate_memory_changes.slurm
git diff --check
```

The selected tests use tiny fixtures; substantial testing/benchmark work stays
in SLURM. The complete suite has not yet been run against this patch.

## Prepared validation

Existing completed baseline (unchanged):
`results/memory-profile/2026-09-21-ricetelib-cov30x-rep1`.

New candidate (63 checksummed source/config/adapter/helper/baseline files):
`results/memory-profile/2026-09-23-memory-optimized-cov30x-rep1`.
The manifest records the base Git revision plus the uncommitted source diff;
the snapshot contains the actual updated source, not just the base revision.

One job performs three gates in order:

1. Full test suite against the candidate source snapshot (test files from
   the current repository), with JUnit output at `full-tests.xml`.
2. Independent baseline and candidate trim-only processes using the same saved
   TE BAMs and original FASTQs. Measure each separately and require byte-for-byte
   SHA-256 equality for all eight trim artifacts against the completed replay,
   including empty files, original qualities, and unclassified TE-hit names.
   Reports go under `trim-validation/{baseline,candidate}/`.
3. Complete one-sample candidate profile through BLAT, trimming, genome
   alignment, insertion calling, and characterization. This measures the BLAT
   buffering change as well as trimming and requires normalized calls and all
   scored tables to match the saved benchmark. Standard profile reports and
   `.profile_complete` are written only on success.

The trim-only comparison avoids alignment cost when measuring the dictionary
lifetime change. The subsequent full replay is still necessary to evaluate
the separate BLAT sequence-buffering change in context. This is one sample,
not the entire six-caller benchmark panel. No old outputs are removed.

Resource request: `epyc`, eight CPUs, 64 GB, 24 hours. The baseline full replay
took 11h41m, plus approximately 23 minutes in its trimming interval; the request
leaves headroom for two isolated trims, tests, and file hashing. Actual time
is not guaranteed. Cached local SLURM configuration confirms `epyc`
`MaxTime=30-00:00`, so the requested wall time is compatible.

## Submission and next action

The controller query timed out. Submission using `timeout 20s sbatch --parsable`
also timed out (exit 124) without a numeric job ID; no `validation-started`
claim existed at inspection. The job is not confirmed queued/running.
Check the queue in an ordinary cluster terminal before retrying:

```bash
squeue -u "$USER" -o '%.18i %.40j %.10T %.10M %.30R'
```

If no matching validation job exists, submit from the RelocaTE3 repository root:

```bash
sbatch --parsable scripts/validate_memory_changes.slurm \
  results/memory-profile/2026-09-23-memory-optimized-cov30x-rep1 \
  results/memory-profile/2026-09-21-ricetelib-cov30x-rep1
```

From `relocate-benchmark`, first `cd ../../RelocaTE3_jason/RelocaTE3`.
Do not rerun preparation or delete the baseline. An exclusive validation-start
directory and per-stage output guards reject duplicate/incomplete reruns.
Failure logs remain available; use a new candidate directory for an intentional
retry after investigating failures. Do not update the shared tool environment
or repository tests while validation runs.

Next action: obtain a confirmed job ID and review both memory profiles and
the final call-equality checks before accepting the memory changes for release.
