# Memory optimization validation: passed

Recorded: 2026-09-25, 10:40 PDT (America/Los_Angeles).
Project: RelocaTE3, branch `main`, following merged PR #49.
Base revision: `942a5ea266bd000a47e27df7ef9e2536c877acbd`, with uncommitted
memory changes in `aligners.py` and `librelocate.py`.
Validation job: `29071085`, node `r33`, eight CPUs, 64 GB allocated.
Status: full suite, baseline/candidate trimming, and full single-sample replay
all passed. This review did not change calling code, submit jobs, or create tags.

## Outcome

The memory patch is supported by this validation: peak full-workflow memory
decreased 40.9%, with unchanged normalized insertion calls and truth scoring.
All 267 tests passed, with zero skips, failures, or errors (JUnit time 83.428 s).
The current working-tree Python source is identical to the validated candidate
source snapshot. Changes remain uncommitted; no stable release was published.

The job ran September 24, 21:57:29 through September 25, 08:51:36 PDT, including
tests, file hashing, two trim-only runs, and the full replay. The `.complete`
trim-validation sentinel and final `.profile_complete` sentinel both exist.

## Resource comparison

GNU time maximum RSS is converted from KiB to GiB using 1,048,576 KiB/GiB.
Full-workflow baseline is the completed September 22 profiling job `29005979`,
not the older full benchmark. Isolated trimming compares two fresh processes
in the same validation allocation, using the same saved TE BAMs and FASTQs.

| Measurement | Baseline | Updated | Change |
|---|---:|---:|---:|
| Full-workflow peak RSS | 34.45 GiB | 20.36 GiB | −40.9% |
| Isolated trimming peak RSS | 33.77 GiB | 19.03 GiB | −43.6% |
| Full-workflow wall time | 11:41:05 | 10:13:16 | −12.5% observed |
| Isolated trimming wall time | 18:04.90 | 18:56.02 | +4.7% observed |

Exact RSS values: full baseline 36,121,384 KiB; full updated 21,345,556 KiB;
trim baseline 35,408,588 KiB; trim updated 19,958,000 KiB. No swaps were reported.

The paired trim experiment directly supports the benefit of releasing each
mate's dictionary before parsing the next. Runtime is less conclusive: the
full runs used different nodes (`r26` versus `r33`) at different times, while
the sequential trim runs may differ in cache state. Do not claim a controlled
12.5% algorithmic speedup, or that every stage became faster.

## Biological/output equivalence

The isolated baseline and candidate trims each reproduce **all eight** original
trim artifacts byte-for-byte, as checked by SHA-256 in the compute job. This
includes both flanking FASTQs, both containing-read FASTQs, the selected-TE-hit
names, read-repeat table, and both TE-portion FASTAs. There are no missing,
additional, or differing files in the comparison.

The completed full replay contains the same 368 normalized calls:

- 302 true positives;
- 66 false positives;
- 198 false negatives against 500 truth events.

`normalized_calls`, `matches.tsv`, `precision.tsv`, `correctness.tsv`, and
`false_positive_calls.tsv` all pass the saved baseline comparison. Independent
row-multiset comparisons of those small tables during this review confirmed
the same result. This is not a claim that every intermediate BAM or raw report
is byte-identical.

Input/index SHA-256 manifests and the recorded executable hashes also agree
between baseline and candidate full profiles. The stale editable package
version string remains a separate release-packaging concern; the run logs
verify import from the explicitly frozen source snapshot.

This validates preservation of biology for the representative 30x riceTElib
sample, plus the full automated test suite. It does not constitute a rerun of
all 72 samples per caller, a new real-data validation, or an improvement in
false-positive/false-negative rates. Memory changes intentionally leave those
rates unchanged.

## Remaining memory peak and next candidate

The full candidate's 20.36-GiB peak is in the Python `relocaTE3 run` process,
at elapsed 17,186.362 seconds. External SAM-to-BAM conversion starts around
17,192.633 seconds. The peak therefore occurs during the first mate's
BLAT-result conversion, before trimming, rather than in a BLAT/BWA executable.
The earlier baseline conversion peaks were 20.68/20.80 GiB. The large overall
gain is consequently primarily associated with eliminating the later,
overlapping-mate trimming peak; sequence-buffer changes did not eliminate the
conversion peak.

Further source inspection identified a remaining whole-output collection:
`psl_to_sam` in `src/RelocaTE3/aligners.py` initializes `out = []`, appends every
converted SAM line, and returns the entire list. `_blat_side` then iterates
that list to write SAM. Thus all converted SAM strings coexist with the query
sequence dictionary. This is a concrete retention pattern consistent with
the remaining peak, not a measured per-allocation attribution.

Recommended next targeted change: stream converted SAM records to the writer.
Preserve the existing list-returning interface where compatibility requires
it, using an internal iterator for the production writing path. Test exact
record content/order, reverse strands, CIGAR/NM fields, alignment admission,
empty input, and incremental consumption. This also gives a bounded record
interface suitable for later Rust implementation without changing Nextflow's
file-level contracts.

Set expectations: isolated trimming still requires approximately 19 GiB, so
removing the SAM list alone cannot make the entire workflow fit RelocaTE2's
roughly 6.5-GiB footprint. Larger reductions will require compact or bounded
per-read processing, a separately scoped change. Do not lower production
memory requests solely from this one sample. Keep the false-positive/false-
negative audit separate from behavior-preserving memory changes.

## Evidence and lightweight review

Paths relative to the RelocaTE3 repository:

```text
logs/relocate3-memory-validation.29071085.log
results/memory-profile/2026-09-23-memory-optimized-cov30x-rep1/
  full-tests.xml
  adapter.time-v.txt
  comparison.json
  input_sha256.json
  tools.txt
  profile/summary.json
  profile/process_memory.tsv
  trim-validation/comparison.json
  trim-validation/baseline/time-v.txt
  trim-validation/candidate/time-v.txt
```

```bash
cat logs/relocate3-memory-validation.29071085.log
cat results/memory-profile/2026-09-23-memory-optimized-cov30x-rep1/comparison.json
cat results/memory-profile/2026-09-23-memory-optimized-cov30x-rep1/adapter.time-v.txt
```

Additional lightweight checks parsed the JUnit XML, compared scored-table row
multisets, compared saved input/executable hashes, compared current source with
its snapshot, and streamed the process-memory TSV to locate the peak. No
large BAM/FASTQ scans, alignments, or new tests were run during this review.

Next action: retain the validated memory patch and address the remaining
whole-SAM-list retention with a separately tested incremental writer.
