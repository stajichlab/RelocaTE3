# Streaming SAM validation results

Recorded: 2026-09-28 10:13 PDT (America/Los_Angeles).
Project: RelocaTE3, branch `main`, following merged PR #49.
Status: validation passed; memory changes remain uncommitted. No release pin.

## Completed validation

SLURM job `29093144` ran on `r38` and finished successfully on
2026-09-26 at 01:50:40 PDT. This resolves the unconfirmed submission status
in `docs/2026-09-25-streaming-sam.md`.

The job ran the full test suite, paired isolated trimming comparisons, and
one complete BLAT/bwa-aln replay of riceTElib `cov30x_rep1`:

- 279 tests passed, with no failures, errors, or skips.
- Both isolated trims matched all eight reference artifacts byte-for-byte.
- Normalized calls and all four scored tables matched the previous validated
  run: `matches.tsv`, `correctness.tsv`, `precision.tsv`, and
  `false_positive_calls.tsv`.
- Calls remained 368: 302 true positives, 66 false positives, and 198 false
  negatives against 500 truth events. This patch preserves calling behavior;
  it does not improve accuracy.
- Input and executable hashes matched the previous replay. Current Python
  source matched the tested candidate snapshot.

## Resources

GNU time maximum resident set size, converted from KiB to GiB:

| Version | Peak RSS (GiB) | Adapter wall time |
| --- | ---: | ---: |
| Original implementation | 34.45 | 11:41:05 |
| Earlier memory improvements | 20.36 | 10:13:16 |
| Incremental SAM writer | 19.10 | 11:50:55 |

The incremental improvement is 6.2%; the combined reduction is 44.6%.
The latest replay reported 20,026,056 KiB maximum RSS and exit status zero.
It was 15.9% slower than the preceding replay. Runs used different nodes;
these measurements do not establish a speed improvement or isolate the cause
of the runtime difference. The unchanged trimming step still used about
19.03 GiB in the paired isolated profiles.

This is one representative 30x sample, not a rerun of the entire benchmark.
Do not lower global SLURM memory requests or declare release readiness from
this result alone.

## Evidence and reproducibility

All paths below are relative to the RelocaTE3 repository root:

- `logs/relocate3-memory-validation.29093144.log`
- `results/memory-profile/2026-09-25-streaming-sam-cov30x-rep1/`
- Previous baseline: `results/memory-profile/2026-09-23-memory-optimized-cov30x-rep1/`
- Original baseline: `results/memory-profile/2026-09-21-ricetelib-cov30x-rep1/`

The candidate directory contains the frozen source snapshot, input hashes,
tool metadata, `full-tests.xml`, `comparison.json`, `adapter.time-v.txt`,
trim comparisons, and sampled process memory. Read-only inspection commands:

```bash
tail -n 14 logs/relocate3-memory-validation.29093144.log
cat results/memory-profile/2026-09-25-streaming-sam-cov30x-rep1/comparison.json
cat results/memory-profile/2026-09-25-streaming-sam-cov30x-rep1/adapter.time-v.txt
```

No analysis was rerun and no previous output was removed during this review.
No validation failures were observed. Recommended next action: audit
RelocaTE3-specific false positives and missed detections relative to RelocaTE2,
using the existing full benchmark and read evidence before changing filters.
Keep the validated memory changes separate from any accuracy changes.
