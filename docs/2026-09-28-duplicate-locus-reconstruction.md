# Duplicate-locus reconstruction and regression fixtures

Follow-up: [the release deep dive](2026-09-28-release-deep-dive.md) documents
the implemented pairing candidate, passing regressions, and a newly verified
benchmark coordinate mismatch. The expected failures described below record
the earlier investigation; they have now been replaced by passing tests.
To re-extract these historical fixtures after the production change, use the
frozen pre-patch source via
`PYTHONPATH=results/pairing-replay/2026-09-28/baseline/src` and a new output
directory. The extractor intentionally requires reproducing historical calls.

Recorded: 2026-09-28 14:32 PDT (America/Los_Angeles).
Project: RelocaTE3, branch `main`, base `942a5ea`, following merged PR #49.
Status: two complete local clusters reconstructed; portable fixtures and
focused tests added. No production calling code changed, no benchmark rerun,
no commit or release pin. Existing memory changes were preserved.

## Finding: a coordinate convention changes the paired reads

The original report correctly identified two same-family neighboring calls,
but the mechanism is not duplicate use of junction read names. Their junction
read sets are disjoint. Both clusters have this geometry, with `p` the truth
coordinate and coordinates expressed as 1-based inclusive:

| Fixture | Left reads ending at p | Right reads starting at p | Right reads starting at p+1 |
|---|---:|---:|---:|
| cov30x_rep1, Os2571, p=24775704 | 21 | 1 | 26 |
| cov30x_rep2, Os0980, p=9006303 | 32 | 1 | 32 |

R3 compares left endpoint `p` with right start coordinates. Its nearest match
is the singleton right read at `p`, producing a one-base `A` TSD. The larger
right group at `p+1` is unclaimed and becomes a separate one-sided candidate.
Different starts prevent same-start consolidation. Both candidates pass the
recorded full-read and quality filters and cluster arbitration. The one-sided
candidate has opposite-side support and is converted to `supporting_junction`.
Conversion changes its end, not the start; it does not create the second
candidate. The reconstructed calls match archived R3 family, TSD, coordinates,
and left/right junction counts exactly.

The R2 source uses `record.reference_end + 1` when building the read cluster
(`references/RelocaTE2/scripts/relocaTE_insertionFinder.py`, relative to the
parent directory, lines 1359 and 1667). Pairing uses that stored exclusive
left endpoint (lines 578–601), so its left comparison coordinate is `p+1`.
The nearest right group is now the dominant group at `p+1`; the singleton at
`p` is left over instead. R2 separately subtracts one when computing TSD
geometry (`TSD_len_calculate`, starting at line 852). This distinction means
one must not globally shift R3 output or observation coordinates by one.

A controlled in-memory experiment shifts only the left coordinates presented
to the existing pairing helper, then converts them back before candidate
construction. The dominant paired flanks have zero overlap, not a one-base
TSD, and the current R3 candidate builder rejects them. Only the singleton
right candidate remains. This is a source-guided local experiment, not an
execution of the whole legacy pipeline.

## Sensitivity cost: not a free two-FP improvement

- Os2571: the singleton right alignment is low quality. After the pairing
  experiment it cannot validate a call by itself. The local result has no
  qualifying insertion, agreeing with the absence of an archived R2 call.
  R3 would lose a currently scored true detection as well as its extra FP.
- Os0980: the singleton passes quality and has supporting evidence. It becomes
  `supporting_junction` at 9006303..9006305, matching the archived R2 call.
  Here the local experiment retains the truth detection and removes the FP.

Consequently, a distance-based deduplication rule or universal increase in
support thresholds is not justified by these cases. The coordinate discrepancy
is demonstrated, but the best biological handling of the zero-overlap dominant
flanks is a separate decision. Do not report a full-benchmark improvement from
this local experiment. Both loci are Helitron-labeled; this review does not
establish behavior for other families or all zero-overlap candidates.

## Fixtures and tests

`tests/data/duplicate_loci/{cov30x_rep1,cov30x_rep2}.json` retains actual read
names, sequences, spans, quality classifications, junction family annotations,
supporting reads, local reference boundaries, stage-by-stage candidate records,
filter outcomes, archived R2/R3 calls, and both callers' local SAM records.
Together the fixtures are approximately 230 KiB. They require no benchmark
filesystem access during tests.

The complete accepted cluster spans are 24775400..24775993 and
9005979..9006564. They contain 48/65 junction observations and 1/3 unpaired
support observations, respectively. Indexed queries extend at least one full
1,000-bp clustering allowance beyond each accepted cluster boundary, preventing
truncation of a chained neighboring cluster. Paired non-junction reads are
retained in the SAM snapshot but, like production code, not counted as
independent support.

`scripts/reconstruct_duplicate_loci.py` reproduces extraction and refuses
existing output directories. It uses bounded indexed BAM queries, streaming
selection of relevant read-family rows, and the production cluster/candidate
helpers. Provenance records input table hashes, BAM index hash/size, source
and extractor hashes, and pysam version. Whole BAMs are not hashed. Full-read
filter outcomes were checked against the original indexed sidecar and recorded;
the portable tests do not pretend to rerun that filter without its BAM.

`tests/duplicate_loci_test.py` includes:

- actual breakpoint-count and cluster-completeness checks;
- current candidate/consolidation/arbitration/conversion reproductions;
- disjoint junction-read-set assertions;
- controlled legacy-coordinate and sensitivity-cost experiments;
- two strict expected failures documenting the unresolved pairing discrepancy;
- controls for a genuine one-base TSD, distinct nearby insertions, and a
  standalone one-sided candidate.

The current-output snapshots are characterization tests, not an instruction to
preserve duplicate behavior. Update them intentionally when implementing a
correction. Strict expected failures turn an unexpected pass into a test
failure so the known-issue markers cannot silently persist after a fix.

## Commands and validation

From the RelocaTE3 repository root:

```bash
.pixi/envs/default/bin/python scripts/reconstruct_duplicate_loci.py \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --outdir tests/data/duplicate_loci
.pixi/envs/default/bin/python -m pytest -q tests/duplicate_loci_test.py \
  tests/insertions_test.py::TestPairBreakpoints
git diff --check
```

Validation: **15 passed, 2 strict expected failures**, with no unexpected
failures. `git diff --check` passed. The expected failures identify the known
pairing discrepancy; this is not a claim that a correction is already complete.

Extraction completed successfully. Choose a new output directory to repeat it;
do not delete or overwrite the checked-in fixture candidates. Existing project
Python/pysam was used for library calls, not a new environment or replacement
aligner. The preceding audit had established that module loading was blocked
by the managed session's `/dev/log` restriction. Work here was limited to small
indexed regions, bounded streaming table reads, and lightweight unit tests.
No substantial analysis or alignment was run on the login node.

## Next action

Design a narrowly scoped pairing-coordinate correction with an explicit
zero-overlap policy, then evaluate it against these fixtures and existing
nearby-insertion tests. Before adopting it, use an HPC step-5 replay of stored
alignments to measure false-positive removals and true-detection losses across
the benchmark. This can test the candidate correction without repeating the
expensive TE alignment and trimming stages. Do not silently equate R2 parity
with greater biological sensitivity.
