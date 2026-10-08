# Candidate: retain TE-hit ties and resolve family using distinguishing reads

Recorded: 2026-10-01, after 15:10 PDT, America/Los_Angeles.
Project: RelocaTE3, `main`, following merged PR #49; changes remain uncommitted.
Status: implementation and 83 focused tests complete; 72-sample frozen replay
prepared. Submission timed out without a job ID and is unconfirmed. No release
pin or claim of improved full-benchmark metrics.

## Implemented behavior

The saved LINE fixture has eight reads tying between Os3328 and Os0596.
R3's alphabetical selection converted those ties into eight Os0596 votes.
The candidate retains the same selected alignment and flanking sequence,
but preserves the set of equally best TE families before that selection.

Insertion-family resolution is isolated in `src/RelocaTE3/family.py`, a pure
module taking selected families, best-family sets, and geometry compatibility.
Both the fixed-TSD and variable-TSD calling paths use this module. This small
interface concentrates the evidence rule and its tests, following the
codebase-design skill; the same rule can later be implemented in Rust.

A primary family can resolve tied evidence only when:

1. At least two distinct junction-read names have that single best family.
2. That family has a strict majority among the distinguishing reads.
3. It appears in the best-family sets of a strict majority of all informative
   junction reads at the insertion.
4. Tied alignments have compatible query start/end, strand, mismatch count,
   and TE-end classification, preserving the selected flank geometry.

Otherwise, the deterministic selected-vote primary remains, with explicit
ambiguity. Supporting mates do not resolve or outvote junction families.
Only one primary family is reported. Counts are read evidence, not independent
molecule counts: no new PCR/UMI or mate-pair deduplication model was introduced.
The two-read quorum is a conservative candidate policy, not a biologically
calibrated threshold or established full-benchmark improvement.

When tied-family evidence is present, duplicate read names cannot provide
multiple distinguishing votes in candidate consolidation. Legacy mappings
without tied evidence keep their existing selected-vote behavior, including
their historical consolidation counts.

## Result on the real fixture

The integrated parser → mapping → clustering → insertion test now reports:

- Primary: **Os3328**, matching simulated truth.
- Coordinates: **Chr1:29402710..29402725**, unchanged.
- TSD: **GGAGAGGGAGGTGGCC**, unchanged.
- Junction counts: **L=2, R=9**, unchanged.
- Distinguishing-family support: **Os3328=2, Os0596=1**.
- Best-family candidate support: **Os3328=10, Os0596=9**; these are overlapping
  memberships, not mutually exclusive votes, and their sum can exceed 11.
- Ambiguous reads: **8**; resolution: `unique_junction_consensus`.
- Status: `ambiguous`; confidence: **2/11 ≈0.181818**, the fraction of all
  informative reads independently distinguishing the selected family.

This confidence is not a calibrated probability. It deliberately does not
convert tied reads into strong independent support. Replaying the same
genomic reads with their original three-column mapping still produces Os0596.

## File contracts and compatibility

The first three columns of `read_repeat_name.txt` remain read, selected TE,
and strand. Multi-family ties append one column:

```text
TE_best_families:{"families":["Os0596#LINE/unknown","Os3328#LINE/unknown"],"compatible":true}
```

Unambiguous reads retain the original three-column layout. Older mapping
columns such as R2 chromosome annotations continue to be ignored unless they
have this explicit prefix. Malformed new evidence fails rather than silently
becoming a confident selected vote. A better-scoring hit clears prior weaker
ties; record order does not change the selected alignment or family set.

When ties exist, `TE_family_support` counts distinguishing reads and
`TE_family_confidence` divides distinguishing support by all informative
reads. Legacy evidence keeps its prior support/confidence semantics. Three
append-only fields explain the richer evidence:

- `TE_family_candidate_support`: best-family set membership counts;
- `TE_family_ambiguous_reads`: number of reads with multiple best families;
- `TE_family_resolution`: `selected_votes`, `unresolved_ties`,
  `unique_junction_consensus`, or `unassigned`.

These fields survive raw/fixed-TSD tables, structured insertion TXT/GFF,
insertion-tier GFF conversion, GFF reading, and both characterization writers.
Existing normalized benchmark calls still use the first eight characterization
columns. Cached old mapping files lack the discarded alternatives; rerunning
step 5 alone with those files will not exercise this change.

## Verification

**83 focused tests pass** across:

```bash
.pixi/envs/default/bin/python -m pytest -q \
  tests/family_resolution_test.py \
  tests/te_family_determinism_test.py tests/te_family_metadata_test.py \
  tests/line_family_evidence_test.py tests/trim_reverse_strand_test.py \
  tests/insertions_test.py::TestInsertionFinder tests/characterize_test.py \
  tests/duplicate_loci_test.py tests/insertions_tsd_parity_test.py \
  tests/insertions_tsd_class_parity_test.py tests/memory_optimization_test.py
```

The new regression tests cover the integrated real LINE fixture, all-tied
evidence, a single distinguishing read, conflicting distinguishing families,
incompatible geometry, a candidate absent from most best-family sets, missing
families, duplicate read names, better-score reset, hit-order independence,
malformed mapping evidence, output round trips, and replay read filtering.

Lint on the changed production modules and new code, formatting, shell syntax,
and `git diff --check` pass. Initial checks caught outdated trailing-column
assertions and a legacy consolidation-count regression; both were corrected.
An early broad test selection also included two external-aligner tests that
failed because minimap2 was absent from the session PATH. They were not rerun
on the login node. The final 83 passing checks are focused library tests,
not the entire test suite or a fresh alignment benchmark.

The portable genome fixture is
`tests/data/line_family/TE000172.junctions.json`. It retains eleven archived
genomic SAM records and its extraction hash. The test provides a synthetic
contig header; it does not claim that header is the original genome header.

## Frozen validation across all three panels

Prepared output: `results/family-replay/2026-10-01/`, containing 72 tasks:
nine mPing, nine riceTElib, and 54 diverged riceTElib. The manifest verifies
27 frozen source/runner/scoring files by SHA-256 and records input metadata.

Each compute task:

1. Collects genome-mapped junction names and scans stored TE BAMs for their
   best-family sets using the production parser with a read-name filter.
2. Requires its selected family and strand to reproduce the original mapping
   before appending tied evidence; missing hits/mappings fail explicitly.
3. Calls and characterizes both original and enriched evidence with **identical
   frozen current source**. Thus the difference isolates tied evidence from
   the already pending breakpoint-pairing change.
4. Scores both variants and cached R2 calls with the same `tsd-interval`
   policy and 10-bp tolerance. It saves event-level matches, FPs, correctness,
   precision, normalized calls, and completion markers.
5. Records whether original evidence reproduces historical R3 calls, and
   whether geometry/TSD/strand/genotype remain unchanged between variants.
   Historical disagreement is reported explicitly, not assumed away.

The replay checks all three panels without repeating expensive alignments.
It exercises new family evidence for junction reads; it does not measure the
memory/runtime cost of full trimming across every original TE-hit read. That
cost remains a release-validation requirement. Per-variant elapsed time here
is step 5 plus characterization, not whole-pipeline runtime.

The epyc array requests one CPU, 32 GB, and eight hours per task, with four
tasks concurrently. Cached cluster configuration gives epyc MaxTime=30 days.
The script captures `SLURM_SUBMIT_DIR` and has no `--chdir` directive.
It uses the existing project Python/pysam; no bioinformatics executable or
new environment is required. Large scans are guarded by `SLURM_JOB_ID`.

Preparation command (already completed; do not overwrite the frozen output):

```bash
.pixi/envs/default/bin/python scripts/replay_family_resolution.py prepare \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --output results/family-replay/2026-10-01
```

Submission was attempted with a bounded timeout and returned exit 124 without
a job ID. No second submission was attempted. Check the queue for an existing
`replay_family_resolution.slurm` job before resubmitting from RelocaTE3 root:

```bash
squeue -u "$USER" -o '%.18i %.45j %.10T'
```

If no matching job exists and the frozen output has no task starts, submit:

```bash
sbatch --array=0-71%4 \
  results/family-replay/2026-10-01/replay_family_resolution.slurm \
  results/family-replay/2026-10-01
```

After completion, compare R2, original evidence, and tied evidence across
coverage, divergence, TE group, and somatic fractions, explicitly tracking
lost as well as gained detections. Do not pin a release or infer a net
benchmark improvement from the one corrected LINE fixture.
