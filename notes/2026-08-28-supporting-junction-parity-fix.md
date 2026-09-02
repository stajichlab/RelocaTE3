# RelocaTE2 `supporting_junction` parity fix

Date: 2026-08-28 14:36 PDT

## Purpose

Recover the five riceTElib `cov5x_rep1` truth events still detected only by
RelocaTE2 after the full-read and same-start parity fixes, without reopening
RelocaTE3's large permissive one-sided call population.

## Status

Implemented locally on branch `r2-parity-work`. Focused and related tests pass,
the full non-acceptance suite passes, and a read-only step-5 plus
characterize/scorer replay against the retained benchmark BAM confirms the
expected result. A fresh end-to-end SLURM benchmark has not yet been submitted.

## Root cause

All five missing events were RelocaTE2 `supporting_junction` calls. This is its
only retained one-sided class: one junction side is observed and a supporting
read brackets the side with no junction.

RelocaTE3 deliberately postponed this exception because its supporting counts
were much wider. The remaining implementation difference was in
`align_process`: RelocaTE2 adds a non-junction alignment to
`teSupportingReads` only when `record.is_paired` is false. RelocaTE3 added every
non-junction record, including paired mates of junction reads, causing most
one-sided candidates to appear bracketed.

## Logic

- `_stream_clusters` now admits only unpaired non-junction records as support,
  matching RelocaTE2.
- `_as_supporting_junction` is the shared classification seam used by both
  step-5 entry paths. It requires support on the side missing a junction.
- The sentinel spans three bases, matching RelocaTE2's `UKN` coordinate width:
  a right junction anchors the start directly; a left junction anchors it one
  base after the aligned flank.
- Strict `--require-both-junctions` retains this sentinel but continues to drop
  every other one-sided class.
- `ST:` now reports `SR + SL` rather than zero.

## Validation

Commands run from the RelocaTE3 repository root:

```bash
PATH=".pixi/envs/default/bin:$PATH" \
  .pixi/envs/default/bin/python -m pytest -q \
  tests/insertions_test.py \
  tests/insertions_tsd_class_parity_test.py \
  tests/false_junction_parity_test.py

PATH=".pixi/envs/default/bin:$PATH" \
  .pixi/envs/default/bin/python -m pytest -q \
  tests/insertions_test.py \
  tests/false_junction_parity_test.py \
  tests/pipeline_test.py tests/main_test.py tests/characterize_test.py

/rhome/nmath020/.local/bin/ruff check \
  src/RelocaTE3/insertions.py tests/insertions_test.py
.pixi/envs/default/bin/python -m compileall -q \
  src/RelocaTE3/insertions.py tests/insertions_test.py
git diff --check
```

The focused suite passed 28 tests, including symmetric left- and right-junction
sentinel anchoring; the related suite passed 59 tests. The full non-acceptance
suite collected 230 tests: 229 passed and one optional BLAT test was skipped
because BLAT was not on the direct environment's `PATH`. Ruff lint,
byte-compilation, and whitespace checks passed. Whole-file
`ruff format --check` reports that the two already-modified files would be
reformatted; no bulk formatting was applied so the existing narrow parity
diffs remain reviewable.

## Retained-BAM replay

Step 5 was replayed into
`/tmp/relocate3-supporting-junction.ojRXmb`; benchmark outputs were not changed.
It produced:

- 191 final calls: the existing 185 two-sided calls plus 6
  `supporting_junction` calls;
- 7,837 supporting-only clusters, down from 18,742;
- all five previously missing truth coordinates with RelocaTE2-equivalent
  sentinel coordinates and junction/support counts;
- one additional candidate at `Chr4:14561152`, which has no truth event.

Characterization reused the retained all-read BAM, and the benchmark scorer
measured:

```text
total_calls              191
matched_calls            154
false_positive_calls      37
overall_precision          0.806282722513089
recall                     0.308
F1                         0.445730824891462
```

This recovers `TE000155`, `TE000241`, `TE000291`, `TE000298`, and `TE000407`.
Relative to RelocaTE2 (188 calls, 153 matches, 35 false positives), RelocaTE3
now retains one additional true event (`TE000328`) and two additional false
positive sites (`Chr1:22171544` and `Chr4:14561152`).

## Failures and limitations

- The first targeted trace used a wildcard import, which excludes private
  helpers; it was rerun with explicit imports.
- The first scorer invocation inherited a relative pixi `PATH` after changing
  directories and could not resolve Python. Scoring was rerun successfully with
  `/usr/bin/python3.12`, matching the benchmark runner.
- Secondary supporting-read TE-family suffixes are not yet appended.
- The replay reused retained step-4 and characterization BAMs. Only a fresh
  end-to-end SLURM run can validate the complete benchmark.

## Next step

Delete the previous RelocaTE3 `cov5x_rep1` benchmark artifacts and submit the
same `relocate3-blat-bwaaln` riceTElib sample through SLURM. Retain the
RelocaTE2 baseline.
