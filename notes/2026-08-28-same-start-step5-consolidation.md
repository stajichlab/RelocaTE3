# RelocaTE2-parity pooling of same-start step-5 subcandidates

Date: 2026-08-28 09:11 PDT

## Purpose

Remove residual duplicate calls and restore the primary TE-family vote when
RelocaTE2's nearest-right breakpoint pairing creates several subclusters with
the same inferred TSD start.

## Status

Implemented and validated with focused tests plus a read-only replay of step 5
against the completed `riceTElib/cov5x_rep1` benchmark BAM. A fresh full
benchmark has not yet been submitted.

## Root cause and logic

RelocaTE2 pairs each left breakpoint with its nearest right breakpoint, then
writes every paired read back into `teInsertions[event][tsd_start]`. Pairs that
share a right breakpoint are consequently pooled into one call. The former
RelocaTE3 path returned each pair as an independent `Insertion` and sent it
directly through filtering and arbitration.

This produced two visible benchmark errors:

- `Chr1:299228` was emitted twice as the same Os0149 insertion, with `AT` and
  `ATTCC` TSDs. RelocaTE2 pools five junction-read observations and emits only
  the better-supported `ATTCC` call.
- `Chr1:27848190` was emitted as separate Os1668 and Os0909 calls. Both pairs
  reuse the Os3640 right-junction read, so RelocaTE2 counts that family twice
  and selects Os3640 as the primary family.

The new `_consolidate_same_start` module is the single seam used by both step-5
entry paths. It:

1. groups candidates only when chromosome and TSD start are identical;
2. preserves RelocaTE2's repeated read evidence across nearest-right pairs;
3. selects the TSD with the most junction evidence, preferring a resolved,
   shorter TSD on an exact evidence tie;
4. recomputes the junction-family and orientation votes from the pooled reads;
5. recomputes supporting-read counts for the selected coordinates; and
6. runs before full-read, low-quality, and cluster-arbitration filters, matching
   RelocaTE2's ordering.

Distinct TSD starts are never merged. Supporting-read TE-family suffixes are
not appended yet because RelocaTE3's supporting-read population is still wider
than RelocaTE2's at some clusters; using it prematurely would change correct
primary labels to unrelated tied families.

## Validation commands

From the RelocaTE3 repository root:

```bash
pixi run pytest -q tests/insertions_test.py

PATH=".pixi/envs/default/bin:$PATH" \
  .pixi/envs/default/bin/python -m pytest -q \
  tests --ignore=tests/acceptance_test.py

.pixi/envs/default/bin/python -m pytest -q \
  tests/insertions_test.py::TestInsertionFinder::test_tsd_unknown_pools_same_start_subcandidates \
  tests/insertions_test.py::TestInsertionFinder::test_tsd_unknown_infers_variable_length_tsd \
  tests/false_junction_parity_test.py

ruff check src/RelocaTE3/insertions.py tests/insertions_test.py
.pixi/envs/default/bin/python -m compileall -q \
  src/RelocaTE3/insertions.py tests/insertions_test.py
git diff --check -- src/RelocaTE3/insertions.py tests/insertions_test.py
```

Results: all 14 tests in `tests/insertions_test.py` passed; the focused parity
selection passed 12 tests; and the non-acceptance suite collected 228 tests,
with 227 passing and one optional BLAT test skipped because BLAT was not on the
direct environment's `PATH`. Ruff, byte-compilation, and whitespace validation
passed.

The real-data diagnostic reused the completed step-4 and junction-full-read
BAMs and wrote to a new `mktemp` directory under `/tmp`. It did not overwrite
benchmark results. The replay produced 185 calls and the expected rows:

```text
Os0149_complete#LTR/Gypsy  ATTCC  Chr1  299228..299232  T:5 R:2 L:3
Os3640_complete#LTR/Copia  TA     Chr1  27848190..27848191  T:4 R:2 L:2
```

At `Chr1:16883100`, pooling also restored RelocaTE2's junction counts from
RelocaTE3 `T:6 R:2 L:4` to `T:9 R:4 L:5`; the call identity is unchanged.

## Expected benchmark result

Characterization-dependent status still requires a fresh benchmark, but the
existing statuses at the affected calls make the deterministic scoring
projection:

- 185 total calls;
- 149 matched truth events;
- 36 false positives;
- precision 0.8054054054;
- recall 0.298;
- F1 0.4350364964.

This should recover truth event `TE000124`, eliminate the Os0149 duplicate,
and leave only one scoring-normalized RelocaTE3-only false-positive site
(`Chr1:22171544`) relative to RelocaTE2.

## Failures and limitations

- The HPCC module command could not update `PATH` because its logger call was
  rejected in the restricted shell. Validation used the already-installed
  project pixi environment instead.
- `sacct` was unavailable during the preceding benchmark comparison because
  `slurmdb1` could not be reached. This implementation turn did not submit a
  SLURM job.
- The diagnostic predicts scoring from unchanged characterization statuses; a
  fresh end-to-end benchmark is required to confirm the final metrics.

## Next step

Delete the previous RelocaTE3 `cov5x_rep1` benchmark artifacts, submit only
`relocate3-blat-bwaaln`, and compare the completed 185-call result with the
retained RelocaTE2 baseline.
