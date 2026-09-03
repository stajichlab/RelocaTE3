# Inferred-TSD junction validation parity fix

- Date: 2026-08-29 21:14 PDT
- Branch: `r2-parity-work`
- Status: implemented and locally validated; benchmark rerun pending

## Purpose

Remove the final RelocaTE3-only false positive in the riceTElib `cov5x_rep1`
comparison at `Chr1:22171544` by reproducing RelocaTE2's per-junction
validation after it infers a target-site duplication (TSD) length.

## Logic

The candidate pairs breakpoints at `22171544` and `22171583`, implying a 40 bp
TSD. Its left junction flank is only 29 bp, while its right flank is 88 bp.
RelocaTE2 applies a 40-base wildcard TSD pattern to each flank independently;
the left flank cannot match and is omitted. The surviving one-sided right
junction has an original full read that maps through the breakpoint, so the
existing full-read false-junction filter removes the call.

RelocaTE3 previously converted the implausibly long TSD to `UNK` but retained
both junctions. It now removes too-short junctions before calculating family,
orientation, read names, and side counts. One-sided TSD inference remains
unchanged because RelocaTE2 does not apply this paired-TSD validation to its
one-sided subclusters.

## Validation

Commands run from the RelocaTE3 repository root:

```bash
pixi run pytest -q tests/insertions_tsd_parity_test.py tests/false_junction_parity_test.py
pixi run ruff check src/RelocaTE3/insertions.py tests/insertions_tsd_parity_test.py
pixi run pytest -q
```

All tests passed. Two tests were skipped because BLAT was not present on the
login-node `PATH`. A read-only replay using the completed benchmark BAMs
produced a one-sided `Os0029#MITE/Stow` candidate (`L:0`, `R:1`) and confirmed
that `_fullread_false_junction` rejects it.

The first Ruff invocation exposed imports below test functions in the edited
test module. The imports were moved to the module header and Ruff then passed.

## Next step

Delete only the `cov5x_rep1` RelocaTE3 run and per-sample scoring directories,
then resubmit that riceTElib sample through SLURM. RelocaTE2 can retain its
completed output and will be skipped by its sentinel.
