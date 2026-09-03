# Junction-majority TSD selection

- Date: 2026-08-30
- Branch: `r2-parity-work`
- Status: implemented and locally validated; benchmark rerun pending

## Purpose

Match RelocaTE2's read-derived TSD selection and correct the shared riceTElib
call at `Chr1:9554116`, where the truth and three of four junction reads support
`GCA` but RelocaTE3 previously reported `TCA`.

## Logic

RelocaTE2 counts the TSD sequence captured by every junction read and selects
the most-supported sequence. RelocaTE3 previously returned the first valid
right-junction capture. It now collects captures from all valid right and left
junctions and returns the most frequent sequence.

Right reads remain ordered before left reads, so tied sequences retain the
previous deterministic first-capture result. The existing 20 bp plausible-TSD
cap is unchanged; longer inferred spans still report `UNK`.

## Validation

Focused tests cover the observed 3-to-1 `GCA` vote, tied sequences, the TSD
length cap, and the existing inferred-TSD behavior.

Commands run from the RelocaTE3 project root:

```bash
pixi run ruff check src/RelocaTE3/insertions.py tests/tsd_plausibility_test.py
.pixi/envs/default/bin/pytest -q \
  tests/tsd_plausibility_test.py \
  tests/insertions_tsd_parity_test.py \
  tests/false_junction_parity_test.py
PATH="$PWD/.pixi/envs/default/bin:$PATH" .pixi/envs/default/bin/pytest -q
git diff --check
```

Results:

- Ruff passed.
- All 31 focused tests passed.
- The complete suite passed with two expected skips because BLAT is provided
  by the separate `blat` pixi environment.
- Replaying the caller against the existing `cov5x_rep1` BAM and read-family
  intermediates returned `Chr1 9554116 9554118 GCA 1 3`, confirming that the
  target call changes from `TCA` to `GCA` without changing its coordinates or
  junction counts.
- `git diff --check` passed.

The first chained `pixi run` validation attempt stalled while starting its
second pixi process. A direct full-suite attempt without activating or adding
the pixi environment to `PATH` then failed eight external-tool tests because
`minimap2`, `samtools`, and related programs were not visible. These were
environment-resolution failures, not assertion failures; prepending
`.pixi/envs/default/bin` produced the successful complete-suite result above.

## Next step

Rerun riceTElib `cov5x_rep1` through SLURM and confirm that exact-TSD accuracy
increases without changing call count, true matches, or false positives.
