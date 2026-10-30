# `--tsd` accepts only `UNK` (2026-10-27)

Branch `fix/fixed-tsd-r2-parity`, from `release/v0.2.0` at f7c93f0.

## Why

- RelocaTE2 2.0.1 has no `--tsd` option. Its driver writes `UNK` into
  regex.txt on every run (`relocaTE2.py:346`).
- Its fixed-TSD branch cannot run. `TSD_check` takes 14 parameters but is
  called with 13 (`relocaTE_insertionFinder.py:923, 941`) or 10 (`:1603, :1611`).
- RelocaTE3's motif (`TTA`) and wildcard (`...`) modes were therefore never
  validated. They partly ported that dead code and added logic of their own:
  `_merge_offset_starts`, plus a one-sided TSD rule (4b0624a). Supporting-read
  counts were always 0.
- The first real-rice harness ran `--tsd "..."` and differed from RelocaTE2 by
  about 1,270 calls (`docs/2026-10-21-release-v0.2.0-validation.md`).

## What changed

- `insertions.require_unknown_tsd` accepts `UNK`/`UKN`/`unknown`, the
  case-insensitive regex of `relocaTE_insertionFinder.py:1337`. It raises
  `ValueError` otherwise.
- `--tsd` on `run-all` and `find-insertions` validates with it. Other values
  exit with code 2 and an explanation.
- `InsertionFinder.find_insertions` always takes the `UNK` path. These methods
  were removed: `_cluster_reads`, `_assign_cluster`, `_tsd_check`,
  `_merge_offset_starts`, `_write_output`, `_write_event_start`,
  `_insertion_family(_evidence)`.
- Tests:
  - Removed the fixed-path tests: `insertions_tsd_class_parity_test.py`,
    and in `insertions_test.py` `test_known_tsd_call`,
    `test_offset_split_merges_to_one_call` and `test_distinct_sites_not_merged`.
  - Switched `test_target_all_preserves_per_cluster_chrom` to `UNK`.
  - Added `tests/tsd_unknown_only_test.py`.

## Verify

```bash
E=<pixi blat env bin>
PATH=$E:$PATH PYTHONPATH=$PWD/src $E/python -m pytest tests --junitxml=junit.xml
```

Results:
- Before: 319 passed.
- After: 336 passed, 0 failed, 0 skipped.
- New test file, red first: 14 failed and 8 passed (the 8 already accepted `UNK`).
- `UNK` output is unchanged: that code path was not edited.

Note: tests need the env's `bin` on PATH. Without it, 8 fail with "minimap2 not found".

## Next

- Tag `v0.2.0rc2` and rerun `validation/release/fresh_clone_check.slurm`.
- Not touched here: the README "Parameter" notes are stale:
  - "`--size` not implemented": it exists.
  - `--min-mapq` default 1: it is 0.
  - "both junctions … a deliberate divergence": it is RelocaTE2 logic.
