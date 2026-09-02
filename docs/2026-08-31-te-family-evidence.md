# Primary TE family with transparent evidence

- Date: 2026-08-31
- Branch: `r2-parity-work`
- Pull request: none
- Status: implemented and benchmark validated; extended with separate
  supporting-mate evidence on 2026-09-01

## Purpose

Keep one biologically interpretable primary TE family per insertion while
exposing conflicting junction-read assignments. This deliberately does not
reproduce RelocaTE2's slash-joined compound labels, which combine unrelated TE
groups at 13 shared `cov5x_rep1` sites.

At the eleven true-positive sites among those thirteen, RelocaTE3's primary
family exactly matches the simulated source family. The remaining two sites are
false positives. Replacing the primary label with RelocaTE2's compound string
would therefore improve literal string parity while reducing biological
specificity.

## Output fields

The existing primary family remains in `te_name`, `Name`, or `TE`, depending on
the output format. Three append-only fields carry its junction-read evidence:

- `TE_family_support`: deterministic comma-separated `family=count` votes.
- `TE_family_confidence`: primary-family votes divided by all informative
  family votes.
- `TE_family_status`: `unique`, `dominant`, `ambiguous`, or `unassigned`.

`dominant` requires an absolute majority. A deterministic primary label is
still retained when no family has a majority so legacy consumers continue to
receive one family, but `TE_family_status=ambiguous` makes the limitation
explicit.

Supporting-mate family assignments are weaker, indirect evidence and are now
reported separately through `TE_supporting_family_*` and
`TE_family_concordance`; they do not change the primary family.

The metadata is propagated through variable- and fixed-TSD step-5 tables,
tiered GFF outputs, step-7 characterized TXT/GFF, and the object-based Python
TXT/GFF writers. Existing columns are not reordered, and older tables without
the appended fields remain readable.

## Existing benchmark replay

A read-only replay against the completed riceTElib `cov5x_rep1` intermediates
recovered 12 of the 13 inspected compound-label coordinates before downstream
arbitration. Eleven had one family with confidence 1.0. At `Chr1:27848190`, the
correct Copia primary received two of four votes while two different MITE
families received one vote each; this is now explicitly `ambiguous` with
confidence 0.5 instead of becoming a compound name.

## Validation

Commands are run from the RelocaTE3 project root:

```bash
pixi run ruff check src/RelocaTE3/models.py src/RelocaTE3/insertions.py \
  src/RelocaTE3/characterize.py tests/te_family_determinism_test.py \
  tests/insertions_tsd_parity_test.py tests/characterize_test.py \
  tests/output_tiers_test.py tests/te_family_metadata_test.py
PATH="$PWD/.pixi/envs/default/bin:$PATH" .pixi/envs/default/bin/pytest -q \
  tests/te_family_determinism_test.py tests/insertions_tsd_parity_test.py \
  tests/characterize_test.py tests/output_tiers_test.py \
  tests/te_family_metadata_test.py tests/false_junction_parity_test.py \
  tests/tsd_plausibility_test.py
PATH="$PWD/.pixi/envs/default/bin:$PATH" .pixi/envs/default/bin/pytest -q
pixi run docs
```

Ruff reported that all checks passed. Pixi remained alive after printing the
successful result and was interrupted during cleanup. All 63 focused tests
passed when pytest was run directly from the pinned pixi environment. The
complete suite also passed with two expected BLAT-dependent skips.

The Sphinx HTML build completed successfully. It reported four non-code
warnings: the optional `_static` directory is absent, and the three external
intersphinx inventories could not be reached from the restricted compute
environment. Pixi again remained alive after reporting build success and was
interrupted during cleanup.

## Next step

Rerun riceTElib `cov5x_rep1` through SLURM for the supporting-mate metadata
extension and confirm unchanged detection and primary-family calls.
