# Separate supporting-mate TE-family evidence

- Date: 2026-09-01 19:13 PDT
- Project: RelocaTE3
- Branch: `r2-parity-work`
- Pull request: none
- Status: implemented and locally validated; benchmark rerun pending

## Purpose

Expose the supporting-mate family evidence that RelocaTE2 encodes by appending
a second family to compound labels, without allowing weaker indirect evidence
to replace the family supported by breakpoint-spanning junction reads.

The preceding read-level audit showed that all 13 riceTElib `cov5x_rep1`
family-string differences are junction-versus-supporting disagreements.
RelocaTE3 retained the supporting assignments in `read_repeat_name`, but final
calls retained only the number and side of supporting reads.

## Implementation

The primary `te_name` remains the primary junction-read family. Existing fields
now have an explicit junction-read scope:

- `TE_family_support`
- `TE_family_confidence`
- `TE_family_status`

Four append-only fields expose the indirect evidence:

- `TE_supporting_family_support`: ordered `family=count` supporting-mate votes.
- `TE_supporting_family_confidence`: the supporting primary's fraction of all
  informative supporting votes.
- `TE_supporting_family_status`: `unique`, `dominant`, `ambiguous`, or
  `unassigned`.
- `TE_family_concordance`: `concordant`, `discordant`, `junction_only`,
  `supporting_only`, or `unassigned`.

Supporting-family votes follow RelocaTE2's mate lookup order: unsuffixed name,
`/1`, `/2`, `.f`, then `.r`. Repeated genome alignments with the same read name
retain their legacy contribution to `ST`, `SR`, and `SL` counts but cast only
one family vote. Supporting-only calls retain `te_name=NA`; indirect evidence
is reported but never promoted to the primary family.

The fields are propagated through:

- Step-5 legacy insertion tables.
- Tiered GFF3 outputs.
- Structured TXT/GFF3 writers and readers.
- Step-7 characterized TXT/GFF3 outputs.
- The `Insertion` Python object.

Older GFF3 and Step-5 tables without these columns remain readable and receive
empty, zero, and `unassigned` defaults.

## Benchmark replay

A read-only replay reconstructed the 13 family-discordant shared calls from the
completed riceTElib genome BAM and `read_repeat_name` table using the updated
code. All 13:

- retained the same RelocaTE3 primary family;
- recovered the exact RelocaTE2 supporting primary;
- had `TE_supporting_family_status=unique`;
- had `TE_family_concordance=discordant`.

This replay did not rewrite benchmark outputs or measure detection metrics.

## Validation

Commands were run from the RelocaTE3 project root:

```bash
PATH="$PWD/.pixi/envs/default/bin:$PATH" \
  .pixi/envs/default/bin/pytest -q \
  tests/te_family_metadata_test.py \
  tests/supporting_reads_test.py \
  tests/characterize_test.py \
  tests/output_tiers_test.py \
  tests/insertions_tsd_parity_test.py

pixi run ruff check \
  src/RelocaTE3/models.py src/RelocaTE3/insertions.py \
  src/RelocaTE3/characterize.py tests/te_family_metadata_test.py \
  tests/supporting_reads_test.py tests/characterize_test.py \
  tests/output_tiers_test.py scripts/audit_te_family_discordance.py

PATH="$PWD/.pixi/envs/default/bin:$PATH" \
  .pixi/envs/default/bin/pytest -q

pixi run docs
```

The focused suite passed. Ruff reported that all checks passed. The complete
suite passed with two expected skips because BLAT was not available on the
login-node `PATH`. The Sphinx HTML build succeeded with four expected
environment warnings: the optional `_static` source directory is absent and
three external intersphinx inventories were unreachable from the restricted
environment.

## Failures

No code or test failures remain. Pixi stayed alive after Ruff printed its
successful result and was interrupted during cleanup; Ruff itself completed
successfully.

## Next step

Rerun riceTElib `cov5x_rep1` through SLURM. Detection metrics and primary-family
calls should remain unchanged, while all 13 formerly compound-label sites
should expose their supporting family and `discordant` concordance.
