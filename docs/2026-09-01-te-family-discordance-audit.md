# Junction versus supporting-read TE-family audit

- Date: 2026-09-01 18:57 PDT
- Project: RelocaTE3
- Branch: `r2-parity-work`
- Pull request: none
- Status: diagnostic complete; recommended metadata propagation implemented

## Purpose

Determine why RelocaTE2 reports slash-joined TE-family labels at 13 shared
riceTElib `cov5x_rep1` insertion sites while RelocaTE3 reports one primary
family. The initial hypothesis was that RelocaTE3 discarded co-optimal
TE-library alignments before insertion aggregation.

## Inputs

The audit used only retained outputs from the completed benchmark:

- RelocaTE2 and RelocaTE3 normalized call tables.
- RelocaTE2 `Chr1.repeat.reads.list` and chromosome-split
  `read_repeat_name` table.
- RelocaTE3 `read_repeat_name` table and coordinate-sorted genome BAM
  containing junction flanks and supporting mates.
- RelocaTE2's retained reference implementation in
  `../references/RelocaTE2/scripts/relocaTE_insertionFinder.py`.

No alignment or benchmark job was rerun.

## Method

For every shared site with a different normalized family label, the audit:

1. Recovered the junction and supporting read names.
2. Mapped each junction read back to its selected TE family.
3. Applied RelocaTE2's supporting-read lookup order to recover the TE family of
   the TE-containing mate.
4. Reconstructed the corresponding RelocaTE3 insertion candidate from the
   retained genome BAM.
5. Counted junction-family and supporting-family evidence separately.
6. Tested whether RelocaTE2's compound string was exactly
   `<junction primary>/<supporting primary>` and whether RelocaTE3 retained the
   same supporting primary.

The reproducible command, run from the RelocaTE3 project root, was:

```bash
PYTHONPATH=src .pixi/envs/default/bin/python \
  scripts/audit_te_family_discordance.py \
  --benchmark-root ../../relocate_benchmark/relocate-benchmark \
  --output results/te-family-discordance-audit.tsv
```

The script refuses to overwrite an existing output.

## Result

All 13 family-string differences have the same explanation:

- RelocaTE2's first family is the primary family among breakpoint-spanning
  junction reads.
- RelocaTE2's second family is the primary family among bracketing supporting
  mates.
- RelocaTE2 concatenates those two primary names when they disagree.
- RelocaTE3 retains both evidence groups, but its family metadata currently
  summarizes only the junction reads.

The result table confirms:

- `r2_compound_explained=true` at 13 of 13 sites.
- `r3_retains_r2_supporting_primary=true` at 13 of 13 sites.
- Twelve sites have unanimous junction-family assignments.
- `Chr1:27848190` is the only site with mixed junction-family evidence.

Therefore, co-optimal TE alignment loss is not the cause of this parity gap.
Preserving co-optimal hits might be useful for a separate analysis, but it
would not reproduce or explain these RelocaTE2 compound labels.

## Example: Chr1:27848190

RelocaTE2 reports:

```text
Os3640_complete#LTR/Copia/Os0765#DNAnona/CACTA
```

The evidence is:

```text
junction reads:
  Os3640_complete#LTR/Copia = 2
  Os0909#MITE/Stow         = 1
  Os1668#MITE/Stow         = 1

supporting mates:
  Os0765#DNAnona/CACTA     = 2 in RelocaTE2
  Os0765#DNAnona/CACTA     = 3 in RelocaTE3
```

RelocaTE2 selects Copia as the junction primary, selects CACTA as the
supporting primary, and joins the two names. RelocaTE3 correctly keeps Copia as
the primary family and reports the mixed junction evidence as `ambiguous`, but
does not yet expose the retained CACTA supporting evidence.

## Exact information-loss point

RelocaTE3 has not lost the supporting-read family assignment from the raw data.
The loss happens while constructing the final `Insertion` object:

- `_make_insertion` calculates `TE_family_support`, confidence, and status from
  junction observations only.
- `_count_support` counts bracketing supporting reads by side but discards their
  names and family assignments.
- The final output therefore contains `ST`, `SR`, and `SL` counts without the
  corresponding supporting-family evidence.

RelocaTE2 preserves the distinction internally, calculates separate junction
and supporting primaries, and then collapses that distinction into a compound
name. RelocaTE3 should preserve the distinction explicitly instead.

## Biological interpretation

Junction reads directly span a TE/genome breakpoint and provide the strongest
family-to-locus linkage. Supporting evidence is indirect: one mate is assigned
to a TE while its partner maps near the insertion. A disagreement between the
two evidence classes does not imply that two unrelated TE families inserted at
the same location.

Keeping the junction primary as the reported TE family is therefore the
appropriate default. The supporting-family disagreement should be exposed as
lower-tier evidence rather than appended to the primary name or allowed to
outvote breakpoint-spanning reads.

## Recommended implementation

Extend the insertion evidence model without changing the primary family or
call filtering:

1. Keep `TE_family_support`, `TE_family_confidence`, and `TE_family_status` as
   junction-read evidence and document that scope explicitly.
2. Change support counting so it also retains supporting-read names and maps
   their TE-containing mates through `read_repeat_name`.
3. Append separate supporting evidence fields, for example:
   `TE_supporting_family_support`, `TE_supporting_family_confidence`, and
   `TE_supporting_family_status`.
4. Append `TE_family_concordance` with values such as `concordant`,
   `discordant`, `junction_only`, `supporting_only`, or `unassigned`.
5. Propagate the fields through TXT, GFF3, characterization, and the modern
   object API.
6. Add tests ensuring supporting reads cannot change the primary junction
   family, each read contributes once, and old inputs without the new columns
   remain readable.
7. Rerun the riceTElib benchmark and require unchanged calls and accuracy while
   all 13 sites report the retained supporting-family disagreement.

## Failures and corrected assumptions

The initial co-optimal-alignment hypothesis was rejected by the retained source
and read-level evidence. The discrepancy is generated after read alignment by
combining two evidence classes, not by joining equally scoring alignments from
one read. No computational failure occurred during the final audit.

## Next step

The recommended metadata change is implemented. Rerun riceTElib `cov5x_rep1`
and confirm unchanged calls and primary-family metrics while all 13 sites expose
their supporting-family disagreement.
