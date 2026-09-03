# RelocaTE2 TE-hit mate-state parity fix

**Date:** 2026-08-21 PDT

**Status:** implemented and tested; riceTElib benchmark pending

## Purpose

Port the RelocaTE2 step-3/step-4 distinction between a classified TE read and
any read with an admitted TE hit. This is the second major parity fix after the
BLAT PSL gap gate.

## Logic implemented

Step 3 now writes `te_containing/<sample>.te_hit_names.txt`, containing one
canonical mate-qualified name for every selected, aligner-admitted TE hit. This
population includes reads that do not satisfy a junction or middle trimming
branch. `ContainingReads.fq` was restored to the same population, matching
RelocaTE2; classified reads retain their junction/middle tag in that FASTQ.

Step 4 now has one alignment-input planner whose interface consumes:

- the classified flanking FASTQs;
- the all-TE-hit name set; and
- the original R1/R2 FASTQs.

The planner distinguishes junction (J), middle (M), unclassified admitted hit
(U), and no admitted hit (N), and implements the RelocaTE2 actions:

| Pair state | Output |
|---|---|
| J/J | both trimmed junctions paired with `sampe` |
| J/M or J/U | junction single-end; mate suppressed |
| J/N | trimmed junction paired with original genomic mate |
| M/M or M/U | no step-4 evidence |
| M/N | original genomic mate single-end as support |

The rules are symmetric. Paired outputs preserve original R1/R2 ordering, so a
junction originating in R2 remains read 2 rather than being rewritten as read
1. Step 4 fails clearly when the new all-hit artifact is missing or inconsistent
instead of silently applying the previous incorrect inference.

The old `build_flank_pairs` and `recover_support_mates` implementations were
replaced rather than layered: their separate interfaces could not express the
four-state decision consistently.

## Validation

Commands run from the RelocaTE3 repository root:

```bash
pixi run pytest -q \
  tests/genome_align_test.py \
  tests/read_name_conventions_test.py \
  tests/trim_reverse_strand_test.py \
  tests/aligners_test.py \
  tests/trim_test.py \
  tests/pipeline_test.py \
  tests/main_test.py

pixi run test

pixi run ruff check \
  src/RelocaTE3/librelocate.py \
  src/RelocaTE3/genome_align.py \
  src/RelocaTE3/cli.py \
  tests/genome_align_test.py \
  tests/read_name_conventions_test.py \
  tests/trim_reverse_strand_test.py \
  tests/pipeline_test.py \
  tests/main_test.py

git diff --check
```

Results:

- focused tests passed; one BLAT executable integration test skipped;
- full repository test task passed; the BLAT acceptance and executable tests
  skipped because BLAT is not on the login-node `PATH`;
- Ruff and `git diff --check` passed.

No benchmark or computationally intensive analysis was run on the login node.

## Benchmark checkpoint

Use a fresh `riceTElib/cov5x_rep1` RelocaTE2 versus RelocaTE3 BLAT + `bwa aln`
run. First verify that step-3 classified counts remain at the post-PSL-gate
parity level. Then compare step-4 total, junction, untagged, and paired BAM
records. Before this fix the corresponding counts were:

| Population | RelocaTE2 | RelocaTE3 |
|---|---:|---:|
| Total | 97,496 | 198,226 |
| Junction | 50,413 | 50,432 |
| Untagged support/mates | 47,083 | 147,794 |
| Paired | 37,108 | 87,164 |

The expected signal is stable junction counts with a large contraction of the
untagged, paired, and total RelocaTE3 populations toward RelocaTE2.

## Deliberately deferred

The full-read false-junction control remains a separate major change so this
benchmark can isolate mate-state recovery. TE-family tie behavior, endpoint
conventions, orientation ties, insert-size wiring, and the narrow RelocaTE2
`supporting_junction` exception also remain deferred.
