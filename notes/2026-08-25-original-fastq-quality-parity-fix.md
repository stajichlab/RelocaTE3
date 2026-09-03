# Original FASTQ quality restoration for BLAT parity

Date: 2026-08-25 16:08 PDT

## Purpose

RelocaTE2 obtains sequence and base qualities from the original FASTQ while it
uses BLAT output only for TE-alignment coordinates. RelocaTE3 previously rebuilt
the reads from a PSL-derived BAM. PSL contains no quality scores, so every BLAT
read received a synthetic uniform quality score before the trimmed flanks were
aligned to the genome with `bwa aln`.

On the riceTElib 5x replicate, approximately 2,500 of 50,400 junction reads had
different genomic placements between RelocaTE2 and RelocaTE3 even though their
alignment-record totals were nearly identical. The synthetic qualities were a
demonstrated upstream difference that could cause those placement changes.

## Change

- `RelocaTE.identify_TE_reads` now passes the original FASTQ files into the trim
  stage.
- The trim stage streams each corresponding FASTQ and restores the selected
  TE-hit read's original sequence and quality before applying the alignment
  coordinates.
- Paired reads are matched using the source side, so `/1` and `/2`, Illumina
  comment-field names, and unsuffixed SRA-style names remain compatible.
- The standalone `relocaTE3 trim` command accepts `--fastq` inputs in the same
  order as `--bam`, preserving the staged interface needed by workflow engines.
- Calls that omit the original FASTQ remain backward compatible and retain the
  BAM-derived behavior.

The implementation streams the FASTQ and only modifies reads already selected
from the TE-alignment BAM. It does not build another full-library read mapping.

## Validation

Command run from the RelocaTE3 project root:

```bash
PATH="$(pwd)/.pixi/envs/blat/bin:$PATH" pytest -q tests
```

Status: 228 tests collected; 226 passed and two optional `seqtk` tests were
skipped because `seqtk` was unavailable. `git diff --check` and Python
byte-compilation also passed.

Focused regression coverage verifies that:

- a reverse-strand BLAT-like BAM read receives the correctly sliced original
  flank quality;
- `ContainingReads` retains the complete original quality;
- an unsuffixed paired FASTQ name resolves to the correct canonical mate; and
- the standalone `trim --fastq` arguments are accepted.

## Next step

Rerun `riceTElib/cov5x_rep1` for `relocate3-blat-bwaaln`, retain the existing
RelocaTE2 result, and compare alignment placement parity and final precision.
