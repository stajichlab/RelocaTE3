# Full-read false-junction parity wiring

Date: 2026-08-26 PDT

## Purpose

RelocaTE2 rejects an insertion candidate when at least 30% of the junction
reads on each side, in their original untrimmed form, map straight through the
proposed breakpoint. RelocaTE3 already ports this rule in `insertions.py`, but
the staged CLI and benchmark did not provide its junction-only full-read BAM to
step 5.

## Status

Implemented locally on branch `r2-parity-work`. Focused and non-acceptance
tests pass; scheduler validation and the `riceTElib/cov5x_rep1` benchmark are
pending because both Slurm controllers were unavailable during the first
submission attempt.

## Logic

- Step 4 identifies junction membership from the tagged flanking FASTQs.
  `read_repeat_name.txt` cannot provide this membership because it deliberately
  stores untagged original names.
- The original versions of only those junction reads are streamed from R1/R2,
  preserving their sequence and base qualities.
- They are aligned single-end with the configured genome aligner and written to
  `genome_aln/<sample>.fullreads.genome.bam`.
- `run-all` and the benchmark pass that BAM to
  `find-insertions --fullreads-bam` before characterization.
- The separate all-read minimap2 BAM remains characterization input only.

## Expected result

A read-only diagnostic on the retained pre-fix `riceTElib/cov5x_rep1` artifacts
reconstructed all 295 RelocaTE3 calls. Applying the existing filter with
RelocaTE2's junction-only `bwa aln` BAM retained 187 calls and all 148 matched
truth events, reducing false positives from 147 to 39. A fresh benchmark must
confirm this estimate using a RelocaTE3-generated BAM.

## Validation commands

```bash
PATH="$(pwd)/.pixi/envs/blat/bin:$PATH" \
  .pixi/envs/blat/bin/python -m pytest -q \
  tests/genome_align_test.py tests/false_junction_parity_test.py tests/main_test.py

PATH="$(pwd)/.pixi/envs/blat/bin:$PATH" \
  .pixi/envs/blat/bin/python -m pytest -q tests \
  --ignore=tests/acceptance_test.py

ruff format --check \
  src/RelocaTE3/genome_align.py src/RelocaTE3/cli.py \
  tests/genome_align_test.py tests/main_test.py

ruff check \
  src/RelocaTE3/genome_align.py src/RelocaTE3/cli.py \
  tests/genome_align_test.py tests/main_test.py

git diff --check
```

Results:

- Focused parity/CLI suite: 39 passed.
- Non-acceptance suite: 227 collected, 225 passed, 2 optional `seqtk` tests
  skipped.
- Ruff format/check, byte-compilation, Bash syntax, and whitespace validation:
  passed.
- Benchmark repository suite: 112 passed.
- Full acceptance-suite Slurm submission: unconfirmed; `sbatch` returned no
  numeric ID and `scontrol ping` reported both controllers down.

## Next step

Run focused and full tests, then submit the same RelocaTE3
`riceTElib/cov5x_rep1` benchmark without rerunning the retained RelocaTE2
baseline.
