# Paired full-read anchoring parity fix

Date: 2026-08-29 08:44 PDT

## Purpose

Remove the riceTElib `cov5x_rep1` false-positive `supporting_junction` call at
`Chr4:14561152` without losing the five truth events recovered by the preceding
RelocaTE2 supporting-junction parity fix.

## Status

Implemented locally on branch `r2-parity-work`. Focused BWA-MEM/BWA-ALN
integration tests, the full non-acceptance suite, lint, byte-compilation, and
whitespace checks pass. A retained-artifact step-5 replay removes only the
target false positive. A fresh end-to-end SLURM benchmark is still required.

## Root cause

RelocaTE2 and RelocaTE3 both:

- assigned `cov5x_rep1:clone40:h1:s1792083413:Chr4-74248/2` to
  `Os1933#DNAauto/CACTA`;
- trimmed the same 39-base TE portion with two mismatches; and
- aligned the 111-base junction flank at `Chr4:14561152`.

The untrimmed read has two equally good genome placements. RelocaTE2 aligned
the full read together with its mate and anchored it at `Chr4:14561113`, where
it spans the proposed breakpoint and causes the false junction to be rejected.
RelocaTE3 aligned the full junction read single-end, selected the alternate
placement at `Chr4:15267684`, and therefore missed the rejection.

## Logic

- For paired input, `collect_junction_fullreads` now retrieves both original
  ends of every pair containing at least one junction read.
- `align_junction_fullreads` maps the resulting R1/R2 FASTQs as pairs. Single-end
  input retains the existing single-end behavior.
- Full-read lookup reconstructs `/1` or `/2` from paired BAM flags because BWA
  removes the suffix from paired SAM query names.
- The lookup retains a suffix-free fallback for existing single-end sidecars,
  but does not let a spanning non-junction mate stand in for the junction end.
- Missing original pair ends raise an explicit error instead of silently
  weakening false-junction filtering.

The changes remain behind the shared step-4 full-read evidence interface and
the shared step-5 false-junction filter.

## Validation

Commands run from the RelocaTE3 repository root:

```bash
module load bwa/0.7.17 minimap2/2.24
pixi run pytest -q \
  tests/genome_align_test.py \
  tests/false_junction_parity_test.py \
  tests/insertions_test.py

module load bwa/0.7.17
pixi run pytest -q tests/genome_align_test.py \
  -k anchors_ambiguous -vv

PATH=".pixi/envs/default/bin:$PATH" \
  .pixi/envs/default/bin/python -m pytest -q \
  --ignore=tests/acceptance_test.py

PATH=".pixi/envs/default/bin:$PATH" \
  .pixi/envs/default/bin/python -m compileall -q \
  src/RelocaTE3/genome_align.py \
  src/RelocaTE3/insertions.py \
  tests/genome_align_test.py \
  tests/false_junction_parity_test.py

PATH=".pixi/envs/default/bin:$PATH" \
  pixi run ruff check \
  src/RelocaTE3/genome_align.py \
  src/RelocaTE3/insertions.py \
  tests/genome_align_test.py \
  tests/false_junction_parity_test.py

git diff --check
```

The ambiguity regression passes with both `bwa` and `bwaaln`. The full
non-acceptance suite collected 233 tests: 232 passed and one optional BLAT test
was skipped because BLAT was unavailable. Ruff, byte-compilation, and
whitespace checks passed.

## Retained-artifact replay

Step 5 was replayed into
`/tmp/relocate3-paired-fullread-replay.UCgOtm`; benchmark outputs were not
changed. The retained all-read paired BAM was used as a proxy for the corrected
paired junction-only sidecar so the existing expensive upstream stages did not
run on the login node.

The replay produced 190 calls. Relative to the fresh 191-call benchmark, the
only removed row was:

```text
Os1933#DNAauto/CACTA  supporting_junction  Chr4  14561152..14561154
```

All five recovered RelocaTE2-only truth events remain. If characterization is
otherwise unchanged, the expected fresh score is:

```text
total_calls              190
matched_calls            154
false_positive_calls      36
overall_precision          0.810526315789474
recall                     0.308
F1                         0.446376811594203
```

RelocaTE3 would then differ from RelocaTE2 by one additional true event
(`TE000328`) and one additional false-positive site (`Chr1:22171544`).

## Failures and limitations

- The restricted shell's HPCC module hook reported `/dev/log` permission errors
  and did not retain the requested tool paths. The installed project pixi
  environment was used as the documented fallback after module loading was
  attempted.
- The first broad suite command used a marker expression, but the acceptance
  file is not marked and was therefore included; it was rerun with an explicit
  `--ignore` and passed.
- A one-pair BWA-MEM fixture is paired and mate-anchored but is not flagged as a
  proper pair because BWA cannot estimate a library distribution from one
  pair. The regression asserts the relevant invariant: paired placement at the
  mate-selected repeat copy.
- The replay uses the retained paired all-read alignment rather than a newly
  generated paired junction-only sidecar. The fresh benchmark must verify the
  final end-to-end count and characterization.

## Next step

Delete only the existing RelocaTE3 `cov5x_rep1` benchmark output, submit the
same `relocate3-blat-bwaaln` riceTElib SLURM job, and retain the RelocaTE2
baseline for comparison.
