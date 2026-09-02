# RelocaTE2 PSL gap-filter parity fix

**Date:** 2026-08-21 18:51 PDT

**Status:** implemented and unit-tested; riceTElib benchmark pending

## Purpose

Implement the highest-priority issue identified in
`notes/2026-08-21-ricetelib-relocate2-relocate3-parity-diagnosis.md` while
keeping the change isolated enough to measure in the next small benchmark.

## Change

`src/RelocaTE3/aligners.py::psl_to_sam` now rejects a BLAT PSL alignment when
any of these RelocaTE2 conditions is true:

```text
qNumInsert > 1
qBaseInsert > 3
tNumInsert > 1
tBaseInsert > 3
blockCount >= 3
```

The gate is applied before PSL-to-SAM conversion. Consequently, rejected
alignments are absent from the TE BAM and cannot participate in RelocaTE3's
per-read best-hit selection. This is the same location in the logical flow as
the filter in `../references/RelocaTE2/scripts/relocaTE_trim.py::parse_align_blat`.

This change intentionally does **not** alter:

- BLAT command-line sensitivity (`-minScore=10 -tileSize=7`);
- best-hit tie-breaking;
- reconstruction of PSL `matches` from SAM/CIGAR;
- TE-hit mate recovery in genome-alignment step 4;
- the full-read false-junction filter; or
- insertion clustering/calling.

Those behaviors remain separate parity changes so the benchmark can attribute
any improvement specifically to the PSL admission gate.

## Tests

`tests/aligners_test.py::TestPslToSam` now verifies that:

- the exact R2 limits (one query gap, three query inserted bases, one target
  gap, three target inserted bases, and two blocks) are accepted; and
- exceeding each of the five limits independently rejects the alignment.

Command run from the RelocaTE3 repository root:

```bash
pixi run pytest -q \
  tests/aligners_test.py \
  tests/relocate_test.py \
  tests/trim_test.py \
  tests/trim_reverse_strand_test.py
```

Result: all executed tests passed. One pre-existing BLAT executable integration
test was skipped because `blat` was not on the login-node `PATH`; the PSL
conversion and new filter tests do not require the BLAT executable and ran
successfully.

No full benchmark was run on the login node.

## Expected benchmark signal

The pre-fix riceTElib `cov5x_rep1` artifacts contained 101,226 RelocaTE3-only
classified reads, of which 101,198 (99.97%) had a selected alignment violating
at least one of these conditions. That population included 9,571 junction and
91,655 middle-read classifications.

The next isolated benchmark should therefore compare at least:

1. classified read names and junction/middle counts after step 3;
2. total, tagged, untagged, and paired records in the step-4 genome BAM;
3. raw and characterized call counts;
4. matched calls, false positives, recall, and precision; and
5. runtime for the BLAT/trim stage.

The counts need not fall by exactly the pre-fix diagnostic totals because a
read can have several PSL alignments: rejecting its former best hit may expose
another admissible hit. The critical check is that rejected alignments never
become the selected TE match.

## Benchmark handoff

Use only the BLAT + `bwa aln` RelocaTE3 caller that resolves to this development
checkout. In the current benchmark troubleshooting configuration,
`config/benchmark.bisect.toml`, that caller is
`relocate3-blat-bwaaln-mainonly`. The configuration must be supplied explicitly
because the benchmark submission wrapper otherwise uses `config/benchmark.toml`.

A suitable single-sample submission shape is:

```bash
CONFIG=config/benchmark.bisect.toml \
  bash pipeline/submit_benchmark.sh \
  --dataset ricetelib \
  --caller relocate3-blat-bwaaln-mainonly \
  --sample cov5x_rep1 \
  --no-aggregate
```

Before submission, use a fresh run/report destination or deliberately archive
the previous task outputs and completion sentinel. The caller is idempotent and
will otherwise reuse completed output. Confirm that `sbatch` returns a numeric
job ID and record the job ID, RelocaTE3 commit/dirty state, configuration hash,
and resulting SLURM log with the benchmark results.

## Next step

Do not start the TE-hit mate-recovery fix until this isolated benchmark is
scored and its stage counts are recorded. If the step-3 population contracts as
expected but the final false-positive burden remains high, proceed to the
second primary issue: port RelocaTE2's distinction between any admitted TE hit
and a classified junction/middle read into step 4.
