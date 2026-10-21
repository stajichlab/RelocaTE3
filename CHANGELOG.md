# Changelog

## 0.2.0 — unreleased (release candidate `v0.2.0rc1`, 2026-10-21)

First release aimed at replacing RelocaTE2 for paired-end short-read data.

### Validated

- RelocaTE2 parity on the simulated rice benchmark (mPing, riceTElib, and
  riceTElib at 0–20% sequence divergence; 72 samples, 5–30x). RelocaTE3
  keeps 99.8% or more of RelocaTE2's true detections and has slightly higher
  pooled recall in every panel, with a few more false positives in the
  multi-TE panels (+56 true positives, +33 false positives overall).
  Benchmark: `relocate-benchmark`, profile `checkpoint-213a0df`.
- Results are reproducible: a full rerun matched the previous calls in 72/72
  samples.

### Changed

- Default aligners now match RelocaTE2: `blat` for the TE-library search and
  `bwa aln` for genome placement. Each stage can be switched with
  `--te-aligner` / `--genome-aligner`.
- RelocaTE2's two-sided junction gate, TSD-length estimator and emission
  rules, read admission and `high_conf` filter are ported as in RelocaTE2.
- TE-family labels are resolved from all supporting reads at a locus.
- `--require-both-junctions` is on by default, as in RelocaTE2.
- TE-search memory reduced by streaming BLAT output; BLAT runs in chunks.

### Added

- `relocaTE3 run-all` (whole pipeline, including reference insertions and
  genotyping) and batch mode over a FASTQ directory or sample sheet.
- Selectable aligners: `minimap2`, `bwa`, `bwa-mem2`, `bwa aln`, `bowtie2`,
  `blat`.
- CRAM input for genotyping; TSD inference for unknown TSDs.

### Known limitations

- Peak memory at 30x on a multi-TE library is about 3x RelocaTE2's
  (19–20 GiB vs 6.5 GiB). Request 64 GB for comparable runs.
- Somatic insertions are detected but are rarely labelled somatic; genotype
  thresholds are RelocaTE2's count rules, not a calibrated model.
- Helitrons and elements more than ~5% diverged from the library are poorly
  detected, as in RelocaTE2.
- Validated only for paired-end short reads on rice. Single-end input,
  reference/absence calls and excision output (`characterize -x`) are not
  covered by the benchmark.
- Excision mode needs `bcftools`, which is not installed by the pixi
  environments.
- The default TE aligner, BLAT, is free for academic, nonprofit and personal
  use only, and its bioconda package is linux-64 only.

## 0.1.0 — 2026-06-12

Initial Python port of RelocaTE2.
