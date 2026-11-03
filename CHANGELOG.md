# Changelog

## Unreleased

### Added

- Nextflow pipeline, milestone 1:
  - `main.nf` in nf-core layout, with parameter and sample-sheet checks
    (nf-schema) and `pixi` and `test` profiles.
  - One task per `run-all` step; the genome is indexed once; lanes are merged;
    genotyping runs in parallel.
  - Per-sample output is identical to `run-all` on the Chr3 test data
    (`tests/nextflow/parity_chr3.sh`).
  - Needs Nextflow ≥ 25.10.
- `relocaTE3 align-reads`: the genotyping alignment as its own step. `run-all`
  uses it.
- New pixi environment `nf` (Nextflow, nf-test).

### Changed

- `index-genome --genome-aligner` (default `bwaaln`) also builds that aligner's
  genome index. `run-all` now always runs `index-genome`, which skips indexes
  already present.

## 0.2.0 — unreleased (release candidates `v0.2.0rc1` 2026-10-21, `v0.2.0rc2` 2026-10-29)

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
- The missing-BLAT error no longer suggests `bowtie2` as a near-equal
  substitute.
- BLAT TE search and trimming run one chunk at a time, as in RelocaTE2's
  `--split`: peak memory at 30x on riceTElib fell from 20.1 GB to 3.4 GB, and
  genome placement from 8.9 GB to 3.0 GB, with identical outputs (rc2).
- `--tsd` accepts only `UNK` (rc2) (also `UKN`/`unknown`); other values exit with an
  explanation. The motif (`TTA`) and wildcard (`...`) modes are removed: RelocaTE2
  always runs `UNK` (`relocaTE2.py:346`) and its fixed-TSD code cannot run, so
  they had no RelocaTE2 behavior to match. `UNK` still reports each TSD.

- README "RelocaTE2 defaults" notes corrected:
  - `-s/--size` is implemented (default 500).
  - `--min-mapq` defaults to 0.
  - `--require-both-junctions` is RelocaTE2's own `characterizer.pl:91` gate,
    not a divergence.
  - A test now checks the notes against the CLI defaults.

### Added

- `relocaTE3 run-all` (whole pipeline, including reference insertions and
  genotyping) and batch mode over a FASTQ directory or sample sheet.
- Selectable aligners: `minimap2`, `bwa`, `bwa-mem2`, `bwa aln`, `bowtie2`,
  `blat`.
- CRAM input for genotyping; TSD inference for unknown TSDs.

### Known limitations

- BLAT trimming needs about 6–7 GB of temporary disk per read side at 30x
  (under the output directory; deleted when the side finishes).
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
