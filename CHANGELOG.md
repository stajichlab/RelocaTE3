# Changelog

## Unreleased

### Added

- Nextflow pipeline, milestone 4: the BLAT TE search spread over jobs.
  - New subcommands, which together give exactly `relocaTE3 run`'s outputs
    with BLAT:
    - `relocaTE3 te-split`: chunk FASTAs of `--chunk-reads` reads (default
      200,000, as `run`) and a `manifest.tsv`;
    - `relocaTE3 te-search-chunk`: BLAT any subset of the chunks (mates from
      the manifest);
    - `relocaTE3 te-merge-trim`: merge every chunk BAM in manifest order and
      trim. It refuses missing, extra or repeated chunk BAMs.
  - `relocaTE3 run` runs the same three steps in one process; its outputs are
    unchanged.
  - The pipeline uses them by default with `--te_aligner blat`:
    - `RELOCATE3_TE_SPLIT`;
    - `RELOCATE3_TE_SEARCH_CHUNK`, one task per `--blat_chunks_per_task`
      chunks (default 8);
    - `RELOCATE3_TE_MERGE_TRIM`, once per sample.
    The published files are identical. `--blat_scatter false` restores the
    single `RELOCATE3_TE_SEARCH` task.
  - `tests/nextflow/parity_chr3.sh` accepts extra pipeline params
    (`NXF_EXTRA_PARAMS`). CI runs the Chr3 parity gate in the image both with
    the defaults and with small chunks (forced scatter).
  - Measured on the riceTElib 30x simulated sample:
    - the TE search took **1 h 21 min** across 47 chunk jobs, against
      13 h 13 min as one job;
    - TE-search outputs are identical to the pre-M4 single-job `run`;
    - chunk jobs peaked at 0.36 GB of memory.

- Nextflow pipeline, milestone 3:
  - New `relocaTE3 cohort` writes:
    - `cohort.characTErized.txt`: every sample's characterized rows, unchanged;
    - `cohort.genotype_matrix.tsv`: exact sites × samples, cell = status,
      `.` = no call;
    - `cohort.summary_mqc.tsv`.
  - `run-batch` and the pipeline write it to `<outdir>/cohort/`.
  - New `relocaTE3 sample-stats`: per-sample read and call counts for the
    report.
  - The pipeline ends with a MultiQC report, `multiqc/relocate3_report.html`.
    New settings `--skip_multiqc` and `--multiqc_title`.
  - MultiQC 1.35 joins the image through a new pixi environment, `pipeline`
    (`blat` + MultiQC, solved together, so every tool build is identical). The
    pixi profile and the image now use `pipeline`.
  - `run-batch` builds the cohort only when it genotypes (`--genotype` or
    `--genotype-bam`); otherwise it logs "cohort skipped: run without
    --genotype" and warns if an old `<outdir>/cohort/` is left in place.
    With genotyping, a sample counts as done only when its characterized
    table exists and is not older than its nonref table.
  - `run-batch` now exits non-zero if the cohort step fails (the error and
    any failed samples are logged).
  - `run-batch` rejects the sample name `cohort` and names containing `=`; the
    pipeline's sample sheet rejects `cohort`, `multiqc` and `pipeline_info`.
  - `relocaTE3 cohort` writes its three files together (a failed run leaves the
    previous outputs unchanged). The TSD is not part of the matrix site key.
  - The MultiQC report comment no longer holds the run name and start time, so
    an unchanged `-resume` caches the report.
  - Checked on the 10 real-rice samples: the cohort holds all 5,775 calls,
    identical to the per-sample tables. The matrix has 2,847 sites, 1,056 of
    them called in at least 2 samples (`validation/real_rice/cohort_report.sh`).

- Nextflow pipeline, milestone 1:
  - `main.nf` in nf-core layout, with parameter and sample-sheet checks
    (nf-schema) and `pixi` and `test` profiles.
  - One task per `run-all` step; the genome is indexed once; lanes are merged;
    genotyping runs in parallel.
  - Per-sample output is identical to `run-all` on the Chr3 test data
    (`tests/nextflow/parity_chr3.sh`).
  - Needs Nextflow ≥ 26.04 (nf-schema 3.0.0).
  - `--help` prints the parameters and exits 0; `--help <parameter>` describes
    one.
  - Sample names, `--te_name` and `--target` are limited to `[A-Za-z0-9._-]`.
  - A sample with a `genotype_bam` is genotyped even with `--genotype false`,
    as `run-all` does with `--genotype-bam`.
  - New `--bcftools` setting (as `run-all --bcftools`).
- Nextflow pipeline, milestone 2:
  - Container image `ghcr.io/stajichlab/relocate3`, built by GitHub Actions
    from `pixi.lock` (`containers/Dockerfile`). Tags:
    - the version on release tags;
    - `dev` for `main`;
    - `pr-<n>` and `sha-<commit>`.
  - Profiles `apptainer`, `singularity` and `docker` (`--container_tag`), and
    `hpcc` for UCR SLURM (`--slurm_partition`).
  - CI runs pytest, nf-test, `nextflow lint` and the Chr3 parity check inside
    the image on every pull request.
  - Validation:
    - Chr3 parity passes with pixi, Docker and Apptainer.
    - The 10 real-rice samples run through `-profile hpcc,apptainer` are
      identical to v0.2.0rc2.
- `relocaTE3 align-reads`: the genotyping alignment as its own step. `run-all`
  uses it.
- New pixi environment `nf` (Nextflow ≥ 26.04, nf-test).

### Changed

- `index-genome --genome-aligner` (default `bwaaln`) also builds that aligner's
  genome index. `run-all` now always runs `index-genome`, which skips indexes
  already present.
- `bcftools` (needed by excision mode, `characterize -x` / `--excision`) is now
  in the pixi environments and the container image. Before, excision output
  was silently empty there.
- `biopython` and `pybedtools` are no longer dependencies (nothing used them).
  Build requirements are `setuptools>=61,<85` and `versioningit>=3,<4`.
- Pipeline `-profile test` now runs from any directory, also as
  `nextflow run stajichlab/RelocaTE3 -profile test,<engine>`.
- Pipeline `-profile pixi` stops at start if `--pixi_env` has no `relocaTE3`.
- Pipeline runs write `pipeline_info/software_environment.txt` (container
  image or pixi environment).
- `validation/real_rice/run_pipeline.sh` refuses the mutable image tags `dev`
  and `pr-<n>` unless `ALLOW_MUTABLE_TAG=1`, and records the image each run
  used in `OUTDIR/container_image.txt`.

### Fixed

- Pipeline: boolean pipeline settings given as `false` on the command line
  were treated as true. Nextflow >= 26.04 (the pipeline needs it) passes an
  untyped `--genotype false` as the text "false", which counts as true. So
  `--require_both_junctions false` still passed `--require-both-junctions`
  (different calls), `--excision false` turned excision on,
  `--skip_multiqc false` skipped MultiQC and the sample stats,
  `--help false` printed the help and exited, and `--genotype false` and
  `--blat_scatter false` had no effect. Affected every pipeline version so far
  (M1 onward); `relocaTE3` itself and settings given in a config file,
  `-params-file` or a profile were not affected. The five settings are now
  typed `Boolean` in `main.nf`'s `params {}` block; values other than
  true/false (any case) stop the run.

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
