# Nextflow M3 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** after all samples finish, the pipeline writes:
- a cohort call table;
- a sites × samples genotype matrix;
- a MultiQC summary report.

The cohort rows must match the per-sample files exactly. Design:
`plans/2026-11-02-nextflow-pipeline-design.md` (milestone M3, CLI addition 3).

**Architecture:**
- The logic lives in a new tested Python subcommand, `relocaTE3 cohort`. The
  workflow contains none.
- `run-batch` calls the same subcommand, so non-Nextflow users get identical
  cohort files.
- MultiQC renders "custom content" tables that `cohort` and the per-sample
  tasks write. It runs only on summary files and never touches calls.

**Branch / worktree:** `feat/nextflow-m3` in `../RelocaTE3-nextflow-m3`
(created with `--no-track`), from `origin/main` 583cd59, which merges PR #64.
The maintainer pushes with `git push -u origin feat/nextflow-m3`.

---

## Conventions (as M1/M2)

**Test commands:**
- pytest: `E=$PWD/.pixi/envs/blat/bin; PATH=$E:$PATH PYTHONPATH=$PWD/src $E/python -m pytest tests -p no:cacheprovider --junitxml=<scratch>/x.xml`.
  Read the counts from the XML.
- `pixi run -e nf nf-test test`
- `pixi run -e nf nextflow lint .`
- actionlint, if workflows change.
- **Baseline at 583cd59:** pytest 533, nf-test 30. Record it in Task 0.

**Process:**
- TDD.
- sbatch for anything heavy: `#!/usr/bin/bash -l`, `set -euo pipefail`,
  `-p intel`, `logs/<name>.%j.log`.
- Never on the login node.
- Commit per task with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Never push.

**Real data for checks:** M2's real-rice pipeline output is at
`../RelocaTE3/results/nextflow-m2/nextflow-real-rice/dev/<sample>/` (10
samples, mPing). The characterized table has a header plus one row per call,
e.g. G_1 has 623 calls. Columns:
`strain TE TSD chromosome.pos strand avg_flankers spanners status` and then the
`TE_family_*` columns. `chromosome.pos` looks like `Chr1:1132975..1132977`.
`status` values include `homozygous`, `heterozygous`, `somatic` and combined
forms such as `homozygous/excision_no_footprint`.

---

### Task 0: Setup
- `pixi install -e blat -e nf`.
- Run the baseline suites and record the counts.

### Task 1: `relocaTE3 cohort` subcommand

**Files:**
- `src/RelocaTE3/cohort.py` (new)
- `src/RelocaTE3/cli.py` (`_menu_cohort`, `cmd_cohort`, registration)
- `tests/cohort_test.py`

**Interface:**
`relocaTE3 cohort -i <characTErized.txt> [<...>] -o <outdir> [--te-name NAME]`.
Inputs may be given in any order. The output is sorted by sample name, so it
is deterministic and independent of the order Nextflow collects files.

**Outputs, in `<outdir>/`:**
1. **`cohort.characTErized.txt`:** one header (the per-sample header,
   unchanged), then every input row unchanged, ordered by sample (`strain`)
   and then by original row order within the sample.
   - The rows equal the concatenation of the per-sample files' data rows. Test
     this as a multiset and in order.
   - Fail if headers differ between inputs, or if two inputs claim the same
     sample.
2. **`cohort.genotype_matrix.tsv`:**
   - Columns: `chromosome start end strand TE` and then one column per sample,
     sorted.
   - One row per **exact** site key `(chromosome, start, end, strand, TE)`,
     parsed from `chromosome.pos` and the `TE`/`strand` columns.
   - Cell = that sample's `status`, or `.` when the sample has no call there.
   - If one sample has two rows with the same key, fail loudly. Check first
     whether this happens in real data. If it does legitimately, decide and
     document a rule: joining with `;` is a candidate.
   - Rows are sorted by chromosome in natural order (Chr1 < Chr2 < Chr10),
     then start, end, strand and TE.
   - The header documents that nearby sites are **not** merged, a decision
     deferred to the VCF design.
3. **`cohort.summary_mqc.tsv`:** MultiQC custom content, a general-stats-style
   table. One row per sample: total characterized calls and counts per status
   class. Split combined statuses on `/` and count the first part as the
   genotype class, e.g. `homozygous/excision_no_footprint` counts as
   homozygous plus an "excision" column. Document the rule.
   - Use MultiQC's documented `_mqc.tsv` header-comment format (`# id:`,
     `# section_name:`, `# plot_type: 'table'`, …). Check it against the
     MultiQC custom content docs, for the MultiQC version pinned in Task 3.

**Tests:**
- Small synthetic tables covering ordering, header mismatch, duplicate
  sample, duplicate site in one sample, `.` cells and natural chromosome sort.
- A test on the real-rice outputs, skipped if the directory is missing.
  `cohort.characTErized.txt` has 1 + Σ(per-sample data rows) lines, and each
  sample's rows equal its file's rows exactly. Report the matrix row count.

Commit: "relocaTE3 cohort: cohort call table, genotype matrix, MultiQC summary".

### Task 2: `run-batch` writes the cohort too

- After all samples succeed, or with `--keep-going` after the run, using the
  succeeded samples, `cmd_run_batch` calls `cmd_cohort` on every sample's
  `results/<target>.<te_name>.all_nonref_insert.characTErized.txt`. It writes
  to `<outdir>/cohort/`.
- This happens only when genotyping produced characterized tables. Otherwise
  it logs that the cohort was skipped.
- Delegation test in the style of `tests/main_test.py`: `rec()` monkeypatch.

Commit: "run-batch writes the cohort table via relocaTE3 cohort".

### Task 3: MultiQC availability

**Decision for the implementer (record it in the docs):**
- **Preferred:** add `multiqc` to the image's environment, so the one-image
  design stays intact. The blat env pins zlib 1.2.11 for bioconda blat, and
  MultiQC's dependencies may not co-solve with it.
  1. Try a new pixi feature `report` with `multiqc` (pinned, e.g. `1.*`
     current), plus an environment `pipeline = blat + report`.
  2. If it solves with `--locked` and does not change any tool version in the
     blat env (diff `pixi list -e blat` against `pixi list -e pipeline` for
     blat, bwa, samtools, bcftools, minimap2, seqtk, bedtools, bowtie2,
     pysam, python), build the image from `pipeline` instead of `blat`.
     Update `containers/build_env.sh`, `Dockerfile` (`ENV_DIR`), `check_image.sh`
     (add multiqc), `dry_run_recipe.sh` and the tests.
  3. Keep `blat` as the pixi profile default for tasks. The pixi profile's
     `pixi_env` can stay blat, with MultiQC from the `pipeline` env, or the
     profile points at `pipeline`. Choose the simpler option and document it.
- **Fallback, only if (1) fails or changes tool versions:** a standalone pixi
  env `report` for MultiQC, and in container profiles a per-process container
  for the MultiQC task, the official
  `quay.io/biocontainers/multiqc:<pinned version>` (verify the tag exists).
  Document the deviation from one image as reporting-only.
- Run the container recipe dry run (`sbatch containers/dry_run_recipe.sh`). It
  must PASS.

Commit: "MultiQC in the pipeline environment" (or the fallback wording).

### Task 4: Per-sample counts for the report
- Each per-sample module, or a new tiny `RELOCATE3_SAMPLE_STATS` process,
  emits `<id>.stats_mqc.tsv` with cheap counts from files the steps already
  produce:
  - TE-containing reads: lines of `read_repeat_name.txt`;
  - flanking reads: FASTQ lines / 4;
  - non-reference calls per tier: lines of each `.txt` tier, minus a header if
    one exists;
  - reference calls: `all_ref_insert.txt` lines.
- Prefer one new process that takes the already-published outputs as inputs,
  so no existing module script changes. Use `wc -l`/`awk` in the image, or a
  `relocaTE3 cohort --sample-stats` mode if that is cleaner and testable in
  Python. **Never** publish these into `<outdir>/<sample>/`: the M1 parity
  file list must stay 33 files. Publish to `<outdir>/multiqc/` inputs or
  `pipeline_info/`.

Commit: "Per-sample stats for the summary report".

### Task 5: Wire COHORT and MULTIQC into the workflow
- **Modules:**
  - `modules/local/relocate3_cohort/main.nf` runs `relocaTE3 cohort` on the
    collected characterized tables. Use `collect()` and sort by sample for
    `-resume` determinism.
  - `modules/local/multiqc/main.nf` runs
    `multiqc --config <assets/multiqc_config.yml> --filename relocate3_report.html .`
    over cohort `_mqc` files, sample stats, `software_versions.yml` (as nf-core
    does, as custom content) and the Nextflow `pipeline_info` if useful.
  - Both get `versions.yml` and stubs.
- **`assets/multiqc_config.yml`:** title "RelocaTE3", report comment with the
  pipeline version and container, section order.
- **Workflow:**
  - COHORT runs once after every CHARACTERIZE, only when at least one sample
    was characterized.
  - MULTIQC runs once at the end.
  - New params `skip_multiqc` (false) and `multiqc_title`, plus schema and
    pytest.
- **Publishing:**
  - `<outdir>/cohort/` gets the 3 cohort files;
  - `<outdir>/multiqc/` gets `relocate3_report.html` and its data dir;
  - per-sample dirs unchanged, still 33 files: assert this in the existing
    nf-test.
- **nf-test:**
  - stub tests: cohort and multiqc run once; `--genotype false` with no
    genotype_bam means no COHORT; `--skip_multiqc` means no MULTIQC;
  - a real test on the Chr3 test profile: the cohort table equals the single
    sample's characterized rows, with 197 lines including the header.
- `tests/nextflow/parity_chr3.sh` must still PASS. Its per-sample comparison
  is unchanged.

Commit: "Workflow: cohort table and MultiQC report".

### Task 6: Real-rice cohort check
- Run `relocaTE3 cohort` on the 10 M2 real-rice characterized tables locally.
  It is small; run it in an sbatch job anyway to stay off the login node.
- Check rows against per-sample files exactly. Report:
  - the matrix shape;
  - the number of sites shared by ≥2 samples;
  - status counts per sample.
- Optionally re-run the real-rice pipeline from the M2 output dir with
  `-resume`, if the tag/branch logic allows, so only COHORT and MULTIQC run.
  Otherwise document the command for the maintainer after CI publishes a new
  image.

### Task 7: Docs and review
- README: cohort outputs, matrix semantics (exact sites, `.` = no call) and
  the MultiQC report.
- CHANGELOG.
- `docs/2026-11-05-nextflow-m3.md`: commands, results, decisions (the MultiQC
  environment choice, the status-class rule, the duplicate-site rule).
- Final whole-branch review as in M1/M2; fix Important findings.
- PR message at `../RelocaTE3/results/pr_messages/nextflow-m3.md`.
