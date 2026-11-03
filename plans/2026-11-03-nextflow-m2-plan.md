# Nextflow M2 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:**
- Run the pipeline from a published container: `-profile apptainer`,
  `singularity` or `docker`.
- Run it on UCR HPCC SLURM: `-profile hpcc,apptainer`.
- Per-sample output identical to `run-all` (Chr3) and to the v0.2.0rc2 results
  (10 real-rice samples).

**Architecture:**
- **Image:** one image built from `pixi.lock`'s `blat` environment, which holds
  the validated tool versions, plus a non-editable install of RelocaTE3. One
  recipe, `containers/Dockerfile`, built only by GitHub Actions.
- **Tags:** CI publishes to `ghcr.io/stajichlab/relocate3`:
  - `pr-<n>` on pull requests from this repo;
  - `dev` on `main`;
  - `<version>` on tags.
- **Cluster:** no Docker or Podman, so the cluster tests that image by pulling
  it with Apptainer 1.5.3.

Design: `plans/2026-11-02-nextflow-pipeline-design.md` (milestone M2).
M1 note: `docs/2026-11-03-nextflow-m1.md`.

**Branch / worktree:** `feat/nextflow-m2` in `../RelocaTE3-nextflow-m2`, from
`origin/main` 1b026ae, which merges PR #62.

**Two phases**, because the image only exists once CI has run:
- **Phase A, local (Tasks 0–5):** recipe, CI workflows, profiles, real-rice
  scripts. Everything checkable without the image.
- **Phase B, after the maintainer pushes (Tasks 6–8):**
  - the maintainer pushes the branch and opens a draft PR, so CI builds
    `pr-<n>`;
  - the maintainer makes the ghcr package public once;
  - then the parity gate and real rice run in the container.

---

## Conventions (same as M1)

**Test commands:**
- pytest: `E=$PWD/.pixi/envs/blat/bin; PATH=$E:$PATH PYTHONPATH=$PWD/src $E/python -m pytest tests -p no:cacheprovider --junitxml=<scratch>/x.xml`.
  Read the counts from the XML.
- `pixi run -e nf nf-test test`
- `pixi run -e nf nextflow lint .`
- **Baseline at 1b026ae:** pytest 393, nf-test 30, lint clean. Record it in Task 0.

**Process:**
- TDD: failing test first, then the fix, then the full suites.
- sbatch scripts: `#!/usr/bin/bash -l`, `set -euo pipefail`,
  `logs/<name>.%j.log`, partition `intel`, date and parameters printed at the
  start.
- Never on the login node. Wait for jobs with a
  `timeout … until sacct …` loop.
- `module load apptainer/1.5.3` only works in a login shell (`bash -lc`) or an
  sbatch script.
- Commit per task with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Never push.

---

## Phase A

### Task 0: Setup

- Run `pixi install -e blat -e nf` in the worktree.
- Run the three suites and record the counts.

### Task 1: Container recipe

**Files:** `containers/Dockerfile`, `.dockerignore`,
`containers/check_image.sh`.

**Requirements** (verify every instruction against the official docs:
pixi.sh deployment/container docs and the `ghcr.io/prefix-dev/pixi` image
tags):

1. **Build stage:**
   - `FROM ghcr.io/prefix-dev/pixi:<pinned version>`.
   - Copy `pixi.toml`, `pixi.lock`, `pyproject.toml`, `src/`, `README.md`, and
     `.git`, which versioningit needs for the version.
   - `pixi install --locked -e blat`.
   - Then install RelocaTE3 **non-editable** into that environment. The env's
     editable install would point at the build directory.
   - Write `pixi shell-hook -e blat` to a file.
2. **Final stage:**
   - A slim base matching the pixi image's glibc, typically
     `debian:<same>-slim` or `ubuntu`.
   - Copy the environment to the **same absolute path**: conda environments
     aren't relocatable.
   - Copy the shell hook.
   - Install `procps`: Nextflow's task metrics need `ps`.
   - Set `PATH` to the env's `bin`, so tools resolve without an entrypoint.
     Nextflow runs `bash .command.run` and ignores ENTRYPOINT under Apptainer.
3. **No `.pixi/`, `work/`, `results/` or `tests/data/` in the build context.**
   That is what `.dockerignore` is for.
4. **`check_image.sh`** runs *inside* the image and prints, then asserts:
   - `relocaTE3 --version` is a real version, not `0.0.0` or the
     `default-tag`;
   - `blat`, `bwa`, `samtools`, `minimap2`, `seqtk`, `bedtools` and `bowtie2`
     run, with versions equal to `pixi.lock`'s `blat` env;
   - `ps` exists;
   - `python -c "import RelocaTE3, pysam"` works.

**Dry run** (no Docker available): run the recipe's commands in an Apptainer
sandbox from `docker://ghcr.io/prefix-dev/pixi:<same tag>`. Use
`apptainer build --sandbox`, then `apptainer exec --writable --fakeroot`, or a
scratch-dir equivalent. Then run `check_image.sh` inside it. This checks the
commands, not the Dockerfile syntax; CI checks that. Run as an sbatch job:
`containers/dry_run_recipe.sh`, committed.

**Test:** a pytest checks `.dockerignore` excludes `.pixi`, `work`, `results`
and `tests/data`, and that the Dockerfile:
- uses `--locked`;
- pins the pixi image tag (no `latest`);
- copies the env to the same path.

Commit: "Container recipe from pixi.lock (blat env) + in-image checks".

### Task 2: GitHub Actions

**Files:** `.github/workflows/container.yml`, `.github/workflows/ci.yml`.

- **`container.yml`** runs on PRs, pushes to `main`, and `v*` tags:
  - Steps: `actions/checkout` with `fetch-depth: 0` for versioningit tags,
    `docker/setup-buildx-action`, `docker/login-action` to `ghcr.io` with
    `GITHUB_TOKEN`, `docker/metadata-action` (tags `pr-<n>`, `dev` on main,
    semver on tags, `sha-<short>`), `docker/build-push-action`.
  - Permissions `packages: write`, `contents: read`.
  - Push only when the PR comes from this repository, not a fork.
  - Run `check_image.sh` in the built image before pushing.
  - Image name `ghcr.io/stajichlab/relocate3`, lower case.
- **`ci.yml`** runs on PRs and on `main`:
  - Job 1: `prefix-dev/setup-pixi` with the `blat` and `nf` envs, then pytest,
    `nf-test test` and `nextflow lint .`.
  - Job 2: after the image job (`workflow_run`, or one workflow with `needs:`;
    pick the simpler, documented way),
    `nextflow run . -profile test,docker --container_tag pr-<n>`, then
    `tests/nextflow/parity_chr3.sh` in container mode (see Task 3).
- Pin every action to a released major version. Use only documented inputs;
  check each action's README.
- **Test:** a pytest parses the YAML (PyYAML is in the env? If not, use a
  minimal check). It asserts `fetch-depth: 0`, `packages: write` and
  `--locked` use, and that `ci.yml` runs pytest, nf-test and lint. Also
  `actionlint` if available through pixi/conda-forge; optional, so note it if
  absent.

Commit: "CI: build/publish image to ghcr; tests on every PR".

### Task 3: Container profiles and image parameter

**Files:** `nextflow.config`, `nextflow_schema.json`,
`tests/nextflow_schema_test.py`, `conf/base.config` or a new
`conf/containers.config`, `tests/nextflow/parity_chr3.sh`.

- **Params:**
  - `container_registry = 'ghcr.io/stajichlab/relocate3'`.
  - `container_tag`: default = the manifest version for releases. For dev
    builds the default is `'dev'`. Decide and justify: a `0.3.0dev` manifest
    should not point at a non-existent `:0.3.0dev` tag.
  - `process.container = "${params.container_registry}:${params.container_tag}"`,
    set only in the container profiles, so `pixi` ignores it.
- **Profiles:**
  - `docker`: `docker.enabled = true`,
    `docker.runOptions = '-u $(id -u):$(id -g)'`.
  - `apptainer`: `apptainer.enabled = true`, `apptainer.autoMounts = true`,
    `apptainer.cacheDir` from the `NXF_APPTAINER_CACHEDIR` env var if set.
  - `singularity`: the same settings with `singularity.*`.
  - Each profile disables the others (nf-core pattern).
  - The `pixi` profile must not set a container.
- **Parity script:** `parity_chr3.sh` gains a mode argument:
  - `pixi`, the default and current behaviour;
  - `apptainer` or `docker`. The `run-all` side then also runs inside the same
    image (`apptainer exec <image> relocaTE3 run-all …`), so both sides use
    identical tools.
  - Keep the existing outputs and the PASS/FAIL contract.
- **Test:** nf-test and pytest check:
  - `nextflow config -profile test,apptainer` sets `process.container`, and
    `apptainer.enabled` is true;
  - `-profile test,pixi` has no container;
  - the schema has the new params.

Commit: "Container profiles (apptainer/singularity/docker), image params".

### Task 4: `hpcc` profile (UCR SLURM)

**File:** `conf/hpcc.config`, included from the `hpcc` profile.

- `process.executor = 'slurm'`; `process.queue = params.slurm_partition`.
  New param, default `'intel'`, in the schema.
- `executor { queueSize = 50; submitRateLimit = '10/1min' }`. Justify the
  values and keep them conservative.
- `process.resourceLimits = [cpus: 64, memory: 440.GB, time: 30.d]`: intel
  nodes have 64 CPUs and about 450 GB (`sinfo -p intel`). Check other
  partitions only if the param allows them; the limits apply to `intel`.
- `apptainer.cacheDir`: the `NXF_APPTAINER_CACHEDIR` env var if set, else
  `${HOME}/.apptainer/nextflow-cache`, with the reason documented. The pull
  happens on the head job's node.
- Document the launch pattern: a small sbatch head job runs
  `module load nextflow apptainer` or `pixi run -e nf nextflow`, and Nextflow
  submits tasks. Never launch on the login node.
- **Test:** `nextflow config -profile hpcc,apptainer` shows executor slurm,
  queue intel and the limits; `--slurm_partition epyc` changes the queue.

Commit: "hpcc profile: SLURM on UCR HPCC".

### Task 5: Real-rice pipeline scripts

**Files:**
- `validation/real_rice/make_samplesheet.py`
- `validation/real_rice/run_pipeline.sh`
- `validation/real_rice/compare_pipeline_to_rc2.sh`

1. **`make_samplesheet.py`** reads `validation/real_rice/config.toml`, falling
   back to `config.example.toml`, and writes `sample,fastq_1,fastq_2,genotype_bam`:
   - absolute paths: `fq_dir/<s>_trim_R1.fastq.gz`, `_R2`;
   - genotype BAM: `genome_aln_dir/<s>/<s>.sorted.markdup.cram`, whose `.crai`
     exists.
   - Test it with pytest on a temporary fake tree.
2. **`run_pipeline.sh`:** an sbatch head job, about 2 CPUs, 8 GB and 3 days.
   - Settings equal to `run_relocate3.sh`'s run-all:
     `--te_library mping.fa --te_name mping --size 250 --mismatch 2 --min_mapq 0 --tsd UNK --target ALL --repeatmasker MSU_r7.fa.out --genome MSU_r7.fa`.
   - Read them from the same config, so they can't drift.
   - `-profile hpcc,apptainer --container_tag <arg>` and `-resume`.
   - `-work-dir` on bigdata, and `--outdir results/nextflow-real-rice/<tag>`.
3. **`compare_pipeline_to_rc2.sh`:** for each sample, run
   `validation/blat_memory/compare_outputs.sh` between
   `../RelocaTE3/results/release-validation/v0.2.0rc2/clone/validation/real_rice/results/<s>`
   and the pipeline's `<outdir>/<s>`. Also compare the file lists. Print one
   summary line per sample and an overall PASS/FAIL.
   - Check what the rc2 directory contains first: run-all with
     `--genotype-bam` writes no `reads.genome.bam`. Explain any expected
     file-set difference in the script header, and never silently exclude a
     differing file.
   - Note: the genome there may already have index files from earlier runs.
     The pipeline builds its own in the work directory, which is fine.

Commit: "Real-rice via the pipeline: sample sheet, sbatch head job, rc2 comparison".

**End of Phase A.**
- Full suites green.
- Update `docs/2026-11-03-nextflow-m2.md` with progress so far.
- **Stop and hand off to the maintainer.** Push the branch, open a **draft**
  PR, wait for `container.yml` to publish `pr-<n>`, then make the ghcr package
  public: GitHub → Packages → relocate3 → Package settings → Change
  visibility. This is the maintainer's action.

---

## Phase B (after the image exists)

### Task 6: Image check and parity in the container
- `apptainer pull` `docker://ghcr.io/stajichlab/relocate3:pr-<n>` in an sbatch
  job, then run `check_image.sh` inside it.
- `sbatch tests/nextflow/parity_chr3.sh apptainer pr-<n>`: must PASS.

### Task 7: Real rice through the pipeline
- `sbatch validation/real_rice/run_pipeline.sh pr-<n>`.
- When it finishes, run `compare_pipeline_to_rc2.sh`. Every sample must be
  identical to rc2, or each difference explained.
  - Thread counts differ: rc2 used 16. bwa aln/sampe should not depend on
    threads. If a difference appears, investigate before concluding anything.
- Record wall time and peak memory per task from the Nextflow trace.

### Task 8: Docs and hand-off
- **README:** container profiles, the `hpcc` launch pattern, and
  `container_tag`.
- **CHANGELOG** and AGENTS.md.
- **`docs/2026-11-03-nextflow-m2.md`:** commands, job IDs, results,
  timings, failures, next steps (M3: cohort table + MultiQC).
- Final whole-branch review, as in M1. Fix Important findings.
- PR message file at `../RelocaTE3/results/pr_messages/nextflow-m2.md`.
