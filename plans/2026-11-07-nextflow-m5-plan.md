# Nextflow M5 Implementation Plan: user docs, real-rice harness, v0.3.0rc1

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** the last Nextflow milestone:
- a user guide for the pipeline;
- the real-rice validation running through the pipeline instead of the
  `run-all` SLURM array;
- the full release check;
- the release candidate **`v0.3.0rc1`**. The maintainer chose to stay on
  release candidates (2026-10-29 and 2026-11-07).

Design: `plans/2026-11-02-nextflow-pipeline-design.md` (M5: usage docs;
real-rice harness on the pipeline; release `v0.3.0`; full release check passes).

**Branch / worktree:** `feat/nextflow-m5` in `../RelocaTE3-nextflow-m5`
(`--no-track`), from `origin/main` 2c9a448, which merges PR #66.
Push with `git push -u origin feat/nextflow-m5`.

**Release mechanics already in place** (M2 review fixes):
- `.github/workflows/container.yml` on a `v*` tag runs
  `containers/check_release.sh`. The build fails unless all of these equal
  `${GITHUB_REF_NAME#v}`:
  - the tag;
  - the image's `relocaTE3 --version`;
  - `manifest.version`;
  - the default `params.container_tag`.
- It then publishes `ghcr.io/stajichlab/relocate3:<PEP 440 version>`.
- `tests/nextflow_schema_test.py` requires a release or rc manifest version to
  default `container_tag` to that same version.

---

## Conventions (as M1–M4)

**Tests:**
- pytest in the `pipeline` env:
  `E=$PWD/.pixi/envs/pipeline/bin; PATH=$E:$PATH PYTHONPATH=$PWD/src $E/python -m pytest tests -p no:cacheprovider --junitxml=<scratch>/x.xml`.
  Read the counts from the XML.
- `pixi run -e nf nf-test test`
- `pixi run -e nf nextflow lint .`
- actionlint
- **Baseline at 2c9a448:** pytest 724, nf-test 51. Record it in Task 0.

**Process:**
- TDD for code changes.
- sbatch for heavy work (`#!/usr/bin/bash -l`, `set -euo pipefail`, `-p intel`,
  `logs/<name>.%j.log`).
- Never on the login node.
- Commits end with:
  ```
  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01CT7G8cnGmxRpYgacJ9uvqg
  ```
- Never push. The maintainer pushes and tags.

**Editable env caveat:** long runs that load this worktree's `src` must not
overlap edits to `src/`. Use a frozen `git archive` copy with
`PYTHONPATH=<copy>/src`, or the container.

---

### Task 0: Setup
- `pixi install -e blat -e pipeline -e nf`.
- Run the baseline suites and record the counts.

### Task 1: User guide

The README has grown to about 640 lines, mixing the CLI, RelocaTE2 notes and
the pipeline.

1. **New Sphinx page `docs/source/pipeline.rst`**, with an entry in the
   `index.rst` toctree:
   - Quick start: `-profile test,apptainer` from GitHub with `-r v0.3.0rc1`.
   - Requirements: Nextflow ≥ 26.04 and one of apptainer, singularity, docker
     or pixi.
   - Sample sheet: columns, lanes, `genotype_bam` with its index, name rules
     and reserved names.
   - Parameters by group. Generate the table from `nextflow_schema.json` with
     a small script or a Sphinx include, so it can't drift, plus a pytest that
     the page lists every schema param.
   - Profiles: `apptainer`, `singularity`, `docker`, `pixi`, `hpcc`, `test`,
     and profile order (the last one wins).
   - UCR HPCC launch: a head-job sbatch template,
     `NXF_APPTAINER_CACHEDIR`, `--slurm_partition` and `-resume`.
   - Outputs: per-sample layout (same as `run-all`), `cohort/`, `multiqc/`,
     `pipeline_info/` (`software_environment.txt`).
   - BLAT scatter: what it does, `--blat_chunks_per_task`, disk cost, and when
     to disable it.
   - Reproducibility: pin `-r` and `--container_tag`, `software_versions.yml`,
     mutable tags.
   - Troubleshooting: the common errors seen during M1–M4 and their fixes.
   - A short "for agents" section: `--help`, the schema, `-params-file`, and
     where provenance is recorded.
2. **Shorten the README's Nextflow section** to a quick start plus a link to
   the guide. Keep the README accurate; don't remove CLI documentation.
3. **Build check:** `pixi run docs` (sphinx). If it fails because of
   environment issues, fix them minimally, or document the failure and the
   command. Add a pytest that the guide's commands use params and profiles
   that exist: parse the `--param` names from the code blocks and check them
   against the schema.

Commit: "Pipeline user guide (docs/source/pipeline.rst); README quick start".

### Task 2: Real-rice harness on the pipeline

- `validation/real_rice/run_validation.sh` currently runs `run_relocate3.sh`, a
  run-all SLURM array, then normalize and compare against the frozen RelocaTE2
  calls.
- **Add a pipeline mode, and make it the default:** run `run_pipeline.sh`
  (`-profile hpcc,apptainer`, with a pinned `--container_tag`), then the same
  normalize and compare steps on the pipeline outdir, which has the same
  per-sample layout.
- Keep the run-all mode as an option, `--mode run-all`, for comparing the two.
- Update `validation/real_rice/README.md` and
  `validation/release/submit_real_rice.sh`, which currently submits the run-all
  array. The release check should use the pipeline path, plus documentation of
  how the rc's image tag is supplied.
- **Tests:** in the style of `tests/real_rice_pipeline_scripts_test.py`, with
  `DRY_RUN` and fake trees. Check the mode switch and that both modes feed the
  same normalize and compare commands.
- **Real check:** run the pipeline mode on the 10 samples with the current
  `dev` image. Use `ALLOW_MUTABLE_TAG=1`, or better, `sha-2c9a448` if CI
  published it; check with `containers/pull_and_check.sh sha-2c9a448`. The
  RelocaTE2 comparison report must equal the rc2 report:
  - characterized: 5,763 matched, 5 RelocaTE2-only, 17 RelocaTE3-only;
  - non-reference GFF: 5,760 shared, 875 RelocaTE2-only, 15 RelocaTE3-only.
  - Diff the report TSVs ignoring path columns, as in
    `docs/2026-10-21-release-v0.2.0-validation.md`.

Commit: "Real-rice validation runs through the pipeline by default".

### Task 3: Release preparation for v0.3.0rc1
- **Version strings:**
  - `nextflow.config` `manifest.version = '0.3.0rc1'`;
  - default `params.container_tag = '0.3.0rc1'`, with the schema default
    updated;
  - `CITATION.cff` version and date (check its current fields);
  - CHANGELOG: turn the "Unreleased" section into `## 0.3.0rc1 — 2026-11-xx`.
    Add a short summary at the top: the Nextflow pipeline, containers, cohort
    and report, BLAT scatter, and the CLI boolean fix. Mention that 0.2.0
    stayed at rc2.
- **What `main` carries after the tag:**
  - The current "dev" convention (`0.3.0dev` with tag `dev`) would sort
    *before* 0.3.0rc1 in PEP 440 (`0.3.0.dev0 < 0.3.0rc1`).
  - Decide the post-release convention, e.g. `0.3.0rc2.dev0` with
    `container_tag 'dev'`. Check that the schema pytest rule ("dev suffix →
    default `dev`") still holds, and that versioningit's dev versions for
    commits after the tag look sensible.
  - Document it in `docs/` and in `AGENTS.md`, under "release steps".
  - Do **not** make the post-release commit in this branch. Write the exact
    steps for the maintainer, or provide a script, `scripts/post_release.sh`
    or similar, with a test.
- **`check_release.sh` dry run:** run it locally against this branch with a
  simulated `GITHUB_REF_NAME=v0.3.0rc1` and a locally built version string, if
  possible. Otherwise, at least run its config checks. Paste the result.
- **Release-steps doc**, for the maintainer, in
  `docs/2026-11-07-release-v0.3.0rc1.md`:
  1. merge this PR;
  2. on `main`, `git tag -a v0.3.0rc1 -m ...` and `git push origin v0.3.0rc1`;
  3. wait for `container.yml` on the tag: the release check, then the push
     of `:0.3.0rc1`;
  4. create the GitHub Release (`gh release create v0.3.0rc1 --prerelease --notes-file ...`);
  5. run the post-release version bump.

Commit: "Release preparation: v0.3.0rc1 version strings, CHANGELOG, steps".

### Task 4: Full release check (on the release commit, before tagging)
- `sbatch validation/release/fresh_clone_check.slurm <release commit sha> 0.3.0rc1`.
  The expected version will be the tag's. Before tagging, the clone reports a
  post-tag dev version, so adapt or document how this check is run pre-tag
  versus post-tag. Prefer to run it pre-tag on the commit, accepting the dev
  version, and post-tag on `v0.3.0rc1`.
- Chr3 parity: pixi mode, plus the forced-scatter variant.
- Real rice through the pipeline: from Task 2.
- Record all job IDs and results in the release doc.

### Task 5: Review and hand-off
- Final whole-branch review; fix Important findings.
- PR message at `../RelocaTE3/results/pr_messages/nextflow-m5.md`.
- Release notes file for `gh release create`, at
  `../RelocaTE3/results/pr_messages/v0.3.0rc1-release-notes.md`.

### After the maintainer tags (post-merge; this session)
- `containers/pull_and_check.sh --force-pull 0.3.0rc1`: PASS.
- `nextflow run stajichlab/RelocaTE3 -r v0.3.0rc1 -profile test,apptainer`
  from a scratch dir, with the module's Nextflow: PASS.
- Record both in the release doc in a follow-up PR, together with the
  post-release version bump.
