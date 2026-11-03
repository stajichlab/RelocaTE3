# Nextflow M1 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (or subagent-driven-development) to implement this plan task-by-task.

**Goal:** Run one or more samples through RelocaTE3 with
`nextflow run . -profile test,pixi`, producing per-sample files identical to
`relocaTE3 run-all` on the Chr3 2 Mb fixture.

**Architecture:**
- Each Nextflow process runs exactly one existing `relocaTE3` subcommand, with
  the arguments `cmd_run_all` passes.
- Two CLI additions make every `run-all` step reachable as a subcommand:
  - `index-genome --genome-aligner`
  - `align-reads`
- `run-all` uses both, so the two paths can't drift apart.
- nf-core layout; nf-schema validation.

Design: `plans/2026-11-02-nextflow-pipeline-design.md`.

**Tech stack:** Python 3.12 (argparse CLI, pytest), Nextflow ≥ 25.04 DSL2,
nf-schema plugin, nf-test, pixi.

**Branch / worktree:** `feat/nextflow-m1` in `../RelocaTE3-nextflow`, from
`origin/main` be68347.

---

## Conventions for every task

**Test command:**
```bash
E=<repo>/.pixi/envs/blat/bin   # after `pixi install -e blat` in this worktree
PATH=$E:$PATH PYTHONPATH=$PWD/src $E/python -m pytest <tests> -p no:cacheprovider --junitxml=<scratch>/x.xml
```
- Read the counts from the XML: `addopts -ra -q` hides the summary.
- Without the env's `bin` on PATH, 8 tests fail with "minimap2 not found".

**Baseline:** the full suite at be68347 is 347 passed. Record it before Task 1.

**Process:**
- TDD: write the test, watch it fail, implement, watch it pass, run the full suite.
- Never run heavy jobs on the login node. This session is a SLURM job, so small
  runs are fine. Long runs use sbatch: `#!/usr/bin/bash -l`,
  `set -euo pipefail`, logs in `logs/<name>.%j.log`, partition `intel`.
- Commit after each task. End messages with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Never push.

---

### Task 0: Environments and the run-all baseline

**Files:** `pixi.toml` (new `nextflow` feature/env); `pixi.lock`.

**Step 1: Install the `blat` env in this worktree.**
- Run `pixi install -e blat`.
- Run the full suite and record the count. Expect 347 passed.

**Step 2: Add a standalone Nextflow env.** It must not co-solve with blat:
bioconda blat pins zlib 1.2.11.

```toml
[feature.nextflow]
platforms = ["linux-64", "osx-64", "osx-arm64"]
[feature.nextflow.dependencies]
nextflow = ">=25.04"
nf-test = ">=0.9"
[environments]
nf = { features = ["nextflow"], no-default-feature = true }
```

- Run `pixi install -e nf`.
- Check: `pixi run -e nf nextflow -version` reports ≥ 25.04, and
  `pixi run -e nf nf-test version` works.

**Step 3: Make the run-all reference output for the Chr3 fixture.** Copy the
genome to scratch first, so indexes are not written into `tests/data/`.

```bash
B=results/nextflow-m1/baseline_run_all; mkdir -p $B/ref
cp tests/data/sim_genome/MSU7.Chr3_2M.fa $B/ref/
PATH=$E:$PATH PYTHONPATH=$PWD/src $E/python -m RelocaTE3 run-all \
  -1 tests/data/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz \
  -2 tests/data/sim_reads/MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz \
  -T tests/data/TE_lib/RiceTE.fa -g $B/ref/MSU7.Chr3_2M.fa \
  -n chr3 -o $B/chr3 --threads 4 --mismatch 2 --tsd UNK --te-name repeat \
  --require-both-junctions \
  --repeatmasker tests/data/sim_genome/MSU7.Chr3_2M.fa.RepeatMasker.out --genotype
(cd $B/chr3 && find . -type f | sort) > $B/file_list.txt
```

- These are the flags of `tests/acceptance_test.py`
  `test_run_all_cli_end_to_end`. Note that `find-reference` also writes
  `existingTE.bed` at the output root.
- `results/` is git-ignored; confirm it with `git check-ignore`.
- `file_list.txt` is the contract for Task 6: which files each step leaves,
  under which subfolder.

**Step 4: Commit** `pixi.toml` and `pixi.lock`:
"pixi: standalone nextflow env (nextflow, nf-test)".

---

### Task 1: `index-genome --genome-aligner`

**Files:**
- Modify: `src/RelocaTE3/cli.py`
  - `_menu_index_genome` (~:635)
  - `cmd_index_genome` (~:1265)
  - `cmd_run_all` step 1 (~:966-981)
- Test: `tests/index_genome_aligner_test.py` (new)

**Why:**
- `align-genome` builds the bwa index next to the genome on demand
  (`backend.index(...)`, cli.py ~:1312/1362).
- In Nextflow, parallel samples must not each build it. One task must build
  it once.

**Step 1: Failing tests.**
```python
import random
from pathlib import Path
import pytest
from RelocaTE3 import cli
from RelocaTE3.cli import build_parser, main
from RelocaTE3.aligners import GENOME_ALIGNERS

def _genome(tmp_path) -> Path:
    rng = random.Random(1)
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">c{i}\n{''.join(rng.choice('ACGT') for _ in range(3000))}\n" for i in range(2)))
    return fa

def test_index_genome_default_aligner_is_bwaaln():
    ns = build_parser().parse_args(["index-genome", "-g", "g.fa"])
    assert ns.genome_aligner == "bwaaln"

def test_index_genome_rejects_blat():
    with pytest.raises(SystemExit):
        build_parser().parse_args(["index-genome", "-g", "g.fa", "--genome-aligner", "blat"])

def test_index_genome_builds_the_genome_aligner_index(tmp_path):
    fa = _genome(tmp_path)
    assert main(["index-genome", "-g", str(fa), "--genome-aligner", "bwaaln"]) == 0
    for ext in (".fai", ".mmi", ".bwt"):
        assert Path(f"{fa}{ext}").is_file(), ext

def test_run_all_indexes_for_its_genome_aligner(monkeypatch, tmp_path):
    """Even when .fai/.mmi exist, run-all must make sure the genome-aligner index exists."""
    # reuse tests/main_test.py's rec()/_fake_stage_outputs pattern; touch ref.fa.fai and ref.fa.mmi,
    # run run-all with --genome-aligner bwaaln, assert the first call is index-genome with
    # genome_aligner == "bwaaln"
```

Run the tests. Expected:
- The first three fail: no `genome_aligner` attribute, and no `.bwt`.
- The run-all test fails: today the step is skipped when `.fai` and `.mmi` exist.

**Step 2: Implement.**
- `_menu_index_genome`: add `--genome-aligner` (`dest="genome_aligner"`,
  `default="bwaaln"`, `choices=list(GENOME_ALIGNERS)`). Copy the help wording
  from `_menu_align_genome`.
- `cmd_index_genome`: after `aln.index_genome(...)`, add
  `get_aligner(args.genome_aligner).index(args.genome_fasta, force=args.force)`.
- `cmd_run_all`: always call `_step("index-genome", ...)` with
  `genome_aligner=args.genome_aligner`, and remove the `.fai`/`.mmi` skip.
  Every index step already skips existing files: faidx is checked,
  `index_minimap` has `force=False`, and the backend's `index` checks its
  `_index_sentinel`.
- In `test_run_all_delegates_to_the_staged_commands`, also check
  `by["index-genome"].genome_aligner == "bwaaln"`.

**Step 3:** new tests pass; full suite green.

**Step 4: Commit:**
"index-genome builds the genome aligner's index; run-all always delegates".

---

### Task 2: `align-reads` subcommand (genotyping alignment)

**Files:**
- Modify: `src/RelocaTE3/cli.py`
  - new `_menu_align_reads` and `cmd_align_reads`
  - register it next to `align-genome` (~:1440-1560)
  - `cmd_run_all` step 7 (~:1072-1090)
- Test: `tests/align_reads_test.py` (new)

**Contract:**
`relocaTE3 align-reads -1 R1 [-2 R2] -g GENOME -n NAME -o OUTDIR [--threads N]`
writes `OUTDIR/genome_aln/NAME.reads.genome.bam` and its `.bai`. This is the
same path and the same code (`Aligner(threads).map_library_to_genome(...)`)
`run-all` uses today.

**Step 1: Failing tests.**
- Parser test: required arguments, `--threads` default 1.
- Functional test:
  - simulate 20 read pairs from the Task 1 random genome: 100 bp reads,
    300 bp insert, reverse-complement R2;
  - run `main(["align-reads", ...])`;
  - assert the BAM exists, is coordinate-sorted (`@HD SO:coordinate`), is
    indexed, and has 40 mapped records.
- Run-all delegation test, copying the `rec()` pattern:
  - `--genotype` without `--genotype-bam` calls `"align-reads"` before
    `"characterize"`;
  - `characterize.bam == [str(outdir/"genome_aln"/"HEG4.reads.genome.bam")]`;
  - with `--genotype-bam x.cram`, `align-reads` is not called.

**Step 2: Implement.**
- Move the inline block from `cmd_run_all` into `cmd_align_reads` unchanged.
- `cmd_run_all` calls `_step("align-reads", cmd_align_reads, ...)`.

**Step 3:** tests pass. Full suite green, including `acceptance_test.py`
`test_run_all_cli_end_to_end`.

**Step 4: Re-run the Task 0 Step 3 command** into
`results/nextflow-m1/run_all_after_cli/`. Compare it with the baseline:
```bash
bash validation/blat_memory/compare_outputs.sh results/nextflow-m1/baseline_run_all/chr3 results/nextflow-m1/run_all_after_cli/chr3
```
Expected: `differing=0`. This proves Tasks 1–2 didn't change run-all output.

**Step 5: Commit:**
"align-reads subcommand for the genotyping alignment; run-all delegates".

---

### Task 3: Pipeline skeleton and parameter validation

**Files:**
- `main.nf`
- `nextflow.config`
- `nextflow_schema.json`
- `assets/schema_input.json`
- `conf/base.config`
- `conf/test.config`
- `nf-test.config`
- `tests/main.nf.test`
- `.gitignore` additions: `.nextflow*`, `work/`, `.nf-test/`

**Step 1: Failing nf-test** (`tests/main.nf.test`), pipeline-level:
- `-profile test,pixi --tsd TTA`: fails, and stdout/stderr names `tsd` and `UNK`.
- `--input missing.csv`: fails during validation, before any process runs.
- `-profile test,pixi -stub`: succeeds. Stubs come in Task 4–5; until then this
  test is expected to fail.

Run `pixi run -e nf nf-test test tests/main.nf.test`. Expected: FAIL (no main.nf).

**Step 2: Implement.**

`nextflow.config`:
- `manifest`:
  - `name = 'stajichlab/RelocaTE3'`
  - `nextflowVersion = '!>=25.04.0'`
  - `mainScript = 'main.nf'`
  - `version` read from `git describe` is not possible in config, so set it to
    the package version string and add a pytest check that it matches
    `RelocaTE3.__version__`'s release part.
- `plugins { id 'nf-schema@2.x.y' }`: pin the newest 2.x; check with
  `nextflow plugin list` or the registry.
- `params`, with the names and defaults from the design doc:
  - `input`, `genome`, `te_library`, `repeatmasker = null`, `outdir`
  - `te_aligner = 'blat'`, `genome_aligner = 'bwaaln'`
  - `mismatch = 2`, `min_match = 10`, `min_trimmed = 10`, `size = 500`,
    `min_mapq = 0`
  - `tsd = 'UNK'`, `require_both_junctions = true`, `target = 'ALL'`
  - `te_name = 'repeat'`
  - `genotype = true`
  - `pixi_env = "${projectDir}/.pixi/envs/blat"`
- Profiles:
  - `pixi`: `process.beforeScript = "export PATH=${params.pixi_env}/bin:\$PATH"`
  - `test`: includes `conf/test.config`, which points at `tests/data/...`,
    with `te_name` and the TE library matching Task 0, a sample sheet in
    `tests/data/samplesheet_chr3.csv`, and resource caps of 4 CPUs and 8 GB.
- `includeConfig 'conf/base.config'`: labels `process_single`, `process_low`,
  `process_medium` and `process_high` (TE search: 8 CPUs, 8 GB, 24 h);
  `errorStrategy` retries on exit 137/140 with memory times the attempt number.

`nextflow_schema.json` (draft 2020-12, nf-schema 2):
- `$defs` groups `input_output_options`, `calling_options`,
  `genotyping_options` and `profile_options`, each with `description` and
  `help_text`.
- `tsd: {"type":"string","enum":["UNK","UKN","unknown"],"default":"UNK"}`
- `te_aligner` and `genome_aligner` enums copied from
  `aligners.TE_ALIGNERS`/`GENOME_ALIGNERS`, minus `blat` for the genome.
- Add a pytest test that the schema enums equal those Python tuples, so they
  can't drift.

`assets/schema_input.json`: fields `sample` (string, no spaces), `fastq_1`
(file, exists, `.f(ast)?q(.gz)?`), `fastq_2` (optional, same) and
`genotype_bam` (optional, `.bam` or `.cram`, exists).

`main.nf`:
- `validateParameters()`
- `samplesheetToList(params.input, "assets/schema_input.json")`
- call `RELOCATE3`, which starts as an empty workflow until Task 6.

**Step 3:** the two validation nf-tests pass; `pixi run -e nf nextflow lint .`
is clean.

**Step 4: Commit:** "Nextflow skeleton: config, profiles, schema validation".

---

### Task 4: Step processes (modules)

**One folder each,** `modules/local/relocate3_<step>/main.nf` plus
`tests/main.nf.test`:

| Process | Command | Main inputs | Outputs (fill exact globs from Task 0 `file_list.txt`) |
|---|---|---|---|
| `RELOCATE3_INDEX_GENOME` | `relocaTE3 index-genome -g $genome --genome-aligner $a` | genome | genome + all index files (`*.fai`, `*.mmi`, aligner index) |
| `CAT_FASTQ` | `cat` (gzip-safe concatenation) | lanes of one sample | `${meta.id}_1.fq.gz`, `_2` |
| `RELOCATE3_TE_SEARCH` | `relocaTE3 run -1 -2 -T -n -o . --threads --te-aligner --min-match --min-trimmed --mismatch` | reads, TE library | `te_containing/`, `flanking/`, `te_portions/`, left/right BAMs |
| `RELOCATE3_ALIGN_GENOME` | `relocaTE3 align-genome -g -f flanking/*.flankingReads.fq -n -o . -1 -2 --genome-aligner --threads` | genome + index, reads, `te_containing/`, `flanking/` | `*.repeat.<aligner>.sorted.bam(.bai)`, `genome_aln/` |
| `RELOCATE3_FIND_INSERTIONS` | `relocaTE3 find-insertions` with every argument `cmd_run_all` passes | genome BAM, `read_repeat_name.txt`, fullreads BAM, RepeatMasker (optional) | `results/` non-reference files |
| `RELOCATE3_FIND_REFERENCE` | `relocaTE3 find-reference` | genome BAM, read_repeat, RepeatMasker | `results/` reference files |
| `RELOCATE3_ALIGN_READS` | `relocaTE3 align-reads` | reads, genome + index | `genome_aln/*.reads.genome.bam(.bai)` |
| `RELOCATE3_CHARACTERIZE` | `relocaTE3 characterize -s results/<target>.<te_name>.all_nonref_insert.txt -b <bam> -g genome -o results` | sites, BAM or CRAM, genome | `results/*characTErized*` |

**Rules:**
- Inputs are staged under the folder names the subcommand expects, e.g.
  `path(te_containing, stageAs: 'te_containing/*')`. Commands run with `-o .`.
- The argument lists are copied from `cmd_run_all`, not retyped from memory.
  Add a pytest that greps each module's `script:` block for each flag
  `cmd_run_all` passes to that step, so a new run-all argument can't be
  silently dropped from the pipeline.
- Every process:
  - writes `versions.yml` (`relocaTE3 --version`, plus the main tool: `blat`,
    `bwa`, `samtools`);
  - has a `stub:` block that touches the declared outputs;
  - has a `tag "$meta.id"`;
  - takes resource labels from `conf/base.config`.
- `publishDir` stays out of the modules; it goes in `conf/modules.config`
  (Task 6).

**Step 1:** a failing nf-test stub test per module, checking that outputs are
emitted with the expected names.

**Step 2:** implement until the stub tests pass.

**Step 3:** one real nf-test, for `RELOCATE3_INDEX_GENOME` on the random
genome from Task 1: builds `.bwt`.

**Step 4: Commit** after every two or three modules.

---

### Task 5: Subworkflows and workflow

**Files:**
- `subworkflows/local/prepare_genome.nf`: one `RELOCATE3_INDEX_GENOME`; emits
  `genome_with_index`, a single value channel of the genome plus its index
  files.
- `subworkflows/local/one_sample.nf`: chains the per-sample steps:
  - `ALIGN_READS` runs from the reads alone, in parallel with `TE_SEARCH`,
    only when `params.genotype` is true and `meta.genotype_bam` is empty;
  - `CHARACTERIZE` takes `genotype_bam` if given, else `ALIGN_READS`' BAM;
  - `FIND_REFERENCE` runs only with `params.repeatmasker`.
- `workflows/relocate3.nf`:
  - group sample-sheet rows by `sample`;
  - more than one row → `CAT_FASTQ`;
  - `PREPARE_GENOME` once;
  - `ONE_SAMPLE` per sample;
  - collect `versions.yml` into `pipeline_info/software_versions.yml`.
- `conf/modules.config`: `publishDir` for each process, to
  `${params.outdir}/${meta.id}/<same relative path as run-all>`, mode `copy`.
  Publish exactly the files in Task 0's `file_list.txt` and nothing else; no
  staged inputs.
- Nextflow `report`, `timeline` and `trace` go to
  `${params.outdir}/pipeline_info/`.

**Step 1:** the stub pipeline test from Task 3 (`-profile test,pixi -stub`)
plus a new stub test:
- a sample sheet with two rows for one sample;
- `CAT_FASTQ` runs once;
- every step runs once per sample.

Expected: fail.

**Step 2:** implement until both pass, then run `nextflow lint .`.

**Step 3: Commit:** "Wire per-sample workflow, shared genome index, lane merging".

---

### Task 6: Parity gate: pipeline vs run-all on Chr3

**Files:** `tests/nextflow/parity_chr3.sh` (sbatch-able), and
`tests/data/samplesheet_chr3.csv` if Task 3 didn't create it.

**Step 1: Write the script.**
- Runs `pixi run -e nf nextflow run . -profile test,pixi --outdir $OUT/nf`.
- Runs `relocaTE3 run-all` with the same settings into `$OUT/run_all/chr3`:
  the genome copied to `$OUT/ref`, and `--genotype`.
- Compares with
  `validation/blat_memory/compare_outputs.sh $OUT/run_all/chr3 $OUT/nf/chr3`.
- Exits non-zero on any difference.
- `#!/usr/bin/bash -l`, `set -euo pipefail`, `#SBATCH -p intel --cpus-per-task=4
  --mem=16G --time=02:00:00 -o logs/nextflow_parity.%j.log`; prints the date
  and all parameters at the start.

**Step 2: Run it.**
- `sbatch tests/nextflow/parity_chr3.sh`.
- Wait with a `timeout … until sacct …` loop; foreground `sleep` is blocked.

**Step 3: Required result:** `identical=N differing=0`, with
N = the number of files in `file_list.txt`, minus indexes.
- On any `DIFF`, find the cause: a missing flag, a path written into a text
  file, or a different input order.
- Fix the module, not the comparison.
- Never exclude a differing file without a documented, reviewed reason.

**Step 4: Commit:** "Parity gate: Nextflow equals run-all on Chr3 (N files identical)".

---

### Task 7: Docs

**Files:**
- `README.md`: new "Running with Nextflow" section. M1 supports only the
  `pixi` profile, so keep it short:
  - install with `pixi install -e blat -e nf`;
  - the sample-sheet format;
  - one command;
  - where outputs land;
  - "container profiles: coming in M2".
- `CHANGELOG.md`: an "Unreleased" section listing `align-reads`,
  `index-genome --genome-aligner` and the Nextflow pipeline (M1).
- `docs/2026-11-02-nextflow-m1.md`: dated note covering what was done, the
  commands, the parity result with job ID, failures, and next steps (M2).
- `AGENTS.md`: "Roadmap context" order: Nextflow → long reads → multiple
  genomes → pangenome → Rust.

**Step 1:** full pytest suite, `nf-test test` (all) and `nextflow lint .` are
green. Record the counts in the dated note.

**Step 2: Commit:** "Docs: Nextflow M1 usage, changelog, roadmap order".

**Step 3: Hand off.**
- Report test counts, the parity job ID and its result.
- Give a ready-to-paste pull-request message file under
  `../RelocaTE3/results/pr_messages/`.
- The maintainer pushes.
