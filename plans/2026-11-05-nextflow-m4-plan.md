# Nextflow M4 Implementation Plan: BLAT scatter across jobs

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** spread the BLAT TE search of one sample across many SLURM jobs, so a
30x multi-family sample takes about an hour instead of about 13 hours.

The outputs must be **identical** to the unsplit `relocaTE3 run`:
- the TE BAMs, record by record;
- `te_containing/`, `flanking/` and `te_portions/`;
- everything downstream.

Design: `plans/2026-11-02-nextflow-pipeline-design.md` (milestone M4: "new
chunk-search and merge-trim subcommands; identical to unsplit; 30x wall time
measured").

**Why this is safe to build:** the BLAT path is already chunked inside one
process (`src/RelocaTE3/aligners.py` `BlatBackend`):
- `_write_query_chunks`: seqtk → FASTA → `split` into `BLAT_CHUNK_SEQS` = 200,000-read chunks;
- `_blat_chunk`: one chunk → sorted, indexed chunk BAM with canonical mate names → `TEChunk(bam, n_reads)`;
- `_blat_side`: a thread pool over the chunks, then `_merge_sorted_bams`;
- trimming in `src/RelocaTE3/librelocate.py`:
  `identify_TE_reads` → `map_te_library_chunked` → `write_trimmed_reads_chunked`,
  which consumes the chunk BAMs in order and needs `n_reads` per chunk plus the
  source FASTQs.

The M2-era memory work proved per-chunk trimming is byte-identical to
whole-side trimming. M4 only moves the chunk boundary between processes.

**Architecture:**
- **Three new subcommands,** thin wrappers over the existing functions:
  - `te-split`: per read side, write chunk FASTAs plus a manifest
    (side, chunk index, file, n_reads).
  - `te-search-chunk`: run `_blat_chunk` on one or more chunk FASTAs, writing
    chunk BAMs.
  - `te-merge-trim`: given every chunk BAM in manifest order plus the original
    reads, merge each side BAM (`_merge_sorted_bams`) and run
    `write_trimmed_reads_chunked`. The outputs are the files `relocaTE3 run`
    writes.
- **Refactor `run`** (`identify_TE_reads`) so the BLAT path calls the same three
  functions in-process. Then `run`, `run-all` and the scatter path share one
  implementation.
- **Pipeline:** when `te_aligner == 'blat'` and scatter is enabled:
  TE_SPLIT → TE_SEARCH_CHUNK (one task per group of `blat_chunks_per_task`
  chunks) → TE_MERGE_TRIM. Otherwise the existing single `RELOCATE3_TE_SEARCH`
  runs. Downstream steps see the same files.
- **New params:**
  - `blat_scatter` (default **true**; decide after measuring in Task 6);
  - `blat_chunks_per_task` (default chosen from measurement);
  - hidden `blat_chunk_reads`, default 200000 = `BLAT_CHUNK_SEQS`, so tests can
    force many chunks on small data.

**Branch / worktree:** `feat/nextflow-m4` in `../RelocaTE3-nextflow-m4`
(`--no-track`), from `origin/main` a3db2bd, which merges PR #65.
Push with `git push -u origin feat/nextflow-m4`.

---

## Conventions (as M1–M3)

**Tests:**
- pytest in the `pipeline` env:
  `E=$PWD/.pixi/envs/pipeline/bin; PATH=$E:$PATH PYTHONPATH=$PWD/src $E/python -m pytest tests -p no:cacheprovider --junitxml=<scratch>/x.xml`.
  Read the counts from the XML.
- `pixi run -e nf nf-test test`
- `pixi run -e nf nextflow lint .`
- **Baseline at a3db2bd:** pytest 614, nf-test 41. Record it in Task 0.

**Process:**
- TDD.
- sbatch for heavy work: `#!/usr/bin/bash -l`, `set -euo pipefail`, `-p intel`,
  `logs/<name>.%j.log`, date and parameters at the start.
- Never on the login node. Wait with a `timeout … until sacct …` loop.
- Commit per task, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Never push.

**Core rule:** no detection-logic change; outputs identical. Any difference is a
bug to find, not a tolerance to add.

**30x data** (simulated riceTElib panel; the paths are in
`validation/blat_memory/run_step3_memcheck.sh`):
`PANEL=/bigdata/wesslerlab/shared/Rice/Nathan/rice/make_simulated_genome/make_simulation_new/results/riceTElib_benchmark`
- reads: `${PANEL}/reads/cov30x_rep1/R{1,2}.fastq.gz`
- TE library: `${PANEL}/selected_te_library.fa`
- Find the genome FASTA in the panel and record it.

---

### Task 0: Setup
- `pixi install -e pipeline -e nf` (plus `-e blat` if tests need it).
- Run the baseline suites and record the counts.

### Task 1: Split the BLAT path into three callable steps (library refactor)

**Files:** `src/RelocaTE3/aligners.py`, `src/RelocaTE3/librelocate.py`, tests.

- **Public functions,** keeping behaviour identical:
  - `split_te_query(read_file, outdir, stem, chunk_reads=BLAT_CHUNK_SEQS) -> list[(path, n_reads)]`
  - `search_te_chunk(te_library, chunk_fa, chunk_bam, mate, opts) -> TEChunk`.
    It is `_blat_chunk`. The SAM header comes from the TE library, so it can be
    computed independently.
  - `merge_trim(...)`: merges each side from its chunk BAMs, then runs
    `write_trimmed_reads_chunked`.
  - All three must use the same BLAT options (`RT2_BLAT_OPTS` plus `te_opts`)
    and the same mate canonicalisation.
- **Manifest:** a small TSV written by `te-split` and read by `te-merge-trim`.
  Columns: side, index, chunk FASTA name, n_reads, mate. Validation:
  - chunk order is preserved by the zero-padded index;
  - `te-merge-trim` refuses a missing or extra chunk BAM;
  - Σ n_reads must equal the reads in the side's FASTQ;
    `write_trimmed_reads_chunked` already checks this.
- `identify_TE_reads` uses these functions for BLAT, still with its in-process
  thread pool, so `relocaTE3 run` is unchanged.
- **Tests:**
  - the existing BLAT tests stay green;
  - a new test: the three functions, run one by one, give byte-identical
    `te_containing/`, `flanking/` and `te_portions/` and record-identical BAMs
    compared with `identify_TE_reads`, on the Chr3 fixture with
    `chunk_reads` = a small number (e.g. 5000) so there are many chunks;
  - a test that `chunk_reads` doesn't change outputs (5000 vs 200000).
- **Gate:** `sbatch tests/nextflow/parity_chr3.sh` still PASSES.

Commit: "BLAT TE search as three callable steps: split, search chunk, merge+trim".

### Task 2: CLI subcommands

`relocaTE3 te-split`, `te-search-chunk` and `te-merge-trim`, in the existing
argparse style.

- **Arguments:** reads `-1/-2`, `-T` TE library, `-n`, `-o .`, `--threads`,
  `--te-opts`, `--min-match`, `--min-trimmed`, `--mismatch`,
  `--chunk-reads` (hidden or documented).
- **Output layout:** `te-merge-trim -o .` must write exactly the files
  `relocaTE3 run -o .` writes (M1 `file_list`): `<name>.{left,right}.bam(.bai)`,
  `te_containing/` (4 files), `flanking/` (2), `te_portions/` (2). No leftover
  chunk directories.
- **Tests:**
  - parser tests;
  - an end-to-end test on Chr3 with small `--chunk-reads`, where the
    subcommands chained in temporary directories match
    `relocaTE3 run`'s outputs exactly (`compare_outputs.sh` style);
  - an error test: a missing chunk BAM is refused.
- `tests/nextflow_modules_args_test.py`-style guard: every argument `cmd_run`
  takes that affects the BLAT search or trim must also be accepted by
  `te-search-chunk`/`te-merge-trim`.

Commit: "te-split / te-search-chunk / te-merge-trim subcommands".

### Task 3: Nextflow modules and wiring

- **Modules:** `RELOCATE3_TE_SPLIT`, `RELOCATE3_TE_SEARCH_CHUNK` and
  `RELOCATE3_TE_MERGE_TRIM`, with versions and stubs.
- **Grouping:** `RELOCATE3_TE_SEARCH_CHUNK` takes a group of chunk FASTAs and
  runs them with `task.cpus` threads. Group by the manifest index in batches of
  `params.blat_chunks_per_task`, per sample and side (`buffer`/`collate`).
  Keep it deterministic for `-resume`, with sorted inputs.
- **Labels:**
  - TE_SEARCH_CHUNK: about 4 CPUs, 4 GB. Measure: one 200k chunk at RT2 BLAT
    options takes about 4 min per CPU at 30x riceTElib (verify in Task 6).
  - TE_MERGE_TRIM: about 4 GB, with trim memory bounded per chunk.
- **Wiring** in `subworkflows/local/one_sample.nf`:
  `params.te_aligner == 'blat' && params.blat_scatter` → the scatter branch;
  otherwise TE_SEARCH. Both branches emit the same channels (`te_containing`,
  `flanking`, `read_repeat`, …) so downstream is untouched.
- **Params** `blat_scatter`, `blat_chunks_per_task` and `blat_chunk_reads`
  (hidden) go into the config, the schema and `tests/nextflow_schema_test.py`.
  These are pipeline-only params.
- `conf/modules.config`: publish TE_MERGE_TRIM's outputs exactly where
  TE_SEARCH's go. The per-sample 33-file invariant must hold. Do not publish
  chunk files.
- **nf-test:**
  - stub: chunk tasks = ceil(chunks / per_task) per side, merge once per sample;
  - `blat_scatter false` gives the old single task;
  - a non-blat aligner never scatters.
- **Parity:** `tests/nextflow/parity_chr3.sh` gains a way to force scatter with
  small chunks (e.g. env `NXF_EXTRA="--blat_chunk_reads 5000 --blat_chunks_per_task 2"`)
  and must PASS both ways.

Commit: "Pipeline: BLAT scatter (split → search chunks → merge+trim)".

### Task 4: Chr3 gates in pixi and in the container
- `sbatch tests/nextflow/parity_chr3.sh` (default), plus the forced-scatter
  variant: both PASS.
- CI: the `container.yml` pipeline job runs the forced-scatter parity too
  (cheap).

### Task 5: Docs (before the long run)
- README: what scatter does, the params, when to turn it off.
- CHANGELOG.
- `docs/2026-11-05-nextflow-m4.md`.

### Task 6: 30x riceTElib gate and timing (long; sbatch)
1. **Unsplit baseline:**
   `relocaTE3 run -1 R1 -2 R2 -T selected_te_library.fa -n cov30x_rep1 -o <base> --threads 8`
   in one sbatch job (about 13 h, 16 GB, 8 CPUs; the M2 measurement was
   13:03:27). Record the wall time.
2. **Scatter:** the pipeline's TE-search part only, or the full pipeline if
   the genome is available, with `-profile hpcc,pixi`. Running the full
   pipeline is better; it also gives the downstream calls. Use default
   `blat_chunks_per_task`.
   - Before choosing the default, try 2–3 values on the first side, or measure
     per-chunk time from the trace and pick a value giving tasks of about
     20–40 min.
   - Record the end-to-end TE-search wall time and the queue/run times from the
     trace.
3. **Compare:** `te_containing/`, `flanking/`, `te_portions/` and the side
   BAMs between the baseline and the scatter run, with `compare_outputs.sh`.
   Required: identical. If the full pipeline ran, also compare against
   `run-all` on the same sample. Running run-all would need another ~14 h, so
   that's optional; the TE-search identity plus the Chr3 full parity suffices.
4. Record everything in the docs note: wall times, task counts, peak memory and
   the chosen defaults.

### Task 7: Review and hand-off
- Final whole-branch review; fix Important findings.
- PR message at `../RelocaTE3/results/pr_messages/nextflow-m4.md`.
