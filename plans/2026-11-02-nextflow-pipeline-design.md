# Nextflow pipeline: design

2026-11-02. Approved in conversation with the maintainer, section by section.
This supersedes the sketch in `plans/FEATURES.md` (Feature 2).

Roadmap order (maintainer, 2026-11-02):
1. Nextflow (this document).
2. Long reads.
3. Multiple reference genomes as separate FASTA files.
4. Pangenome reference.

Rust acceleration moves after these.

## Decisions

| Question | Decision |
|---|---|
| Where the pipeline lives | In the RelocaTE3 repo with nf-core layout and tooling. Not submitted to nf-core. Pipeline version = RelocaTE3 version. |
| Software | One container image built from `pixi.lock` (validated tool versions), at `ghcr.io/stajichlab/relocate3:<version>`. A `pixi` profile for container-free runs. |
| Parallelism | Staged. M1: one task per step per sample, with genome indexing shared. M4: BLAT chunks spread across jobs. Per-chromosome splitting only if needed. |
| Outputs | Per-sample results identical to `run-all`, a cohort table, a genotype matrix and a MultiQC summary. A multi-sample VCF is deferred. |

## Principle

The workflow contains no detection logic. Every process runs one `relocaTE3`
subcommand with the arguments `run-all` uses (`cli.py` `cmd_run_all`), so a
sample run through Nextflow produces the same files as `run-all`. RelocaTE2
parity then carries over.

Lesson from the 0.2.0 release check: the real-rice harness ran the steps by hand
with different settings and lost about 1,270 calls. Every milestone is therefore
gated on identical output against `run-all`.

## Layout

```
main.nf                      entry point
nextflow.config              params, profiles, manifest (version = package version)
nextflow_schema.json         parameter validation + --help (nf-schema plugin)
assets/schema_input.json     sample-sheet validation
workflows/relocate3.nf       sample loop, genome preparation, cohort stage
subworkflows/local/          prepare_genome.nf, one_sample.nf
modules/local/relocate3_*/   one process per subcommand: main.nf + tests/ (nf-test)
conf/                        base.config (labels/resources), test.config, hpcc.config
containers/Dockerfile        image from pixi.lock (`pipeline` env = blat env + MultiQC)
tests/nextflow/              pipeline-vs-run-all parity test
```

## CLI additions

Each is built test-first. `run-all` switches to using them, so the two paths
can't diverge.

1. `index-genome --genome-aligner <a>`. Today it builds only the `.fai` and the
   minimap2 `.mmi` index. It will also build the chosen genome aligner's index,
   which `align-genome` currently builds on demand next to the genome.
2. `align-reads`: the genotyping alignment of all original reads. `run-all`
   currently builds it inline (`Aligner.map_library_to_genome`).
3. `cohort` (M3): combines characterized tables into one TSV with a sample
   column, plus a sites-by-samples genotype matrix keyed by exact chromosome,
   start, end, strand and TE family. Merging nearby sites is deferred to the
   VCF design.

## Data flow

```
sheet -> MERGE_LANES (only when a sample has several rows)
genome -> PREPARE_GENOME (once: faidx, minimap2 index, genome-aligner index)
per sample: TE_SEARCH (run) -> ALIGN_GENOME -> FIND_INSERTIONS -> CHARACTERIZE
                                            \-> FIND_REFERENCE (with --repeatmasker)
            ALIGN_READS (genotyping alignment, if no genotype_bam), in parallel
                        with the chain above
cohort:     COHORT_TABLE -> SUMMARY (MultiQC)
```

- Each process stages the previous task's subfolders (`te_containing/`,
  `flanking/`, `genome_aln/`, `results/`) under their usual names and runs the
  subcommand with `-o .`. No subcommand changes are needed for this.
- Published outputs:
  - `results/<sample>/`: same layout as `run-all -o <sample>`.
  - `results/cohort/`
  - `results/pipeline_info/`: params, image digest, versions, Nextflow
    report, timeline and trace.
- Each process writes `versions.yml` (nf-core convention) and a small counts
  file for MultiQC.
- Resource labels. TE search: 8 CPUs, 8 GB (the 30x peak is 3.4 GB after
  chunked trimming). Retry with more memory on exit 137/140.

## Inputs and parameters

- **Sample sheet:** `sample,fastq_1,fastq_2[,genotype_bam]`. Rows that share a
  sample name are lanes, merged before the TE search.
- **Shared inputs:** `--genome`, `--te_library`, `--repeatmasker` (optional),
  `--outdir`.
- **Calling parameters**, with the same names and defaults as `run-all`:
  `--te_aligner blat`, `--genome_aligner bwaaln`, `--mismatch 2`,
  `--min_match 10`, `--min_trimmed 10`, `--size 500`, `--min_mapq 0`,
  `--tsd UNK` (schema enum: UNK only), `--require_both_junctions true`,
  `--target ALL`.
- **Genotyping:** `--genotype` defaults to **true** in the pipeline, so the
  genotype matrix is filled. `run-all` makes it opt-in.
- Everything can be set with `-params-file`. The schema serves `--help` and
  machine-readable parameter descriptions for agent users.

## Profiles

| Profile | Use |
|---|---|
| `docker` / `apptainer` / `singularity` | run in the image |
| `pixi` | no container; launch inside `pixi run -e pipeline` |
| `test` | Chr3 2 Mb fixture (`tests/data/sim_*`, `TE_lib`), minutes |
| `hpcc` | UCR SLURM: `--slurm_partition` (default `intel`), queue size, Apptainer cache, time limits |

Minimum Nextflow version: **25.10** (changed during M1 Task 3). nf-schema ≥ 2.6
needs it, and earlier releases predate fixes for typed CLI parameters under
strict syntax. The cluster provides 25.10.4 and 26.04.x; the `nf` pixi env
has 26.10.4.

## Testing

1. pytest, failing first, for each CLI addition. `run-all` output must be
   unchanged after switching to them.
2. nf-test per module: a stub run, plus a real run on tiny data.
3. **Parity gate:** on the Chr3 fixture, the pipeline with `-profile test`
   matches `run-all` with the same parameters. Text files match byte for byte;
   BAMs match record by record via `samtools view`.
4. Real rice (M2): 10 samples with `-profile hpcc,apptainer`, identical to the
   v0.2.0rc2 results.
5. CI: pytest, nf-test, `nextflow lint`, and `-profile test,docker` on every
   pull request.

## Milestones

| # | Delivers | Passes when |
|---|---|---|
| M1 | CLI additions 1–2; per-sample processes; `pixi` and `test` profiles | Chr3 pipeline output identical to `run-all` |
| M2 | Image and CI publishing; container and `hpcc` profiles | Chr3 match in the container; 10 real-rice samples identical to rc2 |
| M3 | `cohort` and MultiQC summary | cohort rows match per-sample files exactly |
| M4 | BLAT chunk scatter across jobs (new chunk-search and merge-trim subcommands) | identical to unsplit; 30x wall time measured |
| M5 | Usage docs; real-rice harness on the pipeline; release `v0.3.0` | full release check passes |
