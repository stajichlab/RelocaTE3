# RelocaTE3

RelocaTE3 identifies transposable element (TE) insertion polymorphisms from
short-read resequencing data, at single-base resolution, by comparison to a
reference genome. It is a modern, pure-Python reimplementation of
[RelocaTE2](https://github.com/JinfengChen/RelocaTE2) that depends only on
`minimap2` and `samtools` (plus `pysam`/`biopython`) — no blat, bwa, bowtie2,
seqtk, or Perl.

> Status: reference implementation. The full pipeline (read trimming → genome
> re-alignment → non-reference insertion calling → reference/shared insertions →
> genotyping) runs end to end. On the rice Chr3 2 Mb benchmark it recovers
> ~178/200 simulated insertions (~89% recall, ~90% precision); see
> `tests/acceptance_test.py`. See `PLAN.md` for the roadmap (Rust acceleration and
> Nextflow scatter are planned).

## Installation

With [pixi](https://pixi.sh) (recommended — pins `minimap2`/`samtools`/`bedtools`):

```bash
pixi install
pixi run relocaTE3 --help
```

Or into an existing environment that already provides `minimap2` and `samtools`:

```bash
pip install -e .
```

## Usage

Run the whole pipeline for one paired-end sample:

```bash
relocaTE3 run \
  --r1 reads_1.fq.gz --r2 reads_2.fq.gz \
  --te RiceTE.fa \
  --genome reference.fa \
  --outdir HEG4_out --sample HEG4 \
  --threads 8 --mismatch 2 \
  --repeatmasker reference.fa.RepeatMasker.out \   # optional: reference/shared calls + FP filter
  --genotype                                       # optional: homo/heterozygous/somatic
```

### Individual steps

Each step is also a subcommand, so a workflow engine can scatter them:

| Command | Step | Purpose |
|---------|------|---------|
| `trim` | 3 | Map reads to the TE library, trim the TE portion, emit flanking reads |
| `align-genome` | 4 | Re-align flanking reads + supporting mates to the genome |
| `find-insertions` | 5 | Cluster junction/supporting reads → non-reference insertions |
| `find-reference` | 0/6 | Parse RepeatMasker `.out`, call reference/shared insertions |
| `characterize` | 7 | Genotype insertions (homozygous/heterozygous/somatic) |
| `run` | all | Full single-sample pipeline (the steps above) |

Run `relocaTE3 <command> --help` for options.

## Inputs

- **Reads**: one (single-end) or two (paired-end) FASTQ files (`.fq`/`.fq.gz`).
- **TE library** (`--te`): FASTA of TE/repeat consensus sequences.
- **Genome** (`--genome`): reference genome FASTA.
- **RepeatMasker `.out`** (`--repeatmasker`, optional): reference TE annotation,
  enabling reference/shared insertion calls and false-positive filtering.

## Outputs

Under `--outdir`:

```
flanking/        trimmed flanking reads (FASTQ)
te_portions/     TE-matching read portions (FASTA)
te_containing/   read → TE assignment table
genome_aln/      flanking/support/full reads aligned to the genome (BAM)
existingTE.bed   reference TE annotation (when --repeatmasker given)
results/
  <sample>.all_nonref_insert.gff / .txt          non-reference insertions
  <sample>.all_ref_insert.gff / .txt             reference/shared insertions
  <sample>.all_nonref_insert.characTErized.gff   genotyped insertions (--genotype)
```

The non-reference GFF carries the RelocaTE2 attribute set: `TSD`, `Name` (TE
family), `Note`, and `Left/Right_junction_reads` and `Left/Right_support_reads`.

## Migrating from RelocaTE2

| RelocaTE2 | RelocaTE3 |
|-----------|-----------|
| `--fq_dir` (directory) | `--r1` / `--r2` (explicit files) |
| `--te_fasta` | `--te` |
| `--genome_fasta` | `--genome` |
| `--reference_ins` | `--repeatmasker` |
| `--mismatch` | `--mismatch` (default 0; use 2 to match the RelocaTE2 benchmark) |
| `--aligner blat/bwa/bowtie2` | minimap2 only |
| `characterizer.pl` (Perl) | `characterize` / `--genotype` |

## Development

```bash
pixi run test          # run the test suite
pixi run pytest tests/acceptance_test.py   # the benchmark acceptance gate
```
