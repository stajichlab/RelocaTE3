# RelocaTE3

RelocaTE3 identifies transposable element (TE) insertion polymorphisms from
short-read resequencing data, at single-base resolution, by comparison to a
reference genome. It is a modern, pure-Python reimplementation of
[RelocaTE2](https://github.com/JinfengChen/RelocaTE2) built on `minimap2`
(default) and `samtools` (plus `pysam`/`biopython`) — no seqtk or Perl. The
aligner for each stage is selectable: `minimap2`, `bwa`, `bwa-mem2`, `bowtie2`,
or `blat` (see "Choosing an aligner").

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

### Quick start: one command

`run-all` executes the whole pipeline for a single sample:

```bash
relocaTE3 run-all \
  -1 reads_1.fq.gz -2 reads_2.fq.gz \
  -T RiceTE.fa \
  -g reference.fa \
  -n HEG4 -o HEG4_out \
  --threads 8 --mismatch 2 --tsd UNK \
  --repeatmasker reference.fa.RepeatMasker.out \
  --genotype
```

That indexes the genome (if needed), finds and trims TE-containing reads, places
the flanks, calls non-reference insertions, and — because `--repeatmasker` and
`--genotype` were given — also emits reference/shared calls and genotypes every
site.

`run-all` dispatches the *same* handlers as the staged subcommands below, so its
results match the staged workflow. Pick whichever fits: `run-all` for a laptop or
a single HPC node, the staged commands when a workflow engine needs to scatter
work across samples or chromosomes.

To reproduce the configuration validated against RelocaTE2 on the rice benchmark,
add `--te-aligner blat --genome-aligner bwaaln --min-mapq 1
--require-both-junctions`.

### Staged subcommands

Note that `run` performs **TE-read identification and flank
generation only** (map + trim, steps 2-3) — it is *not* the complete pipeline.
A full single-sample analysis chains the stages:

```bash
# 1. one-time: index the reference genome
relocaTE3 index-genome -g reference.fa

# 2-3. identify TE-containing reads and emit flanking reads (map + trim)
relocaTE3 run \
  --left reads_1.fq.gz --right reads_2.fq.gz \
  --te-library RiceTE.fa \
  --name HEG4 --outdir HEG4_out \
  --threads 8 --mismatch 2

# 4. re-align the trimmed flanking reads to the genome
relocaTE3 align-genome \
  -g reference.fa \
  -f HEG4_out/flanking/HEG4.left.flankingReads.fq HEG4_out/flanking/HEG4.right.flankingReads.fq \
  -n HEG4 -o HEG4_out --threads 8

# 5. cluster junction/supporting reads into non-reference insertions
relocaTE3 find-insertions \
  -b HEG4_out/HEG4.repeat.minimap.sorted.bam \
  --read-repeat HEG4_out/te_containing/HEG4.read_repeat_name.txt \
  --tsd TTA --target ALL --name HEG4 --outdir HEG4_out --te-name mping \
  --reference-ins reference.fa.RepeatMasker.out \
  --mismatch 2 --min-mapq 1

# 7. genotype the insertions from a reads-to-genome BAM/CRAM
relocaTE3 characterize \
  -s HEG4_out/results/ALL.mping.all_nonref_insert.txt \
  -b reads_to_genome.bam -g reference.fa \
  -o HEG4_out/results
```

### Subcommands

| Command | Step | Purpose |
|---------|------|---------|
| `run-all` | 0-7 | **Whole pipeline for one sample in one command** (dispatches the staged handlers below) |
| `index-genome` | 1 | Index the reference genome (`samtools faidx` + minimap2) |
| `map` | 2 | Align reads to the TE library (produces per-side BAMs) |
| `trim` | 3 | Trim the TE portion from TE-library BAMs, emit flanking reads |
| `run` | 2-3 | **Map + trim** in one step (TE-read identification + flank generation) — not the full pipeline |
| `annotate-ref` | 0 | Annotate where reference TE copies *are* (minimap2 → `existingTE.bed`), to filter novel calls |
| `find-reference` | 0/6 | Call reference TEs that are **also present in this sample** → `all_ref_insert.{gff,txt}` |
| `align-genome` | 4 | Re-align flanking reads to the genome |
| `find-insertions` | 5 | Cluster junction/supporting reads → non-reference insertions |
| `characterize` | 7 | Genotype insertions (homozygous/heterozygous/somatic; `-x` for excision) |

Run `relocaTE3 <command> --help` for the full flag list.

## Inputs

- **Reads** (`-l/--left/--r1`, `-r/--right/--r2`): one (single-end) or two
  (paired-end) FASTQ files (`.fq`/`.fq.gz`).
- **TE library** (`-T/--te-library`): FASTA of TE/repeat consensus sequences.
- **Genome** (`-g/--genome-fasta`): reference genome FASTA.
- **Existing-TE annotation** (`find-insertions --reference-ins`, optional): a
  RepeatMasker `.out` or BED of reference TE copies, used to skip/flag reference
  insertions.
- **Characterize alignments** (`characterize -b`): sorted, indexed BAM or CRAM of
  the original reads aligned to the genome. CRAM input also requires `-g` (the
  genome FASTA).

## Outputs

Under `--outdir` (example for `--name HEG4`, and `find-insertions --te-name mping
--target ALL`):

```
HEG4.left.bam / HEG4.right.bam            reads aligned to the TE library (map)
HEG4.repeat.minimap.sorted.bam           flanking reads aligned to the genome (align-genome)
flanking/
  HEG4.left.flankingReads.fq             trimmed flanking reads (5'/3')
  HEG4.right.flankingReads.fq
te_containing/
  HEG4.read_repeat_name.txt              read → TE-family assignment table
  HEG4.left.ContainingReads.fq / HEG4.right.ContainingReads.fq
te_portions/
  HEG4.five_prime.fa / HEG4.three_prime.fa   TE-matching read portions
results/
  ALL.mping.all_nonref_insert.txt                        non-reference insertions (find-insertions)
  ALL.mping.all_nonref_insert.characTErized.gff / .txt   genotyped insertions (characterize)
  HEG4.all_ref_insert.gff / .txt                         reference/shared insertions (find-reference)
existingTE.bed                                           reference TE copies (find-reference/annotate-ref)
```

> **Reproducibility note.** Reads that match several TE families equally well can
> be assigned to different families on repeated runs, so a small number of calls
> may change their `Name`/family label between otherwise identical runs. Positions
> and read counts are stable. This affects the staged and `run-all` paths alike.

The characterized GFF carries the RelocaTE2 attribute set: `TSD`, `Name` (TE
family), `Note`, and `Left/Right_junction_reads` and `Left/Right_support_reads`.

## Migrating from RelocaTE2

| RelocaTE2 | RelocaTE3 |
|-----------|-----------|
| `--fq_dir` (directory) | `--left` / `--right` (explicit files) |
| `--te_fasta` | `-T/--te-library` |
| `--genome_fasta` | `-g/--genome-fasta` |
| `--reference_ins` | `find-insertions --reference-ins` |
| `--mismatch` | `--mismatch` (default 0; use 2 to match the RelocaTE2 benchmark) |
| `--aligner blat/bwa/bowtie2` | `run --te-aligner` and `align-genome --genome-aligner` (see below) |
| `characterizer.pl` (Perl) | `characterize` |
| whole pipeline in one command | `run-all` |
| `all_ref_insert.gff/.txt` | `find-reference` (or `run-all --repeatmasker`) |

`run-all` is the closest equivalent to a single RelocaTE2 invocation. The staged
subcommands remain available for workflow engines.

## Choosing an aligner

The aligner is selectable per stage:

- **TE-library search** (`map` / `run`): `--te-aligner {minimap2,bwa,bwamem2,bowtie2,blat}`
  (default `minimap2`). `--aligner` is a deprecated alias.
- **Genome re-alignment** (`align-genome`): `--genome-aligner {minimap2,bwa,bwamem2,bowtie2}`
  (default `minimap2`). `blat` is TE-search only and is rejected here.

```bash
relocaTE3 run --left r1.fq --right r2.fq -T RiceTE.fa -n HEG4 -o HEG4_out --te-aligner bwa
relocaTE3 align-genome -g reference.fa -f HEG4_out/flanking/*.flankingReads.fq -n HEG4 -o HEG4_out --genome-aligner bwa
```

`bwa`, `bwa-mem2`, and `bowtie2` are pinned by pixi. `blat` is optional — provide it
on `PATH` (its bioconda build conflicts with the pinned plotting stack). Non-minimap2
genome BAMs are named `{name}.repeat.{aligner}.sorted.bam`.

## Development

```bash
pixi run test          # run the test suite
pixi run pytest tests/acceptance_test.py   # the benchmark acceptance gate
```

### UCR HPCC User Development

The following commands run RelocaTE3 on real rice data and compares results
to legacy RelocaTE2 results. See `validation` directory for more information.

```bash
pixi run validate-rice --local --force B_10   # smoke test (one sample, forced re-run)
pixi run validate-rice                        # full 10-sample SLURM array
```
