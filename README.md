# RelocaTE3

RelocaTE3 identifies transposable element (TE) insertion polymorphisms from
short-read resequencing data, at single-base resolution, by comparison to a
reference genome. It is a modern, pure-Python reimplementation of
[RelocaTE2](https://github.com/JinfengChen/RelocaTE2) built on `samtools`
(plus `pysam`/`biopython`) — no seqtk or Perl. By default it uses the same
aligner pair as RelocaTE2, `blat` for TE search and `bwa aln` for genome
placement, but the aligner for each stage is selectable: `minimap2`, `bwa`,
`bwa-mem2`, `bwa aln`, `bowtie2`, or `blat` (see "Choosing an aligner").

> Status: reference implementation. The full pipeline (read trimming → genome
> re-alignment → non-reference insertion calling → reference/shared insertions →
> genotyping) runs end to end. On the rice Chr3 2 Mb benchmark, at the shipped
> defaults, it recovers 194/200 simulated insertions (97% recall, 90%
> precision) — RelocaTE2's published figure on the same data is 196/200. See
> `tests/acceptance_test.py`, which gates the minimap2 configuration
> (~178/200) so it can run without blat. See `PLAN.md` for the roadmap (Rust acceleration and
> Nextflow scatter are planned).

## Installation

With [pixi](https://pixi.sh) (recommended — pins the aligners/`samtools`/`bedtools`):

```bash
pixi install
pixi run relocaTE3 --help
```

The default TE aligner is `blat`, which cannot be installed alongside the
plotting stack (its bioconda build hard-pins `zlib 1.2.11`, and it is linux-64
only). It therefore lives in a separate pixi environment that also carries
RelocaTE3 itself:

```bash
pixi run -e blat relocaTE3 --help     # linux-64; blat + RelocaTE3 together
```

Any `blat` on `PATH` works just as well — the backend simply shells out to it.
Without one, pass `--te-aligner bowtie2` (the closest alternative on the rice
benchmark); RelocaTE3 will say so if `blat` is missing.

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
  --threads 8 \
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

**Every default matches RelocaTE2**, so the command above reproduces a stock
RelocaTE2 run without tuning — see "RelocaTE2 defaults" below for the mapping.

### Many samples

`run-batch` runs the same pipeline across a cohort, one `run-all` per sample into
`<outdir>/<sample>`:

```bash
# a directory of paired FASTQs (RelocaTE2's --fq_dir)
relocaTE3 run-batch --fq-dir reads/ -T RiceTE.fa -g reference.fa \
  -o results --threads 8 --jobs 4 --repeatmasker reference.fa.RepeatMasker.out

# or an explicit sample sheet
relocaTE3 run-batch --samples samples.csv -T RiceTE.fa -g reference.fa -o results
```

`--fq-dir` pairs `<sample>_R1`/`_R2` and `<sample>_1`/`_2`. The sheet is CSV or
TSV with `sample_id,r1_fq[,r2_fq]`, plus optional per-row `te_library`,
`reference_genome` and `repeatmasker` columns that override the global flags:

```csv
sample_id,r1_fq,r2_fq
HEG4,reads/HEG4_R1.fq.gz,reads/HEG4_R2.fq.gz
A123,reads/A123_R1.fq.gz,reads/A123_R2.fq.gz
```

`--jobs` sets how many samples run at once (each still uses `--threads`).
Batches are **resumable**: a sample that already has a results table is skipped
unless `--force`. By default the batch stops at the first failing sample; pass
`--keep-going` to run the rest and get a summary at the end (the exit code is
non-zero either way if anything failed).

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
  --threads 8

# 4. re-align the trimmed flanking reads to the genome
relocaTE3 align-genome \
  -g reference.fa \
  -f HEG4_out/flanking/HEG4.left.flankingReads.fq HEG4_out/flanking/HEG4.right.flankingReads.fq \
  -n HEG4 -o HEG4_out --threads 8

# 5. cluster junction/supporting reads into non-reference insertions
relocaTE3 find-insertions \
  -b HEG4_out/HEG4.repeat.bwaaln.sorted.bam \
  --read-repeat HEG4_out/te_containing/HEG4.read_repeat_name.txt \
  --target ALL --name HEG4 --outdir HEG4_out --te-name mping \
  --reference-ins reference.fa.RepeatMasker.out

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
| `run-batch` | 0-7 | **Whole pipeline for many samples** from a sample sheet or FASTQ directory |
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
  (paired-end) FASTQ files (`.fq`/`.fq.gz`). Any read-name convention works --
  `/1`-`/2` suffixes, modern Illumina (`@A00519:... 1:N:0:INDEX`, where the name
  is identical in both files), or SRA/ENA dumps with no mate marker. Which file
  a read came from determines its mate, so the name does not have to say.
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
HEG4.repeat.bwaaln.sorted.bam            flanking reads aligned to the genome (align-genome;
                                          named for the aligner, minimap2 abbreviates to "minimap")
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

> **Reproducibility.** The same inputs and the same version produce the same
> result tables. Reads matching several TE families equally well, and tied
> cluster-level family votes, are resolved by a fixed rank rather than by
> whichever alignment the aligner emitted first, so results do not shift between
> runs. Note that intermediate BAMs are written by multi-threaded aligners and
> may still differ byte-for-byte between runs; the tables under `results/` are
> the reproducible artifacts.

The characterized GFF carries the RelocaTE2 attribute set: `TSD`, `Name` (TE
family), `Note`, and `Left/Right_junction_reads` and `Left/Right_support_reads`.

## Migrating from RelocaTE2

| RelocaTE2 | RelocaTE3 |
|-----------|-----------|
| `--fq_dir` (directory) | `run-batch --fq-dir`, or `--left`/`--right` for one sample |
| `--te_fasta` | `-T/--te-library` |
| `--genome_fasta` | `-g/--genome-fasta` |
| `--reference_ins` | `find-insertions --reference-ins` |
| `--mismatch` / `--mismatch_junction` | `--mismatch` (one knob covers both; default 2, same as RelocaTE2) |
| `--aligner blat/bwa/bowtie2` | `run --te-aligner` and `align-genome --genome-aligner` — the defaults (`blat`, `bwa aln`) already match RelocaTE2 |
| `characterizer.pl` (Perl) | `characterize` |
| `--mate_1_id` / `--mate_2_id` / `--unpaired_id` | not needed — the mate is taken from which file the read came from |
| whole pipeline in one command | `run-all` |
| `all_ref_insert.gff/.txt` | `find-reference` (or `run-all --repeatmasker`) |

`run-all` is the closest equivalent to a single RelocaTE2 invocation. The staged
subcommands remain available for workflow engines.

## RelocaTE2 defaults

RelocaTE3 ships RelocaTE2's defaults, so an untuned run reproduces RelocaTE2's
behaviour. Values below are cited from RelocaTE2's source.

| Parameter | RelocaTE2 | RelocaTE3 | Source |
|---|---|---|---|
| TE-search aligner | `blat` | `--te-aligner blat` | `relocaTE2.py:204` |
| Genome placement | `bwa aln` | `--genome-aligner bwaaln` | `relocaTE_align.py` |
| Mismatches (reads vs TE) | `--mismatch 2` | `--mismatch 2` | `relocaTE2.py:207` |
| Mismatches (junction reads) | `--mismatch_junction 2` | same `--mismatch` | `relocaTE2.py:208` |
| Match-length cutoff | `--len_cut_match 10` | `--min-match 10` | `relocaTE2.py:205` |
| Trimmed-length cutoff | `--len_cut_trim 10` | `--min-trimmed 10` | `relocaTE2.py:206` |
| TSD | `UNK` (hardcoded) | `--tsd UNK` | `relocaTE2.py:346` |
| Junction reads required | `left >= 1` **or** `right >= 1` | same; `--require-both-junctions` is opt-in | `relocaTE_insertionFinder.py:365,1732-4` |
| Threads | `--cpu 1` | `--threads 1` | `relocaTE2.py:199` |

Two RelocaTE2 options have no RelocaTE3 equivalent because they are no longer
needed: `--mate_1_id`/`--mate_2_id`/`--unpaired_id` (the mate comes from which
file a read is in) and `--split` (chunking for blat/bwa).

`--size` (insert size, default 500) is not implemented: RelocaTE2 uses it only
to estimate the span of insertions supported by mate pairs alone, which
RelocaTE3 does not yet emit as a separate output.

`--min-mapq` (default 1) has no RelocaTE2 counterpart — RelocaTE2 instead
classifies reads below MAPQ 29 as low-quality using bwa's `XM`/`X1`/`XO` tags
(`relocaTE_insertionFinder.py:1523`).

## Choosing an aligner

The aligner is selectable per stage:

- **TE-library search** (`map` / `run`): `--te-aligner {minimap2,bwa,bwamem2,bwaaln,bowtie2,blat}`
  (default `blat`, matching RelocaTE2). `--aligner` is a deprecated alias.
- **Genome re-alignment** (`align-genome`): `--genome-aligner {minimap2,bwa,bwamem2,bwaaln,bowtie2}`
  (default `bwaaln`, matching RelocaTE2). `blat` is TE-search only and is
  rejected here.

Measured on the riceTElib benchmark (9 samples x 3 coverages), matched true
calls and precision at 30x:

| TE search / genome | matched 5x | 15x | 30x | precision 30x |
|---|---|---|---|---|
| `blat` / `bwaaln` (default) | 455 | 781 | 962 | 0.845 |
| `bowtie2` / `bwa` | 428 | 761 | 946 | 0.848 |
| `minimap2` / `minimap2` | 346 | 632 | 798 | 0.827 |
| RelocaTE2 | 449 | 739 | 897 | 0.834 |

`bowtie2` is the best option when `blat` is unavailable — it costs ~2% of the
matched calls and needs no extra environment.

```bash
relocaTE3 run --left r1.fq --right r2.fq -T RiceTE.fa -n HEG4 -o HEG4_out --te-aligner bwa
relocaTE3 align-genome -g reference.fa -f HEG4_out/flanking/*.flankingReads.fq -n HEG4 -o HEG4_out --genome-aligner bwa
```

`minimap2`, `bwa`, `bwa-mem2`, and `bowtie2` are pinned in the default pixi
environment; `blat` is in the separate `blat` environment (see Installation).
Non-minimap2 genome BAMs are named `{name}.repeat.{aligner}.sorted.bam` — with
the default `bwaaln` that is `{name}.repeat.bwaaln.sorted.bam`.

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
