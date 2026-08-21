# BLAT TE search: RelocaTE2 parameters, chunking, and seqtk

2026-08-16

## Why this matters

RelocaTE3 defaults to BLAT for the TE-library search because RelocaTE2 does.
But RelocaTE2 runs BLAT with **sensitivity settings**, not at BLAT's defaults:

    blat -minScore=10 -tileSize=7 <te_fasta> <reads.fa> <out.psl>   # relocaTE2.py:545

BLAT's own defaults are `minScore=30`, `tileSize=11`. RelocaTE3 ran at those
defaults, and it cost real calls. Measured on the mPing benchmark
(`cov30x_rep1`), over the 19 junction reads RelocaTE2 used at sites RelocaTE3
missed outright:

| BLAT settings | reads found, of 19 |
|---|---|
| defaults (minScore=30, tileSize=11) | **0** |
| `-minScore=10 -tileSize=7` (RelocaTE2) | **19** |

End-to-end effect on mPing `cov5x_rep1`:

| | TE reads at trim | calls | TP | FP | precision | recall | F1 |
|---|---|---|---|---|---|---|---|
| before | 3895 | 183 | 183 | 0 | 1.000 | 0.366 | 0.536 |
| **after** | **4128** | **201** | **201** | **0** | **1.000** | **0.402** | **0.573** |
| RelocaTE2 | — | 202 | 202 | 0 | 1.000 | 0.404 | 0.575 |

## The interaction you must not forget

These parameters were added once and reverted (`02e68eb`) because they collapsed
precision on riceTElib. That was real, but the cause was **not** the parameters:
`characterize.py` had ported RelocaTE2's gate (`characterizer.pl:91`) with `or`
where RelocaTE2 has `and`, so every one-sided call the extra sensitivity
produced sailed into the output. With that gate corrected, one-sided calls are
dropped the way RelocaTE2 drops them.

**The BLAT sensitivity and the characterize gate must always be measured
together.** Changing either alone will look like a regression in the other.

## Runtime: why chunking and seqtk are required

The sensitised parameters make BLAT much slower per base, and BLAT is
single-threaded. RelocaTE3 originally ran one BLAT process over the entire read
library, which is not viable: a single mPing 5x *side* ran past 13 minutes,
while RelocaTE2 completes the whole sample in about that.

RelocaTE2 solves this in two ways, and RelocaTE3 now does the same:

1. **Chunking** — split the query into 200k-read chunks and search them
   concurrently (`relocaTE2.py:91,110`, `fastq_split.pl`).
   `BlatBackend.BLAT_CHUNK_SEQS` + a thread pool in `_blat_side`.
2. **seqtk for FASTQ→FASTA** — BLAT cannot read FASTQ. Doing the conversion in
   Python dominated the stage: roughly 4M reads per 35 minutes, meaning hours on
   one 30x side before BLAT even started. `seqtk seq -A -l 0` does it in C and
   emits exactly two lines per record, so plain `split -l` cuts on record
   boundaries.

`seqtk` is now a pixi dependency. A pure-Python fallback remains for
environments without it, so nothing breaks — it is just slow.

Sequences for the PSL→SAM conversion are pulled back with `seqtk subseq` for
**only the reads BLAT matched** (a few thousand), instead of holding the whole
library in memory.

## Operational note

Run the benchmark through SLURM (`pipeline/submit_benchmark.sh`), never on a
login node. The BLAT stage is genuinely compute-heavy: 8 concurrent BLAT
processes per sample for tens of minutes at 5x, longer at 30x. Give the blat
variants `--cpus-per-task` at least equal to the chunk concurrency you want.
