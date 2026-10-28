# BLAT path memory: per-chunk trim (2026-10-27)

Branch `perf/blat-chunked-trim` (from `release/v0.2.0` f7c93f0). Not pushed.

## Purpose

Bring RelocaTE3's peak memory on the BLAT TE-search path down to RelocaTE2's
per-chunk bound (19 GiB at riceTElib 30x before) without changing any output.

## What changed

| Commit | Peak | Change |
|---|---|---|
| f2cf776 | 1. trim dict over all TE-hit reads with seq/qual (~16 GiB) | `RelocaTE.write_trimmed_reads_chunked`: per chunk, unchanged best-hit selection (`_parse_te_bam`) + FASTQ restore from that chunk's slice + unchanged trim branches (`_write_read`), spilled to disk; spills merged on each read's first-hit coordinate key |
| fb7f2fd | 2. merged PSL + whole-side matched-sequence dict + in-process `pysam.sort` | `BlatBackend._blat_chunk`: BLAT -> this chunk's sequences -> SAM -> `samtools sort -m 256M` -> chunk BAM; side BAM = `samtools merge`. `identify_TE_reads` uses `map_te_library_chunked` + chunked trim, then deletes `<outdir>/<name>.<side>.te_chunks/` |
| e4384fb | 3. step-4 set of all `te_hit_names` (~34M names) | `read_te_hit_names(path, bases=flanking_read_bases(flanks))` keeps only names of flanking pairs; same answers for every planner query |

Unchanged on purpose: `map_te_library`, `relocaTE3 trim` (and so the
benchmark replay tier), and the minimap2/bwa/bowtie2 TE paths still use the
whole-BAM trim.

## Decisions / logic

- Chunk-local best hit = global best hit: a read's hits come only from the
  BLAT chunk holding its FASTQ record; mates are separate sides, and names are
  canonicalized (`/1`, `/2`) per side when the chunk SAM is written, exactly
  as `canonicalize_te_bams` did to the merged BAM.
- Output order is kept, not changed to RelocaTE2's FASTQ order: `bwa aln`
  `samse/sampe` pick among equal genome hits with a sequential RNG, so
  flanking-read order changes placements. Order = first TE hit in the
  coordinate-sorted side BAM; `samtools sort` and `merge` are both stable
  (checked: ties keep input order, then file order), so merging chunk spills
  on `(tid, pos, is_reverse)` with chunk order on ties reproduces it.
- The chunked trim checks each chunk consumes exactly its `n_reads` FASTQ
  records and that the FASTQ ends with the last chunk (else `ValueError`).
- Assumption: read names are unique within a side (as in RelocaTE2, which
  also trims per chunk). A duplicate name split across chunks would now be
  trimmed twice instead of merged.
- Disk: the trim spill (~ size of one side's ContainingReads + flanking, about
  6-7 GB per side at riceTElib 30x) sits in `te_containing/.<name>.<side>.trim_chunks.*`
  and is deleted per side; chunk BAMs (~ TE BAM size) live until trim ends.

## Commands

```bash
PIXI_BIN=../RelocaTE3/results/release-validation/v0.2.0rc1/clone/.pixi/envs/blat/bin
PATH=$PIXI_BIN:$PATH PYTHONPATH=$PWD/src $PIXI_BIN/python -m pytest tests -p no:cacheprovider --junitxml=<scratch>/junit.xml

# real-size memory, old vs new (from the repo root)
sbatch validation/blat_memory/run_step3_memcheck.sh ../RelocaTE3-release results/blat_memory/old_release_f7c93f0
sbatch validation/blat_memory/run_step3_memcheck.sh . results/blat_memory/new_fb7f2fd
sbatch validation/blat_memory/run_step4_memcheck.sh ../RelocaTE3-release results/blat_memory/step4_old_release_f7c93f0
sbatch validation/blat_memory/run_step4_memcheck.sh . results/blat_memory/step4_new_e4384fb
```

## Parity evidence

- Tests: 319 baseline -> 327, 0 failures, 0 skipped (blat env). New:
  `tests/blat_chunked_trim_test.py` (synthetic paired 9-chunk fixture: chunked
  trim byte-identical to whole-file trim; one chunk parsed at a time; FASTQ /
  chunk mismatch rejected; real BLAT 4 chunks vs 1 chunk: same BAM records and
  trim outputs), `tests/te_hit_names_filter_test.py`. Rewrote two tests that
  pinned the old merged-PSL internals (`aligners_test` chunking test,
  `memory_optimization_test` sequences test).
- End to end: `run-all` on the Chr3 2 Mb acceptance fixture (6x, RiceTE.fa,
  blat + bwaaln, `--genotype --repeatmasker`), release f7c93f0 vs this branch
  forced to 11 BLAT chunks/side and at the default 1 chunk/side: all 28 output
  files identical (text byte-for-byte, BAMs record-for-record; only `@PG`
  header lines differ). 197-line characTErized GFF in both.

## Memory (riceTElib cov30x_rep1, 8 threads)

Peak RSS two ways: `/usr/bin/time -v` (largest single process) and a 5 s
sampler summing RSS over the whole process tree (python + blat/samtools/bwa
children). Slurm MaxRSS is not used (it counts page cache on this cluster).

| Step | Code | Largest process | Process tree | Wall | Outputs |
|---|---|---|---|---|---|
| 4 `align-genome` (cached step-3 inputs, checkpoint-213a0df) | release f7c93f0 (job 31134363) | 6.40 GB | 8.89 GB | 15:08 | |
| 4 `align-genome` | this branch e4384fb (job 31134364) | 2.48 GB | 2.98 GB | 11:59 | identical (2/2 BAMs, 586,061 genome records) |
| 3 `run` (BLAT + trim) | release f7c93f0 (job 31133962) | pending | pending | | |
| 3 `run` | this branch fb7f2fd (job 31133963) | pending | pending | | |

Step 3 is BLAT-bound (~190 chunks/side, 8 at a time); expect 8+ h per job.
When done: `cat results/blat_memory/{old_release_f7c93f0,new_fb7f2fd}/{time.txt,tree_rss_peak.txt}`
and `bash validation/blat_memory/compare_outputs.sh results/blat_memory/old_release_f7c93f0 results/blat_memory/new_fb7f2fd`
(sbatch it: the TE BAMs and ContainingReads are several GB).

## Failures / issues

- None in tests. Lint: 6 pre-existing ruff findings in `reference_te.py`
  (untouched).

## Next steps

- Read the step-3 job results (31133962 old, 31133963 new); diff their outputs.
- Full benchmark tier (the dev replay tier refuses this change because
  `aligners.py` differs from 213a0df).
