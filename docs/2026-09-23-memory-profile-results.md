# Completed memory profile: baseline reproduced, Python trimming peak localized

Recorded: 2026-09-23, 14:05 PDT (America/Los_Angeles).
Project: RelocaTE3, branch `main`, merged PR #49.
Algorithm revision: `942a5ea266bd000a47e27df7ef9e2536c877acbd`.
Job: `29005979`, node `r26`, eight CPUs, 64 GB requested.
Dataset/caller/sample: riceTElib / BLAT-bwa aln / `cov30x_rep1`.
Status: profiling and baseline comparisons completed successfully; analysis
complete. No algorithm changes or additional compute jobs in this review.

## Completion and preservation of results

The job ran September 22, 11:58:33–23:39:54 PDT. The measured adapter wall
time was 11:41:05, with exit status zero, no swaps, and maximum RSS
36,121,384 KiB = **34.45 GiB**. The historical sample was 33.79 GiB and
11.19 hours; this replay reproduces the high-memory behavior, not an
optimization. Timing differences do not isolate profiling overhead from
hardware/load/cache effects.

The `.profile_complete` sentinel exists. `comparison.json` records five
successful comparisons. An independent reread during this review confirmed
row-multiset equality for normalized calls and all four scored tables:
`matches.tsv`, `precision.tsv`, `correctness.tsv`, and
`false_positive_calls.tsv`. This is normalized/scored equality, not a claim
that every raw output file is byte-identical.

Both runs contain 368 calls: 302 true positives, 66 false positives, and
198 false negatives against 500 truth events. Thus neither sensitivity nor
false-positive behavior changed during the profiling replay.

## Memory localization

| Process/phase | Observed peak RSS (GiB) | Interpretation |
|---|---:|---|
| RelocaTE3 Python map + trim | 34.45 | Dominant peak, during post-alignment trimming |
| Same Python process, left BLAT sequence recovery/conversion interval | 20.68 | Earlier, separate memory-pressure interval |
| Same Python process, right BLAT sequence recovery/conversion interval | 20.80 | Earlier, separate memory-pressure interval |
| External samtools sort for all-read BAM | 6.63 | Secondary cost, not the 34-GiB peak |
| Python align-genome process | 6.17 | Secondary cost; includes in-process pysam operations |
| External minimap2 all-read alignment | 3.62 | Not the dominant cost |
| Python find-insertions process | 0.81 | Not the dominant cost |

These are sampled RSS values, not a sum of independent peaks and not
Python-only allocations. Native libraries inside Python share its RSS.
The sampled peak for the map/trim process agrees with GNU time's recorded
high-water value. The sample interval was two seconds; short peaks in other
processes may be missed.

The dominant process was PID 1493090 (`relocaTE3 run`). Its peak occurred
38,962.647 seconds after monitoring began. The final TE BAM timestamps were
22:20:22 (left) and 22:30:19 (right); the 34.45-GiB peak occurs later, during
the trimming portion, before `align-genome` starts at 22:53:07. The trim
stage reports 2,303,152 flanking reads. Neither BLAT nor BWA executable
memory explains the dominant peak.

For the earlier intervals, the parent's peaks were at elapsed 18,910.865
seconds (left) and 36,374.334 seconds (right), immediately before external
SAM-to-BAM conversion starts at approximately 18,923.5 and 36,387.1 seconds.
This localizes additional pressure to Python's BLAT-result sequence recovery
and SAM-writing path, but does not distinguish individual allocations there.

## First optimization candidate: release completed mate records

In `src/RelocaTE3/librelocate.py`, `write_trimmed_reads` loops over the two
mate BAMs. At line 186 it assigns:

```python
coord = self._parse_te_bam(...)
```

At the end of the first iteration, `coord` still references the first mate's
dictionary even though that mate's output has been written. Python evaluates
the next call to `_parse_te_bam` before replacing `coord`, so the first
dictionary remains live while the second dictionary is built. The lifetime
overlap follows directly from the code; its precise contribution to peak RSS
has not yet been measured by an allocation trace or a before/after experiment.
It is consistent with the observed roughly 18-GiB to 34-GiB increase during
trimming.

Recommended first change: explicitly release the completed mate dictionary
after `_write_direction` returns, before the next parsing iteration, or use
a per-mate helper whose locals end at that boundary. Preserve read-selection,
sequence/quality restoration, iteration order, output content, and thresholds.
Add a regression test that the previous mate's record container is no longer
retained when the next one is parsed, as well as existing output tests.

This is a memory-lifetime change, not stricter evidence filtering. It should
not discard any biological evidence. The expected memory benefit must still
be measured; no numeric reduction is guaranteed yet.

## Second target: buffered sequence recovery

`src/RelocaTE3/aligners.py`, `_query_sequences`, uses
`subprocess.run(..., capture_output=True, text=True)`, then
`proc.stdout.splitlines()`, while constructing a sequence dictionary. This
keeps a complete output string, a line list, and dictionary data alive
together. The no-seqtk fallback similarly uses `fh.read().splitlines()`.
Streaming the FASTA records rather than materializing redundant buffers is
the next candidate; preserve headers, sequences, subprocess error handling,
and both the seqtk and fallback paths.

The run's provenance confirms that seqtk was available, so the fallback
cannot be blamed for this particular replay. The measured 20.80-GiB earlier
peak also means fixing the later dictionary overlap alone will not make
the entire workflow fit in RelocaTE2's roughly 6.5-GiB footprint.

## Efficient verification plan

The replay retained both TE-alignment BAMs and the original FASTQ paths.
For the first lifetime change, rerun only trimming from those BAMs into a
new directory, supplying the original FASTQs so quality restoration remains
identical. The existing `trim` CLI supports this via `--bam` and `--fastq`.
Submit this substantial work through SLURM; do not parse the multi-GB BAMs
on a login node.

Compare every regenerated trimming artifact with the frozen baseline and
measure peak RSS. Reusing the TE BAMs avoids repeating approximately ten
hours of TE searching/conversion. Before accepting the optimization, replay
the downstream stages and compare normalized calls and scoring as well.
If needed, a baseline trim-only run provides an apples-to-apples memory
comparison without retained allocations from the preceding alignment phase.

Only after validating the first change should we address the earlier
sequence-buffering peak. Keep false-positive/false-negative algorithm
changes separate, so memory savings cannot be confused with dropped read
evidence.

## Evidence and lightweight review commands

All profile paths are relative to the RelocaTE3 root:

```text
logs/relocate3-memory.29005979.log
results/memory-profile/2026-09-21-ricetelib-cov30x-rep1/
  .profile_complete
  comparison.json
  adapter.time-v.txt
  tools.txt
  profile/summary.json
  profile/process_memory.tsv
  profile/adapter.log
  baseline/
  replay/
  score/
```

```bash
cat logs/relocate3-memory.29005979.log
cat results/memory-profile/2026-09-21-ricetelib-cov30x-rep1/comparison.json
cat results/memory-profile/2026-09-21-ricetelib-cov30x-rep1/adapter.time-v.txt
```

Additional lightweight checks streamed the 102-MB process-memory table,
selected PID 1493090, and identified its maximum sampled RSS and the earlier
conversion peaks. Only report/log reads and source inspection were needed;
no alignments, BAM scans, or new full tests were run.

Next action: implement and test the per-mate dictionary-lifetime correction,
then measure a trim-only replay before changing other memory paths or filters.
