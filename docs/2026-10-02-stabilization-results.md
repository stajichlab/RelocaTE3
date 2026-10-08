# Stabilization replay: sensitivity restored, precision gap remains

Recorded: 2026-10-02 15:10 PDT, America/Los_Angeles.
Project: RelocaTE3, `main`, base `942a5ea`, following merged PR #49.
Status: array 29349059 completed all 72 tasks; tested pairing restored to
working source. No commit, push, merge, release pin, or dashboard publication.

## What this run establishes

All 72 task completion markers, comparison files, and successful SLURM logs
are present. Frozen source hashes verify. With the historical breakpoint
function restored, original-evidence calls exactly reproduce historical
normalized calls in **72/72 samples**. The previously observed 30 lost truth
detections disappear. The restored code also brings back 11 false positives
that the rejected pairing change had removed; this is not a free precision fix.

The corrected family-evidence implementation gains three truth detections
and removes three FPs relative to that reproduced baseline, with **no lost
truth detections**. Gains are the repeated Os3328 LINE event TE000172 at
5%, 10%, and 15% divergence in replicate 3 at 30x. They are three benchmark
observations of one locus, not three independent biological loci.
The prior somatic hAT loss is absent after preserving the unresolved fallback.
Coordinates, TSD, strand, and genotype remain identical between evidence arms
in all 72 samples.

Production `_pair_breakpoints` now matches the successfully tested snapshot.
All production Python files byte-match their frozen source hashes, including
the family fallback fix. Existing memory work and other changes are preserved.
The local duplicate tests retain the rejected coordinate experiment explicitly
as an experiment; it is no longer the expected production behavior. Those
duplicate loci remain a known limitation requiring a sensitivity-preserving fix.

**106 focused tests pass**, including breakpoint pairing, real duplicate
fixtures, family resolution/determinism/metadata, real LINE evidence, trimming,
insertion/TSD calling, characterization, and memory regressions. This is not
the full external-tool suite or a fresh end-to-end alignment benchmark.

## R2 comparison: same inputs and scoring policy

Both callers use BLAT/bwa aln evidence. R2 calls are archived; R3 downstream
calling and characterization are replayed from stored alignments. Matching
uses the previously documented `tsd-interval` policy with a 10-bp tolerance.
Genotype/status and TSD correctness remain separate metrics.

| Panel | Caller | TP | FP | FN | Precision | Recall |
|---|---|---:|---:|---:|---:|---:|
| mPing | R2 | 2820 | 0 | 1680 | 100% | 62.667% |
| mPing | Stabilized R3 | 2826 | 0 | 1674 | 100% | 62.800% |
| riceTElib | R2 | 2494 | 18 | 2006 | 99.283% | 55.422% |
| riceTElib | Stabilized R3 | 2513 | 27 | 1987 | 98.937% | 55.844% |
| Diverged riceTElib | R2 | 5888 | 77 | 21112 | 98.709% | 21.807% |
| Diverged riceTElib | Stabilized R3 | 5919 | 129 | 21081 | 97.867% | 21.922% |

Relative to R2, paired truth observations gained/lost are 8/2 for mPing,
20/1 for canonical riceTElib, and 46/15 for diverged riceTElib. Pooled recall
is higher in all three panels, but this does not establish universal parity.
R3 has fewer TPs in 2/9 mPing samples and 6/54 divergence samples. It has
more FPs in 5/9 canonical samples and 30/54 divergence samples.

## Somatic insertions and difficult strata

| Panel | Somatic truth | R2 detected | R3 detected | R2 correctly labeled somatic | R3 correctly labeled somatic |
|---|---:|---:|---:|---:|---:|
| mPing | 2700 | 1191 | 1194 | 249 | 250 |
| riceTElib | 2700 | 1033 | 1039 | 181 | 180 |
| Diverged riceTElib | 16200 | 2135 | 2145 | 526 | 525 |

The overall somatic detection gain hides a low-frequency shortfall. At 10%
cellular fraction, canonical R3 detects 188 versus R2's 189; diverged R3
detects 341 versus R2's 345. Detection improves at 20% and 40% fractions.
Correctly labeling a detected insertion as somatic remains weak in both tools.

Canonical SINE detections remain 279 versus R2's 280. Diverged LINE detections
improve with the family fix but remain 670 versus 673. At 15x mPing, R3
detects 987 versus 989. At 15% divergence, R3 detects 30 versus 31 and emits
14 versus 4 FPs. At 20% divergence both detect only 10/4500, while R3 emits
6 versus 2 FPs. Similar performance at that divergence is not strong absolute
sensitivity.

The 0%-divergence panel overlaps canonical riceTElib; coverage and replicate
conditions also reuse loci. Treat these as paired benchmark observations,
not independent biological discoveries. FP counts are not assigned to truth
genotype/fraction strata because unmatched calls have no truth class.

## Reproducibility and output

Result directory: `results/family-replay-analysis/2026-10-02-stabilization/`.
It contains `summary.tsv`, `samples.tsv`, `changed_events.tsv`, `paired.tsv`,
and provenance hashes. Every coverage, divergence, TE-group, biological-class,
and class/fraction stratum reconciles to panel totals.

Command already executed; use a new output directory to rerun:

```bash
.pixi/envs/default/bin/python scripts/summarize_family_replay.py \
  --replay results/family-replay/2026-10-02-stabilization \
  --historical results/coordinate-rescore/2026-09-28 \
  --output results/family-replay-analysis/2026-10-02-stabilization
```

This lightweight analysis reads result tables. No alignment scans, pipeline
reruns, benchmark cleanup, or additional compute submission were needed.

## Decision and next action

Keep the family fix and restored breakpoint behavior. This candidate improves
on historical R3 detection without a new observed detection regression, but
**do not pin a release or promise no sacrifice versus RelocaTE2 yet**.
The remaining precision gap is nine extra FPs in canonical riceTElib and
52 extra FPs in the divergence panel, alongside the paired detection losses
and low-fraction somatic shortfall described above.

Next, classify those remaining R3-only false positives against the R2 source
and read evidence, beginning with a repeated causal category. Any proposed
correction must retain this stabilized sensitivity and include event-level
gain/loss checks. Avoid another broad threshold or coordinate shift justified
only by one or two loci. Full-pipeline memory/runtime, complete tests, and
the R2 utility/feature audit remain separate release gates; this replay does
not certify them. The future Rust/Nextflow/long-read/pangenome roadmap remains
deferred. No additional benchmark submission is needed from the maintainer
for this completed result.
