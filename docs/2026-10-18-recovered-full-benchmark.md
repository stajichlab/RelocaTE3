# Full benchmark after RelocaTE2 recovery

Recorded: 2026-10-18 12:41 PDT, America/Los_Angeles.
RelocaTE3: `main`, `213a0df99af71a57fffad726138593b152fc9a77`, after merged PR #54.
Benchmark: `main`. Status: six baseline recoveries and aggregation verified.
No new caller changes, submissions, commits, or release tag in this review.

Recovery array **30041652** completed all six tasks with exit status zero.
All recovery `.complete` and `.promoted` markers exist, and the installed
chromosome outputs pass the new R2 validator. Finalizer **30041653** archived
the original outputs, promoted the recovered results, and completed full
aggregation on **2026-10-17 at 12:37:52 PDT**. Source/configuration guards pass.
Evidence comes from logs, markers, outputs, and resource exit records; this
review does not claim an independent scheduler-accounting check.

All **432 caller/sample accuracy reports** are now present: six callers across
9 mPing, 9 riceTElib, and 54 divergence samples. The comparisons below use
all 72 samples, superseding the provisional 66-sample comparisons in the
[October 16 report](2026-10-16-full-checkpoint-results.md).

## Reproduction and accuracy

Normalized-call multiset comparison, excluding only the caller-label field,
confirms **72/72 R2 results reproduce the historical successful R2 calls** and
**72/72 R3 BLAT/bwa-aln results reproduce the validated local-filter replay**.
Coordinates, family, TSD, strand, and genotype are included in this comparison.
Thus the large apparent gains caused by the incomplete R2 runs are resolved.

| Panel | Caller | TP | FP | FN | Precision | Recall |
|---|---|---:|---:|---:|---:|---:|
| mPing | R2 | 2,820 | 0 | 1,680 | 100% | 62.667% |
| mPing | R3 BLAT/bwa aln | 2,826 | 0 | 1,674 | 100% | 62.800% |
| riceTElib | R2 | 2,494 | 18 | 2,006 | 99.283% | 55.422% |
| riceTElib | R3 BLAT/bwa aln | 2,513 | 24 | 1,987 | 99.054% | 55.844% |
| Divergence | R2 | 5,888 | 77 | 21,112 | 98.709% | 21.807% |
| Divergence | R3 BLAT/bwa aln | 5,919 | 104 | 21,081 | 98.273% | 21.922% |

R3 has modestly higher pooled recall and F1 in each panel, with reduced
precision in the multi-TE panels. Net differences are **+56 TP, +33 FP,
and −56 FN** across the benchmark observations.

| Panel | Shared detections | R3 only | R2 only | Net R3 gain |
|---|---:|---:|---:|---:|
| mPing | 2,818 | 8 | 2 | +6 |
| riceTElib | 2,493 | 20 | 1 | +19 |
| Divergence | 5,873 | 46 | 15 | +31 |

The 74 R3-only and 18 R2-only observations are not independent loci; coverage
and divergence conditions reuse insertions, and canonical/zero-divergence
panels overlap. No statistical superiority or universal noninferiority claim
is made from these small net differences.

## Somatic and alternative-aligner results

Detected somatic observations increase from 1,191 to 1,194 in mPing, 1,033
to 1,039 in canonical riceTElib, and 2,135 to 2,145 in divergence. At 10%
cellular fraction (expected 5% allele fraction), counts are respectively
205 → 205, 189 → 188, and 345 → 341: the lowest-fraction detection gaps remain.
Correct somatic labels among detected somatic insertions are 249 → 250,
181 → 180, and 526 → 525. Insertion detection and somatic classification
remain distinct capabilities; the shared genotype limitation is not resolved.

Recall across all 72 samples:

| Caller / TE→genome alignment | mPing | Canonical | Divergence |
|---|---:|---:|---:|
| R2 | 62.667% | 55.422% | 21.807% |
| R3 BLAT→bwa aln | 62.800% | 55.844% | 21.922% |
| R3 BLAT→bwa MEM | 37.222% | 34.378% | 17.074% |
| R3 Bowtie2→bwa MEM | 29.778% | 30.178% | 13.815% |
| R3 bwa MEM→bwa MEM | 29.822% | 26.733% | 12.389% |
| R3 minimap2→minimap2 | 25.778% | 25.200% | 11.419% |

The BLAT/bwa-aln combination remains the best overall R3 configuration.
Other aligners do not provide comparable detection performance under these
settings. The prior limitations in Helitron and high-divergence detection
remain, consistent with exact reproduction of the validated call sets.

## Resource interpretation and release decision

Accuracy now covers 72 samples per caller, but valid full-pipeline R2 resource
measurements cover only **66**. The six invalid original resource rows were
archived, and recovery timings are kept separate. Combined resource-table
row counts are 52, 52, and 322 across the three panels (426 total), as intended.
Compare resources on the same 66 samples for both tools; do not count recovery
wall time as a full R2 run or interpret the missing rows as zero cost.

The October 16 complete-pair resource comparison therefore remains applicable:
at 30× multi-TE coverage, R3 BLAT/bwa-aln uses roughly 19–20 GiB peak RSS versus
6.5 GiB for R2. These are GNU time measurements from different nodes/dates,
not controlled speedup estimates or simultaneous aggregate job memory.

This is a stable, reproducible checkpoint with close accuracy parity and
modest net detection gains. It does not yet support the promised blanket
"same or better with no performance sacrifice": individual missed detections,
extra FPs, and high-depth memory costs remain.

Recommended next action: audit the **18 R2-only observations and the R3
false-positive differences**, prioritizing the low-fraction somatic losses.
Group repeated observations by locus and failure mechanism before considering
another small scientific change. Preserve the now-complete benchmark as the
baseline and avoid another full alignment sweep merely to inspect these cases.

## Evidence and review commands

From the benchmark root:

```bash
cat logs/checkpoint-213a0df-r2-recovery.submission
tail -n 5 logs/r2-recovery-finalize.30041653.log
python3.12 callers/relocate3-checkpoint-213a0df/verify_source.py \
  config/checkpoint-213a0df.source.json
```

Inputs reviewed: `reports/checkpoint-213a0df/datasets/*/{precision,correctness,resources}.tsv`,
per-sample R2/R3 `matches.tsv`, and normalized calls under
`runs/checkpoint-213a0df/`. Historical counterpart paths come from
`results/intact-fullreads-replay/2026-10-03-local/manifest.json` in RelocaTE3;
R3 comparison files are each task's `intact_candidate/calls.normalized.tsv`.
Paired matching reconciles 500 unique truth event IDs in both callers for
each sample. Counts are summed from integer numerators, not averaged percentages.

The old analysis script still assumes complete full-pipeline resource rows
and original resource logs for every task; it must be adapted before reuse
on the recovered cohort. No outputs were rewritten in this review. Aggregated
checkpoint reports were refreshed by the finalizer, but this review does not
claim the dashboard's active report-root selection was changed or verified.
