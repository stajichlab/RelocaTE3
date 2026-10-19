# Full benchmark: merged RelocaTE3 checkpoint

Follow-up: the six R2 recoveries and aggregation succeeded on October 17.
Use [the October 18 full-cohort assessment](2026-10-18-recovered-full-benchmark.md)
for the final 72-sample accuracy comparison. The provisional results below
remain preserved; the complete-pair resource comparison still applies.

Recorded: 2026-10-16, America/Los_Angeles, after 20:20 PDT.
Project: RelocaTE3, branch `main`, revision
`213a0df99af71a57fffad726138593b152fc9a77`, following merged PR #54.
Status: analysis complete; no caller changes, job submissions, commits, or release tag.

The merged BLAT/bwa-aln implementation exactly reproduces the previously
validated normalized calls in **all 72 samples**, including coordinates,
family, strand, TSD, and genotype (caller label excluded). This is successful
end-to-end confirmation of the earlier downstream replay. On the 66 samples
with complete chromosome outputs from both tools, RelocaTE3 has slightly
higher pooled detection recall and F1 in all three panels, but more false
positives in the two multi-TE panels. Individual R2-only detections and the
high-coverage memory gap remain. This does not establish universal
equal-or-better performance or justify an unconditional release claim.

## Completion and a RelocaTE2 baseline problem

Recovery array **29932898** produced all 360 R3 results. Original array
**29610466** produced 72 R2 reports. Aggregation **29932899** completed on
October 16 at 19:49:41 PDT. All 432 tasks have their expected completion log
line, run sentinel, report tables, and adapter `time -v` exit status zero.
The pinned source/configuration verifier passes. These are filesystem/log
checks; scheduler accounting was not needed to establish report completion.

Task success did not guarantee successful internal R2 workers. Six R2 runs
lack chromosome-specific non-reference output files:

| Panel | Sample | Array index | Missing chromosomes |
|---|---|---:|---|
| mPing | cov30x_rep2 | 7 | 5, 6, 9 |
| mPing | cov30x_rep3 | 8 | 1, 2, 3, 4, 5, 6, 9 |
| riceTElib | cov15x_rep1 | 57 | 1, 2, 3, 4, 5, 6, 7 |
| riceTElib | cov15x_rep3 | 59 | 9 |
| Diverged riceTElib | div010_rep02_cov5x | 138 | 4, 5, 6, 7 |
| Diverged riceTElib | div010_rep03_cov5x | 141 | 5, 6, 7 |

The two runs missing Chr1 report zero calls, losing 403 and 281 historical
true detections respectively. Their logs show concurrent `createdir(result)`
failures (`OSError: [Errno 17] File exists`) in R2 insertionFinder. The R2
source uses `if not exists: os.mkdir`, which races between chromosome workers.
The parent workflow continued, producing a header-only characterized table;
the adapter's nonempty-file check accepted that header as success.

The other 70 R2 normalized call sets equal their historical counterparts,
including the four runs missing only background chromosomes. Nevertheless,
all six incomplete samples are conservatively excluded from **every caller**
in the primary comparison: missing background outputs could conceal FPs.
That leaves **7 mPing, 7 canonical, and 52 divergence samples**. Original
reports are preserved. Their unfiltered totals should not be used as the
release comparison; they artificially inflate R3's apparent advantage.

## Primary comparison: complete paired samples

R3 means BLAT for TE reads and **bwa aln** for genome alignment. The profile's
`relocate3-blat-bwa` instead uses bwa MEM and is materially different.
Scoring uses one-to-one family-compatible matching to the truth TSD interval,
with 10 bp tolerance outside it. TSD sequence and genotype are assessed
separately. Counts are sample-by-event observations, not unique biological loci.

| Panel | Caller | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| mPing, 7 samples | R2 | 2,011 | 0 | 1,489 | 100% | 57.457% | 72.981% |
| mPing, 7 samples | R3 | 2,015 | 0 | 1,485 | 100% | 57.571% | 73.073% |
| riceTElib, 7 samples | R2 | 1,906 | 14 | 1,594 | 99.271% | 54.457% | 70.332% |
| riceTElib, 7 samples | R3 | 1,919 | 19 | 1,581 | 99.020% | 54.829% | 70.577% |
| Divergence, 52 samples | R2 | 5,877 | 75 | 20,123 | 98.740% | 22.604% | 36.786% |
| Divergence, 52 samples | R3 | 5,908 | 103 | 20,092 | 98.286% | 22.723% | 36.912% |

Paired event identities show why net gains are not universal noninferiority:

| Panel | Found by both | R3 only | R2 only | Net R3 TP gain |
|---|---:|---:|---:|---:|
| mPing | 2,009 | 6 | 2 | +4 |
| riceTElib | 1,905 | 14 | 1 | +13 |
| Divergence | 5,863 | 45 | 14 | +31 |

There are 65 R3-only and 17 R2-only observations. Some repeat the same
insertion across conditions, and canonical riceTElib overlaps the 0%
divergence panel. These are descriptive comparisons, not independent evidence
for statistical superiority. Exclusions also unbalance coverage, so panel
recalls cannot be directly compared to earlier 9/9/54-sample pooled recalls.

For completeness, the **full 72-sample R3** totals are unchanged from replay:
mPing 2,826 TP / 0 FP / 1,674 FN; riceTElib 2,513 / 24 / 1,987;
divergence 5,919 / 104 / 21,081. Historical successful R2 totals were
2,820 / 0 / 1,680; 2,494 / 18 / 2,006; and 5,888 / 77 / 21,112.
Those historical figures are context, not substituted fresh-run results.

## Coverage, somatic insertions, and genotype

Detection recall on complete paired samples, R2 → R3:

| Panel | 5× | 15× | 30× |
|---|---:|---:|---:|
| mPing | 40.667 → 41.067% | 65.933 → 65.800% | 82.400 → 82.400% |
| riceTElib | 36.533 → 36.667% | 58.000 → 58.400% | 71.200 → 71.800% |
| Divergence | 13.713 → 13.737% | 22.889 → 23.011% | 30.222 → 30.422% |

Only one complete mPing 30× replicate and one canonical 15× replicate remain.
The divergence 5× cohort has 16 samples; other divergence depths have 18.

| Panel | Homozygous recall R2 → R3 | Heterozygous recall R2 → R3 | Somatic recall R2 → R3 |
|---|---:|---:|---:|
| mPing | 95.000 → 95.000% | 82.857 → 83.286% | 36.476 → 36.524% |
| riceTElib | 84.143 → 85.143% | 75.714 → 76.000% | 37.476 → 37.667% |
| Divergence | 38.519 → 38.750% | 33.462 → 33.635% | 13.679 → 13.744% |

Somatic detection has small net gains (+1, +4, +10 observations), but the
lowest cellular fraction still exposes losses:

| Panel | 10% cells: TP R2 → R3 | 20% cells | 40% cells |
|---|---:|---:|---:|
| mPing | 113 → 113 | 255 → 254 | 398 → 400 |
| riceTElib | 155 → 154 | 241 → 242 | 391 → 395 |
| Divergence | 345 → 341 | 639 → 645 | 1,150 → 1,158 |

These simulated insertions occupy one haplotype in the affected cells;
10%, 20%, and 40% cellular fractions correspond to expected allele fractions
of 5%, 10%, and 20%. Detecting the insertion and labeling it somatic are
different tasks. Among detected somatic insertions, the correct somatic label
occurs in only 12.92 → 12.78%, 21.35 → 21.11%, and 24.65 → 24.49% for the
three paired panels. Most are labeled heterozygous by both tools. This is a
shared characterization limitation, not a new broad R3 regression.

On exactly the same jointly detected events, correct genotype counts are
1,139 → 1,139; 1,104 → 1,102; and 3,674 → 3,673. Exact TSD counts on those
events are 1,969 → 1,967; 1,847 → 1,853; and 5,077 → 5,094. Genotype
accuracy is close and TSD correctness is mixed/slightly improved overall.
The benchmark adapters use different all-read characterization alignments
(R2 bwa MEM, R3 minimap2), so genotype/runtime differences are workflow-level
observations, not isolated measurements of Python implementation changes.

## TE groups and divergence

TP differences on complete paired samples:

| TE group | Canonical R3 minus R2 | Divergence R3 minus R2 |
|---|---:|---:|
| CACTA | +1 | +2 |
| Helitron | +4 | +9 |
| LINE | 0 | −2 |
| LTR Copia | 0 | 0 |
| LTR Gypsy | +2 | +7 |
| MULE | +1 | +1 |
| PIF/Harbinger | 0 | +3 |
| SINE | −1 | 0 |
| Tc1/Mariner | +5 | +8 |
| hAT | +1 | +3 |

Helitron remains poorly detected by both tools: canonical recall 1.43% versus
2.57%, diverged 3.00% versus 3.35%. These data do not support broad sensitivity
claims across TE types. The R2-only SINE event is Os3874/TE000386 at 10% cells
and 30×, repeated in canonical and zero-divergence data. Several LINE losses
involve Os3328/TE000172 across divergence conditions. All lost/gained event
identities are retained in `changed_events.tsv`.

| Divergence | Paired truth observations | R2 recall | R3 recall |
|---|---:|---:|---:|
| 0% | 4,500 | 55.422% | 55.844% |
| 2% | 4,500 | 46.111% | 46.200% |
| 5% | 4,500 | 23.289% | 23.422% |
| 10% | 3,500 | 6.257% | 6.343% |
| 15% | 4,500 | 0.689% | 0.667% |
| 20% | 4,500 | 0.222% | 0.222% |

R3 closely follows R2's severe sensitivity decline with divergence. It does
not solve high-divergence detection. At 10% and 20%, BLAT/bwa MEM has slightly
more TPs than BLAT/bwa aln (227 versus 222; 11 versus 10), but this does not
offset its substantially worse overall precision and recall.

## Other aligner combinations

Recall on the same complete paired samples:

| Caller / TE→genome aligners | mPing | Canonical | Divergence |
|---|---:|---:|---:|
| R2 | 57.457% | 54.457% | 22.604% |
| R3 BLAT→bwa aln | 57.571% | 54.829% | 22.723% |
| R3 BLAT→bwa MEM | 35.257% | 33.486% | 17.688% |
| R3 Bowtie2→bwa MEM | 28.514% | 29.514% | 14.327% |
| R3 bwa MEM→bwa MEM | 28.486% | 26.029% | 12.838% |
| R3 minimap2→minimap2 | 24.657% | 24.486% | 11.831% |

All configurations have zero mPing FPs. Canonical/divergence FP counts are
19/103 for BLAT/aln, 69/421 for BLAT/MEM, 29/124 for Bowtie2/MEM, 21/74 for
MEM/MEM, and 12/52 for minimap2. Lower FP counts for the latter methods come
with far fewer true detections; they are not equivalent replacements for
the BLAT/aln workflow. BLAT/aln has the best pooled F1 in every panel.

## Runtime and memory

Median adapter wall time and GNU `time -v` peak RSS for complete paired
samples. R2 and R3 ran on different dates/nodes; these are observed workflow
costs, not controlled speedup estimates. Peak RSS is the metric reported by
GNU time, not independently measured aggregate simultaneous job memory.

| Panel | Coverage | Samples | Hours R2 → R3 | GiB R2 → R3 |
|---|---:|---:|---:|---:|
| mPing | 5× | 3 | 0.58 → 0.14 | 6.48 → 5.27 |
| mPing | 15× | 3 | 1.19 → 0.52 | 6.48 → 6.72 |
| mPing | 30× | 1 | 2.36 → 0.93 | 6.49 → 6.72 |
| Canonical | 5× | 3 | 2.31 → 1.95 | 6.48 → 5.27 |
| Canonical | 15× | 1 | 5.96 → 5.90 | 6.48 → 9.85 |
| Canonical | 30× | 3 | 11.71 → 12.53 | 6.49 → 19.21 |
| Divergence | 5× | 16 | 2.20 → 2.11 | 6.48 → 5.39 |
| Divergence | 15× | 18 | 6.16 → 5.98 | 6.48 → 10.14 |
| Divergence | 30× | 18 | 12.06 → 12.28 | 6.49 → 19.74 |

R3 is faster in these mPing runs. Multi-TE BLAT/aln runtimes are broadly
similar, with higher observed R3 times at 30×. Median within-sample R3/R2
time ratios at 30× are 1.054 for canonical and 1.064 for divergence (ratios
of medians are a different statistic). Memory is approximately threefold
higher at 30× multi-TE depth, despite the earlier improvement from roughly
34 GiB to 19 GiB. Maximum observed R3 BLAT/aln RSS is about 19.77 GiB.

## Release assessment and next action

This establishes reproducibility of the merged R3 scientific changes and
close accuracy parity for the BLAT/bwa-aln configuration. It does not establish
that users sacrifice nothing: residual FPs, individual missed insertions,
low-fraction somatic losses, and higher high-depth memory remain. The
simulation also cannot certify utility on all real samples or reference-TE,
excision, and other broader feature behavior.

The immediate next action is to repair the **benchmark R2 execution checks**
and recover the six incomplete R2 runs before changing R3's scientific rules
or pinning a release. Pre-create the shared results directory before parallel
chromosome workers and detect failed/missing expected chromosome outputs;
do not require nonzero call counts, because genuine zero-call runs are valid.
Prefer recovering insertion calling/characterization from retained alignments
where validated. Preserve all existing output and rerun aggregation only
after recovery succeeds. None of that repair or resubmission was performed
as part of this analysis request.

## Reproduce this audit

Run from the RelocaTE3 root. This reads small tables, log tails, and file
metadata; it does not rerun alignments. Use a new output directory to repeat:

```bash
python3.12 scripts/audit_full_checkpoint.py \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --replay results/intact-fullreads-replay/2026-10-03-local \
  --stabilization results/family-replay/2026-10-02-stabilization \
  --exclude-pair mping/cov30x_rep2 \
  --exclude-pair mping/cov30x_rep3 \
  --exclude-pair ricetelib/cov15x_rep1 \
  --exclude-pair ricetelib/cov15x_rep3 \
  --exclude-pair ricetelib_divergence/div010_rep02_cov5x \
  --exclude-pair ricetelib_divergence/div010_rep03_cov5x \
  --r2-array 29610466 --r3-array 29932898 \
  --output results/full-checkpoint-analysis/2026-10-16-complete-chromosomes
```

The command above has already completed. Its output contains `summary.tsv`,
`samples.tsv`, `strata.tsv`, `paired.tsv`, `paired_strata.tsv`,
`changed_events.tsv`, `status_confusion.tsv`, `resources.tsv`,
`historical_reproduction.tsv`, `task_audit.tsv`, and `provenance.json`.
Provenance records 946 input hashes and the analysis-script hash; it does not
claim to hash all FASTQs, BAMs, or external executables. The script checks
sample uniqueness/completeness, 500 truth events per sample, TP reconciliation,
coordinate policy, event identities, adapter exits, and log completions.
Ruff formatting/checks pass. No new computationally intensive tests were run.

The earlier `results/full-checkpoint-analysis/2026-10-16/` directory is a
preserved preliminary audit excluding only the two Chr1 failures; use the
`-complete-chromosomes` directory for the primary assessment in this document.
The automated benchmark reports and dashboard inputs were not rewritten or
rewired. Any display of unfiltered checkpoint results retains the R2 failures
until recovery and aggregation are completed.
