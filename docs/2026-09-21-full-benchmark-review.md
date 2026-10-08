# Full benchmark review and dashboard audit

Recorded: 2026-09-21, approximately 16:45 America/Los_Angeles.
Project: RelocaTE3, branch `main`; parity implementation merged in PR #49.
Benchmark repository: `../../relocate_benchmark/relocate-benchmark`, branch `main`.
Purpose: assess the six-caller full benchmark and verify dashboard presentation.

## Completion and scope

Array `28028811`: all 432 expected task logs end with the correct dataset,
caller, and sample completion marker; all 432 `.run_complete` files exist.
Aggregation `28028812` completed on September 7 at 04:48 PDT.
Configuration: `config/benchmark.full-aligners.toml` in the benchmark repository.
Reports: `reports/datasets/{mping,ricetelib,ricetelib_divergence}/`.

There are six callers per dataset, with 9 mPing samples, 9 riceTElib samples,
and 54 divergence samples (72 samples per caller). The dashboard's schema
validator accepts all reports. Precision and resources have 54, 54, and 324
rows respectively; precision has no duplicate caller/sample keys.
No new caller jobs or scoring runs were submitted during this review.

## Detection results

These values pool counts across samples within each dataset. Precision is
sum(TP)/sum(calls), recall is sum(detected truth)/sum(truth), and F1 is their
harmonic mean. Each caller has 4,500 truth-event opportunities in each ordinary
panel and 27,000 in the divergence panel. These include repeated simulated
events across coverage/replicate conditions, not that many independent loci.

| Dataset | Caller (TE search / genome alignment) | TP | FP | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|
| mPing | RelocaTE2 | 2819 | 1 | 0.9996 | 0.6264 | 0.7702 |
| mPing | R3 BLAT / bwa aln | 2826 | 0 | 1.0000 | 0.6280 | 0.7715 |
| mPing | R3 BLAT / bwa-MEM | 2785 | 0 | 1.0000 | 0.6189 | 0.7646 |
| mPing | R3 Bowtie2 / bwa-MEM | 2647 | 0 | 1.0000 | 0.5882 | 0.7407 |
| mPing | R3 bwa-MEM / bwa-MEM | 2643 | 0 | 1.0000 | 0.5873 | 0.7400 |
| mPing | R3 minimap2 / minimap2 | 2599 | 0 | 1.0000 | 0.5776 | 0.7322 |
| riceTElib | RelocaTE2 | 2085 | 427 | 0.8300 | 0.4633 | 0.5947 |
| riceTElib | R3 BLAT / bwa aln | 2104 | 436 | 0.8283 | 0.4676 | 0.5977 |
| riceTElib | R3 BLAT / bwa-MEM | 2079 | 509 | 0.8033 | 0.4620 | 0.5866 |
| riceTElib | R3 Bowtie2 / bwa-MEM | 2093 | 432 | 0.8289 | 0.4651 | 0.5959 |
| riceTElib | R3 bwa-MEM / bwa-MEM | 1815 | 420 | 0.8121 | 0.4033 | 0.5390 |
| riceTElib | R3 minimap2 / minimap2 | 1782 | 397 | 0.8178 | 0.3960 | 0.5336 |
| Divergence | RelocaTE2 | 4928 | 1037 | 0.8262 | 0.1825 | 0.2990 |
| Divergence | R3 BLAT / bwa aln | 4962 | 1086 | 0.8204 | 0.1838 | 0.3003 |
| Divergence | R3 BLAT / bwa-MEM | 4971 | 1518 | 0.7661 | 0.1841 | 0.2969 |
| Divergence | R3 Bowtie2 / bwa-MEM | 4667 | 922 | 0.8350 | 0.1729 | 0.2864 |
| Divergence | R3 bwa-MEM / bwa-MEM | 4163 | 863 | 0.8283 | 0.1542 | 0.2600 |
| Divergence | R3 minimap2 / minimap2 | 3892 | 828 | 0.8246 | 0.1441 | 0.2454 |

BLAT/bwa aln preserves close aggregate parity across all panels. Relative to
RelocaTE2 it adds 7 TP and removes 1 FP on mPing; adds 19 TP and 9 FP on
riceTElib; and adds 34 TP and 49 FP on divergence. These aggregate comparisons
do not establish exact locus-by-locus family, coordinate, or status parity.
Nor do they quantify improvement against an archived pre-update full run.

## Accuracy and computational tradeoffs

On riceTElib, exact TSD accuracy among true detections is 96.4% for RelocaTE2,
96.2% for R3 BLAT/bwa aln, 95.9% for BLAT/bwa-MEM, 91.6% for
Bowtie2/bwa-MEM, 84.0% for bwa-MEM/bwa-MEM, and 87.1% for minimap2/minimap2.
On divergence the corresponding values are 85.7%, 85.7%, 85.1%, 74.6%,
62.7%, and 68.1%.

RiceTElib median wall times are 6.75 h for RelocaTE2, 6.19 h for R3
BLAT/bwa aln, 5.12 h for BLAT/bwa-MEM, 39.9 min for Bowtie2/bwa-MEM,
44.0 min for bwa-MEM/bwa-MEM, and 54.5 min for minimap2/minimap2.
R3 BLAT/bwa aln uses about 17.1 GiB median peak RSS on riceTElib versus
6.5 GiB for RelocaTE2 and 6.7 GiB for Bowtie2/bwa-MEM.
These are whole-adapter resource measurements, not isolated aligner timings.

Bowtie2/bwa-MEM is a promising speed/accuracy option on canonical riceTElib,
but it loses sensitivity on mPing and diverged elements and has lower exact
TSD accuracy. BLAT/bwa aln remains the strongest choice for RelocaTE2 parity.
All configurations lose most sensitivity at high divergence: at 10% divergence,
recall is 4.2% for RelocaTE2 and 4.4% for R3 BLAT/bwa aln; at 20%, both
are approximately 0.2%.

Somatic classification remains a shared weakness. On riceTElib, RelocaTE2
detects 873/2,700 somatic-event opportunities and R3 BLAT/bwa aln detects
879/2,700, but only 153/873 (17.5%) and 152/879 (17.3%) of those detections,
respectively, receive the correct somatic status. Close parity therefore does
not imply high biological classification accuracy.

## Dashboard validation and remaining correction

The real `reports/datasets.tsv` loads all three datasets through
`dashboard.data.loaders.load_report_suite`. All eight Streamlit pages were
executed with `streamlit.testing.v1.AppTest` against the actual reports.
They rendered without exceptions. Switching to riceTElib and divergence
also rendered without exceptions, with six callers selected on data pages.
The cross-dataset page exposes all five R3 variants and defaults to BLAT/bwa aln.

The dashboard deliberately uses mean per-sample precision, pooled recall,
and their harmonic mean for F1. It does not display the pooled precision/F1
used above. This matters especially for divergence: the dashboard shows F1
0.2927 for RelocaTE2 and 0.2905 for R3 BLAT/bwa aln, whereas pooled F1 is
0.2990 and 0.3003. This is an aggregation difference, not stale data. Sparse
high-divergence samples carry equal weight in mean per-sample precision.

One actual wiring defect remains: `dashboard/pages/05_provenance.py` hardcodes
`config/benchmark.toml`, so it displays the two-caller configuration instead of
the six-caller configuration used for these results. The pinned environment
manifest it displays also does not by itself prove the exact editable source
revision used in each completed run.

The managed filesystem rejected editing the sibling benchmark repository.
A correction is supplied in `2026-09-21-dashboard-config.patch` beside this
document. `git apply --check` succeeds in the benchmark repository. The patch
lets Provenance use `CONFIG`, retaining the default if unset, and labels the
displayed configuration as a current file rather than an archived run snapshot.
An in-memory Streamlit test of the configuration-selection change rendered
without exceptions and displayed all six enabled callers.

From the benchmark repository, apply and launch:

```bash
git apply --check ../../RelocaTE3_jason/RelocaTE3/docs/2026-09-21-dashboard-config.patch
git apply ../../RelocaTE3_jason/RelocaTE3/docs/2026-09-21-dashboard-config.patch
CONFIG=config/benchmark.full-aligners.toml bash pipeline/run_dashboard.sh --report-dir reports
```

Status: analysis and current dashboard rendering verified; provenance fix
prepared but not applied due to workspace permissions. No outputs deleted.
Next step: apply the provenance patch and launch the dashboard with the full
configuration explicitly selected.
