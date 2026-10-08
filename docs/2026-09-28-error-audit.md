# Full-benchmark error audit: RelocaTE2 versus RelocaTE3

Subsequent finding: [the release deep dive](2026-09-28-release-deep-dive.md)
verified a truth/caller coordinate mismatch that explains most scored LINE/SINE
FP/FN pairs. Preserve the historical counts below, but do not interpret all of
them as biological caller errors or use them to tune filters before rescoring.

Recorded: 2026-09-28, America/Los_Angeles (analysis begun after 10:13 PDT).
Project: RelocaTE3, `main`, following merged PR #49; memory improvements remain
uncommitted. Status: report-level audit complete, selected indexed alignments
inspected; no production algorithm changes or new benchmark submissions.

## Scope and method

Compare `relocate2` with `relocate3-blat-bwaaln` across all 72 paired sample
runs. These are the existing full-benchmark results, not a full rerun of the
memory-optimized source. The memory changes have separately preserved calls
on one representative sample.

Truth detections are joined by dataset, sample, and event ID. Scoring requires
the family name to agree within 10 bp. The benchmark strips the suffix after
the first `#`, lowercases, and removes underscores. For RelocaTE2 compound
labels, this effectively selects the first family, not every listed family.
Thus family disagreements can produce a false positive and false negative
even when a caller locates an insertion correctly.

FP calls are paired one-to-one between tools by normalized family and position
within 10 bp, exact positions first, then increasing distance. This is a
descriptive greedy correspondence, not biological adjudication or rescoring.
Do not interpret family-specific unmatched FPs as necessarily distinct loci.
Nearest-call/raw-table annotations include their distance and are not evidence
that a distant record caused the error. Raw table evidence is reported read
support, not a fresh reconstruction from every BAM.

## Counts

| Panel | Shared TP | R2-only TP | R3-only TP | Shared FP | R2-only FP | R3-only FP |
|---|---:|---:|---:|---:|---:|---:|
| mPing | 2,818 | 1 | 8 | 0 | 1 | 0 |
| riceTElib | 2,084 | 1 | 20 | 422 | 5 | 14 |
| Divergence | 4,919 | 9 | 43 | 1,000 | 37 | 86 |

Both tools miss 1,673, 2,395, and 22,029 truth-event opportunities,
respectively. Most errors are therefore shared, not new regressions in R3.
The divergence panel includes the canonical zero-divergence experiment;
counts across panels are not independent biological replication.

Six of the 14 riceTElib R3-only family-specific FPs, and 38 of 86 divergence
FPs, have some R2 call within 10 bp regardless of family. The latter includes
both R2 true detections and R2 false positives. They are not all new loci.

## 1. Most divergence losses are family ambiguity, not missing breakpoints

Seven of the nine divergence R2-only detections have an R3 call within 10 bp
of truth with a different primary family. All seven have equal top junction
votes and `TE_family_status:ambiguous`. R3's deterministic alphabetical
tie-break is implemented in `_te_family_evidence` in `insertions.py`.

| Sample | Event | Truth family | R3 primary family | Junction votes | Offset from truth |
|---|---|---|---|---|---:|
| div002_rep02_cov15x | TE000215 | Os1415 | Os0303_complete | 3:3 | -1 |
| div002_rep02_cov30x | TE000093 | Os1403 | Os1100 | 1:1 | -8 |
| div005_rep01_cov30x | TE000484 | Os3576 | Os0029 | 1:1 | +1 |
| div005_rep02_cov15x | TE000030 | Os2507 | Os0596 | 1:1 | -6 |
| div005_rep03_cov5x | TE000028 | Os0269 | Os0043 | 1:1 | -1 |
| div010_rep02_cov15x | TE000277 | Os2224 | Os1218 | 1:1 | +1 |
| div015_rep03_cov30x | TE000212 | Os3391 | Os0596 | 1:1 | -9 |

Supporting-mate votes are not a universal solution: they favor the truth
family Os2224 in one case but the non-truth family Os0303_complete in another.
Changing tie order to match these truth labels would not demonstrate a better
biological classifier. Retain the single primary family plus ambiguity metadata;
do not restore compound labels just to improve the benchmark score.

The other two divergence losses are the repeated zero-divergence SINE loss
below, and `div002_rep03_cov30x` TE000085: R3 reports Os0029 at Chr1:22171544,
39 bp before the Os1279 truth position. The raw candidate spans through the
truth coordinate, with Os0029:11 versus Os1279:2 junction votes and unknown TSD.
This is a mixed-family/coordinate case, not a tied-label case.

## 2. Duplicated neighboring calls are a bounded precision target

Canonical riceTElib contains two R3-only FPs immediately beside an R3 true
detection of the same family:

- `cov30x_rep1`, Os2571: resolved call Chr1:24775704 with L=21/R=1;
  extra `supporting_junction` call Chr1:24775705 with L=0/R=26.
- `cov30x_rep2`, Os0980: resolved call Chr1:9006303 with L=32/R=1;
  extra `supporting_junction` call Chr1:9006304 with L=0/R=32.

Truth matching consumes the resolved call and leaves the neighboring call as
an FP. These have abundant reads and unanimous family votes: a generic minimum
read count or family-confidence threshold would not remove them.

The divergence panel has five same-family, within-10-bp unmatched calls:
the two above at 0%, the same two at 2%, and a resolved Os2881 call at
Chr1:240679 in `div002_rep03_cov30x`, beside truth at 240678.

There are also extra one-sided Os2781 calls at Chr1:5936699 in canonical
30x replicates 2 and 3, adjacent to correctly detected Os1422 at
Chr1:5936696. These are different-family cases and should not be blindly
collapsed under a same-family rule.

Recommendation: trace breakpoint pairing, same-start consolidation, and
conversion to `supporting_junction` for the two clear same-family duplicates.
The code currently consolidates equal starts, then arbitrates candidates,
then converts accepted one-sided candidates. The observed adjacent records
motivate inspecting that sequence; this review has not proven which step
created the duplicate or that the calls reuse identical read names.
Do not implement arbitrary distance-based deduplication without checking
evidence overlap and preserving genuinely adjacent insertions.

## 3. Two canonical missed detections have distinct evidence

### mPing: upstream alignment difference

`cov15x_rep2`, TE000497, truth Chr1:28969859. R2 reports TCCA at 28969856
with L=3/R=1. Indexed 300-bp BAM windows show the three left-junction reads
in both tools. R2 additionally retains read
`cov15x_rep2:clone40:h1:s889027347:Chr1-27214/1:start:5`, a 17M alignment at
28969856, flag 99, MAPQ 0, XT:R, X0=2, X1=78. That junction alignment is
absent from the inspected R3 locus. Its mate also has different pairing and
alignment in R3. This is not established as a downstream family-filter bug;
recovering an ambiguous 17-base alignment may carry a precision cost.

### SINE: junction evidence retained, downstream explanation unresolved

`cov30x_rep2`, TE000386, Os3874, truth Chr1:20045776. R2 reports a
`supporting_junction` at 20045777, L=2/R=0, ST=4/SR=4. Both indexed BAM
windows contain the same 66M and 42M junction alignments ending at 20045776,
with MAPQ 60 and 37. R3 does not report the insertion. The 300-bp inspection
does not enumerate all mates in the complete clustering window; missing
support, reference-edge filtering, and candidate arbitration remain to be
distinguished. The matching zero-divergence sample repeats this loss.

## 4. Additional FPs are heterogeneous

The 14 canonical R3-only family-specific FPs comprise six
`supporting_junction`, three `UNK`, and five resolved-TSD records. The 86
divergence records comprise 25, 21, and 40 respectively.

Repeated off-truth-chromosome calls include Os3912 at Chr4:19264325
(L=19/R=1 in canonical 30x replicates 1 and 3) and Os3442 at Chr3:31981544
(L=1/R=28 in canonical 30x replicate 3). Their asymmetric junction counts
motivate checking reference-copy edges and alignment evidence. They are not
explained by low total read support or ambiguous family votes alone.

Do not remove all one-sided calls: even mPing has 52 currently scored true
detections with one-sided raw evidence. No counterfactual full rescoring or
proposed threshold optimization was performed here.

## Reproducibility and limitations

Audit outputs: `results/error-audit/2026-09-28-full-benchmark/`:

- `summary.tsv`: reconciled sample-level totals pooled by dataset.
- `discordant_truth_events.tsv`: all lost/gained truth detections and nearest
  calls, including reported R3 read counts and family votes.
- `false_positive_audit.tsv`: every FP, correspondence category and context.
- `provenance.json`: SHA-256 hashes of all tables read and the audit script.
- `.complete`: successful audit marker.

From the RelocaTE3 repository root, the completed command was:

```bash
python scripts/audit_benchmark_errors.py \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --outdir results/error-audit/2026-09-28-full-benchmark
```

The script refuses an existing output directory. Choose a new directory to
repeat it; no cleanup is necessary. It asserts paired sample/truth identities
and reconciles normalized-call, TP, and FP counts with the precision reports.
Seven lightweight synthetic checks covered empty matching, duplicate handling,
10-bp boundaries, family mismatch/normalization, and empty nearest lookup.
`git diff --check` passed. This small-table analysis and four indexed BAM
queries did not run large tests or realignment on the login node.

Cluster `samtools/1.22.1` was available but module loading failed because its
logging hook could not access `/dev/log` in the managed session. Indexed BAM
inspection used existing project Python/pysam 0.24.0 instead; no environment
was created. The inspected BAM paths are
`runs/<dataset>/relocate2/<sample>/raw/repeat/bwa_aln/MSU_r7.repeat.bwa.sorted.bam`
and `runs/<dataset>/relocate3-blat-bwaaln/<sample>/raw/<sample>.repeat.bwaaln.sorted.bam`
under the benchmark root. Queries used `fetch("Chr1", position-150, position+150)`.

An exploratory exact-start join of characterized calls to raw rows encountered
a coordinate absent from the raw table; no incomplete filter-cost estimate
from that join is used above. Nearest-row annotations retain explicit offsets
and must be checked before interpreting evidence. The scorer and normalizer
were read directly; benchmark outputs and scoring rules were not modified.

## Next action

First reconstruct the two same-family duplicate loci from their complete
local read clusters and build regression fixtures. Then propose a narrow,
evidence-aware correction if the same insertion is being emitted twice.
Validate against true adjacent insertions and single-sided true calls before
an HPC replay. Keep family-tie semantics and the two canonical sensitivity
losses as separate investigations. A release pin remains premature until
these decisions and the previously identified artifact checks are resolved.
