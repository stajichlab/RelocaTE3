# RelocaTE3 release deep dive and pairing-correction candidate

Follow-up: [corrected coordinate rescoring](2026-09-28-coordinate-rescoring-results.md)
is now complete for the historical R2/R3 calls across all three datasets.
Use those corrected accuracy tables for current comparisons; the old-scoring
values below remain a record of the diagnosis, not current biological accuracy.

Recorded: 2026-09-28, afternoon PDT (America/Los_Angeles).
Project: RelocaTE3, `main`, base `942a5ea`, following merged PR #49.
Status: candidate correction implemented; 34 focused tests pass. Historical
full-benchmark analysis complete. HPC submission timed out without a job ID:
**submission is unconfirmed**, not known running or complete. No release tag,
commit, benchmark deletion, or official rescore was performed.

## Decision

**Do not pin a stable release with an unconditional “no performance sacrifice”
claim yet.** Detection parity is excellent, but resource tradeoffs remain, the
pairing change needs validation, and a newly verified scoring-coordinate
mismatch materially distorts absolute accuracy.

Harmonize benchmark coordinates and rescore both tools symmetrically before
tuning filters against apparent false positives. This supersedes the earlier
interpretation that most LINE/SINE errors represented discovery failure.
Somatic detection and somatic labeling also require separate claims.

## Critical finding: long TSDs exceed the scoring window

The canonical simulator extracts the TSD as
`reference[position - tsd_length : position]` and inserts the TE plus a new
copy of the TSD after `reference[:position]`. For nonzero TSD length L, its
original leftmost reference base is **position - L + 1**, 1-based.

The caller normalizer exports the *start* of the reported TSD interval as
`position`. The scorer compares it directly to the simulator's insertion
anchor within 10 bp. A correct insertion with L > 11 can therefore fail solely
because these coordinates describe different edges.

Verified source, under simulator root
`/bigdata/wesslerlab/shared/Rice/Nathan/rice/make_simulated_genome/make_simulation_new`:

- `pipeline/make_riceTElib_benchmark/build_multite_panel.py`, lines 555 and
  678–682: TSD extraction and sequence insertion.
- Divergence builder references this multi-TE builder at line 33.
- Benchmark `lib/calls.py` exports the interval start;
  `scoring/score_calls.py` compares it directly with truth `position`.

Example: Os1605 SINE at truth Chr1:488190 has a 12-base TSD. The call at
488179 is exactly the left TSD boundary, but its 11-bp offset fails scoring.

I counted existing FP/FN pairs with the same normalized family and **exact**
agreement with `truth_position - tsd_length + 1`, one-to-one per sample:

| Panel | Reported FP, R2 / R3 | Exact long-TSD anchor FP/FN pairs, R2 / R3 |
|---|---:|---:|
| mPing | 1 / 0 | 0 / 0 |
| riceTElib | 427 / 436 | 408 / 408 |
| Divergence | 1037 / 1086 | 928 / 922 |

Canonical pairs comprise 231 LINE and 177 SINE opportunities per tool:
95.6% of R2's and 93.6% of R3's reported canonical FPs. This conservative
diagnostic is **not a complete rescore**. Do not publish revised accuracy by
simply subtracting counts: other matches can change and coordinate semantics
for unknown TSDs, one-sided calls, and no-TSD elements must be specified.

These pairs include 160 true somatic opportunities per tool in canonical
riceTElib, including 28 correctly labeled somatic. Divergence has 309/310 such
somatic opportunities, including 76 correctly labeled per tool. Thus the old
scorer also understates somatic-label precision. The shared tendency to label
true somatic events heterozygous remains after recognizing these pairs.

Helitrons are a separate issue: all 50 canonical Helitron truth events have
`tsd=NONE`, `tsd_length=0`. Their zero exact-TSD score is not an ordinary
sequence mismatch. The two new pairing fixtures are Helitron-labeled, so
inventing a positive duplication to rescue them would be inappropriate.

## Historical panel comparison

The following values reproduce the **existing scorer**, not corrected
absolute accuracy. R3 means BLAT / bwa aln, not BLAT / bwa-MEM. There are
9 mPing, 9 riceTElib, and 54 divergence samples per caller: 4,500, 4,500, and
27,000 truth-event opportunities. These reuse loci across conditions; the
zero-divergence panel overlaps canonical riceTElib and is not independent.

| Panel | TP, R2 → R3 | FP, R2 → R3 | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|
| mPing | 2819 → 2826 | 1 → 0 | 99.96% → 100% | 62.64% → 62.80% | .7702 → .7715 |
| riceTElib | 2085 → 2104 | 427 → 436 | 83.00% → 82.83% | 46.33% → 46.76% | .5947 → .5977 |
| Divergence | 4928 → 4962 | 1037 → 1086 | 82.62% → 82.04% | 18.25% → 18.38% | .2990 → .3003 |

Metrics pool counts, rather than mixing mean sample precision with pooled
recall as the dashboard does. Sample-level F1 wins/ties/losses are 5/3/1,
7/1/1, and 20/19/15. Worst F1 decrease is about .0040 at 15% divergence,
replicate 3, 30x. Small pooled gains do not establish uniform superiority.

R3 retains 2818/2819, 2084/2085, and 4919/4928 of R2's scored true detections:
99.96%, 99.95%, and 99.82%. It gains 8, 20, and 43 others. Seven of nine
divergence losses are within-window family-vote ties, not absent breakpoints.
These paired comparisons should also be regenerated after harmonization.

## Somatic detection is not correct somatic classification

Detection recall divides by all somatic truth events. Label accuracy divides
correct somatic labels by detected somatic events. End-to-end labeled recall
divides correct somatic labels by all somatic truth events.

Current-scoring results (R2 → R3):

| Panel | Somatic truth | Detected | Detection recall | Correct label | Label accuracy among detections | End-to-end labeled recall |
|---|---:|---:|---:|---:|---:|---:|
| mPing | 2700 | 1190 → 1194 | 44.07% → 44.22% | 249 → 250 | 20.92% → 20.94% | 9.22% → 9.26% |
| riceTElib | 2700 | 873 → 879 | 32.33% → 32.56% | 153 → 152 | 17.53% → 17.29% | 5.67% → 5.63% |
| Divergence | 16200 | 1824 → 1833 | 11.26% → 11.31% | 448 → 447 | 24.56% → 24.39% | 2.77% → 2.76% |

Canonical R3 labels 716/879 detected somatic events heterozygous, versus
709/873 for R2. R3 labels another 10 homozygous and one
homozygous/excision_no_footprint. This is predominantly shared behavior, not
evidence of a memory-patch regression.

Correct somatic labels divided by *all calls labeled somatic*, including FPs,
are 100% → 100% on mPing, 83.61% → 82.16% on riceTElib, and
81.45% → 79.26% on divergence **under the old scorer**. The latter values are
contaminated by the anchor issue and must not be advertised as validated
biological precision.

### Coverage and cellular fraction

Each cell has 300 somatic truth opportunities. Entries are **detected / correctly
labeled somatic**, R2 → R3. Simulator cellular fractions 10%, 20%, 40%
correspond to expected VAF 5%, 10%, 20%; they are not interchangeable labels.

| Panel / cellular fraction | 5x | 15x | 30x |
|---|---|---|---|
| mPing / 10% | 10/0 → 10/0 | 55/7 → 55/7 | 140/109 → 140/109 |
| mPing / 20% | 31/0 → 31/0 | 141/11 → 141/11 | 231/98 → 233/100 |
| mPing / 40% | 84/0 → 87/0 | 224/5 → 223/4 | 274/19 → 274/19 |
| riceTElib / 10% | 11/0 → 11/0 | 44/3 → 44/3 | 106/76 → 105/76 |
| riceTElib / 20% | 22/0 → 22/0 | 97/6 → 99/6 | 150/56 → 150/56 |
| riceTElib / 40% | 68/0 → 68/0 | 171/4 → 173/4 | 204/8 → 207/7 |

Neither tool correctly labels any of these somatic detections at 5x. At 30x,
R3 detects 207/300 canonical 40%-cellular events but labels only seven somatic.
This is consistent with the shared count-threshold classifier, not a calibrated
cell-fraction model. `Characterizer._classify` includes a rule assigning
heterozygous when both flanker and spanner counts exceed ten. The thresholds
are not posterior probabilities or VAF estimates. No genotype rule was changed.

The benchmark's earlier issue audit also documents that simulator
`observed_vaf` is a junction-support ratio with two insertion-junction
opportunities versus one reference-spanning opportunity. It was not used here
as an unbiased measured VAF.

### Divergence and somatic insertions

Each row has 2,700 somatic truth opportunities:

| Divergence | Detected, R2 → R3 | Correctly labeled, R2 → R3 |
|---|---:|---:|
| 0% | 873 → 879 | 153 → 153 |
| 2% | 642 → 643 | 176 → 174 |
| 5% | 266 → 268 | 98 → 98 |
| 10% | 36 → 37 | 18 → 19 |
| 15% | 5 → 4 | 2 → 2 |
| 20% | 2 → 2 | 1 → 1 |

R3's zero-divergence label count differs by one from its separate canonical
run despite equal detection totals. Do not assume byte-identical
characterization or merge these runs silently. Higher conditional accuracy
on a tiny surviving subset is not improved overall somatic performance.

## Germline, TE groups, divergence, and TSDs

Germline recall under current scoring:

| Panel | Homozygous, R2 → R3 | Heterozygous, R2 → R3 |
|---|---:|---:|
| mPing | 95.44% → 95.44% | 85.56% → 85.89% |
| riceTElib | 70.44% → 71.56% | 64.22% → 64.56% |
| Divergence | 30.69% → 31.00% | 26.80% → 26.94% |

Among detected canonical germline events, correct homozygous labels are
581/634 → 591/644; correct heterozygous labels are 444/578 → 445/581.
The exact status scorer keeps `homozygous/excision_no_footprint` separate;
this label does not itself establish that a biological excision occurred.

Canonical TE-group TP counts, with 450 truth opportunities per row:

| Group | R2 | R3 |
|---|---:|---:|
| CACTA | 253 | 254 |
| Helitron | 7 | 13 |
| LINE | 68 | 68 |
| LTR Copia | 272 | 272 |
| LTR Gypsy | 264 | 266 |
| MULE | 282 | 283 |
| PIF/Harbinger | 280 | 281 |
| SINE | 102 | 101 |
| Tc1/Mariner | 280 | 287 |
| hAT | 277 | 279 |

LINE/SINE totals omit the exact-anchor pairs discussed above. Helitron's
no-TSD limitation is distinct. Divergence R3 has one fewer scored LINE
detection overall (140 versus 141), despite its pooled gain.

All-event recall at divergence 0/2/5/10/15/20% is, for R2,
46.33/38.40/19.87/4.24/.49/.18%, and for R3,
46.76/38.47/20.02/4.36/.49/.18%. Harmonization changes absolute rates but
does not erase the observed decline. At 15%/20%, R3 adds FP without adding
total TP under current scoring.

Exact TSD sequence among scored detections is 98.12% → 97.88% on mPing,
96.40% → 96.20% on riceTElib, and 85.73% → 85.65% on divergence. Exact
counts are 2766 → 2766, 2010 → 2024, and 4225 → 4250: a conditional rate
can fall while absolute correct detections increase. Literal caller-start /
simulator-anchor equality is not single-base accuracy: even mPing has zero
such matches because they describe different edges.

## Runtime and memory

Historical whole-adapter measurements, R2 → R3:

| Panel | Median wall hours | Median peak GiB | Maximum peak GiB |
|---|---:|---:|---:|
| mPing | .603 → .487 | 6.48 → 6.72 | 6.49 → 6.72 |
| riceTElib | 6.75 → 6.19 | 6.48 → 17.15 | 6.49 → 33.81 |
| Divergence | 5.84 → 6.38 | 6.48 → 17.42 | 6.49 → 36.97 |

Median *paired-sample* runtime ratios R3/R2 are .665, .917, and 1.058.
R3 is slower in 0/9, 3/9, and 37/54 runs. The largest observed ratio is
1.435 at 5% divergence, replicate 3, 5x. Node/load differences prevent treating
these as controlled intrinsic speed estimates; the adapter includes external
tools and all processing stages.

Subsequent memory optimization reduced one canonical 30x replay from 34.45 to
19.10 GiB, preserving calls, with 279 tests passing. This is substantial but
still above R2's approximately 6.5 GiB. Its latest runtime was 11h51m, not a
demonstrated speed improvement. Do not extrapolate one sample's savings to all
72 workloads or lower global allocations without validation; retain the
64-GB recommendation for comparable high-coverage full-workflow runs.

## Pairing correction and validation status

`_pair_breakpoints` now uses inclusive-left + 1 for distance/order comparisons,
including the 100-bp threshold, while returning unchanged inclusive coordinates.
No output-coordinate shift, family rule, new threshold, or CLI switch was added.
Non-positive TSD geometry remains rejected; no-TSD discovery is a separate
feature decision. The local Os2571 sensitivity cost remains explicit.

Both real fixtures now pass without expected-failure markers. Controls cover
true one-base TSDs, nearby independent insertions, one-sided candidates, and
100/101-bp boundaries on both sides. **34 focused insertion/parity tests
passed**, and `git diff --check` passed. Full-suite and all-sample candidate
validation are pending.

Prepared replay: `results/pairing-replay/2026-09-28/`, 72 tasks, capped at four
concurrent tasks. Each reuses stored flank, untrimmed-junction, and all-read
BAMs, runs baseline/candidate insertion calling and characterization, and
saves normalized calls and scores. No TE alignment, trimming, or genome
realignment is repeated. Baseline normalized calls must exactly reproduce
historical family/position/TSD/strand/status records before the candidate runs.
A mismatch fails the task so historical drift cannot masquerade as patch effect.

Source, tests, normalizer, scorer, and configuration are frozen with hashes.
Baseline uses the pre-patch HEAD insertion finder; other source is identical
between variants. Large inputs have size/mtime guards, not whole-BAM hashes.
Output guards refuse repeated/partial tasks. The replay retains the old scorer
for continuity; preserved calls can be rescored once the coordinate contract
is fixed.

SLURM: epyc, 32 GB, one CPU, eight hours per paired task. Cached
`/var/spool/slurmd/conf-cache/slurm.conf` gives epyc MaxTime=30-00:00.
`scontrol` timed out. Submission with a 15-second timeout exited 124 without
a numeric ID; queue inspection also could not be confirmed. **Do not assume
the replay started or resubmit without checking from a normal cluster terminal.**

From the RelocaTE3 root:

```bash
squeue -u "$USER" -o '%.18i %.40j %.10T %.10M %.30R'
```

If no matching replay exists:

```bash
sbatch --parsable --array=0-71%4 \
  results/pairing-replay/2026-09-28/replay_pairing_benchmark.slurm \
  results/pairing-replay/2026-09-28
```

Do not delete old outputs. The existing `sbatch scripts/run_full_tests.slurm`
can additionally validate the full suite; record source and inspect failures
and skips, rather than treating submission as success.

## Release gates

1. Define/test a common coordinate contract for truth and both callers,
   including long/unknown TSDs, one-sided calls, and no-TSD elements. Rescore
   archived calls for all datasets and regenerate subgroups/dashboard reports.
   Do not simply widen the window to recover known truth labels.
2. Complete baseline/candidate replay and examine every lost/gained true event,
   including coverage, divergence, TE-group, and somatic strata. A pooled gain
   must not silently compensate for a subgroup regression.
3. Validate a clean installed artifact and end-user CLI, not only editable
   source and a staged adapter. Capture versions, source and tool hashes;
   complete full tests and a representative current real-rice regression.
4. Define acceptable noninferiority margins before asserting equivalence.
   None were prespecified. Three replicates and repeated loci do not justify
   naive event-wise significance tests or universal no-sacrifice claims.

An interim description is “a modern implementation with close RelocaTE2-like
behavior on evaluated rice benchmarks,” with resource/genotyping limitations.
It is not yet “uniformly at least as accurate, fast, and memory-efficient.”

## Reproducibility

From the RelocaTE3 root, completed lightweight analyses:

```bash
python scripts/release_benchmark_audit.py \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --outdir results/release-audit/2026-09-28-full-benchmark
python scripts/audit_truth_anchors.py \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --outdir results/release-audit/2026-09-28-anchor-diagnostic
.pixi/envs/default/bin/python -m pytest -q \
  tests/duplicate_loci_test.py tests/insertions_test.py::TestPairBreakpoints \
  tests/insertions_tsd_parity_test.py tests/insertions_tsd_class_parity_test.py
```

Analysis scripts refuse existing output directories. Outputs include pooled
and sample counts, paired comparisons, coverage/divergence/TE-group/status/
fraction subgroups, full genotype confusion tables, somatic-label denominators,
and SHA-256 provenance. The anchor diagnostic lists every event/call pair.
No benchmark or simulator files were modified; no substantial analysis or
alignment ran on the login node.
