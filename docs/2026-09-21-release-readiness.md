# RelocaTE3 release-readiness assessment

Superseded for release decisions by [the September 28 deep dive](2026-09-28-release-deep-dive.md):
long-TSD coordinate matching distorts absolute benchmark accuracy, and the new
pairing candidate still requires validation. Historical values below are retained
as the record of the earlier assessment, not a current stable-release approval.

Recorded: 2026-09-21, 17:56 PDT (America/Los_Angeles).
Project: RelocaTE3, branch `main`, merged PR #49.
Source assessed: `942a5ea266bd000a47e27df7ef9e2536c877acbd`.
Benchmark repository: `../../relocate_benchmark/relocate-benchmark`, branch
`main`, HEAD `8ad522f` (PR #27), with local adapter/configuration changes.
Purpose: assess the default BLAT / `bwa aln` configuration against RelocaTE2
across the completed full benchmark and recommend a release decision.

## Decision

**Ready to freeze the current source as a release candidate; not yet ready to
recommend a stable public release pin without clean-artifact validation.**

The algorithm meets the project's immediate goal of close RelocaTE2 parity
on these simulated paired-end rice benchmarks. There is no evidence here
that another insertion-calling redesign is necessary before a release
candidate. The remaining release work is bounded: validate the installation
users will receive, capture its provenance, strengthen the default acceptance
gate, and refresh a representative real-rice check on that artifact.

This is a readiness recommendation, not a release operation. No source code,
tags, dependency pins, benchmark outputs, or published packages were changed.

## Scope and scoring

The completed full run used `config/benchmark.full-aligners.toml`, array
`28028811`, and aggregation job `28028812`. The previous review verified all
432 task completion markers and all 432 `.run_complete` files. This assessment
compares `relocate2` with `relocate3-blat-bwaaln`, not the separate
`relocate3-blat-bwa` configuration, which uses bwa-MEM.

Each caller has 72 sample runs: 9 mPing, 9 riceTElib, and 54 riceTElib-divergence
runs. Both ordinary panels cover 5x/15x/30x with three replicates; divergence
adds 0%, 2%, 5%, 10%, 15%, and 20% sequence divergence. The zero-divergence
condition overlaps the canonical riceTElib experiment and is not independent
biological replication.

The benchmark uses family-aware matching within 10 bp. Tables below pool
counts: precision = TP/(TP+FP), recall = TP/truth, F1 = 2TP/(truth+TP+FP).
Truth counts are sample-event opportunities, not independent unique loci.
False-positive counts come from `precision.tsv`, not duplicated class rows
in `correctness.tsv`. Resource measurements cover the full caller adapter.

The dashboard uses mean per-sample precision together with pooled recall,
so its F1 is a different statistic. For divergence it reports approximately
0.2927 for RelocaTE2 versus 0.2905 for RelocaTE3; pooled F1 is 0.2990 versus
0.3003. Consequently, "best performing" must specify the metric: BLAT / bwa
aln is the strongest R3 configuration by pooled F1 across these panels, not
a uniform winner under every aggregation or in every sample.

## Whole-panel results

Values are RelocaTE2 → RelocaTE3. Each ordinary panel contains 4,500 truth
opportunities; divergence contains 27,000.

| Panel | True positives | False positives | Precision | Recall | Pooled F1 |
|---|---:|---:|---:|---:|---:|
| mPing | 2,819 → 2,826 | 1 → 0 | 0.9996 → 1.0000 | 0.6264 → 0.6280 | 0.7702 → 0.7715 |
| riceTElib | 2,085 → 2,104 | 427 → 436 | 0.8300 → 0.8283 | 0.4633 → 0.4676 | 0.5947 → 0.5977 |
| Divergence | 4,928 → 4,962 | 1,037 → 1,086 | 0.8262 → 0.8204 | 0.1825 → 0.1838 | 0.2990 → 0.3003 |

These are small improvements in pooled sensitivity/F1, with a precision cost
in the complex-library panels. They establish close parity, not substantial
biological superiority or a statistically significant improvement. They do
not measure change relative to an archived pre-update full benchmark.

### Ordinary panels by coverage

Each row pools three replicates and 1,500 truth opportunities.

| Panel | Coverage | TP, R2 → R3 | FP, R2 → R3 | Recall, R2 → R3 | F1, R2 → R3 |
|---|---:|---:|---:|---:|---:|
| mPing | 5x | 610 → 616 | 0 → 0 | 0.4067 → 0.4107 | 0.5782 → 0.5822 |
| mPing | 15x | 988 → 987 | 1 → 0 | 0.6587 → 0.6580 | 0.7939 → 0.7937 |
| mPing | 30x | 1,221 → 1,223 | 0 → 0 | 0.8140 → 0.8153 | 0.8975 → 0.8983 |
| riceTElib | 5x | 449 → 451 | 103 → 102 | 0.2993 → 0.3007 | 0.4376 → 0.4394 |
| riceTElib | 15x | 739 → 747 | 145 → 147 | 0.4927 → 0.4980 | 0.6200 → 0.6241 |
| riceTElib | 30x | 897 → 906 | 179 → 187 | 0.5980 → 0.6040 | 0.6964 → 0.6988 |

Across individual samples, R3's F1 wins/ties/losses are 5/3/1 for mPing,
7/1/1 for riceTElib, and 20/19/15 for divergence. The largest sample-level
F1 decrease is approximately 0.0040 (divergence 15%, replicate 3, 30x).
These comparisons are descriptive, not significance tests.

### Divergence by sequence divergence

Each row pools nine samples and 4,500 truth opportunities.

| Divergence | TP, R2 → R3 | FP, R2 → R3 | Precision, R2 → R3 | Recall, R2 → R3 | F1, R2 → R3 |
|---|---:|---:|---:|---:|---:|
| 0% | 2,085 → 2,104 | 427 → 436 | 0.8300 → 0.8283 | 0.4633 → 0.4676 | 0.5947 → 0.5977 |
| 2% | 1,728 → 1,731 | 368 → 380 | 0.8244 → 0.8200 | 0.3840 → 0.3847 | 0.5240 → 0.5237 |
| 5% | 894 → 901 | 174 → 184 | 0.8371 → 0.8304 | 0.1987 → 0.2002 | 0.3211 → 0.3226 |
| 10% | 191 → 196 | 51 → 56 | 0.7893 → 0.7778 | 0.0424 → 0.0436 | 0.0806 → 0.0825 |
| 15% | 22 → 22 | 13 → 22 | 0.6286 → 0.5000 | 0.0049 → 0.0049 | 0.0097 → 0.0097 |
| 20% | 8 → 8 | 4 → 8 | 0.6667 → 0.5000 | 0.0018 → 0.0018 | 0.0035 → 0.0035 |

At 15% and 20%, F1 decreases slightly before rounding. R3 adds false positives
without adding true positives. Both tools have very low sensitivity here;
parity is not evidence of satisfactory divergent-element discovery.

All divergence/coverage cells are shown below as R3 minus R2 counts. These
retain the interactions hidden by pooling coverage.

| Divergence | 5x: ΔTP / ΔFP | 15x: ΔTP / ΔFP | 30x: ΔTP / ΔFP |
|---|---:|---:|---:|
| 0% | +2 / −1 | +8 / +2 | +9 / +8 |
| 2% | +1 / −1 | +2 / +5 | 0 / +8 |
| 5% | −1 / +2 | +4 / +2 | +4 / +6 |
| 10% | +1 / 0 | 0 / +2 | +4 / +3 |
| 15% | 0 / +1 | 0 / +3 | 0 / +5 |
| 20% | 0 / +1 | 0 / +1 | 0 / +2 |

## Event-level parity

Within each sample, join the callers' `matches.tsv` on truth `event_id` and
compare detected truth events. This tests whether small net gains conceal
large exchanges of detected loci.

| Panel | Shared true detections | R2-only | R3-only | Identical position among shared | Identical status among shared | Identical TSD among shared |
|---|---:|---:|---:|---:|---:|---:|
| mPing | 2,818 | 1 | 8 | 2,817 | 2,815 | 2,814 |
| riceTElib | 2,084 | 1 | 20 | 2,084 | 2,082 | 2,073 |
| Divergence | 4,919 | 9 | 43 | 4,918 | 4,917 | 4,860 |

R3 retains 99.96%, 99.95%, and 99.82% of R2's true detections, respectively.
Positions and status agree for more than 99.8% of shared true detections in
each panel. This does not establish byte-identical output, raw family-label
equality, shared false-positive parity, or truth accuracy of agreed statuses.

## Biological limitations and subgroup results

On canonical riceTElib, R3's recall is equal or higher in nine of ten TE
groups. SINE loses one detection (recall 0.2267 → 0.2244). Both tools remain
weak for Helitron (0.0156 → 0.0289) and LINE (0.1511 → 0.1511), so average
performance should not be represented as uniform across TE biology.

Exact TSD accuracy among true detections is approximately 98.1% → 97.9%
for mPing, 96.4% → 96.2% for riceTElib, and 85.7% → 85.7% for divergence.
For R3 it declines to about 70.3% at 5% divergence and 44.9% at 10%.
Caller-to-caller agreement does not remove this shared loss of truth accuracy.

On riceTElib, homozygous-event recall is 0.7044 → 0.7156 and heterozygous
recall is 0.6422 → 0.6456. Somatic recall is 0.3233 → 0.3256, but only
153/873 (17.5%) and 152/879 (17.3%) of detected somatic events receive the
correct status. A release can accurately claim RelocaTE2-like behavior; it
should not claim validated high-accuracy somatic genotyping from this evidence.

## Runtime and memory

RSS conversion assumes the Linux resource report's KB values are KiB;
divide by 1,048,576 for GiB. Values are R2 → R3, across samples in each panel.

| Panel | Median wall time (hours) | Median peak RSS (GiB) | Maximum peak RSS (GiB) |
|---|---:|---:|---:|
| mPing | 0.60 → 0.49 | 6.48 → 6.72 | 6.49 → 6.72 |
| riceTElib | 6.75 → 6.19 | 6.48 → 17.15 | 6.49 → 33.81 |
| Divergence | 5.84 → 6.38 | 6.48 → 17.42 | 6.49 → 36.97 |

R3 is not universally faster and has a substantial complex-library memory
cost. Recommend a 64-GB allocation for workloads comparable to the evaluated
high-coverage riceTElib runs, not as a universal minimum for all datasets.

## Release verification still needed

1. **Validate a clean, non-editable release artifact and its version.** Both
   current environments import the live checkout, but their installed metadata
   reports older versions: default `0.1.0.post77+g0ae78a2`, BLAT
   `0.1.0.post137+ga0d0baf.d20260812`. This is consistent with stale editable
   installation metadata, not proof of stale imported source. Build a wheel
   and source distribution from the candidate, install outside the checkout,
   verify imports/entry points/version, and test the supported default workflow.
   Let the existing versioningit configuration derive the version; do not
   hand-edit it. This assessment did not perform that build or installation.

2. **Make the source/environment/benchmark provenance reproducible.** The
   intended source is current `main`, but the completed run does not archive a
   per-task source hash that proves it. The benchmark adapter prefers the
   editable source environment and falls back to an older `6fe6d5d` pin.
   Preserve/commit the local adapter and full configuration, update the fallback
   for the candidate, and record source revision, clean/dirty state, artifact
   checksum, environment lock, tool versions, commands, and input checksums for
   release validation. Do not retrospectively claim certainty about the exact
   historical source from today's checkout or package version string.
   Also reconcile or explicitly document R2's configured `size=250` versus
   R3's default `--size=500` for mate-only insertion spans; these are the
   current-config results, not proof of parameter-identical implementations.

3. **Turn the default's observed performance into a durable release gate.**
   The full SLURM suite completed successfully on September 2
   (`logs/relocate3-tests.28015204.log`), with two seqtk-dependent skips.
   Exercise the optional seqtk path as well as the Python fallback during
   artifact validation. In `tests/acceptance_test.py`, the stronger quantitative
   gate targets minimap2 (at least 170/200 detections and precision at least
   0.85); the default `run-all` smoke test only requires more than 100/200
   detections. Establish baseline-appropriate precision/recall gates for BLAT /
   bwa aln without borrowing the minimap2 thresholds blindly. Test the installed
   CLI, since the full benchmark uses a staged adapter rather than simply
   invoking the end-user `run-all` command. Submit substantial tests via SLURM.

4. **Refresh a representative real-rice check and release limitations.**
   A real-rice harness and historical results exist under
   `validation/real_rice/`, but the inspected result/report artifacts predate
   PR #49. Run a representative real sample with the candidate artifact and
   compare non-reference calls/status against the existing baseline. This is
   a real-data regression check, not a precision estimate without validated
   truth. Describe the release as validated for the tested Linux paired-end
   rice workflow; document divergence, somatic classification, memory, and
   family-label semantics. Do not infer validation of every optional backend,
   operating system, single-end input, reference-insertion or excision feature.

These are release-validation gates, not evidence of a newly demonstrated
calling bug. No second full six-caller/432-task benchmark is needed merely to
repeat unchanged algorithm comparisons. First validate the exact artifact on
focused representative inputs; expand reruns if its results disagree or if
the calling algorithm/configuration changes.

## Provenance, commands, and review status

Source inputs for this assessment, relative to the benchmark root:

- `reports/datasets/{mping,ricetelib,ricetelib_divergence}/correctness.tsv`
- `reports/datasets/{mping,ricetelib,ricetelib_divergence}/precision.tsv`
- `reports/datasets/{mping,ricetelib,ricetelib_divergence}/resources.tsv`
- `reports/datasets/<dataset>/per_sample/<caller>/<sample>/matches.tsv`
- `config/benchmark.full-aligners.toml`, `callers/relocate3/{run.sh,env.sh,pixi.toml}`

Lightweight review commands, from the RelocaTE3 repository:

```bash
git status --short
git branch --show-current
git rev-parse HEAD
sed -n '1,80p' logs/relocate3-tests.28015204.log
rg -n 'run_all_cli_end_to_end|170|100|0.85' tests/acceptance_test.py
```

Tables were calculated by grouping the report counts by caller/panel,
coverage, divergence, biological class, and TE group; pairing samples for
F1 differences; and joining each sample's truth-event matches as described
above. No alignments, simulations, full test suites, or new jobs were run.

The earlier dashboard review is in `2026-09-21-full-benchmark-review.md`.
At this assessment's final check, the benchmark worktree now contains the
CONFIG-aware provenance-page change supplied in that earlier review; this
turn did not apply it. The page labels its configuration as a current file,
not an archived run snapshot. Dashboard wiring cannot substitute for run-time
source/artifact provenance.

Status: all-panel assessment complete; no newly observed benchmark task
failures. Clean artifact installation and current real-data validation remain
unperformed in this review. Next action: prepare a focused, reproducible
release-candidate validation workflow for the frozen source before choosing
and publishing the stable release tag.
