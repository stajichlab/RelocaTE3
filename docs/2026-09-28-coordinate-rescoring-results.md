# Coordinate-contract correction and three-panel rescoring

Recorded: 2026-09-28 19:54–20:10 PDT (America/Los_Angeles).
Project: RelocaTE3, `main`, following PR #49. Benchmark: sibling
`../../relocate_benchmark/relocate-benchmark`, `main`, with existing user edits.
Status: corrected scorer tested; all 144 historical R2/R3 sample runs rescored
in a separate directory. Benchmark integration patch prepared and checked, but
not applied: that sibling repository is outside this session's writable roots.
No caller reruns, source release pin, old-report deletion, or dashboard replacement.

## Matching contract

New explicit policy: `--coordinate-policy tsd-interval` (version 1).
The scorer's default remains `legacy` for reproducibility of historical runs.

For simulator anchor p and positive truth TSD length L, the matching target
is the inclusive reference interval **[p-L+1, p]**. Distance is zero inside
that interval, otherwise distance to its closest endpoint. A zero-length TSD
uses [p,p]. A missing length is inferred only from literal A/C/G/T/N sequence;
unknown/sentinel values do not invent an interval. Invalid negative lengths
or lengths extending before chromosome position 1 fail explicitly.

This is an **interval-aware detection contract**, not a claim that every point
inside a TSD is the same exact nucleotide coordinate. It accommodates left-TSD
boundary calls and anchor-like one-sided/unknown-TSD calls without making a
wrong predicted TSD sequence automatically count as a missed insertion.
Exact TSD sequence and genotype remain separate accuracy measures. Matching
does not use the predicted TSD string to move a call toward truth.

The 10-bp tolerance is unchanged *outside* the interval; no global enlarged
window was introduced. Effective target width naturally depends on true TSD
length. Family normalization and one-to-one assignment remain unchanged.
Candidates rank by interval distance, then original anchor distance, then
input order. Truth events retain the original greedy coordinate order; this
does not introduce a global maximum-cardinality matching algorithm. Strand is
not a new detection requirement, matching the historical scorer.

For auditability, original truth `position` and `call_position` are retained.
New match fields give `truth_match_start`, `truth_match_end`,
`coordinate_policy`, and `anchor_distance_bp`. Under the new policy,
`distance_bp` means distance to the truth interval, **not** exact-base accuracy.

Synthetic tests cover both interval ends and both tolerance limits, long-TSD
rescue, wrong family, unchanged incorrect TSD/genotype labels, duplicate calls,
one call near two truth events, zero-TSD sentinels, missing lengths, invalid
lengths, and negative windows. The six inherited scorer tests plus 18 policy
tests all pass: **24 passed, no skips or expected failures**.

## Reproduction guard and provenance

Before changing the policy, the runner reproduced each sample's historical
`correctness.tsv`, `matches.tsv`, `false_positive_calls.tsv`, and `precision.tsv`
as row multisets using legacy mode. All 144 samples passed. Thus these changes
are due to matching policy, not new alignments, dropped samples, changed calls,
or a hidden source update. Calls and truth were read from the original run.

The isolated output is `results/coordinate-rescore/2026-09-28/`. It contains
per-sample and pooled report tables, unchanged normalized call copies, copied
historical resource measurements, changed-event and sample-change tables,
a frozen scorer copy, SHA-256 input/source provenance, and a completion marker.
The runner refuses any existing output directory.

These are the **historical full-benchmark calls**. They do not validate the
new breakpoint-pairing candidate or reflect its predicted sensitivity cost.

## Corrected results

R3 means BLAT / bwa aln. Counts pool 9 mPing, 9 riceTElib, and 54 divergence
samples per caller. Metrics are pooled precision/recall/F1, not the dashboard's
mixed aggregation. Truth opportunities are 4,500 / 4,500 / 27,000; loci repeat
across conditions, and zero-divergence overlaps canonical riceTElib.

| Panel | TP, R2 → R3 | FP, R2 → R3 | FN, R2 → R3 | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| mPing | 2820 → 2826 | 0 → 0 | 1680 → 1674 | 100% → 100% | 62.67% → 62.80% | .7705 → .7715 |
| riceTElib | 2494 → 2513 | 18 → 27 | 2006 → 1987 | 99.28% → 98.94% | 55.42% → 55.84% | .7114 → .7139 |
| Divergence | 5888 → 5916 | 77 → 132 | 21112 → 21084 | 98.71% → 97.82% | 21.81% → 21.91% | .3572 → .3580 |

This removes 409 previously counted FP/FN pairs per caller on canonical
riceTElib, and 960 R2 / 954 R3 pairs on divergence. One mPing R2 call is also
rescued. The earlier exact-left-boundary diagnostic counted a narrower subset;
the new contract also handles permitted offsets from that boundary.
There are **no previously matched truth events lost to rescoring and no
changed matched call/status assignments**: every changed-event row is 0→1.

The small relative precision deficit remains and is clearer after removing
the scoring artifact: R3 has nine extra canonical FPs and 55 extra divergence
FPs, versus 19 and 28 additional TPs, respectively. In particular, a small
pooled F1 improvement is not evidence of uniform noninferiority.

### Event-level differences now visible

| Panel | Shared TP | R2-only TP | R3-only TP | R3 sample F1 wins / ties / losses |
|---|---:|---:|---:|---:|
| mPing | 2818 | 2 | 8 | 4 / 3 / 2 |
| riceTElib | 2493 | 1 | 20 | 7 / 1 / 1 |
| Divergence | 5870 | 18 | 46 | 15 / 18 / 21 |

The corrected divergence panel exposes **18 R2-only detections**, not nine.
Previous lost-event lists and family-tie explanations covered only the old
scored subset and must not be treated as a complete current error audit.

## Somatic results after harmonization

These distinguish finding a true somatic insertion from labeling it somatic.

| Panel | Somatic truth | Detected, R2 → R3 | Detection recall | Correctly labeled somatic | Label accuracy among detected somatic events |
|---|---:|---:|---:|---:|---:|
| mPing | 2700 | 1191 → 1194 | 44.11% → 44.22% | 249 → 250 | 20.91% → 20.94% |
| riceTElib | 2700 | 1033 → 1039 | 38.26% → 38.48% | 181 → 180 | 17.52% → 17.32% |
| Divergence | 16200 | 2135 → 2145 | 13.18% → 13.24% | 526 → 525 | 24.64% → 24.48% |

Among *all predictions labeled somatic*, including unmatched calls, corrected
somatic-label precision is 100% → 100%, 98.91% → 97.30%, and
95.64% → 93.09%, respectively. Canonical unmatched somatic-labeled calls fall
to 1 for R2 and 4 for R3; divergence has 12 and 27. The old approximate
82%/79% values for R3 were heavily affected by coordinate mismatch.

However, correct end-to-end somatic-label recall remains low: canonical
181/2700 (6.70%) for R2 versus 180/2700 (6.67%) for R3. Most detected true
somatic events still receive other labels. The shared classifier limitation
has not been fixed by rescoring, and no genotype rule was modified.

## TE-group interpretation changes

Canonical LINE recall is now **299/450 = 66.44% for both callers**, not 15.11%.
SINE recall is **280/450 = 62.22% for R2 and 279/450 = 62.00% for R3**,
not approximately 22%. The earlier apparent severe LINE/SINE discovery deficit
was largely a matching artifact.

Helitron remains 7/450 versus 13/450: its zero-TSD truth is not shifted or
given an artificial duplication. In divergence, corrected LINE detections
are 673 versus 667, SINE 673 versus 673, and Helitron 78 versus 87.
This strengthens the need to inspect the corrected R2-only LINE cases before
claiming users sacrifice no sensitivity.

Coverage, divergence, cellular-fraction, genotype confusion, and resource
tables have been regenerated in
`results/release-audit/2026-09-28-coordinate-corrected/`. Resource values are
copied historical measurements, not new timings of either caller. The known
memory/runtime caveats are unchanged by rescoring.

## Benchmark integration and dashboard status

Patch: `docs/2026-09-28-coordinate-scoring.patch`.

It adds the scorer option/tests, exports `SCORING_COORDINATE_POLICY` from
configuration, passes it in `run_benchmark_array.sh`, and opts the full-aligner
configuration into `tsd-interval`. That configuration's report root becomes
`reports.coordinate-v1`, keeping historical `reports` separate. Other configs
without the policy retain legacy behavior. No SLURM resource or directory
directives were modified. `git apply --check` succeeds against the current
benchmark worktree, preserving its unrelated changes.

From the benchmark repository, install with:

```bash
git apply --check ../../RelocaTE3_jason/RelocaTE3/docs/2026-09-28-coordinate-scoring.patch
git apply ../../RelocaTE3_jason/RelocaTE3/docs/2026-09-28-coordinate-scoring.patch
```

This patch does not move or overwrite existing reports. Applying it alone
does **not** publish this two-caller rescore to the dashboard. The isolated
suite intentionally includes only R2 and R3 BLAT/bwa aln, not all six callers,
and is not presented as a complete dashboard-ready suite. A later dashboard
publication must rescore all displayed callers under the same policy and
preserve dataset metadata; do not mix old and new scores.

The benchmark patch has not been applied from this session because write
access is restricted to the RelocaTE3 workspace. No permission bypass was used.

## Reproduction and checks

From the RelocaTE3 root (existing outputs are protected; use a new path to rerun):

```bash
.pixi/envs/default/bin/python -m pytest -q validation/coordinate_scoring/tests
python scripts/rescore_coordinate_contract.py \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --outdir results/coordinate-rescore/2026-09-28
python scripts/release_benchmark_audit.py \
  --benchmark results/coordinate-rescore/2026-09-28 \
  --outdir results/release-audit/2026-09-28-coordinate-corrected
```

Validation also included shell syntax, patched configuration parsing and
legacy fallback, patch applicability, and `git diff --check`. An initial
configuration assertion assumed unquoted shell values; it was corrected to
parse the actual quoted exports, and both policy branches passed. No production
analysis failed. Work was limited to small report/call tables and lightweight
tests, without BAM scans, alignment, or substantial login-node computation.

## Release decision and next action

The corrected contract substantially improves the interpretation of both
tools, but does not justify a stable “no performance sacrifice” pin. R3 still
has a precision cost, some R2-only true detections, shared somatic-labeling
limitations, and higher memory use. These are rescored historical calls; the
pairing-correction candidate remains unvalidated across all samples.

Next analysis priority: inspect the **18 corrected divergence R2-only
detections and residual R3 false positives**, then assess the frozen pairing
replay under this same scoring contract. Keep benchmark integration, candidate
validation, and clean-artifact release checks distinct.
