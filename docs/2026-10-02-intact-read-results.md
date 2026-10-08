# Intact-read experiment: precision gain with a sensitivity cost

Recorded: 2026-10-02 23:39 PDT, America/Los_Angeles.
Project: RelocaTE3, `main`, base `942a5ea`, after merged PR #49.
Status: all 72 tasks of array 29354331 completed. Candidate not adopted;
production retains the stabilized behavior. No commit, push, or release pin.

## Verified outcome

All completion markers and successful logs are present; frozen hashes verify.
Every task reproduces the stabilized baseline before testing the added veto.
The selected full-read records contain no secondary/supplementary alignments,
so their handling does not explain this comparison.

The added intact-read veto removes **44 false positives and five true
detections**, with no newly added normalized calls. It removes six canonical
and 38 diverged riceTElib FPs without losing a truth detection in either panel.
The five losses are all mPing, where there were no FPs to remove. Consequently
the candidate fails the intended sensitivity-preservation criterion despite
remaining slightly above R2 in pooled mPing recall.

| Panel | Caller | TP | FP | FN | Precision | Recall |
|---|---|---:|---:|---:|---:|---:|
| mPing | R2 | 2820 | 0 | 1680 | 100% | 62.667% |
| mPing | Stabilized R3 | 2826 | 0 | 1674 | 100% | 62.800% |
| mPing | Intact candidate | 2821 | 0 | 1679 | 100% | 62.689% |
| riceTElib | R2 | 2494 | 18 | 2006 | 99.283% | 55.422% |
| riceTElib | Stabilized R3 | 2513 | 27 | 1987 | 98.937% | 55.844% |
| riceTElib | Intact candidate | 2513 | 21 | 1987 | 99.171% | 55.844% |
| Diverged riceTElib | R2 | 5888 | 77 | 21112 | 98.709% | 21.807% |
| Diverged riceTElib | Stabilized R3 | 5919 | 129 | 21081 | 97.867% | 21.922% |
| Diverged riceTElib | Intact candidate | 5919 | 91 | 21081 | 98.486% | 21.922% |

The net divergence FP gap would fall from 52 to 14, and the canonical gap
from nine to three. These are condition-level observations; the canonical
and zero-divergence panels overlap and should not be treated as independent
biological datasets. Scores use the same `tsd-interval` policy and 10-bp
tolerance for all variants.

## Lost detections

| Sample | Truth event | Class | Cellular fraction | R2 detects it? |
|---|---|---|---:|---|
| cov30x_rep2 | TE000316 | Somatic | 0.2 | No |
| cov5x_rep1 | TE000125 | Heterozygous | 1.0 | No |
| cov5x_rep2 | TE000157 | Heterozygous | 1.0 | No |
| cov5x_rep2 | TE000434 | Somatic | 0.4 | No |
| cov5x_rep2 | TE000452 | Somatic | 0.4 | No |

Thus the candidate erases genuine R3 advantages; matching R2's rejection does
not make these losses desirable. Four occur at 5x and three are somatic.
Aggregate mPing somatic detection decreases from 1194 to 1191 and heterozygous
detection from 773 to 771. All event-level gains/losses are retained in the
analysis tables rather than inferred from net TP differences.

## What the read checks support

The prior Os3912 fixture established the useful case: an intact full read
aligns nearby, but its short trimmed flank falsely maps at a reference TE edge.
The experimental veto also accepts intact alignments anywhere in the genome,
following R2's intact-flag logic. That broader scope can oppose a genuine
insertion whose original read maps to another TE copy.

For three losses, bounded indexed inspection of the existing +/-500 bp window
found **zero intact alignments for their candidate read ends locally**, despite
the global intact flags triggering rejection:

| Event | Candidate read ends | Locally observed ends | Local intact ends | Global intact left/right |
|---|---:|---:|---:|---|
| TE000316 | 4 | 1 | 0 | 1 / 1 |
| TE000125 | 4 | 1 | 0 | 1 / 1 |
| TE000434 | 1 | 0 | 0 | 0 / 1 |

The locally observed alignments for the first two are clipped (`39M111S` and
`43S105M2S` respectively). The opposing intact evidence is outside the queried
window. Its precise remote placement has not been extracted, so this is
evidence for nonlocal vetoing, not proof of which repeat copy attracted it.
The other two losses have left-only raw candidates whose starts differ from
their final `supporting_junction` coordinates; their veto records are present,
but they were not included in the three direct-coordinate inspections above.

## Reproduction and checks

```bash
.pixi/envs/default/bin/python scripts/summarize_intact_replay.py \
  --replay results/intact-fullreads-replay/2026-10-02 \
  --stabilization results/family-replay/2026-10-02-stabilization \
  --output results/intact-fullreads-analysis/2026-10-02
```

Already run; choose a new output directory to repeat. Outputs include pooled
and per-sample metrics, biological-class/fraction/TE-group/divergence strata,
paired gain/loss counts against both baseline and R2, and input hashes. The
summarizer checks completion, baseline reproduction, frozen source hashes,
task identity, unique event IDs, common truth sets, and TP count consistency.
Normalized-call multiset comparison independently confirms exactly 49 removals
and no additions. This turn used small-table analysis and bounded BAM queries,
not a new caller run on the login node. A lint-only variable rename followed
the analysis; the provenance retains the script hash actually executed.

## Decision and next action

Do not enable the global intact-read veto in production. Do not special-case
mPing or select thresholds from truth labels. The next controlled hypothesis
is to require intact-read evidence local to the candidate, retaining the
existing spanning-read test and the mate-specific lookup. This could preserve
the demonstrated nearby-reference-copy correction while avoiding remote-copy
vetoes, but its full-panel accuracy has not been tested. It should undergo
the same paired replay before adoption, including all five losses as explicit
regression targets. No new job has been prepared or submitted in this analysis
turn. Release readiness remains unproven.
