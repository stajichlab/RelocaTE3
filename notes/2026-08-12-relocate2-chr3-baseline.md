# RelocaTE2 Chr3 2 Mb baseline — measured, not quoted

**Date:** 2026-08-12

## Purpose

The figure "RelocaTE2 recovers 196/200 on this dataset" appears in
`tests/acceptance_test.py:4`, `tests/acceptance_test.py:69`, `plans/PLAN.md:201`
and `README.md`. Every instance traced back to the same unsourced claim, which
originates in RelocaTE2's own README — it had never been reproduced here. Since
"matches RelocaTE2" is the v1 release claim, it needed measuring.

## Result: the claim is correct

RelocaTE2 run on its own `test_data` with its own defaults (its README passes no
tuning flags), scored with `bedtools window -w 10` semantics against
`MSU7.Chr3_2M.ALL.gff`:

| RelocaTE2 tier | calls | recall | precision | F1 |
|---|---|---|---|---|
| `.raw` / `.all` | 207 | 196/200 | 0.947 | 0.963 |
| **headline `.gff`** | **196** | **196/200 (0.980)** | **1.000** | **0.990** |
| `.high_conf` | 196 | 196/200 | 1.000 | 0.990 |

196 calls, 196 window matches — exactly RelocaTE2's README. Runtime ~20 s.

## Head to head, identical inputs

RelocaTE3 (this branch, shipped defaults) on **the same** `test_data` files:

| tier | calls | recall | precision | F1 |
|---|---|---|---|---|
| RelocaTE2 headline | 196 | 196/200 (0.980) | 1.000 | **0.990** |
| RelocaTE3 headline | 199 | 193/200 (0.965) | 0.970 | 0.967 |
| RelocaTE3 `.high_conf` | 193 | 191/200 (0.955) | 0.990 | **0.972** |

**RelocaTE2 is ahead on this dataset**, F1 0.990 vs 0.972, but by less than it
first appeared. The first measurement had RelocaTE3's `.high_conf` at F1 0.939;
that was a bug in *our* tier filter, not a property of the caller — see below.

Note this is the opposite of the riceTElib benchmark, where RelocaTE3 leads at
every coverage (F1 0.444/0.645/0.729 vs 0.438/0.620/0.696 at 5x/15x/30x).
Chr3 2 Mb is a single easy region at 6X; riceTElib is a harder multi-TE panel.
Neither result generalises on its own.

## Correction: our high_conf filter was stricter than RelocaTE2's

RelocaTE2's `high_conf` step greps out the *literal* strings
`Right_junction_reads=1;Left_junction_reads=0` and its mirror
(clean_false_positive.py:108) — **exactly one junction read against zero**, and
nothing else. A one-sided call backed by several junction reads survives.

Our implementation dropped every call with a zero side. On this fixture that
discarded 16 calls RelocaTE2 keeps:

| RelocaTE3 `.high_conf` rule | calls | recall | precision | F1 |
|---|---|---|---|---|
| drop any zero side (wrong) | 177 | 177/200 | 1.000 | 0.939 |
| drop only 1-vs-0 (RelocaTE2) | 193 | 191/200 | 0.990 | **0.972** |

Note also that RelocaTE2's own headline and `.high_conf` are identical here
(196 = 196) simply because no `1 vs 0` calls survived its class filter — not
because it has no one-sided calls. It keeps 4 of them.

## What is actually left

With the filter corrected, **5 truth sites** are found by RelocaTE2 and not by
RelocaTE3 — down from the 19 the buggy filter implied. At all five, RelocaTE3
makes **no call at all**; none is a one-sided call that could be rescued by
recovering a second junction. So "second-junction recovery" was the wrong
diagnosis.

A separate, larger quality gap was TSD inference: 16 of the 191 shared sites had
RelocaTE2 resolving a concrete TSD where RelocaTE3 reported `UNK`. **Fixed in
50fcbc8** — the depth estimator was fed supporting mates instead of junction
reads, so it returned 0 (see the commit). After the fix, `high_conf`:

| | calls | recall | precision | F1 | UNK |
|---|---|---|---|---|---|
| RelocaTE3 before | 193 | 191/200 | 0.990 | 0.972 | 18 |
| **RelocaTE3 after** | 193 | 193/200 | **1.000** | **0.982** | **3** |
| RelocaTE2 | 196 | 196/200 | 1.000 | 0.990 | 0 |

The F1 gap to RelocaTE2 is now 0.008, from 0.018.

## Incidental finding: the vendored RepeatMasker is synthetic

`tests/data/sim_genome/MSU7.Chr3_2M.fa.RepeatMasker.out` is **not** RelocaTE2's
real annotation. RelocaTE2's is 1255 lines of genuine RepeatMasker output;
ours is 121 lines with round SW scores (999/501/502), `DNA/MITE` classes and
sequential IDs — clearly hand-made.

It materially changes results. Same reads, same TE library, same genome:

| RepeatMasker | R3 raw calls | recall | precision |
|---|---|---|---|
| vendored (synthetic) | 215 | 194/200 | 0.902 |
| RelocaTE2's (real) | 199 | 193/200 | **0.970** |

The synthetic annotation **understates RelocaTE3's precision by ~0.07**, because
the reference-TE filtering has far less to work with. Genome, truth GFF, reads
(uncompressed) and TE library sequences are all identical between the two
fixtures — only the annotation differs.

Worth replacing the vendored file with RelocaTE2's real one, or at least saying
in the fixture that it is synthetic.

## Reproducing

    /rhome/nmath020/bigdata/github/github_tools/RelocaTE/rt2_chr3_baseline/run_rt2.sh

Environment comes from the relocate-benchmark RT2 adapter (digest-pinned
BioContainers); RelocaTE2 is a dead py2.7 package and cannot be installed
directly. Outputs and the R3 comparison run live alongside it under
`rt2_chr3_baseline/`.

## Next

- The 5 sites RelocaTE3 misses entirely (not a second-junction problem) —
  now the whole remaining detection gap.
- The 3 residual `UNK` calls, if they matter.
- Replace or label the synthetic RepeatMasker fixture.
- The 196/200 references in the acceptance test / PLAN / README can now cite
  this measurement rather than RelocaTE2's README.
