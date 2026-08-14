# Recover the insertions RelocaTE2 finds and RelocaTE3 misses on Chr3

Closes the last detection gap against RelocaTE2 on the rice Chr3 2 Mb fixture:
RelocaTE3 now matches it **call-for-call** — 196 calls, 196 of 200, precision
1.000, F1 0.990, with zero calls unique to either caller.

> [!IMPORTANT]
> **Do not merge this PR on its own.** On the riceTElib panel these two commits
> regress F1 by 0.036 / 0.053 / 0.061 at 5x / 15x / 30x — 9 of 9 samples worse.
> The fix is #`<PR 2>` (`fix/tsd-sequence-accuracy`), which stacks on this
> branch. Merge that immediately after, or merge it alone since it contains
> both. Detail in "Known regression, fixed downstream" below.

## Background

Three truth sites were found by RelocaTE2 and missed entirely by RelocaTE3. They
turned out to have two unrelated causes, and neither was the "second-junction
recovery" the earlier analysis predicted.

## `db77dca` — don't merge two adjacent insertions into one impossible call

`_pair_breakpoints` groups breakpoints within `SUBCLUSTER_GAP` (25 bp) and keeps
the dominant left and right of each group. When two real insertions sit close
together their *inner* breakpoints fall inside that gap, so one group holds
both, and the dominant pair becomes insertion A's left breakpoint with insertion
B's right. One call is emitted, at neither site, and both true insertions are
lost.

The bad pairing is geometrically impossible, which is what makes it detectable
with no extra evidence: a left-junction read marks the TSD's *right* edge and a
right-junction read its *left* (see `_make_insertion`), so a genuine pair always
has `left >= right`.

```
truth Chr3:155988   RT2: 155951..155962  and  155984..155988
                    R3 : one call at 155962      (L=155962, R=155984 -> -21)
truth Chr3:672724   RT2: 672693..672695  and  672720..672724
                    R3 : one call at 672695      (L=672695, R=672720 -> -24)
```

Both recovered, with their TSDs resolved (`ACTCT`, `CCATATAT`).

## `8c44dce` — keep two-sided calls that abut a reference TE boundary

RelocaTE3 discarded any call whose edge coincided with a stored reference TE
edge, regardless of support. RelocaTE2 applies that rule only to calls with an
**empty junction side** (`clean_false_positive.py:82`, `Right == 0 or Left == 0`):
reads running off an intact reference copy's edge mimic a novel junction, but a
call with junction reads on both sides is real evidence — transposons do insert
next to transposons.

It cost a genuine call. The mPing insertion at Chr3:257446..257448 has 3 left and
6 right junction reads and resolves TSD `ACG` — the same call RelocaTE2 makes —
but a TEOS1 copy ends at 257444 and the annotation loader stores a window of end
positions around it, so 257445 matched and the call was dropped. It did not even
reach the raw tier; the cluster fell through to the supporting-reads path.

Extracted as `_excluded_by_reference_edge`, gated on junction support, matching
the rule already used by the tier-level boundary filter (`_has_empty_side`).

## Chr3 result

| step | calls | recall | precision | F1 |
|---|---|---|---|---|
| before | 193 | 193/200 | 1.000 | 0.982 |
| + adjacent-insertion split | 195 | 195/200 | 1.000 | 0.987 |
| + gated edge exclusion | **196** | **196/200** | **1.000** | **0.990** |
| RelocaTE2 | 196 | 196/200 | 1.000 | 0.990 |

Zero calls unique to either caller; all 196 shared.

## Known regression, fixed downstream

Validating on riceTElib (9 samples, 3 coverages) showed these commits regress
**9 of 9 samples**:

| coverage | baseline F1 | this branch | RelocaTE2 |
|---|---|---|---|
| 5x | 0.444 | 0.409 | 0.438 |
| 15x | 0.647 | 0.593 | 0.620 |
| 30x | 0.731 | 0.670 | 0.696 |

Cause is `db77dca`: splitting on *any* `left < right` is too aggressive. The
geometry argument holds, but at low coverage a genuine single insertion shows
1–3 bp of jitter in that direction, and splitting those yields two one-sided
half-calls that both miss. The Chr3 fixture cannot expose this — its two real
cases are separated by 22 and 25 bp, far above jitter.

`5f8aa81` in the stacked PR requires separation beyond `SPLIT_MIN_SEPARATION`
(10 bp), which restores riceTElib to at-or-above baseline and above RelocaTE2 at
every coverage while keeping every Chr3 recovery.

A bisect on the same benchmark cleared PR #44: it scored identically to the
pre-merge baseline (183 calls, F1 0.448) on the probe sample, so `main` was never
affected.

## Verification

- Chr3: 196/200, precision 1.000, F1 0.990 — identical to RelocaTE2.
- Suite: 177 passed / 0 skipped in the `blat` env; 175 / 2 (blat-gated) in the
  default env.
- New: `tests/breakpoint_pairing_test.py` — pairing geometry, the split, the
  edge-exclusion gating, each citing the RelocaTE2 source line it mirrors.

Full write-up: `notes/2026-08-12-relocate2-chr3-baseline.md`.
