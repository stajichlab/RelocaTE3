# `--require-both-junctions`: what it does and when to use it

2026-08-13

## Purpose

Decide whether RelocaTE3 reports a non-reference insertion supported by junction
reads on only **one** side of the breakpoint.

- `--require-both-junctions` (default): report a call only when both the 5' and
  3' junctions are observed.
- `--no-require-both-junctions`: also report one-sided calls.

## Current status

Default flipped to **on** in `288cdb1`. RelocaTE2 has no equivalent switch; it
reports one-sided calls, but almost never a wrong one.

## What the measurements say

Two-sided calls are essentially free of false positives in both callers. The
entire behavioural difference lives in the one-sided population.

mPing 30x, `cov30x_rep1`, scored against per-sample truth (window 10bp):

| population | RelocaTE3 | RelocaTE2 |
|---|---|---|
| two-sided | 394 calls, 100% correct | 407 calls, 100% correct |
| one-sided | 78 calls, **45%** correct | 53 calls, **100%** correct |

So turning the flag on costs RelocaTE3 35 true calls and removes 43 false ones.

Across datasets the trade-off inverts:

| dataset | flag on | flag off |
|---|---|---|
| riceTElib (many TE families) | F1 **0.722** | F1 0.144 |
| mPing (single element) | slightly lower recall, precision ~1.00 | slightly higher F1 |

## Guidance for users

**Leave it on** (the default) for any library with more than one TE family. On
riceTElib it is the difference between a usable and an unusable call set: with
it off, one-sided calls from paralogous and nested elements swamp the output and
F1 collapses from 0.722 to 0.144.

**Consider turning it off** only when all of the following hold:

- a single, well-characterised element (an mPing-style screen),
- recall matters more than precision, and
- downstream work will validate calls anyway.

It is best described as a **precision gate**, not a general "high-confidence
mode": it does not rank or score calls, it removes exactly one class of
evidence. On a single-element library the class is mostly signal; on a
multi-family library it is mostly noise.

## Why RelocaTE3's one-sided calls are worse

Open. Five RelocaTE2 mechanisms were ported or tested and **none** closes the
gap:

1. `supporting_junction` classification rule — tested, not the cause.
2. MAPQ / support thresholds — ported (`8e45e61`), not the cause.
3. Reference-copy proximity — tested, not the cause.
4. Low-quality-read call validation (`relocaTE_insertionFinder.py:226-241`) —
   ported (`8e45e61`), not the cause.
5. False-junction filter (`:212-221`) — ported and fixed (`6276472`); with the
   lookup verified live it drops **0 of 472** calls. A junction read carries TE
   sequence and is soft-clipped at the breakpoint, so it cannot map through it;
   the filter targets a failure mode this data does not produce.

Since no downstream filter explains it, RelocaTE2 must not *create* the 43 bad
candidates in the first place. The cause is upstream, in candidate
construction rather than candidate filtering. Deferred to v1.1.

## Remaining known parity gaps

- `ST` is hardcoded to `ST:0` in the RelocaTE3 writer.
- `--mismatch_junction` is collapsed into `--mismatch`.
- `pipeline.run_sample` and `InsertionFinder` are separate code paths; five
  behavioural differences have hidden in that split so far. Collapsing them is
  the highest-value structural cleanup.
