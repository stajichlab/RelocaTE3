# `--require-both-junctions`: what it does and when to use it

2026-08-13

## Purpose

Decide whether RelocaTE3 reports a non-reference insertion supported by junction
reads on only **one** side of the breakpoint.

- `--require-both-junctions` (default): report a call only when both the 5' and
  3' junctions are observed.
- `--no-require-both-junctions`: also report one-sided calls.

## Current status

Default flipped to **on** in `288cdb1`.

> **Correction (2026-08-14).** This document previously said "RelocaTE2 has no
> equivalent switch; it reports one-sided calls". **That is wrong**, and the
> error propagated into the benchmark harness, where it justified running
> RelocaTE3 in the permissive mode for a supposedly like-for-like comparison.
>
> RelocaTE2's *final* call set is two-sided-only, plus one narrow exception.
> The claim came from reading `relocaTE_insertionFinder.py:365`
> (`l_count >= 1 OR r_count >= 1`) in isolation, but that line only admits a row
> into the intermediate `all_nonref_insert.txt`. The branch at `:373` writes a
> real TSD **only** when both sides carry junction reads; every one-sided case
> instead gets a sentinel word in the TSD column (`supporting_junction`,
> `singleton`, `insufficient_data`). Those sentinels are then filtered twice:
>
> - `characterizer.pl:91` —
>   `if ( ($left_count >= 1 and $right_count >= 1) or $TSD eq 'supporting_junction' )`
> - `clean_false_positive.py:99,107` —
>   `grep -v "singleton|insufficient_data|supporting_reads"`
>
> Counted on riceTElib `cov30x_rep1`: RelocaTE2's intermediate table holds 421
> rows (357 real TSD + 36 `insufficient_data` + 25 `singleton` + 3
> `supporting_junction`) and its characterized output holds 360 = 357 + 3.
> mPing `cov30x_rep1`: 407 + 32 + 16 + 5 → 412 = 407 + 5.
>
> So `--require-both-junctions` is **RelocaTE2's own behaviour**, not a
> RelocaTE3 policy layered on top of it. As of 2026-08-28, RelocaTE3 also
> retains RelocaTE2's narrow `supporting_junction` exception after correcting
> the supporting-read population described below.

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
| riceTElib (many TE families), blat + bwa aln | F1 **0.715** | F1 0.143 |
| riceTElib, blat + bwa mem | F1 0.722 | F1 0.722 (no collapse) |
| mPing (single element) | slightly lower recall, precision ~1.00 | slightly higher F1 |

> An earlier version of this table compared riceTElib "flag on 0.722" against
> "flag off 0.144" — those two numbers came from **different genome aligners**
> (bwa mem and bwa aln respectively), so they overstated the flag's effect.
> The collapse is specific to `bwa aln` × permissive: `bwa aln` has no seed
> floor and places the short junction flanks bwa mem drops, which is a recall
> win and, without the gate, a flood of one-sided calls. bwa mem permissive
> does not collapse (precision 0.85).

## Guidance for users

**Leave it on** (the default). It is what RelocaTE2 does, and on a multi-family
library with the default `bwa aln` genome aligner it is the difference between a
usable and an unusable intermediate call set.

**Consider turning it off** only when all of the following hold:

- a single, well-characterised element (an mPing-style screen),
- recall matters more than precision, and
- downstream work will validate calls anyway.

Note that since the characterize gate was corrected (below), turning the flag
off no longer floods the *characterized* output — those calls are dropped at
characterize instead. The flag now mainly determines what appears in the
intermediate `.txt` and the `.raw` tier.

## The gate that enforces this

As of 2026-08-14 the flag is no longer the only thing standing between a
one-sided call and the output. `characterize.py` now implements
`characterizer.pl:91` faithfully:

```python
junction_supported = left_count >= 1 and right_count >= 1
if not (junction_supported or tsd == "supporting_junction"):
    continue
```

It previously read `left_count >= 1 or right_count >= 1` with a
`total_count > 1` guard, which admitted every one-sided multi-read cluster.
That single operator was the whole riceTElib precision collapse. Re-scoring the
existing `cov30x_rep1` permissive run through both gates:

| gate | calls | TP | FP | precision | recall | F1 |
|---|---|---|---|---|---|---|
| old (`or`, `total > 1`) | 5264 | 348 | 4916 | 0.066 | 0.696 | 0.121 |
| **new (`and`) = RelocaTE2** | **391** | **316** | **75** | **0.808** | **0.632** | **0.709** |
| RelocaTE2, measured | 360 | 298 | 62 | 0.828 | 0.596 | 0.693 |

`--require-both-junctions` was masking that bug upstream rather than fixing it.
With the gate corrected the flag is close to redundant: it now controls only
what reaches the intermediate `.txt` and the `.raw` tier, not what is
characterized.

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

### Traced 2026-08-15 — the cause is upstream, and it is two things

Taking the 20 mPing `cov30x_rep1` sites RelocaTE2 reports and RelocaTE3 misses,
and asking what RelocaTE3 has within 200 bp:

- **20 of 20 have a RelocaTE3 candidate; never both sides as separate calls.**
  So this is *not* a breakpoint-pairing problem — RelocaTE3 genuinely lacks the
  second junction's reads. (RelocaTE2 pairs left/right junctions up to 100 bp
  apart in `TSD_from_read_depth`; RelocaTE3's `SUBCLUSTER_GAP` is 25, but that
  difference is not what costs these calls.)
- **12 of 20 are sites where RelocaTE2's weaker junction side has exactly one
  read** and RelocaTE3 has zero there. 5 more are RelocaTE2 `supporting_junction`
  calls (see the gap above).

Tracing RelocaTE2's 59 junction reads at those sites through RelocaTE3's own
intermediates:

| stage | reads |
|---|---|
| present in RelocaTE3's trim output | 40 of 59 |
| present in RelocaTE3's genome BAM | 40 of 59 |
| **lost at trim** (never identified as TE-containing) | **19** |
| lost at genome alignment | 0 |

Two distinct causes, then:

1. **A MAPQ admission gate RelocaTE2 does not have.** `min_mapq` defaulted to 1,
   discarding MAPQ-0 reads outright, where RelocaTE2 admits them and only marks
   them low quality. **Fixed 2026-08-15** — default is now 0; measured +4 true
   calls at `cov30x_rep1` with no new false positives (F1 0.879 → 0.884).
2. **Trim sensitivity** — 19 of 59 reads RelocaTE2 identifies as TE-containing
   never appear in RelocaTE3's trim output at all. This is the larger remaining
   share and is untouched. Nothing is lost at genome alignment.

## Remaining known parity gaps

- **Resolved 2026-08-28:** RelocaTE2 adds only unpaired, non-junction BAM
  records to `teSupportingReads`; RelocaTE3 added every non-junction record.
  Applying the legacy predicate reduced eligible one-sided
  `supporting_junction` candidates from 18 to 6 in riceTElib `cov5x_rep1`,
  while retaining all 185 two-sided calls. The six comprise all five truth
  events RelocaTE2 detected alone plus one additional false positive.
  RelocaTE3 now emits the sentinel class and coordinates, and `ST:` is the real
  sum of `SR:` and `SL:` rather than a hardcoded zero.
- Secondary TE-family labels inferred only from supporting reads are still not
  appended to the primary junction family. This affects exact label strings,
  not the primary family used by the benchmark scorer.
- `--mismatch_junction` is collapsed into `--mismatch`.
- `pipeline.run_sample` and `InsertionFinder` are separate code paths; five
  behavioural differences have hidden in that split so far. Collapsing them is
  the highest-value structural cleanup.
