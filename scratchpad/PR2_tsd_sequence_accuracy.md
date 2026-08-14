# Fix the riceTElib regression, and stop inventing TSDs

Stacks on `fix/chr3-missed-detections` (#`<PR 1>`). Two commits: the first
repairs a regression that branch introduced, the second fixes a long-standing
defect in TSD reporting.

After both, RelocaTE3 beats RelocaTE2 on **detection at every coverage** and on
**TSD accuracy wherever a TSD exists** — while keeping the Chr3 call-for-call
parity PR 1 achieved.

> [!NOTE]
> If PR 1 has not merged yet, merging this branch alone is safe and preferred —
> it contains all four commits. PR 1 on its own regresses riceTElib.

## `5f8aa81` — require real separation before splitting a breakpoint pair

`db77dca` (PR 1) split a breakpoint group whenever its dominant breakpoints were
ordered `left < right`, on the grounds that the arrangement is geometrically
impossible for a single insertion. The reasoning holds; the threshold did not. At
low coverage a genuine single insertion routinely shows 1–3 bp of jitter in that
direction, and splitting those produced two one-sided half-calls that both then
missed.

Caught by riceTElib, where it regressed **9 of 9 samples**:

| coverage | baseline | PR 1 alone | **with this fix** | RelocaTE2 |
|---|---|---|---|---|
| 5x | 0.444 | 0.409 | **0.448** | 0.438 |
| 15x | 0.647 | 0.593 | **0.643** | 0.620 |
| 30x | 0.731 | 0.670 | **0.722** | 0.696 |

A split now requires separation beyond `SPLIT_MIN_SEPARATION` (10 bp): clear of
jitter, clear of the genuine Chr3 cases at 22 and 25 bp, and below
`SUBCLUSTER_GAP`, above which grouping already separates breakpoints. Chr3 is
unchanged at 196/200, precision 1.000, F1 0.990 — the threshold costs none of the
recoveries `db77dca` was written for.

The mechanism confirms the diagnosis: calls went 378 → 362 (PR 1 destroying real
insertions) → 401 (fix restoring them, plus the genuine wide-separation splits).

## `38d2f1f` — report UNK rather than inventing an implausibly long TSD

RelocaTE3 reported the genomic sequence across the inferred TSD span whatever its
length, so when the depth estimator over-reached it emitted 87 bp and 100 bp
"TSDs". Target-site duplications are short; the longest in the riceTElib truth
set is 20 bp.

Splitting truth sites by whether they actually have a TSD makes it precise:

| site class | R3 detects | R3 reports a TSD | longest true TSD |
|---|---|---|---|
| has a real TSD | 1972 | all ≤ 12 bp | 20 bp |
| TSD-less (`TSD=NONE`) | 255 | 140, of which **114 > 20 bp** | — |

Every over-long TSD is on an element with no TSD at all, and no correct call comes
within 8 bp of the cap, so `MAX_PLAUSIBLE_TSD = 20` removes only fabricated
sequence. **The insertion is still called; only the TSD field changes to UNK.**

Verified detection-neutral on all 9 samples — F1 identical to three decimals at
every coverage, real-TSD exactness identical, TSD-less handling 31–53% → 87–95%.

### This also corrects an earlier misreading of the same data

RelocaTE3's TSD accuracy is **better** than RelocaTE2's, not worse. An earlier
analysis reported 86.7% vs 95.9% by scoring the 10% of truth sites that are
TSD-less as TSD errors. On sites that genuinely have a TSD:

| caller | real-TSD exact | TSD-less sites handled |
|---|---|---|
| RelocaTE2 | 2010/2085 (96.4%) | 0/10 (0%) |
| **RelocaTE3** | 1930/1972 (**97.9%**) | 229/255 (**90%**) |

RelocaTE3 detects 255 TSD-less insertions to RelocaTE2's 10 — a large recall win
on that class — which is exactly why RelocaTE3 has to handle them and RelocaTE2
never had to.

## Final position (all four commits)

riceTElib, 9 samples, blat/bwa-aln, ±10 bp matching:

| coverage | RelocaTE2 R/P/F1 | **RelocaTE3** R/P/F1 |
|---|---|---|
| 5x | 0.300 / 0.815 / 0.439 | **0.308 / 0.825 / 0.448** |
| 15x | 0.495 / 0.839 / 0.622 | **0.526 / 0.827 / 0.643** |
| 30x | 0.602 / 0.839 / 0.701 | **0.651 / 0.812 / 0.722** |

Chr3 2 Mb: 196/200, precision 1.000, F1 0.990 — identical to RelocaTE2.

## Verification

- Suite: 186 passed / 2 skipped (blat-gated).
- New: `tests/tsd_plausibility_test.py` — the cap, the genome-fallback path, and
  that the cap sits above the longest real TSD observed.
- Every number above is from a real `run-all` on the benchmark data, scored the
  same way for both callers.

## Known and not addressed

- **30x precision trails RelocaTE2** (0.812 vs 0.839) while F1 leads. Traceable
  to `8c44dce` keeping ~20 extra calls per sample; recovering it means also
  gating `_excluded_by_reference_edge` on minimum junction support.
- **mping dataset**: RelocaTE3 trails RelocaTE2 at 5x (F1 0.546 vs 0.578, both at
  precision 1.000). Verified **pre-existing** — a July 28 run of RelocaTE3 scores
  byte-identically, so no commit here caused or worsened it. mPing is the single
  element RelocaTE2 was tuned around; riceTElib is a 500-family panel where
  RelocaTE3's generality wins.
- 3 residual `UNK` calls on Chr3.
