# RelocaTE2 parity: matching defaults, matching outputs, and a measured baseline

Makes RelocaTE3 a drop-in for RelocaTE2 — same defaults, same output files — and
replaces the long-quoted "RelocaTE2 gets 196/200" with a number we actually
measured. Along the way, two accuracy bugs surfaced and are fixed.

On RelocaTE2's own Chr3 2 Mb test data, identical inputs and identical scoring:

| | calls | recall | precision | F1 |
|---|---|---|---|---|
| RelocaTE2 (headline) | 196 | 196/200 | 1.000 | **0.990** |
| RelocaTE3 before this branch | 215 | 194/200 | 0.902 | 0.935 |
| **RelocaTE3 after** (`high_conf`) | 193 | 193/200 | 1.000 | **0.982** |

Nine commits, each independently reviewable and each carrying its own
measurements. **Reviewing commit-by-commit is much easier than the squashed
diff.**

---

## Part 1 — Parity (5 commits)

### `32dc485` — default to RelocaTE2's aligner pair

`blat` for TE search, `bwa aln` for genome placement, replacing minimap2 for
both. On the riceTElib benchmark that pair is not just RT2-matching but simply
better — matched true calls 455/781/962 at 5x/15x/30x vs 346/632/798 for
minimap2, at higher precision too.

blat can't go in the main pixi environment (its bioconda build hard-pins
zlib 1.2.11, which matplotlib won't solve against, and it's linux-64 only). The
plotting stack — used only by `validation/real_rice/`, never by `src/` — moves to
a `plot` feature, freeing a `blat` environment that carries the core runtime:

```bash
pixi run -e blat relocaTE3 run-all ...
```

`default` still includes plotting, so `pixi run test` and `validate-rice` are
unchanged. A missing blat now raises a message naming both ways out instead of a
bare `FileNotFoundError` from inside the trim step.

### `39f6729` — ship RelocaTE2's parameter defaults

`--mismatch` 0 → 2, and `--tsd` on `find-insertions` goes from *required* to
defaulting `UNK` (RelocaTE2 hardcodes UNK; it never asks the user). A bare
`run-all` with no tuning flags now reproduces a stock RelocaTE2 run.

Two things that look like they should match RelocaTE2 but must not — both
verified in its source, both worth a reviewer's eye:

- **`--require-both-junctions` stays opt-in.** RelocaTE2 keeps a site on
  `l_count >= 1 OR r_count >= 1` (`relocaTE_insertionFinder.py:365`), so it emits
  single-junction calls. Turning this on by default would drop calls RelocaTE2
  reports.
- **`--min-mapq` has no RelocaTE2 counterpart at all.** RelocaTE2 classifies
  reads below MAPQ 29 as low-quality via bwa's `XM`/`X1`/`XO` tags, not a
  threshold.

`tests/relocate2_defaults_test.py` pins every value with the RelocaTE2 source
line cited per assertion, including that the library entry points can't drift
from the CLI.

### `28d3c63` — emit RelocaTE2's tiered call sets, and a GFF at all

RelocaTE2 publishes several call sets, not one, and **its documented 196 is the
filtered file**. RelocaTE3 emitted a single unfiltered table and no GFF at all
from `find-insertions` — so every RelocaTE3-vs-RelocaTE2 comparison anyone had
made was raw-versus-filtered.

Worth reviewing: **the plain `.txt` is deliberately left unfiltered.** RelocaTE2
cleans only its GFF (`relocaTE2.py:704`), concatenates the table unfiltered
(`:703`), and genotypes *that* table (`:707`). An earlier draft filtered the
table too, which would have quietly shrunk what `characterize` sees.

### `bba0492` — RelocaTE2's reference-TE boundary filter, and `.raw`

RelocaTE2's filter is narrower than ours was: it drops a call only when it is
one-sided **and** an endpoint sits within `--distance` (3 bp) of a reference TE
boundary — reads from an intact copy's edge mimicking a novel junction. A
two-sided call at a boundary is kept.

Recorded negative result: on Chr3 this removes nothing (37 one-sided calls, none
within 50 bp of any of 1200 boundaries). Correct to have; not a lever here.

### `1b907d8` — mate-pair-only insertions, `-s/--size`

RelocaTE2 calls sites from supporting mates when no junction read mapped and
files them separately as `all_nonref_supporting.{txt,gff}`.

The gap was bigger than a missing file: `_call_support_only` existed but hangs
off the module-level `find_insertions` that only `run_sample` uses — **unreachable
from `InsertionFinder`**, which is what the CLI and the benchmark run. Every
shipped run emitted zero support-only calls, silently. Fourth time the
`run_sample`/`InsertionFinder` split has hidden a behavioural difference.

These go only to the supporting file, never into `all_nonref_insert` or its
tiers, so they cannot move existing numbers. Second negative result: on Chr3 the
path finds 6 sites and **none is a true insertion**. They're emitted because
RelocaTE2 emits them, and segregated for the same reason.

---

## Part 2 — Baseline and accuracy (4 commits)

### `5d7c2fb` — measure the RelocaTE2 baseline instead of quoting it

"196/200" appeared in four places and every one traced to the same unsourced
claim in RelocaTE2's README. Ran RelocaTE2 itself on its own `test_data`: **the
claim holds** — 196 calls, 196 window matches, precision 1.000, ~20 s.

**Incidental and important:** `tests/data/sim_genome/MSU7.Chr3_2M.fa.RepeatMasker.out`
is **synthetic**. RelocaTE2's real annotation is 1255 lines; ours is 121 with
round SW scores (999/501/502) and sequential IDs. Same reads, TE library and
genome — only the annotation differs:

| RepeatMasker | R3 calls | recall | precision |
|---|---|---|---|
| vendored (synthetic) | 215 | 194/200 | 0.902 |
| RelocaTE2's (real) | 199 | 193/200 | **0.970** |

It understates RelocaTE3's precision by ~0.07. Not changed in this PR —
flagged for a follow-up.

### `a66cff2` — `high_conf` must drop only single-read one-sided calls

A bug in the tier filter added earlier in this same branch. RelocaTE2's
`high_conf` greps out the *literal* strings
`Right_junction_reads=1;Left_junction_reads=0` and its mirror
(`clean_false_positive.py:108`) — exactly one read against zero. Ours dropped
every call with a zero side.

| `.high_conf` rule | calls | recall | precision | F1 |
|---|---|---|---|---|
| drop any zero side (wrong) | 177 | 177/200 | 1.000 | 0.939 |
| drop only 1-vs-0 (RelocaTE2) | 193 | 191/200 | 0.990 | **0.972** |

14 recall for 0.01 precision. Note the *boundary* filter at `:82` genuinely does
use "either side is zero" — a different, broader test that had been conflated
into one helper; now `_has_empty_side` vs `_is_single_read_one_sided`, each
citing its line.

This also **corrects the diagnosis** in `5d7c2fb`: the "19 sites RelocaTE2
resolves two-sided that we don't" was an artefact of our own filter. The real
gap is 5 sites, and at all five RelocaTE3 makes no call at all — so
"second-junction recovery" was the wrong target.

### `50fcbc8` — estimate TSD length from junction reads, not supporting mates

RelocaTE3 reported `TSD=UNK` on 18 of 193 high-confidence Chr3 calls where
RelocaTE2 resolves a concrete TSD.

RelocaTE2 builds its depth pileup from the cluster's **junction** reads, dividing
by their count (`insertionFinder.py:1069-1076` feeding `tsd_finder` at `:843`) —
junction reads all abut the same breakpoint, so the TSD is the run of positions
nearly all of them cover. `_make_insertion` passed `cluster.support` instead: the
supporting *mates*, spread across the library insert. With 13–37 of them no base
is covered by even 60%, so the estimator returned 0 and `_resolve_tsd` fell to
`UNK`.

Diagnosed all 16 disputed sites before writing code: **15 failed at
`tsd_len == 0`**, one at read capture. (The genome fallback in `_resolve_tsd` is
also unreachable from the CLI — `find-insertions` has no genome argument — but
that would have fixed 1 site of 16, so it's left alone.)

| `high_conf` | calls | recall | precision | F1 | UNK |
|---|---|---|---|---|---|
| before | 193 | 191/200 | 0.990 | 0.972 | 18 |
| **after** | 193 | 193/200 | **1.000** | **0.982** | **3** |

Recall and precision improve alongside the labels because `tsd_len` also sets the
call's start/end. F1 gap to RelocaTE2: 0.018 → **0.008**.

### `46e3c35` — record the fix in the baseline note

---

## Verification

- **164 passed / 2 skipped** (default env; both skips blat-gated).
  **0 skipped** in the `blat` env, where the shipped defaults are exercised.
- Acceptance gate unchanged and passing; deliberately pinned to minimap2 so it
  still runs where blat is absent, with `test_run_all_cli_end_to_end` covering
  the shipped defaults.
- Every accuracy number above comes from a real `run-all` on RelocaTE2's test
  data, scored the same way for both callers.

New test files: `relocate2_defaults_test.py`, `output_tiers_test.py`,
`supporting_reads_test.py`, `tsd_depth_input_test.py`.

## Not in this PR

- The synthetic RepeatMasker fixture (understates precision ~0.07).
- The 5 Chr3 sites RelocaTE3 misses entirely — the whole remaining detection gap.
- 3 residual `UNK` calls.
- RelocaTE2 conveniences with no equivalent yet: `-b/--bam` reuse, `--dry_run`,
  `--split`.

Full baseline write-up: `notes/2026-08-12-relocate2-chr3-baseline.md`.
Reproducer for the RelocaTE2 run lives outside both repos at
`rt2_chr3_baseline/run_rt2.sh`.
