# RelocaTE2 vs RelocaTE3: stage-by-stage logic and parameter audit

2026-08-18

Purpose: RelocaTE3 is meant to be a modernised RelocaTE2 — same aligners, same
logic, same parameters. This audit walks the pipeline stage by stage and records
every difference found, separating **measured** findings from **code-read**
ones.

---

## Headline: RelocaTE3 aligns TE-internal ("middle") reads to the genome; RelocaTE2 does not

This is the single largest divergence and it explains the riceTElib collapse.

A read that lies **entirely inside** a TE gets tagged `:middle` by both callers'
trim step (`relocaTE_trim.py:447`, `librelocate._trim_record`). It carries no
flank, so it cannot mark a breakpoint — it is pure transposon sequence.

**RelocaTE2 removes these before genome alignment.** `clean_pairs_memory.py`
states it directly:

> `*.unPaired.fq`: contains trimmed reads that do not have mates found and
> contain the mate pairs of reads that matched to middle of repeat **only if the
> mate pair is not repeat, but not reads themselve as they are part of repeat**

So RelocaTE2 keeps the *genomic mate* of a middle read (a real supporting read)
and discards the middle read itself.

**RelocaTE3 aligns them all.** `genome_align.build_flank_pairs` reads every
record in `*.flankingReads.fq` — middle-tagged included — and emits it for
alignment.

Measured on riceTElib `cov30x_rep1`:

| | RelocaTE2 | RelocaTE3 |
|---|---|---|
| middle reads in files fed to bwa | **0** (`.matched` ×2, `.unPaired.fq`) | — |
| middle-tagged records in the genome BAM | — | **2,546,333 (64.6 %)** |
| flanking FASTQ composition | — | left 186,141 junction / 1,308,723 middle; right 171,612 / 1,242,166 |

### Why this wrecks precision

`insertions._stream_clusters` does two things with every mapped record:

1. **`current.extend(gstart, gend)`** — the record extends the cluster, and
   clusters chain while successive reads fall within `RANGE_ALLOWANCE` (1000 bp).
2. Records that do not match `:(start|end):[53]$` go to **`current.support`**.

Middle reads are TE sequence, so they map to *every reference copy of that
family*. At each copy they pile up, glue clusters together across up to 1000 bp,
and are counted as supporting reads.

That accounts for observations that had no explanation before:

- RelocaTE3's false positives are **66 % within 100 bp of a reference TE copy**,
  versus 38 % of its own true positives and 29 % of RelocaTE2's FPs.
- Of one-sided calls, RelocaTE3 has bracketing support at **3562 of 4873** sites;
  RelocaTE2 at **1 of 64**.
- riceTElib (500 families, 256,637 reference copies) is catastrophic while mPing
  (1 family) is nearly fine.
- RelocaTE3 emits 2430 calls to RelocaTE2's 1076 *from fewer TE-containing reads*
  (2.91 M vs 3.33 M).

**Fix:** exclude `:middle` records from genome alignment, and recover the genomic
mate of each middle read as a supporting read, exactly as
`clean_pairs_memory.py` does.

---

## Stage-by-stage

### Step 0 — reference TE annotation
`existingTE_RM_ALL` → `reference_te.parse_repeatmasker`. No difference found.
Both build a `{chrom: {start:…, end:…}}` boundary table used only to drop
one-sided calls at intact-copy edges.

### Step 2 — FASTQ → FASTA
RelocaTE2 uses `seqtk`. RelocaTE3 did this in Python until 2026-08-16; now also
`seqtk`. **Resolved.**

### Step 3 — TE search + trim
- **Aligner:** both default to `blat`. **Versions differ (v35 vs v36) but were
  tested head-to-head and produce byte-identical PSL** — not a factor.
- **BLAT parameters:** RelocaTE2 runs `-minScore=10 -tileSize=7`
  (`relocaTE2.py:545`). RelocaTE3 ran BLAT defaults until 2026-08-16; now
  matches. **Resolved**, though see the caveat below.
- **Trim branch logic:** identical. `relocaTE_trim.py:381,412` and
  `librelocate._trim_record` apply the same `tStart<=2` / `tEnd>=tLen-3` /
  ends-align / `match+mismatch >= len_cut_match` /
  `length-(match+mismatch) >= len_cut_trim` / `mismatch <= allowance` tests.
- **Difference:** RelocaTE3 collapses `--mismatch_junction` into `--mismatch`.
  Same default (2), so no behavioural change unless a user sets them apart.

### Step 4 — genome alignment (**the big one**)
RelocaTE2 (`relocaTE_align.py` + `clean_pairs_memory.py`) builds **three**
populations:

1. `*.matched` — flank paired with its mate (mate taken from the other trimmed
   flank file *or* from the original FASTQ) → `bwa aln` ×2 + **`sampe`**
2. `*.unPaired.fq` — flanks with no mate, **plus genomic mates of middle reads**
   → `bwa aln` + **`samse`**
3. `*.fullreads` — untrimmed junction reads + pairs → `sampe`, used as the
   false-junction control

RelocaTE3 (`genome_align.build_flank_pairs` / `align_flanks_anchored`) builds:

1. flank + genomic mate pairs (only when the mate did **not** match a TE)
2. flanks without a mate → single-end
3. fullreads

Differences:
- **Middle reads are aligned by RelocaTE3, excluded by RelocaTE2** (headline).
- **Middle reads' genomic mates**: RelocaTE2 explicitly recovers them as
  supporting reads; RelocaTE3 does not (its `te_ends` test excludes any pair
  whose other end matched a TE, which is exactly the middle-read case).
- **bwa version** 0.6.2 vs 0.7.19 — tested head-to-head on 200,000 flanks:
  identical placements, identical XT/X0/X1, 0 disagreements. **Not a factor.**

### Step 5 — insertion calling
- **Cluster range:** both use 1000 bp (`range_allowance` /
  `_RANGE_ALLOWANCE`). Same.
- **Breakpoint pairing:** RelocaTE2 pairs each left position with its *nearest*
  right within 100 bp (`TSD_from_read_depth:610-645`). RelocaTE3 chains
  positions into groups by `SUBCLUSTER_GAP = 25` and takes the most-supported
  position per side. **Different algorithm** — RelocaTE3 changed this
  deliberately; a chain of ≤25 bp steps can span far more than 25 bp.
- **TSD length:** RelocaTE2 computes depth-based (`tsd_finder`) **and** geometric
  (`TSD_len_calculate`), takes the geometric on disagreement, and emits nothing
  when it is not positive (`:800-818`). RelocaTE3 had only the depth estimate and
  emitted a collapsed call when geometry failed; the reconciliation + emission
  gate were ported 2026-08-17. **Resolved, but measured F1-neutral on riceTElib.**
- **Read admission:** RelocaTE2 has **no MAPQ gate** — MAPQ<29 only *marks* a
  read low quality (`:1523,1539`) and calls resting solely on such reads are
  dropped later (`:226-241`). RelocaTE3 defaulted `min_mapq=1`. **Resolved
  2026-08-16** (default now 0).
- **One-sided classes:** RelocaTE2 labels one-sided calls
  `supporting_junction` / `singleton` / `insufficient_data` and keeps only
  `supporting_junction` downstream. RelocaTE3 emits **zero**
  `supporting_junction` rows — the class is unreachable because
  `require_both_junctions` drops one-sided calls before the sentinel is
  assigned, and its supporting-read counts are not R2-comparable anyway
  (see `todo/supporting-read-counts-not-r2-equivalent.md`). **Open.**
- **`ST:` field** hardcoded to 0 in RelocaTE3. **Open.**

### Step 6 — reference/shared insertions
`relocaTE_absenceFinder.py` implements true *absence* detection (TE present in
the reference, deleted in the sample). RelocaTE3 does not. **Open, known,
deliberate** (`plans/PLAN.md` Phase 4 "deferred").

### Step 7 — genotyping and filtering
- **characterizer gate:** RelocaTE2 requires
  `($left >= 1 and $right >= 1) or $TSD eq 'supporting_junction'`
  (`characterizer.pl:91`). RelocaTE3 ported this with `or` instead of `and`.
  **Resolved 2026-08-14.**
- **Tiering:** `clean_false_positive.py` greps out
  `singleton|insufficient_data|supporting_reads` for the headline GFF and the
  literal `1-vs-0` junction pattern for `high_conf`. RelocaTE3 matches.
- **Excision/VCF footprint analysis** from `characterizer.pl` — not ported.
  **Open, known.**

---

## Parameter comparison

| RelocaTE2 | default | RelocaTE3 | default | status |
|---|---|---|---|---|
| `--size` | 500 | `-s` / `insert_size` | 500 | same |
| `--len_cut_match` | 10 | `--min-match` | 10 | same |
| `--len_cut_trim` | 10 | `--min-trimmed` | 10 | same |
| `--mismatch` | 2 | `--mismatch` | 2 | same |
| `--mismatch_junction` | 2 | *(collapsed into `--mismatch`)* | — | **merged** |
| `--aligner` | blat | `--te-aligner` | blat | same |
| *(bwa aln, hardcoded)* | — | `--genome-aligner` | bwaaln | same |
| `--cpu` | 1 | `--threads` | 1 | same |
| `-d/--distance` (clean_false_positive) | 3 | `--distance` | 3 | same |
| *(no MAPQ admission gate)* | — | `--min-mapq` | 0 | same since 2026-08-16 |
| BLAT `-minScore` / `-tileSize` | 10 / 7 | same | 10 / 7 | same since 2026-08-16 |
| *(none)* | — | `--require-both-junctions` | True | **R3-only flag**, but equivalent to R2's characterizer gate |
| `--mate_1_id` / `--mate_2_id` | `_1` / `_2` | inferred from filenames | — | different mechanism |
| *(cluster range 1000)* | 1000 | `_RANGE_ALLOWANCE` | 1000 | same |
| *(pair within 100 bp)* | 100 | `SUBCLUSTER_GAP` | 25 | **different** |

---

## Ranked list of remaining work

1. ~~**Exclude `:middle` reads from genome alignment; recover their genomic mates
   as supporting reads.**~~ **DONE 2026-08-18.** `build_flank_pairs` now skips any
   record without a `:(start|end):[53]` tag. Mate recovery needed no change:
   `recover_support_mates` already reads the same `read_repeat` table and pulls
   the genomic mate of every TE-containing read (middle reads included) whose
   mate did not itself match a TE — it simply never saw those mates before,
   because the middle read had already claimed them as pair partners. Tests:
   `tests/genome_align_test.py::test_build_flank_pairs_skips_te_internal_middle_reads`
   and `::test_middle_read_mate_is_still_recovered_as_support`. Not yet measured
   on a benchmark.
2. ~~**Breakpoint pairing**~~ **DONE 2026-08-19.** `_pair_breakpoints` is now a
   direct port of `TSD_from_read_depth` (:603-770): walk left positions
   ascending, pair each with its *nearest* right within `PAIR_MAX_DISTANCE`
   (100 bp), emit unclaimed rights one-sided, plus RelocaTE2's several-on-one-side
   and one-and-one branches. `SUBCLUSTER_GAP` (25) and `SPLIT_MIN_SEPARATION`
   (10) are gone — both were RelocaTE3 inventions. The adjacent-insertion
   protection they provided now comes from where RelocaTE2 gets it: the pair is
   made, then suppressed at emission when `TSD_len_calculate` is non-positive
   (:818).
   **Measured 2026-08-19: worse on its own.** mPing precision 1.000 → 0.992/0.978/0.976
   (0 → 5/22/30 FPs), F1 0.572/0.784/0.893 → 0.571/0.778/0.883; riceTElib 5x
   F1 0.409 → 0.360 (calls 671 → 970, FP 227 → 525). **Recall unchanged on both
   panels** — the port purely added candidates. Cause: see item 2b.

2b. **Cluster-level candidate arbitration in `write_output` (:257-330) — NOT
   PORTED.** RelocaTE2 generates many candidate starts per cluster (which is why
   the pairing port produces more) and then *arbitrates among them*:
   - all candidates two-sided → keep all
   - some two-sided (`start_both_junction > 0`) → keep two-sided candidates, plus
     one-sided ones with **>= 3 junction reads**; discard the rest, reasoning
     "if we found both junction for one insertion we should find both junction
     for others too"
   - none two-sided → keep only candidates not overlapping an existing TE, and
     only those with >= 3 junction reads

   RelocaTE3 emits every sub-insertion independently with no arbitration. The
   pairing port supplied RelocaTE2's *generator* without its *filter*, which is
   exactly why false positives rose while recall did not move. Porting this is
   the direct completion of item 2.

   **DONE 2026-08-19 (`_arbitrate_cluster`), and measured 2026-08-20: it is a
   NO-OP in the shipped configuration.** Output was byte-identical on both
   panels. Every branch arbitration affects concerns *one-sided* candidates, and
   `--require-both-junctions` (on by default) already removes those before they
   are written. Arbitration only changes behaviour with
   `--no-require-both-junctions`. The code is correct and RelocaTE2-faithful; it
   simply has nothing to do here.

   Consequence: **the remaining riceTElib gap is entirely in two-sided call
   generation**, which neither arbitration nor the junction policy touches.
   Measured at riceTElib cov5x_rep1: raw call tables are close (RelocaTE2 308,
   RelocaTE3 327) but RelocaTE2's characterize gate drops ~40 % of its rows
   (one-sided sentinels) while RelocaTE3's drops ~1 % (all rows are already
   two-sided) -- so RelocaTE3 produces roughly 1.8x more *two-sided* calls.
   RelocaTE3 also places 2.4x more reads on the genome (232,748 vs 97,496,
   excluding the fullreads control).
3. **Supporting-read counts** — make `SR:`/`SL:`/`ST:` mean what RelocaTE2 means,
   then re-land the `supporting_junction` exemption.
4. Absence detection (step 6) and the excision/VCF analysis (step 7).
5. `--mismatch_junction` as a separate knob.

## Benchmark fairness: how each caller is actually invoked

Audited 2026-08-19 against `relocate-benchmark/callers/*/run.sh` and
`config/benchmark.toml`. Two asymmetries found — the comparison is not
parameter-controlled.

### 1. Library insert size differs, and neither value matches the data

| | value | source |
|---|---|---|
| RelocaTE2 | **250** | `config/benchmark.toml` `[callers.relocate2] size = 250` → `--size` |
| RelocaTE3 | **500** | adapter passes nothing → CLI default |
| the simulated data | **300** | `run_metadata.json` `"fragment_size": 300` (both panels, `fragment_std` 30) |

RelocaTE2's own default is 500, so the 250 is a benchmark choice with no
recorded rationale (it arrived in the initial config commit, `1db8057`).

This parameter sets the span of a mate-only supporting cluster —
`insert_size * 1.2` in both callers (RelocaTE2 :446,
RelocaTE3 `call_support_only`). RelocaTE3 is therefore using a 600 bp window
where RelocaTE2 uses 300 bp. Those calls land in `all_nonref_supporting`
rather than the scored set, so the direct effect on recall/precision is
probably small — but it is an uncontrolled 2x difference and should be fixed.

**Fix:** pass the same value to both. 300 matches the simulation; if the goal is
strictly "RelocaTE2's defaults", use 500 for both. Either is defensible; the
current split is not.

### 2. The genotyping BAM is built with a different aligner for each caller

| | aligner |
|---|---|
| RelocaTE2 | `bwa mem` (bwa 0.7.19 container), passed via `--bam` |
| RelocaTE3 | `minimap2 -a -x sr`, passed to `characterize -b` |

This BAM is what spanner counting reads, so it drives the
homozygous/heterozygous/somatic call. It does not affect detection
recall/precision, but every genotype-status comparison between the callers is
confounded by it. **Fix:** build it once per sample with one aligner and give
the same file to both adapters.

### Verified equal

`--mismatch` 2 both; `len_cut_match`/`len_cut_trim` 10/10 both (RelocaTE2 uses
its defaults, RelocaTE3 is passed 10/10 explicitly); TE aligner `blat` both;
genome placement `bwa aln` both; `--mismatch_junction` 2 both (RelocaTE2 default,
RelocaTE3 folded into `--mismatch`); threads/cpu the same; RelocaTE2 runs
`--step 1234567` (all steps); both receive the same `TE_LIBRARY`, `REFERENCE` and
`REPEATMASKER` paths from the dataset env; and both chunk BLAT input at 200,000
reads (RelocaTE2 `split_fq -s 200000`, RelocaTE3 `BLAT_CHUNK_SEQS`).

## Caveat on the BLAT sensitivity setting

RelocaTE2's `-minScore=10 -tileSize=7` is correct parity and is required for
mPing recall (it closed ~80 % of that gap). It was measured as catastrophic on
riceTElib — but that measurement was taken **with middle reads still being
aligned**. Since middle-read volume scales with TE-library size, item 1 above is
the most likely explanation for why the sensitised setting looked unsafe.
**Re-measure the BLAT setting after fixing item 1 before reverting it.**
