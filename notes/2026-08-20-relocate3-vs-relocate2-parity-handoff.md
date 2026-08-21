# Handoff: RelocaTE3 does not match RelocaTE2 on the multi-TE benchmark

**Date:** 2026-08-20
**Status:** open, root cause not found
**Audience:** an agent picking this up cold

---

## 1. The goal (first principles)

RelocaTE3 is meant to be a **modernised RelocaTE2**, not an improved one. Same
aligners, same parameters, same algorithms — modern implementation. Parity on
the benchmarks is the acceptance bar. Only after parity is demonstrated should
anyone consider changing behaviour to score better.

This matters because several past "improvements" in RelocaTE3 turned out to be
undocumented divergences from RelocaTE2 that were individually defensible and
collectively made the two callers incomparable. **A benchmark gain produced by
logic RelocaTE2 does not have should be treated as a red flag, not a win.**

---

## 2. The issue in one paragraph

On a **single-TE-family** dataset (mPing) RelocaTE3 is close to RelocaTE2. On a
**multi-family** dataset (riceTElib, 500 TE families) RelocaTE3 emits roughly
twice as many calls as RelocaTE2 at similar recall, so its precision is far
worse — 0.46 vs 0.81 at 5x coverage. Recall is essentially identical between the
two callers; **the entire gap is excess false-positive calls.** The excess is in
**two-sided** calls (junction reads on both sides of the breakpoint), which is
the population that no junction-policy or one-sided filter touches.

---

## 3. Measurements

### Metric definitions

Scored by `scoring/score_calls.py` in the benchmark repo: greedy nearest-match,
strictly 1:1, truth-driven. For each truth event take the closest unused call
within **±10 bp** whose `te_family` matches after normalisation (strip
`#Class/Family`, lowercase, collapse the somatic vocabulary). Unmatched calls are
false positives.

- **recall** = truth events detected / truth events
- **precision** = matched calls / total calls (from `precision.tsv`; this is the
  trustworthy one)
- **F1** = harmonic mean of the two

Note the 1:1 rule: a duplicated breakpoint for a real insertion scores as a false
positive, not as a second hit.

### Current numbers (2026-08-20 runs)

**mPing** (single family) — RelocaTE3 close but slightly behind:

| caller | cov | calls | FP | recall | precision | F1 |
|---|---|---|---|---|---|---|
| relocate2 | 5x | 610 | 0 | 0.407 | 1.000 | 0.578 |
| relocate3-blat-bwaaln | 5x | 606 | 5 | 0.401 | 0.992 | 0.571 |
| relocate2 | 15x | 989 | 1 | 0.659 | 0.999 | 0.794 |
| relocate3-blat-bwaaln | 15x | 990 | 22 | 0.645 | 0.978 | 0.778 |
| relocate2 | 30x | 1221 | 0 | 0.814 | 1.000 | 0.897 |
| relocate3-blat-bwaaln | 30x | 1239 | 30 | 0.806 | 0.976 | 0.883 |

**riceTElib** (500 families) — the real problem. **Only 5x is current**; the 15x
and 30x rows in the report are stale (their run directories were deleted and only
coverage 5 was re-run — verify before quoting):

| caller | cov | calls | FP | recall | precision | F1 |
|---|---|---|---|---|---|---|
| relocate2 | 5x | 552 | 103 | 0.299 | 0.813 | 0.438 |
| relocate3-blat-bwaaln | 5x | **970** | **525** | 0.297 | **0.459** | 0.360 |

Recall 0.297 vs 0.299 — statistically the same. Calls 970 vs 552.

### Where the excess sits (riceTElib cov5x_rep1, one sample)

| | RelocaTE2 | RelocaTE3 |
|---|---|---|
| reads placed on the genome (excl. fullreads control) | 97,496 | **232,748** |
| raw `all_nonref_insert.txt` rows | 308 | 327 |
| rows surviving to the characterized table | ~184 | ~323 |

The raw tables are nearly the same size, but RelocaTE2's characterize gate drops
~40 % of its rows (they are one-sided, carrying sentinel strings in the TSD
column) while RelocaTE3's drops ~1 % because its rows are already all two-sided.
So **RelocaTE3 produces ~1.8x more two-sided calls from ~2.4x more placed
reads.** Both of those ratios are unexplained.

---

## 4. Benchmark details

### Repository

`/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/relocate_benchmark/relocate-benchmark`
(HEAD `333b2c4`, branch `feat/configurable-junction-policy`)

Runs callers and scores them; it does **not** generate data. Simulated reads,
genomes and truth come from a separate project and are referenced read-only by
`config/benchmark.toml` → `[datasets.<key>].panel_root`.

### Datasets

| key | truth events | samples | TE library | notes |
|---|---|---|---|---|
| `mping` | 500 | 9 | **1 seq, 430 bp** | single element, easy case |
| `ricetelib` | 500 across 10 TE groups | 9 | **500 seqs, 1.27 Mb** | the discriminating one |
| `ricetelib_divergence` | 9000 | 54 | same | TE divergence 0–20 %; both callers collapse past 10 % |

All panels: Chr1 of MSU_r7, coverage 5x/15x/30x x 3 replicates, match window
±10 bp. Truth composition per panel is 100 each of homozygous, heterozygous, and
somatic at cellular fraction 0.1 / 0.2 / 0.4 (expected VAF 0.05 / 0.10 / 0.20).
Simulated fragment size is **300** (`fragment_std` 30) per `run_metadata.json`.

The reference carries **256,637** annotated TE copies (`existingTE.bed`), which
matters for anything involving reference-TE proximity.

### Key paths

```
config/benchmark.toml                      dataset paths, caller registry, scoring params
callers/relocate2/{env.sh,run.sh,normalize.py}
callers/relocate3/{env.sh,run.sh,normalize.py}
truth/<dataset>/{truth.tsv,samples.tsv,per_sample/<sample>.tsv}
runs/<dataset>/<caller>/<sample>/          raw caller output (GITIGNORED)
reports/datasets/<dataset>/precision.tsv   per-sample precision  <- most trustworthy
reports/datasets/<dataset>/correctness.tsv per-class recall/status/TSD
reports/datasets/<dataset>/per_sample/<caller>/<sample>/{matches,false_positive_calls}.tsv
scoring/score_calls.py                     the scorer
```

RelocaTE2 intermediates worth knowing (under `runs/<ds>/relocate2/<sample>/raw/repeat/`):
`blat_output/`, `te_containing_fq/`, `flanking_seq/`, `bwa_aln/`, `results/`.
`results/<chrom>.repeat.reads.list` gives **the exact read names RelocaTE2
assigned to each called site** — the most useful file for a per-locus diff.

RelocaTE3 intermediates (under `runs/<ds>/<caller>/<sample>/raw/`):
`te_containing/<sample>.read_repeat_name.txt`, `flanking/`,
`<sample>.repeat.bwaaln.sorted.bam`, `results/`.

### Running it

```bash
cd <benchmark repo>
bash pipeline/submit_benchmark.sh --dataset mping
bash pipeline/submit_benchmark.sh --dataset ricetelib --caller relocate3-blat-bwaaln --coverage 5
```

**Always run through SLURM.** Never on the login node — a previous session did
this and made the cluster unusable for other users. The BLAT stage is genuinely
heavy: riceTElib is ~2 h at 5x, ~6.5 h at 15x, ~13 h at 30x per sample, peak RSS
up to 37 GB. mPing is 8 / 21 / 49 min.

### Two traps that have wasted time

1. **`run.sh` exits early when `runs/<ds>/<caller>/<sample>/.run_complete`
   exists**, so a "re-run" silently re-scores stale output. Delete the run
   directory, then verify freshness afterwards (check `.run_complete` timestamps
   and the trim read counts).
2. **Report tables combine per-caller directories of different vintages.** A row
   in `correctness.tsv` is not necessarily from the same code as its neighbour.
   Check `reports/datasets/<ds>/per_sample/<caller>/` mtimes before quoting a
   head-to-head.

### Environments

- RelocaTE2 runs from a **digest-pinned BioContainer** via apptainer shims
  (`callers/relocate2/images.txt`); its bioconda package is a dead py2.7 build.
  `relocaTE2.py` executes *inside* the container, so it uses the container's own
  blat 35 and bwa 0.6.2 — the host shims are not visible to it.
- RelocaTE3 runs from the **dev repo's editable pixi env** (`RT3_REPO` in
  `config/benchmark.toml`). This means benchmark runs execute the working tree,
  **not** the pinned rev in `callers/relocate3/pixi.toml`. Check the dev repo's
  branch and file mtimes before submitting.

---

## 5. State of the code

RelocaTE3 repo:
`/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/RelocaTE3`
HEAD `16ab499`, branch `main`. **Everything below is uncommitted** — working-tree
changes only, in both repos.

Modified: `src/RelocaTE3/{aligners,characterize,cli,genome_align,insertions}.py`,
`pixi.toml`/`pixi.lock`, and eight test files. New docs:
`docs/relocate2-vs-relocate3-audit.md` (the stage-by-stage audit — read this
first), `docs/blat-te-search-performance.md`.

Benchmark repo: `callers/relocate3/run.sh` and `config/benchmark.toml` modified;
`scripts/compare_bwa_aln_versions.sh` added.

240 tests pass, including the Chr3 2 Mb acceptance gate.

---

## 6. What has been tried

Each of these was a **real, verified divergence from RelocaTE2** and was fixed by
porting RelocaTE2's version. None of them closed the riceTElib gap.

| # | Divergence | RelocaTE2 reference | Result |
|---|---|---|---|
| 1 | characterize gate used `or` where RelocaTE2 has `and` | `characterizer.pl:91` | **Large fix.** riceTElib precision 0.066 → 0.808 in a re-score. Kept. |
| 2 | `min_mapq` defaulted to 1; RelocaTE2 has no MAPQ admission gate (it only *marks* MAPQ<29 low quality) | `relocaTE_insertionFinder.py:1521-1558`, `:226-241` | +4 true calls on mPing 30x, 0 new FPs. Kept. |
| 3 | BLAT ran at defaults; RelocaTE2 uses `-minScore=10 -tileSize=7` | `relocaTE2.py:545` | Closed ~80 % of the mPing recall gap (0.536 → 0.572 F1 at 5x). Required adding chunked/parallel BLAT + seqtk to stay tractable. Kept — **but see the warning below.** |

> **The TE-search sensitivity is the change most strongly correlated with the
> riceTElib gap, and this is not explained.** In the 2026-08-17 riceTElib run,
> with all other code identical, the variants whose TE search used RelocaTE2's
> sensitised BLAT scored F1 ≈ 0.45 at 30x (`relocate3-blat-bwaaln` 0.456,
> `relocate3-blat-bwa` 0.453) while `relocate3-bowtie2-bwa`, which uses a
> different and less sensitive TE search, scored **0.719** — above RelocaTE2's
> 0.696. Before the BLAT change, blat+bwa-aln at BLAT's own defaults scored
> 0.441 / 0.632 / 0.715.
>
> The obvious reading is "revert the BLAT parameters", but that would be
> reverting to something RelocaTE2 does not do. **RelocaTE2 runs these exact
> parameters and achieves precision 0.834.** So the parameters are not the
> defect; something downstream in RelocaTE3 amplifies the extra BLAT hits into
> calls in a way RelocaTE2's pipeline does not. Finding that amplifier is
> probably the same problem as open questions 1 and 2 below.
>
> A previous hypothesis — that `:middle` reads were the amplifier — was tested
> and **disproven** (fix #5: it removed 64.6 % of the genome BAM and changed the
> 30x call set by one call).
| 4 | TSD length: only the depth estimator was ported, not the geometric one nor its emission gate | `:800-818` (`TSD_len_calculate`, `if TSD_len > 0`) | F1-neutral on riceTElib (0.459 → 0.456). Kept as correct. |
| 5 | `:middle` (TE-internal) reads were aligned to the genome; RelocaTE2 excludes them and keeps only their genomic mates | `clean_pairs_memory.py` docstring | Removed 2,546,333 of 3,942,639 BAM records (64.6 %). **Changed the 30x call set by one call.** Slight recall loss at 5x/15x. Kept (R2-faithful, 2.8x smaller BAM). |
| 6 | breakpoint pairing used a RelocaTE3-invented scheme (chain by 25 bp, take most-supported per side) instead of nearest-right-within-100bp | `TSD_from_read_depth:603-770` | **Made things worse alone**: mPing precision 1.000 → 0.976, riceTElib 5x F1 0.409 → 0.360, recall unchanged. |
| 7 | no cluster-level arbitration among a cluster's candidates | `write_output:257-330` | **No-op.** Byte-identical output. Every branch concerns one-sided candidates, which `--require-both-junctions` already removes. |

### Ruled out by direct experiment

- **BLAT version** (RelocaTE2 v35 vs RelocaTE3 v36): run head-to-head with
  identical parameters on the same 1000-read query — **byte-identical PSL**
  (same md5).
- **bwa version** (RelocaTE2 0.6.2 vs RelocaTE3 0.7.19, both `aln`+`samse`): run
  head-to-head on 200,000 flanking reads (SLURM job 27546533,
  `scripts/compare_bwa_aln_versions.sh`) — 199,847/199,847 mapped to the **same
  position**, 0 XT-tag disagreements, identical MAPQ/X0/X1 distributions.
- **Trim over-calling**: RelocaTE2 identifies *more* TE-containing reads than
  RelocaTE3 (3.33M vs 2.91M at riceTElib 30x), so trim is not producing excess.
- **Trim branch logic**: verified line-by-line identical to
  `relocaTE_trim.py:381,412`.

### Verified-equal parameters

`--mismatch` 2, `len_cut_match`/`len_cut_trim` 10/10, TE aligner blat, genome
placement bwa aln, `--mismatch_junction` 2, threads, RelocaTE2 runs all seven
steps, identical TE library / reference / RepeatMasker inputs, BLAT input chunked
at 200,000 reads on both sides.

### Known-unequal parameters (NOT yet fixed — fix before drawing conclusions)

1. **Library insert size.** RelocaTE2 is given `--size 250`
   (`config/benchmark.toml`); RelocaTE3 is passed nothing and uses its default
   **500**; the simulated data is **300**. RelocaTE2's own default is 500. This
   sets the mate-only supporting-cluster span (`insert_size * 1.2` in both). Those
   calls go to `all_nonref_supporting` rather than the scored set, so the scoring
   impact is probably small — but it is an uncontrolled 2x difference.
   RelocaTE3 exposes the knob (`find-insertions -s/--size`); the adapter simply
   never passes it.
2. **The genotyping BAM uses a different aligner per caller** — `bwa mem` for
   RelocaTE2, `minimap2 -a -x sr` for RelocaTE3. Feeds spanner counting, so it
   confounds hom/het/somatic comparisons. Does not affect detection.

---

## 7. Open questions

These are the facts that remain unexplained. They are stated deliberately without
a preferred hypothesis.

1. Why does RelocaTE3 place **2.4x more reads** on the genome than RelocaTE2
   (232,748 vs 97,496 at riceTElib cov5x_rep1, excluding the fullreads control),
   given that RelocaTE2 identifies *more* TE-containing reads at trim?
2. Why does RelocaTE3 form **~1.8x more two-sided calls** than RelocaTE2 from
   those reads?
3. Why is the discrepancy specific to the multi-family library? mPing (1 family)
   is nearly at parity; riceTElib (500 families, 256,637 reference copies) is not.
   Candidate axes: library size, reference-copy density, TE family diversity,
   paralogy between families.
4. Earlier (before fixes 5-7) RelocaTE3's false positives were measured at
   **66 % within 100 bp of a reference TE copy**, versus 38 % of its own true
   positives and 29 % of RelocaTE2's false positives. That measurement has not
   been repeated since. It may or may not still hold.

### A suggested next step (not the only one)

RelocaTE2 writes `results/<chrom>.repeat.reads.list`, mapping each called site to
the exact read names it used. Picking ~10 loci where RelocaTE2 emits one
two-sided call and RelocaTE3 emits several, and diffing both callers' state at
those exact positions — which reads each placed, at which breakpoints, with what
counts — would distinguish "extra placed reads" from "extra breakpoint positions"
from "extra pairs per breakpoint set". Those imply different fixes. This has not
been done.

**Caution for whoever picks this up:** five plausible mechanisms have now been
identified from reading code, each confirmed as a genuine divergence, and none
resolved the gap. Prefer direct measurement over another code-reading hypothesis,
and re-measure claims in this document rather than trusting them — several
earlier conclusions in this project's notes turned out to be wrong or to have
been measured under a confound.
