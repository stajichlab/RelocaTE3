# Performance Roadmap — R2 Parity & Simulated-Data Benchmarking

> **For Claude / a future maintainer:** This is a **roadmap** document. It sets the priority order for work that keeps R3 correct and measurable. Do NOT attempt to execute this file task-by-task. When a bullet here is ready to work on, write a dated implementation plan alongside it (e.g. `plans/2026-07-XX-fragmentation-fix.md`) with failing-test → fix → measurement structure, using `superpowers:writing-plans`.

## Where things stand (as of the merge of `fix-genotyping-mismatch`)

Rice 10-sample validation against the RelocaTE2 harness:

| Metric | Value | Notes |
|---|---|---|
| TSD strict agreement | 98.29 % | Essentially maxed on this benchmark |
| Status agreement | 93.67 % → 97.32 % | Diagonal grew from 4764 → 5416 |
| Recall vs R2 | 96.6 % | Matched 5565 / 5763 R2 calls |
| Precision vs R2 | 85.3 % | 962 R3-only calls — see caveat below |
| Low-recall sites (r2_T ≥ 2 × r3_T) | 1888 | Down from 5535 |

**Caveat on the precision number:** the R3-only investigation at the end of that PR found that ~90 % of the 962 R3-only calls carry canonical mping TSDs (`TTA`/`TAA`) and most sit at sites R2 saw but tagged `insufficient_data` in its nonref TXT (so R2's characterize gate dropped them). R3 is *emitting real calls R2 could not fully characterize*, not junk — but the confusion matrix cannot distinguish these from false positives without a ground-truth reference. **Any further R2-parity work is measured against a benchmark with known false negatives, and precision numbers below 90 % from now on are more likely benchmark artefacts than R3 defects.**

That caveat is why the top priority below is *not* the fragmentation fix, even though it is the smallest scoped win. It is the simulated-data panel.

---

## Priority 1: simulated-data benchmark panel

Move from measuring R3 against R2 (imperfect benchmark) to measuring both against ground truth (known-placed insertions).

### Motivation
- Absolute recall / precision numbers R3 can quote directly (paper, thesis, releases).
- Distinguishes "R3 found something R2 missed" from "R3 emitted a false positive". Every future parity fix becomes objectively measurable.
- Un-couples R3 quality gates from R2's idiosyncrasies (`insufficient_data` filtering, blat-specific behaviour, etc.).
- Enables regression testing under conditions we cannot get from real rice data (variable coverage, TE family diversity, insertion age, host cultivar).

### Infrastructure that already exists
- `/bigdata/stajichlab/nmath020/github/github_tools/data_sim/simulate_data` — modular CLI with TE insertion, SV placement, and Illumina/ONT/PacBio read modules; SLURM support; tests; plans directory. **Ready to use**, not a build-from-scratch cost.
- `RelocaTE3/tests/acceptance_test.py` — the existing Chr3 2 Mb single-dataset acceptance gate. This is the seed of the panel harness; the panel replaces it with a matrix of datasets.

### Panel design (minimum useful shape)
A 12-sample matrix that stresses the recall/precision tradeoff along axes we care about:

| Axis | Levels |
|---|---|
| Coverage | 10 ×, 30 ×, 60 × |
| TE family | mping (3 bp TSD), one longer element (Ping / Tos17 / RIRE1) |
| Cultivar | Nipponbare + 1 other rice genome |
| Insertions per sample | ~200 spread across all chromosomes |
| Replicates | 1 per cell to start; add later for statistical power |

Total: 3 × 2 × 2 = 12 samples. Small enough for a single SLURM array; large enough to see recall floors per family and coverage.

### Design considerations
- **Placement bias.** Random uniform placement is easy but unrealistic — real TEs land preferentially in AT-rich intergenic regions. Consider mimicking known mping insertion-site preferences so the benchmark measures R3 on the distribution the tool will actually see.
- **Insertion age.** Freshly placed TEs have identical sequence to the library; older TEs accumulate mutations. Simulating a mismatch spectrum tests the mismatch-tolerance parameters in trim (default `--mismatch 2`).
- **TSD variability.** Real mping TSDs are `TTA`/`TAA` but noisy. If the sim always inserts a perfect TSD, R3's TSD-agreement rate is inflated. Consider injecting per-base error at the TSD position at rates matching the sequencing platform.
- **Ground-truth format.** The truth set should be a BED with `chrom start end TE_family TSD strand` at minimum. That schema lines up with what `compare_char.py` already consumes for R2, so a small `normalize_truth.py` can slot into the existing harness.
- **Comparison window.** The current `compare_char.py` uses ± 100 bp for matching; for a ground-truth comparison drop that to ± 10 bp (exact placement is known). Makes recall stricter and more meaningful.
- **Both R2 and R3.** Run RelocaTE2 (blat aligner, as configured in the rice harness) on the SAME simulated inputs. Then both tools are ranked on absolute recall/precision, and every "R3-only" call becomes classifiable as true positive (in truth set) or false positive.

### Milestones
1. **M1 — one dataset end-to-end.** Simulate one 30× mping panel on Nipponbare, run R2 + R3, produce a comparison matrix against truth. Validates the wiring.
2. **M2 — panel matrix.** Extend to the 12-sample grid. Produce a summary table with recall/precision per axis level for both tools.
3. **M3 — panel becomes the acceptance gate.** Replace or supplement `tests/acceptance_test.py` with a lightweight subset of M2 that runs in < 5 min so PRs can gate on it.
4. **M4 — panel becomes the R3-vs-R2 headline.** The comparison matrix from M2 goes into the release notes and any paper prep. R3-only calls that hit the truth set are reported as *R3 advantages*, not false positives.

### Open questions
- Which second TE family? Ping shares mping's TSD but is longer; Tos17 is retrotransposon with a different TSD; RIRE1 is retrotransposon with a long TSD. Different pick, different test coverage.
- Do we simulate paired-end R1/R2 correctly (matching real insert-size distribution) or just single reads? Paired matters for R3's `_call_support_only` path.
- How do we handle nested / adjacent insertions in the truth set? R2 collapses; R3 fragments (see Priority 2). Truth set design determines how each is scored.

---

## Priority 2: insertion fragmentation fix

> **COMPLETED (2026-07-20).** Shipped on branch `fix-insertion-fragmentation`. Plan:
> `plans/2026-07-20-insertion-fragmentation-fix.md`. Verified panel-wide on the
> simulated benchmark: **R3 false positives 693 → 2 across the 9 samples; precision
> 0.70–0.90 → ~1.0 (parity with R2); recall improved or unchanged in all 15
> class×coverage cells (no regressions).** The historical diagnosis below was
> wrong about the function — kept for the record with the correction inline.

### Corrected diagnosis (what actually shipped)
The earlier diagnosis named `_pair_breakpoints`. **That is the wrong function** — it
belongs to the module-level `pipeline.py` path, which the benchmark and the
`find-insertions` CLI do not run. The live path is the `InsertionFinder` class, and
`_write_output` emits **one row per `tsd_start` in a cluster**. With the fixed 3 bp
wildcard TSD (`tsd="..."`), the right junction records `tsd_start = start` while the
left junction records `tsd_start = end - len(tsd)`; when the true TSD is 4–5 bp these
diverge by `true_len - 3` = 1–2 bp, splitting one insertion into two single-sided
rows. Measured: 691/693 (99.7 %) of benchmark FPs were exactly 1–2 bp from a matched
call.

### Fix
`InsertionFinder._merge_offset_starts` (`src/RelocaTE3/insertions.py`): within each
cluster, before `_write_output`, collapse `tsd_start` sub-buckets whose TSD spans
overlap (`gap < captured TSD length`), pooling counts/reads. Split single-sided
fragments recombine into one both-junction call. Distinct insertions have
non-overlapping TSDs and are never merged. R2-parity-preserving: RelocaTE2 avoids the
split by using the correct TSD length, so the merge only reconstructs what the
wildcard mode fragmented.

### Not addressed here (see Priority 5)
The 2 residual FPs are the far-apart complementary-junction case (e.g. L:13 / R:9 tens
of bp apart) that needs RelocaTE2's `TSD_from_read_depth` sub-cluster pairing.

---

## Priority 3: precision investigation (R3-only low-flanker calls)

The 551 low-flanker (avg_flankers < 2) R3-only calls need a second look. Some are real weak calls; some may be spurious junctions from soft-clipped noise. With the simulated panel from Priority 1 this becomes trivial:

- Bucket R3-only by junction count → measure precision per bucket against truth.
- If the precision floor is below ~50 % for avg_flankers < 2, add a minimum-junction filter to `find-insertions` (either a hard `--min-junctions 2` or a MAPQ-aware filter).
- If precision is high even for low counts, leave the calls alone.

**Do not add a filter without the panel.** Filtering against R2 as ground truth would drop real calls.

---

## Priority 4: two-CLI-file consolidation

> **COMPLETED (2026-07-11).** Shipped in two steps: consolidation (PR #19) merged the
> two front-ends, then the canonical CLI was relocated into `src/RelocaTE3/cli.py` with
> a thin `__main__.py` launcher (Option A below is what shipped). **Current state: one
> CLI module, `cli.py`; entry point `relocaTE3 = "RelocaTE3.cli:main"`; `python -m
> RelocaTE3` delegates to `cli.main`.** See `AGENTS.md` → "The CLI lives in one file".
> A follow-up remains: the dispatch cleanup (`build_parser()`, pass `Namespace` instead
> of `**vars(parsed)`, drop per-handler `**kwargs`) — see
> `todo/cmd-run-drops-trim-thresholds.md`. The historical analysis below is kept for record.

`src/RelocaTE3/__main__.py` and `src/RelocaTE3/cli.py` have parallel implementations of the same subcommands. The registered entry point is `__main__:main`, so the harness runs through `__main__.py`; `cli.py` is documented but not wired. This has produced **two null-result SLURM runs** during recent parity work — each time I edited the wrong path (`map_reads_to_genome` instead of `map_genome_minimap`; function-based `find_insertions` instead of `InsertionFinder` class). CONTEXT.md warns about it but the warning is easy to miss.

### The right consolidation
Pick one file as the authoritative source. Two options:

**Option A — `cli.py` becomes the only file.** Rewrite `pyproject.toml`'s entry point to `RelocaTE3.cli:main`; delete `__main__.py`. Cleaner but changes the module invocation semantics (`python -m RelocaTE3` breaks unless the entry moves).

**Option B — `__main__.py` becomes the only file, `cli.py` deleted.** Keep the current entry point; move the docstrings and any unique subcommand handling from `cli.py` into `__main__.py`. Lower migration risk.

Either way, each subcommand's *implementation* (the `cmd_foo` body) should call ONE shared function in an implementation module. No parallel `map_reads_to_genome` / `map_genome_minimap` duplication.

### Prerequisites
- Priority 1 completed. Otherwise a consolidation regression looks like a parity regression and you cannot distinguish them.
- All existing tests passing. The consolidation is behaviour-preserving; the test suite is the gate.

---

## Priority 5: `supporting_junction` sub-cluster port

The pre-existing plan `plans/2026-06-25-tsd-supporting-junction-port.md` targeted ~960 cells; that count is now down to ~40 after the trim-recall fix. The remaining sites are cases where R2 recovers a TSD via `TSD_from_read_depth` sub-cluster splitting that R3's class path doesn't do.

40 sites out of 5565 matched = 0.7 %. Consider whether it is worth doing at all. The written plan can be closed out with a "resolved incidentally" note, or executed for completeness. The simulated panel from Priority 1 will tell you which of those 40 are real biological insertions vs harness artefacts.

---

## Priority 6: acceptance-test breadth

Once M3 from Priority 1 lands, retire the single Chr3 2 Mb `tests/acceptance_test.py` in favour of a fast subset of the simulated panel. Recall/precision floors are asserted per family × coverage combination. A PR that regresses any cell fails CI.

### Considerations
- Runtime budget. Real acceptance tests need to finish in under 10 min or PR authors skip them locally. The subset used for CI must be a strict subset of the M2 panel — probably one coverage level, one family, ~50 insertions, one cultivar.
- Full-panel runs happen nightly (or per-release), not per-PR.

---

## Priority 7: low-quality / MAPQ read filtering

Still deferred from `PLAN.md`. Whether it matters is a question the simulated panel answers directly. If R3's precision floor is bottlenecked by soft-clipped low-MAPQ junction reads, this becomes a real fix. If not, deprioritize.

---

## Open questions across the whole roadmap

- Do we want R3 to *match* R2 or *exceed* it? Every parity gate implicitly answers "match"; the simulated panel implicitly answers "exceed if truth agrees". The two goals occasionally conflict.
- What is the release cadence? If R3 ships to lab users before Priority 1 lands, they will report bugs that the simulated panel would have caught. Consider gating a `v0.2` release on Priority 1.
- Ground truth for real biology? The simulated panel is measurable but synthetic. A small hand-curated set of ~20 loci per rice sample (verified by long-read data or PCR) would triangulate the simulated numbers with reality. Do NOT let this become an ongoing project — it is a one-off validation.

---

## Cross-references

- Historical narrative & completed work: `plans/PLAN.md`
- Long-term capability roadmap (Rust, Nextflow, long reads, pangenomes): `plans/FEATURES.md`
- Completed parity plans (executed, in-tree for the audit trail):
  - `plans/2026-06-24-tsd-depth-inference.md`
  - `plans/2026-06-25-tsd-class-path-port.md`
  - `plans/2026-06-25-tsd-supporting-junction-port.md`
  - `plans/2026-07-02-trim-recall-parity.md`
- Retired dead-end (kept as documentation): `plans/2026-06-26-genotype-status-parity.md`
