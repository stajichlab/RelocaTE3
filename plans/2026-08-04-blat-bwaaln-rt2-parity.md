# BLAT + bwa-aln Parity with RelocaTE2

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.
> **For a future session:** Read the "Diagnosis" section in full before touching code. The `bwaaln` genome backend already matches RelocaTE2 (PR #39, commit `60082c7`). The remaining, previously-undocumented gap is on the **BLAT TE-search side**: RelocaTE3 runs BLAT at its insensitive defaults while RelocaTE2 runs it with sensitized `-minScore=10 -tileSize=7`. Chunking (`--split`) is a red herring for correctness — see Diagnosis §B.

**Goal:** Make the `--te-aligner blat --genome-aligner bwaaln` combination in RelocaTE3 reproduce RelocaTE2's TE-read recall, by matching RelocaTE2's BLAT sensitivity parameters (and confirming nothing else in the BLAT→trim→bwa-aln path diverges). Achieve R3 ≈ R2 for this aligner combo on the multi-TE and divergence benchmark panels before further R3 development.

**Non-goal:** Do not tune the classifier, TSD logic, or the genome (`bwaaln`) backend — that stage is already at parity. Do not add fastq chunking unless §B's verification proves it changes results (it should not).

**Tech stack:** Python 3.12, pysam, BLAT (unpinned, provided via `callers/relocate3/blat-env` in the benchmark), pytest, RelocaTE3 `src/RelocaTE3/aligners.py`.

---

## Diagnosis

RelocaTE2's default caller path (`--aligner blat`) is: split reads → BLAT reads against the TE library → trim TE portion off junction reads → place trimmed flanks on the genome with `bwa aln` → cluster into insertions. RelocaTE3's `blat` TE backend + `bwaaln` genome backend is meant to reproduce this. Direct source comparison found three differences; only §A is an open correctness gap.

### A. BLAT sensitivity parameters — PRIMARY GAP (undocumented)

RelocaTE2 (`RelocaTE2/scripts/relocaTE2.py:545`):
```python
blatcmd = '%s -minScore=10 -tileSize=7 %s %s %s ...' % (blat, te_fasta, fa, blatout, ...)
```

RelocaTE3 (`src/RelocaTE3/aligners.py:513-516`, `BlatBackend._blat_side`):
```python
subprocess.run(
    ["blat", str(te_library), query_fa, "-noHead", "-out=psl", psl],
    check=True,
)
```

RelocaTE3 passes **no** `-minScore`/`-tileSize`, so BLAT runs at its defaults (verified from `blat` help): **`minScore=30`, `tileSize=11`** (and `minIdentity=90`). RelocaTE2 lowers these to **`minScore=10`, `tileSize=7`**.

Consequences:
- **`tileSize=11 → 7`**: BLAT only seeds an alignment when it finds perfect tiles of the tile size. An 11-mer perfect seed is far rarer than a 7-mer in a short, mismatch-bearing junction read, so R3 fails to *seed* many TE-junction reads that R2 finds.
- **`minScore=30 → 10`**: even when seeded, R3 discards any hit scoring < 30 (matches − mismatches − gap penalty). Short junction overlaps (a read that is mostly genomic flank with only ~20–30 bp of TE) and divergent copies score below 30 and are dropped.
- `psl_to_sam` (`aligners.py:428`) applies **no** additional score filter — every BLAT-emitted PSL row becomes a SAM record — so the entire sensitivity difference lives in the BLAT command line. This makes the fix self-contained.

This is expected to be the dominant driver of any `blat-bwaaln` recall gap, and it is most severe on the **divergence panel** (reads drawn from TE copies 2–20 % diverged from the canonical library BLAT searches): fewer, shorter perfect tiles per read.

### B. Chunking / `--split` — NOT a correctness difference (verify, don't port)

RelocaTE2 optionally splits each fastq into 200,000-read chunks (`relocaTE2.py:91`, `fastq_split.pl`) and BLATs each chunk; the benchmark runs it **with** `--split`. RelocaTE3 does not chunk. This is a parallelization strategy, not a sensitivity lever:
- BLAT scores each *query read* independently against the same fixed TE-library database; the query set's size does not change any single read's alignment or score.
- RelocaTE2's own usage note (`relocaTE2.py:27`) says the *maximum-sensitivity* mode is **"do not split"** + multi-CPU `pblat`. Splitting is the fallback when a single run is too large.

Therefore chunking should be reproduced only in the sense of "no chunking = the sensitive mode." Task 2 verifies this empirically rather than assuming it, because the user specifically flagged chunking.

### C. Trim flank-length filter — separate, already-planned gap

`plans/2026-07-02-trim-recall-parity.md` documents that R3's trim step drops junction reads whose trimmed flank fails `len(trimmed_seq) >= len_cutoff_l`, even when the TE aligner *did* find them. That gap is orthogonal to §A (it acts *after* a read is found) but interacts with it: lowering BLAT's thresholds surfaces more candidate junction reads that then must pass the same trim gates as R2. Do not re-solve it here; Task 3 measures the residual after §A and confirms whether the 2026-07-02 plan still needs execution.

---

## Task 1: Parameterize BlatBackend with RelocaTE2's sensitivity flags

**Files:**
- Modify: `src/RelocaTE3/aligners.py` (`BlatBackend`, ~line 488-546)
- Test: `tests/aligners_test.py` (already covers `BlatBackend` / `psl_to_sam`)

**Step 1 — Write the failing test.** Assert the BLAT command `BlatBackend` builds contains `-minScore=10` and `-tileSize=7`. Prefer refactoring the command into a small builder method (e.g. `_blat_cmd(te_library, query_fa, psl)`) so the test asserts on the returned list without running BLAT:
```python
def test_blat_command_uses_relocate2_sensitivity_params():
    backend = BlatBackend()
    cmd = backend._blat_cmd("te.fa", "query.fa", "aln.psl")
    assert "-minScore=10" in cmd
    assert "-tileSize=7" in cmd
    assert cmd[:2] == ["blat", "te.fa"]  # db then query, unchanged orientation
```

**Step 2 — Run it, watch it fail** (`_blat_cmd` does not exist / flags absent).

**Step 3 — Implement.** Extract the command into `_blat_cmd` and add the two flags. Keep db-then-query order and `-noHead -out=psl`:
```python
def _blat_cmd(self, te_library, query_fa, psl):
    # RelocaTE2 parity: sensitize BLAT for short/divergent TE-junction reads.
    # Defaults (minScore=30, tileSize=11) miss reads with only a short or
    # mismatch-bearing TE overlap; RelocaTE2 uses -minScore=10 -tileSize=7.
    return [
        "blat", str(te_library), str(query_fa),
        "-minScore=10", "-tileSize=7", "-noHead", "-out=psl", str(psl),
    ]
```
Call it from `_blat_side`.

**Step 4 — Run tests, confirm green** (new test + existing BLAT/aligner tests unchanged).

**Step 5 — Commit** `feat(blat): sensitize BLAT to RelocaTE2 params (-minScore=10 -tileSize=7)`.

## Task 2: Empirically confirm chunking is not a sensitivity factor

**Goal:** Refute (or, if wrong, quantify) the chunking hypothesis so we can justify *not* adding fastq splitting.

**Steps (validation, not production code):**
1. Take one divergence-panel read set (e.g. `div000_rep01/reads/cov30x_rep1`, the same sample the smoke run used).
2. Run `BlatBackend._blat_side` (post-Task-1) on the whole read file → count PSL hit rows and distinct query reads with ≥1 hit.
3. Split the same reads into 200k chunks (`split`/`seqtk`), BLAT each chunk with the identical Task-1 command, concatenate PSL → count the same two metrics.
4. Assert the distinct-read-with-hit sets are **identical** (order aside). Record the counts in `validation/` and in this plan's "Validation record".

**Expected:** identical read sets → confirms §B, no chunking needed. If they differ, stop and investigate BLAT `-repMatch`/`.ooc` query-size effects before proceeding.

**Commit:** the validation script + a short note under `validation/`.

## Task 3: End-to-end blat-bwaaln parity measurement (R3 vs R2)

**Files:** benchmark repo (`relocate-benchmark`), not R3 source. No code change expected unless a residual gap points back at trim (§C).

**Steps:**
1. In the benchmark, run the `relocate3-blat-bwaaln` variant and `relocate2` on a shared slice: multi-TE panel (0 % divergence) at 5×/15×/30×, plus divergence panel at 0 % and (for the divergence-specific view) 2 %/5 %.
2. Compare pooled detection recall and per-TE-group recall R3 vs R2 for the blat-bwaaln combo.
3. Acceptance: R3 `blat-bwaaln` detection recall ≥ 0.95 × R2 on the 0 % multi-TE slice (parity band), and the R3−R2 gap on divergence 0/2/5 % is within noise of the multi-TE result.
4. If a residual gap remains after Task 1, trace one disagreement site (reuse the method in `plans/2026-07-02-trim-recall-parity.md` "Diagnosis") to attribute it to trim (§C) vs another cause, and hand off to that plan.

**Commit/record:** parity numbers in this plan's "Validation record" and in the benchmark docs.

---

## Acceptance criteria

- [ ] `BlatBackend` emits `-minScore=10 -tileSize=7`; unit test guards it.
- [ ] Chunking verified irrelevant to BLAT hit sets (Task 2 record).
- [ ] R3 `blat-bwaaln` detection recall within the parity band of R2 on the 0 % multi-TE slice (Task 3).
- [ ] Any residual gap is attributed (trim §C or documented new cause), not left unexplained.

## Open questions / risks

- **minIdentity:** both R2 and R3 leave BLAT's `minIdentity=90` default, capping usable divergence near ~10 %. This is a *shared* limit (not a parity gap) and explains why *both* callers collapse by 20 % divergence. Matching R2 is the goal here; raising minIdentity is a separate future experiment, not part of parity.
- **Runtime:** `tileSize=7` + `minScore=10` makes BLAT slower and emits more hits. Confirm the benchmark's blat variant walltime stays within its SLURM budget; if not, that is the real (performance) motivation for optional chunking — revisit §B as a *speed* feature only.
- **psl_to_sam volume:** more low-scoring hits flow into `psl_to_sam`/trim; confirm no downstream assumption breaks on multi-hit reads (existing multi-copy handling should cover it).

## References

- RelocaTE2 BLAT command: `RelocaTE2/scripts/relocaTE2.py:545`
- RelocaTE2 bwa-aln genome placement: `RelocaTE2/scripts/relocaTE_align.py:269-284`
- RelocaTE2 split: `relocaTE2.py:72-95` (`split_fq`), `fastq_split.pl`
- RelocaTE3 BLAT backend: `src/RelocaTE3/aligners.py:488-553`
- RelocaTE3 bwaaln backend (already parity): `src/RelocaTE3/aligners.py:272-336`
- Related trim gap: `plans/2026-07-02-trim-recall-parity.md`
