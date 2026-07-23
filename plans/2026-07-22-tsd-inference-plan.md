# Variable-Length TSD Inference Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Port RelocaTE2's read-depth TSD inference into RelocaTE3 so `--tsd UNK` yields variable-length TSD calls that match RT2, replacing RT3's fixed-length-only behavior.

**Architecture:** Add a TSD-unknown path to `InsertionFinder` (`src/RelocaTE3/insertions.py`) that infers TSD length from junction-cluster overlap (RT2 `TSD_len_calculate`) with a read-depth cross-check (RT2 `tsd_finder`), extracts the TSD sequence at the inferred span, and feeds the existing `te_insertions[event][tsd_start][tsd_seq]` structure. The fixed-motif path is untouched; `UNK` selects inference.

**Tech Stack:** Python 3.10+, pysam, pytest (unittest style, `tests/*_test.py`). Env: `pixi run --manifest-path pixi.toml test`. Reference implementation: `references/RelocaTE2/scripts/relocaTE_insertionFinder.py` (Python 2). Design: `plans/2026-07-22-tsd-inference-design.md`.

**Conventions:** Tests are `tests/<name>_test.py` (unittest.TestCase). Run one with `pixi run --manifest-path pixi.toml pytest tests/insertions_tsd_infer_test.py -v`. Commit after each task. **This is a faithful PORT** — cite the RT2 line you are porting in each commit/comment; validate against RT2 output, do not reinvent.

---

## Task 0: Capture RT2 ground truth (reference behavior)

Before porting, record exactly what RT2 produces so later tasks can assert parity.

**Files:** Create `plans/notes/rt2-tsd-depth-trace.md` (scratch notes, committed).

**Steps:**
1. Read `references/RelocaTE2/scripts/relocaTE_insertionFinder.py` in full for the
   `UNK` path: `main()` (L1726) → `find_insertion_cluster_bam` (L1335) →
   `TSD_from_read_depth` (L551) → `TSD_len_calculate` (L852) / `tsd_finder` (L843).
   Write down, precisely: (a) how left/right junction reads are identified (the
   `start:5|3` / `end:5|3` suffix + strand logic), (b) the exact boundary
   coordinates pushed into `TSD_left` / `TSD_right`, (c) the length formula
   (`right_max - left_max + 1`), (d) **where the TSD sequence comes from** (read
   bases vs reference bases — grep the depth path for genome/`seq[` slicing), (e)
   the `tsd_depth` threshold default.
2. Record answers in the notes file. The (d) answer decides whether Task 4 needs a
   genome FASTA argument.
3. Commit: `git add plans/notes/rt2-tsd-depth-trace.md && git commit -m "docs(tsd): trace RT2 depth-mode TSD inference for the port"`.

---

## Task 1: Port `TSD_len_calculate` as a pure function

**Files:**
- Modify: `src/RelocaTE3/insertions.py`
- Test: `tests/insertions_tsd_infer_test.py` (create)

**Step 1: Write the failing test.** Model the junction reads as small records
(name with `:start:5`/`:end:3` suffix, `start`, `end`, `strand`) placed so the
left/right boundaries overlap by a known length.

```python
import unittest
from RelocaTE3.insertions import InsertionFinder, _infer_tsd_length

class TestTSDLength(unittest.TestCase):
    def _reads(self, left_end, right_start):
        # right-boundary reads: name '...:start:5', strand '+', start=right_start
        # left-boundary reads:  name '...:end:3',   strand '+', end=left_end
        return [
            {"name": "rA:start:5", "start": right_start, "end": right_start + 40, "strand": "+"},
            {"name": "rB:end:3",   "start": left_end - 40, "end": left_end,       "strand": "+"},
        ]

    def test_tsd_len_overlap_5bp(self):
        # right cluster starts at 100, left cluster ends at 104 -> overlap 5 bp
        n = _infer_tsd_length(self._reads(left_end=104, right_start=100))
        self.assertEqual(n, 5)

    def test_tsd_len_overlap_3bp(self):
        n = _infer_tsd_length(self._reads(left_end=102, right_start=100))
        self.assertEqual(n, 3)
```

**Step 2: Run to verify it fails** (`ImportError: _infer_tsd_length`).
Run: `pixi run --manifest-path pixi.toml pytest tests/insertions_tsd_infer_test.py -k tsd_len -v`

**Step 3: Implement** `_infer_tsd_length(reads)` in `insertions.py`, porting
`TSD_len_calculate` (RT2 L852): build `tsd_left` (start coords of `:start` reads,
per strand rules) and `tsd_right` (end-1 coords of `:end` reads), take the
max-support coord on each side, return `right_max - left_max + 1` (guard: 0 if
either side empty). Match the strand handling from RT2 L874–893 exactly.

**Step 4: Run to verify it passes.**

**Step 5: Commit** `feat(tsd): port RT2 TSD_len_calculate (overlap length inference)`.

---

## Task 2: Port `tsd_finder` (read-depth boundary cross-check)

**Files:** Modify `src/RelocaTE3/insertions.py`; Test `tests/insertions_tsd_infer_test.py`.

**Step 1: Failing test** — synthetic per-position depth map where positions with
depth ≥ `tsd_depth × read_total` span a known length; assert `_tsd_finder(...)`
returns that length. (Port RT2 `tsd_finder` L843–850: count positions meeting the
depth threshold.)

**Step 2–4:** fail → implement → pass, citing RT2 L843.

**Step 5: Commit** `feat(tsd): port RT2 tsd_finder (depth-based TSD length)`.

---

## Task 3: TSD sequence extraction at the inferred span

**Files:** Modify `src/RelocaTE3/insertions.py`; Test `tests/insertions_tsd_infer_test.py`.
Depends on Task 0(d): source = reads or reference.

**Step 1: Failing test** — given the inferred `(tsd_start, tsd_len)` and the
cluster's reads (or a small reference), assert the extracted `tsd_seq` equals the
expected substring for a 4-bp and a 5-bp case.

**Step 2–4:** fail → implement `_extract_tsd_seq(...)` faithful to RT2 → pass.
If Task 0(d) shows reference bases are required, add a `genome_fasta` parameter to
`_extract_tsd_seq` here (the plumbing lands in Task 4).

**Step 5: Commit** `feat(tsd): extract inferred-length TSD sequence`.

---

## Task 4: Wire the `UNK` path into `find_insertions`

**Files:**
- Modify: `src/RelocaTE3/insertions.py` (`find_insertions`, `_cluster_reads`/`_assign_cluster` as needed)
- Modify (only if Task 0(d) requires genome): `src/RelocaTE3/insertions.py` signature + `src/RelocaTE3/cli.py::cmd_find_insertions` + `src/RelocaTE3/pipeline.py` caller
- Test: `tests/insertions_tsd_infer_test.py`

**Step 1: Failing test** — call `InsertionFinder().find_insertions(..., tsd="UNK", ...)`
on a tiny hand-built genome-aligned BAM (fixture) containing junction reads for one
insertion with a 5-bp TSD; assert the written table has that insertion with the
correct 5-bp `TSD` field (today this raises `NotImplementedError`).

**Step 2: Run to verify it fails** (NotImplementedError).

**Step 3: Implement.** Replace the `NotImplementedError` block (L71-80) with a call
into a new `_infer_insertions_depth_mode(...)` that walks clusters (reuse
`_cluster_reads`), and for each candidate site calls `_infer_tsd_length` (+
`_tsd_finder` cross-check), `_extract_tsd_seq`, and records into
`te_insertions[event][tsd_start][tsd_seq]` — the same structure `_write_output`
consumes. Port the pairing/orchestration from RT2 `TSD_from_read_depth` (L551) and
`find_insertion_cluster_bam` (L1335). If a genome is required, thread it through
(see Files).

**Step 4: Run to verify it passes.**

**Step 5: Commit** `feat(tsd): UNK read-depth inference path in find_insertions`.

---

## Task 5: Parity test against RelocaTE2

**Files:** Test `tests/insertions_tsd_parity_infer_test.py` (create). (There is a
precedent: `tests/insertions_tsd_parity_test.py`.)

**Step 1:** Build (or reuse) a small genome-aligned junction-read BAM with several
insertions of TSD length 3/4/5. Run the RT2 reference `relocaTE_insertionFinder.py`
UNK path on it once, capture the `all_nonref_insert.txt` TSD column as the golden
expectation (store as a fixture). Assert RT3 `--tsd UNK` produces the same TSD
sequence per site.

**Step 2:** Mark `@unittest.skipUnless` if RT2/py2 isn't runnable in CI; keep the
captured golden fixture so the test runs from RT3 alone.

**Step 3: Commit** `test(tsd): RT2-parity for UNK TSD inference`.

---

## Task 6: Full suite + docs

**Steps:**
1. `pixi run --manifest-path pixi.toml test` — all green, no regressions in the
   fixed-motif path (`insertions_test.py`, `insertions_tsd_parity_test.py`).
2. Update `README.md` / `plans/FEATURES.md`: `--tsd UNK` now infers variable-length
   TSDs (default-style behavior matching RT2).
3. Commit `docs: document --tsd UNK variable-length inference`.
4. Finish the branch (superpowers:finishing-a-development-branch) → PR; user pushes/merges.

---

## Follow-up (separate, downstream — NOT in this plan)

After merge, in relocate-benchmark: set `tsd = "UNK"` for the RelocaTE3 callers in
`config/benchmark.toml`, re-pin `callers/relocate3/pixi.toml` to the new RT3 rev,
clear + rerun the RT3 variants, and confirm TSD-exact rises from ~0.30 toward RT2's
~0.98.

---

## Notes / risks

- **Faithful port, not reinvention.** RT2's depth path (L551–910) is Python-2 and
  subtle (strand/boundary off-by-ones). Task 0 + the Task 5 parity test are the
  guardrails — port line-by-line and compare output.
- **Genome dependency:** decided in Task 0(d). If needed, it is a small signature
  addition to `find_insertions` + its two callers; the fixed-motif path keeps its
  current signature (genome optional / unused there).
- **`_merge_offset_starts`** exists to collapse fixed-wildcard splits; confirm it
  still behaves for variable-length inferred starts (it should be a no-op when the
  inference emits a single correct start per site).
