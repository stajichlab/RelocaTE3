# TSD Depth-Inference & R2 Literal Parity Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Close the TSD-string gap between RelocaTE3 and legacy RelocaTE2 by (a) porting R2's read-depth TSD-length estimator into `insertions.py`, (b) capturing the TSD characters from a representative junction read (R2's literal behavior, not the forward-strand reference), and (c) removing the `MAX_TSD` cap so R3 mirrors R2's behavior even on long/odd TSDs.

**Architecture:**
- Extend the junction/cluster data model so reads carry their sequence into TSD inference.
- Replace the all-or-nothing `_fetch_tsd(genome, ...)` branch in `_make_insertion` with: depth-based length estimator → literal read-derived TSD capture → reference-genome fallback only when no read sequence is available.
- Same fix in the legacy class-based `InsertionFinder._emit` for `supporting_junction` calls.
- No threshold added to `acceptance_test.py`. TSD agreement stays a soft metric in `validation/real_rice/report/characterized/summary.txt`.

**Tech Stack:** Python ≥3.10, pysam, pytest (via `pixi run …`). Code lives in `src/RelocaTE3/insertions.py`; tests in `tests/insertions_test.py` (new fixtures alongside the existing ones).

**Reference (R2):** `references/RelocaTE2/scripts/relocaTE_insertionFinder.py`
- `TSD_from_read_depth` line 551 — orchestrates depth-based inference
- `tsd_finder` line 843 — depth-threshold length estimator (1.0, 0.8, 0.6 of cluster reads)
- `TSD_len_calculate` line 852 — position-based length cross-check
- `TSD_check_cluster` line 1249 — captures TSD chars from read with regex `^(.{TSD_len})` / `(.{TSD_len})$`
- Reconciliation rule: prefer position-based length when they disagree (line 814-816)

---

## Task 1: Establish failing TSD-parity test fixtures

**Files:**
- Test: `tests/insertions_tsd_parity_test.py` (new)
- Reference data: a minimal pysam-built BAM written in-test (no large binary fixtures)

**Step 1: Write a unit test for the depth-length estimator (does not exist yet)**

```python
# tests/insertions_tsd_parity_test.py
"""Tests for the depth-based TSD length/sequence inference (R2 parity)."""
from RelocaTE3.insertions import _estimate_tsd_length_from_depth


def test_depth_estimator_returns_zero_when_no_reads():
    assert _estimate_tsd_length_from_depth([], breakpoint=100) == 0


def test_depth_estimator_matches_overlap_width():
    # Five reads all overlapping positions 100..102 inclusive (3bp), then tapering.
    reads = [(95, 102), (96, 102), (97, 102), (98, 102), (100, 110)]
    # depth at 100,101,102 == 5; threshold 1.0 yields width 3
    assert _estimate_tsd_length_from_depth(reads, breakpoint=100) == 3
```

**Step 2: Run the test to confirm it fails (import error)**

Run: `pixi run pytest tests/insertions_tsd_parity_test.py -v`
Expected: `ImportError: cannot import name '_estimate_tsd_length_from_depth'`

**Step 3: Commit the failing test**

```bash
cd RelocaTE3
git add tests/insertions_tsd_parity_test.py
git commit -m "test(insertions): failing parity test for depth-based TSD inference"
```

---

## Task 2: Implement `_estimate_tsd_length_from_depth` (R2 `tsd_finder` port)

**Files:**
- Modify: `src/RelocaTE3/insertions.py` (add helper near `_fetch_tsd` ≈ line 655)
- Test: `tests/insertions_tsd_parity_test.py`

**Step 1: Add the helper**

```python
def _estimate_tsd_length_from_depth(
    spans: list[tuple[int, int]],
    breakpoint: int,
    thresholds: tuple[float, ...] = (1.0, 0.8, 0.6),
) -> int:
    """Estimate TSD length from read-depth overlap near ``breakpoint``.

    Port of RelocaTE2 ``tsd_finder`` (relocaTE_insertionFinder.py:843). Builds a
    per-base depth pileup from ``spans`` (1-based inclusive ``(start, end)``
    tuples), then for each fractional threshold (in order) counts contiguous
    positions whose depth >= ``threshold * len(spans)``. Returns the first
    non-zero length, or 0 if none qualify.
    """
    if not spans:
        return 0
    depth: dict[int, int] = {}
    for s, e in spans:
        for p in range(s, e + 1):
            depth[p] = depth.get(p, 0) + 1
    total = len(spans)
    for frac in thresholds:
        cutoff = frac * total
        length = sum(1 for p, d in depth.items() if d >= cutoff)
        if length:
            return length
    return 0
```

**Step 2: Run the test, confirm it passes**

Run: `pixi run pytest tests/insertions_tsd_parity_test.py -v`
Expected: both tests PASS.

**Step 3: Commit**

```bash
cd RelocaTE3
git add src/RelocaTE3/insertions.py tests/insertions_tsd_parity_test.py
git commit -m "feat(insertions): port RelocaTE2 tsd_finder as _estimate_tsd_length_from_depth"
```

---

## Task 3: Carry read sequences into the cluster/junction model

**Files:**
- Modify: `src/RelocaTE3/models.py` (add `seq` field to `JunctionObservation`)
- Modify: `src/RelocaTE3/insertions.py` (lines 491-550, `_Cluster` + `_stream_clusters`)

**Why:** R2's TSD-character capture reads the literal first/last `TSD_len` bases of the read sequence; that data is currently dropped in `_stream_clusters`.

**Step 1: Read the model**

Run: `grep -nE "class JunctionObservation|@dataclass" src/RelocaTE3/models.py`
Then read the matched dataclass to see existing fields.

**Step 2: Add `seq: str` field with default `""` so old callers/tests still work**

```python
# in models.py JunctionObservation dataclass
seq: str = ""
```

**Step 3: Plumb `seq` into junction and support records in `_stream_clusters`**

In `insertions.py:498-500` extend `_Cluster.support` to `(name, gstart, gend, strand, seq)`:

```python
# supporting reads: (name, gstart, gend, strand, seq)
self.support: list[tuple[str, int, int, str, str]] = []
```

In `_stream_clusters` (around line 527-548):

```python
seq = rec.query_sequence or ""
# ...
current.junctions.append(
    JunctionObservation(
        name, side, pos, strand, _te_family(read_repeat, name), te_end, seq
    )
)
# ...
current.support.append((name, gstart, gend, strand, seq))
```

**Step 4: Update `_count_support` (line 664-673) to unpack 5-tuple**

```python
for _name, gstart, gend, strand, _seq in cluster.support:
```

**Step 5: Update `_call_support_only` (line 687-690) to unpack 5-tuple**

```python
plus_ends = [gend for _n, _s, gend, strand, _seq in cluster.support if strand == "+"]
minus_starts = [
    gstart for _n, gstart, _e, strand, _seq in cluster.support if strand == "-"
]
```

**Step 6: Run the full existing test suite to confirm no regressions**

Run: `pixi run pytest tests/insertions_test.py tests/pipeline_test.py -v`
Expected: all green (the `seq` field defaults to `""` so old fixtures pass).

**Step 7: Commit**

```bash
cd RelocaTE3
git add src/RelocaTE3/models.py src/RelocaTE3/insertions.py
git commit -m "refactor(insertions): carry read sequence through cluster/junction model"
```

---

## Task 4: Implement `_capture_tsd_from_read` (R2 `TSD_check_cluster` port)

**Files:**
- Modify: `src/RelocaTE3/insertions.py`
- Test: `tests/insertions_tsd_parity_test.py`

**Step 1: Add failing test for character capture**

Append to `tests/insertions_tsd_parity_test.py`:

```python
from RelocaTE3.insertions import _capture_tsd_from_read


def test_capture_tsd_from_left_read_uses_last_bases():
    # left-junction read: TSD is the last `length` chars of the seq
    assert _capture_tsd_from_read("AAAAATAA", side="left", length=3) == "TAA"


def test_capture_tsd_from_right_read_uses_first_bases():
    # right-junction read: TSD is the first `length` chars of the seq
    assert _capture_tsd_from_read("TTAGGGCC", side="right", length=3) == "TTA"


def test_capture_tsd_short_seq_returns_empty():
    assert _capture_tsd_from_read("AT", side="right", length=3) == ""
```

**Step 2: Implement**

In `insertions.py` near the new depth helper:

```python
def _capture_tsd_from_read(seq: str, side: str, length: int) -> str:
    """Return the literal TSD characters from a junction read.

    Mirrors RelocaTE2 ``TSD_check_cluster`` (line 1249): a *right*-side junction
    read carries the TSD at the start of the read, a *left*-side read at the end.
    Returns ``""`` when the read is too short or ``length <= 0``.
    """
    if length <= 0 or len(seq) < length:
        return ""
    return (seq[:length] if side == "right" else seq[-length:]).upper()
```

**Step 3: Run, confirm pass**

Run: `pixi run pytest tests/insertions_tsd_parity_test.py -v`
Expected: PASS.

**Step 4: Commit**

```bash
cd RelocaTE3
git add src/RelocaTE3/insertions.py tests/insertions_tsd_parity_test.py
git commit -m "feat(insertions): port RelocaTE2 TSD_check_cluster as _capture_tsd_from_read"
```

---

## Task 5: Remove `MAX_TSD` cutoff and rewrite `_make_insertion` to use depth + read capture

**Files:**
- Modify: `src/RelocaTE3/insertions.py:454-457, 598-635`
- Test: `tests/insertions_tsd_parity_test.py`

**Step 1: Add failing integration test (`_make_insertion` directly)**

This builds a `_Cluster` by hand, calls `_make_insertion`, and asserts the TSD is captured from the read rather than UNK. Append:

```python
from unittest.mock import MagicMock

from RelocaTE3.insertions import _Cluster, _make_insertion
from RelocaTE3.models import JunctionObservation


def _mk_junction(name, side, pos, strand, seq, te_end="5"):
    return JunctionObservation(
        read_name=name, side=side, position=pos, te_orientation=strand,
        te_name="mPing", te_end=te_end, seq=seq,
    )


def test_make_insertion_one_sided_recovers_tsd_from_depth_and_read():
    chrom = "Chr1"
    # Three right-junction reads at position 100, all sequences start with "TTA"
    cluster = _Cluster(chrom)
    cluster.junctions = [
        _mk_junction(f"r{i}", "right", 100, "+", "TTAGGGCCAAAA") for i in range(3)
    ]
    # supporting reads with 3bp overlap at the breakpoint => depth length 3
    cluster.support = [
        ("s1", 98, 102, "+", "NNNNN"),
        ("s2", 98, 102, "+", "NNNNN"),
        ("s3", 98, 102, "+", "NNNNN"),
    ]
    cluster.extend(98, 110)

    genome = MagicMock()  # should NOT be consulted when read capture succeeds
    ins = _make_insertion(chrom, left_reads=[], right_reads=cluster.junctions,
                          genome=genome, cluster=cluster)
    assert ins.tsd == "TTA"
    genome.fetch.assert_not_called()
```

**Step 2: Run to confirm failure**

Run: `pixi run pytest tests/insertions_tsd_parity_test.py::test_make_insertion_one_sided_recovers_tsd_from_depth_and_read -v`
Expected: FAIL — currently returns `"UNK"`.

**Step 3: Remove `MAX_TSD` and add `cluster` parameter to `_make_insertion`**

- Delete the `MAX_TSD = 20` constant (line 456-457) and the matching `MAX_TSD` mention in the comment above.
- Change `_make_insertion(chrom, left_reads, right_reads, genome)` to `_make_insertion(chrom, left_reads, right_reads, genome, cluster)`.
- Replace the body (lines 612-623) with:

```python
    spans = [(s, e) for _n, s, e, _strand, _seq in cluster.support]

    if left_reads and right_reads:
        i_end = left_reads[0].position
        i_start = right_reads[0].position
        if i_end - i_start + 1 > 0:
            tsd_len = i_end - i_start + 1
        else:
            tsd_len = _estimate_tsd_length_from_depth(spans, i_start)
            i_start = i_end = min(i_start, i_end)
    else:
        present = left_reads or right_reads
        i_start = i_end = present[0].position
        tsd_len = _estimate_tsd_length_from_depth(spans, i_start)
        if tsd_len > 0:
            if right_reads:
                i_end = i_start + tsd_len - 1
            else:
                i_start = i_end - tsd_len + 1

    tsd = _resolve_tsd(left_reads, right_reads, chrom, i_start, i_end, tsd_len, genome)
```

**Step 4: Add the `_resolve_tsd` helper (read-first, genome-fallback)**

```python
def _resolve_tsd(
    left_reads: list[JunctionObservation],
    right_reads: list[JunctionObservation],
    chrom: str,
    i_start: int,
    i_end: int,
    tsd_len: int,
    genome: pysam.FastaFile,
) -> str:
    """Capture TSD from a junction read (R2 parity) or fall back to the genome."""
    if tsd_len <= 0:
        return "UNK"
    for obs in right_reads:
        captured = _capture_tsd_from_read(obs.seq, "right", tsd_len)
        if captured:
            return captured
    for obs in left_reads:
        captured = _capture_tsd_from_read(obs.seq, "left", tsd_len)
        if captured:
            return captured
    fetched = _fetch_tsd(genome, chrom, i_start, i_end)
    return fetched or "UNK"
```

**Step 5: Update `_call_insertions` (line 638-652) to pass `cluster` through**

```python
ins = _make_insertion(cluster.chrom, left_reads, right_reads, genome, cluster)
```

**Step 6: Run the new test, confirm PASS**

Run: `pixi run pytest tests/insertions_tsd_parity_test.py -v`
Expected: all PASS.

**Step 7: Run the full unit test suite for regressions**

Run: `pixi run pytest tests/insertions_test.py tests/pipeline_test.py tests/acceptance_test.py -v`
Expected: existing tests still pass. The acceptance test may legitimately shift TSD strings (now read-derived); the recall/precision gate is on position match, not TSD string, so it should be unaffected.

**Step 8: Commit**

```bash
cd RelocaTE3
git add src/RelocaTE3/insertions.py tests/insertions_tsd_parity_test.py
git commit -m "feat(insertions): R2-parity TSD inference (depth length + read capture, drop MAX_TSD)"
```

---

## Task 6: Apply the same parity fix to the legacy `InsertionFinder._emit`

**Files:**
- Modify: `src/RelocaTE3/insertions.py:386-432` (the `_emit` method)

**Why:** The class-based path emits `tsd_field = "supporting_junction"` (line 428) when no junction reads paired. R2 still resolves a TSD by depth in that case. Mirror Task 5's logic here so both code paths produce comparable output.

**Step 1: Trace `_emit` and verify it has access to the equivalent of `cluster.support`**

Run: `grep -nE "supporting|_emit|tsd_field" src/RelocaTE3/insertions.py | head -20`
Read lines 365-430 to confirm the data available at the call site.

**Step 2: Where `tsd_field = "supporting_junction"` is set, attempt depth-based recovery first**

Replace the relevant block (around line 426-428):

```python
elif te_insertions[event][tsd_start][found_tsd].get("supporting"):
    tsd_field = "supporting_junction"  # ← drop this literal
```

…with a call into `_estimate_tsd_length_from_depth` + `_capture_tsd_from_read` (mirroring Task 5). If both return empty, keep `"supporting_junction"` as the literal fallback so downstream consumers still see a sentinel.

**Step 3: Run unit tests**

Run: `pixi run pytest tests/insertions_test.py -v`
Expected: PASS.

**Step 4: Commit**

```bash
cd RelocaTE3
git add src/RelocaTE3/insertions.py
git commit -m "feat(insertions): apply depth-based TSD recovery in InsertionFinder._emit"
```

---

## Task 7: Update `plans/PLAN.md` and the soft validation metric

**Files:**
- Modify: `plans/PLAN.md:280-283`
- Modify: `validation/real_rice/compare_char.py` (add canonical TSD-match column)

**Step 1: Move the deferred TSD bullet**

In `plans/PLAN.md`, take the *"Still deferred (lower value): read-depth TSD length refinement (`tsd_finder`/`TSD_len_calculate` — currently TSD = the genomic overlap span)"* line out of the deferred list. Add a one-liner under Phase 3 ✅ noting depth-based length + literal read capture landed via `_estimate_tsd_length_from_depth` / `_capture_tsd_from_read` (insertions.py), `MAX_TSD` removed.

**Step 2: Report canonical TSD-match rate alongside strict match in `summary.txt`**

In `validation/real_rice/compare_char.py` (around the summary lines 379-401), add a `canonical_tsd_match` counter that treats `tsd == revcomp(tsd)` as a match, and print both rates. Keep all behavior soft — no test assertion.

**Step 3: Commit**

```bash
cd RelocaTE3
git add plans/PLAN.md validation/real_rice/compare_char.py
git commit -m "docs(plans): mark depth-based TSD inference complete; report canonical TSD match"
```

---

## Task 8: Full validation re-run (manual, user-initiated)

This task is **not auto-run** — the rice 10-sample validation is a SLURM job. Document the expected commands and outcomes; the user kicks it off.

**Commands (from `RelocaTE3/`):**

```bash
pixi run validate-rice --local B_10   # quick smoke
pixi run validate-rice                 # full 10-sample
```

**Expected delta in `report/characterized/tsd_confusion.tsv`:**

- `TTA/TTA` cell grows substantially.
- `TAA → UNK` cell shrinks toward 0 (now resolved via depth + read capture).
- `* → supporting_junction` cells shrink (Task 6).
- `TTA/TAA` may appear as a new cell where R3's read happened to be on the opposite strand from R2's — this is expected with literal capture; canonical-match rate should still be ≥ 0.95.

---

## Conventions

- One commit per task. No bundled refactors.
- Run `pre-commit run --all-files` before each commit; ruff/black/codespell must pass.
- All commands run from `RelocaTE3/` (the git repo) via `pixi run …`.
- No `mkdir` of result/log directories in scripts that aren't already creating them.
