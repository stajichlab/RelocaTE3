# Insertion Fragmentation Fix (TSD-start merge) Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Stop `InsertionFinder` from emitting one biological insertion as two adjacent single-sided calls (1–2 bp apart), which inflates false positives on the ground-truth benchmark.

**Architecture:** Add a per-cluster merge step in `InsertionFinder` that collapses `tsd_start` sub-buckets whose TSD coordinate spans overlap (`|Δstart| < captured_TSD_len`), pooling their junction reads/counts, before `_write_output`. Two split single-sided fragments recombine into one both-junction call. This restores RelocaTE2's effective behavior for the fixed-length wildcard TSD mode without porting full variable-length TSD support.

**Tech Stack:** Python 3.10+, pysam, pytest/unittest, pixi (pinned minimap2/samtools). Run everything under `pixi run`.

---

## Background (read before starting)

**This is the LIVE path.** The CLI `find-insertions` and the benchmark run
`InsertionFinder` (class) in `src/RelocaTE3/insertions.py`, whose emission logic
is `_write_output` → `_write_event_start`. Do **not** edit `_pair_breakpoints` /
`_call_insertions` / the module-level `find_insertions` — those belong to the
separate `pipeline.py` path and are not exercised by the benchmark. (This
corrects `plans/PERFORMANCE.md` Priority 2, which named the wrong function.)

**Root cause.** `_write_output` emits one row per distinct `tsd_start` in
`te_insertions[event][tsd_start][tsd_seq]`. For one insertion, the right junction
sets `tsd_start = start` (TSD left edge) while the left junction sets
`tsd_start = end - len(tsd)` using the **wildcard pattern width** (3 for
`tsd="..."`). When the true TSD is 4–5 bp, the two sides land `true_len - 3`
= 1–2 bp apart → two rows (one `R:n L:0`, one `R:0 L:m`). 99.7 % of benchmark
false positives are exactly this 1–2 bp split (measured: 691/693 across 9 samples).

**Why the merge is safe.** Distinct insertions cannot have overlapping TSDs, so
their `tsd_start` values are ≥ TSD-length apart and never merge. The merge window
is exactly the captured TSD length.

**Parity note.** Keep RelocaTE2's classifier semantics faithful (see memory
`feedback_r2_logic_preserving`). RelocaTE2 avoids this split by using the correct
TSD length; RelocaTE2's own multi-start arbitration
(`references/RelocaTE2/scripts/relocaTE_insertionFinder.py:258–333`) is a
*different* concern (multiple genuine starts in one cluster) and is deliberately
**out of scope** here — see "Deferred" at the end.

## Key data structures (already built by `_cluster_reads`)

```
te_insertions[event][tsd_start][tsd_seq]        # defaultdict -> {"count","left","right","+","-"} : int
te_insertions_reads[event][tsd_start][tsd_seq]  # defaultdict -> {"read","left_read","right_read"} : list[str]
```
`event` = int cluster id, `tsd_start` = int genomic coord, `tsd_seq` = captured TSD string.

## Insertion point

`InsertionFinder.find_insertions` (`src/RelocaTE3/insertions.py:99–121`): call the
new merge method **after** `_cluster_reads` populates the dicts and **before**
`_write_output`.

---

## Task 1: Add the failing unit test for the offset-split merge

**Files:**
- Test: `tests/insertions_test.py` (add a method to `class TestInsertionFinder`)

**Step 1: Write the failing test**

Reproduces the benchmark pattern: `tsd="..."` (3 bp wildcard), true 4 bp TSD
`GCAT` at Chr1:1000–1003. Right junction captures `GCA` at start 1000; left
junction captures `CAT` and lands at `end - 3 = 1001` (the 1 bp offset). Today
this yields two rows; after the fix it must yield one both-junction row.

```python
    def test_offset_split_merges_to_one_call(self):
        """A wildcard-TSD 1 bp offset between the two junction sides must not
        fragment one insertion into two single-sided calls (issue: TSD-start
        divergence)."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "flank.bam")
            # True 4 bp TSD "GCAT" at genomic 1000..1003; wildcard tsd="..." (3 bp).
            # Right (":start:5", +): 1-based start = 1000 -> captures "GCA", tsd_start=1000.
            # Left  (":end:5",  +): seq ends "CAT", ref_end+1 = 1004 -> tsd_start = 1001.
            right = {"name": "rR:start:5", "seq": "GCA" + "A" * 37, "start0": 999}
            left = {"name": "rL:end:5", "seq": "A" * 37 + "CAT", "start0": 963}
            _write_junction_bam(bam_path, "Chr1", 5000, [left, right])

            read_repeat = os.path.join(workdir, "read_repeat_name.txt")
            with open(read_repeat, "w") as fh:
                fh.write("rR\tmping\t+\n")
                fh.write("rL\tmping\t+\n")

            finder = InsertionFinder(mismatch_allow=0, min_mapq=1)
            out_txt = finder.find_insertions(
                bam_file=Path(bam_path),
                read_repeat_file=Path(read_repeat),
                tsd="...",
                target="Chr1",
                sample="HEG4",
                outdir=Path(workdir),
            )
            rows = [ln for ln in out_txt.read_text().splitlines() if ln.strip()]
            self.assertEqual(len(rows), 1, f"expected 1 merged row, got {rows}")
            cols = rows[0].split("\t")
            self.assertEqual(cols[6], "T:2")   # pooled total
            self.assertEqual(cols[7], "R:1")   # right junction retained
            self.assertEqual(cols[8], "L:1")   # left junction retained
```

**Step 2: Run it to confirm it fails**

Run: `pixi run pytest tests/insertions_test.py::TestInsertionFinder::test_offset_split_merges_to_one_call -v`
Expected: FAIL — currently emits 2 rows (`AssertionError: expected 1 merged row, got [...]`).

**Step 3: Commit the failing test**

```bash
git add tests/insertions_test.py
git commit -m "test(insertions): failing test for TSD-start offset fragmentation"
```

---

## Task 2: Add a no-over-merge guard test (distinct nearby sites stay separate)

**Files:**
- Test: `tests/insertions_test.py`

**Step 1: Write the test**

Two genuine insertions 50 bp apart (TSDs cannot overlap) must remain two rows.

```python
    def test_distinct_sites_not_merged(self):
        """Two real insertions farther apart than the TSD length must stay two
        separate calls (merge must not collapse distinct sites)."""
        with tempfile.TemporaryDirectory() as workdir:
            bam_path = os.path.join(workdir, "flank.bam")
            # Site A at tsd_start=1000, Site B at tsd_start=1050 (>3 bp apart).
            reads = [
                {"name": "aR:start:5", "seq": "TTA" + "A" * 37, "start0": 999},
                {"name": "aL:end:5", "seq": "A" * 37 + "TTA", "start0": 962},
                {"name": "bR:start:5", "seq": "TTA" + "A" * 37, "start0": 1049},
                {"name": "bL:end:5", "seq": "A" * 37 + "TTA", "start0": 1012},
            ]
            _write_junction_bam(bam_path, "Chr1", 5000, reads)
            read_repeat = os.path.join(workdir, "read_repeat_name.txt")
            with open(read_repeat, "w") as fh:
                for n in ("aR", "aL", "bR", "bL"):
                    fh.write(f"{n}\tmping\t+\n")
            finder = InsertionFinder(mismatch_allow=0, min_mapq=1)
            out_txt = finder.find_insertions(
                bam_file=Path(bam_path),
                read_repeat_file=Path(read_repeat),
                tsd="TTA",
                target="Chr1",
                sample="HEG4",
                outdir=Path(workdir),
            )
            rows = [ln for ln in out_txt.read_text().splitlines() if ln.strip()]
            self.assertEqual(len(rows), 2, f"expected 2 distinct rows, got {rows}")
```

**Step 2: Run it**

Run: `pixi run pytest tests/insertions_test.py::TestInsertionFinder::test_distinct_sites_not_merged -v`
Expected: PASS today (they are already separate) — this is a regression guard so the fix in Task 3 does not over-merge. If it errors for another reason, fix the test before proceeding.

**Step 3: Commit**

```bash
git add tests/insertions_test.py
git commit -m "test(insertions): guard against over-merging distinct nearby sites"
```

---

## Task 3: Implement `_merge_offset_starts` and wire it in

**Files:**
- Modify: `src/RelocaTE3/insertions.py` — add method to `InsertionFinder`; call it in `find_insertions` between `_cluster_reads` and `_write_output` (around line 108–113).

**Step 1: Add the merge method**

Add inside `class InsertionFinder` (e.g. just after `_cluster_reads`/`_assign_cluster`, before the "output" section). Complete code:

```python
    # ------------------------------------------------------------------
    # sub-cluster reconciliation
    # ------------------------------------------------------------------
    def _merge_offset_starts(self, te_insertions, te_insertions_reads):
        """Merge tsd_start sub-buckets split by the fixed-length wildcard TSD.

        The right junction records ``tsd_start = start`` (TSD left edge) while the
        left junction records ``tsd_start = end - len(tsd)``. When the true TSD is
        longer than the wildcard pattern the two sides diverge by
        ``true_len - wildcard_len`` (1-2 bp), fragmenting one insertion into two
        single-sided rows. Within each cluster, collapse tsd_start entries whose
        TSD coordinate spans overlap (gap < captured TSD length) into a single
        canonical start, pooling counts and reads. Distinct insertions have
        non-overlapping TSDs (>= TSD length apart) and are never merged.

        Mutates ``te_insertions`` and ``te_insertions_reads`` in place.
        """
        for event in list(te_insertions.keys()):
            starts = te_insertions[event]
            if len(starts) < 2:
                continue
            # captured TSD length for this cluster (wildcard captures are uniform);
            # fall back to 1 (no merge) if nothing was captured.
            tsd_len = 0
            for pos_map in starts.values():
                for seq in pos_map:
                    tsd_len = max(tsd_len, len(seq))
            if tsd_len < 2:
                continue

            ordered = sorted(starts.keys(), key=int)
            # chain consecutive starts whose spans overlap (gap < tsd_len)
            groups: list[list[int]] = []
            current = [ordered[0]]
            for s in ordered[1:]:
                if int(s) - int(current[-1]) < tsd_len:
                    current.append(s)
                else:
                    groups.append(current)
                    current = [s]
            groups.append(current)

            for group in groups:
                if len(group) < 2:
                    continue
                # canonical start = highest total count; tie -> smallest coordinate
                def _group_total(st):
                    return sum(b["count"] for b in starts[st].values())

                canonical = min(group, key=lambda st: (-_group_total(st), int(st)))
                for st in group:
                    if st == canonical:
                        continue
                    for seq, bucket in starts[st].items():
                        dst = te_insertions[event][canonical][seq]
                        for k, v in bucket.items():
                            dst[k] += v
                        rsrc = te_insertions_reads[event][st][seq]
                        rdst = te_insertions_reads[event][canonical][seq]
                        for k, v in rsrc.items():
                            rdst[k].extend(v)
                    del te_insertions[event][st]
                    if te_insertions_reads.get(event, {}).get(st) is not None:
                        del te_insertions_reads[event][st]
```

**Step 2: Wire it into `find_insertions`**

In `src/RelocaTE3/insertions.py`, after the `self._cluster_reads(...)` call
(ends ~line 108) and before building `result_dir`/`_write_output`, insert:

```python
        self._merge_offset_starts(te_insertions, te_insertions_reads)
```

**Step 3: Run the Task 1 test — expect PASS**

Run: `pixi run pytest tests/insertions_test.py::TestInsertionFinder::test_offset_split_merges_to_one_call -v`
Expected: PASS (1 row, `T:2 R:1 L:1`).

**Step 4: Run the Task 2 guard — expect PASS**

Run: `pixi run pytest tests/insertions_test.py::TestInsertionFinder::test_distinct_sites_not_merged -v`
Expected: PASS (2 rows).

**Step 5: Commit**

```bash
git add src/RelocaTE3/insertions.py
git commit -m "fix(insertions): merge wildcard-TSD offset-split starts into one call"
```

---

## Task 4: Run the full insertions + parity suite (no regressions)

**Files:** none (verification only).

**Step 1: Run the insertions-related tests**

Run:
```bash
pixi run pytest tests/insertions_test.py tests/insertions_tsd_parity_test.py \
  tests/insertions_tsd_class_parity_test.py -v
```
Expected: all PASS. In particular `test_known_tsd_call` (a correct-length TSD
where both sides already share `tsd_start`) must still emit exactly one row —
confirms the merge is a no-op when there is nothing to merge.

**Step 2: Run the full unit suite**

Run: `pixi run test`
Expected: all PASS (acceptance_test may be slow / shell out to minimap2+samtools;
it must not regress).

**Step 3: Commit (only if any test needed adjustment)**

```bash
git add -A && git commit -m "test(insertions): keep suite green after offset-merge fix"
```

---

## Task 5: Measure against the ground-truth benchmark

**Files:** none (measurement). Benchmark repo:
`/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/relocate_benchmark/relocate-benchmark`.

**Step 1: Sync the benchmark's RelocaTE3 env to this code**

RelocaTE3 is installed editable in the benchmark's pixi env, so edits are live,
but sync to be safe:
```bash
cd /bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/relocate_benchmark/relocate-benchmark
bash pipeline/update_relocate3.sh
```
Note the printed RelocaTE3 git commit — it must be the commit from Task 3.

**Step 2: Re-run RelocaTE3 on the 30x samples**

```bash
bash pipeline/submit_benchmark.sh --caller relocate3 --coverage 30
```
Wait for the SLURM array + dependent aggregation to finish (see README §Quickstart).

**Step 3: Check false positives collapsed**

```bash
column -t -s$'\t' reports/precision.tsv | grep relocate3
```
Expected: `relocate3 cov30x_rep*` `false_positive_calls` drop from ~127–142 to
≤ ~2; `overall_precision` rises from ~0.70 toward ~1.0. `matched_calls` must be
stable or higher (recall not harmed).

**Step 4: Spot-check the known locus**

```bash
awk -F'\t' '$5 ~ /^30817[0-9]/' runs/relocate3/cov30x_rep1/raw/results/ALL.mPing.all_nonref_insert.txt
```
Expected: a single row near 308175–308178 with both `R:` and `L:` > 0 (not two
single-sided rows).

**Step 5: Record results**

Update `plans/PERFORMANCE.md` Priority 2 with the before/after precision and FP
counts, and correct the function attribution (`InsertionFinder._write_output`,
not `_pair_breakpoints`). Add a `.living/learnings.md` entry per the repo's
post-action protocol. Commit:
```bash
git add plans/PERFORMANCE.md
git commit -m "docs(perf): record offset-merge precision fix results"
```

---

## Deferred (explicitly out of scope for this plan)

- **RelocaTE2 multi-start arbitration** (`relocaTE_insertionFinder.py:258–333`):
  drops weak single-sided starts when a cluster has multiple genuine starts.
  Not needed for the 99.7 % offset-split class and risks lowering recall; revisit
  only if post-fix precision is still below the truth-derived floor.
- **`TSD_from_read_depth` far-apart junction pairing** (PERFORMANCE.md Priority 5):
  handles the ~2 residual distant FPs (complementary junctions tens of bp apart).
  Separate, larger port.
- **Full variable-length TSD support**: the underlying `tsd="..."` limitation.
  This fix makes the wildcard mode behave correctly without it.

## Risks / watch-items

- **Canonical-coordinate choice** shifts the reported position by ≤ 2 bp; the
  benchmark match window is ±10 bp so recall is unaffected, and exact-TSD is
  already caveated for the 3 bp wildcard. Document the choice (highest count,
  tie → smallest coordinate).
- **Over-merge**: guarded by Task 2 and the `< tsd_len` window; distinct
  insertions are ≥ TSD-length apart.
- **`defaultdict` side effects**: `_group_total` and the merge read
  `starts[st]` / `te_insertions_reads[event][st]` — iterate over `list(...)`
  keys and delete after pooling (as written) to avoid mutating during iteration.
