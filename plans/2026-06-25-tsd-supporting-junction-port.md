# Close the Remaining `supporting_junction` Gap (Single-Sided Junctions)

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.
> **For a future maintainer (or a fresh Claude session):** Read the "Where this fits" section in full before touching code. The prior two plans (`2026-06-24-tsd-depth-inference.md`, `2026-06-25-tsd-class-path-port.md`) closed most of the TSD-string gap with RelocaTE2 but left a structural ~16% mismatch caused by the legacy class path discarding read-captured TSDs whenever junction reads land on only one side of the insertion.

**Goal:** Drop the remaining ~960 `TTA/TAA → supporting_junction` cells in the rice validation TSD confusion matrix toward zero. Two options exist (B1 surgical, B2 faithful R2 parity); plan executes **B1** by default with B2 sketched as the alternative.

---

## Where this fits

### Current state (as of 2026-06-25, branch `fix-tsd-bug`)
- Function-based `find_insertions` (`src/RelocaTE3/insertions.py:888+`, used by `pipeline.run_sample` / `cli.py`) has full depth-mode TSD inference (`_estimate_tsd_length_from_depth`, `_capture_tsd_from_read`, `_resolve_tsd`, rewritten `_make_insertion`). See `plans/2026-06-24-tsd-depth-inference.md`.
- Legacy `InsertionFinder` class (`src/RelocaTE3/insertions.py:27+`, used by `__main__.py`'s `find-insertions` subcommand and therefore by the validation harness) gets read-captured TSDs via a wildcard regex (`tsd = "..."` in `validation/real_rice/config.example.toml`). See `plans/2026-06-25-tsd-class-path-port.md`.
- Rice 10-sample validation: strict TSD agreement 80.4% (4088/5086). The remaining ~16% are dominated by R2-`TTA`/`TAA` ↔ R3-`supporting_junction` cells: **524** + **434** + 8 + 5 + 3 ≈ **974 mismatches**.

### Why those rows are `supporting_junction` in R3

Spot-check of `validation/real_rice/results/B_10/results/ALL.mping.all_nonref_insert.txt` after the wildcard run:

```
mPing  supporting_junction  B_10  Chr1  3487911..3487913  -  T:2  R:2  L:0  ST:0  SR:0  SL:0
mPing  supporting_junction  B_10  Chr1  4803952..4803954  -  T:7  R:0  L:7  ST:0  SR:0  SL:0
mPing  supporting_junction  B_10  Chr1  4928689..4928691  +  T:6  R:0  L:6  ST:0  SR:0  SL:0
...
```

`R:N L:0` or `R:0 L:N` — every row has junction reads on exactly one side. `top_tsd` IS a real read-captured 3-mer (`TTA`, `TAA`, etc., because of the wildcard fix). The class path's `_write_event_start` discards it:

```python
# src/RelocaTE3/insertions.py:422-428
if left_count > 0 and right_count > 0:
    tsd_field = top_tsd
elif total_count == 1:
    tsd_field = "singleton"
else:
    tsd_field = "supporting_junction"   # ← single-sided lands here, top_tsd lost
```

### Why R2 emits `TTA`/`TAA` for the same sites

`references/RelocaTE2/scripts/relocaTE_insertionFinder.py:376-388` has the *same* gate: single-sided junctions are written as `supporting_junction`. **But R2 also runs `TSD_from_read_depth` (line 551)** which, for each cluster with one-sided junctions, walks the read-depth pileup of *all* reads in the cluster (junction + supporting + low-quality), splits the cluster into sub-clusters, and calls `TSD_check_cluster` (line 1249) on those reads with a `.{TSD_len}` pattern. The synthetic registrations land in `teInsertions[event][TSD_start][TSD_seq]` *as if they were junction reads*, populating the missing side. After the depth pass the both-sided branch fires and emits the real TSD. Only when the depth pass *also* fails does R2 write the `supporting_junction` sentinel — which is why the matrix still has 42 `supporting_junction/supporting_junction` agreement cells.

### Two paths to close the gap

| | **B1 — trust the wildcard-captured single-sided TSD** | **B2 — port R2's `TSD_from_read_depth`** |
|---|---|---|
| Code touched | 1 hunk in `_write_event_start` (~5 LOC) | New `_synth_missing_side_from_depth` helper + new pass over closed clusters (~150 LOC) |
| Faithfulness | Closer to R2's *output* than R2's *logic* (R2 still labels these `supporting_junction` when depth synthesis fails) | Faithful R2 port; matches R2's exact emission decisions |
| Expected impact on matrix | The 974 mismatch cells collapse to `TTA/TTA` + `TAA/TAA`; strict agreement → ≥ 95% | Same, plus the 42 `supporting_junction/supporting_junction` cells stay aligned in edge cases where both tools fail to recover |
| Risk | `top_tsd` confidence is lower on single-sided sites; mis-captured 3-mers (e.g. soft-clipped reads) could leak in | Larger surface area, two-pass refactor of `_cluster_reads` lifecycle |

**Plan executes B1.** B2 is sketched in the "If B1 is insufficient" section at the bottom — open a new plan if/when the ~5% residual matters.

---

## B1 — Emit `top_tsd` for single-sided junctions when wildcard mode is active

### Task 1: Failing test for the single-sided emission

**Files:**
- Test: `tests/insertions_tsd_class_parity_test.py` (extend the existing file)

**Step 1: Add a fixture with only a right-junction read carrying TAA**

Append to `tests/insertions_tsd_class_parity_test.py`:

```python
def _write_single_sided_bam(tmp_path: Path) -> Path:
    """One right-junction read framing a TAA TSD at Chr1:100..102; no left side."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"LN": 1000, "SN": "Chr1"}],
    }
    raw = tmp_path / "syn1s.raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        # two right-junction reads so total_count > 1 (avoid the singleton branch)
        for i in range(2):
            r = pysam.AlignedSegment(bam.header)
            r.query_name = f"r{i}:start:5"
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 99  # 1-based 100
            r.mapping_quality = 60
            r.cigartuples = [(0, 10)]
            r.query_sequence = "TAAGGGCCAA"
            r.query_qualities = pysam.qualitystring_to_array("I" * 10)
            r.set_tag("NM", 0)
            bam.write(r)
    sorted_bam = tmp_path / "syn1s.bam"
    pysam.sort("-o", str(sorted_bam), str(raw))
    pysam.index(str(sorted_bam))
    return sorted_bam


def test_class_path_emits_captured_tsd_for_single_sided_junction(tmp_path):
    """One-sided junction with wildcard tsd should emit top_tsd, not 'supporting_junction'."""
    bam = _write_single_sided_bam(tmp_path)
    repeat = tmp_path / "read_repeat.txt"
    repeat.write_text("r0\tmPing\t+\nr1\tmPing\t+\n")
    outdir = tmp_path / "out"
    finder = InsertionFinder(mismatch_allow=2, min_mapq=1)
    out_txt = finder.find_insertions(
        bam_file=bam,
        read_repeat_file=repeat,
        tsd="...",
        target="Chr1",
        sample="syn",
        outdir=outdir,
        te_name="mPing",
    )
    rows = [
        line.split("\t")
        for line in Path(out_txt).read_text().splitlines()
        if line and not line.lower().startswith("strain")
    ]
    assert len(rows) == 1, rows
    # _emit format: TE, TSD, sample, chrom, coord, strand, T:N, R:N, L:N, ST, SR, SL
    assert rows[0][0] == "mPing"
    assert rows[0][1] == "TAA", rows  # currently fails: emits "supporting_junction"
    assert rows[0][7] == "R:2"
    assert rows[0][8] == "L:0"
```

**Step 2: Run, confirm failure**

```
pixi run pytest tests/insertions_tsd_class_parity_test.py::test_class_path_emits_captured_tsd_for_single_sided_junction -v
```

Expected: FAIL on `assert rows[0][1] == "TAA"`, actual `"supporting_junction"`.

**Step 3: Commit**

```bash
cd RelocaTE3
git add tests/insertions_tsd_class_parity_test.py plans/2026-06-25-tsd-supporting-junction-port.md
git commit -m "test(insertions): failing test for single-sided wildcard TSD emission"
```

### Task 2: Implement — emit `top_tsd` whenever it's a real read-captured string

**Files:**
- Modify: `src/RelocaTE3/insertions.py:422-428` (the `_write_event_start` gate)

**Step 1: Read the current logic and the read context**

The relevant block is the one identified above. `top_tsd` is `max(tsd_count.items(), key=lambda kv: kv[1])[0]` (line 414) and is only ever `"UNK"` when the caller passed a literal TSD motif (e.g. `"TTA"`) AND no junction read matched it. Under the harness config (`tsd = "..."` since `plans/2026-06-25-tsd-class-path-port.md`), `top_tsd` is always a captured 3-mer.

**Step 2: Replace the gate with a "trust real captures" rule**

```python
# Emit the read-captured TSD whenever one was inferred. Mirrors RelocaTE2's
# behavior post-`TSD_from_read_depth`: that pass synthesizes the missing-side
# junction reads from the depth pileup, lifting single-sided clusters into
# the both-sided emission. We achieve the same emission directly when the
# wildcard TSD mode (validation/real_rice/config.example.toml -> tsd="...")
# has filled top_tsd with a real read-captured 3-mer.
real_capture = top_tsd and top_tsd not in {"UNK", "UKN"}
if real_capture and (left_count > 0 and right_count > 0):
    tsd_field = top_tsd
elif real_capture and total_count > 1:
    # single-sided junction with multiple reads; trust the wildcard capture.
    tsd_field = top_tsd
elif total_count == 1:
    tsd_field = "singleton"
else:
    tsd_field = "supporting_junction"
```

**Why two separate `real_capture` branches and not a flat `if real_capture: ...`?**
Keep the existing both-sided emission intact (familiar for future readers), make the new single-sided branch explicit so its `total_count > 1` guard reads as a deliberate "don't promote singletons" rule.

**Why not collapse to my reverted Task 6 wording?**
The earlier revert (`f241599`) was done under the *literal-TSD* config (`tsd = "TTA"`); back then `top_tsd == "UNK"` was common and the change misrouted those calls to `supporting_junction`. With the wildcard config now in place, `top_tsd == "UNK"` only happens if the caller explicitly opted into literal-TSD mode — the `real_capture` guard makes that case keep the old behavior.

**Step 3: Run the new test, confirm PASS**

```
pixi run pytest tests/insertions_tsd_class_parity_test.py -v
```

Expected: both tests in the file PASS.

**Step 4: Full regression sweep**

```
pixi run pytest -q
```

Expected: green. Pay special attention to `tests/insertions_test.py::TestInsertionFinder` and `tests/main_test.py` — those exercise the legacy class path under the literal-TSD config and must not regress.

**Step 5: Commit**

```bash
cd RelocaTE3
git add src/RelocaTE3/insertions.py
git commit -m "feat(insertions): legacy _emit keeps captured TSD for single-sided junctions

Only fires when the captured top_tsd is a real read-derived string (i.e. the
caller is in wildcard TSD mode); literal-TSD callers keep the old behavior
where 'UNK' captures route to 'supporting_junction'."
```

### Task 3: Validate on B_10 first (no SLURM)

**Files:** none modified. Runtime check.

**Step 1: Delete the stale characterize sentinel for B_10**

```bash
find validation/real_rice/results/B_10 -name '*.all_nonref_insert.characTErized.txt' -delete
```

**Step 2: Re-run B_10 locally**

```
pixi run validate-rice --local B_10
```

**Step 3: Verify the supporting_junction count dropped**

```
grep -c supporting_junction validation/real_rice/results/B_10/results/ALL.mping.all_nonref_insert.txt
```

Expected: much lower than the pre-Task-2 count (which was 263 lines in B_10).

```
awk -F'\t' '{print $2}' validation/real_rice/results/B_10/results/ALL.mping.all_nonref_insert.txt | sort | uniq -c | sort -rn
```

Expected: `TTA` + `TAA` dominate; `supporting_junction` is residual (sites with `T:1` singletons or zero-junction clusters).

### Task 4: Full 10-sample re-run

**Files:** none modified.

**Step 1: Delete all 10 stale characterize sentinels**

```bash
find validation/real_rice/results -name '*.all_nonref_insert.characTErized.txt' -delete
```

**Step 2: Submit**

```
pixi run validate-rice
```

(Lets the harness submit the SLURM array with `--wait`.)

**Step 3: Read the new matrix**

```
cat validation/real_rice/report/characterized/tsd_confusion.tsv
cat validation/real_rice/report/characterized/summary.txt
```

Expected:
- `TTA → supporting_junction` cell (was 524) drops toward 0; those calls migrate to `TTA/TTA`.
- `TAA → supporting_junction` cell (was 434) drops toward 0; those calls migrate to `TAA/TAA`.
- Strict TSD agreement crosses 95%.
- `supporting_junction/supporting_junction` cell (was 42) may grow slightly as R2's depth-synthesis failures stay aligned with R3's true zero-junction sites. That's fine.

**Step 4: Update `plans/PLAN.md`**

Append a one-liner under the Phase 3 "Depth-based TSD inference" bullet noting this plan landed, with the new strict-match rate.

### Task 5: Final sanity — does anything downstream care about the new emission?

**Files:**
- Check: `src/RelocaTE3/characterize.py:152` (the parser of the `all_nonref_insert.txt`)
- Check: `validation/real_rice/normalize_relocate3_char.py`

**Step 1: Confirm `characterize` accepts any TSD string in column 2**

`grep -nE "supporting_junction|tsd ==" src/RelocaTE3/characterize.py`

Expected: a couple of branches that check `tsd == "supporting_junction"` to skip TSD-aware genotype calls; these still fire for the residual cluster of true zero-junction sites. No code change required.

**Step 2: Commit any incidental clarifying comments**

If `characterize.py` has hard-coded assumptions about which sentinels exist, document them next to the assertion. Do not change behavior in this plan.

---

## If B1 is insufficient (residual ≥ 5%)

Open a new plan implementing **B2 — port `TSD_from_read_depth`**. Sketch:

1. **New cluster lifecycle in `_cluster_reads`**: open clusters in pass 1 (collect every read's name/start/end/seq/strand into the cluster, junction OR supporting), close clusters in pass 2 (run TSD synthesis before `_write_output`).
2. **Per-cluster depth pileup**: port `tsd_finder` (`relocaTE_insertionFinder.py:843`) and `TSD_len_calculate` (line 852). Reuse the function-path's `_estimate_tsd_length_from_depth` (`insertions.py:673+`) if its signature can be widened.
3. **Sub-cluster split**: port `TSD_from_read_depth` (line 551). The hard part: pairs of left/right read pileups that overlap by ≤ TSD_len constitute one synthesized TSD; isolated pileups are dropped.
4. **Synthesize junction registrations**: for each sub-cluster, call a new `_synth_junction_register(...)` that mimics `TSD_check_cluster` (line 1249) — captures first/last `tsd_len` bases of one representative read into `te_insertions[event][tsd_start][tsd_seq]` with `pos` and `te_orient` set from the read's name suffix and strand.
5. **Existing `_write_event_start` then emits naturally** — both sides will be populated for sites the depth pass rescued, supporting_junction will be reserved for true failures.

Expected size: ~150 LOC + new fixtures. Test by hand-building a BAM with one cluster that has only right-junction reads + bracketing supporting reads, and asserting the emitted row has `R:N L:0` upgraded to `R:N L:0' (synthetic)` or the depth-derived TSD.

---

## Conventions

- One commit per task. Pre-commit (ruff check/format) must pass.
- All `pixi run …` commands from `RelocaTE3/`.
- Do not pass `--force` to `validate-rice`; deleting the characterize sentinels is enough.
- If the report at Task 4 shows strict-match still < 95%, **stop** and characterize the residual cells before iterating — there may be a small population of sites with non-3-bp captured TSDs (e.g. soft-clipped reads producing 4-mers) that need a separate fix.

## Repository pointers for a fresh session

- Branch: `fix-tsd-bug` (multiple commits leading up to this plan).
- Most-recent matrix: `RelocaTE3/validation/real_rice/report/characterized/tsd_confusion.tsv`.
- Most-recent summary: `RelocaTE3/validation/real_rice/report/characterized/summary.txt`.
- Class path code: `RelocaTE3/src/RelocaTE3/insertions.py:27-447` (`InsertionFinder`).
- Function path code (already has depth inference): `RelocaTE3/src/RelocaTE3/insertions.py:448+`.
- R2 reference: `references/RelocaTE2/scripts/relocaTE_insertionFinder.py` (line numbers above are R2's).
- CONTEXT.md warns: "The two CLI files (`__main__.py` vs `cli.py`) have diverged — touch both when changing subcommands". This plan only touches the class path used by `__main__.py`; the function path is already correct.
