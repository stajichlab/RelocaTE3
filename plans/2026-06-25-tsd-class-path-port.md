# TSD Parity for the Legacy `InsertionFinder` Class Path

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Close the TSD-string gap between RelocaTE3 and legacy RelocaTE2 on the rice validation harness, which runs through the legacy `InsertionFinder` class (not the function-based `find_insertions` already fixed in `2026-06-24-tsd-depth-inference.md`).

## Why this exists

The prior plan ported depth-based TSD inference (`_estimate_tsd_length_from_depth`, `_capture_tsd_from_read`, the rewritten `_make_insertion`) into the **function-based** path in `insertions.py:888`. That path is wired up by `cli.py:218` and `pipeline.py:124`. Acceptance tests pass.

But `pyproject.toml` registers the `relocaTE3` entry point as `__main__:main`, and `__main__.py:596` invokes the **legacy class** `InsertionFinder` for the `find-insertions` subcommand. The validation harness (`validation/real_rice/run_relocate3.sh:237`) calls `relocaTE3 find-insertions`, so all of the previous plan's depth work is bypassed at runtime.

Empirical confirmation (post-run, 2026-06-25):
- Strict TSD agreement moved 41.9% → 46.8% (only Task 6's tiny tweak fired; reverted in this branch).
- Canonical (rev-comp-tolerant) TSD agreement 46.8% — confirms there is no strand-canonicalization bug to fix.
- Total matched calls dropped 4990 → 4444 (regression caused by Task 6's misclassification of `UNK` calls as `supporting_junction`).

R2 default behavior: `relocaTE2.py:346` writes `regex.txt` containing the sentinel `UNK`, which triggers `TSD_from_read_depth` → `TSD_check_cluster` in `relocaTE_insertionFinder.py`. There the TSD regex becomes `'.'*TSD_len`, so the captured `TSD_seq` is whatever the read's first/last `TSD_len` bases are (TTA, TAA, or anything else). The TSD *length* is estimated per-site by the depth pileup.

## Architecture decision: two real options

**Option D — surgical "wildcard TSD" config flip (recommended for mPing-only validation).**
The class path's `_tsd_check` already uses `r5_tsd = re.compile(rf"^({tsd})")` and assigns `tsd_seq = m.group(1)` — i.e., it captures whatever the regex matched. If we pass `tsd = "..."` (three-char wildcard) instead of `tsd = "TTA"`, every junction read literally captures its first/last 3 bases. This mirrors R2's depth-mode behavior for the special case of a known, fixed-length TSD (mPing = 3 bp). Implementation cost: ~1 line in the harness `config.toml`, plus a hardening tweak in `__main__.py`'s entry-point check so wildcard patterns aren't reflexively rejected.

Pros: tiny change, no two-pass refactor, the class path stays single-pass. Covers ~99% of the matrix delta on this dataset.
Cons: hardcodes TSD length to 3 for any TE family the harness runs on. Long/odd R2 TSDs (`TTAGGACGT`, `TTAGCAGGTCGTTCAAAGAGCGT`, …) won't match — but those are <20 cases and R2 itself only produces them by overestimating length.

**Option C — full depth-inference port into the class (general but invasive).**
Add per-cluster depth-based length estimation to `InsertionFinder`, mirroring R2's `tsd_finder`/`TSD_len_calculate`. Requires:
- Lifting the `UNK`-rejection at `insertions.py:71-74`.
- Teaching `_cluster_reads` to also collect non-junction reads (currently discarded) and to operate in two passes per cluster: pass 1 collects all reads, pass 2 estimates TSD length and runs `_tsd_check` with a dynamic `.{N}` pattern.
- New cluster lifecycle hook (cluster close/flush) because depth estimation needs the cluster *complete* before length can be computed.

Pros: faithful R2 parity for all TE families and TSD lengths.
Cons: ~200 LOC refactor, harder to keep alongside the function path, duplicates work already done in the function path.

**Plan executes Option D.** If/when the validation needs TE families with variable TSD length we'll revisit Option C.

---

## Task 1: Reproduce the regression baseline with a tiny dedicated test

**Files:**
- Test: `tests/insertions_tsd_class_parity_test.py` (new)

**Step 1: Write a failing test that exercises `InsertionFinder.find_insertions` with a synthetic BAM where every junction read has a 3-bp TAA prefix/suffix**

```python
"""Class-path TSD parity test: wildcard TSD mode emits read-captured 3-mers."""
from pathlib import Path

import pysam

from RelocaTE3.insertions import InsertionFinder


def _write_bam(tmp_path: Path) -> Path:
    """Build a 1-chrom BAM with two junction reads framing one TAA TSD."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"LN": 1000, "SN": "Chr1"}],
    }
    bam_path = tmp_path / "syn.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        # right-junction read, flank seq starts with TAA at position 100
        r1 = pysam.AlignedSegment(bam.header)
        r1.query_name = "read1:start:5"
        r1.flag = 0  # forward strand
        r1.reference_id = 0
        r1.reference_start = 99  # 0-based; reports as 100 (1-based)
        r1.mapping_quality = 60
        r1.cigartuples = [(0, 10)]
        r1.query_sequence = "TAAGGGCCAA"
        r1.query_qualities = pysam.qualitystring_to_array("I" * 10)
        r1.set_tag("NM", 0)
        bam.write(r1)
        # left-junction read, flank seq ends with TAA at position 102 (1-based)
        r2 = pysam.AlignedSegment(bam.header)
        r2.query_name = "read2:end:3"
        r2.flag = 0
        r2.reference_id = 0
        r2.reference_start = 92  # spans 93..102 1-based
        r2.mapping_quality = 60
        r2.cigartuples = [(0, 10)]
        r2.query_sequence = "AAAAAA TAA".replace(" ", "")
        r2.query_qualities = pysam.qualitystring_to_array("I" * 10)
        r2.set_tag("NM", 0)
        bam.write(r2)
    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path


def _write_read_repeat(tmp_path: Path) -> Path:
    p = tmp_path / "read_repeat.txt"
    p.write_text("read1\tmPing\t+\nread2\tmPing\t-\n")
    return p


def test_class_path_wildcard_tsd_captures_read_bases(tmp_path):
    bam = _write_bam(tmp_path)
    read_repeat = _write_read_repeat(tmp_path)
    outdir = tmp_path / "out"
    finder = InsertionFinder(mismatch_allow=2, min_mapq=1)
    out_txt = finder.find_insertions(
        bam_file=bam,
        read_repeat_file=read_repeat,
        tsd="...",  # wildcard: capture any 3 bases from each junction read
        target="Chr1",
        sample="syn",
        outdir=outdir,
        te_name="mPing",
    )
    rows = [
        line.split("\t")
        for line in Path(out_txt).read_text().splitlines()
        if line and not line.startswith("strain")
    ]
    assert len(rows) == 1, rows
    # TSD column (index 1) must be the read-captured "TAA", not "UNK" or "..."
    assert rows[0][1] == "TAA", rows
```

**Step 2: Run, confirm failure**

Run: `pixi run pytest tests/insertions_tsd_class_parity_test.py -v`
Expected: FAIL at the entry-point rejection — `tsd="..."` may be accepted (no `UNK` substring), but if the test fails for a different reason (e.g. cluster never closes, or TSD pattern is treated as literal), capture the failure message; that drives Task 2.

**Step 3: Commit**

```bash
cd RelocaTE3
git add tests/insertions_tsd_class_parity_test.py plans/2026-06-25-tsd-class-path-port.md
git commit -m "test(insertions): failing class-path TSD parity test (wildcard TSD)"
```

---

## Task 2: Make the class-path entry point accept regex TSD patterns

**Files:**
- Modify: `src/RelocaTE3/insertions.py:71-74` (the UNK rejection)

**Step 1: Inspect the failure from Task 1, then loosen the entry point**

Current guard:
```python
if re.search(r"UNK|UKN|unknown", tsd, re.IGNORECASE):
    raise NotImplementedError(
        "TSD-unknown (read-depth) inference is not yet ported; provide a TSD motif."
    )
```

Change to:
```python
if re.search(r"UNK|UKN|unknown", tsd, re.IGNORECASE):
    # The literal "UNK" sentinel triggers R2's full depth-mode pipeline, which
    # is not yet ported. Regex-friendly TSD patterns (e.g. "..." for a 3 bp
    # wildcard) are accepted and flow through TSD_check naturally.
    raise NotImplementedError(
        "TSD-unknown (read-depth) inference is not yet ported; provide a TSD "
        'motif or a fixed-length wildcard regex (e.g. "..." for a 3 bp TSD).'
    )
```

(The rejection itself is unchanged; only the error message is updated so a future maintainer knows wildcards are legal.)

**Step 2: Run the test from Task 1**

Run: `pixi run pytest tests/insertions_tsd_class_parity_test.py -v`
Expected: now PASS, because `"..."` was already accepted — the rejection only fires on the literal `UNK`/`UKN`/`unknown` strings. If the test still fails, the failure is in `_tsd_check`'s `tsd_start = end - len(tsd)` arithmetic (treating regex chars as length characters) — which is correct for `"..."` because `len("...") == 3`.

**Step 3: Commit**

```bash
cd RelocaTE3
git add src/RelocaTE3/insertions.py
git commit -m "docs(insertions): clarify that wildcard TSD regex (e.g. '...') is accepted"
```

---

## Task 3: Make the harness pass a wildcard TSD

**Files:**
- Modify: `validation/real_rice/config.toml` (the `[relocate3]` section)

**Step 1: Read the current value**

Run: `grep -nE 'tsd\s*=' validation/real_rice/config.toml`
Expected: `tsd = "TTA"  # find-insertions --tsd (mPing canonical TSD)`

**Step 2: Change to the wildcard pattern**

```toml
# find-insertions --tsd. "..." is a 3-bp wildcard so R3's class path captures
# the literal TSD bases from each junction read (matching R2's depth-mode
# behavior on the mPing dataset). Use a known motif ("TTA") only when you
# want to require a literal match.
tsd = "..."
```

**Step 3: Sanity-check that the harness picks it up**

Run: `pixi run python -c "from validation.real_rice._config import load_config; print(load_config('validation/real_rice/config.toml')['relocate3']['tsd'])"`
Expected: `...`

**Step 4: Commit**

```bash
cd RelocaTE3
git add validation/real_rice/config.toml
git commit -m "validation: pass 3-bp wildcard TSD to find-insertions for R2 parity"
```

---

## Task 4: Validate on one sample (B_10) before the full SLURM run

**Files:** none modified. This is a runtime check.

**Step 1: Delete the stale characterize sentinel for B_10**

Run:
```bash
find validation/real_rice/results/B_10 -name '*.all_nonref_insert.characTErized.txt' -delete
```

**Step 2: Re-run B_10 locally**

Run: `pixi run validate-rice --local B_10`
Expected: completes without SLURM, regenerates `results/B_10/...all_nonref_insert.txt` and `.characTErized.txt`.

**Step 3: Spot-check the nonref TXT**

Run:
```bash
awk -F'\t' '{print $2}' validation/real_rice/results/B_10/results/ALL.mping.all_nonref_insert.txt | sort | uniq -c | sort -rn
```
Expected: the dominant TSD strings should be `TTA` and `TAA` in roughly equal counts, with few or no `UNK` / `supporting_junction` entries (the few remaining `supporting_junction` are sites where no junction read mapped at all).

**Step 4: Spot-check the per-sample summary**

Open `validation/real_rice/report/characterized/summary.txt` (Stage B will have re-run as part of `validate-rice`). The strict-TSD agreement should jump from ~47% toward 95%+.

**Step 5: Commit any incidental fixes if needed**

If Task 4 reveals follow-on issues (e.g. characterize.py treating `TAA` differently from `TTA`), open a sibling task and address before the full run.

---

## Task 5: Full 10-sample re-run + report

**Files:** none modified.

**Step 1: Delete all 10 stale characterize sentinels**

Run:
```bash
find validation/real_rice/results -name '*.all_nonref_insert.characTErized.txt' -delete
```

**Step 2: Submit the SLURM array**

Run: `pixi run validate-rice`
(The harness submits `sbatch --wait`; let it complete before the report stage.)

**Step 3: Verify the report regenerated**

Run: `stat -c '%y  %n' validation/real_rice/report/characterized/summary.txt`
Expected: today's date.

**Step 4: Read the new confusion matrix**

Run: `cat validation/real_rice/report/characterized/tsd_confusion.tsv`
Expected outcome:
- `TTA/TTA` and `TAA/TAA` cells dominate.
- The `TAA → UNK` cell is gone (UNK no longer emitted for these sites).
- `* → supporting_junction` cells shrink dramatically.
- Strict-match rate ≥ 0.95; canonical-match rate within a hair of strict.

**Step 5: Update PLAN.md once the matrix collapses**

Append a one-liner under Phase 3 noting the class-path parity landed via this plan; cross-reference both plan files.

---

## Conventions

- One commit per task. Pre-commit (ruff check/format) must pass.
- All `pixi run …` commands from `RelocaTE3/`.
- Do not pass `--force` to `validate-rice`; deleting the characterize sentinels is enough to re-trigger steps 5+7 without re-running expensive trim/align.
- If Task 4's spot-check still shows a large `* → supporting_junction` cell, **stop** and reassess — that means clusters are missing junction reads entirely (an upstream alignment issue, not a TSD-capture issue).
