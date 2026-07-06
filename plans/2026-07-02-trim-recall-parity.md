# Trim-Step Recall Parity with RelocaTE2

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.
> **For a future session:** Read the "Diagnosis" section in full before touching code. The prior plan (`2026-06-26-genotype-status-parity.md`) targeted the genome-align step and moved zero mismatches after two SLURM runs because the recall gap is not in align-genome. It is in the trim step. This plan documents the site-level trace that led to that conclusion.

**Goal:** Close the residual 318-call genotype-status-confusion gap (206 `hom/excision → het`, 112 `het → somatic`) by recovering the junction reads RelocaTE3's trim step drops. Preserve R2's classifier logic; do not tune `_classify`.

---

## Diagnosis

### Symptom (from validation run 2026-07-02, `-k 11 -w 5` on `map_genome_minimap`)

- Recall diff: **5535 sites where R2 `T:N` ≥ 2× R3 `T:N`** across the 10 rice samples (unchanged from pre-tuning baseline).
- Status confusion: **206** `hom/excision → het`, **112** `het → somatic` (unchanged).
- Full site-level trace for one representative disagreement site (B_10 Chr1:25680936):
  - R2 finds **17** junction reads at this site.
  - R3 finds **6**.
  - Of the 11 missing reads:
    - **7** are present in R3's `te_containing/*.ContainingReads.fq` (trim SAW them and classified them as junction reads) but **do not survive to `flanking/*.flankingReads.fq`**.
    - **4** are absent from `te_containing` too (the R3 TE-library aligner didn't find a usable match at all).

Command that produced the trace:
```bash
# 1. Extract R2's junction reads at the site
samtools view validation_data/real_rice/relocate2_results/B_10/repeat/bwa_aln/MSU_r7.repeat.bwa.sorted.bam \
  Chr1:25680930-25680945 | awk '$1 ~ /:(start|end):[35]$/ {print $1}'

# 2. Check whether each missing read is in R3's flanking / te_containing
grep -c '^@lh00134:.*:2165:38651:9491:' \
  validation/real_rice/results/B_10/flanking/*.fq \
  validation/real_rice/results/B_10/te_containing/*.fq
```

### Root cause: R3 classifies reads as junctions but produces an empty trimmed flank

For read `lh00134:...2165:38651:9491`:
- **R2 emits** it as `:end:3` at Chr1:25680824 with cigar `115M` — the trimmed flanking sequence is 115 bp.
- **R3's `left.ContainingReads.fq`** has `...9491:end:3` (trim's 3' branch was entered).
- **R3's `left.flankingReads.fq`** does NOT have it — the trimmed sequence failed the `len(trimmed_seq) >= len_cutoff_l` filter.

Trace through `src/RelocaTE3/librelocate.py:_trim_record` for this read:
- Best TE-library alignment (from `_parse_te_bam`, tiebreak on boundary + match):
  - `flag=16` (reverse strand), `POS=395`, `CIGAR=36M115S`, `seqlen=151`, TE length ≈ 430.
  - pysam `query_alignment_start = 0`, `query_alignment_end = 36`.
  - `qstart = 0`, `qend = 35`, `qlen = 151`, `strand = "-"`.
  - `tstart = 394`, `tend = 429` (`tend >= tlen - 3` ⇒ 3' junction).
- The 3' branch (`src/RelocaTE3/librelocate.py:281-296`):
  ```python
  if strand == "-":
      te_subseq = reverse_complement(te_subseq)
      trimmed_seq, trimmed_qual = seq[0:start], qual[0:start]   # <-- seq[0:0] = ""
      header = f"{rl_name}:end:3"
  ```
- Because `qstart = 0`, `trimmed_seq = seq[0:0] = ""`. The subsequent `if len(trimmed_seq) >= len_cutoff_l` gate at line 290 rejects it.

### Why the coordinates come out this way

`_parse_te_bam` builds `seq` by `_original_orientation` — reverse-complementing `record.query_sequence` for `flag=16` records to recover the FASTQ-frame sequence. But `qstart`/`qend` come from `record.query_alignment_start` / `_end`, which are indices into the **stored** BAM sequence, not the original. For reverse-strand alignments, these two frames are flipped:

```
FASTQ (original)   frame:  [ flank_115bp ][ te_36bp ]
BAM  (stored, RC'd) frame: [ te_36bp     ][ flank_115bp_RC ]
                            ^ qstart=0    qend=35
```

R3 then indexes the FASTQ-frame `seq` with stored-frame `qstart` = 0, producing `seq[0:0] = ""`. The **correct** trimmed_seq for R2's `strand="-"` 3' rule is `seq[0:115]` (the first 115 bp of the FASTQ-frame sequence).

The fix: for reverse-strand alignments, translate the stored-frame `qstart`/`qend` into FASTQ-frame indices before slicing:

```python
qstart_fastq = qlen - qend - 1     # first FASTQ-frame position of the alignment
qend_fastq   = qlen - qstart - 1   # last FASTQ-frame position of the alignment
```

For this read: `qstart_fastq = 151 - 35 - 1 = 115`, `qend_fastq = 151 - 0 - 1 = 150`. Then `trimmed_seq = seq[0:qstart_fastq] = seq[0:115]` — the 115-bp flank R2 sees. Matches R2's 115M cigar.

### Why R2's trim doesn't hit this bug

RelocaTE2 uses BLAT (via `--aligner blat` in the harness — see `validation_data/real_rice/example_relocate2_pipeline/01_relocate_native_cram.sh`) for the TE-library step. BLAT's PSL output reports query coordinates already in the original strand's frame regardless of alignment direction, so `parse_align_blat` never needs the flip. RelocaTE2's `parse_align_bwa` code path DOES have the same latent bug — it's just not exercised in the rice validation.

### Not confirmed but plausible

The 4 reads missing from `te_containing` (aligner-side loss) are a separate bucket. Fixing the coordinate flip likely closes ≥ 60% of the recall gap; the remaining ~40% is genuine aligner-sensitivity differences (BLAT vs minimap2) that would need a different plan (Option D from `2026-06-26-genotype-status-parity.md`).

---

## Strategy

**Option A — Fix the FASTQ-frame coordinate translation in `_parse_te_bam` / `_trim_record`.** Small, surgical change; preserves aligner choice. Should recover the ~7 of 11 missing reads per site that are already in `te_containing`.

**Option B — Add BLAT as an optional TE-library aligner backend.** Faithful R2 parity but a new dependency and larger surface area. Defer to a future plan unless Option A leaves > 5 % status-agreement gap.

**Plan executes Option A.** If matrix improvement plateaus above the target, open a follow-up plan for Option B.

---

## Task 1: Failing test that reproduces the coord-frame bug on a synthetic BAM

**Files:**
- Test: `tests/trim_reverse_strand_test.py` (new)

**Step 1: Write a synthetic BAM with one reverse-strand read whose alignment covers the last 36 bp of a 430-bp mock TE, and 115 bp of flank at the read start (FASTQ frame)**

```python
"""Reverse-strand 3' junction: coord-frame translation regression guard.

Build a 1-record BAM whose alignment cigar/flag mirror the real-world case
that reveals the trim bug (see plans/2026-07-02-trim-recall-parity.md):

  - 151-bp read
  - reverse-strand alignment (flag=16) at TE position 395 (0-based 394)
  - cigar 36M115S; alignment covers TE positions 394..429 (TE length 430)
  - the 115-bp flank lives at the FASTQ-frame START of the read

Expected: RelocaTE.trim_TE_reads emits one flanking record whose sequence is
the 115-bp flank, tagged ``:end:3``. Before the fix, the trimmed sequence is
empty and the record is silently dropped.
"""
from pathlib import Path

import pysam

from RelocaTE3.librelocate import RelocaTE


TE_LEN = 430
READ_LEN = 151
MATCH_LEN = 36
FLANK_LEN = READ_LEN - MATCH_LEN  # 115


def _write_bam(tmp_path: Path) -> Path:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"LN": TE_LEN, "SN": "mping"}]}
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        r = pysam.AlignedSegment(bam.header)
        r.query_name = "read1"
        r.flag = 16  # reverse strand
        r.reference_id = 0
        r.reference_start = 394  # 0-based; 1-based POS=395
        r.mapping_quality = 60
        r.cigartuples = [(0, MATCH_LEN), (4, FLANK_LEN)]  # 36M115S
        # BAM stores the reverse-complemented read: TE portion first, flank RC after
        r.query_sequence = "A" * MATCH_LEN + "C" * FLANK_LEN
        r.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
        r.set_tag("NM", 0)
        bam.write(r)
    sorted_bam = tmp_path / "syn.bam"
    pysam.sort("-o", str(sorted_bam), str(raw))
    pysam.index(str(sorted_bam))
    return sorted_bam


def test_reverse_strand_three_prime_junction_emits_flank(tmp_path):
    bam = _write_bam(tmp_path)
    outdir = tmp_path / "out"
    rt = RelocaTE()
    rt.write_trimmed_reads(
        name="syn",
        direction_bams=[("left", bam)],
        outdir=outdir,
        minimum_match_length=10,
        minimum_trimmed_length=10,
        mismatch_allowance=2,
    )
    flank_fq = outdir / "flanking" / "syn.left.flankingReads.fq"
    lines = [l for l in flank_fq.read_text().splitlines() if l]
    assert len(lines) == 4, lines           # 1 record: header, seq, +, qual
    assert lines[0].endswith(":end:3"), lines
    # FASTQ-frame flank: RC of the last FLANK_LEN bp of the stored seq = GGGGG...
    assert len(lines[1]) == FLANK_LEN, len(lines[1])
    assert lines[1] == reverse_complement("C" * FLANK_LEN)


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
```

**Step 2: Run, confirm failure**

Run: `pixi run pytest tests/trim_reverse_strand_test.py -v`
Expected: FAIL — the flanking FASTQ is empty (only a header if any), or the test asserts a non-empty file that has zero data lines.

**Step 3: Commit the failing test**

```bash
cd /bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/RelocaTE3
git add tests/trim_reverse_strand_test.py plans/2026-07-02-trim-recall-parity.md
git commit -m "test(trim): failing regression for reverse-strand coord-frame bug"
```

---

## Task 2: Translate stored-frame coords to FASTQ-frame in `_parse_te_bam` (Option A)

**Files:**
- Modify: `src/RelocaTE3/librelocate.py:352-432` (`_parse_te_bam`)
- Modify: `src/RelocaTE3/trim.py:56-124` (`parse_te_alignments`) — same fix, mirrored for the standalone CLI

**Step 1: Add the coord-frame flip immediately after the pysam extractions**

In `_parse_te_bam`, after computing `qstart`, `qend`, `qlen`, `strand`, `record.is_reverse`, insert:

```python
                # pysam reports qstart/qend in the STORED sequence frame. For
                # reverse-strand alignments the stored seq is RC of the FASTQ
                # sequence, so we flip the coords to the FASTQ frame here —
                # every downstream slice (_trim_record etc.) uses the FASTQ
                # sequence returned by _original_orientation. See
                # plans/2026-07-02-trim-recall-parity.md for the bug that
                # this closes.
                if record.is_reverse and qlen > 0:
                    qstart, qend = qlen - qend - 1, qlen - qstart - 1
```

Then the record still stores `start=qstart, end=qend` in FASTQ frame.

Apply the same 3 lines in `trim.py:parse_te_alignments` between line 100 and 102 (before `_boundary_score` uses `qstart`/`qend`).

**Step 2: Verify boundary_score still identifies junction ends correctly**

The boundary score checks `qstart <= 2` (touches read left end) and `qend >= qlen - 3` (touches read right end). After the flip, these still describe FASTQ-frame ends, which is what the classifier logic needs. Reverse-strand alignments where the TE was at the read's stored-start (FASTQ-end) now correctly boundary-score on `qend >= qlen - 3`.

**Step 3: Run the failing test — should now PASS**

```
pixi run pytest tests/trim_reverse_strand_test.py -v
```

Expected: PASS.

**Step 4: Full regression sweep**

```
pixi run test
```

Expected: 45 green (44 pre-existing + 1 new). If the acceptance test regresses, STOP — the coord flip may be inadvertently changing junction-tag suffixes across the whole benchmark. Diagnose which reads changed tags and whether the change is correct (matches R2) or a bug.

**Step 5: Commit**

```bash
cd RelocaTE3
git add src/RelocaTE3/librelocate.py src/RelocaTE3/trim.py
git commit -m "fix(trim): translate reverse-strand coords to FASTQ frame

pysam's query_alignment_start/_end are indexes into the stored BAM
sequence, which is the reverse complement of the FASTQ sequence for
flag=16 records. RelocaTE3's _parse_te_bam / parse_te_alignments then
sliced the reconstructed FASTQ-frame sequence with those stored-frame
coordinates, producing empty trimmed flanks for reverse-strand 3'-end
junctions (and mislabeled boundaries elsewhere).

Fix: flip qstart/qend to the FASTQ frame immediately after extraction,
so every downstream consumer indexes consistently. Adds a synthetic
regression test that reproduces the exact real-world case (151-bp read,
cigar 36M115S, TE 3' end match) that traced this bug.

Expected impact on the rice validation:
  - Recovers ~7 of 11 missing junction reads per representative site.
  - Closes >= 60 percent of the 318 status-confusion mismatches
    (206 hom/excision -> het, 112 het -> somatic).

See plans/2026-07-02-trim-recall-parity.md."
```

---

## Task 3: Rerun the SLURM validation and measure

**Files:** none modified.

**Step 1: Ask the user (do NOT delete on your own) to clear the stale outputs**

Prompt: *"The fix invalidates every downstream output from trim on. Delete `flanking/`, `te_containing/`, `te_portions/`, the `*.left.bam` / `*.right.bam` TE-library BAMs, the `*.repeat.minimap.sorted.bam` genome BAM, and the per-sample `results/*.all_nonref_insert*.txt` files? (Yes / No)"*

If yes, run:
```bash
for d in flanking te_containing te_portions; do
  find validation/real_rice/results -type d -name "$d" -exec rm -rf {} + 2>/dev/null
done
find validation/real_rice/results -maxdepth 2 -name '*.bam*' -delete
find validation/real_rice/results -name '*.all_nonref_insert*.txt' -delete
```

**Step 2: User submits `pixi run validate-rice`**

Wait for their confirmation that the SLURM array finished.

**Step 3: Rerun the recall-diff**

```bash
cd validation/real_rice
python3 junction_recall_diff.py --out post_taskB_recall_diff.tsv
```

Compare `tail -n +2 post_taskB_recall_diff.tsv | wc -l` to the baseline 5535. Target: ≤ 2500 (≥ 55 % reduction).

**Step 4: Read the new status matrix**

```bash
cat validation/real_rice/report/characterized/status_confusion.tsv
cat validation/real_rice/report/characterized/summary.txt
```

Expected cells:
| Cell | Baseline | Target after Task 2 |
|---|---|---|
| `hom/excision → het` | 206 | ≤ 80 |
| `het → somatic` | 112 | ≤ 50 |
| Status agreement | 93.67 % | ≥ 96 % |

Also expected: the TSD confusion matrix should stay essentially unchanged (the trim fix affects which reads reach genome-align but not the TSD-capture rule inside insertions.py).

**Step 5: If Task 2 closes ≥ 60 % of the 318 mismatches, write a brief commit note and stop**

Update `plans/PLAN.md` Phase 3 bullet with the new numbers and cross-reference this plan.

**Step 6: If Task 2 closes < 60 %, open Task 4 (Option B: BLAT backend) as a new plan**

See "Follow-up options" below.

---

## Follow-up options (do not implement in this plan)

**Option B — Add BLAT as TE-library aligner.** R2's actual configuration. Adds `blat` to `pixi.toml`, a `--te-aligner=blat` switch, a `parse_blat_output` parser. ~300 LOC + tests. Only worth doing if Option A leaves > 5 % status-agreement gap.

**Option C — Secondary-alignment fallback in `_parse_te_bam`.** When the current "best by boundary" alignment produces an empty trimmed flank, fall back to the second-best. Small change, may recover a handful more reads but does not close the aligner-sensitivity floor.

---

## Repository pointers for a fresh session

- Branch: `fix-genotyping-mismatch` (contains all TSD work + the misdirected Task 2 minimap2 tuning from the prior plan + this plan).
- Prior plan (align-genome dead end): `plans/2026-06-26-genotype-status-parity.md`.
- Diagnostic tool: `validation/real_rice/junction_recall_diff.py`.
- R2 reference (all Perl/Python 2): `references/RelocaTE2/scripts/relocaTE_trim.py`.
- R3 code to fix: `src/RelocaTE3/librelocate.py:352-432` and `src/RelocaTE3/trim.py:56-124`.
- R2's actual TE-library aligner in this validation: BLAT (per `validation_data/real_rice/example_relocate2_pipeline/01_relocate_native_cram.sh:25 aligner=blat`).
- CONTEXT.md's warning about diverged CLI files: still holds — `__main__.py` and `cli.py` have different `find-insertions` / `align-genome` paths, and both use the same `librelocate.py:_parse_te_bam` for trim, so the fix here lands on the harness path automatically.

## Conventions

- One commit per task. Pre-commit (ruff check/format) must pass.
- All `pixi run …` commands from `RelocaTE3/`.
- Do **not** touch `_classify` in `characterize.py`.
- If Task 3 reveals a new confusion cell (e.g. `het → hom` appears), STOP — the coord flip may be misclassifying reads elsewhere. The synthetic test in Task 1 covers one specific case; a real-world regression means a second failing test is needed before iterating.
