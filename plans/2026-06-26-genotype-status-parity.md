# Genotype-Status Parity with RelocaTE2

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.
> **For a future maintainer (or fresh Claude session):** Read the "Diagnosis" section in full before touching code. The classifier in `src/RelocaTE3/characterize.py:_classify` is **already byte-identical to RelocaTE2's `characterizer.pl`**. Do not tune the thresholds.

**Goal:** Drop the residual ~320 status-confusion mismatches in the rice validation matrix (`homozygous/excision_no_footprint → heterozygous`: 206; `heterozygous → somatic_insertion`: 112; misc: <5) by closing the upstream junction-read recall gap that feeds the classifier the wrong inputs. Preserves R2's classifier logic exactly.

---

## Diagnosis

### Current status (validation run 2026-06-26 on branch `fix-tsd-bug`)

`validation/real_rice/report/characterized/status_confusion.tsv`:

| R2 \ R3 | hom | het | hom/excision | somatic |
|---|---|---|---|---|
| **homozygous** | 4317 | – | – | – |
| **heterozygous** | – | 260 | 2 | **112** |
| **hom/excision** | 1 | **206** | 153 | – |
| **somatic** | – | 1 | – | 34 |

Status agreement: 4764 / 5086 = **93.67%**. Diagonals are otherwise clean; the two off-diagonal masses (206 + 112 = 318 mismatches) account for the entire gap.

### Confirmed: the classifier is not the bug

`grep -nE "_classify" src/RelocaTE3/characterize.py` → `_classify(average_flankers, spanners)` at `characterize.py:275`. Side-by-side with `characterizer.pl:149-185` it is a faithful port — same cascade, same thresholds, same precedence. **Do not edit `_classify`.**

### Confirmed: spanners match exactly

Sample of disagreement rows from `matched_calls.tsv`:

```
B_10 Chr1:25680936  R2: flank=8.5  span=3   R3: flank=3    span=3   (hom/excision -> het)
B_10 Chr1:36270511  R2: flank=10.5 span=3   R3: flank=1.5  span=3   (hom/excision -> het)
B_10 Chr1:39658243  R2: flank=9    span=2   R3: flank=2    span=2   (hom/excision -> het)
B_10 Chr5:1779800   R2: flank=9    span=2   R3: flank=2    span=2   (hom/excision -> het)
B_10 Chr8:25327897  R2: flank=5.5  span=1   R3: flank=2    span=1   (hom/excision -> het)
```

Spanner column is identical for every row. The flanker column is consistently lower for R3.

### Confirmed: R3 finds fewer junction reads upstream

`grep 25680936 results/B_10/results/ALL.mping.all_nonref_insert.txt` and the equivalent in `validation_data/real_rice/relocate2_results/B_10/repeat/results/ALL.all_nonref_insert.txt`:

```
R3: mPing  TTA  B_10  Chr1  25680936..25680938  -  T:6   R:2  L:4   ST:0  SR:0  SL:0
R2: mPing  TTA  B_10  Chr1  25680936..25680938  -  T:17  R:8  L:9   ST:0  SR:0  SL:0
```

`average_flankers = T:N / 2`, so R3 feeds 3 to the classifier and R2 feeds 8.5. Plugging both into the (identical) decision cascade:

- R2 (flank=8.5, span=3) → `avg_flankers >= 5 AND spanners < 5` → **homozygous/excision_no_footprint** ✓
- R3 (flank=3, span=3) → `avg_flankers >= 5` fails → falls through to `abs(avg - span) <= 5` → **heterozygous** ✓

Both classifications are correct *given their inputs*. The inputs differ because R3's genome BAM has fewer junction reads at this site.

### Root cause: known minimap2-vs-blat short-flank gap

PLAN.md (line 260-262, 277-283) already flags this: *"short flanks (<~15 bp) don't map uniquely under minimap2 `-x sr` (the blat-sensitivity tradeoff noted in §7)"*. RelocaTE2 uses `bwa mem -k 15 -T 10` for the genome re-align; RelocaTE3 uses `minimap2 -x sr -k 13 -w 6` (`src/RelocaTE3/align.py:279-290`). For short flanks (10-25 bp), bwa-mem's seed-and-extend places reads that minimap2 leaves unmapped or as low-quality multimappers.

The same site that genotypes correctly under R2 (with T:17) lands in the wrong classifier band under R3 (with T:6) **purely because of upstream alignment recall**. Touching `characterize.py` to compensate (e.g. counting supporting reads as flankers, or padding the flanker estimate) would diverge from R2's logic and is **explicitly out of scope** for this plan.

---

## Strategy

Three options, in increasing blast radius and faithfulness:

### Option A — Tighten minimap2 short-read sensitivity in `genome_align.py`

Reduce `-k` and `-w` further, enable secondary alignments (`--secondary=yes -N`), drop the chaining penalty (`-p 0.5`). Already partially applied at the trim/library-mapping step (PLAN.md line 221: *"the TE-read step now maps with `-k 11 -w 5 -N 20 -p 0.5`"*), but the genome re-align step is more conservative (`-k 13 -w 6`, no `-N`, no `-p`). Bringing those two steps to parity is a one-line change with no new dependencies.

**Pros:** smallest blast radius, no new tool, no schema change.
**Cons:** minimap2 may have a fundamental sensitivity floor on 10-15 bp flanks that no tuning crosses; bwa-mem outperforms it on this regime by design.
**Expected impact:** halves or two-thirds of the 318 mismatches resolve; the residual is the regime where minimap2 can't recover the placement.

### Option B — Map flanks as proper pairs in `genome_align.py`

Currently `map_reads_to_genome` (`align.py:252+`) concatenates all flanking FASTQs and maps single-end (`combined_fq`, no `-I` insert-size hint). PLAN.md line 260 calls this out: *"reads are mapped single-end rather than as proper pairs"*. With paired-end mapping, the flank's mate provides an anchor for the short flank to be placed uniquely.

**Pros:** uses mate-pair information already present in the FASTQ; aligned with how R2 operates.
**Cons:** requires preserving R1/R2 separation through trim → align; some flanking reads have no mate (orphans).
**Expected impact:** closes some of the residual after Option A; not a silver bullet (paired-end helps when one mate has a unique anchor and the flank-bearing mate has a short flank).

### Option C — Add a bwa-mem alignment path mirroring R2 exactly

Optional aligner backend (`Aligner.default_aligner = "bwa-mem"`). `bwa mem -k 15 -T 10 ref flank.fq` matches R2's invocation byte-for-byte. Either replace minimap2 for the genome re-align step entirely, or expose as a `--genome-aligner` flag.

**Pros:** faithful R2 parity; closes the recall gap fully on this benchmark.
**Cons:** new external dependency (bwa or bwa-mem2); diverges from PLAN.md's "minimap2-only" goal at §7. PLAN.md line 124-126 says: *"No bwa/blat/bowtie2/seqtk/perl — consistent with the intent to standardize on minimap2 + pysam."*
**Expected impact:** closes the gap to ≤ 2% (residual is genuine biological ambiguity).

**Plan executes Option A first; Option B if Option A doesn't close ≥ 60% of the 318; Option C only with explicit user sign-off given the dependency trade-off.**

---

## Task 1: Build a per-site recall diff to drive iteration

**Files:**
- New: `validation/real_rice/junction_recall_diff.py`

**Step 1: Write the diff script**

The script reads both R2 and R3 nonref TXTs, joins by `(sample, chrom, pos)` within ±5 bp, and prints rows where R2 has ≥ 2× the junction-read count R3 has. This becomes the iteration target: each tuning change should shrink this list.

```python
"""Per-site junction-read recall diff (R3 minimap2 vs R2 bwa-mem).

Joins R2 and R3 ``all_nonref_insert.txt`` files by (sample, chrom, ±5 bp) and
prints sites where R2 has substantially more junction reads. Output:

  sample  chrom  pos  r2_T  r3_T  ratio  r2_status  r3_status

Used by ``plans/2026-06-26-genotype-status-parity.md`` to measure whether a
genome-align tuning change closes the recall gap that drives the residual
hom/excision -> het and het -> somatic mismatches.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

_COUNT = re.compile(r"T:(\d+).*R:(\d+).*L:(\d+)")


def parse_nonref(path: Path):
    """Yield (chrom, start_pos, total) tuples from an all_nonref_insert.txt."""
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or not cols[3].startswith("Chr"):
                continue
            chrom = cols[3]
            m_coord = re.match(r"(\d+)\.\.(\d+)", cols[4])
            if not m_coord:
                continue
            pos = int(m_coord.group(1))
            m_counts = _COUNT.search("\t".join(cols[6:9]))
            if not m_counts:
                continue
            total = int(m_counts.group(1))
            yield chrom, pos, total


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r2-root", required=True, help="validation_data/real_rice/relocate2_results")
    ap.add_argument("--r3-root", required=True, help="validation/real_rice/results")
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--window", type=int, default=5)
    ap.add_argument("--min-ratio", type=float, default=2.0)
    ap.add_argument("--out", default="-")
    args = ap.parse_args(argv)

    rows = []
    for sample in args.samples:
        r2 = Path(args.r2_root) / sample / "repeat/results/ALL.all_nonref_insert.txt"
        r3 = Path(args.r3_root) / sample / "results/ALL.mping.all_nonref_insert.txt"
        if not (r2.exists() and r3.exists()):
            print(f"WARN: missing for {sample}: r2={r2.exists()} r3={r3.exists()}", file=sys.stderr)
            continue
        r2_by_chrom = defaultdict(list)
        for chrom, pos, total in parse_nonref(r2):
            r2_by_chrom[chrom].append((pos, total))
        r3_by_chrom = defaultdict(list)
        for chrom, pos, total in parse_nonref(r3):
            r3_by_chrom[chrom].append((pos, total))
        for chrom, items in r2_by_chrom.items():
            r3_index = sorted(r3_by_chrom.get(chrom, []))
            for r2_pos, r2_total in items:
                r3_match = None
                for r3_pos, r3_total in r3_index:
                    if abs(r3_pos - r2_pos) <= args.window:
                        if r3_match is None or abs(r3_pos - r2_pos) < abs(r3_match[0] - r2_pos):
                            r3_match = (r3_pos, r3_total)
                if r3_match is None:
                    continue
                r3_total = r3_match[1]
                if r3_total == 0:
                    ratio = float("inf")
                else:
                    ratio = r2_total / r3_total
                if ratio >= args.min_ratio:
                    rows.append({
                        "sample": sample, "chrom": chrom, "pos": r2_pos,
                        "r2_T": r2_total, "r3_T": r3_total, "ratio": f"{ratio:.2f}",
                    })

    out = sys.stdout if args.out == "-" else open(args.out, "w", newline="")
    w = csv.DictWriter(out, fieldnames=["sample", "chrom", "pos", "r2_T", "r3_T", "ratio"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)
    if args.out != "-":
        out.close()
    print(f"sites with R2/R3 T-ratio >= {args.min_ratio}: {len(rows)}", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
```

**Step 2: Run it for the current (baseline) state**

```bash
cd validation/real_rice
SAMPLES=$(cut -d, -f1 ../../validation_data/real_rice/sample_file/samples.csv | tail -n +2 | tr '\n' ' ')
python3 junction_recall_diff.py \
  --r2-root ../../validation_data/real_rice/relocate2_results \
  --r3-root results \
  --samples $SAMPLES \
  --out baseline_recall_diff.tsv
```

Capture stderr line: `sites with R2/R3 T-ratio >= 2.0: N`. Record N as the baseline. **Save this number in the next commit message** — every later iteration's improvement is measured against it.

**Step 3: Commit the diff script**

```bash
cd /bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/RelocaTE3
git add validation/real_rice/junction_recall_diff.py plans/2026-06-26-genotype-status-parity.md
git commit -m "validation: junction-read recall diff tool for genotype-status parity work

Per-site comparison of R2 vs R3 T:N junction-read counts. Iteration target
for closing the status-confusion mismatches (hom/excision -> het, 206;
het -> somatic, 112) that originate in upstream alignment recall, not in
the classifier. See plans/2026-06-26-genotype-status-parity.md."
```

---

## Task 2: Option A — tune minimap2 in `map_reads_to_genome`

**Files:**
- Modify: `src/RelocaTE3/align.py:279-290` (the `cmd` list in `map_reads_to_genome`)

**Step 1: Read the current invocation**

`sed -n '279,295p' src/RelocaTE3/align.py` — note `-k 13 -w 6`, no `--secondary`, no `-N`, no `-p`. Compare to the TE-mapping step at `align.py:220-240` (which already uses `-k 11 -w 5 -N 20 -p 0.5` per PLAN.md line 221).

**Step 2: Mirror the TE-step parameters on the genome step**

Change the `cmd` list to:

```python
cmd = [
    self.minimap,
    "-t",
    str(cpu_threads),
    "-a",
    "-x",
    "sr",
    "-k",
    "11",         # was 13; smaller seed -> shorter flanks find anchors
    "-w",
    "5",          # was 6; finer minimizer window
    "--secondary=yes",
    "-N",
    "20",         # report up to 20 secondaries so a short flank that has
                  # multiple equally-good placements doesn't get dropped
    "-p",
    "0.5",        # lower chaining penalty -> retain weaker chains
    "-o",
    temp_sam,
    str(genome),
    combined_fq,
]
```

**Step 3: Run the acceptance test to confirm Phase-3 recall isn't lost**

```bash
pixi run pytest tests/acceptance_test.py -v
```

Expected: still ≥ 178/200 (the recall target from PLAN.md line 218-220).

**Step 4: Commit the tuning**

```bash
git add src/RelocaTE3/align.py
git commit -m "feat(align): tighten minimap2 short-read sensitivity in map_reads_to_genome

Match the TE-mapping step's parameters (-k 11 -w 5 -N 20 -p 0.5) on the
genome re-align step. Targets the residual hom/excision -> het and
het -> somatic genotype-status mismatches that originate in low junction-
read recall, not in the classifier. See
plans/2026-06-26-genotype-status-parity.md task 2."
```

---

## Task 3: Measure Task 2's impact on the validation matrix

**Files:** none modified.

**Step 1: Clear the find-insertions + characterize sentinels for all 10 samples**

Find-insertions has no skip-if-exists guard, but characterize does:

```bash
find validation/real_rice/results -name '*.all_nonref_insert.characTErized.txt' -delete
```

The genome BAM is also stale (it's the input to find-insertions). We need to re-run from align-genome forward. Inspect `run_relocate3.sh:230` — that step skips if `$GENOME_BAM` exists. To force a re-align:

```bash
find validation/real_rice/results -name '*.flank.genome.sorted.bam*' -delete
find validation/real_rice/results -name '*.fullreads.genome.sorted.bam*' -delete
```

(If the actual sentinel file names differ, grep `run_relocate3.sh:225-232` for the exact `$GENOME_BAM` path and clear that. Do **not** delete the trim output — that step is the expensive one and is unaffected by minimap2 genome-align changes.)

**Step 2: Submit the SLURM array**

```bash
pixi run validate-rice
```

**Step 3: Re-run the recall diff**

```bash
cd validation/real_rice
python3 junction_recall_diff.py \
  --r2-root ../../validation_data/real_rice/relocate2_results \
  --r3-root results \
  --samples $SAMPLES \
  --out post_taskA_recall_diff.tsv
```

Compare `wc -l baseline_recall_diff.tsv` and `wc -l post_taskA_recall_diff.tsv`. Expected: post-Task-A is at least 30% smaller.

**Step 4: Read the new status-confusion matrix**

```bash
cat validation/real_rice/report/characterized/status_confusion.tsv
cat validation/real_rice/report/characterized/summary.txt
```

Expected change in matrix cells:

| Cell | Baseline | Target after Task 2 |
|---|---|---|
| `hom/excision → het` | 206 | ≤ 100 |
| `het → somatic` | 112 | ≤ 60 |
| Status agreement | 93.67% | ≥ 96% |

**Step 5: If Task 2 closes ≥ 60% of the gap, write a brief commit note and stop**

```bash
git commit --allow-empty -m "validation: Task 2 closed N/318 status mismatches (X% -> Y% status agreement)"
```

Then update `plans/PLAN.md` Phase 3 bullet: *"Genotype-status parity (DONE / partial 2026-MM-DD)"* with the new agreement number and a cross-reference to this plan.

**Step 6: If Task 2 closes < 60% of the gap, open Task 4 (Option B: paired-end mapping)**

See the sketched Task 4 below.

---

## Task 4 (conditional): Option B — paired-end mapping for flanks

**Trigger:** Task 3 measures < 60% closure of the 318 baseline mismatches.

**Files:**
- Modify: `src/RelocaTE3/align.py:map_reads_to_genome`
- Modify: `src/RelocaTE3/genome_align.py` (the caller)
- Modify: `src/RelocaTE3/trim.py` (preserve R1/R2 separation through the flank FASTQ emission)

**Architecture sketch**

Currently:
1. `trim.py` emits a single `flank.fq.gz` per sample (R1+R2 mixed).
2. `genome_align.py` calls `Aligner.map_reads_to_genome(flank.fq.gz)`.
3. `align.py:map_reads_to_genome` concatenates inputs and runs minimap2 single-end.

Proposed:
1. `trim.py` emits `flank_R1.fq.gz` + `flank_R2.fq.gz`, preserving mate pairing.
2. `genome_align.py` calls a new `Aligner.map_paired_to_genome(r1, r2, orphans)` for the paired flanks + a single-end pass for orphans.
3. `align.py:map_paired_to_genome` runs `minimap2 -ax sr ... ref r1 r2` (minimap2 already supports paired short reads when both files are passed positionally).

**Why it helps:** for a flank-bearing mate too short to anchor uniquely, the OTHER mate (full-length, untrimmed because the TE wasn't in it) provides the anchor. R2 gets this for free via bwa-mem's paired-end mode.

**Risks**
- Trim output schema change ripples through the cluster cache.
- Orphan flanks (where the mate didn't survive trim) still need a single-end pass.
- Run-time may rise modestly.

**Concrete tasks** (only flesh out when triggered):
- 4a. Failing test: `tests/genome_align_paired_test.py` builds a tiny BAM with one flank read whose mate has a unique anchor; asserts the flank lands at the anchor.
- 4b. `trim.py` change: emit paired FASTQs.
- 4c. `align.py` change: `map_paired_to_genome`.
- 4d. `genome_align.py` change: invoke paired + orphan single-end.
- 4e. Re-run validation (delete genome BAM + trim outputs); re-run recall diff.

---

## Task 5 (conditional, with sign-off): Option C — bwa-mem genome aligner

**Trigger:** Tasks 2+4 together still leave > 5% of the 318 baseline mismatches.

**Open a separate plan.** Decision points: add `bwa-mem` to `pixi.toml` (PLAN.md's minimap2-only goal is broken), expose `--genome-aligner=bwa-mem`, mirror R2's exact `bwa mem -k 15 -T 10` invocation. Do **not** start this task without an explicit user "yes" — the dependency change is a real policy shift.

---

## Conventions

- One commit per task. Pre-commit (ruff check/format) must pass.
- All `pixi run …` commands from `RelocaTE3/`.
- The validation harness has no `--force` flag on align-genome; force re-runs by deleting the per-sample `*.flank.genome.sorted.bam*` files (and the characterize TXT downstream).
- Do **not** touch `_classify` in `characterize.py`. The thresholds are R2-faithful and any change there is a divergence, not a fix.
- If Task 2's matrix shows a *new* cell appearing (e.g. `hom → het`), **stop** and investigate before proceeding. That would indicate the tuned minimap2 is now over-aligning and producing spurious junction reads, not just recovering missing ones.

## Repository pointers for a fresh session

- Branch: `fix-tsd-bug` (the TSD work landed in this branch; the genotype-status work continues here or on a follow-up branch — your call).
- Baseline status matrix: `RelocaTE3/validation/real_rice/report/characterized/status_confusion.tsv` (post the 2026-06-26 wildcard-TSD + single-sided-emission run).
- Classifier: `RelocaTE3/src/RelocaTE3/characterize.py:275` (`_classify`) — DO NOT EDIT.
- Spanner counter: `RelocaTE3/src/RelocaTE3/characterize.py:225` (`_count_spanners`) — verified equivalent to R2.
- Genome re-align step: `RelocaTE3/src/RelocaTE3/align.py:252-298` (`map_reads_to_genome`).
- R2 reference: `references/RelocaTE2/scripts/characterizer.pl` (149-185 for the classifier; 91-145 for the spanner counter).
- R2 alignment: `references/RelocaTE2/scripts/relocaTE_align.py` — uses `bwa mem -k 15 -T 10`.
- Prior PLAN.md notes on the recall gap: lines 206-212, 260-262, 277-283.
