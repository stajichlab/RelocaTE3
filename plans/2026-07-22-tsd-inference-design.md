# Variable-Length TSD Inference (TSD-unknown mode) — Design

Date: 2026-07-22
Status: design (for review); implementation plan in `2026-07-22-tsd-inference-plan.md`
Scope: RelocaTE3 (the tool). A follow-up one-line benchmark change (RT3 callers
use `--tsd UNK`) is downstream and out of scope here.

## Goal

Let RelocaTE3 **infer the target-site duplication (TSD) length and sequence from
the read data**, instead of requiring a fixed-length TSD motif. This matches
RelocaTE2's default behavior and fixes RT3's TSD accuracy on TEs whose TSD length
varies per insertion.

## Problem / evidence

In the relocate-benchmark (simulated mPing panel, TSDs 3–5 bp), TSD-exact accuracy
is:

| caller | TSD-exact / detected |
|---|---|
| relocate2 | 2766/2819 = **0.981** |
| relocate3-minimap2/minimap2 | 829/2767 = 0.300 |
| relocate3-bwa/bwa | 833/2814 = 0.296 |
| relocate3-bowtie2/bwa | 836/2821 = 0.296 |

The ~0.30 is identical across aligners, so it is a RelocaTE3 core issue, not an
alignment effect. Per-event strings show why: truth TSDs are variable length
(`AAGGG`, `GCAT`, `TGCAA`, `TGT`, …); RT2 recovers the full sequence exactly,
while RT3 emits only a fixed 3 bp (`AAGGG`→`GGG`, `GCAT`→`GCA`), matching exactly
only when the truth happens to be 3 bp (~30% of cases).

Root cause: the benchmark runs RT3 with `--tsd "..."` (a fixed 3-bp wildcard), and
RT3 **explicitly does not support TSD-unknown inference** — `InsertionFinder.find_insertions`
raises `NotImplementedError` for `UNK` (`src/RelocaTE3/insertions.py:71-80`). RT2's
default is the opposite: it writes `regex.txt` with the TSD field = `UNK`
(`references/RelocaTE2/scripts/relocaTE2.py:346`), triggering its read-depth TSD
inference — which is why RT2 gets variable-length TSDs right.

## How RelocaTE2 infers the TSD (reference: `references/RelocaTE2/scripts/relocaTE_insertionFinder.py`)

RT2 runs the depth path when TSD == `UNK`:

- **`TSD_from_read_depth` (line 551):** iterates positional read clusters; within a
  cluster, splits junction reads into left-side and right-side groups and pairs
  them into candidate insertions.
- **`TSD_len_calculate` (line 852):** the core. For a cluster it builds two
  histograms over genomic coordinates:
  - `TSD_left[start]`  — start coords of reads whose trimmed TE end marks the
    **right** boundary of the site (RT2's `pos == 'right'`).
  - `TSD_right[end-1]` — end coords of reads marking the **left** boundary
    (`pos == 'left'`).
  It takes the max-support coordinate on each side and computes
  **`TSD_len = TSD_right_max - TSD_left_max + 1`** (lines 899–905): the number of
  overlapping bases between the two junction clusters = the TSD length. This is
  data-driven and needs no motif.
- **`tsd_finder` (line 843):** a read-depth refinement — counts contiguous
  positions whose depth ≥ `tsd_depth × read_count`, giving a depth-based TSD length
  as a cross-check/alternative.
- The TSD **sequence** is then taken from the overlapping read/reference bases at
  the inferred coordinate span (contrast the fixed-pattern path `TSD_check`
  (line 1165), which captures the motif match from the read).

## What RelocaTE3 has today (`src/RelocaTE3/insertions.py`)

- `find_insertions(...)` (line 45): refuses `UNK` (line 71), else clusters via
  `_cluster_reads` → `_assign_cluster` → `_tsd_check`, populating
  `te_insertions[event][tsd_start][tsd_seq]`, then `_merge_offset_starts`, then
  `_write_output`.
- `_tsd_check` (line 279): a faithful port of RT2 `TSD_check` (fixed-motif path).
- The junction clustering machinery (`_cluster_reads`, `_assign_cluster`, per-read
  start/end/strand/`:start|end:5|3` suffix parsing) already exists and mirrors the
  inputs RT2's depth functions consume.

**Missing:** the `UNK` read-depth path (`TSD_from_read_depth` + `TSD_len_calculate`
+ `tsd_finder`).

## Approach

Port RT2's depth-mode TSD inference into RT3 as the `UNK` path, reusing RT3's
existing cluster walk:

1. **Trigger:** `--tsd UNK` (RT3 already recognizes the `UNK/UKN/unknown`
   sentinel; just replace the `NotImplementedError` with the inference path). No
   new CLI flag.
2. **Length inference:** port `TSD_len_calculate` — per cluster, build the
   left/right boundary histograms from the already-parsed junction reads and take
   `TSD_len = right_max - left_max + 1`. Port `tsd_finder` as the depth-based
   cross-check.
3. **Sequence + position:** extract the TSD sequence at the inferred coordinate
   span, faithful to RT2 (Task 1 of the plan pins down RT2's exact source —
   read-derived vs reference-derived; if reference bases are needed, thread the
   genome FASTA into `InsertionFinder.find_insertions`, whose signature currently
   lacks it).
4. **Emit** into the same `te_insertions[event][tsd_start][tsd_seq]` structure so
   `_write_output` is unchanged and variable-length TSDs flow through.
5. Keep the fixed-motif path exactly as-is; `UNK` selects inference.

## Key decisions

- Reuse `UNK` as the trigger (RT2 parity) rather than a new flag.
- Faithful port validated against RT2 output on identical inputs (not a
  reimplementation from first principles) — RT2's depth code is intricate and its
  exact boundary/sequence conventions must be reproduced to match its accuracy.
- Benchmark rewire (`tsd = "UNK"` for the RT3 callers + re-pin) is a **separate,
  downstream** change, done after this merges.

## Testing strategy

- Unit: `TSD_len_calculate` and `tsd_finder` ports on synthetic clusters with
  known TSD lengths (3/4/5 bp) → exact expected lengths.
- Golden/parity: run RT2 (the reference tree) and RT3-`UNK` on the same small
  genome-aligned BAM; assert RT3's per-site TSD sequences match RT2's.
- Integration: `find-insertions --tsd UNK` on `tests/data` produces variable-length
  TSDs and a non-empty table.
- Acceptance (post-merge, benchmark): RT3 TSD-exact rises from ~0.30 toward RT2's
  ~0.98.

## Risks

- RT2's depth path is Python-2, print-heavy, and subtle (boundary off-by-ones,
  strand handling). The parity test against RT2 output is the guardrail.
- If the TSD sequence needs reference bases, `InsertionFinder.find_insertions` and
  its CLI caller gain a genome argument — a small signature change to plan for.
