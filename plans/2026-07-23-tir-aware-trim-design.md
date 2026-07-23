# TIR-Aware Trim (recover both insertion sides for TIR TEs) — Design

Date: 2026-07-23
Status: design (for review)
Scope: RelocaTE3 trim step (`src/RelocaTE3/librelocate.py`). Downstream
clustering / TSD inference (`insertions.py`) is unchanged.
Related: `todo/tsd-single-sided-tir-junctions.md`

## Goal

Recover the second junction side of insertions into TE elements with **terminal
inverted repeats** (TIRs), so RelocaTE3 measures the TSD instead of emitting
`UNK`. Target: TSD-exact on the mPing sim panel rises from ~0.81 to ~0.98,
matching RelocaTE2.

## Problem (summary; full evidence in the todo)

mPing has a ~15 bp TIR, so a read crossing one TE junction aligns to **both** TE
ends. RelocaTE3's trim keeps a **single best TE alignment per read**
(`_parse_te_bam` ~L87, `trim_TE_reads` L341–343, tie-broken by `_is_better` L347:
boundary contact then longer match) and then classifies that one alignment as 5′
**or** 3′ (`_trim_record` L269 vs L285, mutually exclusive). When the longer hit
is to the "wrong" TE end, the correct-side junction flank (`:start:3` / `:end:3`)
is never emitted → the insertion is one-sided → `UNK` TSD. RelocaTE2 emits the
correct flank because it does not collapse to one alignment per read the same way.

## Current data flow (what to change)

```
map_te_library ─► TE-library BAM(s)  (aligner reports BOTH TIR hits per read)
   └─ _parse_te_bam(bam)        -> {read: best rec}       ← collapse #1 (per BAM)
   └─ trim_TE_reads(bams)       -> {read: best rec}        ← collapse #2 (across BAMs)
        └─ _write_direction     -> for read, rec: _trim_record(rec)  ← one flank max
             └─ _trim_record    -> 5' (t_start<=2) XOR 3' (t_end>=t_len-3)
```

Both collapses use `_is_better`. The correct other-end alignment is discarded
**before** `_trim_record` ever sees it.

## Design decisions

- **Keep multiple TE-end alignments per read through junction extraction.**
  Instead of `{read: best rec}`, carry `{read: [recs]}` (or at least the best hit
  to *each distinct TE-end contact*). Then run the junction classification per
  retained alignment and emit a flank for each *valid* junction.
- **Validity gate stays as-is per alignment.** A retained alignment only produces
  a flank if it passes the existing branch conditions (`ends_align`,
  `match_span >= len_cutoff_m`, `flank_len >= len_cutoff_l`,
  `mismatch <= allowance`, and `t_start<=2` or `t_end>=t_len-3`). This naturally
  filters spurious TIR hits that don't form a real read-end/TE-end junction.
- **De-duplicate emitted flanks** by (read_name, side/label, trimmed-flank
  sequence) so a read isn't written twice for the same junction.
- **Faithful to R2, validated by parity** — R2 (`relocaTE_trim.py`) is the
  reference; the acceptance is matching R2's `:start:5/:end:5/:end:3/:start:3`
  emissions on the sim data, not a first-principles reimplementation.
- **Do not touch** `_make_insertion` / `_resolve_tsd` / clustering — the TSD math
  is correct once both sides are present.

## Open question to resolve during implementation (via R2 parity)

For a genuine 3′-junction read, the spurious 5′-TIR alignment shares the same
genomic flank sequence. We must confirm whether R2 emits **one** correct-side flank
or effectively both (and whether the downstream cluster de-dups). The parity test
(below) settles this empirically before finalizing the de-dup rule — e.g. keep,
per read, only junctions whose TE-end contact is at the read end *opposite* the
flank being written, so each read contributes at most one flank per physical side.

## Approach (incremental)

1. **Retain alternates:** change `_parse_te_bam` / `trim_TE_reads` to keep, per
   read, the best alignment to *each* TE end (5′-contact and 3′-contact) rather
   than one overall best. Minimal surface: return a small list/dict keyed by
   TE-end contact.
2. **Classify per alignment:** `_write_direction` iterates the retained
   alignments; `_trim_record` is called per alignment (it already returns one
   labelled flank or `None`).
3. **De-dup + write** as above.
4. Everything downstream (`read_repeat_name`, flanking FASTQs, genome re-align,
   `find-insertions`) is unchanged — it just now receives both sides.

## Testing

- **Unit (trim):** a synthetic TIR read with two TE-end alignments (long 5′-TIR
  hit + short correct 3′ hit) → assert the correct `:start:3`/`:end:3` flank is
  emitted (fails today).
- **Unit (no regression):** a normal single-end-contact read still yields exactly
  its one correct flank; a fully-internal read still `:middle`.
- **Parity / integration:** run trim on a small TIR fixture and compare emitted
  junction reads to RelocaTE2's for the same input; then, on the sim panel,
  confirm single-sided-`UNK` insertions drop and TSD-exact → ~0.98 with no
  precision regression. (The relocate-benchmark `matches.tsv` R2-vs-R3 comparison
  used in the diagnosis is the ready-made check.)

## Risks

- **Spurious double-emission.** Retaining both TE-end alignments risks writing two
  flanks for one physical junction. The validity gate + de-dup + the R2 parity
  test are the guardrails; resolve the exact de-dup rule against R2 output.
- **Over-calling / precision.** More junction reads could, in principle, create
  new low-support clusters. Watch precision/FDR in the parity run (currently
  ~1.0); it should hold.
- **Faithful-port discipline.** R3's trim is a port of R2's `relocaTE_trim.py`;
  the safest fix mirrors how R2 handles multi-hit reads rather than inventing new
  logic.

## Out of scope

- Changes to the TSD inference, clustering, or pairing.
- The blat backend read-chunking issue (separate; `relocate3-blat-bwa`).
- Any benchmark-side change (the benchmark already passes `--tsd UNK`).
