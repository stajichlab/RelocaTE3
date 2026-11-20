# Long-read L3: TE annotation, filters, first calls (design)

**Date:** 2026-11-19 (America/Los_Angeles)
**Status:** approved by the maintainer 2026-11-19; implementation plan next.
**Parent:** `plans/2026-11-09-long-read-algorithm-proposal.md`, milestone L3
(steps 4 and 6). Builds on L2 (`docs/2026-11-18-long-read-l2.md`).

## Purpose

Turn L2's candidate sites (clusters) into the first TE insertion calls:
annotate every read's inserted sequence against the TE library, vote a
family per site, and apply filters. Spanning insertions only are gated;
exact breakpoints, TSDs and clip-only placement are L4.

## Decisions (2026-11-19)

| Question | Decision |
|---|---|
| Gate | At 20x on ONT-HQ and HiFi: recall >= 0.95 and precision >= 0.95. A match is a call within 50 bp of a true site with the correct family (greedy 1:1). Strand accuracy and 5x are reported, not gated |
| Clip-only clusters | Annotated and written as tier `clip_only`; not scored by the gate; their hits and misses reported separately |
| TE matcher | minimap2 (existing dependency), one run per sample, flags tuned for short queries; a unit test proves a 145 bp insert is hit |
| Structure | Approach A: two subcommands, `lr-annotate` and `lr-call`. Rejected: one `lr-call` that runs minimap2 internally (no per-read evidence table, refiltering reruns alignment); one representative read per cluster (cannot apply the 80% rule) |
| Minimum support | 2 reads for `hifi` and `ont-hq`; 3 for `ont` and `clr` |

## Data flow

```
lr-align -> lr-signatures -> lr-cluster          (L2, unchanged except reference id)
                 |               |
                 v               |
           lr-annotate           |   signatures + te_library.fa -> <sample>.te_hits.tsv.gz
                 |               |
                 +-----> lr-call <+  clusters + signatures + hits -> <sample>.calls.tsv.gz
```

- New modules: `longread/annotate.py`, `longread/call.py`.
- New table schemas in `longread/tables.py`: `lr-te-hits`, `lr-calls`
  (versioned, `.partial` then rename, refuse to overwrite).
- Both stages write a `timing.tsv` row via `stage_timer`.
- Short-read code is not touched.
- L2 known limit fixed here: `lr-align` records the reference id in the BAM
  header (`@CO`); `lr-signatures` refuses a mismatch with `--reference-id`.
- Deferred: read groups / multi-sample BAMs; short-read format outputs
  (`all_nonref_insert.txt`, etc.) stay in L5.

## `lr-annotate`

- Each signature's sequence (inserted bases; clipped bases for CLIP) goes to
  a temporary FASTA; minimap2 aligns it to the TE library with `-c`.
- Flags: start from the read-type preset, lower k-mer size for short
  queries. The final flags and why are recorded in the results note.
- One row per hit: `signature_id`, `read_id`, `family`, `te_strand`
  (TE orientation on the reference-forward strand), `te_start`, `te_end`,
  `te_length`, `query_start`, `query_end`, `identity`.
- A signature with no hit gets one row with `family=.` ("looked, nothing
  found" is distinct from "not looked").

## `lr-call`

Per signature (computed from hits):
- best family = most aligned bases;
- insert TE coverage = union of hit query spans / insert length;
- `nested` = more than one family covering disjoint parts of the insert
  (flag only; structure classes are later).

Per cluster:
1. Family vote: >= 80% of reads have a TE hit; the top family holds >= 80%
   of those reads; strand by majority, agreement fraction reported.
2. Tier: `spanning` if >= 1 INS or SPLIT_INS signature, else `clip_only`.
3. Filters: minimum reads (per read type, above); median insert TE
   coverage >= 0.8 for spanning calls; minimum hit identity per read type
   (tuned once on the panel, recorded).
4. Merge calls of the same family within the cluster window into one call
   (fixes the L2 limit where length splitting leaves two clusters at one
   site). Different families nearby stay separate.
5. Nothing is dropped: every cluster yields a row with `filter=PASS` or a
   fixed-order list of reasons (`no_te_hit`, `family_discordant`,
   `low_support`, `low_te_fraction`).
6. Position = cluster median breakpoint (exact breakpoints are L4).

Calls columns: `call_id`, `reference_id`, `sequence_id`, `position`, `tier`,
`filter`, `family`, `strand`, `family_fraction`, `strand_fraction`,
`n_reads`, `n_ins`, `n_split_ins`, `n_clip_left`, `n_clip_right`,
`median_insert_length`, `median_te_coverage`, `median_identity`, `nested`,
`cluster_ids`.

## Testing

Unit (TDD, hand-built, seconds):
- annotate: 145 bp mGing insert hit; reverse-strand insert -> `te_strand=-`;
  random insert -> `family=.`; two-family insert -> `nested`; CLIP annotated
  from clipped bases.
- call: 80% vote pass and discordant; min support per read type; TE-content
  filter; tiering; same-family merge vs different-family separation; every
  cluster yields a row; reason order.
- tables: round trip for both schemas; wrong schema/version refused.
- reference id: mismatch refused.
- CLI: both subcommands registered, refuse overwrite, write timing.

Panel: one end-to-end test on `tests/data/longread_small` (HiFi 5x, four
stages, calls carry families); skips when the panel is not linked.

## Gate: `validation/longread_l3/`

- `run_gate.slurm` (sbatch only): L2 stages (skipped when output exists),
  then `lr-annotate` and `lr-call`, on 4 samples + 4 reference-only controls.
- `score_calls.py`: greedy 1:1 within 50 bp with correct family; recall,
  precision, F1, strand accuracy, per-family table, `clip_only` hits and
  misses, PASS calls on controls (pure false positives).
- Pass: 20x ONT-HQ and HiFi recall >= 0.95 and precision >= 0.95.
- Results note: `docs/2026-11-xx-long-read-l3.md` (gate table, timing,
  minimap2 flags and why, known limits).

## Execution

Worktree `RelocaTE3-lr-l3`, branch `feat/long-read-l3`; one subagent per
task with reviews between; one PR.
