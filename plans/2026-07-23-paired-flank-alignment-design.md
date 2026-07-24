# Paired-End Junction-Flank Re-Alignment — Design

Date: 2026-07-23
Status: design (root cause confirmed by controlled trace)
Scope: RelocaTE3 step 4 genome re-alignment (`src/RelocaTE3/genome_align.py`).
Clustering / TSD inference (`insertions.py`) and trim (`librelocate.py`) unchanged.
Related: `todo/tsd-single-sided-ambiguous-flanks.md`

> **Supersedes `2026-07-23-tir-aware-trim-design.md`** (TIR + trim-collapse), which
> a controlled read-level trace disproved.

## Goal

Recover the second junction side of insertions whose flank maps ambiguously, so
RelocaTE3 measures the TSD instead of emitting `UNK`. Target: TSD-exact on the
mPing sim panel ~0.81 → ~0.98, matching RelocaTE2.

## Confirmed root cause

RT3 re-aligns trimmed junction flanks **single-end** (`align_to_genome(...,
paired=False)`). Junction flanks are short and often map **ambiguously (MAPQ=0)**
to several equal-scoring loci; single-end, RT3 scatters them off the true site, so
the true insertion loses one side → single-sided → `UNK` TSD. RelocaTE2 aligns
each flank **paired with its mate**; the mate maps **uniquely** at the true site
and anchors the flank (proper-pair rescue). Evidence (Chr1:11761374, a true
homozygous mPing) is in the todo; key: RT3 flanks are `flag=16` MAPQ=0 at non-true
loci, RT2's are proper pairs (`flag=83/163`, MAPQ 29/37) at the true locus. The
effect is aligner-invariant, consistent with a step-4 (not aligner) cause.

## Current data flow (`genome_align.py::align_to_genome`)

```
flanking/*.flankingReads.fq  ─┐
recover_support_mates()       ─┤─► ONE fastq ─► backend.map_genome(paired=False) ─► genome.bam
                                     (flanks + mates aligned SE, side-by-side)
```

The mates are recovered but aligned as separate SE reads — never *paired* with the
flank, so they cannot anchor it. (The module docstring makes this an explicit,
now-wrong, design choice.)

## Design

Align each junction flank **paired-end with its genomic mate**:

```
for each junction flank read F (base name B, mate end e):
    mate M = the /~e read of B  (the genomic mate; already found by recover_support_mates)
    if M exists:  emit (F, M) into paired R1/R2 files (matched order/name)
    else:         emit F into an SE fallback file
map_genome(paired R1/R2, paired=True)  +  map_genome(SE fallback, paired=False)  -> merge
```

- **Name matching:** strip the junction tag from F so its base matches M; write R1/R2
  in the same order so the aligner pairs them (bwa-mem/minimap2 pair by order+name).
- **Flank stays the read whose breakpoint we want:** F keeps its `:start:3` tag in
  the BAM (the insertion finder reads the tag); only its *placement* is anchored by M.
- **Fallback SE** for flanks whose mate also matched the TE (both ends are junctions)
  or is absent — unchanged behavior for those.
- **Merge** the paired BAM + SE BAM into the single coordinate-sorted `genome.bam`
  the finder consumes. `fullreads_bam` (false-junction filtering) is unchanged.

Downstream: the insertion finder clusters by tag + genomic position — it does **not**
need proper-pair flags — so it simply now finds the flank at the correct locus. No
`insertions.py` change.

## Testing (TDD)

- **RED (integration):** build a tiny synthetic genome with a short duplicated
  region (so a flank maps to 2 loci at MAPQ 0) at only one of which sits the
  insertion; give the flank a mate that maps uniquely at the true locus. Assert
  `align_to_genome` places the flank at the true locus (today: it scatters).
- **Unit:** the pairing builder emits matched R1/R2 for flanks-with-mates and SE
  for flanks-without-mates, with correct names/tags.
- **No-regression:** existing `genome_align` / `pipeline` tests stay green; a
  flank with no mate still aligns (SE) as before.
- **Acceptance (benchmark, out of band):** single-sided-`UNK` insertions drop;
  TSD-exact → ~0.98; precision (~1.0) holds.

## Risks / open questions

- **Pairing correctness:** flanks and mates must be matched exactly by base name and
  written in lockstep; a mismatch would mis-anchor. The unit test guards this.
- **Which reads get a mate:** confirm against RT2 that the set of flanks RT2 pairs
  matches what `recover_support_mates` yields (it recovers the mate only when the
  other end did not match the TE — the correct set for anchoring).
- **Insert-size expectations:** map_genome paired mode should not over-constrain
  insert size (junction flank + mate can be short); verify the aligner flags don't
  drop valid proper pairs. Compare placement to RT2 on the sim data.
- Secondary (not this fix): minimap2 TE-library sensitivity misses ~3/10 junction
  reads blat catches — a separate item.

## Out of scope

- Trim / TSD inference / clustering changes.
- The blat backend read-chunking issue.
- Benchmark-side changes (already passes `--tsd UNK`).
