# Long-read L2: alignment, signatures, clustering (design)

**Date:** 2026-11-17 (America/Los_Angeles)
**Status:** approved by the maintainer 2026-11-17; implementation plan next.
**Parent:** `plans/2026-11-09-long-read-algorithm-proposal.md`, milestone L2
(steps 1-3, and the "Designing for multiple references" rules).

## Purpose

First long-read code in RelocaTE3: align long reads to the genome, collect
insertion evidence per read, and cluster it into candidate sites. No calls
yet (TE annotation and filters are L3).

## Decisions (2026-11-17)

| Question | Decision |
|---|---|
| Where `--read-type` lives in L2 | Only on the new subcommands `lr-align`, `lr-signatures`, `lr-cluster`. `run-all`, `run-batch` and the sample sheet are untouched until L3 |
| Deletion signatures (DEL / Split-DEL) | Deferred to L5 (reference-TE absence). The schema's `kind` field leaves room |
| Gate | Site recall >= 95 % at 20x on both ONT-HQ and HiFi, a cluster within 50 bp of the true position. 5x recall and false clusters are reported, not gated |
| Structure | Approach A: a `RelocaTE3.longread` subpackage, three subcommands, versioned TSV.gz intermediates. Rejected: JSON Lines intermediates; one streaming command with no intermediates |

## Architecture: `src/RelocaTE3/longread/`

| Module | Holds |
|---|---|
| `model.py` | Frozen dataclasses `Locus`, `Reference`, `PlacementRecord`, `Signature`, `Cluster`; enums `ReadType`, `SignatureKind`, `ClipSide` |
| `reference.py` | `ReferenceSpace` interface (`sequences()`, `fetch(locus)`); `LinearFastaReferenceSpace` (pysam `FastaFile`) |
| `align.py` | Step 1: minimap2 preset per read type -> sorted, indexed BAM; checks for a supplied BAM/CRAM |
| `placements.py` | BAM adapter: alignments -> `PlacementRecord`. The only module touching pysam alignment objects |
| `signatures.py` | Step 2: placement records -> INS, SPLIT_INS, CLIP signatures |
| `cluster.py` | Step 3: signatures -> clusters |
| `tables.py` | Versioned TSV.gz readers/writers for signatures and clusters |
| `timing.py` | Stage timer: one row per stage in `<outdir>/timing.tsv` |
| `cli.py` | `lr-align`, `lr-signatures`, `lr-cluster`; registered from `RelocaTE3/cli.py` by one call |

No existing short-read module is edited except the one registration call in
`RelocaTE3/cli.py`.

## Data model

- `Locus(reference_id, sequence_id, start, end, strand)`: 0-based half-open.
  Every signature and cluster carries one. No dict keyed by chromosome alone.
- `Reference(id, fasta, index)`: the CLI takes one (`--reference-id`,
  default `ref`); internals take a list or a `ReferenceSpace`.
- `PlacementRecord`: read id, read length, query start/end (in read
  orientation, hard clips included), target `Locus`, `is_supplementary`,
  `is_reverse`, MAPQ, the large insertions (reference position, query
  interval) and soft clips (side, length). No pysam objects.
- `Signature`: `kind`, `Locus`, insert length, inserted or clipped sequence,
  ~200 bp read flank each side, read id, query coordinates, clip side.
- `Cluster`: explicit `cluster_id`, `Locus` (median breakpoint), counts per
  kind, insert lengths, left/right clip counts, member signature ids,
  distinct read count.

## Stages

**lr-align.** `--reads FASTQ...` or `--bam BAM/CRAM`; `--genome`,
`--read-type`, `--threads`, `-o`. Presets: `ont` map-ont, `ont-hq` lr:hq,
`hifi` map-hifi, `clr` map-pb; flags `-a -Y --MD --secondary=no`; then
`samtools sort` + index. A supplied BAM must be coordinate-sorted, indexed,
and its `@SQ` names/lengths must match the genome FASTA, else an error naming
the `samtools` fix. minimap2 2.31 in the env supports all four presets.

**lr-signatures.** Primary and supplementary records with MAPQ >=
`--min-mapq` (20).
- INS: CIGAR insertion >= `--min-ins-len` (50); insertions within
  `--ins-merge-dist` (20 bp) on the reference in one read are merged first.
- SPLIT_INS: consecutive alignments of one read in read order, same
  sequence and strand, reference gap/overlap <= `--split-max-ref-gap` (50),
  read gap >= `--min-ins-len`.
- CLIP: soft clip >= `--min-clip-len` (100), with side.
- Hard clips: coordinates exact; sequence from the read's primary record;
  the count of hard-clipped records is reported.
- Output header carries schema version, read type and reference id.

**lr-cluster.** Sort by sequence and position; join within
`--cluster-window` (default from the header's read type: 50 for hifi/ont-hq,
100 for ont/clr); split when insert lengths disagree by more than 30 %;
attach CLIPs to the nearest cluster, left/right counted separately; each read
counts once per cluster. Splitting by sequence dissimilarity is deferred to
L4 (consensus).

**Shared.** Outputs are written to `.partial` and renamed; existing outputs
are refused. Each stage appends a row to `<outdir>/timing.tsv` (stage, wall
time, peak RSS, records in/out, RelocaTE3 version). Missing tools, wrong read
types, and schema-version mismatches between stages are clear errors.

## Testing

- **Tier 0 (TDD):** hand-built `PlacementRecord`s for spanning insertion,
  split read, clip-only, reverse strand, merged ONT insertion, nested
  insert, hard clip; small pysam-built BAMs for the adapter; a tiny synthetic
  genome for each minimap2 preset.
- **Tier 1 (gate):** `validation/longread_l2/run_gate.slurm` runs the three
  stages on the eight panel read sets (via `tests/data/longread_small`);
  `validation/longread_l2/score_sites.py` scores site recall (cluster within
  50 bp of truth). Pass: >= 95 % at 20x on ONT-HQ and HiFi. Reported: 5x
  recall, unmatched clusters per sample, and all clusters on the
  reference-only controls (pure false positives). The scorer only reads
  results; nothing is simulated in RelocaTE3.

## Not in L2

DEL/Split-DEL (L5); TE annotation and filters (L3); consensus, sequence-based
cluster splitting, TSDs (L4); `run-all`/`run-batch`/sample-sheet
`--read-type` (L3); GAF input (later); CRAM-specific tuning.
