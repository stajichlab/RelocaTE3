# Long-read TE detection: algorithm proposal

2026-11-09. Status: **proposal; approach B approved, decisions 2–6 open**. Roadmap item after the
Nextflow pipeline (`plans/2026-11-02-nextflow-pipeline-design.md`). It
replaces the sketch in `plans/FEATURES.md`, Feature 3.

## Purpose

Choose the algorithm RelocaTE3 uses to find TE insertions in long reads (ONT
and PacBio). This document compares the possible approaches, recommends one,
and lists the decisions the maintainer needs to make before an implementation
plan is written.

## What changes with long reads

RelocaTE2, and RelocaTE3's short-read path, exist because a 100–150 bp read
cannot contain a TE plus enough flank to place it. The method therefore:
1. finds reads that touch a TE;
2. trims the TE part off;
3. places the short flank on the genome;
4. infers the insertion from junction pairs.

A long read removes that limit. One 10–20 kb read usually holds
flank + TSD + whole TE + TSD + flank. Three read situations occur:

| Read situation | What the alignment to the reference shows | Typical for |
|---|---|---|
| **Spanning**: the read crosses the whole insertion | one alignment with a large insertion (`I`) in the CIGAR, or two alignments of one read with unaligned read sequence between them | MITEs, SINEs, most DNA TEs, solo LTRs, short LTR elements |
| **Partial**: the read ends inside the insertion | a long soft clip at one end, the clip being TE sequence | LTR retrotransposons of 5–15 kb, low-N50 libraries |
| **Reference TE absent**: the sample lacks a TE the reference has | a deletion (`D`), or a split alignment, matching an annotated TE | excisions, reference-only insertions |

**The first-principle rule changes here.** "Port RelocaTE2's logic exactly"
cannot apply: RelocaTE2 has no long-read logic. For this feature, parity is
replaced by three gates (see "Validation"):
- truth on a simulated panel;
- agreement with RelocaTE3's own short-read calls on the same samples;
- comparison with published long-read callers.

The short-read path is not touched.

## Approaches considered

| | Approach | How it works | Used by | Strengths | Weaknesses |
|---|---|---|---|---|---|
| A | **TE-first** (keep RelocaTE's shape) | Align reads to the TE library, cut out TE parts, align the flanks to the genome, pair junctions | RelocaTE2/3 short-read path | Reuses the existing code and outputs | Most TE-containing long reads are *reference* copies, so nearly every read in a TE-rich genome becomes a candidate. It throws away the read's own proof of the insertion (the contiguous alignment), and nested TEs break the one-TE-per-read trim |
| B | **Genome-first signatures** | Align reads to the genome once; collect insertion, clip and split signatures; cluster; test the inserted sequence against the TE library; build a consensus; genotype | xTea (long-read mode), TLDR, PALMER, sTELLeR, LOCATE | Uses the read as direct evidence; one framework for all three read situations; observes the TSD and full TE sequence; cost scales with *non-reference* sites, not with TE content | Depends on the aligner representing insertions well in repeats; needs its own clustering and consensus code |
| C | **SV-caller-first** | Run a general SV caller (Sniffles2, cuteSV), then keep insertions whose sequence matches the TE library; local assembly per site | TELR / STELR, GraffiTE (read mode) | Least new code; mature callers | Callers are tuned for germline SVs ≥ 50 bp with several supporting reads. Low-frequency (somatic) insertions and partial-read evidence are filtered before we see them. Two extra heavy dependencies; TSD and TE structure still need our own step |
| D | **Assembly-first** | Assemble the sample, align the assembly to the reference, call insertions from the assembly | GraffiTE (assembly mode), pangenome pipelines | Best sequence resolution; natural for pangenomes | Needs ≥ 20–30× and an assembly per sample; loses read-level allele frequency and somatic calls; belongs with the later pangenome milestone, not here |

Evidence from published benchmarks:
- **Systematic benchmark, bioRxiv 2025.09.29.679192:** it compared PALMER,
  TLDR, xTea and sTELLeR. All four are approach B. On simulated reads, F1 was
  96.3–97.8% for xTea, sTELLeR and PALMER, and 94.4% for TLDR. xTea was best:
  95.8% recall, 99.9% precision. Accuracy rose with coverage up to 40× for
  xTea and sTELLeR.
- **LOCATE, bioRxiv 2025.02.26.640385:** also approach B. It reports beating
  nine other callers, including TELR, TrEMOLO and GraffiTE (approach C or D),
  and assembles each insertion.
- **Caveat:** these benchmarks are human (Alu, L1, SVA). None covers plants or
  MITEs. TLDR's default minimum insertion is 200 bp; mPing is 430 bp, and many
  MITEs are shorter than 200 bp. This is the main reason not to simply wrap an
  existing tool.

## Recommendation: approach B, with a per-site consensus

**Genome-first signature calling, TE-library annotation, consensus per site,
read-count genotyping.** In outline:

```
reads ──► 1 align to genome (minimap2, technology preset)
          2 collect signatures per read: large insertion, long clip, split, deletion
          3 cluster signatures by breakpoint
          4 annotate each cluster's sequence against the TE library
          5 consensus + realign  → exact breakpoint, TSD, TE sequence
          6 filter
          7 genotype from supporting vs reference-spanning reads
          8 write the usual RelocaTE3 tables (+ consensus FASTA, later VCF)
```

### Step 1. Align reads to the reference
- **Aligner:** minimap2 (already a dependency).
- **Preset per `--read-type`:**
  - `ont`: `map-ont`;
  - `ont-hq` (R10 / Q20+): `lr:hq`;
  - `hifi`: `map-hifi`;
  - `clr`: `map-pb`.
- **Flags:** `-Y` (soft clips on supplementary records, so clipped sequence is
  kept), `--MD`, `--secondary=no`. Output a sorted BAM/CRAM.
- **Existing alignments:** accept an existing long-read BAM/CRAM instead of
  FASTQ. Most long-read users already have one.

### Step 2. Collect signatures (one pass over the BAM)
Each primary or supplementary alignment with MAPQ ≥ a threshold contributes
some of these signatures:
- **INS:** a CIGAR insertion ≥ `--min-ins-len` (default 50 bp), with its
  reference position and the inserted bases.
- **Split-INS:** two alignments of one read, same chromosome and strand,
  adjacent on the reference (gap or overlap ≤ about 50 bp) but separated on
  the read by ≥ `--min-ins-len`. The read bases between them are the insert.
  minimap2 reports long insertions this way, especially in ONT.
- **CLIP:** a soft clip ≥ `--min-clip-len` (default 100 bp) at a read end,
  with its reference position, side and clipped bases.
- **DEL / Split-DEL:** a deletion ≥ 50 bp, for reference-TE absence.
- **REF-SPAN:** the read aligns across a position with ≥ `--min-flank`
  (default 100 bp) on both sides and no signature there. These are counted
  later, per candidate site, not stored per base.

### Step 3. Cluster
- Sort signatures by chromosome and position.
- Join INS and Split-INS whose positions are within `--cluster-window`
  (default 50 bp; wider for `ont`, because homopolymer errors shift
  breakpoints).
- Attach CLIP signatures within the window. Right-clips and left-clips are
  tracked separately, as the two junction sides.
- Split a cluster when insert lengths disagree by more than a tolerance (two
  different insertions at nearby positions) or when sequences don't align to
  each other.

### Step 4. Annotate against the TE library
- Align each signature's inserted or clipped sequence to the TE library
  (minimap2, same preset family; BLAT is not needed at these lengths).
- Per read, record family, strand, TE start and end, identity, and the
  fraction of the insert covered by TE.
- A cluster becomes a **candidate TE insertion** when ≥ 80% of its reads have
  a TE hit (LOCATE's rule) and a consistent family and strand.
- Classify structure from TE coordinates: full length, 5′-truncated, solo
  LTR, or **nested / composite** (more than one family in the insert,
  reported as an ordered list).
- **Mismatch gate:** identity is a percentage with a per-technology default.
  The short-read "≤ 2 mismatches" rule is meaningless here.

### Step 5. Consensus and exact placement
- **Spanning clusters:**
  1. Take each supporting read's segment (insert ± about 500 bp of flank).
  2. Build a partial-order-alignment consensus (abPOA or spoa; both are on
     Bioconda and have Python bindings).
  3. Realign the consensus to the reference.
  - Breakpoint: where the consensus leaves and rejoins the reference.
  - TSD: the reference bases present on **both** sides of the insert in the
    consensus. It is observed, not estimated.
  - TE sequence: the consensus insert, written to a FASTA.
  - HiFi reads are accurate enough to skip the consensus for 1–2 read sites;
    ONT reads are not.
- **Partial clusters** (clips only):
  1. Build a consensus for each side.
  2. Take the breakpoint from the clip positions.
  3. If left and right junctions overlap, that overlap is the TSD.
  - This is the same geometry as the short-read two-sided call. It is the
    only evidence for elements longer than the reads.

### Step 6. Filters
- **Minimum supporting reads:** default 2 for `hifi`, 3 for `ont`, scaled up
  with depth. `--min-support 1` is allowed for somatic work and labelled as
  such.
- **TE content of the insert:** ≥ 80% of the insert must be TE. This drops
  tandem-repeat expansions and non-TE insertions.
- **Flank quality:** flanking MAPQ and flank-length minimums on at least one
  read per side.
- **Two-sided rule for partial clusters:** both sides must have clips, the
  long-read counterpart of `--require-both-junctions`. One-sided clusters are
  reported in a lower tier and never dropped silently.
- **Low-complexity context:** flagged, not removed.
- **Insertions inside reference TEs** are kept and annotated. They are real
  and common in plants; hiding them would be a regression against the
  short-read path.

### Step 7. Genotype
- For each call, count **alt** reads (insertion or clip support) and **ref**
  reads (REF-SPAN at that position).
- Allele frequency = alt / (alt + ref), with a binomial confidence interval.
- Class by genotype likelihood:
  - homozygous;
  - heterozygous;
  - `somatic_insertion` (low frequency).

  The class names are the ones `characterize` already writes, so cohort
  tables and reports work unchanged.
- **Reference-TE absence:** for DEL clusters that match an annotated
  reference TE (boundaries within a tolerance), the same counts give
  present / absent / heterozygous. This is the long-read counterpart of
  `find-reference` plus excision detection.

### Step 8. Outputs
- The existing per-sample files, so `cohort`, MultiQC and the Nextflow
  pipeline need no new formats:
  - `all_nonref_insert.txt/.gff`;
  - tiers;
  - `characTErized.txt/.gff`.
- New columns:
  - support counts by signature type;
  - allele frequency;
  - TE structure class;
  - insert length;
  - consensus identity.
- New files:
  - `<sample>.te_insertions.fa`: consensus insert sequences;
  - the reference-absence table.
- VCF stays deferred to the cohort / VCF design, as agreed for the pipeline.

## Why this is the best fit

1. **It uses the strongest evidence available.** A spanning read shows the
   insertion, both junctions, the TSD and the element in one molecule.
   Approach A cuts that read into pieces and re-infers what was already
   visible.
2. **Its cost follows the answer, not the genome.** Work is proportional to
   the number of non-reference signatures. Approach A's work is proportional
   to the number of reads touching any TE, which is most reads in rice or
   maize.
3. **It is the approach the best-scoring published tools share** (xTea,
   sTELLeR, PALMER, LOCATE; F1 96–98% simulated). Writing it in RelocaTE3 lets
   us set limits for plant TEs: small MITEs, long LTR elements, nested
   insertions.
4. **One framework covers all three read situations,** including elements
   longer than the reads and reference-TE absence.
5. **It keeps allele frequency at read level,** which SV-caller-first and
   assembly-first lose. The lab's somatic mPing work needs it.
6. **It fits what exists.** Only minimap2, samtools/pysam and one small
   consensus library are needed. The pipeline, cohort and report layers stay
   as they are; a long-read sample becomes another kind of row.

What it costs: new clustering, consensus and genotyping code with no
RelocaTE2 reference to check against. That is why the validation gates below
are part of the proposal, not an afterthought.

## Validation (replaces RelocaTE2 parity for this feature)

| Gate | Data | Passes when |
|---|---|---|
| 1. Truth | The simulated rice panel (mPing, riceTElib) re-simulated as ONT and HiFi reads at 5×, 15×, 30× from the same genomes and truth sets. The simulator needs a long-read mode, e.g. pbsim3 or badread; check what the panel generator supports | Precision and recall against truth, by coverage, family and length. Targets are set after a baseline: at least the short-read numbers at equal coverage, with TSDs exact for spanning calls |
| 2. Self-consistency | The same simulated samples through the short-read path | Long-read calls contain the short-read true positives; every disagreement is examined |
| 3. Other tools | TELR/STELR and a general SV caller + TE filter on the panel. Human-tuned tools (xTea, PALMER, TLDR) only if they run on rice | RelocaTE3 at or above them on the plant panel; losses explained |
| 4. Real data | A rice sample with both short and long reads, if the lab has one (to confirm with the maintainer) | Agreement with the sample's short-read calls; manual check of discordant sites |

The benchmark repository (`relocate-benchmark`) is the place for gates 1–3,
as with the short-read parity work.

## Proposed milestones

| # | Delivers | Gate |
|---|---|---|
| L1 | Long-read simulated panel (ONT + HiFi; 5×/15×/30×) with truth; scoring in the benchmark | panel reproducible; a general SV caller + TE filter gives a baseline number |
| L2 | `--read-type`, genome alignment presets, signature collection, clustering (steps 1–3) as tested subcommands | signatures recover the simulated sites (recall of *sites*, before any filter) |
| L3 | TE annotation + filters (steps 4, 6): first complete calls for spanning insertions | gate 1 on spanning calls |
| L4 | Consensus, exact breakpoints, TSDs, insert FASTA (step 5); partial-read (clip) calls | TSD exactness; long-LTR recall |
| L5 | Genotyping and reference-TE absence (step 7); outputs and cohort compatibility (step 8) | allele-frequency accuracy on simulated heterozygous and low-frequency sites |
| L6 | Nextflow: long-read rows in the sample sheet, modules, test profile; gates 2–4; docs | pipeline output equals the CLI's; release candidate |

Each milestone is one branch and pull request, TDD, as with the pipeline.

## Decisions for the maintainer

1. **Approach. DECIDED 2026-11-09: approach B**, genome-first signatures
   with consensus. The question asked was: accept approach B?
   The alternative worth discussing is C, wrapping Sniffles2 or cuteSV. It is
   faster to a first result, but weaker for somatic and partial-read calls,
   with two more dependencies.
2. **Technologies in the first release.** ONT (R10 / Q20+) and PacBio HiFi
   both from the start, with older ONT R9 and PacBio CLR as "supported, less
   tested"? Or one technology first? Recommendation: ONT-HQ and HiFi
   together; they share everything but presets and thresholds.
3. **Consensus engine.** abPOA or spoa through their Python bindings (small,
   fast, on Bioconda)? Or a heavier polisher (racon/medaka) for old ONT data?
   Recommendation: abPOA/spoa; revisit only if R9 data matters to the lab.
4. **Real data for gate 4.** Which rice samples have long reads, and are there
   matched short reads?
5. **Mixed short + long reads for one sample.** Out of scope for the first
   release (recommended). Each read type gives its own call set, and a later
   step can merge them.
6. **Command shape.** One CLI with `--read-type` (recommended; the user-facing
   surface stays one tool) and separate internal modules for the long-read
   path, so the validated short-read code is not edited.

## Sources

- Systematic comparative benchmarking of long-read TE detection (PALMER, TLDR, xTea, sTELLeR): https://www.biorxiv.org/content/10.1101/2025.09.29.679192.full.pdf
- LOCATE: https://www.biorxiv.org/content/10.1101/2025.02.26.640385.full.pdf
- sTELLeR (Bioinformatics 2024): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11601167/
- TELR (NAR 2022): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9757076/ ; STELR: https://github.com/bergmanlab/STELR
- TLDR (Ewing et al., Mol Cell 2020): https://nanoporetech.com/resource-centre/nanopore-sequencing-enables-comprehensive-transposable-element-epigenomic-profiling
- McClintock 2 (reproducible TE-caller evaluation): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10347736/

Not verified for this proposal, to check during L1:
- the exact default thresholds of each tool beyond those quoted above;
- which of the human-focused tools run on a plant genome and TE library;
- the long-read mode of the lab's panel simulator.
