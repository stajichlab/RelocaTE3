# Long-read TE detection: algorithm proposal

2026-11-09, updated 2026-11-13. Status: **accepted design; all six decisions made (2026-11-13)**. Roadmap item after the
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
  2. Build a partial-order-alignment consensus with abPOA (`pyabpoa`,
     decision 3). The call sits behind one function,
     `consensus(sequences) -> sequence`, so spoa (`pyspoa`) can replace it
     if the simulated panel shows a difference.
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

## Designing for multiple references and pangenomes

Requested by the maintainer on 2026-11-09. The two roadmap items after long
reads are (1) several reference genomes given as separate FASTA files and
(2) a pangenome graph as the reference. The long-read feature ships with one
linear reference, but its data model must not assume there is only one. These
rules apply from milestone L2:

**1. A position always names its reference.**
- Use a `Locus` record, never a bare chromosome string: `reference_id`,
  `sequence_id` (contig, or later a graph path), `start`, `end`, `strand`.
- Every signature, cluster and call carries a `Locus`.
- In a graph, the same record holds a path name plus offset, or a node walk.
- No dictionary is keyed by chromosome name alone.

**2. References are a list, behind one interface.**
- Internally: `references: list[Reference]`, each with `id`, FASTA, TE
  annotation and aligner index. The CLI accepts one for now.
- A `ReferenceSpace` interface provides:
  - `sequences()`;
  - `fetch(locus)`;
  - `te_annotation(locus)`;
  - later, `project(locus, to_reference)` (liftover between FASTAs, or
    projection onto a graph path).
- A linear FASTA is the first implementation; a graph is a later one.

**3. Evidence is separated from the alignment format.**
- Signature collection reads **placement records**, not pysam objects:
  read id, query interval, target `Locus`, strand, mapping quality, identity,
  and the edit operations that matter (insertions, deletions, clips).
- A BAM adapter produces them now. A GAF adapter (minigraph, vg giraffe,
  GraphAligner) can produce the same records later.
- Steps 3–7 never see the file format.

**4. Signatures keep the read-side evidence.**
- Store the inserted sequence, about 200 bp of flanking read sequence on each
  side, the read name and the query coordinates.
- A site can then be re-placed on another reference by aligning its flanks.
  The whole sample does not need re-aligning.

**5. "Non-reference" is a property relative to a reference, not a kind of
call.**
- An insertion absent from reference A can be present in reference B.
- The internal call (`TEAllele`) holds:
  - the element: family, strand, structure, consensus sequence, TSD;
  - its flanking anchor sequences;
  - `placements: {reference_id: Locus}`;
  - `state: {reference_id: present | absent}`;
  - the sample genotype with its read counts.
- Today's file names (`all_nonref_insert`, `all_ref_insert`) become views of
  that model for one chosen reference.

**6. Site identity is explicit.**
- Each call gets a stable `site_id`, built from reference, locus, strand and
  family.
- The cohort matrix groups by `site_id`, not by re-deriving a key from
  columns.
- Multi-reference support can then define equivalence classes of `site_id`s
  (the same insertion seen against different references, or nearby
  breakpoints) without changing every consumer.
- This is also where the opposite-strand and nearby-site merging deferred in
  the pipeline milestones belongs.

**7. Outputs carry the reference.**
- New long-read columns include `reference_id` (one value for now).
- Intermediate files (signatures, clusters, calls) have a written schema with
  a version field.
- The existing short-read tables keep their columns unchanged.

**8. No single-genome global state.**
- Nothing in the long-read modules reads "the genome" from module-level
  state.
- A function that needs a reference takes a `Reference` or `ReferenceSpace`
  argument.

What this costs now: a few small dataclasses and one adapter layer. What it
avoids later: rewriting clustering, consensus and genotyping when the
reference stops being one FASTA.

## Candidates for a later Rust rewrite

Requested by the maintainer on 2026-11-09. Rust comes after the multi-reference
work on the roadmap, and nothing here is to be done now. This section records
where a rewrite is likely to pay off, and how to write the Python so the port
is cheap.

**Long-read path** (new code; estimates only until profiled on the 30× panel):

| Stage | Why it is likely hot | Rust route |
|---|---|---|
| Signature collection (step 2) | One pass over every alignment, parsing CIGAR strings and supplementary-alignment tags in a Python loop. This is the largest loop in the design | `noodles` or `rust-htslib`; one record in, signature records out |
| Reference-spanning read counts (step 7) | A second pass over the alignments around every candidate | Fold into the same Rust pass as signature collection, with an interval index |
| Clustering (step 3) | Sort and sweep over millions of signatures; cheap per item, large in total | Straightforward once signatures are plain records |
| TE annotation (step 4) | Many short alignments of inserts to the TE library. The cost is process start-up and file round trips, not the aligner | `minimap2-rs` bindings, in process |
| Consensus (step 5) | Partial-order alignment per site. abPOA/spoa are already C/C++; the Python cost is batching and conversion | Call the library from Rust; parallel over sites with `rayon` |

**Short-read path** (existing code; these were the measured or suspected
hotspots during the parity and memory work):

| Code | Evidence |
|---|---|
| `librelocate.py` `_parse_te_bam`, `write_trimmed_reads_chunked` (trim) | The 20 GB memory peak fixed by chunking; still a Python object per TE hit |
| `aligners.py` `_iter_psl_to_sam`, `_query_sequences` (BLAT PSL → SAM) | About 6 GB and minutes per side before chunking |
| `insertions.py` `_stream_clusters`, `_call_insertions`, `_estimate_tsd_length_from_depth` | Per-record loops over the flank BAM; named in `plans/FEATURES.md`, Feature 1 |
| `characterize.py` | One `samtools view` subprocess per insertion site |
| FASTQ → FASTA and chunking | Already moved to `seqtk` + `split`; a Rust reader would remove the temporary whole-side FASTA |

**Rules that keep the port cheap** (apply while writing the long-read code):
1. **Stage boundaries are files with a written schema** (see rule 7 above).
   A Rust binary can then replace one stage, for example
   `relocaTE3 lr-signatures`, and must produce the identical file.
2. **Each stage has golden-output tests.** Those tests become the parity gate
   for the Rust version, exactly as `run-all` output was the gate for the
   Nextflow pipeline.
3. **Inner loops are pure functions over small typed records.** That means
   dataclasses with primitive fields, with no pysam objects or nested dicts
   passed between stages.
4. **Measure before porting.** Each stage logs its wall time, peak memory and
   record counts. The pipeline trace already records per-task time and
   memory.
5. **Leave the choice of binding open:** a PyO3 extension module (the
   in-process `--backend rust` idea of `plans/FEATURES.md`) or standalone
   binaries behind the same subcommands. File-level stage boundaries make
   either possible.

## Validation (replaces RelocaTE2 parity for this feature)

| Gate | Data | Passes when |
|---|---|---|
| 1. Truth | The simulated rice panel (mPing, riceTElib) re-simulated as ONT and HiFi reads at 5×, 15×, 30× from the same genomes and truth sets. The simulator needs a long-read mode, e.g. pbsim3 or badread; check what the panel generator supports | Precision and recall against truth, by coverage, family and length. Targets are set after a baseline: at least the short-read numbers at equal coverage, with TSDs exact for spanning calls |
| 2. Self-consistency | The same simulated samples through the short-read path | Long-read calls contain the short-read true positives; every disagreement is examined |
| 3. Other tools | TELR/STELR and a general SV caller + TE filter on the panel. Human-tuned tools (xTea, PALMER, TLDR) only if they run on rice | RelocaTE3 at or above them on the plant panel; losses explained |
| 4. Real data | A rice sample with both short and long reads, if the lab has one (to confirm with the maintainer) | Agreement with the sample's short-read calls; manual check of discordant sites |

The benchmark repository (`relocate-benchmark`) is the place for gates 1–3,
as with the short-read parity work.

### Test data: what exists (checked 2026-11-13)

Four locations, given by the maintainer. Nothing in them was changed.

| Location | What it holds | State |
|---|---|---|
| `/rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data` | The simulation toolkit. Long-read simulation (PBSIM3, plus `ccs` for HiFi) is implemented in `longread.py` | Branch `feat/refactor-replan` with uncommitted work (molecule-level support accounting, 2026-10-12). Also holds an earlier plan, `plans/2026-09-02-relocate3-long-read-support.md` |
| `/rhome/nmath020/wessler_bigdata/rice/make_simulated_genome/make_simulation_new` | Simulation runs. `results/riceTElib_lr_benchmark` is the long-read panel: the 500-event riceTElib truth set, identical to the short-read panel, as ONT-HQ and HiFi at 5×/15×/30× × 3 replicates, each with a matched reference-only control | **15 of 18 tasks complete.** All ONT-HQ and HiFi 5×/15× are done. The three HiFi 30× tasks hit the 48-hour limit (partial FASTQ, empty control, no completion marker). Support tables predate the 2026-10-12 accounting fix. 134 GB |
| `/rhome/nmath020/bigdata/github/github_tools/RelocaTE/relocate_benchmark/relocate-benchmark` | The benchmark harness (callers, scoring, reports) | Short-read only: the manifest reader requires `r1`/`r2`, and no long-read dataset is registered |
| `/rhome/nmath020/bigdata/nanopore_data` | Real ONT reads, `all_fq_final/concat_fq/<line>.ont.fq.gz`, 52 lines. Read headers show `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`, so `--read-type ont-hq` | 33 `B_*`/`G_*` RILs have matched Illumina reads in `/bigdata/stajichlab/nmath020/pop2_and_pop3_analysis/find_TE/fq_input/trim`. Parents A119, A123 and HEG4 are present. EG4 and NB (Nipponbare, the reference line) are in `UCD_03.16.2026_EG4_NB/concat_fq`, not yet in `all_fq_final`. No real HiFi data |

### Test tiers (smallest first)

Each tier must be credible before the next is run. All computation goes
through `sbatch`, never the login node.

| Tier | Data | Run time target | Used for |
|---|---|---|---|
| 0. Unit | Hand-built reads and small BAMs, a few records each: spanning read, split read, clip-only read, reverse strand, nested insert, hard clip | seconds | every TDD cycle |
| 1. Small panel | **To build in L1.** The repository's Chr3 2 Mb fixture (`tests/data/sim_genome`, 200 insertions with truth) re-simulated as ONT-HQ and HiFi reads with the `simulate_data` toolkit. Stored outside git if it exceeds a few tens of MB; a script regenerates it | minutes per sample | the long-read acceptance test, the counterpart of `tests/acceptance_test.py`; fast iteration |
| 2. Full panel | `riceTElib_lr_benchmark` through `relocate-benchmark`, in cost order: one 5× ONT-HQ replicate, one 5× HiFi, all 5×, then 15× and 30× | hours | gates 1–3 |
| 3. Real | The lab's ONT lines | hours | gate 4 |

### Gate 4 samples (decision 4)

- **NB (Nipponbare):** the sample is the reference line, so non-reference
  calls should be near zero. A negative control on real reads.
- **Parents HEG4, EG4, A119, A123:** known mPing content.
- **A few RILs with matched Illumina reads,** for agreement with the
  short-read calls. Chosen in L6 by coverage.
- ONT only. HiFi is validated on simulated data until real HiFi exists.

### Points adopted from the earlier plan (2026-09-02)

The earlier plan agrees with approach B. These items from it are added to
this design:
- **FASTQ and BAM/CRAM are two input adapters behind one interface.** The
  same sample run both ways must give the same calls. A hard-clipped record
  at a junction is an error or a marked reduced-capability result, never a
  silent loss.
- **The benchmark sample key includes the platform**
  (`ont-hq_cov5x_rep1`), because the panel repeats sample names across
  platforms.
- **An evidence table per read** is written before calls are made, so a call
  can be traced.
- **Reference / non-reference / absent classes** are scored on the pangenome
  SampleA dataset. The riceTElib panel tests non-reference calls and false
  positives only.
- **Genotype thresholds** are calibrated on the panel's homozygous,
  heterozygous and 10%/20%/40% mixtures, not copied from short reads.
- **Failed HiFi 30× outputs are archived, not overwritten,** before a rerun.

## Proposed milestones

| # | Delivers | Gate |
|---|---|---|
| L1 | Test data. (a) Build the small Chr3 2 Mb long-read panel (tier 1). (b) Teach `relocate-benchmark` to read long-read manifests and register the existing full panel, ONT-HQ and HiFi 5×/15× first. (c) Repair the three HiFi 30× tasks; this does not block L2 | small panel reproducible from one script; existing short-read benchmark scores unchanged; a general SV caller + TE filter gives a baseline number |
| L2 | The data model (`Locus`, `Reference`, placement records, signature schema; see "Designing for multiple references"); `--read-type`, genome alignment presets, signature collection, clustering (steps 1–3) as tested subcommands with per-stage timing | signatures recover the simulated sites (recall of *sites*, before any filter) |
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
2. **Technologies in the first release. DECIDED 2026-11-09:** ONT (R10 /
   Q20+) and PacBio HiFi are both validated from the start. Older ONT R9 and
   PacBio CLR run through `--read-type` as "supported, less tested".
3. **Consensus engine. DECIDED 2026-11-13:** a small partial-order-alignment
   library through its Python binding. Start with abPOA, behind one function
   so spoa can be swapped in. Packages checked on conda-forge/bioconda on
   2026-11-13: `pyabpoa` 1.5.7, `pyspoa` 0.3.2. Rejected: racon/medaka
   (large dependencies, chemistry-specific models, slower per site; revisit
   only if R9 data matters to the lab) and single-best-read (ONT breakpoints
   and TSDs would carry read errors).
4. **Real data for gate 4. DECIDED 2026-11-13:** the lab's ONT R10.4.1 data
   in `/rhome/nmath020/bigdata/nanopore_data`; see "Gate 4 samples". The
   maintainer also asked for a small test set for fast iteration (tier 1)
   and named `relocate-benchmark` as the place for the large test (tier 2).
5. **Mixed short + long reads for one sample. DECIDED 2026-11-13:** out of
   scope for the first release. Each read type gives its own call set, which
   keeps the short-read path untouched and keeps gate 2 and gate 4 an
   independent cross-check. A later step can merge call sets by `site_id`.
6. **Command shape. DECIDED 2026-11-13:** one tool. `run-all`, `run-batch`
   and the pipeline sample sheet gain `--read-type`
   (`short` default, `ont`, `ont-hq`, `hifi`, `clr`); `short` is today's
   behaviour unchanged. Long-read input is a FASTQ or an existing BAM/CRAM.
   Long-read stages are their own subcommands in separate modules, so the
   validated short-read code is not edited. An option given with the wrong
   read type is rejected with a clear error. Rejected: a separate
   `long-run` command.

All six decisions are made. Next: one implementation plan per milestone,
starting with L1.

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
- the long-read mode of the panel simulator: confirmed 2026-11-13 (PBSIM3 + ccs in `simulate_data`); read realism against the real ONT data is still unmeasured.
