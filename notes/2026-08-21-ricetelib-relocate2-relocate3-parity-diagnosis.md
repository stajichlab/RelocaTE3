# riceTElib RelocaTE2/RelocaTE3 parity diagnosis

**Date:** 2026-08-21 18:41 PDT

**Status:** root-cause candidates identified and measured; no code changed and no
new benchmark run submitted

**Purpose:** explain why the BLAT + `bwa aln` RelocaTE3 benchmark has RelocaTE2-like
recall on riceTElib but many more false-positive calls, and define the shortest
path to behavioral parity.

## Executive conclusion

The sensitized BLAT command is not the defect. RelocaTE3 does not preserve two
downstream RelocaTE2 contracts that make `-minScore=10 -tileSize=7` safe:

1. **RelocaTE3 dropped RelocaTE2's PSL gap-complexity filter.** RelocaTE2 rejects
   an alignment before best-hit selection when it contains too many query/target
   gaps, too many inserted bases, or at least three alignment blocks.
   RelocaTE3 converts every PSL row to SAM and never applies this filter. In
   riceTElib `cov5x_rep1`, this accounts for **101,198 of the 101,226 reads that
   RelocaTE3 classifies as TE reads but RelocaTE2 does not**. The excess includes
   **9,571 junction reads**, which can directly generate extra breakpoints.

2. **RelocaTE3 confuses “did not become a junction/middle read” with “did not hit
   a TE.”** RelocaTE2 writes every read with an admissible selected BLAT hit to
   `ContainingReads.fq`, including hits that do not satisfy a trimming branch.
   `clean_pairs_memory.py` uses that population to prevent such reads from being
   recovered as genomic mates. RelocaTE3 writes `ContainingReads.fq` only after
   a trim branch succeeds and, more importantly, step 4 uses only
   `read_repeat_name.txt` to decide whether a mate hit a TE. In the same sample,
   RelocaTE3 recovers 173,368 mates as if they were genome-only; **at least
   123,798 of them have a BLAT alignment that passes RelocaTE2's gap gate and
   therefore would be suppressed by RelocaTE2**.

These two errors explain the previously mysterious step-4 expansion. RelocaTE3
has 232,748 mapped step-4 records versus RelocaTE2's 97,496. The 123,798
wrongly admitted mates plus 9,571 extra junctions account for 133,369 records,
or **98.6% of the 135,252-record excess** before accounting for interactions
between mate pairs. This is why BLAT sensitivity strongly correlates with the
precision collapse even though RelocaTE2 uses the same BLAT parameters.

There is also an active benchmark wiring difference: RelocaTE2 always applies
its untrimmed-full-read false-junction filter, whereas the RelocaTE3 benchmark
does not pass a full-read BAM to `find-insertions`. That omission can retain
additional false calls, but it does not explain the upstream 2.4-fold BAM
expansion.

The first parity work should therefore be in steps 3 and 4, not in TSD pairing,
reference-TE filtering, or BLAT sensitivity.

## Scope and repository state

This analysis read the current versions of:

- the handoff at `notes/2026-08-20-relocate3-vs-relocate2-parity-handoff.md`;
- RelocaTE2 under `references/RelocaTE2`, especially `relocaTE2.py`,
  `relocaTE_trim.py`, `clean_pairs_memory.py`, `relocaTE_align.py`,
  `relocaTE_insertionFinder.py`, and `characterizer.pl`;
- RelocaTE3 step 3 through step 7 code under `src/RelocaTE3`;
- the benchmark configuration, caller adapters, normalizers, array driver, and
  scorer in `../../relocate_benchmark/relocate-benchmark`.

State observed during the audit:

- RelocaTE3: branch `r2-parity-work`, commit `5bd6f83` (also current `main` and
  `origin/main`), clean before this note was added.
- Benchmark: branch `main`, commit `8ad522f`; two pre-existing untracked local
  configs, `config/benchmark.bisect.toml` and
  `config/benchmark.chr3fixes.toml`, were not modified.
- The current riceTElib RelocaTE3 5x runs completed on 2026-08-20. The retained
  RelocaTE2 5x runs completed on 2026-07-28 and were normalized/rescored on
  2026-08-11. The 15x and 30x report rows are from older code and must not be
  treated as a current head-to-head result.

No computational pipeline was run on the login node. The analysis only streamed
existing FASTQ/BAM/tables and inspected code. Loading the HPCC samtools module
failed in this restricted agent environment because the module logger could not
open `/dev/log`; read-only BAM counts therefore used the already installed
RelocaTE3 pixi `samtools`/`pysam`.

## Benchmark behavior verified

The benchmark does the following for each dataset/caller/sample tuple:

1. `pipeline/submit_benchmark.sh` submits the substantial work to the `epyc`
   SLURM partition through `pipeline/run_benchmark_array.sh`.
2. The RelocaTE2 adapter builds a `bwa mem` whole-read BAM for characterization,
   then runs all RelocaTE2 steps with BLAT, `bwa aln`, `--mismatch 2`, and
   `--size 250`.
3. The RelocaTE3 parity adapter runs map/trim with BLAT, step 4 with `bwa aln`,
   step 5 with `--mismatch 2 --min-mapq 0 --require-both-junctions`, and then
   builds a minimap2 whole-read BAM for characterization.
4. Both characterized tables are parsed by the same
   `lib.calls.parse_characterized_txt` implementation.
5. `scoring/score_calls.py` greedily matches one unused call to each truth event
   within +/-10 bp and requires normalized TE-family equality. Every unused call
   is a false positive.

Current riceTElib 5x aggregate:

| Caller | Calls | Matched | False positives | Recall | Precision |
|---|---:|---:|---:|---:|---:|
| RelocaTE2 | 552 | 449 | 103 | 0.299 | 0.813 |
| RelocaTE3 BLAT/`bwa aln` | 970 | 445 | 525 | 0.297 | 0.459 |

The nearly identical matched counts confirm that the parity problem is excess
call generation, not a different truth set or scorer.

## Stage-level evidence from riceTElib `cov5x_rep1`

The following counts were taken from the retained current artifacts. “Junction”
means a FASTQ/BAM name ending in `:(start|end):[53]`; “support” means an untagged
record in the step-4 genome BAM.

| Population | RelocaTE2 | RelocaTE3 |
|---|---:|---:|
| Step-3 flanking FASTQ: junction | 50,431 | 59,545 |
| Step-3 flanking FASTQ: middle | 333,822 | 425,809 |
| `read_repeat_name` rows | 384,253 | 485,354 |
| `ContainingReads.fq` rows | 5,804,562 | 497,796 |
| `ContainingReads.fq` rows without a junction/middle tag | 5,408,944 | 0 |
| Step-4 mapped junction records | 50,413 | 59,536 |
| Step-4 mapped untagged support/mate records | 47,083 | 173,212 |
| Step-4 mapped total | 97,496 | 232,748 |
| Step-4 paired records | 37,108 | 100,098 |
| Raw non-reference rows | 308 | 327 |
| Characterized/scored calls | 188 | 327 |

RelocaTE2 and RelocaTE3 share 384,128 classified read names. RelocaTE2 has only
125 classified names absent from RelocaTE3, while RelocaTE3 has 101,226 absent
from RelocaTE2. Among the shared reads, classification differs for only 338
reads. Thus the large count difference is a specific R3-only population, not a
general inability to match read names between implementations.

## Primary issue 1: the RelocaTE2 PSL gate was not ported

### RelocaTE2 behavior

`references/RelocaTE2/scripts/relocaTE_trim.py:63-66` rejects a PSL alignment
when any of the following is true:

```text
qNumInsert > 1
qBaseInsert > 3
tNumInsert > 1
tBaseInsert > 3
blockCount >= 3
```

Only surviving alignments participate in the per-read best-hit comparison at
lines 87-116. The ranking is boundary score first and PSL `matches` second.

### RelocaTE3 behavior

`src/RelocaTE3/aligners.py:521-578` converts every PSL record to SAM. It retains
gap information in the CIGAR but discards the original PSL admission decision.
`src/RelocaTE3/librelocate.py:384-474` then selects a best record without an
equivalent gap/block filter.

This is not merely a theoretical difference. For each R3-only classified read,
I reconstructed the five RelocaTE2 conditions from the selected record's CIGAR,
which is generated directly from the PSL block/query-gap/target-gap arrays.

| R3-only selected alignment property | Reads |
|---|---:|
| R3-only classified reads | 101,226 |
| Violates at least one R2 gap rule | **101,198** |
| `blockCount >= 3` | 56,536 |
| `qBaseInsert > 3` | 57,821 |
| `qNumInsert > 1` | 25,542 |
| `tBaseInsert > 3` | 80,052 |
| `tNumInsert > 1` | 39,974 |

The rule counts overlap. Of the R3-only classified reads, 9,571 are junctions
and 91,655 are middle reads. The selected-alignment violation rate is 99.97%.

This finding explains why reverting to BLAT defaults appeared to improve
riceTElib precision: fewer complex low-score hits reached R3's missing gate. The
correct parity fix is to restore the gate, not to use different BLAT parameters.

### Required fix

Apply the five RelocaTE2 PSL conditions before PSL-to-SAM conversion and before
best-hit selection. Preserve the distinction between PSL `matches`,
`misMatches`, and gap bases; do not infer a subtly different score later from
`M` and `NM` if exact BLAT parity is the goal.

For the BLAT backend, equal boundary/match hits should initially preserve
RelocaTE2's first-surviving-PSL-hit behavior. RelocaTE3 currently adds a
lexicographic TE/coordinate total-order tie-break after coordinate-sorting the
converted BAM. That makes runs deterministic, but it is not RelocaTE2's logic.
A stable original-PSL ordinal would allow both determinism and faithful
first-hit behavior.

## Primary issue 2: step 4 uses the wrong definition of a TE-hit mate

### The three populations must remain distinct

RelocaTE2 distinguishes:

1. a read with any selected, gap-admissible TE alignment;
2. a read that becomes a 5'/3' junction or middle classification;
3. a read with no admitted TE alignment, which may be a genomic mate.

The distinction is explicit in `relocaTE_trim.py`:

- `read_repeat_name.txt` is written only in a successful junction/middle branch
  (`:407`, `:435`, `:449`);
- `ContainingReads.fq` is written for every read in the post-filter best-hit
  coordinate table, even if no branch succeeded (`:454-455`).

`clean_pairs_memory.py:163-221` checks the mate against `ContainingReads.fq`,
not against `read_repeat_name.txt`. A mate with an admissible TE hit is not a
genomic support/anchor read merely because it failed the trim thresholds.

RelocaTE3 collapses populations 1 and 2:

- `src/RelocaTE3/librelocate.py:217-224` continues before writing
  `ContainingReads.fq` when `_trim_record` returns `None`, so its file no longer
  means “any selected TE hit.”
- `src/RelocaTE3/genome_align.py:109-123` and `:153-195` construct `te_ends`
  only from `read_repeat_name.txt` and assume any absent mate is genomic.

The comments and unit tests currently encode the same false equivalence. For
example, `tests/genome_align_test.py:73-87` says absence from `read_repeat`
means the mate “did NOT match the TE.” That is not what the R2 artifact means.

### Direct measurement

The R3 classified FASTQs contain 329,361 read-pair bases:

| R3 pair state | Count |
|---|---:|
| Both ends classified | 155,993 |
| Only one end classified | 173,368 |

For the 173,368 absent mates, the current R3 logic creates:

- 50,061 junction + presumed-genomic-mate pair candidates;
- 123,307 recovered support mates for middle reads.

The resulting BAM contains almost exactly those records: 50,055 paired junction
records, 50,043 paired untagged mates, and 123,169 other untagged records.

Streaming the TE-alignment BAMs for just those absent mate names showed:

| Source classified end | Recovered mate candidates | Any BLAT hit | At least one R2-gap-admissible hit |
|---|---:|---:|---:|
| Junction | 50,061 | 30,559 | **29,042** |
| Middle | 123,307 | 103,373 | **94,756** |
| Total | 173,368 | 133,932 | **123,798** |

Thus 123,798 mates are admitted by R3 that R2's `ContainingReads` check would
exclude. This alone is nearly the entire unexplained step-4 count difference.

### Why this generates two-sided false calls

The effect is larger than simply adding support counts:

- A junction plus a falsely presumed genomic mate is sent through `bwa sampe`.
  If `sampe` marks it proper, `InsertionFinder._passes_quality` accepts it under
  the paired branch and does not require the junction to have `XT:A:U`. The same
  junction mapped single-end in R2 would have to pass the unpaired uniqueness
  gate.
- The untagged mate is itself added to cluster support.
- Every admitted support record extends the 1,000 bp streaming cluster. Dense
  TE-derived placements can therefore join breakpoint populations that would
  otherwise be separate.
- The additional 9,571 gap-complex junctions from issue 1 supply real left/right
  breakpoint tags. The falsely admitted mates help place and cluster them.

This combination explains why the final excess is two-sided even though the
largest numerical inflation is in untagged mates.

### Required fix

Step 3 needs an explicit, reproducible artifact containing the names of **all
post-PSL-filter selected TE-hit reads**, regardless of trim outcome. It may be a
RelocaTE2-compatible `ContainingReads.fq` or a compact names table, but its
contract must be distinct from `read_repeat_name.txt`.

Step 4 must use this all-hit membership for the exact R2 truth table:

| Classified read | Mate state | RelocaTE2 action |
|---|---|---|
| Junction | Mate is a junction | align the two trimmed flanks together with `sampe` |
| Junction | Mate is middle or has another admitted TE hit | align the junction flank single-end; do not use the mate |
| Junction | Mate has no admitted TE hit | align flank + original genomic mate with `sampe` |
| Middle | Mate has no admitted TE hit | discard middle read; align original mate as support |
| Middle | Mate has an admitted TE hit | discard both for step 4 |

The all-hit table should be a declared step-3 output and step-4 input, which also
fits the intended future Nextflow process boundary.

## Secondary active issue: the benchmark bypasses false-junction filtering

RelocaTE2 always constructs
`<target>.repeat.fullreads.bwa.sorted.bam` in `relocaTE_align.py:337-350`.
`relocaTE_insertionFinder.py:1804-1817` always loads it, associates untrimmed
reads with junction reads, and `write_output:212-221` removes a candidate when
at least 30% of the junction reads on both sides map through the putative
breakpoint.

RelocaTE3 now has `_fullread_false_junction` and exposes
`find-insertions --fullreads-bam`. However, the benchmark adapter currently:

1. calls `relocaTE3 align-genome`, whose direct CLI path writes only the flank
   BAM;
2. calls `relocaTE3 find-insertions` without `--fullreads-bam`;
3. only afterward builds `original_reads.sorted.bam` with minimap2 for
   characterization.

Therefore `_fullread_false_junction(None, ins)` returns false for every candidate
and the R3 benchmark applies no false-junction filter. The library
`pipeline.run_sample` path does wire a full-read BAM, so unit and acceptance
tests of that path do not prove that the benchmarked CLI sequence uses the
filter.

For exact parity, do not simply pass the later minimap2 whole-read BAM and call
the issue closed. RelocaTE2 extracts the corresponding original read pairs and
aligns them with `bwa aln`/`sampe`; RelocaTE3's current
`collect_junction_fullreads` extracts junction reads only and aligns them
single-end. The parity implementation should reproduce the R2 population and
pairing, expose it as a step-4 output, and require it as a step-5 input for the
BLAT/`bwa aln` parity mode.

This issue has not been end-to-end measured on the current 327 R3 calls, so its
incremental call reduction is unknown.

## Additional code-level parity differences

These are real differences but are lower-priority explanations for the large
riceTElib false-positive excess.

### Both mates are junctions: R3 maps them separately

`clean_pairs_memory.py:120-127` writes two non-middle trimmed mates to matched
FASTQs and RelocaTE2 maps them with `sampe`. In
`genome_align.build_flank_pairs`, the other end's presence in `te_ends` makes
`mate_name` remain `None`, so both junction flanks go to the single-end file.
There are 888 such junction ends (444 read pairs) in `cov5x_rep1`. This is not
the dominant count difference but can alter repetitive short-flank placement.

### TE-family ties do not follow R2

R3 intentionally chooses lexicographically among equal read-to-TE hits and
among tied cluster family votes. R2 keeps the first equal boundary/match PSL hit
and its Python 2 cluster vote sorts by count only. On the 384,128 commonly
classified reads in `cov5x_rep1`, **8,073 (2.1%) have different selected TE
families**. Family differences turn otherwise positional matches into false
positives under this benchmark's family-aware scorer, although they do not
explain the excess call count.

Parity mode should reproduce R2's read-level first-hit behavior. A separately
documented deterministic policy can remain for future non-parity modes.

### One-base genome-end convention differs

RelocaTE2 assigns `start = reference_start + 1` and
`end = reference_end + 1` in `relocaTE_insertionFinder.py:1351-1359`.
RelocaTE3 `_stream_clusters` uses `gstart = reference_start + 1` but
`gend = reference_end` at `insertions.py:826-829`. This changes left
breakpoints/TSD geometry by one base. The scorer's +/-10 bp window masks most
of the impact, but exact output parity requires choosing and testing the R2
convention.

### Orientation ties differ

RelocaTE2 reports `+` only when forward votes are strictly greater than reverse
votes; ties become `-` (`relocaTE_insertionFinder.py:202`). RelocaTE3 uses `+`
on a tie (`insertions.py:1163-1164`). Strand is not used by the current scorer,
but this is a direct output-parity difference.

### One-sided `supporting_junction` behavior remains unported

The benchmark asks R3 to require both junctions. That matches most of R2's final
characterizer gate, but R2 also retains a small `supporting_junction` class.
R3 cannot safely restore that exception until step-4 mate admission and support
counts match. The mate-membership bug identified above is likely a major reason
the existing `SR`/`SL` counts were dramatically inflated. Fix steps 3/4 first,
then remeasure and port the narrow exception.

### Insert size and characterization BAM are uncontrolled

- The benchmark passes `--size 250` to R2 but nothing to R3, which uses 500.
  This mainly affects support-only intervals and is unlikely to create the
  scored two-sided excess, but the values must be equal for a parity claim.
- R2 characterization uses a `bwa mem` BAM; R3 uses minimap2. Detection
  precision/recall is decided before this BAM, but status comparisons are
  confounded. Build one shared characterization BAM or use the same command for
  both.

### R3-only TSD policy remains in parity mode

`MAX_PLAUSIBLE_TSD = 20` changes long inferred TSD strings to `UNK`. This may be
a sensible future policy and is not a source of the current call-count excess,
but it is not R2 logic. Keep improvements out of the parity baseline until R2
equivalence is established.

## Recommended implementation and validation order

### 1. Restore the BLAT PSL contract

- Filter PSL rows with the exact five R2 conditions before conversion.
- Rank surviving BLAT records as R2 does: boundary, then PSL matches, then first
  surviving PSL record on an exact tie.
- Add focused tests for every threshold (`1` versus `2` gap openings, `3` versus
  `4` inserted bases, and `2` versus `3` blocks).
- Add a multi-family test where an invalid complex hit currently outranks or
  coexists with a valid simple hit.

Acceptance checkpoint on existing `cov5x_rep1` input: the 101,226 R3-only
classified population should collapse to approximately zero, subject to the
125 R2-only and tie-order differences already measured.

### 2. Port the complete `clean_pairs_memory.py` state machine

- Emit all-hit membership separately from `read_repeat_name`.
- Use it when deciding whether an original mate is genomic.
- Pair two junction flanks with each other.
- Preserve the junction/middle/weak-hit/no-hit cases in the table above.
- Add unit tests where the other mate has a TE hit that fails trimming; current
  tests lack this essential case.

Acceptance checkpoint: for riceTElib `cov5x_rep1`, step-4 mapped records should
fall from 232,748 to the same order as R2's 97,496, with comparable
paired/unpaired and junction/support partitions.

### 3. Make the R2 full-read control an explicit step output/input

- Generate the original-read pair subset with R2-equivalent naming and
  `bwa aln`/`sampe` behavior.
- Make the benchmark pass it to `find-insertions --fullreads-bam`.
- Log the number of candidates removed by the 30% rule.

### 4. Re-run one sample through SLURM before touching step-5 algorithms

Run only `ricetelib`, `relocate3-blat-bwaaln`, `cov5x_rep1` through the benchmark
submission script after removing/moving the old run directory. Record, in order:

1. all post-gap-filter TE-hit names;
2. junction/middle classifications;
3. paired-junction, junction+genomic-mate, single-junction, and support-mate
   step-4 inputs;
4. mapped/proper/unpaired records admitted by step 5;
5. clusters, two-sided candidates, full-read removals, raw rows, and
   characterized rows.

Expected result: the BAM and two-sided call counts should move sharply toward
R2 without reverting BLAT sensitivity. If they do not, perform a per-locus read
diff using R2's `Chr1.repeat.reads.list` and a new equivalent R3 read-list
artifact before changing breakpoint pairing again.

### 5. Expand only after the single-sample stage counts agree

1. all three riceTElib 5x replicates;
2. mPing 5x/15x/30x to confirm its recovered sensitivity is retained;
3. riceTElib 15x/30x;
4. only then the divergence panel.

Use fresh `.run_complete` timestamps and per-sample report timestamps. Do not
mix old 15x/30x reports with new 5x results.

### 6. Close remaining exact-output differences

After the false-positive gap closes, align insert size, full-read and genotype
aligners, endpoint convention, TE-family/orientation ties, and the
`supporting_junction` exception. At that point parity should be asserted at each
step, not only through final precision/recall.

## Final assessment

The riceTElib failure is best understood as a missing pair of filters around
sensitive TE search:

```text
sensitized BLAT
    -> R3 retains complex PSL hits that R2 rejects
    -> R3 classifies ~101k extra reads, including ~9.6k junctions
    -> R3 uses only classified reads to decide whether a mate hit a TE
    -> ~124k TE-hit mates are incorrectly treated as genome-only
    -> sampe anchoring + extra support place/cluster those reads
    -> excess two-sided calls
```

This chain fits all central observations from the handoff: mPing is nearly at
parity, riceTElib fails after sensitized BLAT, R3 has fewer/all-different TE-hit
semantics yet many more genome placements, and the final excess survives the
both-junction gate. It also identifies specific R2 logic to port rather than a
new heuristic intended to improve precision.
