# LINE family assignment: confirmed tie-breaking effect

Recorded: 2026-10-01, after 14:43 PDT (America/Los_Angeles).
Project: RelocaTE3, branch `main`, following merged PR #49.
Status: extraction and analysis complete; portable fixture and five passing
tests. No production calling change, new benchmark run, or release pin.

## Result

For TE000172 in `div005_rep03_cov30x`, **eight junction reads have equally
scoring Os3328 and Os0596 TE hits**. RelocaTE2 keeps the first hit in PSL
order and chooses Os3328 for those reads. RelocaTE3 resolves the tie by
target name and chooses Os0596. This alone changes the locus's family vote
from Os3328 10:1 to Os0596 9:2, producing a false negative and false positive
under the family-aware benchmark score.

The per-read TE evidence is the same in both callers. This case does not
support recovering more alignments or changing breakpoint pairing to fix its
family assignment. It demonstrates that deterministic assignment of ambiguous
reads can create a misleading apparent family majority.

## Completion and inputs

SLURM job **29331469** ran from 13:38:07 to 13:43:18 PDT, producing
`results/error-audit/2026-10-01-line-family-trace/family_evidence.json` and
`.complete`. `logs/line-family-trace.29331469.log` records successful completion.
The extraction retained eleven read-family mapping rows for each caller,
78 original R2 PSL hits, and 77 archived R3 TE BAM alignments (15 from left
reads and 62 from right reads). Input table hashes and BAM metadata are saved.

The analysis operates on the 239-KiB extracted file and tiny local BAMs. It
does not scan the large original mappings or BAMs again.

## Causal checks

1. Apply the documented R2 gap-admission thresholds to the extracted PSL
   records: one complex alignment is rejected, leaving 77 eligible hits.
2. Convert those PSL hits using the current R3 converter and compare record
   multisets with the archived R3 SAM records. All **77 signatures agree**:
   read name, flag, TE target, position, mapping quality, CIGAR, and NM tag.
   No tool-only eligible alignment remains. This comparison excludes sequence
   and quality strings; it does not claim byte identity of whole BAMs.
3. Select each read's best R2 hit by boundary contact and match count, retaining
   the first hit on exact ties. All eleven selections reproduce archived R2
   family mappings, including agreement between its chromosome mapping and
   original chunk mappings.
4. Replay the 77 archived SAM records through production R3 `_parse_te_bam`
   and merge its left/right outputs with `_is_better`. All eleven selections
   reproduce archived R3 family mappings.
5. Apply only R3's deterministic ranking to the **same eligible R2 PSL hits**.
   All eleven family choices now reproduce R3. Eight change; three do not.

The one rejected complex hit belongs to an unchanged read, so admission
filtering does not explain the eight family switches. There is no missing
read-name mapping or inconsistent original-versus-chromosome R2 family row
in this extracted set.

## What the eleven reads actually support

| Read evidence | Reads | R2 assigned family | R3 assigned family |
|---|---:|---|---|
| Os3328 is the only family at the best primary score | 2 | Os3328 | Os3328 |
| Os0596 is the only family at the best primary score | 1 | Os0596 | Os0596 |
| Os3328 and Os0596 have the same best primary score | 8 | Os3328 | Os0596 |
| Total assigned votes | 11 | Os3328:10, Os0596:1 | Os3328:2, Os0596:9 |

For every switched read, the two leading hits agree in boundary score,
matching-base count, mismatch count, and query start/end. Both are terminal
TE alignments. Matching-base counts range from 27 to 56. R3's family-name
tie-break adds reproducibility; it supplies no biological evidence that one
of these two families is the correct origin.

Counting only reads with a single best family gives **Os3328:2, Os0596:1**.
That supports investigating whether unresolved read hits should be evaluated
using the locus's independent evidence. It is a diagnostic count, not a
validated replacement classifier. The observed reads can include mapping
ambiguity and background support, and this one locus cannot establish how a
general resolution rule behaves across all TE families or low-frequency
insertions.

R3 currently reports `TE_family_status:dominant` and confidence 9/11 ≈0.818
because these fields summarize already selected per-read labels. They do
not capture the eight unresolved per-read ties. The label is correct for the
existing vote definition, but users could reasonably interpret it as stronger
family evidence than this trace supports.

## Fixture and tests

`tests/data/line_family/TE000172.json` preserves original PSL records, mapping
rows, and all selected SAM records, with headers reduced to their relevant TE
targets. It records the full extraction's SHA-256 hash. It is a portable
characterization fixture: it reproduces historical behavior without asserting
that the historically selected R3 family is biologically correct.

`tests/line_family_evidence_test.py` covers:

- equality of eligible hit sets and reproduction of all archived assignments;
- equal primary scores, mismatches, and query spans on all eight switched reads;
- stability of deterministic selection when hit order is reversed, while
  documenting the order dependence of R2's first-on-tie selection.

Together with the two extraction-helper checks in `line_family_trace_test.py`,
**five tests pass**. The fixture requires no original benchmark files or
external aligner. Local BAM indexing uses pysam on only the extracted records.
`git diff --check` passes. These are focused validation checks, not a new full
test suite or whole-benchmark performance measurement.

## Reproduction

From the RelocaTE3 repository root:

```bash
.pixi/envs/default/bin/python scripts/analyze_line_family_trace.py \
  --evidence results/error-audit/2026-10-01-line-family-trace/family_evidence.json \
  --outdir results/error-audit/2026-10-01-line-family-comparison-v2
.pixi/envs/default/bin/python -m pytest -q \
  tests/line_family_evidence_test.py tests/line_family_trace_test.py
```

The analysis directory already exists; choose a new name to repeat it.
It contains per-read scores/selections, pooled vote counts, evidence and
source hashes, and a completion marker. The initial
`2026-10-01-line-family-comparison` directory records an earlier complete
comparison without the subsequently added hit-set and unique-family summaries;
use `-v2` for the final findings. Neither directory was overwritten.

## Scope of the conclusion and next action

The preceding audit found eight benchmark losses repeating at TE000172 across
coverage/divergence conditions. This detailed trace establishes the cause for
**one representative 30x sample**, not the per-read history of all eight runs.
The other losses and false-positive loci retain their separate explanations.
No full-benchmark precision, recall, somatic classification, or memory metric
has changed as a result of this investigation.

Restoring first-arrival tie-breaking would reproduce this R2 label but also
restore ordering dependence. It is not an evidence-based general family fix.

Recommended next implementation: retain each read's equally best TE families
alongside its deterministic primary alignment, propagate that ambiguity to
insertion evidence, and evaluate locus-level family resolution using reads
that distinguish families. Preserve a single reported primary family and
explicit uncertainty. First demonstrate the proposed evidence rule on this
fixture and controls for conflicting families and sparse somatic evidence;
then evaluate it across all three benchmark panels before adopting changed
primary labels. The known LINE ambiguity and release decision cannot be
resolved by a one-locus truth-directed relabeling.
