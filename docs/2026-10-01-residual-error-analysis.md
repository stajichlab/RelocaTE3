# Residual errors after coordinate-corrected scoring

Recorded: 2026-10-01, 13:03–13:15 PDT (America/Los_Angeles).
Project: RelocaTE3, branch `main`, following merged PR #49.
Status: completed read-only benchmark analysis; new analysis scripts and tests
only. No production caller changes, benchmark submission, release pin, or
dashboard publication. Existing uncommitted caller changes were preserved.

## Outcome

The main new finding is **family attribution, not missing junction reads**.
Sixteen of the 18 corrected divergence-panel R2-only detections have an R3
call inside the truth TSD interval or its 10-bp tolerance, with the wrong
primary family. Eight losses repeat at one LINE insertion. A bounded indexed
inspection of a representative sample found the same 11 junction-read names
in both callers at the insertion's two breakpoints.

This does not justify changing the primary family to whichever label matches
truth. The next investigation should trace TE-family assignment for those
same reads through trimming/alignment, retaining one primary family plus
ambiguity metadata. Broad minimum-support or one-sided-call exclusions would
discard many true detections.

These are **historical full-benchmark calls**, rescored on September 28 using
the explicit TSD-interval contract. They are not a full benchmark of the
uncommitted memory or breakpoint-pairing changes. Release remains on hold.

## Scope and reconciled counts

All 72 paired samples were checked: nine mPing, nine riceTElib, and 54 diverged
riceTElib samples. R3 here means BLAT/bwa aln. Truth identity, normalized-call
counts, matched-call counts, and false-positive counts reconcile with the
corrected reports. The audit uses the scored truth interval for proximity;
caller-to-caller FP correspondence still uses family and point distance ≤10 bp.
That correspondence is descriptive, not a second truth-scoring policy.

| Panel | Shared TP | R2-only TP | R3-only TP | Shared FP | R2-only FP | R3-only FP |
|---|---:|---:|---:|---:|---:|---:|
| mPing | 2818 | 2 | 8 | 0 | 0 | 0 |
| riceTElib | 2493 | 1 | 20 | 13 | 5 | 14 |
| Diverged riceTElib | 5870 | 18 | 46 | 49 | 28 | 83 |

Thus the divergence excess is **83 R3-only minus 28 R2-only = 55 FPs**, not
83 additional FPs. The 18 lost sample/event observations represent 11 distinct
event IDs. Conditions repeat loci, and zero-divergence repeats canonical
riceTElib: these are not independent biological replicates.

## The 18 R2-only divergence detections

| Explanation supported by saved tables | Observations |
|---|---:|
| Correct location, ambiguous primary family | 9 |
| Correct location, wrong dominant family | 6 |
| Correct location, wrong unique junction family | 1 |
| No normalized R3 call in truth interval ±10 bp | 2 |

### Recurring LINE: TE000172, Chr1:29402725, Os3328

R3 instead reports Os0596 at **Chr1:29402710**, exactly the left truth TSD
boundary. This is not a coordinate-scoring failure. The eight observations are:

| Divergence / replicate / coverage | R3 Os0596 : Os3328 junction votes |
|---|---:|
| 5% / 3 / 15x | 4 : 1 |
| 5% / 3 / 30x | 9 : 2 |
| 10% / 2 / 15x | 2 : 0 |
| 10% / 3 / 5x | 1 : 1 |
| 10% / 3 / 15x | 3 : 1 |
| 10% / 3 / 30x | 7 : 4 |
| 15% / 3 / 15x | 4 : 1 |
| 15% / 3 / 30x | 9 : 2 |

All eight are homozygous truth events. Six have a wrong absolute-majority
junction family, one a tie, and one only the wrong family among junction
reads. Alphabetical tie-breaking therefore cannot explain this whole group.
Supporting reads favor Os3328 in three of these observations, but are absent
in the other five; universally substituting supporting-read family is not
an established solution.

For `div005_rep03_cov30x`, archived R2 and R3 calls both have L=2, R=9,
coordinates 29402710..29402725, and TSD `GGAGAGGGAGGTGGCC`. R2 reports Os3328;
R3 reports Os0596 with votes 9:2. The R2 junction read list names 11 reads.
Indexed BAM queries of Chr1:[29402500,29402900) find exactly those same 11
names at the corresponding two breakpoints in **both** tools (no tool-only
names). Local SAM records and BAM-index hashes are preserved.

The inspected RelocaTE2 `insertion_family` implementation chooses the most
frequent junction-read family, and R3 `_te_family_evidence` also chooses the
most frequent family. R2 can append a differing supporting family; its
archived label here is a single Os3328. These observations prioritize auditing
per-read family assignment and read-name lookup upstream. They do not yet
prove which alignment/ranking/lookup operation differs: the large original
read-family and alignment tables have **not** been traced in this turn.

### Eight other family disagreements

Each has tied leading junction votes with a different deterministic R3 primary:

| Sample | Event | Truth → R3 primary |
|---|---|---|
| div002_rep02_cov15x | TE000215 | Os1415 → Os0303_complete |
| div002_rep02_cov15x | TE000382 | Os3724 → Os1100 |
| div002_rep02_cov30x | TE000093 | Os1403 → Os1100 |
| div005_rep01_cov30x | TE000484 | Os3576 → Os0029 |
| div005_rep02_cov15x | TE000030 | Os2507 → Os0596 |
| div005_rep03_cov5x | TE000028 | Os0269 → Os0043 |
| div010_rep02_cov15x | TE000277 | Os2224 → Os1218 |
| div015_rep03_cov30x | TE000212 | Os3391 → Os0596 |

Do not treat these ambiguous calls as confident family identification or
change tie ordering to favor this benchmark's truth labels.

### Two remaining losses

- `div000_rep02_cov30x`, TE000386, Os3874: repeats the canonical SINE loss.
  Prior bounded inspection found its two junction alignments in both tools;
  the complete candidate/support/filter path remains unresolved.
- `div002_rep03_cov30x`, TE000085, Os1279: R3 emits Os0029 at 22171544,
  37 bp before the truth interval's start (22171581). Its raw record extends
  through the truth anchor and has Os0029:11 versus Os1279:2 junction votes.
  This is mixed-family and breakpoint geometry, not just a missing call.

The 18 losses comprise ten homozygous, two heterozygous, and six somatic truth
observations. Both mPing R2-only detections and the one canonical riceTElib
R2-only detection are somatic. mPing's newly exposed TE000318 at
Chr1:21597869 is an R2 call at the outer tolerance boundary (21597857);
there is no nearby normalized R3 call. It has not been reconstructed here.

## Residual false positives

Of the **83 R3-only** divergence FPs:

- 18 are near a different-family truth event (including the 16 family losses);
- five are near a same-family truth already consumed by another call;
- 60 are outside every truth interval's tolerance, including 45 on chromosomes
  other than Chr1, which has all truth insertions in this benchmark.

The five same-family extras are the two previously reconstructed Helitron
duplicates at 0% and 2%, plus Os2881 at Chr1:240679 in
`div002_rep03_cov30x`. The pending pairing candidate is relevant to this
bounded subset, not the entire 55-FP net excess. Its known local sensitivity
cost must still be evaluated in the paired replay.

The 83 records comprise 36 unique-family, 30 ambiguous-family, and 17
dominant-family calls; 24 supporting-junction, 20 unknown-TSD, and 39
resolved-TSD calls. Seventeen are labeled somatic. They are heterogeneous.

Repeated off-truth-chromosome R3-only calls include:

| Locus / family | Divergence observations |
|---|---:|
| Chr4:19264325 / Os3912 | 10 |
| Chr5:5342237 / Os0915 | 6 |
| Chr6:14377462 / Os0798 | 6 |
| Chr10:10678443 / Os0868 | 5 |
| Chr3:31981544 / Os3442 | 5 |

Some have many junction reads but strongly asymmetric sides, e.g. Os3912.
Reference-copy edges, flank ambiguity, and candidate acceptance deserve
inspection; off-chromosome location establishes an error against this truth
set, not its mechanism. A Chr1-only filter would be benchmark overfitting.

Canonical riceTElib's 14 R3-only FPs include two same-family extras, two
different-family-near-truth calls, and ten outside tolerance. Seven are off
Chr1; three are labeled somatic.

## Why broad exclusions are not a no-sacrifice fix

Each rule below was evaluated against **all R3 calls**, using an exact
chromosome/start/family/TSD join to the corresponding `.raw.txt` evidence.
Every call had exactly one matching raw record. Entries count currently
matched calls and FPs that would be excluded; this is diagnostic exposure,
not a rerun or rematched estimate of final recall. Rules overlap.

| Proposed exclusion | riceTElib matched / FP excluded | Divergence matched / FP excluded |
|---|---:|---:|
| Ambiguous family | 3 / 6 | 15 / 43 |
| Fewer than three junction reads | 257 / 6 | 899 / 59 |
| One-sided junction evidence | 57 / 6 | 500 / 31 |
| Junction/supporting primary-family disagreement | 274 / 5 | 764 / 29 |

mPing has no FPs: excluding low-junction-count calls would affect 295 matched
calls; excluding one-sided calls would affect 52. Even the more selective
ambiguity exclusion sacrifices matched detections. An optional confidence
tier may be useful later, but should not silently replace the primary output.

One provenance correction matters: the earlier audit used `.all.txt` for
nearest raw evidence. A characterized Chr2:19805726 call was absent there but
present in `.raw.txt` and `.txt`. The new exposure analysis uses exact
`.raw.txt` joins, not nearest rows; it accounts for this call. The nearest-row
fields in the first audit remain descriptive and retain their offsets.

## Reproduction, outputs, validation

From the RelocaTE3 repository root:

```bash
python scripts/audit_benchmark_errors.py \
  --benchmark results/coordinate-rescore/2026-09-28 \
  --raw-benchmark ../../relocate_benchmark/relocate-benchmark \
  --outdir results/error-audit/2026-10-01-coordinate-corrected
python scripts/summarize_residual_errors.py \
  --audit results/error-audit/2026-10-01-coordinate-corrected \
  --reports results/coordinate-rescore/2026-09-28 \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --outdir results/error-audit/2026-10-01-residual-diagnosis
.pixi/envs/default/bin/python scripts/inspect_line_junctions.py \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --outdir results/error-audit/2026-10-01-line-junctions
.pixi/envs/default/bin/python -m pytest -q tests/error_audit_test.py
```

These commands have completed; choose new output directories to repeat them.
Scripts reject existing output directories. Results include reconciled
summary tables, all discordant events, all FPs, 21 R2-only event diagnoses,
159 exact-annotated R3 FPs, filter exposure, and local LINE SAM evidence.
Input table and script SHA-256 hashes are recorded. The indexed inspection
records BAM size and index hash, not a whole-BAM checksum.

Five lightweight audit tests passed, covering interval boundaries, legacy
fallback, chromosome separation, family normalization, and one-to-one FP
matching. Table assertions and the indexed read-list comparison passed;
`git diff --check` passed. CLI help was checked before script execution.
No production tests or aligners were run. Existing project Python/pysam was
used for two bounded 400-bp queries; no module-provided binary was replaced.
No BAM or large per-read alignment table was scanned. An exploratory file
listing was overly broad and truncated; subsequent checks used specific files.

## Recommended next action

Trace the **same 11 junction reads at TE000172** through RelocaTE2 and
RelocaTE3's TE-family mapping and TE-hit ranking. Use a bounded, reproducible
extraction job on SLURM for the large tables (the R3 mapping alone is 191 MiB),
save a regression fixture, and distinguish different TE-hit selection from a
read-name lookup or evidence-aggregation discrepancy before changing code.
Then test any general correction across all three panels, including somatic
events and correct family assignments. Keep the already prepared pairing
replay as a separate validation; do not assume it fixes family attribution.

No user decision is needed for this diagnosis. The concrete next work item
is the per-read family-provenance trace, not a release pin or a broad filter.
