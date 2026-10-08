# Recurring false positives: intact original-read evidence

Completed follow-up: [intact-read results](2026-10-02-intact-read-results.md).
Array 29354331 removed 44 FPs but lost five true mPing detections; the global
intact-read candidate was not adopted. Submission notes below are historical.

Recorded: 2026-10-02 15:34 PDT, America/Los_Angeles.
Project: RelocaTE3, `main`, base `942a5ea`, after merged PR #49.
Status: source/read diagnosis and 19 focused tests complete; frozen 72-sample
candidate comparison prepared. Production caller/filter code is unchanged.

## Updated residual audit

The completed stabilization replay contains 27 canonical and 129 divergence
R3 FPs, versus R2's 18 and 77. One-to-one caller FP correspondence (same
normalized family within 10 bp) gives:

| Panel | Shared FP | R2-only FP | R3-only FP | Net R3 excess |
|---|---:|---:|---:|---:|
| mPing | 0 | 0 | 0 | 0 |
| riceTElib | 13 | 5 | 14 | 9 |
| Diverged riceTElib | 49 | 28 | 80 | 52 |

This correspondence is descriptive and does not change truth scoring.
Of the 80 divergence R3-only FPs, 60 are outside truth windows, 15 are near
different-family truth, and five are same-family extras near truth. Canonical
counts are 10, two, and two respectively. Each R3 FP joins exactly one raw
call by chromosome/start/family/TSD, with complete support/ambiguity metadata.

The most repeated divergence-only locus is Os3912 at Chr4:19264325 (10
observations). Other repeated loci include Os0915 at Chr5:5342237 (six),
Os0798 at Chr6:14377462 (six), and Os0868 at Chr10:10678443 (five).
These counts are repeated conditions, not independent biological loci.

Reproducible audit:

```bash
.pixi/envs/default/bin/python scripts/rank_replay_false_positives.py \
  --replay results/family-replay/2026-10-02-stabilization \
  --output results/error-audit/2026-10-02-stabilized
```

Already run; outputs are protected against overwrite. Tables and input hashes
are in that directory. An initial duplicate `sample` keyword in the reporting
script failed before creating outputs; corrected and rerun successfully.

## Concrete example: Os3912 reference SINE edge

In canonical `cov30x_rep1`, R3 reports Os3912 at
Chr4:19264325..19264335, TSD `CTTAGTCGATA`, with 19 left and one right junction.
RepeatMasker annotates a reference SINE at 19264336..19264629, directly adjacent
to the left-junction endpoints. This proximity alone is not a sufficient
reason to reject a call: genuine insertions can occur beside reference TEs.

Bounded indexed queries of Chr4:[19263825,19264825) preserve 20 junction SAM
records and 20 matching full-read SAM records for each caller. The key right
read is `cov30x_rep1:baseline:h1:s547481152:Chr4-402230/1:start:5`:

| Evidence | Alignment |
|---|---|
| Trimmed R3 flank | 19264325..19264336, `12M` |
| Original R3 read, same mate | 19264493..19264642, `150M` |

The original read aligns fully at a different nearby reference position. It
does not span the candidate breakpoint, so R3 does not count it as evidence
against the junction. The 19 left original reads also align intact. Under
R2's intact-read rule the candidate has 19/19 and 1/1 opposing full-read
observations, exceeding the existing 30% threshold on both sides. The same
right read end aligns intact in the archived R2 full-read BAM.

## Source difference, precisely

In `../references/RelocaTE2/scripts/relocaTE_insertionFinder.py`:

- `read_junction_reads_align` (around lines 1645–1702) marks a mapped read
  intact when CIGAR `M + I >= query length - 10`. It does not require a
  proper-pair flag or a MAPQ threshold, despite nearby comments saying
  "mapped properly". Insertions count as aligned query bases; deletions do not.
- `find_insertion_cluster_bam` (around lines 1382–1428) treats that intact
  flag as opposing junction evidence without requiring breakpoint overlap.
  Its later fallback checks matching alignment endpoints and extension.
- `junction_full_reads` plus `write_output` reject when at least 30% of
  junction observations on **each** side have opposing full-read evidence.

R3 `_fullread_false_junction` currently stores local alignment spans and tests
whether those spans cross the candidate breakpoint with a margin. It does
not retain the intact flag. The representative displaced right read therefore
prevents rejection even though both tools have its full-read alignment.

This establishes a missing evidence rule at this fixture, not that all 52
net extra FPs have this cause, nor that applying the rule is harmless to real
insertions. Intact alignments to other copies can also occur for real events.

## Portable fixture and controls

`tests/data/intact_fullreads/Os3912.json` (about 86 KB) contains actual SAM
records, original headers, bounded query coordinates, input paths/sizes,
index hashes, and extractor hash. No invented or shifted reads are used.

```bash
.pixi/envs/default/bin/python scripts/extract_intact_fixture.py \
  --replay results/family-replay/2026-10-02-stabilization \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --output tests/data/intact_fullreads/Os3912.json
.pixi/envs/default/bin/python -m pytest -q \
  tests/intact_fullreads_test.py tests/false_junction_parity_test.py
```

Extraction and **19 tests passed**. Controls cover the ten-base cutoff,
insertions versus deletions, extended CIGAR operations, unmapped records,
mate identity, one-sided and empty candidates, existing spanning-read
behavior, and the actual displaced full-read fixture. Lint and shell syntax
checks pass. The extractor's name-set initialization was made explicit after
lint; the stored fixture hash identifies the extractor version actually run.

## Frozen diagnostic comparison

Prepared `results/intact-fullreads-replay/2026-10-02/` from the successful
stabilization snapshot. No production filter was changed. The experimental
runner adds the intact-read veto by replacing the filter only in its own
compute process. The original spanning-read veto remains: this is an additive
candidate, not a full reimplementation of every R2 filtering detail.

Each of 72 tasks:

1. Verifies frozen code and input metadata/hashes.
2. Scans saved junction/full-read BAMs on SLURM, retaining only matching
   read/mate keys and intact flags, not full sequences or alignment objects.
   This finds intact alignments outside the local breakpoint window too.
   Per-key updates follow input order, including mapped/unmapped replacement;
   secondary/supplementary counts are reported for interpretation.
3. Replays baseline calling/characterization and requires exact reproduction
   of the stabilized normalized calls.
4. Replays with the added intact-read veto, keeping all other settings fixed.
5. Scores both outputs with `tsd-interval` and 10-bp tolerance, saving every
   truth match, FP, summary, and additional veto with its read names/counts.

This tests downstream arbitration effects, not just removal of final rows.
Both gains and losses must be reviewed across TE groups, coverage, divergence,
and somatic fractions before any production change. It does not repeat
alignments or establish full-pipeline runtime/memory performance.

Preparation already completed:

```bash
.pixi/envs/default/bin/python scripts/replay_intact_fullreads.py prepare \
  --baseline results/family-replay/2026-10-02-stabilization \
  --output results/intact-fullreads-replay/2026-10-02
```

The epyc request is 32 GB, one CPU, eight hours, four concurrent tasks;
the cached partition maximum is 30 days. The job captures SLURM_SUBMIT_DIR
and uses relative paths, with no --chdir directive. Submission/status is
recorded below. No cleanup of previous benchmark output is needed.

Submission with a bounded 15-second `sbatch --parsable` call returned exit 124
without a numeric job ID. It is **unconfirmed**, not reported as running.
From the RelocaTE3 root, check for an existing `r3-intact-read` job first:

```bash
squeue -u "$USER" -o '%.18i %.45j %.10T'
```

If none exists, submit the already prepared snapshot:

```bash
REPLAY=results/intact-fullreads-replay/2026-10-02
sbatch --job-name=r3-intact-read --array=0-71%4 \
  "$REPLAY/replay_intact_fullreads.slurm" "$REPLAY"
```

Next action: inspect the completed paired scores and additional-veto records;
retain the current production behavior until sensitivity costs are known.
No release pin, commit, push, or broad filtering change was made.
