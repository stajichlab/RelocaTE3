# Completed family replay and stabilization decision

Follow-up: [the successful stabilization results](2026-10-02-stabilization-results.md)
verify array 29349059, restoration of all 30 detections, and three family-related
gains without losses. The failed submission notes below remain historical.

Recorded: 2026-10-02, America/Los_Angeles (analysis session began 01:12 PDT).
Project: RelocaTE3, `main`, base `942a5ea`, after merged PR #49.
Status: 72/72 replay tasks complete; release gate remains open.

## Verified inputs and comparison

SLURM array 29333789 produced all 72 completion markers and comparison files.
All frozen source hashes verify. The comparison covers 9 mPing, 9 riceTElib,
and 54 diverged riceTElib samples. Both replay arms use identical calling code,
including the pending September breakpoint adjustment. Only tied-family
evidence differs. R2 uses its archived calls. All results below use the same
10-bp `tsd-interval` scoring policy, with genotype and TSD assessed separately.

The historical R3 column is the original full benchmark rescored under that
same policy. It is not the replay baseline: baseline calls match history in
only 52/72 samples. Family-only changes leave coordinates, TSD, strand, and
genotype unchanged in all 72 samples.

| Panel | Version | TP | FP | FN | Precision | Recall |
|---|---|---:|---:|---:|---:|---:|
| mPing | R2 | 2820 | 0 | 1680 | 100% | 62.667% |
| mPing | Historical R3 | 2826 | 0 | 1674 | 100% | 62.800% |
| mPing | Replay original / tied | 2826 | 0 | 1674 | 100% | 62.800% |
| riceTElib | R2 | 2494 | 18 | 2006 | 99.283% | 55.422% |
| riceTElib | Historical R3 | 2513 | 27 | 1987 | 98.937% | 55.844% |
| riceTElib | Replay original / tied | 2499 | 23 | 2001 | 99.088% | 55.533% |
| Diverged riceTElib | R2 | 5888 | 77 | 21112 | 98.709% | 21.807% |
| Diverged riceTElib | Historical R3 | 5916 | 132 | 21084 | 97.817% | 21.911% |
| Diverged riceTElib | Replay original | 5900 | 125 | 21100 | 97.925% | 21.852% |
| Diverged riceTElib | Replay tied | 5902 | 123 | 21098 | 97.959% | 21.859% |

The family candidate gains three detections and loses one. It does not justify
a no-sacrifice release claim. Versus R2, the tied candidate still misses
2 mPing, 2 canonical, and 16 diverged detections that R2 finds, while gaining
8, 7, and 30 respectively. Net counts hide those losses.

## Regression mechanisms and next test

Historical R3 → replay original loses 14 canonical and 16 diverged detections,
with no gains. This removes only 4 and 7 FPs respectively. All 30 losses are
Helitron (13) or Tc1/Mariner (17); four are somatic insertion observations.
The pending exclusive-endpoint breakpoint change is the leading cause, but
the completed replay did not isolate that function from all other source
changes. A frozen follow-up restores only `_pair_breakpoints` from HEAD and
requires original-evidence calls to reproduce history before completing.
Do not yet attribute all 30 losses causally to pairing or alter production
pairing solely from the pooled counts.

R2 source inspection confirms exclusive left coordinates during pairing and
the subtraction during TSD inference (`TSD_from_read_depth` and
`TSD_len_calculate` in `../references/RelocaTE2/scripts/relocaTE_insertionFinder.py`).
Copying the former coordinate convention into R3 without validating the whole
candidate lifecycle was insufficient. Two local duplicate fixtures did not
establish full-panel noninferiority.

Family-only gains are the Os3328 LINE event TE000172 in div005_rep03_cov30x,
div010_rep03_cov30x, and div015_rep03_cov30x. The loss is somatic hAT event
TE000496, div002_rep02_cov30x, at Chr1:960139 (cell fraction 0.4).
Its raw candidate changes Os1524 → Os3293 despite `unresolved_ties`:
old selected votes are 9:7; deduplicated selected reads are 5:7.
Independent family evidence is Os1524=4, Os3293=7, with one ambiguous read.
The fallback was recomputed after deduplication, silently changing the label
even though geometry prevents resolution. This is an implementation regression,
not evidence that the consensus quorum should be tuned to simulated truth.

The small production fix preserves the pre-deduplication selected primary
when ties remain unresolved. Distinguishing evidence stays deduplicated;
repeated names cannot satisfy the quorum. Compatible consensus may still
override the fallback. A regression test reproduces the 9:7 versus 5:7
multiplicity pattern without benchmark-specific thresholds or family rules.
The full-panel effect of this fix remains pending the follow-up replay.

## Somatic and difficult-condition performance

| Panel | Truth somatic events | R2 detected | Tied R3 detected | R2 correctly labeled somatic | Tied R3 correctly labeled somatic |
|---|---:|---:|---:|---:|---:|
| mPing | 2700 | 1191 | 1194 | 249 | 250 |
| riceTElib | 2700 | 1033 | 1037 | 181 | 180 |
| Diverged riceTElib | 16200 | 2135 | 2142 | 526 | 525 |

Detection is slightly higher overall, but genotype/status is a separate
limitation. Both callers miss many low-frequency events and frequently label
detected somatic insertions as another genotype. No long-read/genotype changes
are warranted as part of this small stabilization patch.

At 15x mPing coverage R3 detects 987 versus R2's 989, despite the pooled mPing
advantage. At 15% library divergence R3 detects 30 versus R2's 31 and emits
14 versus 4 FPs. At 20% divergence both detect only 10/4500, with R3 6 versus
R2 2 FPs. Low pooled divergence recall therefore cannot be described as good
absolute sensitivity merely because the tools are close to each other.

All coverage, divergence, TE-group, biological-class, and somatic-fraction
tables are in `results/family-replay-analysis/2026-10-02/summary.tsv`.
Truth-stratum rows intentionally omit FP/precision: unmatched calls cannot be
assigned a truth genotype/fraction. `changed_events.tsv` records every gain
and loss; `paired.tsv` preserves sample-level paired counts. Observations
across coverage/replicates include repeated simulated loci, not independent
discoveries. Runtime in the replay covers downstream calling only, not a
new whole-pipeline runtime or memory validation.

## Reproduction

From the RelocaTE3 repository root (existing outputs are protected):

```bash
.pixi/envs/default/bin/python scripts/summarize_family_replay.py \
  --replay results/family-replay/2026-10-01 \
  --historical results/coordinate-rescore/2026-09-28 \
  --output results/family-replay-analysis/2026-10-02
```

The summarizer checks completion, frozen source hashes, task identities, unique
event IDs, consistent truth sets, and TP agreement between event and precision
tables. Its provenance records hashes of all scored tables and the script.
It reads small result tables only; BAM scans and caller execution stay on SLURM.

## Release decision

Do not pin or promise R2 noninferiority yet. Keep the validated memory work,
fix the unresolved fallback, and validate restoration of the lost sensitivity
before introducing another calling heuristic. Precision, feature coverage,
and full-pipeline resource/testing gates remain. Future direction is recorded
in `2026-10-02-release-first-direction.md`; no Rust, Nextflow, long-read, or
pangenome implementation is included in this change.

## Follow-up prepared at 01:22 PDT

The fallback fix passes **84 focused tests** spanning family resolution,
determinism/metadata, trimming, insertion/TSD behavior, characterization,
duplicate fixtures, and memory regressions. Lint, shell syntax, and
`git diff --check` pass. Every analysis stratum reconciles with panel totals.
This is not the full external-tool integration suite.

Prepared `results/family-replay/2026-10-02-stabilization/` with 72 tasks:

```bash
.pixi/envs/default/bin/python scripts/replay_family_resolution.py prepare \
  --benchmark ../../relocate_benchmark/relocate-benchmark \
  --output results/family-replay/2026-10-02-stabilization \
  --cached-evidence results/family-replay/2026-10-01 \
  --historical-pairing
```

The snapshot includes the fallback fix and restores historical pairing only
inside the frozen validation source. Production pairing remains untouched
pending results. AST comparison confirms that its pairing function exactly
matches HEAD and the rest of the insertion module matches current source.
All frozen file hashes verify. The replay refuses completion if original
evidence does not reproduce historical normalized calls. It then tests tied
evidence against that reproduced baseline, using the same restored pairing.

Enriched mappings come from the completed replay; TE alignments are not
rescanned. Preparation verifies unchanged TE parser source, prior completion,
sample order, and original input size/mtime. Each task verifies frozen source
hashes, input sizes/mtimes, and small truth/call hashes. Cached enriched mapping
files are checked by size/mtime, not newly SHA-hashed in this preparation.
Calling and characterization execute on SLURM. The existing epyc eight-hour
request is below the cached partition limit of 30 days.

Submission was attempted once with `timeout -k 2s 15s sbatch --parsable ...`;
it returned exit 124 with no job ID. Submission is **unconfirmed**. Do not infer
that the follow-up is running. From the RelocaTE3 root, check before submitting
again to avoid a duplicate:

```bash
squeue -u "$USER" -o '%.18i %.45j %.10T'
```

If no `r3-stabilization` job exists and no new task outputs have appeared:

```bash
sbatch --job-name=r3-stabilization --array=0-71%4 \
  results/family-replay/2026-10-02-stabilization/replay_family_resolution.slurm \
  results/family-replay/2026-10-02-stabilization
```

After completion, verify historical reproduction and measure gains/losses
again before deciding whether to restore production pairing. No commit,
push, merge, release tag, benchmark deletion, or dashboard publication occurred.

## Failed submission diagnosis, 2026-10-02 13:50 PDT

The maintainer's array 29340746 finished, but all 72 logs show the same
pre-execution failure: the supplied directory ends in `2026-10-02-stabilizatio`
(missing the final `n`). Python could not open `replay_family_resolution.py`.
No task directory, comparison, or completion marker was created. This run
provides no new accuracy results and does not test the candidate code.

Verified all 72 error logs, all frozen source hashes, and every recorded input
size/mtime. The intended `2026-10-02-stabilization` directory is intact and
requires no cleanup or re-preparation. Rechecked SLURM syntax and the epyc
partition's 30-day maximum against the eight-hour request.

Attempted submission with the complete directory argument once from this
session. The bounded `sbatch` call again returned exit 124 without a job ID;
submission remains unconfirmed. Check for `r3-stabilization` in the queue
before repeating the corrected submission command above. No production code
was changed during this diagnosis.
