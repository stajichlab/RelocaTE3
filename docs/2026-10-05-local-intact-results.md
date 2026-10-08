# Local intact-read filter: precision improvement without observed detection loss

Checkpoint follow-up, 2026-10-05 14:27 PDT: SLURM job **29419250** passed
the full suite: **331 tests, zero failures, zero errors, zero skips**.
All 89 recorded source/test-control hashes still match the current files.
The full-suite gate described as pending below is now satisfied; release-wide
performance and utility limitations remain. See
[the checkpoint handoff](2026-10-05-checkpoint-handoff.md).

Recorded: 2026-10-05, after 13:54 PDT, America/Los_Angeles.
Project: RelocaTE3, `main`, base `942a5ea`, following merged PR #49.
Status: 72-sample local experiment complete; local rule integrated into
production with 127 focused tests passing. Full checkpoint suite pending.
No commits, push, merge, release tag, or dashboard publication performed.

## Evidence supporting adoption

Array 29380561 completed all 72 tasks on October 3 at 16:44 PDT. All tasks
reproduced the stabilized baseline; frozen source hashes verify. The paired
event audit finds **no lost or gained truth detections** relative to baseline.
Normalized-call multiset comparison finds exactly **28 removals and no
additions**: three canonical and 25 diverged riceTElib FPs. Retained calls
preserve coordinates, family, TSD, strand, and genotype.

All five mPing regression targets from the global experiment remain detected:
TE000316/cov30x_rep2, TE000125/cov5x_rep1, and TE000157, TE000434,
TE000452/cov5x_rep2. Thus the precision improvement did not erase these
R3-only true detections. The local rule is adopted; the global rule is not.

| Panel | Caller | TP | FP | FN | Precision | Recall |
|---|---|---:|---:|---:|---:|---:|
| mPing | R2 | 2820 | 0 | 1680 | 100% | 62.667% |
| mPing | Local R3 candidate | 2826 | 0 | 1674 | 100% | 62.800% |
| riceTElib | R2 | 2494 | 18 | 2006 | 99.283% | 55.422% |
| riceTElib | Local R3 candidate | 2513 | 24 | 1987 | 99.054% | 55.844% |
| Diverged riceTElib | R2 | 5888 | 77 | 21112 | 98.709% | 21.807% |
| Diverged riceTElib | Local R3 candidate | 5919 | 104 | 21081 | 98.273% | 21.922% |

The stabilized R3 baseline had 27 and 129 FPs in the latter two panels.
The remaining net FP gaps against R2 are six and 27, respectively. R3 retains
pooled recall advantages but is not uniformly superior: the previously
documented R2-only events and low-fraction somatic shortfalls remain.
No truth detections changed within any biological class, fraction, TE group,
or divergence stratum. Canonical and zero-divergence observations overlap;
these totals are benchmark observations, not independent biological loci.

## Production integration

`src/RelocaTE3/insertions.py` now collects local intact keys during the existing
bounded full-read fetch. A new pure helper recognizes near-complete query
alignment (M/I/=/X covering at least query length minus ten). The same-mate
lookup and legacy unpaired fallback are preserved.

The final predicate is exactly the tested additive policy: the original
spanning-read test passes its rejection threshold on both sides **OR** the
local intact-read test does so on both sides. It does not pool evidence types
across sides. Local means overlap with the existing same-contig 500-bp fetch
window, not a new family-specific distance. Record-order add/discard semantics
for intact keys follow the experiment. Global intact matches are never loaded.

The experiment used two bounded fetches per candidate; production gathers
both evidence types in one fetch. No entire-BAM read-name map is introduced.
The full benchmark measured the experimental implementation, not a fresh
end-to-end run of this integrated source. Focused tests validate the integration;
the full suite is the next checkpoint gate.

**127 focused tests pass** across family evidence, intact/spanning filters,
trimming, breakpoint pairing, real duplicate/LINE/SINE fixtures, TSD calling,
characterization, and memory regressions. New assertions cover actual nearby
reference-copy rejection, retention of remote/other-contig/clipped/wrong-mate
evidence, and the separate spanning/intact thresholds. Lint and diff checks
pass. Production filtering is no longer the old stabilized snapshot; this
change must be included explicitly in the planned precision commit.

## Reproducible analysis

```bash
.pixi/envs/default/bin/python scripts/summarize_intact_replay.py \
  --replay results/intact-fullreads-replay/2026-10-03-local \
  --stabilization results/family-replay/2026-10-02-stabilization \
  --output results/intact-fullreads-analysis/2026-10-05-local
```

Already run. The output directory includes summary, per-sample, paired,
changed-event, and stratified tables plus input hashes. Choose a new directory
to repeat; no previous result was overwritten. This analysis used small tables,
not new alignments or substantial caller runs on the login node.

## Checkpoint and release decision

This is a reasonable scientific change to include in a **development
checkpoint**, after full tests pass. It is not a basis for promising complete
R2 noninferiority or pinning a release: remaining FPs, paired R2-only detections,
whole-pipeline memory/runtime, and the broader feature audit remain open.
Follow the focused commit grouping in `2026-10-03-local-intact-plan.md`, with
the local filter as its own commit. Do not sweep untracked frozen runs and
large generated results into a blanket commit.

Prepared `scripts/run_checkpoint_tests.slurm` and its Python runner. They run
the full test suite on epyc (32 GB, four CPUs, eight hours; partition maximum
30 days verified), record executable paths and source/test-control hashes,
emit JUnit results, and fail the checkpoint if those inputs change during
execution. They do not claim to hash every large external fixture. Python
libraries use the existing project environment; bioinformatics tools use
available, versioned cluster modules. BLAT has no listed module and uses the
existing isolated project installation. Module availability was inspected;
actual executable/test validation occurs on the compute node.

Result path is `results/checkpoint-tests/<job-id>/`, with `.complete` created
only after successful tests and unchanged recorded inputs. Review any skipped
tests in the `-ra` output before claiming full capability coverage. No broad
alignment rerun is needed merely to repeat the completed candidate comparison.

Submission attempt returned exit 124 from the bounded `sbatch` call with no
numeric job ID. It remains unconfirmed. From the RelocaTE3 root, check for an
existing `r3-checkpoint-tests` job before submitting:

```bash
squeue -u "$USER" -o '%.18i %.45j %.10T'
sbatch --job-name=r3-checkpoint-tests scripts/run_checkpoint_tests.slurm
```

Run the second command only if no matching job exists. Next action: inspect
full-suite results for the recorded production source, then prepare focused
commits and the maintainer push/merge handoff. No release claim is implied.
