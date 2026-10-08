# Local intact-read validation and development checkpoint

Completed follow-up: [local filter results](2026-10-05-local-intact-results.md).
All 72 tasks passed; 28 FPs were removed without losing a truth detection.
The local rule is now integrated, with the full checkpoint suite pending.

Recorded: 2026-10-03 14:36 PDT, America/Los_Angeles.
Project: RelocaTE3, `main`, base `942a5ea`, following merged PR #49.
Status: local candidate prepared; production filtering unchanged. No commits,
push, merge, or release pin performed.

## Decision

Finish this focused comparison before deciding whether the precision change
belongs in production. Then organize a development checkpoint; do not let
unrelated future features extend the uncommitted work indefinitely. A reviewed
development merge is distinct from a release asserting R2 noninferiority.

The global intact-read rule removed 44 FPs but lost five real mPing detections.
All five were R3 advantages over R2. Three bounded inspections found no local
intact-read evidence despite the global veto. Locality is therefore a causal
hypothesis to test, not an established solution or a family-specific exception.

## Experimental rule

`scripts/replay_intact_fullreads.py` now accepts `--scope local` during
preparation and freezes the scope in its manifest. Old manifests retain global
semantics; existing completed snapshots were not modified.

The local rule counts nearly complete full-read alignments only when they
overlap the existing same-contig window:
`[max(0, min(start,end)-500), max(start,end)+500)` in BAM fetch coordinates.
This reuses production's existing FULLREAD_WINDOW; it was not fitted to a
family, truth position, or the five lost events. Overlap with this window is
not a requirement to cross the breakpoint. CIGAR coverage, 30% per-side veto,
and mate identity remain as in the preceding experiment. The existing
spanning-read veto remains active.

Local mode skips the whole-BAM scans and global read-name set used by the
global experiment. Its extra lookup is bounded per candidate and retains
only local intact keys. It currently performs an additional local fetch beside
the existing spanning-read check; consolidating those fetches is a possible
implementation detail only if the policy earns adoption.

The experimental filter replaces a function only within the replay process;
no production caller source was changed. Each task first requires exact
reproduction of the stabilized baseline, then reruns calling/characterization
with the local rule and scores against truth. Candidate names in output remain
`intact_candidate`; the manifest/comparison `scope` distinguishes this run.

## Verification and acceptance

Twenty focused tests pass across intact and existing false-junction behavior.
The actual Os3912 SAM fixture is rejected by the local candidate. An indexed
synthetic BAM verifies that a nearby intact read qualifies, a nearby heavily
clipped read does not, and remote or other-contig alignments
are excluded by the genomic window. Existing tests preserve mate identity and
empty-side handling. Lint, shell syntax, and `git diff --check` pass.

The full 72-sample replay must:

- Reproduce stabilized baseline calls in every task.
- Preserve all five named mPing regression targets (TE000316, TE000125,
  TE000157, TE000434, TE000452 in their previously documented samples).
- Report all other gained/lost detections and FP changes across panels,
  coverage, divergence, TE groups, and somatic fractions.

Passing the fixture alone is insufficient. Even a successful replay would not
establish whole-pipeline runtime/memory or certify the complete R2 utility
surface. The full external-tool suite remains a separate pre-merge check.

## Prepared job

```bash
.pixi/envs/default/bin/python scripts/replay_intact_fullreads.py prepare \
  --baseline results/family-replay/2026-10-02-stabilization \
  --scope local \
  --output results/intact-fullreads-replay/2026-10-03-local
```

Already prepared: all snapshot hashes verified, 72 tasks. Input paths and
size/mtime checks plus small input hashes follow the earlier replay. The epyc
eight-hour request was checked against its cached 30-day maximum. All caller
runs remain on SLURM. No old output deletion is needed.

## Proposed checkpoint scope after the result

Use explicit staging and reviewable commits; the working tree includes many
untracked frozen runs and analysis outputs that should not be swept into a
blanket commit.

1. Validated memory/streaming changes and their focused tests, with resource
   measurements and limitations recorded.
2. Family-evidence preservation and fallback correction, metadata contracts,
   portable fixtures, and tests. Keep historical pairing, whose restoration
   recovered 30 detections. Exclude the rejected global filter from production.
3. Reproducible validation scripts and concise result/roadmap documentation;
   curate this scope separately from large generated benchmark outputs.
4. Only if supported by this comparison: a separate local-filter commit,
   rather than mixing another scientific behavior change into the prior work.

Run the full external-tool tests on compute nodes for the exact proposed
checkpoint. Review active CLI/library behavior, release limitations, and
dependency setup before handing the maintainer push commands. This is a
commit-organization recommendation, not authorization inferred to push or
merge from this session. The maintainer retains the push/merge step.

## Submission status

The single submission attempt returned exit 124 after the bounded timeout,
without a numeric job ID. Submission is unconfirmed. From the RelocaTE3 root,
check for `r3-local-intact` before repeating it:

```bash
squeue -u "$USER" -o '%.18i %.45j %.10T'
```

If no matching job exists:

```bash
REPLAY=results/intact-fullreads-replay/2026-10-03-local
sbatch --job-name=r3-local-intact --array=0-71%4 \
  "$REPLAY/replay_intact_fullreads.slurm" "$REPLAY"
```

Next action: run this frozen comparison and inspect its paired results before
deciding on a production filtering change. No new permission decision is
needed to analyze it when complete.
