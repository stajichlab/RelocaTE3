# Validated development checkpoint

Recorded: 2026-10-05 14:27 PDT, America/Los_Angeles.
Project: RelocaTE3, `main`, base `942a5ea`, after merged PR #49.
Status: full-suite gate passed; changes remain uncommitted. No push, merge,
or release tag was performed.

## Verification

SLURM job 29419250 finished at 14:24:34 PDT on October 5. JUnit records
**331 tests, zero failures, zero errors, zero skips**, with 70.729 seconds
of test-suite time on host r27. `.complete` exists. The runner verified its
recorded input hashes before/after testing, and the subsequent review confirms
all 89 recorded source/test-control files still match. This fingerprint is
not a checksum of all large external biological inputs.

Artifacts:

- `logs/checkpoint-tests.29419250.log`
- `results/checkpoint-tests/29419250/full-tests.xml`
- `results/checkpoint-tests/29419250/provenance.json`

Executable provenance confirms module-provided bwa, bwa-mem2, samtools,
minimap2, bedtools, bowtie2, seqtk, and bcftools. Python/RelocaTE3 use the
existing project environment; BLAT uses its separate existing installation.
`git diff --check` also passes.

## Recommendation

Proceed to focused development commits and maintainer review/push. There is
no reason to hold all validated work uncommitted until every residual benchmark
difference disappears. This is a development checkpoint, not a release pin or
a promise of universal R2 noninferiority.

The 72-sample local-filter experiment removed 28 FPs without changing any truth
detection. The production integration then passed the full suite. These are
separate pieces of evidence: the complete original alignment pipeline was not
rerun across all 72 samples after integration. Remaining precision gaps,
whole-pipeline resource validation, and broader R2 utility coverage remain
release work, as documented in the result reports.

## Proposed commit boundaries

1. **Memory and streaming:** `aligners.py`, memory-related hunks of
   `librelocate.py`, memory/PSL tests, and their concise validation record.
2. **Family evidence:** new `family.py`, remaining family-related changes in
   `librelocate.py`, `models.py`, `insertions.py`, and `characterize.py`, plus
   family tests and portable LINE fixtures. Include analysis helpers imported
   by those tests so the commit works outside this working directory.
3. **Local intact-read precision correction:** the full-read filter/helper
   hunks of `insertions.py`, intact-read tests and SINE fixture, and the replay
   helper imported by those tests. Keep the rejected global rule experimental.
4. **Validation/documentation:** remaining reproducible benchmark auditing
   tools, full-suite runner, selected reports, and future-direction notes.

Some files contain multiple topics; use hunk-level staging or equivalent
reviewed patch splitting. Do not group solely by whole filename. The rejected
exclusive-coordinate pairing change is absent from production and must not
be accidentally restored when splitting patches. Duplicate-locus experiment
tests should remain identified as diagnostic history, not a production fix.

Review `.gitignore` separately. Avoid staging the untracked `results/` trees
or frozen source copies with a blanket add. Preserve those local artifacts;
no deletion is needed. Curate small fixtures and the required test helpers,
and retain provenance references for large benchmark evidence outside Git.

Next action: organize these focused commits, verify the final source remains
the tested source, and provide the maintainer the push/PR handoff. If splitting
requires any source change, rerun the relevant checks before handoff. A release
tag is not part of this checkpoint.
