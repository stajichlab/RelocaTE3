# Fresh full benchmark of merged RelocaTE3 checkpoint

Recorded: 2026-10-07, America/Los_Angeles.
Project: RelocaTE3, `main`, after merged PRs #50–#54.
Pinned revision: `213a0df99af71a57fffad726138593b152fc9a77`.
Status update, 2026-10-07 19:21 PDT: setup installed in the sibling benchmark
after write access became available. Source guard verifies revision 213a0df;
enumeration confirms 432 tasks, six callers, and the new output roots.
Seventeen coordinate tests pass against the installed scorer. No jobs
submitted by this session. Earlier installation instructions below remain
the reproducible procedure; repeating installation is unnecessary.

## Verified configuration and scope

The existing `config/benchmark.full-aligners.toml` already points all R3
callers to the current checkout; its rhome path resolves to this repository.
However, original completed run directories would be skipped, and scoring
still used the legacy anchor-only policy. The prepared checkpoint profile
addresses both issues without deleting the historical benchmark evidence.

- 72 samples: nine mPing, nine riceTElib, 54 diverged riceTElib.
- Six callers: R2; R3 BLAT/bwa aln, BLAT/BWA-MEM, BWA-MEM/BWA-MEM,
  Bowtie2/BWA-MEM, and minimap2/minimap2.
- **432 array tasks**, plus the dependent aggregation job.
- Existing eight-thread, 64-GB, 96-hour array resource requests; epyc's cached
  maximum is 30 days, so the partition/time combination is valid.
- Same 10-bp tolerance and explicit `tsd-interval` scoring for both tools.
  Genotype and exact TSD checks stay separate from positional detection.

References, annotations, TE libraries, and existing genome/TE indexes were
checked. No substantial indexing or read alignment was run on the login node.
The normal submission script performs its existing FASTQ/header preflight.
Existing truth and shared index caches are reused; actual caller pipelines
rerun into fresh directories.

## Prepared bundle

Local bundle:
`results/full-benchmark-setup/2026-10-07-213a0df-ready/`.

Profile to install:
`config/benchmark.checkpoint-213a0df.toml`.

Fresh benchmark output roots:

- `runs/checkpoint-213a0df/`
- `reports/checkpoint-213a0df/`

Original `runs/`, historical dataset reports, and default configuration are
not deleted or replaced. The new report root has its own dataset manifest;
the existing dashboard remains pointed at its previous root until explicitly
reconfigured after inspecting the new results.

The installer updates the scorer/config bridge/array scoring argument, installs
the matching-contract tests, and creates a private checkpoint R3 adapter.
Original caller adapter files are preserved, including their existing local
changes. The checkpoint adapter copies those files and adds source checks.
Legacy profiles retain legacy scoring through the default policy; the new
profile selects the explicit interval contract.

The guard checks the exact Git revision, tracked production/package files,
active benchmark files and configuration, and R2 adapter files. It fails a task
before expensive work if these differ. The R3 adapter sets PYTHONPATH to the
configured checkout and verifies its imported package path after activation.
It rechecks source integrity before creating its completion sentinel.
These checks pin source/settings, not every external FASTQ byte or executable
binary. The existing adapters retain their previously configured environments.
GNU time remains enabled for full adapter runtime and maximum RSS per task.

## Validation

The installer was exercised twice in a temporary benchmark copy to verify
idempotence, backups, source guards, and correct destinations. Python/shell
syntax checks pass; **17 coordinate-contract tests pass**. A deliberately
changed temporary profile was rejected by the guard. Config enumeration
confirms all six callers and 432 tasks. RelocaTE3's CLI help resolves in its
existing environment.

Preflight caught two setup-helper path issues (the prior report-root name and
temporary config-bridge imports); both were corrected before installation.
The successful bundle is the `-ready` directory. Failed preparation artifacts
remain under ignored `results/`; no benchmark files were changed by those
checks. Production scientific source remains unchanged from the merged
checkpoint. New setup scripts are local changes for this validation workflow.

## Commands for the maintainer

The benchmark is outside this session's writable workspace, so its installation
must be executed from the maintainer's shell. From the **relocate-benchmark
repository root**:

```bash
RT3=../../RelocaTE3_jason/RelocaTE3
python3.12 "$RT3/scripts/prepare_full_checkpoint_benchmark.py" install \
  --bundle "$RT3/results/full-benchmark-setup/2026-10-07-213a0df-ready"
CONFIG=config/benchmark.checkpoint-213a0df.toml \
  bash pipeline/submit_benchmark.sh --dataset full
```

Installation checks that relevant benchmark files still match their inspected
versions, backs up replaced pipeline/scoring files under
`cache/checkpoint_setup/<full-revision>/`, and refuses conflicting files or
preexisting checkpoint output roots. It prints its verified profile and the
submission command. If it refuses, inspect the reported difference rather
than deleting historical data or forcing an overwrite.

Submission should print array **0-431**, a numeric array job ID, and the
dependent aggregation job ID. An absent job ID is unconfirmed submission.
Monitor with:

```bash
squeue -u "$USER"
```

Keep the RelocaTE3 revision and benchmark settings unchanged during the run.
Do not commit/pull a newer revision mid-run; even an unrelated revision change
will intentionally fail the pin. After completion, compare all three panels,
somatic fractions, TE groups, coverage/divergence, runtime, and peak memory
against the newly rerun R2 results and the retained historical benchmark.
