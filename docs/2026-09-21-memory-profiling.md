# Memory profiling: preserved baseline and single-sample replay

Recorded: 2026-09-21, 19:05 PDT (America/Los_Angeles).
Project: RelocaTE3, branch `main`, parity work merged in PR #49.
Algorithm source: `942a5ea266bd000a47e27df7ef9e2536c877acbd`.
Status: preparation and lightweight validation complete; SLURM submission
unconfirmed. No algorithm changes, deletion of benchmark outputs, or release
tagging were performed.

## Purpose and selected sample

Locate the stage/process responsible for high memory before changing data
structures or insertion filters. Use the full benchmark's riceTElib
`cov30x_rep1`, caller `relocate3-blat-bwaaln`, with eight threads.

Its completed benchmark used 35,431,884 KiB peak RSS (33.79 GiB) and 40,281
seconds (11.19 hours). The preserved normalized calls contain 302 true
positives and 66 false positives (368 calls), against 500 truth events.
The new replay is not another 432-task benchmark.

Prepared output directory, relative to the RelocaTE3 repository:

```text
results/memory-profile/2026-09-21-ricetelib-cov30x-rep1/
  manifest.json       source revisions, local changes, parameters, snapshot hashes
  baseline/           existing call tables, truth, scoring, panel resource table
  snapshot/           frozen Python source, adapter, scorer, config, environment locks
  scorer-preflight/   saved-call rescoring check; no new alignments
```

The snapshot contains 61 checksummed files. It is protected by hash verification
before replay, not by filesystem read-only permissions. RelocaTE3 imports the
snapshot through an explicit `PYTHONPATH`; the existing benchmark environment
provides dependencies and installed package metadata. The wrapper verifies the
import location and records actual tool paths/hashes/versions on the compute
node. It does not install a new profiler or change alignment settings.
External environment directories are reused, not frozen copies; their lock
files and runtime executable hashes provide provenance. Do not update those
environments while profiling is running.

Source scripts:

- `scripts/profile_benchmark_memory.py`: prepare, monitor, normalize, score,
  and compare one replay.
- `scripts/memory_profile_adapter.sh`: benchmark environment activation,
  tool provenance, frozen adapter execution.
- `scripts/profile_benchmark_memory.slurm`: eight CPUs, 64 GB, 24 hours on
  `epyc`, logs under `logs/`; captures `SLURM_SUBMIT_DIR`.
- `tests/memory_profile_test.py`: six lightweight harness tests.

Generated profiling outputs are ignored by Git under `results/memory-profile/`.
The algorithm under `src/` was not edited. Only the profiling harness, tests,
documentation, and the targeted ignore entry were added.

## Validation performed

```bash
python3.12 scripts/profile_benchmark_memory.py --help
bash -n scripts/memory_profile_adapter.sh scripts/profile_benchmark_memory.slurm
.pixi/envs/default/bin/python -m pytest -q tests/memory_profile_test.py
ruff check scripts/profile_benchmark_memory.py tests/memory_profile_test.py
git diff --check
```

All six harness tests pass: order-independent/duplicate-aware table comparison,
process-session filtering, parent/child memory capture, failure propagation and
overwrite refusal, login-node replay refusal, and exclusive metadata creation.
Shell syntax and the targeted Ruff check pass. No full bioinformatics tests
were run on the login node.

The frozen scorer was also run on the 368 already saved calls. `matches.tsv`,
`precision.tsv`, `correctness.tsv`, and `false_positive_calls.tsv` all match the
existing sample report exactly after ignoring row/column order. This is a
lightweight scorer check, not a new caller result.

Environment preflight resolved Python, relocaTE3, BLAT, bwa, minimap2, samtools,
and seqtk, and confirmed import from the source snapshot. The existing
benchmark activation script emitted a bcftools module logger/socket warning in
this restricted session, but completed successfully. Runtime provenance is
checked again on the compute node. The separate BLAT test environment's earlier
seqtk skips do not imply seqtk is missing from this benchmark environment.

## Submission status and command

`scontrol show partition epyc` failed with a controller connection error.
The readable local configuration at
`/var/spool/slurmd/conf-cache/slurm.conf` specifies `epyc MaxTime=30-00:00`;
the requested 24-hour wall time fits that limit. This cached configuration
does not prove current scheduler availability.

The following submission was attempted with a 20-second timeout and returned
exit 124 without a numeric job ID:

```bash
sbatch --parsable scripts/profile_benchmark_memory.slurm \
  results/memory-profile/2026-09-21-ricetelib-cov30x-rep1
```

No `started.json` or `logs/relocate3-memory.*.log` was present at the final
check. Submission is **unconfirmed**, not running or completed. Check the
queue from an ordinary cluster terminal before repeating the command, since
a client timeout alone cannot prove the scheduler rejected a request:

```bash
squeue -u "$USER" -o '%.18i %.35j %.10T %.10M %.30R'
```

If no matching profiling job exists, run the `sbatch` command above from the
RelocaTE3 repository root. From the sibling `relocate-benchmark` root, first:

```bash
cd ../../RelocaTE3_jason/RelocaTE3
```

Do not remove the previous full benchmark. The prepared directory is already
ready; there is no need to rerun `prepare`. An exclusive `started.json` marker
prevents concurrent duplicate jobs from running the same replay. A failed run
retains its logs and cannot silently overwrite them; use a new output directory
for an intentional retry after investigating the failure.

For a future, separate baseline preparation (choose a new directory):

```bash
python3.12 scripts/profile_benchmark_memory.py prepare \
  --output results/memory-profile/NEW-RUN-NAME
```

## What the job records and how to interpret it

Large input/index SHA-256 checksums are calculated on the compute node, outside
the measured adapter runtime. The actual adapter is wrapped with GNU `time -v`
and sampled every two seconds through Linux `/proc`, including child aligners
and other processes in its session.

Expected new outputs:

- `input_sha256.json`, `tools.txt`, `started.json`: runtime provenance.
- `adapter.time-v.txt`: whole-adapter resource usage comparable to the benchmark.
- `profile/adapter.log`: existing stage-start messages and command output.
- `profile/process_memory.tsv`: timestamped per-process RSS and high-water marks.
- `profile/summary.json`: per-process peaks and maximum sampled RSS sum.
- `replay/`, `score/`: isolated new caller outputs and truth scoring.
- `comparison.json`: normalized-call and scored-table equality against baseline.
- `.profile_complete`: created only after successful execution and all comparisons.

Compare process RSS timestamps with the adapter's stage messages/command lines
to distinguish BLAT, Python map/trim, genome alignment, insertion finding,
full-read alignment/sorting, and characterization. The first pass does not
attribute allocations to individual Python lines or data structures. Once the
dominant process/stage is measured, a targeted allocation trace can be scoped
to that stage without repeating every expensive alignment.

RSS sums double-count shared memory and are not unique physical memory. A
two-second sampler can miss short-lived processes/peaks; GNU time supplies an
additional kernel-recorded high-water metric. A process's VmHWM can survive
exec, so use sampled RSS, not a newly named command's inherited HWM alone, to
attribute memory to a command. Resource measurements include profiling-run
conditions and should not be treated as identical hardware/cache conditions
to the historical benchmark.

If calls differ, the job preserves the differences and fails the final gate;
do not attribute changes to a memory optimization that has not occurred.
Investigate environment, ordering, or historical-source differences first.

## Hypotheses, not conclusions

Static inspection found several plausible memory consumers: the BLAT
sequence-recovery path buffers sequence output, TE parsing retains selected
per-read dictionaries, and genome-alignment planning retains flank/original
read records. These are investigation targets, not established causes of the
33.79-GiB peak. In particular, the no-seqtk path reads the complete query FASTA
into memory, but seqtk is present in the currently resolved benchmark
environment; that fallback must not be assumed to explain this benchmark.

Next action: obtain a confirmed SLURM job ID and complete this isolated
profile; then select a measured, behavior-preserving memory optimization.
False-positive/false-negative algorithm changes remain deferred until the
memory baseline has been measured and preserved.
