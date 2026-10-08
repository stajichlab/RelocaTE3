# LINE read-family provenance extraction

Update: SLURM job 29331469 completed the extraction successfully. The
[completed comparison](2026-10-01-line-family-trace-results.md) confirms that
eight equally scoring read hits switch family solely through tie-breaking.
The submission uncertainty below records the earlier session state.

Recorded: 2026-10-01 13:35 PDT, America/Los_Angeles.
Project: RelocaTE3, `main`, following merged PR #49.
Status: extraction implemented and lightly tested; cluster submission unconfirmed.
No production caller changes or new biological conclusions.

## Purpose and evidence

Trace the eleven shared junction reads for TE000172 in
`div005_rep03_cov30x`. Both callers report the same breakpoint and counts,
but R2 reports Os3328 while R3 reports Os0596, with R3 junction votes 9:2.

Source inspection identifies a candidate explanation: R2
`relocaTE_trim.py:parse_align_blat` ranks boundary contact then match count,
keeping the first alignment on exact ties. R3 `librelocate.py:_match_rank`
uses those primary signals followed by deterministic target-name/coordinate
tie-breaks. This could produce systematic family attribution differences at
an ambiguous LINE, but **the actual selected hits have not yet been extracted**.
Do not change tie-breaking based on this hypothesis alone.

## Implementation

`scripts/trace_line_families.py` reads the existing indexed-junction evidence,
extracts R2 chromosome and original chunk family mappings, selects all PSL
hits from relevant R2 chunks, extracts R3 mappings, and scans the archived R3
left/right TE-alignment BAMs for all eleven reads. It preserves original row
order and line numbers, SAM headers/records, input text hashes, BAM size/mtime,
script hash, and pysam version. It checks coverage of all eleven read names
in both callers' mappings and hits before writing `.complete`.

R2 PSL files use `.blatout`, not `.psl`. R3 temporary PSL files are absent;
its retained TE BAMs are the available evidence and may not preserve original
PSL ordering. Both R3 BAMs are approximately 2.6 GB. Large scans therefore
run on SLURM, not on the login node. The script uses streaming selection,
refuses existing output directories, and never writes benchmark inputs.
Failures may leave an incomplete output directory; inspect it and choose a
new output name for a retry rather than overwriting it.

The eventual `family_evidence.json` is a portable candidate regression
fixture. It is **not yet generated**, and no biological regression is claimed
to pass from these unextracted data. The next analysis should compare exact
per-read hit scores, eligible hit sets, emitted family labels, and name lookup;
then test the causal explanation using the small extracted records.

## SLURM and validation

`scripts/trace_line_families.slurm`: short partition, one CPU, 4 GB, one hour.
Cached `/var/spool/slurmd/conf-cache/slurm.conf` states short MaxTime=120
minutes, validating this time request. Runtime captures `SLURM_SUBMIT_DIR`;
there is no `--chdir` directive or script-path-derived output path.
Existing project Python provides pysam; no aligner or module-provided binary
is invoked. SLURM output goes to `logs/line-family-trace.<jobid>.log`.

CLI help, shell syntax, and `git diff --check` passed. Two synthetic tests
passed (exact-name matching, duplicate-hit preservation, text hashing,
PSL column selection, and header handling):

```bash
.pixi/envs/default/bin/python -m pytest -q tests/line_family_trace_test.py
```

`sinfo` timed out; cached partition configuration was used. A bounded
`sbatch --parsable scripts/trace_line_families.slurm` attempt timed out
(exit 124) without a job ID. Queue verification also failed to return promptly;
no matching SLURM log was present at inspection. Submission is **unconfirmed**,
not reported as running or completed. No second submission was attempted.

## Next action

From the RelocaTE3 repository root, first inspect the queue for any already
accepted `trace_line_families.slurm` job to avoid duplicate submission:

```bash
squeue -u "$USER" -o '%.18i %.40j %.10T'
```

If none exists and the default output directory has not been created, submit:

```bash
sbatch scripts/trace_line_families.slurm
```

After completion, inspect `logs/line-family-trace.<jobid>.log` and
`results/error-audit/2026-10-01-line-family-trace/.complete`. Resume the
per-read comparison from `family_evidence.json`. No benchmark rerun is needed.
