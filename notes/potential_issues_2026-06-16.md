# Potential RelocaTE3 issues noticed while building the real-rice validation pipeline

Date: 2026-06-16
Context: Discovered while implementing `validation/real_rice/` (a reproducible
regression test that compares RelocaTE3 output to a frozen RelocaTE2 run on
the same 10 rice samples). None of the items below were the cause of the
validation-script failures we hit during that work — those were pipeline-
script bugs (sbatch cwd handling, pixi env hand-off to compute nodes). The
items here are RelocaTE3 issues to consider fixing in a future session.

---

## 1. `src/RelocaTE3/cli.py` and parts of `src/RelocaTE3/pipeline.py` are dead code

> **RESOLVED / OUTDATED (2026-07-11) — do not act on the description below.** The
> two CLI front-ends were consolidated (PR #19) and the canonical CLI was then
> relocated into `cli.py`. **Current reality is the opposite of what this item
> says:** `src/RelocaTE3/cli.py` is now the *single canonical CLI*, the console
> entry point is `relocaTE3 = "RelocaTE3.cli:main"`, and `src/RelocaTE3/__main__.py`
> is a thin `python -m` launcher. The subcommands/flags an agent should trust are
> the ones in `cli.py` (see `AGENTS.md` → "The CLI lives in one file"). Text kept
> for historical record only.

**Symptom.** `cli.py` registers subcommands (`run`, `find-insertions`,
`find-reference`, `characterize`) with a *different* CLI shape than what
the installed `relocaTE3` binary actually exposes. For example, `cli.py`'s
`run` advertises `--sample`, `--te`, `--genome`, `--repeatmasker`,
`--genotype`, `--len-cut-match`, `--len-cut-trim` and runs the entire
pipeline including reference-TE filtering and characterization. The real
binary (entry point goes through `src/RelocaTE3/__main__.py`) uses
`-n/--name`, `-T/--te-library`, `--min-match`, `--min-trimmed`,
`--mismatch`, etc., and `run` only chains map + trim — not align-genome
or find-insertions.

**Evidence.**
- `pyproject.toml` console script points at `RelocaTE3.__main__:main`, not
  `RelocaTE3.cli:main`.
- `pixi run python -c "from RelocaTE3 import cli; ..."` resolves to
  `src/RelocaTE3/cli.py`, but `relocaTE3 --help` shows the `__main__.py`
  subcommand set (`map`, `trim`, `run`, `annotate-ref`, `index-genome`,
  `align-genome`, `find-insertions`, `characterize`).
- `src/RelocaTE3/pipeline.py` imports `find_insertions`, `read_insertions_gff`,
  `write_insertions_gff`, `find_reference_insertions`,
  `write_existing_te_bed_from_rm`, `read_read_repeat`,
  `characterize_insertions`, `write_characterized` — those names don't exist
  in the current `insertions.py` / `reference_te.py` / `characterize.py` / 
  `genome_align.py` modules (which define `InsertionFinder`,
  `ReferenceTEAnnotator`, `Characterizer` instead).

**Why this matters.** Anyone reading the source to wire RelocaTE3 into a
new workflow (as I did to build the validation script) will read `cli.py`
first and infer the wrong API. It cost me an iteration in the validation
pipeline because I wrote `relocaTE3 run --sample ... --te ... --genome ...`
based on `cli.py`.

**Suggested fix.** Either delete `cli.py` and the now-unreachable
`run_sample`/`run_samples` in `pipeline.py`, or rewire the console script
to use them and bring them up to date. Pick one source of truth.

---

## 2. `find-insertions` does not emit a GFF compatible with RelocaTE2

**Symptom.** `find-insertions` writes only
`<outdir>/results/<target>.<te_name>.all_nonref_insert.txt` (custom
tab-separated columns). RelocaTE2 wrote
`<outdir>/repeat/results/ALL.all_nonref_insert.gff` as well, with the
attribute set:

```
ID=repeat_<chrom>_<start>_<end>;Name=<TE>;TSD=<TSD>;Note=...;
Right_junction_reads=<N>;Left_junction_reads=<N>;
Right_support_reads=<N>;Left_support_reads=<N>;
```

`characterize` does emit a GFF, but with a different schema
(`avg_flankers=`, `spanners=`, `type=`, `TE=`, `TSD=`) and only after step
7. There is currently **no RelocaTE3 step that produces the pre-genotyped,
RelocaTE2-compatible non-reference-insertion GFF**.

**Why this matters.**
- The validation pipeline's `normalize_relocate3.py` had to parse the
  custom `.txt` instead of a GFF. Anyone downstream (workflow engines,
  comparison scripts, mPing-tracking tools the lab already runs) expecting
  the RelocaTE2 layout will need a similar shim.
- The legacy schema is what the field has standardized on and what the
  bundled `validation_data/real_rice/relocate2_results/*/repeat/results/ALL.all_nonref_insert.gff`
  uses for the regression baseline.

**Suggested fix.** Have `find-insertions` also write a parallel
`<outdir>/results/<target>.<te_name>.all_nonref_insert.gff` in the
RelocaTE2 attribute schema (the data is already computed — see
`_write_event_start` in `insertions.py`, which has `right_count`,
`left_count`, `top_tsd`, `te_orient`, `repeat_family`). Adding the GFF
writer next to the existing TXT writer is mechanical.

---

## 3. `find-insertions --tsd` is required and global

**Symptom.** The `find-insertions` subcommand requires `--tsd MOTIF`
(e.g. `TTA`). Help text reads: `"TSD-unknown mode is not yet supported"`.

**Why this matters.**
- The legacy `ALL.all_nonref_insert.gff` includes per-call TSDs like `TAA`,
  `TTA`, `TGA`, `insufficient_data` — that's because mPing junctions in
  real data don't all share one motif. Requiring a single global motif
  means RelocaTE3 will either reject or mislabel non-canonical TSDs.
- A run can only handle one TE family per invocation. The legacy
  pipeline could batch TE libraries.

**Suggested fix.** Implement TSD auto-detection from the bracketing reads,
already half-supported by the dominant-TSD selection logic in
`_write_event_start`. At minimum, allow `--tsd auto` to pick per-cluster
TSDs and treat overlap-less or ambiguous calls as `insufficient_data`
(matching the RelocaTE2 label).

---

## 4. Subcommand naming: `run` does only map + trim

**Symptom.** `relocaTE3 run` chains `map` + `trim` and stops. Producing
insertion calls requires the user to also chain `align-genome` →
`find-insertions` (and optionally `characterize`) by hand.

**Why this matters.** Most users will read "run" as "run the pipeline."
A new collaborator setting up a single-sample test will assume `run`
suffices and be surprised when no insertion table appears. (This is
exactly what I did before reading the help carefully.)

**Suggested fix.** One of:
- Rename `run` → `run-trim` (or `map-trim`) so the partial scope is
  explicit, and add a new `pipeline` / `all` subcommand that runs every
  step end-to-end.
- Or keep `run` but extend it to chain through `find-insertions` when the
  required inputs (`--genome`, `--reference-ins`, etc.) are supplied.

---

## 5. No single-deliverable entry point for one sample

**Symptom.** Producing the per-sample `all_nonref_insert` file requires
four CLI calls (`index-genome`, `run`, `align-genome`, `find-insertions`),
each with their own paths to thread through. The bundled
`validation/real_rice/run_relocate3.sh` exists almost entirely to glue
these together for one sample.

**Why this matters.** Same usability concern as #4. Workflow engines
(Nextflow / Snakemake) like the step granularity; humans running a quick
test do not.

**Suggested fix.** Tied to #4 — a `pipeline` subcommand whose flags are
just `--left/--right/--te-library/--genome/--reference-ins/--name/
--outdir/--tsd/--threads` and which orchestrates all four internally.
Would let the validation script collapse to a single `relocaTE3 pipeline
...` call per sample.

---

## 6. `characterize` GFF schema diverges from `find-insertions` output

**Symptom.** `characterize` writes
`<stem>.characTErized.gff` with attributes
`ID=<chrom>.<pos>.spanners;avg_flankers=...;spanners=...;type=...;TE=...;TSD=...`,
while `find-insertions` (if/when it gets a GFF — see #2) would naturally
want to use the RelocaTE2 attribute set. Tools reading both for the same
sample will see two different schemas.

**Suggested fix.** Once #2 is done, extend the same attribute set in
`characterize` (keep `avg_flankers`, `spanners`, `type` but add the
`Left_/Right_junction_reads` and `Note=` fields the legacy tools expect).

---

## 7. `pipeline.py` provenance helpers reference a nonexistent module API

**Symptom.** `src/RelocaTE3/pipeline.py` imports
`from RelocaTE3.align import Aligner` (fine) but also
`from RelocaTE3.genome_align import align_to_genome, read_read_repeat` and
`from RelocaTE3.insertions import find_insertions, write_insertions_gff,
write_insertions_txt` — none of those module-level symbols exist in the
current codebase (the modules expose classes `Aligner`, `InsertionFinder`,
etc., not module-level functions of those names).

**Why this matters.** `import RelocaTE3.pipeline` will raise
`ImportError` at import time, so anything that touches the module breaks.
Tied to #1: this whole file appears to be leftover from a previous API
design.

**Suggested fix.** Delete the file or port it to the current
class-based API.

---

## 8. `map` step appears to be racy when the output dir is concurrently mutated

**Symptom.** During a local `pixi run validate-rice --local --force B_10`
run, the `map` step produced `B_10.right.bam` (667 KB, indexed OK) but
failed with:

```
[E::hts_open_format] Failed to open file "<outdir>/B_10.left.bam" : No such file or directory
samtools index: failed to open "<outdir>/B_10.left.bam": No such file or directory
RelocaTE3 ERROR Command '['samtools', 'index', '.../B_10.left.bam']' returned non-zero exit status 1.
```

i.e. the left-direction BAM was either never written or was removed before
`samtools index` ran on it, while the right-direction BAM came through fine.

**Why this matters.** The trim/find-insertions steps require BOTH 5'/3'
flanking FASTQs, so a half-map failure silently breaks the downstream
pipeline.

**Suggested investigation.** Look at `src/RelocaTE3/align.py`
(`map_minimap_library` / its direction loop) for:
- Whether the two directions write through a temp path that could collide.
- Whether the indexing step assumes the previous step's output without
  verifying.
- Whether the failure should abort the whole `map` step or be retried.

A repro should be straightforward by running `relocaTE3 run -l R1 -r R2 -T
te.fa -n NAME -o OUT --threads N` against the rice validation FASTQs and
watching whether one direction races the other.

---

## 9. `find-insertions` writes the `--target` filter value into the chromosome column

**Severity.** High — it silently destroys the chromosome of every insertion
in the output file. Any downstream comparison keyed on chromosome (the
real-rice validation pipeline, any RelocaTE2-compat tooling) sees zero
overlap with the legacy calls because every RelocaTE3 row's "chrom" is the
string `"ALL"` (or whatever was passed to `--target`).

**Symptom.** Raw output at
`<outdir>/results/<target>.<te_name>.all_nonref_insert.txt` has the literal
target value in column 4 instead of the chromosome:

```
mPing  TTA        B_10  ALL  1041521..1041523  + T:2 R:1 L:1 ST:0 SR:0 SL:0
mPing  singleton  B_10  ALL  1623995..1623997  - T:1 R:1 L:0 ST:0 SR:0 SL:0
```

Reproduced on `B_10` from the real-rice validation set:
- R2 calls: 650 mPing insertions, `chrom ∈ {Chr1..Chr12}`.
- R3 calls: 617 mPing insertions, `chrom = "ALL"` on every row.
- `compare_calls.py` matches on `(chrom, te_name)` → 0 shared even though
  the call counts (650 vs 617) are within ~5% of each other and the actual
  coordinates almost certainly do overlap heavily.

**Root cause.** `src/RelocaTE3/insertions.py:421` in `_write_event_start`:

```python
out.write(
    f"{repeat_family}\t{tsd_field}\t{sample}\t{target}\t{coor_start}..{coor}\t"
    f"{te_orient}\tT:{total_count}\tR:{right_count}\tL:{left_count}\t"
    f"ST:0\tSR:0\tSL:0\n"
)
```

`target` is the CLI `--target` *filter* (`"ALL"` or a single chromosome
name), not the cluster's actual chromosome. The actual chromosome IS
tracked by `cluster_chrom: dict[int, str]` (built at
`insertions.py:91` and populated at `insertions.py:191` with
`cluster_chrom[count] = chro`), and `cluster_chrom` is even passed into
`_write_output` (insertions.py:358) — it just isn't threaded into
`_write_event_start`, so the row ends up using the filter as a stand-in
for chrom.

**Fix.**

1. Thread the chromosome into the writer:
   - In `InsertionFinder._write_output` (around `insertions.py:350`),
     look up `cluster_chrom[event]` for each event and pass it to
     `_write_event_start` as a new `chrom` argument.
   - In `InsertionFinder._write_event_start`, replace `{target}` in the
     f-string at line 421 with `{chrom}`. Drop the `target` parameter
     from `_write_event_start` if it's no longer used.

2. Sanity check after the fix:
   - Re-run `relocaTE3 find-insertions ... --target ALL` on `B_10`; the
     4th column should now be `Chr1`/`Chr2`/... .
   - Re-run `relocaTE3 find-insertions ... --target Chr10`; the filter
     should still restrict to Chr10 clusters but the output should now
     literally print `Chr10` per row (it already did, but for the wrong
     reason — verify it still does post-fix).

3. Once #2 (RelocaTE2-compatible GFF) is implemented, make sure the GFF
   writer uses the same `chrom` lookup, not `target`.

**Why we are fixing this.**
- The real-rice validation pipeline (`validation/real_rice/`) cannot
  measure recall, precision, or Jaccard against the legacy RelocaTE2
  baseline until R3's chrom column is correct. Right now shared=0 by
  construction.
- Any user running `--target ALL` (the most common case — whole-genome
  insertion scan) gets output that's effectively unusable without
  reprocessing the BAM.
- The bug is cheap to fix and has no behavioral knock-on: `cluster_chrom`
  is already populated correctly, we're only changing what gets written.

**Notes for the future debugging session.**
- Tests: `tests/` does not currently assert the chrom column value of the
  `all_nonref_insert.txt` writer; add a regression test that runs
  `find-insertions --target ALL` against a tiny synthetic BAM with reads
  on at least two different reference contigs, then asserts the output
  has both contig names in column 4 (and no `"ALL"` rows).
- Related code locations:
  - `src/RelocaTE3/insertions.py:91` — `cluster_chrom` declared.
  - `src/RelocaTE3/insertions.py:191` — `cluster_chrom[count] = chro` is
    where the real chromosome is recorded.
  - `src/RelocaTE3/insertions.py:354..374` — `_write_output` signature
    + call into `_write_event_start`.
  - `src/RelocaTE3/insertions.py:376..424` — `_write_event_start` body.
- After fixing, `validation/real_rice/normalize_relocate3.py` requires
  no change — it already reads column 4 as `chrom`. Just re-run
  `pixi run validate-rice --skip-run` to regenerate
  `report/relocate3_calls.tsv` and the comparison.
- This is independent of issue #2 (no RelocaTE2-compatible GFF emitted).
  Fix #9 first; #2 can re-use the same `chrom` plumbing when a GFF
  writer is added.

---

=======

## Lower priority observations

- `relocaTE3` always prints a `pixi` netfs-redirect WARN at startup
  (rattler cache redirected to `$SCRATCH`). Not a RelocaTE3 bug — it's
  pixi telling us to set `[cache.repodata]` or `PIXI_CACHE_DIR`. Worth
  documenting in `README.md`/pixi config so SLURM logs don't include the
  noise on every run.
- `align-genome` writes `<name>.repeat.minimap.sorted.bam` — the literal
  string `repeat` is hardcoded. If a user wants to track multiple TE
  families in one project they'll get filename collisions. Possibly a
  knob to add.
- `relocaTE3 --version` and `relocaTE3 -V` both work, but the global
  `--verbose` flag exists only on subcommands (via `_add_common_args`)
  rather than the top-level parser, so `relocaTE3 -v <cmd>` errors. Minor.

---

## What was NOT a RelocaTE3 bug

Recording these so a future debugger doesn't re-investigate them:

- `ModuleNotFoundError: No module named '_config'` inside a SLURM job —
  that was the validation script using `BASH_SOURCE[0]` to find
  `_config.py`. Under `sbatch`, the script is copied to
  `/var/spool/slurmd/...`. Fix lives in
  `validation/real_rice/run_relocate3.sh`: anchor on the config path's
  directory instead.
- `relocaTE3: command not found` on compute nodes — pixi env not
  inherited by SLURM jobs. Fix lives in the same script: re-exec via
  `pixi run --manifest-path <pixi.toml> -- bash $0 "$@"` when the binary
  isn't already on PATH.
- The first round of "wrong flags" — I trusted `src/RelocaTE3/cli.py`
  rather than `__main__.py`. See item #1.
