# User-Selectable Outputs — Design

Date: 2026-07-22
Status: approved (design); implementation plan to follow
Scope: RelocaTE3 (the tool). The relocate-benchmark repo is downstream and is
NOT changed by this work; it keeps using the granular subcommands and only picks
this up if/when it is re-pinned.

## Goal

Let a RelocaTE3 user control which outputs a run retains — keep only the final
insertion calls, or keep specific intermediate artifacts (e.g. the reads that
contain the TE) — instead of always writing everything to disk. Motivating case:
a user who only wants "how many reads contain the TE," not the final calls.

This is deliberately **not** tied to the `--verbose` flag (which is debug logging
only, and partly stubbed). Output selection is its own explicit system.

## Decisions (locked)

- **Selection model:** a named-category **keep-list**. `--keep cat1,cat2,...`.
- **Default (omitted `--keep`):** `final` only.
- **Legacy behavior:** `--keep all` keeps everything (today's behavior).
- **Drop semantics:** the pipeline runs normally and writes every artifact; on
  **successful** completion, the files for any category NOT in the keep-set are
  deleted. (Later steps still find their inputs mid-run; on error nothing is
  deleted, so failures stay debuggable.)
- **Where it lives:** a NEW end-to-end `relocaTE3 pipeline` command. The existing
  8 subcommands are unchanged. (Rationale below.)
- **Always-on summary:** a small `<name>.summary.tsv` is always written and is
  exempt from cleanup, regardless of `--keep`.
- **Aligner passthrough:** the new command supports `--te-aligner` /
  `--genome-aligner`, which requires threading those params through
  `run_sample()` (they are currently only wired into the `run` / `align-genome`
  subcommands).

## Why a new command (not a flag on existing subcommands)

RelocaTE3's CLI is 8 granular per-step subcommands (`map`, `trim`, `run`,
`annotate-ref`, `index-genome`, `align-genome`, `find-insertions`,
`characterize`). Users — and the benchmark — chain them, and each step's output
is the next step's input. There is **no single end-to-end CLI command today**.
`run`/`trim` produce ONLY intermediates (flanking, te_containing) and have no
"final" output of their own, so a "keep final only" default cannot sit on them
(it would delete the flanking reads that are the entire product of `run`, which
the benchmark's next step consumes).

A keep-list with a "final-only" default only has a coherent meaning where one run
produces both intermediates and the final result. So it belongs on a new
one-shot command, which also fills a genuine usability gap. `pipeline.run_sample()`
already orchestrates steps 3→5 (+ optional 0/6 and 7), so the command is mostly
wiring.

## The `relocaTE3 pipeline` command

```
relocaTE3 pipeline --left R1.fq [--right R2.fq] --te-library te.fa \
    --genome ref.fa --name SAMPLE --outdir OUT \
    [--repeatmasker rm.out] [--genotype] [--threads N] \
    [--te-aligner {minimap2,bwa,bwamem2,bowtie2,blat}] \
    [--genome-aligner {minimap2,bwa,bwamem2,bowtie2}] \
    [--minimum-match-length N] [--minimum-trimmed-length N] [--mismatch-allowance N] \
    [--keep final|te_containing|flanking|read_repeat|te_portions|alignments|all]
```

Runs, in one process: step 3 (identify/trim), step 4 (genome align), step 5
(find insertions), optional steps 0/6 (reference annotation) when
`--repeatmasker` is given, optional step 7 (characterize) when `--genotype` is
set. Writes `<name>.summary.tsv`, then applies cleanup.

## Output categories

Comma-separated; validated against this set (unknown name → hard error). Special
tokens: `final` (default), `all` (keep everything).

| Category        | On-disk (under `<outdir>/`)                     |
|-----------------|--------------------------------------------------|
| `final`         | `results/` — nonref/ref insert GFF+txt, characterized |
| `flanking`      | `flanking/*.flankingReads.fq`                    |
| `te_containing` | `te_containing/*.ContainingReads.fq`             |
| `read_repeat`   | `te_containing/*.read_repeat_name.txt`           |
| `te_portions`   | `te_portions/*.fa`                               |
| `alignments`    | `genome_aln/*.bam`                               |

`te_containing` and `read_repeat` share the `te_containing/` directory, so
cleanup must delete **per-category file globs**, not whole directories, and must
not remove a shared dir while another kept category still needs it.

## Summary file

`<outdir>/<name>.summary.tsv` — top-level (outside `results/`, so it survives
even when `final` is dropped), always written, never cleaned. Counts are captured
from **in-memory** pipeline results during the run (not by re-reading files that
cleanup may delete):

| metric                    | source                                    |
|---------------------------|-------------------------------------------|
| `reads_te_containing`     | count of TE-containing reads (step 3)     |
| `flanking_reads_written`  | `identify_TE_reads` return value          |
| `nonref_insertions`       | len(insertions) from step 5               |
| `nonref_insertions_<class>` | per-class counts (hom/het/somatic/excision) when genotyped |

(`reads_te_library_mapped` included if cheaply available from the map step;
otherwise omitted rather than adding a second pass.)

## Cleanup mechanics & safety

- Resolve keep-set from `--keep` (default `{final}`; `all` short-circuits, no
  deletion).
- Run only after `run_sample()` returns successfully.
- For each non-kept category, delete its glob(s). Idempotent; guard missing
  paths; never touch the summary file.
- On any exception, skip cleanup entirely.

## Backward compatibility / benchmark impact

- Existing subcommands: unchanged. The benchmark composes `run` /
  `align-genome` / `find-insertions` and is therefore **unaffected**.
- The one edit to existing code is extending `run_sample()` to accept and pass
  `te_aligner` / `genome_aligner`; its current callers keep working via defaults.
- After merge, the benchmark can be re-pinned to the new rev to gain the command;
  no benchmark behavior change is required.

## Testing

- Unit: `--keep` parsing (default `final`, `all`, comma-list, unknown → error);
  cleanup over a populated fixture dir (correct globs removed, kept categories and
  the summary survive, shared `te_containing/` handled); summary contents/format.
- Integration: run `pipeline` on the tiny `tests/data` set with two `--keep`
  values; assert kept vs dropped files, a present+valid `summary.tsv`, and a
  non-empty nonref insertion call.

## Out of scope (YAGNI)

- Peak-disk reduction (skip-writing intermediates mid-run) — cleanup addresses
  retained footprint; peak is a separate concern.
- Per-subcommand `--keep`.
- A `--clean/--drop` inverse flag.
- Migrating the benchmark to the new command.
