# Release readiness: one-command pipeline, reference insertions, reproducibility, read-name handling

Four fixes for a stable release, each independently reviewable. The goals were
the two stated for the release: **easy to run**, and **results comparable to
RelocaTE2**. Two of these are silent-wrong-answer bugs, not conveniences.

Reviewing commit-by-commit is easier than the squashed diff — each commit
carries its own measurements.

---

## `372bcf6` — `run-all` and `find-reference` (easy to run)

Two capabilities existed in the library but were unreachable from the installed
CLI, so users had neither a single-command pipeline nor RelocaTE2's
`all_ref_insert` output.

```bash
relocaTE3 run-all -1 reads_1.fq.gz -2 reads_2.fq.gz -T RiceTE.fa -g ref.fa \
  -n HEG4 -o HEG4_out --threads 8 --repeatmasker ref.fa.RepeatMasker.out --genotype
```

`run-all` chains `index-genome → run → align-genome → find-insertions`, plus
`find-reference` when `--repeatmasker` is given and `characterize` when
`--genotype` is.

**Worth reviewing:** it is a thin orchestrator dispatching the *same* `cmd_*`
handlers as the staged subcommands, **not** a wrapper around
`pipeline.run_sample` as the TODO proposed. `run_sample` calls the module-level
`insertions.find_insertions` while the staged CLI uses `InsertionFinder`; the two
are not equivalent, so wrapping it would have shipped a one-command mode
producing different, unvalidated numbers than the benchmark validated. Verified
empirically on the Chr3 2 Mb fixture: `run-all` and the staged workflow produce
identical `all_nonref_insert.txt`.

`find-reference` emits `results/<name>.all_ref_insert.{gff,txt}` + `existingTE.bed`.
Distinct from `annotate-ref`, which only records where reference TEs are in order
to filter novel calls.

Tests cover dispatch order and per-step arguments (so `run-all` cannot silently
drift from the staged path), optional-stage skipping, the between-step guards, and
an end-to-end CLI run on the real 6X fixture.

Also corrects two stale TODO entries verified against the code and benchmark.

## `8d2fcb4` — TE-family assignment is reproducible run-to-run

Two identical runs on identical inputs could label the same call with different
TE families — Chr3:672695 came out `RIRE3` or `mGing`, with byte-identical
positions and read counts. Two independent order-dependent tie-breaks, both now
total orders:

1. `RelocaTE._is_better` returned `False` on an exact tie, so a read's TE match
   was whichever alignment the aligner emitted first — an order that is not
   stable across threads for multi-mapping reads. Replaced by `_match_rank`.
2. `_make_insertion` voted on the cluster family with `max(set(names),
   key=names.count)`. **`set` iteration order over strings varies between
   processes** because CPython randomises string hashing, so a tied vote
   resolved differently on each run. This was the demonstrated cause; extracted
   as `_majority_te_name`, ranking by `(-count, name)`.

Note that (2) is invisible to a single-process test — it only appears across
process boundaries. Which name wins among equals is arbitrary; that it is always
the same one is the point.

**Verified:** five `run-all` runs under `PYTHONHASHSEED=1..5` give byte-identical
`read_repeat_name.txt`, `all_nonref_insert.txt` and `all_ref_insert.txt`.
Accuracy unchanged — against the pre-fix run the output differs by exactly one
line (the tied label): same 158 calls, same positions and counts, zero
`read_repeat` changes. Acceptance thresholds still pass.

## `6c06d74` — docs: correct the reproducibility note

The note added alongside `run-all` described the pre-fix behaviour. Now states
that `results/` tables are the reproducible artifacts, while intermediate BAMs
from multi-threaded aligners may still differ byte-for-byte.

## `b033d57` — derive read mates from the source file, not the read name

The highest-impact fix here. Mate pairing string-matched a trailing `/1`,`/2` on
read names. **Modern Illumina FASTQs carry the mate in a separate field, so the
name is identical in R1 and R2; SRA/ENA dumps often have no marker at all.** On
those inputs `split_mate` returned `mate=''`, `recover_support_mates` skipped
every read, and no junction flank was mate-anchored — **with exit 0 and no
warning.** RelocaTE2 made this configurable via `--mate_1_id`/`--mate_2_id`.

Measured on the Chr3 2 Mb fixture with only the suffixes stripped:

| input | supporting reads | calls with support | recall |
|---|---|---|---|
| with `/1`,`/2` | 3955 | 158/158 | 158/200 |
| without | **17** | **4/154** | 154/200 |
| without, fixed | 3986 | 160/160 | 160/200 |

Which file a read came from settles the mate, so that is now authoritative:

- `canonical_name(raw, mate)` strips any comment field and any existing suffix,
  then applies the side's mate.
- `canonicalize_te_bams` stamps the mate onto the per-side TE-library BAMs for
  **every** backend. Only `bwa` did this before, so minimap2/bowtie2/blat left
  names bare. Single-end libraries are untouched.
- Both FASTQ lookups match on the canonical name rather than the raw one.

**Backward compatible:** on `/1`-suffixed input the output is byte-identical to
the previous code, so the RelocaTE2 parity numbers are unaffected.

One honest caveat: a suffix-less run is still not bit-identical to a suffixed one
(158 vs 160 calls). That is **minimap2, not RelocaTE3** — it seeds its tie-break
among equally-scoring multi-mappings from the read name. Verified directly: the
same read sequence under two names maps to `Bajie` pos=654 vs pos=4421 under an
identical command line. Not fixable here.

---

## Docs

README gains the `run-all` quick start, `run-all`/`find-reference` rows in the
subcommand table, an Outputs section, a read-name-convention note in Inputs, a
reproducibility note, and two new migration-table rows
(`--mate_1_id`/`--mate_2_id` → not needed; whole pipeline in one command →
`run-all`).

## Verification

Full suite passes (1 skip: `blat` not on PATH), including the acceptance gate.
New test files: `tests/te_family_determinism_test.py` (7),
`tests/read_name_conventions_test.py` (10), plus `run-all`/`find-reference`
coverage in `tests/main_test.py` and `tests/acceptance_test.py`.

## Not in this PR

`run-batch` (multi-sample / `--fq_dir`) is stacked on top in
`feat/cli-run-batch`. Still open and non-blocking: `--dry_run`, `--bam` reuse,
and the ~0.02 TSD-exactness gap vs RelocaTE2.
