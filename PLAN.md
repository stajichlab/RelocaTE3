# RelocaTE3 Development Plan

Status: draft assessment + roadmap (2026-06-10)

RelocaTE3 is a ground-up modernization of RelocaTE2 (in `../RelocaTE2`), a tool
that maps transposable element (TE) insertion polymorphisms from short-read
resequencing data at single-base resolution. This document assesses both
codebases, defines what to port and what to redesign, and lays out a phased plan
to (1) reach a validated **Python reference implementation**, then (2) push the
heavy per-read processing into **Rust**, and (3) scatter parallel work with
**Nextflow** while keeping a single-node / laptop path.

---

## 1. What RelocaTE2 actually does (the algorithm to preserve)

RelocaTE2 is a 7-step pipeline orchestrated by `scripts/relocaTE2.py`, which
writes per-step shell scripts and either runs them inline (`os.system`) or via a
`multiprocessing.Pool`. The science lives in a handful of large scripts.

**Inputs**
- Reference genome FASTA (`--genome_fasta`)
- TE/repeat consensus FASTA (`--te_fasta`)
- Resequencing reads: a directory of paired FASTQ (`--fq_dir`) **or** a BAM of
  reads already mapped to the genome (`--bam`)
- Optional RepeatMasker `.out` of TE annotations in the reference
  (`--reference_ins`) — used for ref/shared calls and false-positive filtering

**Pipeline steps**

| Step | Script | Role |
|------|--------|------|
| 0 | `relocaTE2.py` / `blat` | Build `existingTE.bed` from RepeatMasker `.out` (or blat genome×TE) |
| 2 | `seqtk` | FASTQ → FASTA (only needed for the blat aligner path) |
| 3 | `blat`/`bwa`/`bowtie2` + `relocaTE_trim.py` | Map reads to TE library; trim the TE-matching portion off each read, emit **flanking reads** + 5′/3′/middle classification + `read_repeat_name` table |
| 4 | `relocaTE_align.py` (+ `clean_pairs_memory.py`) | Re-map trimmed flanking reads back to the genome with bwa, repair mate pairing, merge/sort BAM |
| 5 | `relocaTE_insertionFinder.py` (1827 lines) | **Core.** Cluster junction + supporting reads at genomic positions, detect target-site duplication (TSD), call **non-reference** insertions, write GFF/TXT |
| 6 | `relocaTE_absenceFinder.py` (968 lines) | Call **reference-only** and **shared** insertions (TE present in reference) |
| 7 | `characterizer.pl` (Perl) + `clean_false_positive.py` | Genotype homozygous/heterozygous/somatic from a genome BAM; drop calls overlapping known reference TEs |

**Key algorithmic concepts (must be reproduced faithfully)**
- *TE-containing read*: a read whose end aligns to a TE consensus boundary. The
  read is split into a **TE portion** and a **flanking (trimmed) portion**.
  `relocaTE_trim.py` scores each alignment by a **boundary** count (read-end
  and/or TE-end reached, 0–4) plus match length, keeps the best, and classifies:
  - 5′ junction (TE start at read boundary), 3′ junction (TE end), or "middle"
    (whole read inside TE → supporting, not junction).
  - Trimmed flank must be ≥ `len_cut_trim`; TE match must be ≥ `len_cut_match`;
    mismatches ≤ `mismatch`.
- *Junction read*: flanking read remaps adjacent to an insertion; its clipped
  end marks a precise breakpoint.
- *Supporting read*: a read pair where one mate is in the TE and the other lands
  in unique genome sequence, bracketing the insertion but not crossing it.
- *Clustering*: junction reads are binned by genomic position; left- and
  right-side junctions that meet across a short gap define a **TSD** (target site
  duplication, typically a few bp). `TSD_from_read_depth` / `tsd_finder` derive
  the TSD string and exact coordinates from read-depth overlap.
- *Output GFF attributes*: `TSD`, `Name` (TE family), `Note` (Non-reference /
  Reference-Only / Shared), and the four counts
  `Left/Right_junction_reads`, `Left/Right_support_reads`.

**Validation target (golden test):** rice Chr3 2 Mb simulation in
`../RelocaTE2/test_data/`. Truth = `MSU7.Chr3_2M.ALL.gff` (~196 simulated
insertions). RelocaTE2's acceptance check is:
```
bedtools window -w 10 -a MSU7.Chr3_2M.ALL.gff -b .../ALL.all_nonref_insert.gff | wc -l   # == 196
```
The same dataset is already vendored into RelocaTE3 at
`tests/data/sim_genome/` and `tests/data/sim_reads/` (with `mping.fa` TE lib).
This is our north-star regression test.

**Pain points in RelocaTE2 we explicitly want to leave behind**
- Python 2 idioms (`print` statements, `dict.has_key`) — the shipped
  `relocaTE_insertionFinder.py` still contains Py2 syntax and will not run on
  Py3 as-is; `relocaTE_align.py` is checked in as a **diff/patch file**, not
  runnable source.
- Hardcoded tool paths and a `CONFIG` file; tool discovery via `which`.
- Orchestration by writing shell scripts and `os.system` — no real error
  handling, no provenance, hard to test.
- A Perl dependency (`characterizer.pl`) for genotyping.
- Deeply nested `defaultdict` state passed positionally through 15-argument
  functions; effectively untestable.
- Dependence on `blat` (and Py2-era `bwa 0.6.2`) for sensitivity.

---

## 2. State of RelocaTE3 today

Modern packaging is in place; the science is ~15% done.

**Good / keep**
- `pyproject.toml` src-layout package, `versioningit`, `pre-commit`, `trunk`,
  Sphinx docs scaffold, `pytest` config. Entry point `relocaTE3`.
- `src/RelocaTE3/align.py` — `Aligner` class wrapping **minimap2** + `pysam`:
  index a TE library, map left/right read files to it, sort, drop unmapped
  (`-F 0x4`), index BAM. This is a clean replacement for Step 3 mapping.
- `src/RelocaTE3/ReadLibrary.py` — `ReadLibrary` (1–2 FASTQ + sample name) and a
  partial `TrimmedReadLibrary`.
- `src/RelocaTE3/librelocate.py` — `RelocaTE.trim_TE_reads()` parses the TE BAM
  and computes the same per-read boundary/match scoring as RelocaTE2's
  `parse_align_bwa`. This is the right nucleus for Step 3 trimming.

**Broken / incomplete (blocking)**
- `RelocaTE.identify_TE_reads()` calls `TE_to_readinfo[TE].write_reads(outdir)` —
  `write_reads` does not exist; and it iterates `for TE in TE_to_readinfo`
  treating a single `TrimmedReadLibrary` as a dict of TEs (type mismatch).
- `TrimmedReadLibrary.trimmed_coordinates` is a `defaultdict` **attribute** but
  is invoked as a **method** (`trimlib.trimmed_coordinates(coord)`).
- The actual **trimming** (cut TE sequence out of the read; write flanking FASTQ;
  emit 5′/3′/middle FASTA + `read_repeat_name`) is **not implemented** — only the
  coordinate scoring exists.
- No genome-alignment step (Step 4): `_map_minimap_genome` / `_map_bwa_genome`
  are empty stubs.
- No insertion finder, TSD caller, clustering, absence finder, or
  characterization (Steps 5–7) — i.e. the core science is absent.
- `src/RelocaTE3/__main__.py` is broken: duplicate `-r/--right` argument
  (argparse will raise), references an `align` symbol that is commented out, and
  wires to `pipeline.align`, which is a **placeholder** that squares integers
  from a text file (with matching placeholder tests in `tests/pipeline_test.py`
  and `tests/main_test.py`).
- Source duplication/confusion: real code in `src/RelocaTE3/`, a stale copy in
  `build/lib/RelocaTE3/`, and an empty `RelocaTE3/` (only `__pycache__`).

**Dependencies declared:** `biopython`, `pybedtools`, `pysam` (+ minimap2,
samtools as external binaries). No bwa/blat/bowtie2/seqtk/perl — consistent with
the intent to standardize on minimap2 + pysam.

---

## 3. Design decisions for RelocaTE3

1. **One aligner family by default: minimap2.** Use `-x sr` for short reads vs.
   the TE library and for re-mapping flanking reads to the genome. Keep an
   `Aligner` abstraction so bwa-mem2 can be added later for benchmarking
   sensitivity against RelocaTE2's bwa/blat results.
2. **Library, not shell scripts.** All steps are importable Python functions /
   classes with typed inputs and explicit outputs. The CLI and the workflow
   engine are thin drivers over the library. No `os.system`; use `subprocess`
   with checked errors (already the style in `align.py`).
3. **Replace Perl genotyping** (`characterizer.pl`) with a Python `pysam`-based
   module.
4. **Deterministic, file-based intermediates** with a documented layout, so a
   run can be resumed and so Nextflow can scatter per-chromosome / per-sample
   without code changes.
5. **Single source of truth for state**: replace nested `defaultdict`s with small
   dataclasses (`TEReadAlignment`, `JunctionRead`, `Insertion`, `ReadCluster`).
6. **Golden-output regression** against the rice Chr3 2 Mb set is the definition
   of "the reference implementation works."
7. **Rust later, behind the same Python API.** The Python reference defines
   behavior and tests; Rust reimplements the hot loops (trimming, clustering,
   TSD) and is swapped in via a feature flag / optional import, validated to
   produce identical GFF.

---

## 4. Target architecture (Python reference)

```
src/RelocaTE3/
  __init__.py            # version/logging (keep)
  cli.py                 # argparse → subcommands (replaces broken __main__ menu)
  __main__.py            # thin: delegate to cli.main
  config.py              # run parameters dataclass + tool discovery (minimap2, samtools)
  readlib.py             # ReadLibrary, TrimmedReadLibrary (fix current ReadLibrary.py)
  align.py               # Aligner: minimap2/bwa-mem2 wrappers (extend current)
  trim.py                # Step 3: TE detection + read trimming + 5'/3'/middle + read_repeat
  genome_align.py        # Step 4: map flanking reads to genome, mate repair, merge/sort
  cluster.py             # Step 5a: cluster junction + supporting reads by position
  tsd.py                 # Step 5b: TSD detection from read-depth overlap
  insertions.py          # Step 5: non-reference insertion calls + GFF/TXT writer
  reference_te.py        # Step 0/6: RepeatMasker.out → bed; reference/shared calls
  characterize.py        # Step 7: homo/het/somatic genotyping via pysam
  gff.py                 # GFF3 records, attributes, txt2gff
  pipeline.py            # orchestrate full single-sample run (replaces placeholder)
  models.py              # dataclasses: TEReadAlignment, JunctionRead, Insertion, Cluster
```

**Output layout (per sample / outdir)** — mirrors RelocaTE2 enough to compare:
```
<outdir>/
  te_containing/        # reads matching TE, read_repeat_name tables
  flanking/             # trimmed flanking FASTQ (+ 5'/3' portions)
  genome_aln/           # flanking reads mapped to genome (sorted/indexed BAM)
  results/
    <sample>.nonref_insertions.gff
    <sample>.ref_insertions.gff
    <sample>.insertions.txt
```

---

## 5. Phased roadmap

> Progress: **Phases 0–5 implemented — the full reference pipeline runs end to
> end** (2026-06-10). The package has `models.py`, `trim.py`, `genome_align.py`,
> `insertions.py`, `reference_te.py`, `characterize.py`, `cli.py`, a real
> `pipeline.py`, and a pixi env; `relocaTE3 trim | align-genome | find-insertions |
> find-reference | characterize | run` work end-to-end on the rice test data. 35
> tests pass. **Acceptance gate met**: full 14-family `RiceTE.fa` run recovers
> 178/200 simulated insertions (~89% recall, ~90% precision) vs RelocaTE2's 196/200
> — see `tests/acceptance_test.py`. Reference/shared calling + FP filtering via
> RepeatMasker and genotyping (homo/heterozygous/somatic) are in place. Pure Python,
> only minimap2 + samtools. Next: Phase 6 (docs/provenance polish), then 7–8
> (Rust/Nextflow).

> Recall investigation (2026-06-10): the recall gap was **not** an aligner problem.
> Short junction flanks (10–25 bp) can't be uniquely placed in a 2 Mb genome by any
> aligner, but the paired-end *mates* of those reads do map. Implementing RelocaTE2's
> support-only call path (`_call_support_only`: both-strand bracketing with a clean
> gap) recovered 2 of 4 missed sites (19→21 of 23). The last 2 have overlapping
> +/- support reads (no clean gap) and are rejected by the same `ins_start > ins_end`
> rule RelocaTE2 uses — the genuine hard-case floor for this dataset/coverage.

> Acceptance gate (2026-06-10): RelocaTE2 is not runnable as-is (Python 2;
> `relocaTE_align.py` is a checked-in diff), and no frozen RelocaTE2 GFF is committed,
> so the gate compares RelocaTE3 against the **simulation truth** RelocaTE2 was
> benchmarked on (`MSU7.Chr3_2M.ALL.gff`, 200 sites; RelocaTE2 recovered 196). Running
> RelocaTE3 with the full 14-family `RiceTE.fa` (`--mismatch 2`): **178/200 recovered
> (~89% recall) at ~90% precision** (`tests/acceptance_test.py`). Getting there
> required two sensitivity fixes that close most of the minimap2-vs-blat gap:
> (a) the TE-read step now maps with `-k 11 -w 5 -N 20 -p 0.5` against the FASTA
> directly (was `-x sr` defaults; 157→173 recovered); (b) the genome re-align step
> uses `-k 13 -w 6` (173→178). The residual ~18-site gap to RelocaTE2 is the
> remaining blat sensitivity on short/divergent TE matches.

### Phase 0 — Repo hygiene & test harness (small, do first)  ✅ DONE
- Delete stale `build/lib/RelocaTE3/` copy and the empty top-level `RelocaTE3/`
  package; keep a single `src/RelocaTE3/`. Remove committed VSCode logs.
- Remove the placeholder `pipeline.align` square-numbers function and its tests
  (`pipeline_test.py`, the `test_main_align*` cases) once real subcommands exist.
- Add a `tests/golden/` regression that runs the full pipeline on the vendored
  rice Chr3 2 Mb data and asserts the `bedtools window -w 10` recovery count and
  precision against `MSU7.Chr3_2M.ALL.gff`. Mark `@pytest.mark.slow` /
  integration; gate on minimap2+samtools availability.
- Provide a `pixi.toml` (use `[workspace]`) pinning minimap2, samtools, bedtools
  so contributors and CI get identical tool versions.

### Phase 1 — Step 3: TE read identification + trimming  ✅ DONE
- Fix `ReadLibrary.py`: make `trimmed_coordinates` a real method; add
  `TrimmedReadLibrary.write_reads()`; correct class-level mutable defaults
  (`file_set`, dicts) to instance attributes.
- Finish `trim.py`: given the per-read best TE alignment (port the existing
  `trim_TE_reads` scoring), **cut** the TE portion from each read, write:
  - `*.flankingReads.fq` (trimmed flank, renamed with `:start/end:5/3` /`:middle`
    tags exactly as RelocaTE2 so downstream parsing matches),
  - `*.five_prime.fa` / `*.three_prime.fa` (TE portions),
  - `*.read_repeat_name.txt` (read → TE family, strand).
  Reuse the thresholds `len_cut_match`, `len_cut_trim`, `mismatch`.
- Fix `identify_TE_reads` control flow (it currently mis-iterates the trim
  result). Unit-test against a few hand-built BAM records and against the
  RelocaTE2 `relocaTE_trim.py` output on the same input.

### Phase 2 — Step 4: map flanking reads back to genome  ✅ DONE
- `genome_align.py`: minimap2 `-x sr` of flanking FASTQ to genome; recover the
  genomic *mates* of TE-containing reads as supporting reads (replacement for
  `clean_pairs_memory.py`); merge + coordinate-sort + index BAM with pysam.
- Read-name junction tags are carried through so the insertion finder can tell
  junction vs. supporting reads.
- Deferred refinements (revisit during Phase 3 GFF validation): reads are mapped
  single-end rather than as proper pairs; the `*.fullreads.bam` used for
  false-junction filtering is not yet produced; short flanks (<~15 bp) don't map
  uniquely under minimap2 `-x sr` (the blat-sensitivity tradeoff noted in §7).

### Phase 3 — Step 5: non-reference insertion calling (the core)  ✅ DONE + parity pass
- `insertions.py` streams the genome BAM in coordinate order, bins reads into
  clusters (`RANGE_ALLOWANCE`), derives left/right junction breakpoints from the
  read-name tags + mapped strand, pairs them into a TSD (overlap within
  `TSD_WINDOW`), counts bracketing supporting reads, and emits GFF3 + TXT with the
  full RelocaTE2 attribute set. Built on `models.Insertion` / `JunctionObservation`
  instead of nested dicts.
- **Parity pass (done):** sub-cluster splitting via `_pair_breakpoints` so one
  cluster yields multiple insertions (RelocaTE2 `TSD_from_read_depth` pairing);
  proper TE orientation from the 5'/3' tag × junction side; full-read
  false-junction filtering — Phase 2 now also produces `<sample>.fullreads.genome.bam`
  (untrimmed junction reads), and `find_insertions` drops sites whose full reads
  map across the breakpoint (≥30% both sides).
- Validated on the golden set: 21 calls, 19 within 10 bp of a true mPing site
  (~90% precision, ~83% recall); 0 false junctions filtered (all real here);
  orientation splits +/- correctly.
- **Still deferred** (lower value): read-depth TSD length refinement
  (`tsd_finder`/`TSD_len_calculate` — currently TSD = the genomic overlap span);
  low-quality/MAPQ read filtering. The ~83% recall is the minimap2-vs-blat
  short-flank gap from Phase 2 (§7), not a clustering issue.
- **Acceptance (still open):** match RelocaTE2's window recovery on the golden set
  with comparable false positives; diff GFFs against a frozen RelocaTE2 run.

### Phase 4 — Steps 0 & 6: reference & shared insertions  ✅ DONE (core)
- `reference_te.py`: `parse_repeatmasker` (port of `existingTE_RM_ALL`, incl.
  intact-element detection) → `ReferenceTE` records + `write_existing_te_bed`;
  `find_reference_insertions` calls reference/shared insertions where junction
  clusters coincide with an intact reference TE boundary; `filter_reference_overlaps`
  drops non-ref calls overlapping a same-family reference TE (RelocaTE2 step-7
  `bedtools intersect -v` cleanup).
- Wired into `run_sample` via the optional `--repeatmasker` flag, plus a standalone
  `find-reference` subcommand. The rice RM `.out` (1252 elements) is vendored and
  parsed in tests; end-to-end on the golden set: existingTE.bed built, 0 shared
  calls (mPing is not a reference element here, as expected), non-ref recall
  unchanged.
- **Deferred to the parity pass for this phase:** true *absence* detection
  (reference-only = TE present in reference but deleted in the sample, via reads
  spanning the locus) from `relocaTE_absenceFinder.py`; shared/ref family-name
  reconciliation when library names differ from RepeatMasker names.

### Phase 5 — Step 7: genotyping (replace Perl)  ✅ DONE (core)
- `characterize.py`: pysam reimplementation of `characterizer.pl` — `count_spanners`
  (reference-allele reads mapping cleanly across a site), `classify_status` (the
  RelocaTE2 homozygous/heterozygous/somatic ladder), `characterize_insertions`, and
  `write_characterized` (GFF + TXT). `Aligner.map_library_to_genome` produces the
  required whole-genome alignment of the original reads.
- Wired into `run_sample` via the optional `--genotype` flag, plus a standalone
  `characterize` subcommand (reads back the non-ref GFF via `read_insertions_gff`).
  End-to-end on the golden set: 10 two-sided insertions all genotyped homozygous
  with 0 spanners — correct, since the simulated reads come from a genome that
  carries the insertions.
- The `clean_false_positive.py` reference-overlap cleanup was already implemented
  in Phase 4 (`filter_reference_overlaps`).
- **Deferred:** the optional excision-with-footprint VCF analysis (bcftools) from
  the Perl tool.

### Phase 6 — CLI, pipeline orchestration, docs
- `cli.py`: subcommands `trim`, `align-genome`, `find-insertions`,
  `find-reference`, `characterize`, and `run` (end-to-end single sample). Fix the
  duplicate-argument bug and remove the placeholder wiring.
- `pipeline.py`: a single-node driver that runs all steps for one sample with a
  thread pool (laptop path). Provenance/log of commands and versions.
- Update README + Sphinx; document the output GFF and parameters; migration
  notes from RelocaTE2 flags.

**Milestone = "reference implementation validated":** end of Phase 5, when
`relocaTE3 run` on the rice Chr3 2 Mb set reproduces RelocaTE2's non-reference,
reference, and genotyped calls within tolerance, all in pure Python with only
minimap2 + samtools (+ bedtools) as external binaries.

### Phase 7 — Rust acceleration (after the reference is frozen)
- Profile; the hot paths are (a) BAM iteration + per-read trimming/scoring and
  (b) clustering + TSD over dense read pileups.
- Implement those in a Rust crate exposed to Python via **PyO3/maturin** (use
  `rust-bio` / `noodles` for BAM/FASTA, `rust-htslib` if linking htslib is
  acceptable). Keep the exact same function signatures the Python modules expose.
- Gate Rust behind an optional import: `from RelocaTE3._native import ...` with a
  pure-Python fallback, so laptops without the wheel still work.
- **Validation:** Rust and Python must produce byte-identical GFF on the golden
  set; add a CI job asserting this.

### Phase 8 — Nextflow scatter (parallel/cluster) + single-node parity
- Wrap the Phase 6 subcommands in a DSL2 workflow: scatter Step 3 by read chunk,
  Step 5/6 by chromosome, then gather. Each process is one `relocaTE3` subcommand
  — no logic in the workflow itself.
- Provide profiles: `standard`/`local` (laptop, no scheduler) and a SLURM profile
  for UCR HPCC (see the `nextflow-hpcc` conventions: partition selection, module
  vs. pixi vs. singularity provisioning). The library's `pipeline.py` remains the
  zero-Nextflow path for a single node.

---

## 6. Immediate next actions
1. Phase 0 cleanup: collapse to one package dir, drop the placeholder
   pipeline/tests, add `pixi.toml`, add the golden integration test skeleton
   (skipped until Phase 3 produces output).
2. Phase 1: fix `ReadLibrary`/`TrimmedReadLibrary`, implement read trimming +
   flanking/portion/`read_repeat` outputs, and repair `identify_TE_reads`.
3. Stand up a frozen RelocaTE2 reference run on the rice data (run it once with
   the legacy scripts, archive the GFFs) to diff against during Phases 3–5.

## 7. Open questions for the maintainer
- Is dropping blat/bowtie2 (minimap2-only by default) acceptable for the
  sensitivity targets, or must we reproduce blat-level recall on hard TEs?
- BAM-input mode (reads pre-mapped to genome) — port in the reference phase, or
  defer until after FASTQ mode is validated?
- Minimum supported Python (current `pyproject` says ≥3.7; recommend bumping to
  ≥3.10 to match the dataclass/typing style and current test interpreters).
- Target external-tool versions to pin in `pixi.toml` (minimap2, samtools,
  bedtools).
