# Characterize-stage validation comparison — design

**Date:** 2026-06-19
**Scope:** `RelocaTE3/validation/real_rice/` (with one small upstream change in
`src/RelocaTE3/characterize.py` to support CRAM input).
**Status:** Approved; ready for implementation per Section 9 order.

> **Note (2026-07-11):** references below to `__main__.py` as "the registered entry
> point" / two CLI front-ends are historical. The CLI now lives in `cli.py` (entry
> point `RelocaTE3.cli:main`); `__main__.py` is a thin launcher. See `AGENTS.md`.
**Update 2026-06-19:** CRAM support has landed upstream in
`src/RelocaTE3/characterize.py` (`_open_alignment` detects `.cram` by suffix
and passes `reference_filename`; `Characterizer.characterize` takes a
`genome_fasta` arg) and in `src/RelocaTE3/__main__.py::_menu_characterize`
(`-g/--genome-fasta` flag). Step 4 of the implementation order is therefore
already done — the validation harness can call `relocaTE3 characterize` with
CRAM input today.

## 1. Goals & scope

Extend the real-rice validation harness so it compares RelocaTE3's
**characterized** (genotyped) output against RelocaTE2's
`ALL.all_nonref_insert.characTErized.{gff,txt}`, alongside the existing raw
non-reference comparison. The harness must:

1. Actually run RelocaTE3's Step 7 `characterize` on each sample, using the
   pre-existing CRAMs at
   `validation_data/real_rice/aln_input/<sample>/<sample>.sorted.markdup.cram`
   as the reads-to-genome alignment, with a future config switch to instead
   have RelocaTE3 align reads itself.
2. Produce **two independent comparisons**:
   - `report/nonref/` — the existing raw-call comparison (preserved; refactored
     to share helpers).
   - `report/characterized/` — a new comparison over genotyped calls that
     scores position agreement, TSD agreement, and genotype/status agreement.
3. Emit summary stats, per-sample TSVs, and plots (venn for shared/unique,
   plus dedicated plots for TSD and status agreement).

**Out of scope:** changing RelocaTE3's characterize *logic*, changing
`run_validation.sh`'s SLURM model, regenerating `validation_data/`, or
consolidating the two CLI front-ends (`__main__.py` vs `cli.py`).

## 2. Why this is needed

The validation harness currently compares only raw `ALL.all_nonref_insert.gff`
files — it never exercises the Step 7 genotype call. `run_relocate3.sh` even
has a placeholder block that explicitly skips characterize because no
reads-to-genome BAM is produced upstream. As a result:

- The single file RelocaTE3 produces per sample
  (`results/<sample>/results/ALL.mping.all_nonref_insert.txt`) carries no
  status (homozygous / heterozygous / somatic), so we cannot tell whether
  RelocaTE3's Step 7 is faithful to RelocaTE2 even when position calling
  agrees.
- We have a ready-made source of truth: RelocaTE2 already produced
  `ALL.all_nonref_insert.characTErized.{txt,gff}` for the same 10 rice
  samples using the bundled CRAMs.

## 3. Directory & file layout after the change

```
RelocaTE3/validation/real_rice/
├── README.md                       # updated to document two-stage compare
├── config.example.toml             # new keys for char + alignment dir
├── config.toml                     # gitignored (collaborators copy from example)
├── _config.py                      # adds resolve_genome_aln() helper
├── _common.py                      # NEW: shared TSV io, greedy nearest-neighbor
│                                   #      match, filter, venn helpers
│                                   #      (extracted from current compare_calls.py)
├── run_relocate3.sh                # updated: actually runs `characterize` when
│                                   #          [relocate3].characterize = true
├── run_validation.sh               # updated: orchestrates both compare stages
│
├── normalize_relocate2_nonref.py   # = current normalize_relocate2.py renamed
├── normalize_relocate3_nonref.py   # = current normalize_relocate3.py renamed
├── compare_nonref.py               # = current compare_calls.py renamed,
│                                   #   thinned to use _common.py
│
├── normalize_relocate2_char.py     # NEW: parse ALL.all_nonref_insert.characTErized.txt
├── normalize_relocate3_char.py     # NEW: parse <name>.all_nonref_insert.characTErized.txt
├── compare_char.py                 # NEW: position+family match, then TSD and status agreement
│
├── results/<sample>/results/       # RelocaTE3 outputs (existing) — now also includes
│                                   # <name>.all_nonref_insert.characTErized.{txt,gff}
├── logs/                           # unchanged
└── report/
    ├── nonref/                     # existing comparison's outputs land here
    │   ├── relocate2_calls.tsv
    │   ├── relocate3_calls.tsv
    │   ├── summary.{tsv,txt}
    │   ├── matched_calls.tsv
    │   ├── relocate{2,3}_only.tsv
    │   ├── venn_total.png
    │   └── venn/venn_<sample>.png
    └── characterized/              # NEW
        ├── relocate2_calls.tsv
        ├── relocate3_calls.tsv
        ├── summary.{tsv,txt}
        ├── matched_calls.tsv       # adds tsd_match, status_match columns
        ├── relocate{2,3}_only.tsv
        ├── tsd_confusion.tsv       # NEW: contingency table of R2 tsd × R3 tsd
        ├── status_confusion.tsv    # NEW: contingency of R2 status × R3 status
        ├── venn_total.png          # shared/unique by position
        ├── venn/venn_<sample>.png
        ├── tsd_agreement.png       # NEW: bar — % of matched pairs with TSD match per sample
        ├── status_agreement.png    # NEW: stacked bar — status concordance per sample
        └── status_confusion.png    # NEW: heatmap of status_confusion.tsv
```

The rename of the three existing scripts is purely cosmetic — their logic
moves into `_common.py` and a thin wrapper. The README adds a one-paragraph
migration note for any user/script still referencing the old names.

## 4. Config schema additions

Append to `config.example.toml`:

```toml
[paths]
# (existing keys unchanged)

# Per-sample reads-to-genome alignments used by Step 7 characterize.
# Format string: ${sample} is substituted per sample.
# Default points at the bundled CRAMs the legacy RelocaTE2 characterize run used.
genome_aln_dir     = "${paths.validation_root}/aln_input"
genome_aln_pattern = "${sample}/${sample}.sorted.markdup.cram"

[relocate3]
# (existing keys unchanged)
characterize        = true          # was false; flipping this on is the main behavior change
# How to source the BAM/CRAM for characterize:
#   "external"  -> use paths.genome_aln_dir / paths.genome_aln_pattern
#   "relocate3" -> have RelocaTE3 align reads to the genome itself (slower; future work)
characterize_source = "external"

[compare_char]
# Same window as compare (defaults to inherit). Status comparison is reported as a
# separate metric — pairs are matched on position, then status/TSD agreement scored.
position_window     = 100
te_family           = "mPing"
# Optional: normalize R2 status strings before scoring (e.g. drop "/excision_with_footprint" suffix).
# Off by default — we want to see those distinctions.
collapse_excision_suffix = false
```

`_config.py` gets one small addition: a `resolve_genome_aln(cfg, sample)`
helper that substitutes `${sample}` into `genome_aln_pattern` and joins under
`genome_aln_dir`, raising if the file is missing. This is the only logic
change to `_config.py`.

The existing `[compare]` block is unchanged; the new `[compare_char]` block is
parallel so the two stages can have independent windows/filters in the future.

## 5. `run_relocate3.sh` changes

The end-of-script `if [[ "$DO_CHARACTERIZE" == "1" ]]` block is currently a
no-op placeholder. Replace it with a real invocation:

1. Resolve the per-sample alignment via `_config.py::resolve_genome_aln`. If
   `characterize_source == "external"`, the resolver returns the path under
   `genome_aln_dir` and the script `exit 1`s with a clear error if the file
   is missing.
2. Build the expected output path:
   `${SAMPLE_OUTDIR}/results/${TARGET}.${TE_NAME}.all_nonref_insert.characTErized.txt`.
3. Idempotency: skip the step if that file already exists and `--force` was
   not passed (the existing `--force` wipes the whole sample dir, so this is
   automatic).
4. Run:
   ```bash
   relocaTE3 characterize \
     -s "$NONREF_TXT" \
     -b "$GENOME_ALN" \
     -g "$GENOME" \
     -o "${SAMPLE_OUTDIR}/results" \
     --samtools samtools --bcftools bcftools
   ```
   `-g/--genome-fasta` doubles as the CRAM reference (the upstream
   `_open_alignment` helper detects `.cram` by suffix and passes it as
   `reference_filename` to pysam). No `--excision` for now; can be exposed as
   a config flag later. `bcftools` is only needed for excision, so leaving it
   off avoids a hidden dependency if/when we drop it from defaults.
5. The `characterize_source == "relocate3"` branch is a stub: it logs that
   this mode is not yet implemented and exits non-zero. The plan calls out
   the future hook into `pipeline.run_sample(genotype=True)` but doesn't
   build it.

**CRAM input note.** Already handled upstream. `_open_alignment` in
`characterize.py` detects CRAM by suffix and opens with
`reference_filename=<genome_fasta>`; `Characterizer.characterize` takes a
`genome_fasta` arg and `_menu_characterize` exposes `-g/--genome-fasta`. No
upstream change is needed from this harness work.

## 6. Normalize scripts for characterized output

### `normalize_relocate2_char.py`

Parses
`validation_data/real_rice/relocate2_results/<sample>/repeat/results/ALL.all_nonref_insert.characTErized.txt`.
Format:

```
strain  TE  TSD  chromosome.pos  strand  avg_flankers  spanners  status
B_10    mPing  TTA  Chr1:1041521..1041523  +  4  0  homozygous
```

`chromosome.pos` field is `<chrom>:<start>..<end>`. Output schema (a flat TSV
mirroring the nonref one's spirit but with the extra fields):

```
sample, chrom, start, end, strand, te_name, tsd, avg_flankers, spanners, status, source_file
```

Skip the header row; skip blank/comment lines; warn on missing files
(mirroring `normalize_relocate2.py`).

### `normalize_relocate3_char.py`

Parses RelocaTE3's
`<outdir>/<sample>/results/<target>.<te_name>.all_nonref_insert.characTErized.txt`.
Same 8-column schema as RelocaTE2's output (confirmed by inspection; both
files share an identical header line). Emits the same flat TSV.

### Shared columns module

Both characterize-normalize scripts import `CHAR_COLUMNS` from `_common.py`
(same pattern the current `normalize_relocate3.py` uses by importing
`COLUMNS` from `normalize_relocate2.py`).

### Filtering

The characterize step only writes rows with two-sided junction support (see
`characterize.py::_score_sites`), so a `min_junction_reads` filter is
redundant here. The `te_family` filter from `_common.py` is reused
(case-insensitive `mping` ↔ `mPing`).

## 7. `compare_char.py` and the new metrics

`compare_char.py` reuses `_common.py` for: TSV io, the `te_family` filter,
and the greedy nearest-neighbor matcher (same `(chrom, te_name.lower())` key,
same `position_window` rule). It does **not** also match on status — per the
brainstorming decision, status is reported as a separate metric on top of
position-matched pairs.

### Per-sample output (`report/characterized/summary.tsv`)

Columns:

```
sample
relocate2_total
relocate3_total
shared                       # position-matched pairs
relocate2_only
relocate3_only
recall_vs_r2
precision_vs_r2
jaccard
# new metrics — denominators are the matched pairs
tsd_match_n
tsd_match_rate               # tsd_match_n / shared
status_match_n
status_match_rate            # status_match_n / shared
status_match_n_collapsed     # if collapse_excision_suffix=true; else same as status_match_n
```

### Pair-level output (`report/characterized/matched_calls.tsv`)

Same columns as the nonref `matched_calls.tsv`, plus:

```
r2_tsd, r3_tsd, tsd_match
r2_status, r3_status, status_match
r2_avg_flankers, r3_avg_flankers
r2_spanners, r3_spanners
```

### Confusion tables

- `tsd_confusion.tsv` — 2-D contingency of distinct R2 TSDs × distinct R3
  TSDs across all matched pairs.
- `status_confusion.tsv` — same for status. This is the right place to see
  "RelocaTE2 calls X homozygous; RelocaTE3 calls them heterozygous" patterns.

### Plots (matplotlib + matplotlib-venn, both already in `pixi.toml`)

- `venn_total.png` and `venn/venn_<sample>.png` — shared vs. unique by
  position (same as nonref).
- `tsd_agreement.png` — per-sample bar chart of `tsd_match_rate`, with the
  total across samples as a final bar.
- `status_agreement.png` — per-sample stacked bar showing matched-pair counts
  split into "status match" vs. each kind of disagreement.
- `status_confusion.png` — heatmap rendering of `status_confusion.tsv`.

All plot generation guarded by an import check, same pattern as `_load_venn`
— if matplotlib is missing, the TSVs are still written.

### Upstream `characterize.py` change — done

CRAM support has already landed: `_open_alignment` detects `.cram` by suffix
and opens with `reference_filename`; `Characterizer.characterize` accepts a
`genome_fasta` arg; `_menu_characterize` exposes `-g/--genome-fasta`. A unit
test fixture covering CRAM input is still worth adding (see Section 9), but
no further upstream code changes are required by this plan.

`cli.py` is not touched by this plan — the two-front-ends consolidation is a
separate task per `AGENTS.md`.

## 8. Orchestrator (`run_validation.sh`) changes

The script keeps its current flags (`--local`, `--slurm`, `--skip-run`,
`--force`, `--config`, positional samples). Changes:

1. **After per-sample RelocaTE3 runs finish**, the existing single
   normalize+compare block is replaced with two parallel blocks:

   ```bash
   # Stage A: raw non-reference calls (unchanged behavior)
   ( cd "$SCRIPT_DIR" && python3 normalize_relocate2_nonref.py --config "$CONFIG" "${SAMPLE_FLAGS[@]}" )
   ( cd "$SCRIPT_DIR" && python3 normalize_relocate3_nonref.py --config "$CONFIG" "${SAMPLE_FLAGS[@]}" )
   ( cd "$SCRIPT_DIR" && python3 compare_nonref.py             --config "$CONFIG" )

   # Stage B: characterized (genotyped) calls
   ( cd "$SCRIPT_DIR" && python3 normalize_relocate2_char.py --config "$CONFIG" "${SAMPLE_FLAGS[@]}" )
   ( cd "$SCRIPT_DIR" && python3 normalize_relocate3_char.py --config "$CONFIG" "${SAMPLE_FLAGS[@]}" )
   ( cd "$SCRIPT_DIR" && python3 compare_char.py             --config "$CONFIG" )
   ```

2. **New flags** for selective re-runs (mirroring `--skip-run`):
   - `--skip-nonref` — skip Stage A.
   - `--skip-char` — skip Stage B (useful while iterating on Stage A).
   - Default: both stages run.

3. **Existing skip check** (currently looks for the raw GFF) is updated to
   require the **characterized** TXT when `[relocate3].characterize = true`,
   so a partial sample (raw done, characterize never ran) is re-attempted
   instead of silently skipped:

   ```bash
   if [[ "$DO_CHAR" == "1" ]]; then
     DONE_MARKER="$OUTROOT/$s/results/${TARGET}.${TE_NAME}.all_nonref_insert.characTErized.txt"
   else
     DONE_MARKER="$OUTROOT/$s/results/${s}.all_nonref_insert.gff"
   fi
   ```

4. **Console output** at end-of-run prints both `report/nonref/summary.txt`
   and `report/characterized/summary.txt`.

## 9. Testing & verification plan

This is a validation harness rather than a library, so most checks are
end-to-end against the bundled rice data. Small unit tests where they add
real coverage.

### Unit tests (new, under `RelocaTE3/tests/`)

1. **`characterize_test.py` — CRAM input.** Test that calls
   `Characterizer.characterize` against a tiny CRAM fixture + matching
   reference FASTA in `tests/data/`, asserting `(txt_path, gff_path)` are
   produced and contain the expected status string. Catches the
   `reference_filename` plumbing.
2. **`validation/test_normalize_char.py`.** Parses a fixture RelocaTE2
   characterized TXT (8 columns) and asserts the normalized rows; same for a
   fixture RelocaTE3 characterized TXT. Confirms both produce identical
   schemas.
3. **`validation/test_compare_char.py`.** Two synthetic TSV inputs covering:
   an exact position+status match (counts to `status_match_n`), a position
   match with status disagreement (counts to `shared` but not
   `status_match_n`), an R2-only call, an R3-only call, and a TSD
   disagreement. Asserts the summary numbers and the contents of
   `tsd_confusion.tsv` / `status_confusion.tsv`.

### Smoke test (manual)

`pixi run validate-rice --local B_10` — single sample on the login node.
Confirms:

- `results/B_10/results/ALL.mping.all_nonref_insert.characTErized.{txt,gff}`
  is produced.
- `report/nonref/summary.txt` and `report/characterized/summary.txt` both
  exist.
- The characterized summary shows non-zero `tsd_match_n` and
  `status_match_n` (sanity: RelocaTE2's TSD is overwhelmingly `TTA` for
  mPing, so TSD match rate should be very high; status agreement on B_10 is
  the headline number we care about).

### Full run

`pixi run validate-rice` (SLURM array, 10 samples). Acceptance bar to call
the plan "done":

- Pipeline completes without manual intervention.
- `report/characterized/summary.tsv` populated for all 10 samples.
- `report/characterized/{tsd,status}_confusion.tsv` and the three new PNGs
  are produced.
- Spot-check: pick one matched pair where R2 says `homozygous` and R3
  disagrees; confirm the `avg_flankers`/`spanners` numbers are visible in
  `matched_calls.tsv` so the disagreement can be triaged.

### Regression bar

The existing acceptance test (`tests/acceptance_test.py`) and the nonref
comparison numbers must not change. The Stage A outputs in `report/nonref/`
should be byte-identical to a baseline captured before the refactor (modulo
`source_file` paths that move with the directory). Baseline captured as
step 1 of Section 10.

## 10. Implementation order

A linear order so each step leaves the harness in a runnable state.

1. **Capture baseline.** Run `pixi run validate-rice --local B_10` against
   the current code; copy `report/` somewhere outside the tree. Reference
   for the no-regression check at the end.
2. **Extract `_common.py`.** Move TSV io, `_filter`, `_match_sample`,
   `_write_pairs`, `_write_rows`, `_load_venn`, `_draw_one_venn` out of
   `compare_calls.py`. Confirm `pixi run validate-rice --local B_10` still
   produces identical `report/` outputs.
3. **Rename trio.** `compare_calls.py` → `compare_nonref.py`;
   `normalize_relocate{2,3}.py` → `normalize_relocate{2,3}_nonref.py`.
   Update `run_validation.sh` and the README. Re-run smoke test; confirm
   `report/nonref/` matches the baseline (move output dir at this step).
4. **`Characterizer` CRAM support.** *Done upstream* — `_open_alignment`
   handles CRAM via `reference_filename`, `Characterizer.characterize` takes
   `genome_fasta`, and `_menu_characterize` exposes `-g/--genome-fasta`. Only
   remaining sub-task is the tiny CRAM fixture + unit test from Section 9.
5. **`_config.py::resolve_genome_aln` + config knobs.** Add
   `genome_aln_dir`, `genome_aln_pattern`, `characterize_source`, and the
   `[compare_char]` section to `config.example.toml`. Update `_config.py`
   minimally.
6. **`run_relocate3.sh` characterize block.** Replace the placeholder with a
   real `relocaTE3 characterize` invocation that reads the CRAM path from
   the config helper. Smoke-test on B_10: confirm `*.characTErized.{txt,gff}`
   are produced.
7. **`normalize_relocate{2,3}_char.py`.** Write the two scripts +
   fixture-based unit tests.
8. **`compare_char.py`.** Write the comparison, the two confusion tables,
   and the summary; fixture-based unit tests for the metrics.
9. **Plots.** Add `tsd_agreement.png`, `status_agreement.png`,
   `status_confusion.png` rendering — guarded by the existing matplotlib
   import check so missing deps still emit TSVs.
10. **Wire into `run_validation.sh`.** Add the two stage blocks, the
    `--skip-nonref` / `--skip-char` flags, and the smarter `DONE_MARKER`.
    README gets a section explaining the two-stage output layout.
11. **Full validation.** `pixi run validate-rice` on all 10 samples. Verify
    against the Section 9 acceptance bar; tweak only if a metric is
    obviously broken (don't touch RelocaTE3's characterize logic in this
    plan).

Each step is a logical commit. Steps 1–3 are pure refactor; steps 4–6 unlock
the characterize step; steps 7–10 add the new comparison; step 11 is
verification.

## 11. Risks and known unknowns

- **CRAM with a different reference.** `pysam` needs the exact FASTA the
  CRAM was compressed against. The bundled CRAMs were aligned against
  `validation_data/real_rice/ref_genome/MSU_r7.fa`, which is the same FASTA
  the config already uses for `paths.genome`. If a collaborator overrides
  `paths.genome` to a renamed/copied FASTA, the CRAM read will fail loudly
  — acceptable, but worth a one-line README note.
- **Spanner-count divergence from RelocaTE2.** RelocaTE2 used BWA-aligned
  CRAMs; `characterize.py`'s spanner-counting logic is faithful to
  `characterizer.pl`, and it counts spanners off those same BWA alignments,
  so the *spanner counts themselves should match*. If they don't, the
  divergence is in our spanner-counting code (real bug to fix) rather than
  upstream aligner choice. The matched-pairs TSV makes this directly
  inspectable.
- **`relocate3` characterize-source mode is a stub.** Section 5 leaves it
  as `exit 1`. Real future work, intentionally not in this plan.
- **CLI front-end divergence.** The CRAM `--reference` flag goes on
  `__main__.py::_menu_characterize` since that's the registered entry
  point. `cli.py` is not touched by this plan; consolidation is separate.
