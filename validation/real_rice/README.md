# Real-rice validation pipeline

Reproducible regression test that runs the current RelocaTE3 on the bundled
rice resequencing dataset and compares calls against a frozen RelocaTE2 run
of the same samples.

The large inputs (FASTQs, CRAMs, genome, RepeatMasker output, legacy GFFs)
live outside this repo under `../../../validation_data/real_rice/` and are
**not** tracked in git. Edit `config.toml` to point at them.

## Layout

```
RelocaTE3/validation/real_rice/
├── README.md                  # this file
├── config.example.toml        # template — copy to config.toml and edit
├── _config.py                 # shared TOML loader (interpolation + path resolution)
├── run_relocate3.sh           # SLURM-aware per-sample driver for RelocaTE3
├── normalize_relocate2.py     # legacy GFFs   -> report/relocate2_calls.tsv
├── normalize_relocate3.py     # RelocaTE3 GFFs -> report/relocate3_calls.tsv
├── compare_calls.py           # diff + stats + venns -> report/summary.{tsv,txt}, report/venn_total.png, report/venn/venn_<sample>.png
├── results/                   # created at runtime: per-sample RelocaTE3 outputs
├── logs/                      # created at runtime: SLURM + python logs
└── report/                    # created at runtime: normalized TSVs + summary + venn
```

## Quick start

From the **`RelocaTE3/` working directory**, one command runs the whole pipeline:

```bash
pixi run validate-rice
```

Equivalent if you don't use pixi:

```bash
bash validation/real_rice/run_validation.sh
```

What it does:

1. Creates `validation/real_rice/config.toml` from `config.example.toml` on
   first run (collaborators on the same HPCC usually don't need to edit it).
2. Submits a SLURM array (one task per sample, partition/cpus/mem/time pulled
   from `[slurm]` in the config) with `sbatch --wait`, so the script blocks
   until every sample's RelocaTE3 run finishes. Per-sample outputs land in
   `validation/real_rice/results/<sample>/`.
3. Runs `normalize_relocate2.py` and `normalize_relocate3.py` to flatten both
   call sets into `validation/real_rice/report/relocate{2,3}_calls.tsv`.
4. Runs `compare_calls.py` to compute per-sample shared / R2-only / R3-only
   counts plus recall, precision, and Jaccard, and to draw the Venn.

Useful flags (pass through the pixi task or directly to the shell script):

```bash
pixi run validate-rice --local B_10            # one sample, no SLURM (fast smoke test)
pixi run validate-rice B_10 B_42               # specific samples via SLURM (one job each)
pixi run validate-rice --local                 # all samples, sequentially in this shell
pixi run validate-rice --skip-run              # only re-normalize + re-compare
pixi run validate-rice --force B_10            # wipe results/B_10/ before re-running
pixi run validate-rice --config /path/to/custom.toml
```

Positional sample names override the config's sample CSV; when given, the
report only reflects those samples (so the Venn / Jaccard numbers are scoped
to what you actually ran).

Outputs land under `validation/real_rice/report/`:

- `summary.tsv` — per-sample shared / R2-only / R3-only counts plus
  recall, precision, and Jaccard index.
- `summary.txt` — overall plain-text summary (printed to stdout too).
- `matched_calls.tsv` — every RelocaTE2 call paired with its nearest
  RelocaTE3 match (or blank if none within `compare.position_window`).
- `relocate2_only.tsv`, `relocate3_only.tsv` — caller-specific calls.
- `venn_total.png` — overall Venn diagram across every sample compared.
- `venn/venn_<sample>.png` — one Venn per sample, labeled with the sample
  name. Both venn outputs require `matplotlib_venn` (already declared as a
  pixi dep); the script silently skips them if the import fails.

## For collaborators on this HPCC

The bundled validation dataset is owned by **nmath020** and lives at
`/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/validation_data/real_rice`.
If you have read access to that path (group membership / file permissions),
the defaults work as-is — `pixi run validate-rice --local B_10` will pick up
the dataset from there and write its own per-clone outputs under
`validation/real_rice/{results,report,logs}/`.

If the dataset has moved or your copy lives elsewhere, override the root
**without editing `config.toml`**:

```bash
export RELOCATE3_VALIDATION_ROOT=/path/to/your/validation_data/real_rice
pixi run validate-rice --local B_10
```

The env var wins over `[paths].validation_root` in the config. Useful for
mirrored copies, scratch staging, or CI runs.

`config.toml`, `results/`, `report/`, and `logs/` are all in this
directory's `.gitignore`, so each collaborator's local state stays out of
the shared repo.

## Manual / per-step usage

You can still drive the individual steps directly (useful for debugging):

```bash
cd RelocaTE3/validation/real_rice
./run_relocate3.sh config.toml B_10          # single sample
python3 normalize_relocate2.py --config config.toml
python3 normalize_relocate3.py --config config.toml
python3 compare_calls.py --config config.toml
```

## Match criterion

Two insertions are "shared" when they:

- come from the same sample,
- sit on the same chromosome,
- name the same `te_name` (e.g. `mPing`),
- and have start positions within `compare.position_window` bp (default 100).

Matching is greedy nearest-neighbor: each RelocaTE2 call is paired with the
closest unmatched RelocaTE3 call within the window. Tune the window in
`config.toml` if you want a stricter or looser definition of agreement.

## Parameter tracking

The defaults in `[relocate3]` mirror the legacy RelocaTE2 run in
`validation_data/real_rice/example_relocate2_pipeline/01_relocate_native_cram.sh`
(`--mismatch 2`, `--min-match 10`, `--min-trimmed 10`). The aligner differs
— RelocaTE2 used BLAT, RelocaTE3 uses minimap2 — so calls will not be
bit-identical, but the recall/Jaccard numbers in `summary.txt` should stay
high. `tsd = "TTA"` matches mPing's canonical TSD; change it if you swap
the TE library. `te_name` is just the label baked into the output filename
`results/<target>.<te_name>.all_nonref_insert.txt` and is compared
case-insensitively against RelocaTE2's `Name=` GFF attribute.

## Pipeline steps

`run_relocate3.sh` chains the underlying RelocaTE3 subcommands per sample
(see `relocaTE3 --help`):

1. `index-genome` — ensure `.fai` + `.mmi` exist (skipped if present).
2. `run` — map reads to the TE library and trim TE sequence (steps 2+3).
3. `align-genome` — align the trimmed flanking FASTQs back to the genome.
4. `find-insertions` — cluster junctions and write
   `<sample>/results/<target>.<te_name>.all_nonref_insert.txt`.

Each step is idempotent: re-runs skip stages whose outputs already exist
(use `--force` from `run_validation.sh` to wipe a sample's directory).

## Dependencies

Activate the project's pixi env (`pixi shell` from `RelocaTE3/`) so
`relocaTE3` is on PATH. The reporting step additionally wants:

```bash
pip install matplotlib matplotlib-venn
```

(They're already declared in `pixi.toml` — `pixi install` from `RelocaTE3/`
covers it.) If those are absent, `compare_calls.py` still writes the
TSV/TXT outputs and just skips the venn PNGs.
