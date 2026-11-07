# Real-rice validation pipeline

Reproducible regression test that runs the current RelocaTE3 on the bundled
rice resequencing dataset and compares calls against a frozen RelocaTE2 run
of the same samples.

The large inputs (FASTQs, CRAMs, genome, RepeatMasker output, legacy GFFs)
live outside this repo under `../../../validation_data/real_rice/` and are
**not** tracked in git. The mPing TE FASTA is tracked at
`tests/data/TE_lib/mping.fa`. Edit `config.toml` to point at the external
validation inputs if needed.

## Layout

```
RelocaTE3/validation/real_rice/
├── README.md                    # this file
├── config.example.toml          # template — copy to config.toml and edit
├── _config.py                   # shared TOML loader (interpolation, path resolution, env overrides)
├── run_validation.sh            # driver: run (pipeline or run-all) + normalize + compare
├── run_pipeline.sh              # pipeline mode: Nextflow head job (sbatch)
├── make_samplesheet.py          # config -> pipeline sample sheet + params
├── run_relocate3.sh             # run-all mode: per-sample relocaTE3 run-all (SLURM array task)
├── normalize_relocate{2,3}_nonref.py  # non-reference GFFs  -> report/nonref/relocate{2,3}_calls.tsv
├── normalize_relocate{2,3}_char.py    # characterized tables -> report/characterized/relocate{2,3}_calls.tsv
├── compare_nonref.py, compare_char.py # diff + stats + venns -> report/{nonref,characterized}/
├── compare_pipeline_to_rc2.sh   # pipeline outdir vs v0.2.0rc2 run-all results, file by file
├── cohort_report.sh             # cohort + MultiQC report over a pipeline outdir
├── results/                     # created at runtime: run-all per-sample outputs
├── logs/                        # created at runtime: SLURM + python logs
└── report/                      # created at runtime: normalized TSVs + summary + venn
```

## Quick start

From the **`RelocaTE3/` working directory**:

```bash
# Default: through the Nextflow pipeline, with a pinned image tag
pixi run validate-rice --container-tag sha-2c9a448
# or without pixi
bash validation/real_rice/run_validation.sh --container-tag sha-2c9a448

# The run-all SLURM array instead (relocaTE3 run-all per sample)
pixi run validate-rice --mode run-all
```

Two modes run RelocaTE3; the report steps after them are the same commands.

**`--mode pipeline` (default).** Runs `run_pipeline.sh`
(`-profile hpcc,apptainer --container_tag <tag>`, settings mapped from this
config by `make_samplesheet.py`) as a head job with `sbatch --wait`; Nextflow
submits one SLURM job per task. Before submitting, `run_pipeline.sh`'s own
`DRY_RUN=1` check reads the config and checks every input.

- `--container-tag TAG` (or env `CONTAINER_TAG`) is required: `sha-<short>`
  (the image of a commit on `main`) or a release version such as `0.3.0rc1`.
  There is no default. `dev` and `pr-<n>` are refused, because CI re-pushes
  them and `-resume` could mix images, unless `ALLOW_MUTABLE_TAG=1`.
- Outputs go to `--pipeline-outdir DIR`, by default
  `results/nextflow-real-rice/<tag>` under the repository root (as
  `run_pipeline.sh`). The pipeline writes the same per-sample layout as run-all
  (`<outdir>/<sample>/results/...`).
- The report reads that outdir: `run_validation.sh` exports
  `RELOCATE3_VALIDATION_OUTDIR=<outdir>`, which `_config.py` uses in place of
  `[paths].relocate3_outdir`. `config.toml` is not edited.
- The report goes to `<[paths].report_dir>/<tag>` (with `--skip-run` and no
  tag: `<[paths].report_dir>/<basename of --pipeline-outdir>`), so reports of
  different images, and the run-all report in `[paths].report_dir`, do not
  overwrite each other.
- `run_validation.sh` refuses to start if `RELOCATE3_VALIDATION_OUTDIR` or
  `RELOCATE3_VALIDATION_REPORT_DIR` is already set in your environment (an
  inherited value would silently redirect the run and the report); use
  `--pipeline-outdir` / `--report-dir`.
- Rerunning the same command resumes (`-resume`). `--force`, `--local` and
  sample names are run-all only; the pipeline always runs the whole sheet.

**`--mode run-all`.** Today's behaviour: a SLURM array (one task per sample,
partition/cpus/mem/time from `[slurm]`) running `run_relocate3.sh` with
`sbatch --wait`; per-sample outputs in `[paths].relocate3_outdir`
(`validation/real_rice/results/<sample>/`).

Both modes then:

1. Normalize both call sets: `normalize_relocate{2,3}_nonref.py` (Stage A,
   non-reference GFFs) and `normalize_relocate{2,3}_char.py` (Stage B,
   characterized tables).
2. Compare: `compare_nonref.py`, `compare_char.py` — per-sample shared /
   R2-only / R3-only counts plus recall, precision, Jaccard, and Venns.

Useful flags:

```bash
# re-report on a finished pipeline run (no tag needed), into its own report dir
bash validation/real_rice/run_validation.sh --skip-run \
    --pipeline-outdir results/nextflow-real-rice/sha-2c9a448 --report-dir /path/to/report
DRY_RUN=1 bash validation/real_rice/run_validation.sh --container-tag sha-2c9a448   # print the commands

pixi run validate-rice --mode run-all --local B_10         # one sample, no SLURM (fast smoke test)
pixi run validate-rice --mode run-all B_10 B_42            # specific samples via SLURM (one job each)
pixi run validate-rice --mode run-all --skip-run           # only re-normalize + re-compare
pixi run validate-rice --mode run-all --force B_10         # wipe results/B_10/ before re-running
pixi run validate-rice --config /path/to/custom.toml ...
```

- `--report-dir DIR` reports into DIR instead of the default
  (`[paths].report_dir` for run-all, `[paths].report_dir/<tag>` for pipeline).
- `--skip-nonref` / `--skip-char` skip Stage A / Stage B.
- `DRY_RUN=1` prints every command (`+ ...`) and creates nothing.
- Positional sample names override the config's sample CSV in run-all mode;
  the report then only covers those samples. In pipeline mode they are
  accepted with `--skip-run` only, to scope the report.

Outputs land under the report dir (`validation/real_rice/report/` for run-all,
`validation/real_rice/report/<tag>/` for pipeline mode, by default),
in `nonref/` and `characterized/`:

- `summary.tsv` — per-sample shared / R2-only / R3-only counts plus
  recall, precision, and Jaccard index.
- `summary.txt` — overall plain-text summary (printed to stdout too).
- `matched_calls.tsv` — every RelocaTE2 call paired with its nearest
  RelocaTE3 match (or blank if none within the position window).
- `relocate2_only.tsv`, `relocate3_only.tsv` — caller-specific calls.
- `venn_total.png`, `venn/venn_<sample>.png` — Venn diagrams; they need
  `matplotlib_venn` (a pixi dep) and are skipped if it does not import.

## Release check

`validation/release/submit_real_rice.sh` runs this harness from the fresh
clone that `validation/release/fresh_clone_check.slurm` made, as two chained
jobs (run, then report with `--skip-run`):

```bash
# pipeline mode (default): the release image tag is required
bash validation/release/submit_real_rice.sh v0.3.0rc1 --container-tag 0.3.0rc1 [--after <release-check job>]
# run-all mode: the clone's locked blat env
bash validation/release/submit_real_rice.sh v0.3.0rc1 --mode run-all [--after <jobid>] [--partition epyc]
```

The image tag is supplied explicitly with `--container-tag`: the release
version (`0.3.0rc1`) once `container.yml` has pushed it for the tag, or the
release commit's `sha-<short>` image before tagging. Mutable tags are refused.
The image cache defaults to this checkout's `results/apptainer-cache`
(`NXF_APPTAINER_CACHEDIR`). `DRY_RUN=1` prints the two `sbatch` commands.

## For collaborators on this HPCC

The bundled validation dataset is owned by **nmath020** and lives at
`/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/validation_data/real_rice`.
If you have read access to that path (group membership / file permissions),
the defaults work as-is — `pixi run validate-rice --mode run-all --local B_10` will pick up
the dataset from there and write its own per-clone outputs under
`validation/real_rice/{results,report,logs}/`.

If the dataset has moved or your copy lives elsewhere, override the root
**without editing `config.toml`**:

```bash
export RELOCATE3_VALIDATION_ROOT=/path/to/your/validation_data/real_rice
pixi run validate-rice --mode run-all --local B_10
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
python3 normalize_relocate2_nonref.py --config config.toml
python3 normalize_relocate3_nonref.py --config config.toml
python3 compare_nonref.py --config config.toml
# on a pipeline outdir instead of [paths].relocate3_outdir:
RELOCATE3_VALIDATION_OUTDIR=/path/to/pipeline_outdir python3 normalize_relocate3_nonref.py --config config.toml
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

`[relocate3]` matches how RelocaTE2 was run on these samples
(`relocate2_results/<sample>/shellscripts/step_6/`, `regex.txt`): BLAT +
`bwa aln`, `--mismatch 2`, `--min-match 10`, `--min-trimmed 10`, TSD `UNK`,
library size 250, no MAPQ gate. Keep `tsd = "UNK"`: a fixed wildcard such as
`"..."` goes through a different writer that RelocaTE2 has no equivalent of.
`te_name` is the label in RelocaTE3's output filenames and is compared
case-insensitively against RelocaTE2's `Name=` GFF attribute.

The non-reference comparison uses the same tier from both tools: RelocaTE2's
`ALL.all_nonref_insert.gff` and RelocaTE3's
`<target>.<te_name>.all_nonref_insert.gff` (both after RelocaTE2's
clean_false_positive.py filters). The unfiltered `.txt` is not comparable.

## Pipeline steps

`run_relocate3.sh` runs the shipped pipeline, `relocaTE3 run-all`, once per
sample with both mates (`-1/-2`) and `--genotype --genotype-bam` pointing at
the external alignment. Running the subcommands one by one, as earlier
versions of this script did, skipped run-all's mate-anchored genome alignment
and RelocaTE2's full-read false-junction filter.

A re-run skips a sample whose `.gff` (and characterized table, if requested)
already exists (use `--force` from `run_validation.sh` to wipe a sample's
directory).

The pipeline (`run_pipeline.sh`) runs the same run-all steps as Nextflow
processes, with the same settings: `make_samplesheet.py --print-params` maps
every run-all option `run_relocate3.sh` passes to its pipeline param
(`tests/real_rice_pipeline_scripts_test.py` guards the mapping), and each
sample's `genotype_bam` is the same external CRAM. `compare_pipeline_to_rc2.sh`
checks a pipeline outdir file by file against the v0.2.0rc2 run-all results.

## Dependencies

Activate the project's pixi env (`pixi shell` from `RelocaTE3/`) so
`relocaTE3` is on PATH. The reporting step additionally wants:

```bash
pip install matplotlib matplotlib-venn
```

(They're already declared in `pixi.toml` — `pixi install` from `RelocaTE3/`
covers it.) If those are absent, `compare_nonref.py` still writes the
TSV/TXT outputs and just skips the venn PNGs. Pipeline mode additionally
needs the `nf` env (`pixi install -e nf`, used by `run_pipeline.sh`), the
`blat` env's Python for `make_samplesheet.py`, and apptainer on the compute
nodes.
