# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

RelocaTE3 is a pure-Python reimplementation of [RelocaTE2](https://github.com/JinfengChen/RelocaTE2) that calls transposable-element (TE) insertion polymorphisms from short-read resequencing data at single-base resolution. The shipped defaults match RelocaTE2: `blat` for the TE-library search (`--te-aligner`) and `bwa aln` for genome placement (`--genome-aligner bwaaln`); it also calls `samtools`, `minimap2` (genome index), `seqtk` (FASTQ to FASTA for BLAT) and the Python library `pysam` (plus `bcftools` for excision mode). `blat` lives in its own pixi environment (`-e blat`; `pipeline` = blat + MultiQC), not in `default`. RelocaTE2 parity was reached in v0.2.0rc2. Roadmap order, set 2026-11-02: Nextflow pipeline (`plans/2026-11-02-nextflow-pipeline-design.md`) → long reads → multiple reference genomes as separate FASTAs → pangenome reference → Rust acceleration.

Status: at the defaults the rice Chr3 2 Mb acceptance benchmark recovers 196/200 simulated insertions at precision 1.000 (F1 0.990), the same as RelocaTE2, gated by `tests/acceptance_test.py`. The Nextflow pipeline (`main.nf`, below) runs many samples with per-sample output identical to `run-all --genotype`; v0.3.0rc1 is its first release.

## Common commands

The project uses [pixi](https://pixi.sh) to pin tool versions and the editable install (environments in `pixi.toml`: `default` = core tools + plotting, no blat; `blat`; `pipeline` = blat + MultiQC; `nf` = Nextflow + nf-test). Prefer pixi over a bare venv so the binary tools resolve consistently.

```bash
pixi install                                # create/refresh the environment
pixi run -e blat relocaTE3 --help          # entry point; -e blat for the default TE aligner
pixi run test                               # pytest -ra -q (all tests)
pixi run -e blat pytest tests/acceptance_test.py   # the benchmark acceptance gate (needs blat)
pixi run pytest tests/trim_test.py::test_x  # single test
pixi run docs                               # sphinx build -> docs/_build/html
pixi run validate-rice --mode run-all --local B_10   # smoke test against legacy RelocaTE2
pixi run validate-rice --container-tag sha-<short>  # 10 samples through the pipeline (default mode)
pixi run validate-rice --mode run-all               # 10-sample run-all SLURM array
pre-commit run --all-files                  # ruff, black, codespell, pydocstyle, pyupgrade
```

Without pixi, `pip install -e .` works if `blat`, `bwa`, `samtools`, `minimap2`, `seqtk` (and `bcftools` for excision) are already on PATH.

## Architecture

The pipeline mirrors the RelocaTE2 step numbering. Each step is a library function in `src/RelocaTE3/` and is also exposed as a CLI subcommand so a workflow engine can scatter them. The full single-sample pipeline is `pipeline.run_sample`.

Step → module → CLI subcommand:

| Step | Module | Subcommand | Role |
|------|--------|------------|------|
| 0/6 | `reference_te.py` | `annotate-ref` (installed); `find-reference` (planned) | Align TE library to genome → `existingTE.bed`; reference/shared insertion calls are library-only until `find-reference` is exposed |
| 1 | `align.py` (`Aligner.index_genome`) | `index-genome` | `samtools faidx` + minimap2 index |
| 3 | `librelocate.py` (`RelocaTE`), `trim.py` | `trim` | Map reads to TE library, score TE alignments by boundary+match, classify 5′/3′/middle, emit flanking FASTQs and the `read_repeat_name` table |
| 4 | `genome_align.py`, `align.py` | `align-genome` | Re-align flanking reads + supporting mates to the genome (sorted, indexed BAM) |
| 5 | `insertions.py` (`InsertionFinder`) | `find-insertions` | Cluster junction + supporting reads, detect TSD, call non-reference insertions, write `all_nonref_insert.{gff,txt}` |
| 7 | `characterize.py` (`Characterizer`) | `characterize` | Genotype each insertion (homozygous/heterozygous/somatic) and optionally detect excision from genome BAMs |
| all | `pipeline.py` (`run_sample`) | `run` | Single-sample driver wiring 3 → 4 → 5 → 6 → 7 |

Supporting modules:
- `ReadLibrary.py` — abstracts SE/PE FASTQ inputs and read-group metadata.
- `models.py` — shared dataclasses (insertions, junction reads, etc.).

### The CLI lives in one file

`src/RelocaTE3/cli.py` is the single canonical CLI (argparse parser + `_menu_*`
builders + `cmd_*` handlers + `main()`) and the registered entry point
(`relocaTE3 = "RelocaTE3.cli:main"` in `pyproject.toml`). Subcommands: `map`,
`trim`, `run` (map + trim only — TE-read identification + flank generation, not
the full pipeline), `annotate-ref`, `index-genome`, `align-genome`,
`find-insertions`, `characterize`.

`src/RelocaTE3/__main__.py` is a thin launcher for `python -m RelocaTE3` that
delegates to `cli.main` (`raise SystemExit(main())`, so `python -m` propagates
the CLI's exit code). Don't put logic there.

**Access points (use these, in order):** the installed console script
`relocaTE3` (regenerated from the `RelocaTE3.cli:main` entry point on
install/reinstall); `python -m RelocaTE3`; or, in tests/library code,
`from RelocaTE3.cli import main`. All three reach the *same* `cli.main`.

**Anti-regression:** there is exactly **one** CLI module (`cli.py`). Never
reintroduce a second parser/front-end or move `main` out of `cli.py`. A past
`__main__.py`-vs-`cli.py` split silently diverged and caused wrong-file edits and
two null-result validation runs.

If you add or rename a subcommand, edit `cli.py` and update the README
subcommand table + `docs/source/usage.rst`. The acceptance test imports library
functions directly, so it does not catch CLI drift — the subprocess smoke tests
in `tests/main_test.py` cover both entry points.

### Nextflow pipeline

`main.nf` → `workflows/relocate3.nf` → `subworkflows/local/` → `modules/local/<process>/` (one `relocaTE3` subcommand each); config in `nextflow.config` and `conf/`. The core rule: no detection logic in Nextflow. Each process passes exactly the arguments `cmd_run_all` (`cli.py`) gives that step, so per-sample output equals `run-all --genotype`; `tests/nextflow_modules_args_test.py` and `tests/nextflow_schema_test.py` pin modules, schema and config to `cli.py`.

```bash
pixi run -e nf nf-test test               # pipeline + module tests
pixi run -e nf nextflow lint .
sbatch tests/nextflow/parity_chr3.sh      # parity gate vs run-all on Chr3; must print RESULT: PASS
```

## Testing

- `tests/` mirrors the module layout (`trim_test.py`, `insertions_test.py`, …) plus `acceptance_test.py` (the recall/precision gate on the rice Chr3 2 Mb fixture) and `pipeline_test.py` (end-to-end on a small fixture).
- Test data lives in `tests/data/`; transient outputs go under `tests/results/` (gitignored).
- The acceptance test runs the real pipeline: one gate pinned to minimap2 (runs in `default`) and the RelocaTE2-parity gate on the shipped defaults (blat + bwa aln), which skips without `blat` on PATH. Run it with `pixi run -e blat pytest tests/acceptance_test.py` to get both with the pinned versions.
- The `validation/real_rice/` harness runs RelocaTE3 against legacy RelocaTE2 calls on real rice samples; the `validate-rice` pixi task is the driver.

## Repository conventions

- Python ≥ 3.10; targets `pyupgrade --py310-plus`.
- Formatting/lint via pre-commit: ruff (check + format), black (py311), codespell, pydocstyle (`--convention=google`, src only).
- Google-style docstrings. Tests are exempt from pydocstyle and pyupgrade.
- `versioningit` derives `__version__` from git tags (`default-tag = "0.1.0"`); commits after tag `vX` get `X.postN+g<sha>`. The Python version is never hand-edited. The pipeline's `manifest.version` and default `container_tag` are (Release steps below).
- `docs/source/` is a Sphinx project; build with `pixi run docs`.

## Release steps

The release commit sets `manifest.version`, the default `container_tag` (`nextflow.config` + `nextflow_schema.json`, then `python docs/source/_gen_params.py`), `CITATION.cff` and the CHANGELOG heading to the version; `containers/check_release.sh` gates this on the `v*` tag in `container.yml`. After the tag, `bash scripts/post_release.sh <next>.devN` on a new branch returns to a development version that sorts after the release in PEP 440 (after `0.3.0rc1`: `0.3.0rc2.dev0`, image `dev`; `0.3.0dev` sorts before it). Bump the `-r v<version>` / `--container_tag <version>` examples in `README.md` and `docs/source/pipeline.rst` to the new version in the same commit.

Tag the merge commit of the release PR explicitly, never whatever `HEAD` is:

```bash
MERGE_SHA=$(gh pr view <n> --json mergeCommit -q .mergeCommit.oid)   # or: git log --merges -1 --format=%H origin/main
git tag -a v<version> "$MERGE_SHA" -m "RelocaTE3 <version>"
git push origin v<version>
```

Full checklist: `docs/2026-11-07-release-v0.3.0rc1.md`.

## Roadmap context

`plans/PLAN.md` is the source of truth for what's intentionally unimplemented vs. a bug. Before adding a feature that feels missing (e.g. TSD-unknown mode, Rust hotspots, Nextflow), check there first. `notes/` holds in-progress investigation notes (e.g. `potential_issues_*.md`).
