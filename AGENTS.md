# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

RelocaTE3 is a pure-Python reimplementation of [RelocaTE2](https://github.com/JinfengChen/RelocaTE2) that calls transposable-element (TE) insertion polymorphisms from short-read resequencing data at single-base resolution. External dependencies are minimized to `minimap2`, `samtools`, `bedtools`, and the Python libraries `pysam`/`biopython`/`pybedtools`. The roadmap in `plans/PLAN.md` targets a validated Python reference, then Rust acceleration of per-read hotspots, then Nextflow scatter.

Status: end-to-end pipeline runs; rice Chr3 2 Mb acceptance benchmark recovers ~178/200 simulated insertions (~89% recall, ~90% precision) and is gated by `tests/acceptance_test.py`.

## Common commands

The project uses [pixi](https://pixi.sh) to pin tool versions (`minimap2`, `samtools`, `bedtools`) and the editable install. Prefer pixi over a bare venv so the binary tools resolve consistently.

```bash
pixi install                                # create/refresh the environment
pixi run relocaTE3 --help                   # entry point (script in pyproject.toml)
pixi run test                               # pytest -ra -q (all tests)
pixi run pytest tests/acceptance_test.py    # the benchmark acceptance gate
pixi run pytest tests/trim_test.py::test_x  # single test
pixi run docs                               # sphinx build -> docs/_build/html
pixi run validate-rice --local B_10         # smoke test against legacy RelocaTE2
pixi run validate-rice                      # full 10-sample SLURM array
pre-commit run --all-files                  # ruff, black, codespell, pydocstyle, pyupgrade
```

Without pixi, `pip install -e .` works if `minimap2`/`samtools`/`bedtools` are already on PATH.

## Architecture

The pipeline mirrors the RelocaTE2 step numbering. Each step is a library function in `src/RelocaTE3/` and is also exposed as a CLI subcommand so a workflow engine can scatter them. The full single-sample pipeline is `pipeline.run_sample`.

Step → module → CLI subcommand:

| Step | Module | Subcommand | Role |
|------|--------|------------|------|
| 0/6 | `reference_te.py` | `find-reference` / `annotate-ref` | Parse RepeatMasker `.out` (or align TE library to genome) → `existingTE.bed` and reference/shared insertion calls |
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

`src/RelocaTE3/__main__.py` is the single canonical CLI and the registered entry
point (`relocaTE3 = "RelocaTE3.__main__:main"` in `pyproject.toml`). Subcommands:
`map`, `trim`, `run` (map + trim only — TE-read identification + flank
generation, not the full pipeline), `annotate-ref`, `index-genome`,
`align-genome`, `find-insertions`, `characterize`.

`src/RelocaTE3/cli.py` is now a thin **compatibility shim** that re-exports
`main` from `__main__` (so `from RelocaTE3.cli import main` still works). Do not
add subcommands there.

If you add or rename a subcommand, edit `__main__.py` and update the README
subcommand table + `docs/source/usage.rst`. The acceptance test imports library
functions directly, so it does not catch CLI drift — the subprocess smoke tests
in `tests/main_test.py` cover the entry points.

## Testing

- `tests/` mirrors the module layout (`trim_test.py`, `insertions_test.py`, …) plus `acceptance_test.py` (the recall/precision gate on the rice Chr3 2 Mb fixture) and `pipeline_test.py` (end-to-end on a small fixture).
- Test data lives in `tests/data/`; transient outputs go under `tests/results/` (gitignored).
- The acceptance test runs the real pipeline and shells out to `minimap2`/`samtools` — it will be skipped or fail loudly if those aren't on PATH. Run it under `pixi run` to get the pinned versions.
- The `validation/real_rice/` harness runs RelocaTE3 against legacy RelocaTE2 calls on real rice samples; the `validate-rice` pixi task is the driver.

## Repository conventions

- Python ≥ 3.10; targets `pyupgrade --py310-plus`.
- Formatting/lint via pre-commit: ruff (check + format), black (py311), codespell, pydocstyle (`--convention=google`, src only).
- Google-style docstrings. Tests are exempt from pydocstyle and pyupgrade.
- `versioningit` derives `__version__` from git tags (`default-tag = "0.1.0"`). Don't hand-edit a version.
- `docs/source/` is a Sphinx project; build with `pixi run docs`.

## Roadmap context

`plans/PLAN.md` is the source of truth for what's intentionally unimplemented vs. a bug. Before adding a feature that feels missing (e.g. TSD-unknown mode, Rust hotspots, Nextflow), check there first. `notes/` holds in-progress investigation notes (e.g. `potential_issues_*.md`).
