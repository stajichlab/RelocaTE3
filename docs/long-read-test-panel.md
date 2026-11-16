# Small long-read test panel (`chr3_lr_small`)

Recorded: 2026-11-16.

## What it is

A small, fast long-read test set for long-read development:

- the Chr3 2 Mb fixture with 200 homozygous TE insertions;
- ONT-HQ and HiFi reads at 5x and 20x (`cov5x_rep1`, `cov20x_rep1`);
- a reference-only control read set for each sample.

`tests/longread_small_panel_test.py` checks the panel's shape: the manifest,
200 homozygous truth events, and the truth copy in `config/`.

## Where it is built

The panel is not built in this repository.

- Toolkit: `simulate_data`, branch `feat/chr3-lr-small-panel`. Subcommands:
  `simulate-data te-panel-from-gff`, `simulate-data longread-panel-reads`,
  `simulate-data check-longread-panel`.
- Workspace pipeline: `make_simulation_new/pipeline/make_chr3_lr_small/`
  (see its `README.md`).
- Output: `make_simulation_new/results/chr3_lr_small/`. A finished build has
  `.complete` and `check_result PASS` in `build_info.tsv`.

Build, from the workspace root:

```bash
bash pipeline/make_chr3_lr_small/submit.sh
```

## Use from RelocaTE3

The path below is specific to the HPCC cluster. On another system, point
to your own build.

From the repository root:

```bash
ln -s /rhome/nmath020/wessler_bigdata/rice/make_simulated_genome/make_simulation_new/results/chr3_lr_small tests/data/longread_small
pixi run python -m pytest tests/longread_small_panel_test.py -o addopts="" -q
```

`tests/data/longread_small` is gitignored. Or skip the link and set the path:

```bash
RELOCATE3_LR_SMALL_PANEL=/path/to/chr3_lr_small \
  pixi run python -m pytest tests/longread_small_panel_test.py -o addopts="" -q
```

Without the link or the variable, the 4 panel tests skip and print where they
looked.

## Provenance

The panel was first built by code inside RelocaTE3 (PR #70, closed without
merging). That code moved to the `simulate_data` toolkit. The rebuilt panel
was verified identical to the RelocaTE3-built one. Design context:
`plans/2026-11-09-long-read-algorithm-proposal.md`; the old build plan
`plans/2026-11-13-long-read-l1-small-panel-plan.md` is superseded.
