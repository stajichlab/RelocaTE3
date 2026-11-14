# Small Chr3 long-read panel

## Purpose

A small ONT-HQ + HiFi read panel over the Chr3 2 Mb fixture
(`tests/data/sim_genome/`, 200 known insertions), so long-read code can be
tested in minutes. Reads are simulated with the `simulate_data` toolkit's
long-read functions, the same ones that built the full
`riceTElib_lr_benchmark` panel.

## Build

Submit from the repository root (never run on the login node):

```bash
mkdir -p logs
export RELOCATE3_LR_SMALL_PANEL=/path/outside/git/longread_small/panel_v1   # optional
export SIM_PROJECT_DIR=/path/to/simulate_data                               # required for anyone but nmath020
sbatch validation/longread_small/build_panel.slurm
```

- Without `RELOCATE3_LR_SMALL_PANEL` the panel is written to
  `validation/longread_small/panel_v1/` (gitignored).
- `SIM_PROJECT_DIR` defaults to a hard-coded home path
  (`/rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data`).
  Other users must set it.
- Run time: about 27 minutes on `epyc`, 8 CPUs, 16 GB (job 32073062,
  17:45 to 18:11 on 2026-11-13). Log: `logs/lr_small_panel.<jobid>.log`.
- Success ends with `RESULT: PASS (8 FASTQs checked)` and writes
  `<panel>/.complete`.
- Rerun-safe. Finished steps are skipped and incomplete output is refused. A
  panel with `.complete` ends the job at once.
- Size on disk: 175 MB.

### Where the existing panel is

The first panel was built **outside the repository**:

```bash
export RELOCATE3_LR_SMALL_PANEL=/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/validation_data/longread_small/panel_v1
```

The in-repo default `validation/longread_small/panel_v1/` does not exist, so
without this export the panel test skips.

## Layout

```text
panel_v1/
├── build_info.tsv           RelocaTE3 + toolkit commits, dirty counts, job id
├── config/                  backbone.fa, te_library.fa, panel.json,
│                            te_library_changes.tsv, truth_events.tsv
├── truth/                   the te-panel bundle (genomes/, events.tsv, ...)
├── truth_events.tsv         one row per insertion, reference coordinates
├── panel_manifest.tsv       sample platform coverage replicate reads control_reads
├── reads/<platform>/<sample>/reads.fastq.gz + simulation.tsv + .complete
├── reads/<platform>/<sample>_reference_control/reads.fastq.gz + simulation.tsv + .complete
├── qc/panel_check.tsv
└── .complete
```

- Platforms `ont-hq` and `hifi`; coverage 5x and 20x; one replicate
  (samples `cov5x_rep1`, `cov20x_rep1`).
- Base seed 916 (the full panel's seed); each FASTQ gets its own seed
  (916 to 923, recorded in `simulation.tsv`).
- Paths in `panel_manifest.tsv` are relative to the panel root.
- Each sample has a reference control: reads simulated from the unmodified
  Chr3 2 Mb reference at the same coverage.

## Truth conventions

- `truth_events.tsv` columns: `event_id chrom position te_family strand tsd
  tsd_length biological_class`. Coordinates are on the reference.
- `position` = column 4 of `MSU7.Chr3_2M.ALL.gff`.
- TSD = `reference[position-L:position]` (0-based slice, `L` =
  `tsd_length`). This holds for all 200 rows of the fixture.
- `te-panel` accepts only A/C/G/T, so the simulation uses a normalized copy
  of `tests/data/TE_lib/RiceTE.fa`. The three changed bases (all in
  `spmlike`) are listed in `config/te_library_changes.tsv`. Callers should
  still use the original `RiceTE.fa`.
- `truth/` holds the toolkit's own truth (sample genome
  `truth/genomes/sample_h1.fa`, `genome_truth.tsv` in sample-genome
  coordinates).

## Use from tests

`tests/longread_small_panel_test.py` reads the panel from
`RELOCATE3_LR_SMALL_PANEL`, or `validation/longread_small/panel_v1/` when it
is unset. It skips unless `<panel>/.complete` exists.

```bash
RELOCATE3_LR_SMALL_PANEL=<panel> pixi run pytest tests/longread_small_panel_test.py
```

## Limits

- Homozygous insertions only (one sample haplotype). Heterozygous and
  low-frequency sites come in `panel_v2`, with genotyping.
- Simulated error profiles are not yet compared with the lab's real R10.4.1
  reads.
- The toolkit was run from an uncommitted working tree (21 dirty paths at
  build time). Its commit is in `<panel>/build_info.tsv` and the job log, but
  that commit alone does not reproduce the build.
- The toolkit's editable install is broken (no `simulate_data` import, no
  `simulate-data` launcher). `build_panel.slurm` therefore runs
  `python -m simulate_data.cli` with `PYTHONPATH=${SIM_PROJECT_DIR}/src`. This
  depends on the toolkit's source layout.
- `reads/*/cov*_rep1/simulation.tsv` records the sample genome as an
  absolute path when the panel path is absolute. Reference-control rows hold
  the relative `tests/data/...` path. The `genome_sha256` column, not the
  path, is the portable identity.
