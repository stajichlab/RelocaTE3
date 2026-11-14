# Long-read L1a: small Chr3 long-read panel

Plan: `plans/2026-11-13-long-read-l1-small-panel-plan.md`. Branch
`feat/long-read-l1-small-panel`. Usage: `validation/longread_small/README.md`.

## 2026-11-14 08:00 PST: panel built and checked

### Purpose
A small ONT-HQ + HiFi panel (5x, 20x) over the Chr3 2 Mb fixture, 200
homozygous insertions, so long-read code (L2 onward) can be tested in minutes.

### Current status
- Panel built and checked: SLURM job 32073062, `RESULT: PASS (8 FASTQs
  checked)`.
- The panel is **outside the repository**:

  ```bash
  export RELOCATE3_LR_SMALL_PANEL=/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/validation_data/longread_small/panel_v1
  ```

  The in-repo default `validation/longread_small/panel_v1/` does not exist.
  Without the export, `tests/longread_small_panel_test.py` skips.
- Added: `validation/longread_small/` (`make_panel_config.py`,
  `simulate_reads.py`, `check_panel.py`, `build_panel.slurm`, `README.md`)
  and tests `tests/longread_small_{config,reads,check,panel}_test.py`.
  No file under `src/` changed.

### Commands

```bash
# build (from repo root, never on the login node)
mkdir -p logs
export RELOCATE3_LR_SMALL_PANEL=<panel dir>
export SIM_PROJECT_DIR=<simulate_data checkout>   # default is nmath020's home path
sbatch validation/longread_small/build_panel.slurm

# result
tail -3 logs/lr_small_panel.<jobid>.log
cat "${RELOCATE3_LR_SMALL_PANEL}/qc/panel_check.tsv"

# tests
RELOCATE3_LR_SMALL_PANEL=<panel dir> pixi run test
```

### Build
- Job 32073062 on `epyc`, 8 CPUs, 16 GB. Start 17:45:20, end 18:11:55 PST
  on 2026-11-13: about 27 minutes.
- Log: `<panel>/logs/lr_small_panel.32073062.log` (copied from the
  worktree's `logs/`, which is deleted with the tree).
- Built with RelocaTE3 `3d462f8`, toolkit `1e19882` (21 uncommitted paths).
  Both are in `<panel>/build_info.tsv`.
- Size: 175 MB (the plan estimated about 130 MB). Reads are 163 MB of it;
  the sample genome is 2.53 Mb, not 2.0 Mb, so sample reads at the same
  coverage carry more bases.

### Numbers (`<panel>/qc/panel_check.tsv`)

| platform | sample | kind | reads | bases | coverage | N50 | max length |
|---|---|---|---:|---:|---:|---:|---:|
| ont-hq | cov5x_rep1 | reads | 1012 | 12651833 | 5.00 | 17251 | 75755 |
| ont-hq | cov5x_rep1 | control_reads | 846 | 10000089 | 5.00 | 16128 | 90843 |
| ont-hq | cov20x_rep1 | reads | 4215 | 50607236 | 20.00 | 16377 | 61285 |
| ont-hq | cov20x_rep1 | control_reads | 3247 | 40000096 | 20.00 | 17037 | 69899 |
| hifi | cov5x_rep1 | reads | 841 | 12572108 | 4.97 | 15210 | 21698 |
| hifi | cov5x_rep1 | control_reads | 666 | 9972882 | 4.99 | 15223 | 21425 |
| hifi | cov20x_rep1 | reads | 3349 | 50324905 | 19.89 | 15211 | 24134 |
| hifi | cov20x_rep1 | control_reads | 2659 | 39764297 | 19.88 | 15128 | 24134 |

- All 8 rows `status ok`.
- Short-read fit of the rebuilt genome: `unexplained_short_reads_sample` 6
  vs `unexplained_short_reads_reference` 30029 (0.02%; limit 5%). The
  rebuilt sample genome explains the original fixture's short reads.
- HiFi `ccs` yield (from each `simulation.tsv`; floor 0.80):

  | sample | ccs_yield |
  |---|---:|
  | cov5x_rep1 | 0.9964 |
  | cov5x_rep1_reference_control | 1.0000 |
  | cov20x_rep1 | 0.9973 |
  | cov20x_rep1_reference_control | 0.9966 |

- Sample reads' coverage is measured on the sample genome (2,530,357 bp);
  control reads on the reference (2,000,001 bp).
- `truth_events.tsv`: 200 rows, all `biological_class` `homozygous`.

### Tests
- New: `tests/longread_small_panel_test.py`, 4 tests. Skip without the
  panel; pass with `RELOCATE3_LR_SMALL_PANEL` set to the panel above.
- `tests/longread_small_*_test.py`: 96 tests, all pass with the panel
  (config 41, check 31, reads 20, panel 4). The plan expected 19; the rest
  were added by review fixes in Tasks 1–5.
- Full suite, `RELOCATE3_LR_SMALL_PANEL=<panel> pixi run test` on compute
  node c17 (2026-11-14, at `ec787a9`): 844 passed, 50 skipped, 0 failed
  (`pixi run python -m pytest -o addopts="" -q`).
  Skips are environment-only: no `blat`/`seqtk`, no `nf` or `pipeline` env,
  no `yaml`. Without the export, the 4 panel tests also skip.
- `pixi run test` prints no summary line (the task's `-q` plus `addopts`
  `-q` give `-qq`); counts were taken from the progress line and
  `pytest --co`.
- Task 0 baseline count was not recorded.
- Lint: pre-commit passes on the longread_small files. The default
  `pixi exec --spec pre-commit` picks Python 3.15, where pyupgrade crashes
  (`TypeError: cannot use a bytes pattern`) and black's `python3.11` hook
  cannot build. Working command:

  ```bash
  PRE_COMMIT_HOME=<scratch dir> pixi exec --spec pre-commit --spec python=3.11 \
    pre-commit run --files validation/longread_small/*.py tests/longread_small_*_test.py
  ```

### Failures / issues
- The toolkit's editable install is broken: no `simulate_data` import and
  no `simulate-data` launcher in its pixi env. `build_panel.slurm` runs
  `python -m simulate_data.cli` with `PYTHONPATH=${SIM_PROJECT_DIR}/src`.
  This depends on the toolkit's source layout.
- The QC column is `unexplained_short_reads_*`, not the plan's
  `clipped_short_reads_*`. Commit `f9513fa` renamed it when the check began
  counting unmapped reads as well as clipped ones (an empty or wrong genome
  used to pass as 0 vs 0).

### Review follow-ups (not fixed)
- `reads/*/cov*_rep1/simulation.tsv` records the sample genome as an
  absolute path, because the panel path was given as absolute.
  Reference-control rows are relative (`tests/data/...`). The panel's
  provenance is therefore not portable; use `genome_sha256` to identify the
  genome.
- `SIM_PROJECT_DIR` defaults to a hard-coded home path. Other users must
  set it.
- The toolkit ran from an uncommitted working tree (21 dirty paths). Its
  commit alone does not reproduce the build.
- Commit `b8cff7b` changed `build_panel.slurm` **after** the build:
  provenance captured after the input checks, a RelocaTE3 dirty count
  (`relocate3_dirty_paths`), `pixi run --locked` for RelocaTE3-env calls, an
  atomic `truth_events.tsv` copy, `completed_relocate3_commit` and
  `completed_slurm_job_id` rows appended on success, and more parameters
  echoed. The panel was built with `3d462f8`, so its `build_info.tsv` lacks
  these rows. The `b8cff7b` script, and the log copy added in `ec787a9`, have not been
  run end to end (the panel's log copy was made by hand). If a job is killed
  after `check_panel.py` writes `.complete` but before the `completed_*`
  rows and log copy, reruns exit early and those steps never happen; add
  them by hand.

### Decisions
- Panel written outside git (175 MB, of which 163 MB is FASTQ), selected by
  `RELOCATE3_LR_SMALL_PANEL`.
- The mutated genome is rebuilt from the fixture's truth GFF with
  `te-panel`, and checked against the original short reads (above).
- `RiceTE.fa`'s three ambiguity codes are replaced for simulation only;
  changes in `config/te_library_changes.tsv`.
- Homozygous only for `panel_v1`.

### Next steps
1. L2 plan: first long-read detection step, tested on this panel.
2. Benchmark long-read manifest plan (`relocate-benchmark`), next to L3.
3. HiFi 30x repair plan (`make_simulation_new`).
4. Optional: rebuild the panel with the `b8cff7b` script and a relative
   panel path to get full provenance rows and portable `simulation.tsv`.
