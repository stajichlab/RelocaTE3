# CLI Canonical Relocation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Move the canonical CLI implementation from `__main__.py` into `cli.py`, leaving `__main__.py` as a tiny `python -m` launcher — the conventional modern Python layout — with zero change to command behavior.

**Architecture:** A literal relocation. `git mv src/RelocaTE3/__main__.py src/RelocaTE3/cli.py` (after removing the current shim), then a new 4-line `__main__.py` that imports and calls `cli.main`. Repoint the console-script entry point to `RelocaTE3.cli:main` and reinstall. The parser, handlers, flags, defaults, and dispatch are moved verbatim — no logic edits. A follow-up branch (out of scope here) will then do the dispatch cleanup (`build_parser()`, pass `Namespace` instead of `**vars(parsed)`, drop per-handler `**kwargs`).

**Tech Stack:** Python 3.10+, argparse, pytest, pixi (editable install + pinned tools), `importlib.metadata` entry-point resolution.

**Note on location:** Per repo convention, plans live in `plans/` (not the skill's default `docs/plans/`), matching the prior `plans/2026-07-10-cli-consolidation-*.md`.

**Reference skills:** @superpowers:executing-plans to drive the plan; @superpowers:test-driven-development discipline (here: keep the suite green across the refactor rather than red→green, since this is a behavior-preserving move).

---

## Background the implementer needs

The prior consolidation (PR #19, on `main`) made `src/RelocaTE3/__main__.py` the single
canonical CLI and left `src/RelocaTE3/cli.py` as a compatibility shim
(`from RelocaTE3.__main__ import main`). The console entry point is
`relocaTE3 = "RelocaTE3.__main__:main"` (`pyproject.toml`).

This branch inverts *where the code lives* (not what it does): `cli.py` becomes
the canonical module and `__main__.py` becomes the thin launcher. This is the
conventional layout (pip and many tools do it): `__main__.py` exists only for the
`python -m package` protocol; the substance lives in a normally-importable module.

**Two load-bearing details:**
1. `tests/main_test.py` does `from RelocaTE3 import __main__` and
   `monkeypatch.setattr(__main__, "cmd_map"/"cmd_characterize", ...)`. After the
   move those handlers live in `cli`, so the patch target must become `cli`.
2. `main()` sets `prog = __entry_points__.get(__name__, "relocaTE3")`.
   `__entry_points__` is built in `__init__.py` from installed console-script
   metadata (`{module: script_name}`). After repointing the entry point to
   `RelocaTE3.cli:main` and reinstalling, it becomes `{"RelocaTE3.cli": "relocaTE3"}`,
   and inside the moved `main()` `__name__ == "RelocaTE3.cli"`, so `prog` resolves
   to `"relocaTE3"`. The `"relocaTE3"` fallback also covers the not-yet-reinstalled
   and `python -m` (`__name__ == "__main__"`) cases. **No code change needed for
   prog** — but a reinstall is required so the `relocaTE3` shell wrapper regenerates.

---

## Task 1: Relocate the canonical CLI into `cli.py` + thin `__main__.py`

**Files:**
- Remove: `src/RelocaTE3/cli.py` (the shim)
- Move: `src/RelocaTE3/__main__.py` → `src/RelocaTE3/cli.py` (verbatim)
- Create: `src/RelocaTE3/__main__.py` (new thin launcher)
- Modify: `pyproject.toml:22` (entry point)
- Modify: `tests/main_test.py` (patch target `__main__` → `cli`)
- Modify: `tests/pipeline_test.py` (rename the now-misnamed "shim" test)

**Step 1: Move the files (preserve git history)**

```bash
git rm src/RelocaTE3/cli.py
git mv src/RelocaTE3/__main__.py src/RelocaTE3/cli.py
```
Git records this as a rename of the big file into `cli.py` (history follows).

**Step 2: Trim the moved `cli.py`**

In the moved `src/RelocaTE3/cli.py`, update the module docstring (line 1) to:
```python
"""RelocaTE3 command-line interface (canonical implementation)."""
```
and REMOVE the trailing module-exec guard at the end of the file:
```python
if __name__ == "__main__":
    main()
```
(The new `__main__.py` owns `python -m` execution.) Do **not** touch any
`_menu_*`, `cmd_*`, `main()`, flags, defaults, or the
`prog = __entry_points__.get(__name__, "relocaTE3")` line.

**Step 3: Create the thin launcher `src/RelocaTE3/__main__.py`**

```python
"""Module entry point so ``python -m RelocaTE3`` runs the canonical CLI."""

from __future__ import annotations

from RelocaTE3.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
```
Note: `raise SystemExit(main())` makes `python -m RelocaTE3` propagate the exit
code (the old guard called `main()` and ignored its return, so a failing command
exited 0). This aligns `python -m` with the setuptools console-script wrapper,
which already does `sys.exit(main())`. The validated harness uses the console
script, so this does not affect validated behavior.

**Step 4: Repoint the console-script entry point**

In `pyproject.toml`, change line 22:
```toml
relocaTE3 = "RelocaTE3.__main__:main"
```
to:
```toml
relocaTE3 = "RelocaTE3.cli:main"
```

**Step 5: Retarget the monkeypatches in `tests/main_test.py`**

Change the import (line 11-12):
```python
from RelocaTE3 import __main__, __version__
from RelocaTE3.__main__ import main
```
to:
```python
from RelocaTE3 import __version__, cli
from RelocaTE3.cli import main
```
Change the three patch targets (currently `__main__`):
```python
monkeypatch.setattr(cli, "cmd_map", mockreturn)          # was __main__ (x2)
monkeypatch.setattr(cli, "cmd_characterize", mockreturn) # was __main__
```
Leave the subprocess smoke tests (`python -m RelocaTE3 --version`,
`relocaTE3 --version`) unchanged.

**Step 6: Rename the now-misnamed test in `tests/pipeline_test.py`**

`test_cli_shim_reaches_canonical_parser` is no longer a "shim" test (cli IS
canonical). Rename to `test_cli_main_handles_installed_subcommand` and update its
docstring to: "RelocaTE3.cli.main is the canonical parser and handles installed
subcommands (index-genome exists only on it)." Keep the assertion
`assert main(["index-genome", "--help"]) == 0`. The import
`from RelocaTE3.cli import main` is unchanged (still valid).

**Step 7: Reinstall so the console script regenerates against the new entry point**

Run: `pixi run pip install -e . --no-deps --force-reinstall`
Then verify the wrapper repointed:
```bash
pixi run python -c "import importlib.metadata as m; print([ (e.name,e.value) for e in m.entry_points(group='console_scripts') if e.dist.name=='RelocaTE3'])"
```
Expected: `[('relocaTE3', 'RelocaTE3.cli:main')]`.

**Step 8: Verify wiring end to end (no behavior change)**

```bash
pixi run relocaTE3 --version              # prints version, exit 0
pixi run relocaTE3 --help                 # prog shows as "relocaTE3", lists all 8 subcommands
pixi run relocaTE3 run --help             # same flags as before (-l/-r/-T/-n/-o/--threads/--aligner/--min-match/--min-trimmed/--mismatch)
pixi run python -m RelocaTE3 --version    # prints version, exit 0
```
Expected: all succeed; `--help` prog is `relocaTE3` (not `cli.py` or `__main__.py`).

**Step 9: Full test suite green**

Run: `pixi run test`
Expected: all pass (same count as main — 48). In particular `main_test.py`'s
monkeypatch tests and the subprocess smoke tests pass.

**Step 10: Lint**

Run: `~/.local/bin/ruff check src/RelocaTE3/cli.py src/RelocaTE3/__main__.py tests/main_test.py tests/pipeline_test.py`
Run: `~/.local/bin/ruff format --check src/RelocaTE3/__main__.py tests/main_test.py tests/pipeline_test.py`
Expected: clean (format the new `__main__.py` if flagged).

**Step 11: Commit**

```bash
git add -A
git commit -m "refactor(cli): relocate canonical CLI to cli.py; __main__.py is a thin launcher"
```

---

## Task 2: Update docs to the new layout

**Files:**
- Modify: `AGENTS.md` (the "The CLI lives in one file" section from PR #19)
- Check: `README.md`, `docs/source/*.rst` for any `__main__`-as-canonical wording

**Step 1: Rewrite the AGENTS.md CLI section**

Update it to state: the canonical CLI is `src/RelocaTE3/cli.py`; the console entry
point is `relocaTE3 = "RelocaTE3.cli:main"`; `src/RelocaTE3/__main__.py` is a thin
launcher for `python -m RelocaTE3` that delegates to `cli.main`. New/renamed
subcommands are added in `cli.py` (+ README table). Keep the note that the
subprocess smoke tests in `tests/main_test.py` cover both entry points.

**Step 2: Grep for stale references**

Run: `grep -rniE "RelocaTE3\.__main__:main|__main__\.py.*canonical|canonical.*__main__" --include="*.md" --include="*.rst" --include="*.toml" . | grep -v .pixi | grep -v plans/`
Expected: no remaining claims that `__main__.py` is the canonical CLI (plans/ are
historical and may keep their wording). Fix any that are user-facing.

**Step 3: Commit**

```bash
git add AGENTS.md README.md docs/
git commit -m "docs: canonical CLI is cli.py; __main__.py is the python -m launcher"
```

---

## Task 3: Validation gate

The relocation is a literal move, so the per-read algorithm cannot change; the
only risk is import/entry-point wiring, which Task 1 Steps 7-9 already exercise
(monkeypatch CLI tests, subprocess entry-point tests, `relocaTE3`/`python -m`
smoke). That is the **required** gate.

**Step 1 (required): confirm the required gate is green**

Re-confirm from Task 1: `pixi run test` all pass; `relocaTE3 --help`,
`relocaTE3 run --help`, `python -m RelocaTE3 --version` all succeed; entry point
metadata is `RelocaTE3.cli:main`.

**Step 2 (recommended end-to-end proof): B_10 parity**

Because this touches the console-script path the real-rice harness depends on, an
end-to-end B_10 check is worth doing for sign-off (it re-confirms the console
script → `cli:main` → handlers path produces identical calls):

Baseline on `main` (the parent), before the move — or reuse if still available:
```bash
git stash list   # ensure clean; baseline must be the pre-move code
# (from main, or note that main == pre-move behavior)
pixi run validate-rice --local --force B_10
# copy report/ + results/B_10/results OUT to a scratch dir
```
After the move (on this branch):
```bash
pixi run validate-rice --local --force B_10
```
Compare **sorted** normalized records (not checksums): the report `*.tsv`,
`summary.txt`, and `results/ALL.mping.all_nonref_insert*.{txt,gff}` must be
identical; recall/precision/Jaccard/TSD/status unchanged. Any diff is a
regression.

**Step 3: If B_10 run, record the result** in this plan file and commit.

---

## Done criteria

- `src/RelocaTE3/cli.py` is the canonical CLI (parser + handlers + `main`), moved
  verbatim; `src/RelocaTE3/__main__.py` is a 4-line launcher.
- `pyproject.toml` entry point is `RelocaTE3.cli:main`; reinstall confirmed.
- `pixi run test` green (48); both entry-point smoke tests pass; `relocaTE3 --help`
  prog is `relocaTE3`.
- Docs describe the new layout.
- (Recommended) B_10 output identical before/after.
- Dispatch cleanup (`build_parser()`, `Namespace` dispatch, drop `**kwargs`) is
  explicitly NOT in this branch — it is the next one.

---

## Follow-up (next branch, not here)

Dispatch cleanup on `cli.py`: extract `build_parser()`, change
`parsed.func(**vars(parsed))` to pass the `Namespace` (or explicit args), and drop
the per-handler `**kwargs` sinks. This also fixes the filed defect where `cmd_run`
silently drops `--min-match/--min-trimmed/--mismatch` (see
`todo/cmd-run-drops-trim-thresholds.md`).
