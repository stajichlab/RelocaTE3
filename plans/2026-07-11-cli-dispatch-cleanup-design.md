# CLI Dispatch Cleanup Design

Date: 2026-07-11
Author: Nathan Mathieu (with Claude Code)
Status: approved — ready for implementation planning
Branch: `refactor-cli-dispatch`

## Purpose

Modernize the internals of the single canonical CLI (`src/RelocaTE3/cli.py`)
without changing any command behavior. Today `main()` builds the argparse tree
inline and dispatches with `parsed.func(**vars(parsed))`, which forces every one
of the 8 `cmd_*` handlers to end in a `**kwargs` sink (a pure smell — none of
them read it) and to duplicate the flag list in their signatures. This is the
"actually well-designed" follow-on to the CLI relocation.

This is **branch 1 of 3** on the CLI internals: (1) this behavior-preserving
dispatch refactor; (2) fix `cmd_run` dropping `--min-match/--min-trimmed/--mismatch`
(behavior change, separate validation); (3) resolve `--aligner bwa` being accepted
but unsupported. Branches 2 and 3 are OUT OF SCOPE here.

## Current state (ground truth, `cli.py`)

- 8 `_menu_*(parser)` builders, each ending `parser.set_defaults(func=cmd_*)`.
- `main()` assembles the subparsers inline (≈ lines 637–724) and dispatches at
  `parsed.func(**vars(parsed))` (≈ line 726).
- 8 `cmd_*` handlers, each `cmd_x(a, b, …, verbose, **kwargs)`. **All 8 take
  `**kwargs` purely to swallow `command`/`func`; none read it.**

## Approach (chosen: A — Namespace dispatch)

- **A (chosen)** — handlers become `cmd_x(args: argparse.Namespace)` and read
  `args.foo`; dispatch becomes `return parsed.func(parsed)`. Conventional argparse
  idiom (matches the old `cli.py` `_run_*(args)` handlers); no magic; each handler
  becomes self-documenting about exactly which args it reads.
- B — keep explicit params, drop `**kwargs`, filter `vars(parsed)` per handler via
  `inspect.signature`. Rejected: introspection magic in the dispatcher.
- C — dataclass per command. Rejected: YAGNI.

## Design

Three changes, all inside `src/RelocaTE3/cli.py`:

1. **Extract `build_parser() -> argparse.ArgumentParser`.** Move the subparser
   assembly out of `main()` into `build_parser()`. `main()` becomes: build parser
   → `parse_args` → handle `--version` / no-args (print help) / `verbose` setup →
   `return parsed.func(parsed)`, keeping the existing try/except (KeyboardInterrupt,
   SystemExit, Exception) guard exactly as-is. The `prog = __entry_points__.get(...)`
   logic stays.

2. **Namespace dispatch.** `parsed.func(**vars(parsed))` → `parsed.func(parsed)`.
   Rewrite each of the 8 `cmd_*` signatures to `cmd_x(args)` and prefix each param
   read in the body with `args.` (mechanical; the `_menu_*` `dest` names are the
   `args` attribute names). Remove all 8 `**kwargs` sinks.

3. **No `_menu_*` / flag / default / choice changes.** Parser *definitions* are
   untouched — only the location of `main()`'s assembly code moves.

### Behavior-preserving guarantee (the crux)

`cmd_run` continues to read exactly what it reads today
(`args.left/right/te_library/name/outdir/threads/aligner/verbose`) and **still does
not forward `--min-match/--min-trimmed/--mismatch`** to `identify_TE_reads`. The
bug is preserved on purpose so B_10 stays byte-identical. Add an explicit guard
comment at that call site so no reviewer/agent "helpfully" fixes it and breaks the
gate:

```python
# BUG preserved intentionally (see todo/cmd-run-drops-trim-thresholds.md):
# --min-match/--min-trimmed/--mismatch are parsed but not forwarded here.
# Fixing changes trim output; do it in its own validated branch, not this refactor.
```

## Testing / validation

- **All 48 existing tests stay green.** `main_test.py` monkeypatches
  `cli.cmd_map`/`cli.cmd_characterize`; since dispatch now passes a `Namespace`
  positionally, update the test's `mockreturn(**kwargs)` to `mockreturn(args)` and
  assert on `args.name` / `args.left` / `args.te_library` / `args.sites_file` /
  `args.bam` (test-only change; no product impact).
- **Add one test:** `build_parser()` returns a parser that parses a known
  subcommand (e.g. `build_parser().parse_args(["index-genome", "-g", "x"]).func`
  is `cmd_index_genome`) — locks in the extraction.
- **B_10 `--local --force` gate: byte-identical to the pre-refactor baseline**
  (sorted normalized records + summaries + result files; recall/precision/Jaccard/
  TSD/status unchanged). This is the proof the refactor changed nothing.

## Out of scope (separate branches)

- `cmd_run` threshold forwarding fix (`todo/cmd-run-drops-trim-thresholds.md`).
- `--aligner bwa` accepted-but-unsupported (`todo/run-aligner-bwa-unsupported.md`).

## Files touched

- `src/RelocaTE3/cli.py` (build_parser extraction, Namespace dispatch, drop `**kwargs`)
- `tests/main_test.py` (mock accepts a `Namespace`; add `build_parser` test)
