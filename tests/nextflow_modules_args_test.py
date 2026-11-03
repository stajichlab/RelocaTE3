"""Each Nextflow module passes every argument ``run-all`` passes to its step.

The pipeline's per-sample output must equal ``relocaTE3 run-all``'s, so each
process (``modules/local/relocate3_<step>/main.nf``) runs one subcommand with
the arguments ``cmd_run_all`` gives that step's handler. Nothing here is a
hand-kept list of flags: the expected flags are derived from the source.

1. ``cmd_run_all``'s ``_step("<subcommand>", handler, key=value, ...)`` calls
   are read with ``ast``: the step label is the subcommand, the keywords are
   argparse ``dest`` names.
2. Each ``dest`` is resolved to its option strings through ``build_parser()``'s
   subparser for that subcommand.
3. One of those option strings must appear in the module's ``script:`` block.

So adding an argument to a ``run-all`` step (or a new step) fails this test
until the matching module passes it too.
"""

from __future__ import annotations

import argparse
import ast
import inspect
import re
from pathlib import Path

import pytest

from RelocaTE3 import cli

ROOT = Path(__file__).resolve().parent.parent
MODULES = ROOT / "modules" / "local"

# run-all step (subcommand) -> module folder.
STEP_MODULE = {
    "index-genome": "relocate3_index_genome",
    "run": "relocate3_te_search",
    "align-genome": "relocate3_align_genome",
    "find-insertions": "relocate3_find_insertions",
    "find-reference": "relocate3_find_reference",
    "align-reads": "relocate3_align_reads",
    "characterize": "relocate3_characterize",
}


def _run_all_steps() -> dict[str, dict[str, ast.expr]]:
    """``{subcommand: {dest: value expression}}`` from cmd_run_all's _step calls."""
    tree = ast.parse(inspect.getsource(cli.cmd_run_all))
    steps: dict[str, dict[str, ast.expr]] = {}
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "_step"
        ):
            label = node.args[0]
            assert isinstance(label, ast.Constant), ast.dump(label)
            assert label.value not in steps, f"duplicate step {label.value}"
            steps[label.value] = {kw.arg: kw.value for kw in node.keywords}
    return steps


def _subparser(name: str) -> argparse.ArgumentParser:
    parser = cli.build_parser()
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action.choices[name]
    raise AssertionError("build_parser() has no subcommands")


def _passed_as_false(action: argparse.Action, value: ast.expr) -> bool:
    """A store_true flag given the constant False is simply not on the command line."""
    return (
        isinstance(action, argparse._StoreTrueAction)
        and isinstance(value, ast.Constant)
        and value.value is False
    )


def _script_block(module: str) -> str:
    text = (MODULES / module / "main.nf").read_text()
    m = re.search(r"^\s*script:\s*$(.*?)^\s*stub:\s*$", text, re.S | re.M)
    assert m, f"{module}/main.nf: no script: block followed by a stub: block"
    return m.group(1)


def _has_option(script: str, option: str) -> bool:
    return re.search(rf"(?<![\w-]){re.escape(option)}(?![\w-])", script) is not None


def _expected(step: str) -> list[tuple[str, tuple[str, ...]]]:
    """``[(dest, option strings)]`` that must appear in the step's module."""
    sub = _subparser(step)
    by_dest = {a.dest: a for a in sub._actions}
    out = []
    for dest, value in _run_all_steps()[step].items():
        assert dest in by_dest, f"run-all passes {dest}= to {step}, which has no such option"
        action = by_dest[dest]
        if _passed_as_false(action, value):
            continue
        assert action.option_strings, f"{step} {dest} is positional"
        out.append((dest, tuple(action.option_strings)))
    return out


def test_every_run_all_step_has_a_module():
    assert set(_run_all_steps()) == set(STEP_MODULE)


def test_the_introspection_sees_known_arguments():
    """Guard against the ast walk silently finding nothing."""
    expected = dict(_expected("find-insertions"))
    for dest in ("fullreads_bam", "insert_size", "distance", "require_both_junctions"):
        assert dest in expected
    assert "paired" not in dict(_expected("align-genome"))  # paired=False
    assert "force" not in dict(_expected("index-genome"))  # force=False


@pytest.mark.parametrize("step", sorted(STEP_MODULE))
def test_module_runs_its_subcommand(step):
    script = _script_block(STEP_MODULE[step])
    assert re.search(rf"\brelocaTE3\s+{re.escape(step)}(?![\w-])", script), step


@pytest.mark.parametrize("step", sorted(STEP_MODULE))
def test_module_passes_every_run_all_argument(step):
    script = _script_block(STEP_MODULE[step])
    missing = [
        f"{dest} ({'/'.join(opts)})"
        for dest, opts in _expected(step)
        if not any(_has_option(script, o) for o in opts)
    ]
    assert not missing, f"{STEP_MODULE[step]} does not pass: {', '.join(missing)}"


@pytest.mark.parametrize("step", sorted(STEP_MODULE))
def test_module_keeps_temporary_files_in_the_task_directory(step):
    """relocaTE3 writes temp files (e.g. align-reads' whole-genome SAM) via
    Python's tempfile, which honours TMPDIR. Pointing TMPDIR at the task
    directory keeps them off a small node-local /tmp and cleans up with the
    work directory; the Python code is unchanged."""
    script = _script_block(STEP_MODULE[step])
    lines = [ln.strip() for ln in script.splitlines()]
    assert 'export TMPDIR="\\${PWD}"' in lines, STEP_MODULE[step]
    assert lines.index('export TMPDIR="\\${PWD}"') < next(
        i for i, ln in enumerate(lines) if ln.startswith("relocaTE3 ")
    )


# Value-level checks (each flag carries its own param, quoting survives the
# shell) need the rendered .command.sh of a real run: see the nf-test
# "every setting reaches its flag in the right task" in tests/main.nf.test.
