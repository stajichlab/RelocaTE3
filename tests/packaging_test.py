"""Packaging metadata: runtime dependencies are used, build requirements are
bounded and minimal.

* Every requirement in ``requirements.txt`` (the package's runtime
  dependencies, pyproject ``[tool.setuptools.dynamic]``) is imported by
  ``src/RelocaTE3``; an unused one only costs install time (pybedtools is an
  sdist that needs a C++ compiler in the image build).
* ``[build-system].requires`` holds only what builds the wheel (setuptools +
  versioningit), each with explicit bounds.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src" / "RelocaTE3"

# distribution name -> import name, where they differ.
IMPORT_NAME = {"biopython": "Bio", "pybedtools": "pybedtools", "pysam": "pysam"}


def _requirements() -> list[str]:
    names = []
    for line in (ROOT / "requirements.txt").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            names.append(re.split(r"[<>=!~\[; ]", line, maxsplit=1)[0].lower())
    return names


def _imported_top_level_modules() -> set[str]:
    modules = set()
    for path in SRC.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                modules.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                modules.add(node.module.split(".")[0])
    return modules


def test_every_runtime_requirement_is_imported():
    imported = _imported_top_level_modules()
    reqs = _requirements()
    assert reqs, "requirements.txt is empty"
    for name in reqs:
        assert IMPORT_NAME.get(name, name) in imported, f"{name} is required but never imported"


def test_pixi_core_has_no_unused_python_packages():
    deps = tomllib.loads((ROOT / "pixi.toml").read_text())["dependencies"]
    for name in ("biopython", "pybedtools"):
        assert name not in deps, f"pixi.toml [dependencies] still lists unused {name}"


def test_build_requirements_are_bounded_and_minimal():
    requires = tomllib.loads((ROOT / "pyproject.toml").read_text())["build-system"]["requires"]
    names = {re.split(r"[<>=!~ ]", r, maxsplit=1)[0].lower(): r for r in requires}
    assert set(names) == {"setuptools", "versioningit"}, requires
    for name, req in names.items():
        assert ">=" in req and "<" in req.replace("<=", ""), f"{req}: needs a lower and an upper bound"
