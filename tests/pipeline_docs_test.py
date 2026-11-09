"""The pipeline user guide must not drift from the pipeline.

- docs/source/pipeline.rst's parameter tables are generated from
  nextflow_schema.json by docs/source/_gen_params.py and committed; the
  committed text must equal a fresh generation, and list every schema param.
- Every ``--param`` and ``-profile`` name in a ``nextflow run`` command of the
  guide's code blocks, or of the README's "Nextflow pipeline" section, must
  exist in the schema / nextflow.config. Single-dash Nextflow options must be
  real ``nextflow run`` options. Keys of YAML blocks (params files) must be
  schema params.
"""

from __future__ import annotations

import importlib.util
import json
import re
import shlex
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
GUIDE = ROOT / "docs" / "source" / "pipeline.rst"
INDEX = ROOT / "docs" / "source" / "index.rst"
README = ROOT / "README.md"
SCHEMA = ROOT / "nextflow_schema.json"
CONFIG = ROOT / "nextflow.config"
GEN = ROOT / "docs" / "source" / "_gen_params.py"

# `nextflow run` options used in the docs, all listed by `nextflow run -h`
# (Nextflow 26.10.4). Extend only with options that command lists.
NEXTFLOW_RUN_OPTIONS = {
    "-r", "-resume", "-profile", "-params-file", "-work-dir", "-w",
    "-preview", "-latest", "-with-report", "-with-trace", "-with-timeline",
    "-with-dag", "-name", "-stub-run",
}

# Options of `nextflow` itself (before `run`), from `nextflow -h`.
NEXTFLOW_TOP_OPTIONS = {"-c", "-config", "-log"}


def _gen_module():
    spec = importlib.util.spec_from_file_location("_gen_params", GEN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def schema_params() -> set[str]:
    schema = json.loads(SCHEMA.read_text())
    return {p for group in schema["$defs"].values() for p in group["properties"]}


def config_profiles() -> set[str]:
    text = CONFIG.read_text()
    block = text[text.index("\nprofiles {"):]
    return set(re.findall(r"^    (\w+) \{", block, flags=re.M))


def rst_code_blocks(text: str) -> list[tuple[str, str]]:
    """(language, body) of every ``.. code-block::`` / ``::`` literal block."""
    blocks = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.match(r"^(\s*)\.\. code-block::\s*(\S*)", line)
        literal = line.rstrip().endswith("::") and not line.lstrip().startswith("..")
        if m or literal:
            indent = len(m.group(1)) if m else len(line) - len(line.lstrip())
            lang = m.group(2) if m else ""
            i += 1
            body = []
            while i < len(lines):
                cur = lines[i]
                if cur.strip() and len(cur) - len(cur.lstrip()) <= indent:
                    break
                if not (m and not body and cur.strip().startswith(":")):
                    body.append(cur)
                i += 1
            blocks.append((lang, "\n".join(body)))
            continue
        i += 1
    return blocks


def readme_nextflow_section() -> str:
    text = README.read_text()
    start = text.index("### Nextflow pipeline")
    end = text.index("\n### ", start + 1)
    return text[start:end]


def md_code_blocks(text: str) -> list[tuple[str, str]]:
    return re.findall(r"^```(\w*)\n(.*?)^```", text, flags=re.M | re.S)


def nextflow_run_commands(body: str) -> list[list[str]]:
    """Token lists after `nextflow run`, backslash continuations joined."""
    joined = re.sub(r"\\\n", " ", body)
    commands = []
    for line in joined.splitlines():
        line = line.split(" #", 1)[0]
        # `nextflow run ...`, or `nextflow -c custom.config run ...`
        m = re.search(r"\bnextflow((?:\s+-\S+\s+\S+)*?)\s+run\s(.*)", line)
        if not m:
            continue
        top = shlex.split(m.group(1))
        assert all(w in NEXTFLOW_TOP_OPTIONS for w in top[::2]), f"unknown nextflow option in {top}"
        commands.append(shlex.split(m.group(2), comments=True))
    return commands


def doc_commands() -> list[tuple[str, list[str]]]:
    found = []
    guide = GUIDE.read_text() if GUIDE.is_file() else ""
    for lang, body in rst_code_blocks(guide):
        for cmd in nextflow_run_commands(body):
            found.append(("pipeline.rst", cmd))
    for lang, body in md_code_blocks(readme_nextflow_section()):
        for cmd in nextflow_run_commands(body):
            found.append(("README", cmd))
    return found


def test_guide_exists_and_is_in_the_toctree():
    assert GUIDE.is_file()
    toctree = INDEX.read_text().split(".. toctree::", 1)[1]
    assert re.search(r"^\s+pipeline\s*$", toctree, flags=re.M)


def test_parameter_tables_equal_a_fresh_generation():
    gen = _gen_module()
    text = GUIDE.read_text()
    assert gen.update_text(text, json.loads(SCHEMA.read_text())) == text, (
        "docs/source/pipeline.rst parameter tables are stale: run "
        "`python docs/source/_gen_params.py`"
    )


def test_every_schema_param_is_in_the_parameter_tables():
    gen = _gen_module()
    region = gen.generated_region(GUIDE.read_text())
    listed = set(re.findall(r"^\s+\* - ``--(\w+)``", region, flags=re.M))
    assert listed == schema_params()


def test_generator_check_mode_flags_a_stale_table(tmp_path):
    gen = _gen_module()
    stale = tmp_path / "pipeline.rst"
    stale.write_text(GUIDE.read_text().replace("* - ``--genome``", "* - ``--genomes``", 1))
    assert gen.main(["--check", str(stale)]) == 1
    assert gen.main(["--check", str(GUIDE)]) == 0


def test_docs_have_nextflow_commands():
    sources = {src for src, _ in doc_commands()}
    assert sources == {"pipeline.rst", "README"}


@pytest.mark.parametrize("source,cmd", doc_commands(), ids=lambda v: v if isinstance(v, str) else " ".join(v)[:60])
def test_command_params_and_profiles_exist(source, cmd):
    params = schema_params()
    profiles = config_profiles()
    i = 1  # cmd[0] is the pipeline (. or stajichlab/RelocaTE3)
    while i < len(cmd):
        word = cmd[i]
        if word.startswith("--"):
            name = word[2:].split("=", 1)[0]
            assert name in params, f"{source}: --{name} is not in nextflow_schema.json"
        elif word.startswith("-"):
            assert word in NEXTFLOW_RUN_OPTIONS, f"{source}: {word} is not a known nextflow run option"
            if word == "-profile":
                for p in cmd[i + 1].split(","):
                    assert p in profiles, f"{source}: profile {p} is not in nextflow.config"
        i += 1


def test_yaml_params_file_keys_are_schema_params():
    params = schema_params()
    yaml_blocks = [b for lang, b in rst_code_blocks(GUIDE.read_text()) if lang == "yaml"]
    assert yaml_blocks, "the guide shows a -params-file example"
    for body in yaml_blocks:
        for key in re.findall(r"^\s*(\w+):", body, flags=re.M):
            assert key in params, f"params file key {key} is not in nextflow_schema.json"


def test_readme_links_to_the_guide():
    assert "(docs/source/pipeline.rst)" in readme_nextflow_section()
