"""The Nextflow pipeline's parameters cannot drift from the Python CLI.

``nextflow_schema.json`` (nf-schema validation and ``--help``) and the
``params`` block of ``nextflow.config`` both restate values the package already
defines: the aligner choices, the ``--tsd`` modes, and the ``run-all`` defaults.
These tests pin each copy to its source of truth.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
from packaging.version import Version

import RelocaTE3
from RelocaTE3.aligners import GENOME_ALIGNERS, TE_ALIGNERS
from RelocaTE3.cli import build_parser
from RelocaTE3.insertions import require_unknown_tsd

ROOT = Path(__file__).resolve().parent.parent
SCHEMA = ROOT / "nextflow_schema.json"
CONFIG = ROOT / "nextflow.config"

# Pipeline param -> run-all argparse dest. Same meaning, same default.
RUN_ALL_DEST = {
    "te_aligner": "te_aligner",
    "genome_aligner": "genome_aligner",
    "mismatch": "mismatch_allowance",
    "min_match": "minimum_match_length",
    "min_trimmed": "minimum_trimmed_length",
    "size": "insert_size",
    "min_mapq": "min_mapq",
    "tsd": "tsd",
    "require_both_junctions": "require_both_junctions",
    "target": "target",
    "te_name": "te_name",
    "repeatmasker": "repeatmasker",
    "te_opts": "te_opts",
    "genome_opts": "genome_opts",
    "distance": "distance",
    "excision": "excision",
    "bcftools": "bcftools",
}


def _schema_properties() -> dict[str, dict]:
    schema = json.loads(SCHEMA.read_text())
    props: dict[str, dict] = {}
    for group in schema["$defs"].values():
        props.update(group["properties"])
    return props


def _config_block(text: str, name: str) -> str:
    """Body of the top-level ``name { ... }`` block (no nested braces expected)."""
    m = re.search(rf"^{name}\s*\{{(.*?)^\}}", text, re.S | re.M)
    assert m, f"no top-level {name} block in nextflow.config"
    return m.group(1)


def _literal(raw: str):
    raw = raw.strip()
    if raw in ("true", "false"):
        return raw == "true"
    if raw == "null":
        return None
    if re.fullmatch(r"-?\d+", raw):
        return int(raw)
    m = re.fullmatch(r"""'([^']*)'|"([^"]*)\"""", raw)
    assert m, f"unsupported literal in nextflow.config params: {raw}"
    return m.group(1) if m.group(1) is not None else m.group(2)


def _config_params() -> dict:
    body = _config_block(CONFIG.read_text(), "params")
    params = {}
    for line in body.splitlines():
        line = re.sub(r"(^|\s)//.*$", "", line).strip()
        if not line:
            continue
        key, _, value = line.partition("=")
        params[key.strip()] = _literal(value)
    return params


def _run_all_defaults():
    return build_parser().parse_args(
        ["run-all", "-1", "r1.fq", "-T", "te.fa", "-g", "g.fa", "-n", "s"]
    )


def test_aligner_enums_match_the_python_tuples():
    props = _schema_properties()
    assert props["te_aligner"]["enum"] == list(TE_ALIGNERS)
    assert props["genome_aligner"]["enum"] == list(GENOME_ALIGNERS)
    assert "blat" not in props["genome_aligner"]["enum"]


def test_tsd_enum_is_the_unknown_mode_only():
    tsd = _schema_properties()["tsd"]
    assert tsd["enum"] == ["UNK", "UKN", "unknown"]
    assert tsd["default"] == "UNK"
    for value in tsd["enum"]:
        assert require_unknown_tsd(value) == value


def test_schema_and_config_declare_the_same_params():
    assert set(_schema_properties()) == set(_config_params())


def test_schema_defaults_equal_config_defaults():
    props = _schema_properties()
    config = _config_params()
    for name, value in config.items():
        if value is None or name == "pixi_env":
            # null means "no default"; pixi_env is built from ${projectDir}.
            assert "default" not in props[name], name
            continue
        assert props[name].get("default") == value, name


def test_config_defaults_equal_run_all_defaults():
    config = _config_params()
    run_all = _run_all_defaults()
    for param, dest in RUN_ALL_DEST.items():
        assert config[param] == getattr(run_all, dest), param


def test_genotype_is_on_in_the_pipeline_and_opt_in_for_run_all():
    """Deliberate difference (design doc): the pipeline fills the genotype matrix."""
    assert _config_params()["genotype"] is True
    assert _run_all_defaults().genotype is False


def test_manifest_version_never_lags_the_package():
    """The pipeline version is the RelocaTE3 version it ships with.

    nextflow.config cannot read ``git describe``, so the manifest holds the
    version under development. Its release part may be ahead of the installed
    package (e.g. ``0.3.0rc2.dev0`` while the last tag is 0.2.0rc2) but never behind,
    and on a tagged final release the two must be equal.
    """
    manifest = _config_block(CONFIG.read_text(), "manifest")
    m = re.search(r"""^\s*version\s*=\s*'([^']+)'""", manifest, re.M)
    assert m, "manifest.version not set"
    pipeline = Version(m.group(1))
    package = Version(RelocaTE3.__version__)
    assert pipeline.release >= package.release
    if not (package.is_prerelease or package.is_postrelease or package.local):
        assert pipeline == package


@pytest.mark.parametrize("name", ["input", "genome", "te_library", "outdir"])
def test_shared_inputs_are_required(name):
    schema = json.loads(SCHEMA.read_text())
    required = {r for g in schema["$defs"].values() for r in g.get("required", [])}
    assert name in required


# ------------------------------------------------------- container params

IMAGE = "ghcr.io/stajichlab/relocate3"


def _manifest_version() -> Version:
    manifest = _config_block(CONFIG.read_text(), "manifest")
    m = re.search(r"""^\s*version\s*=\s*'([^']+)'""", manifest, re.M)
    assert m, "manifest.version not set"
    return Version(m.group(1))


def test_container_params_are_declared():
    props = _schema_properties()
    config = _config_params()
    assert config["container_registry"] == IMAGE
    assert props["container_registry"]["default"] == IMAGE
    assert props["container_tag"]["type"] == "string"


def test_container_tag_default_names_a_published_image():
    """The default tag must exist in ghcr (.github/workflows/container.yml).

    CI publishes ``dev`` (main), ``pr-<n>``, ``sha-<short>`` and, on ``v*``
    tags, the PEP 440 version. A development manifest (``0.3.0rc2.dev0``) has no
    image of its own, so it defaults to ``dev``; a release or release
    candidate manifest defaults to its own version, which its tag publishes.
    """
    tag = _config_params()["container_tag"]
    version = _manifest_version()
    if version.is_devrelease:
        assert tag == "dev"
    else:
        assert tag == str(version)



def test_slurm_partition_param():
    props = _schema_properties()
    assert _config_params()["slurm_partition"] == "intel"
    assert props["slurm_partition"]["type"] == "string"


# ------------------------------------------------------- report params


def test_report_params_are_declared():
    """Pipeline-only: run-all has no MultiQC report (run-batch writes the cohort only)."""
    props = _schema_properties()
    config = _config_params()
    assert config["skip_multiqc"] is False
    assert props["skip_multiqc"]["type"] == "boolean"
    assert config["multiqc_title"] is None
    assert props["multiqc_title"]["type"] == "string"
    run_all = vars(_run_all_defaults())
    assert "skip_multiqc" not in run_all and "multiqc_title" not in run_all


def test_report_params_have_their_own_schema_group():
    schema = json.loads(SCHEMA.read_text())
    group = schema["$defs"]["report_options"]
    assert set(group["properties"]) == {"skip_multiqc", "multiqc_title"}
    assert {"$ref": "#/$defs/report_options"} in schema["allOf"]


# ------------------------------------------------------- BLAT scatter params


def test_blat_scatter_params_are_declared():
    """Pipeline-only: how the BLAT TE search is spread over tasks. Outputs do
    not depend on them, so run-all has no such options."""
    props = _schema_properties()
    config = _config_params()
    assert config["blat_scatter"] is True
    assert props["blat_scatter"]["type"] == "boolean"
    assert config["blat_chunks_per_task"] == 8
    assert props["blat_chunks_per_task"]["type"] == "integer"
    assert props["blat_chunks_per_task"]["minimum"] == 1
    assert props["blat_chunk_reads"]["type"] == "integer"
    assert props["blat_chunk_reads"]["minimum"] == 1
    assert props["blat_chunk_reads"].get("hidden") is True
    run_all = vars(_run_all_defaults())
    for name in ("blat_scatter", "blat_chunks_per_task", "blat_chunk_reads"):
        assert name not in run_all


def test_blat_chunk_reads_default_is_the_unsplit_chunk_size():
    """Chunks of RelocaTE2's 200,000 reads, as `relocaTE3 run` and te-split
    cut them (aligners.BLAT_CHUNK_SEQS)."""
    from RelocaTE3.aligners import BLAT_CHUNK_SEQS

    assert _config_params()["blat_chunk_reads"] == BLAT_CHUNK_SEQS
    te_split = build_parser().parse_args(["te-split", "-1", "r1.fq", "-n", "s"])
    assert te_split.chunk_reads == BLAT_CHUNK_SEQS


def test_blat_scatter_params_have_their_own_schema_group():
    schema = json.loads(SCHEMA.read_text())
    group = schema["$defs"]["scatter_options"]
    assert set(group["properties"]) == {"blat_scatter", "blat_chunks_per_task", "blat_chunk_reads"}
    assert {"$ref": "#/$defs/scatter_options"} in schema["allOf"]


def _main_nf_typed_params() -> dict[str, tuple[str, object]]:
    """main.nf's typed ``params {}`` block: name -> (type, default)."""
    body = _config_block((ROOT / "main.nf").read_text(), "params")
    typed = {}
    for line in body.splitlines():
        line = re.sub(r"(^|\s)//.*$", "", line).strip()
        if not line:
            continue
        m = re.fullmatch(r"(\w+)\s*:\s*(\w+)\s*=\s*(.+)", line)
        assert m, f"unsupported line in main.nf params block: {line}"
        typed[m.group(1)] = (m.group(2), _literal(m.group(3)))
    return typed


def test_every_boolean_param_is_typed_in_main_nf():
    """Nextflow >= 26.04 (strict parser) passes an untyped `--genotype false`
    as the String "false", which is true in Groovy. Only a param typed Boolean
    in main.nf's params block arrives as a Boolean. `help` is boolean or a
    parameter name, so it stays untyped (main.nf handles both)."""
    booleans = {name for name, prop in _schema_properties().items()
                if prop.get("type") == "boolean"}
    typed = _main_nf_typed_params()
    config = _config_params()
    assert booleans == {name for name, (kind, _) in typed.items() if kind == "Boolean"}
    for name in booleans:
        assert typed[name][1] == config[name], name
