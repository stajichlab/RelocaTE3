"""Shared config loader for the real-rice validation pipeline.

Reads a TOML file, expands ${section.key} interpolations, and resolves any
relative output paths against the directory containing the config file so
collaborators can run the pipeline from anywhere on the cluster after a fresh
``git pull``.
"""

from __future__ import annotations

import csv
import os
import re
import sys
from pathlib import Path
from typing import Any

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - Py3.10 fallback
    import tomli as tomllib

_INTERP = re.compile(r"\$\{([a-zA-Z_][a-zA-Z0-9_]*)\.([a-zA-Z_][a-zA-Z0-9_]*)\}")

# Output keys under [paths] that are resolved relative to the config dir if not absolute.
_RELATIVE_OUTPUT_KEYS = ("relocate3_outdir", "report_dir", "log_dir")


def _expand(value: str, cfg: dict[str, Any]) -> str:
    """Resolve ${section.key} references; raise on unknown keys."""
    seen = 0
    while True:
        m = _INTERP.search(value)
        if not m:
            return value
        section, key = m.group(1), m.group(2)
        try:
            replacement = cfg[section][key]
        except KeyError as exc:
            raise KeyError(
                f"Config interpolation '${{{section}.{key}}}' refers to missing key"
            ) from exc
        if not isinstance(replacement, str):
            raise TypeError(
                f"Cannot interpolate non-string value at {section}.{key}"
            )
        value = value[: m.start()] + replacement + value[m.end() :]
        seen += 1
        if seen > 64:
            raise RuntimeError(f"Interpolation loop in value: {value!r}")


def _expand_all(cfg: dict[str, Any]) -> dict[str, Any]:
    for section, body in cfg.items():
        if not isinstance(body, dict):
            continue
        for key, val in list(body.items()):
            if isinstance(val, str):
                body[key] = _expand(val, cfg)
    return cfg


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a validation TOML config and resolve interpolations + relative paths.

    Honors a few environment variable overrides so collaborators don't need
    to edit ``config.toml`` just because their copy of the bundled
    ``validation_data/`` lives in a different absolute path:

    - ``RELOCATE3_VALIDATION_ROOT`` overrides ``[paths].validation_root``.
    """
    config_path = Path(path).resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "rb") as fh:
        cfg = tomllib.load(fh)
    env_root = os.environ.get("RELOCATE3_VALIDATION_ROOT")
    if env_root:
        cfg.setdefault("paths", {})["validation_root"] = env_root
    cfg = _expand_all(cfg)

    config_dir = config_path.parent
    paths = cfg.setdefault("paths", {})
    for key in _RELATIVE_OUTPUT_KEYS:
        val = paths.get(key)
        if val and not os.path.isabs(val):
            paths[key] = str((config_dir / val).resolve())

    cfg["_config_path"] = str(config_path)
    cfg["_config_dir"] = str(config_dir)
    return cfg


def resolve_genome_aln(cfg: dict[str, Any], sample: str) -> Path:
    """Resolve the per-sample reads-to-genome alignment path.

    Substitutes ``${sample}`` into ``paths.genome_aln_pattern`` and joins it
    under ``paths.genome_aln_dir``. Raises if either key is missing or if the
    resulting file does not exist.
    """
    paths = cfg.get("paths", {})
    aln_dir = paths.get("genome_aln_dir")
    pattern = paths.get("genome_aln_pattern")
    if not aln_dir or not pattern:
        raise KeyError(
            "config is missing paths.genome_aln_dir / paths.genome_aln_pattern "
            "(needed when [relocate3].characterize = true and "
            "characterize_source = 'external')"
        )
    rel = pattern.replace("${sample}", sample)
    aln_path = Path(aln_dir) / rel
    if not aln_path.is_file():
        raise FileNotFoundError(
            f"genome alignment not found for sample {sample!r}: {aln_path}"
        )
    return aln_path


def load_samples(csv_path: str | Path) -> list[str]:
    """Read sample names from a one-column CSV with a 'Sample' header."""
    samples: list[str] = []
    with open(csv_path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header or header[0].strip().lower() not in {"sample", "strain"}:
            raise ValueError(
                f"Expected first column header 'Sample' or 'Strain' in {csv_path}"
            )
        for row in reader:
            if row and row[0].strip():
                samples.append(row[0].strip())
    return samples
