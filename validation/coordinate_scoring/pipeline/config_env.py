#!/usr/bin/env python3.12
"""Config-to-shell bridge for the multi-dataset benchmark pipeline.

Canonical task order is dataset, caller, then panel-manifest row. Dataset
selection accepts one key, comma-separated keys, or ``full``. Legacy configs
with a single ``[dataset]`` table remain readable as dataset ``default``.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib.config import load_config
from lib.panel import panel_path, read_panel_manifest

ADAPTER_ENV_MAP: dict[str, dict[str, str]] = {
    "relocate3": {
        "repo": "RT3_REPO",
        "tsd": "TSD_PATTERN",
        "te_aligner": "RT3_TE_ALIGNER",
        "genome_aligner": "RT3_GENOME_ALIGNER",
        "require_both_junctions": "RT3_REQUIRE_BOTH_JUNCTIONS",
    },
    "relocate2": {
        "aligner": "RT2_ALIGNER",
        "size": "RT2_SIZE",
        "mismatch": "RT2_MISMATCH",
    },
}


def _sq(value) -> str:
    return "'" + str(value).replace("'", "'\\''") + "'"


def pretty_caller(key: str) -> str:
    m = re.match(r"^relocate3-([^-]+)-(.+)$", key)
    if m:
        return f"RelocaTE3-{m.group(1)}/{m.group(2)}"
    if key.lower().startswith("relocate"):
        return "RelocaTE" + key[len("relocate") :]
    return key


def _dataset_tables(cfg: dict) -> dict[str, dict]:
    if "datasets" in cfg:
        return cfg["datasets"]
    if "dataset" in cfg:
        return {"default": cfg["dataset"]}
    raise KeyError("config must contain [datasets.<key>] or legacy [dataset]")


def _enabled_datasets(cfg: dict) -> list[str]:
    return sorted(
        key
        for key, table in _dataset_tables(cfg).items()
        if table.get("enabled", True) is True
    )


def _default_dataset(cfg: dict) -> str:
    enabled = _enabled_datasets(cfg)
    if not enabled:
        raise ValueError("no enabled benchmark datasets")
    default = cfg.get("benchmark", {}).get("default_dataset", enabled[0])
    if default not in enabled:
        raise ValueError(
            f"default dataset {default!r} is not enabled; available: {', '.join(enabled)}"
        )
    return default


def _select_datasets(cfg: dict, selection: str | None) -> list[str]:
    enabled = _enabled_datasets(cfg)
    if selection is None:
        return [_default_dataset(cfg)]
    if selection.lower() in {"full", "all"}:
        return enabled
    requested = [item.strip() for item in selection.split(",") if item.strip()]
    unknown = sorted(set(requested) - set(enabled))
    if unknown:
        raise ValueError(
            f"unknown or disabled dataset(s): {', '.join(unknown)}; "
            f"available: {', '.join(enabled)}, full"
        )
    return requested


def _enabled_callers(cfg: dict) -> list[str]:
    return sorted(
        name
        for name, table in cfg.get("callers", {}).items()
        if table.get("enabled") is True
    )


def _manifest_rows(panel_root: Path) -> list[dict]:
    return read_panel_manifest(panel_root)


def _globals(cfg: dict) -> str:
    run = cfg["run"]
    lines = [
        f"WORK_ROOT={_sq(run['work_root'])}",
        f"REPORT_ROOT={_sq(run.get('report_root', 'reports'))}",
        f"TRUTH_ROOT={_sq(run.get('truth_root', 'truth'))}",
        f"THREADS={_sq(run['threads'])}",
        f"MATCH_WINDOW={_sq(cfg['scoring']['match_window'])}",
        f"SCORING_COORDINATE_POLICY={_sq(cfg['scoring'].get('coordinate_policy', 'legacy'))}",
        f"DEFAULT_DATASET={_sq(_default_dataset(cfg))}",
    ]
    # Preserve the old one-table shell contract for external scripts that still
    # use a legacy [dataset] config. Multi-dataset jobs use dataset-env.
    if "datasets" not in cfg:
        ds = cfg["dataset"]
        lines.extend(
            (
                f"REFERENCE={_sq(ds['reference'])}",
                f"TE_LIBRARY={_sq(ds['te_library'])}",
                f"REPEATMASKER={_sq(ds['repeatmasker'])}",
                f"TE_NAME={_sq(ds['te_name'])}",
                f"PANEL_ROOT={_sq(ds['panel_root'])}",
            )
        )
    return "\n".join(lines) + "\n"


def _dataset_env(cfg: dict, name: str) -> str:
    tables = _dataset_tables(cfg)
    if name not in tables:
        raise ValueError(f"unknown dataset: {name}")
    ds = tables[name]
    run = cfg["run"]
    work_root = run["work_root"]
    report_root = run.get("report_root", "reports")
    truth_root = run.get("truth_root", "truth")
    staged_te_library = str(Path("cache") / "te_libraries" / name / "library.fa")
    lines = [
        f"DATASET={_sq(name)}",
        f"DATASET_LABEL={_sq(ds.get('label', name))}",
        f"REFERENCE={_sq(ds['reference'])}",
        f"TE_LIBRARY_SOURCE={_sq(ds['te_library'])}",
        f"TE_LIBRARY={_sq(staged_te_library)}",
        f"REPEATMASKER={_sq(ds['repeatmasker'])}",
        f"REPEATMASKER_GFF={_sq(ds.get('repeatmasker_gff', ''))}",
        f"TE_NAME={_sq(ds['te_name'])}",
        f"PANEL_ROOT={_sq(ds['panel_root'])}",
        f"DATASET_WORK_ROOT={_sq(Path(work_root) / name)}",
        f"DATASET_REPORT_ROOT={_sq(Path(report_root) / 'datasets' / name)}",
        f"DATASET_TRUTH_ROOT={_sq(Path(truth_root) / name)}",
    ]
    return "\n".join(lines) + "\n"


def _datasets(cfg: dict, selection: str | None) -> str:
    tables = _dataset_tables(cfg)
    rows = [
        "\t".join((key, str(tables[key].get("label", key)), str(tables[key]["panel_root"])))
        for key in _select_datasets(cfg, selection)
    ]
    return "".join(row + "\n" for row in rows)


def _callers(cfg: dict) -> str:
    return "".join(f"{name}\n" for name in _enabled_callers(cfg))


def _adapter_key(table: dict) -> str:
    return Path(table.get("adapter", "")).name


def _caller_env(cfg: dict, name: str) -> str:
    table = cfg.get("callers", {}).get(name, {})
    mapping = ADAPTER_ENV_MAP.get(_adapter_key(table))
    if not mapping:
        return ""
    lines = [
        f"{env}={_sq(table[key])}"
        for key, env in mapping.items()
        if key in table
    ]
    return ("\n".join(lines) + "\n") if lines else ""


def _adapter(cfg: dict, name: str) -> str:
    table = cfg.get("callers", {}).get(name, {})
    return (table.get("adapter") or f"callers/{name}") + "\n"


def _labels(cfg: dict) -> str:
    callers = cfg.get("callers", {})
    rows = [
        f"{name}\t{callers.get(name, {}).get('label') or pretty_caller(name)}"
        for name in _enabled_callers(cfg)
    ]
    return "\n".join(rows) + "\n" if rows else ""


def _task_records(cfg: dict, dataset_selection: str | None = None) -> list[dict]:
    records = []
    tables = _dataset_tables(cfg)
    for dataset in _select_datasets(cfg, dataset_selection):
        panel_root = Path(tables[dataset]["panel_root"])
        rows = _manifest_rows(panel_root)
        for caller in _enabled_callers(cfg):
            for row in rows:
                records.append(
                    {
                        "dataset": dataset,
                        "caller": caller,
                        "sample": row["sample"],
                        "coverage": row["coverage"],
                        "replicate": row["replicate"],
                        "divergence": row.get("divergence_percent", ""),
                        "r1": str(panel_path(panel_root, row["r1"])),
                        "r2": str(panel_path(panel_root, row["r2"])),
                    }
                )
    return records


def _tasks(cfg: dict, selection: str | None) -> str:
    fields = (
        ("dataset", "caller", "sample", "coverage", "replicate", "r1", "r2")
        if "datasets" in cfg
        else ("caller", "sample", "coverage", "replicate", "r1", "r2")
    )
    return "".join(
        "\t".join(record[field] for field in fields) + "\n"
        for record in _task_records(cfg, selection)
    )


def _parse_filter(value: str | None) -> set[str] | None:
    if value is None:
        return None
    return {item.strip() for item in value.split(",") if item.strip()}


def _indices(cfg: dict, dataset_selection: str | None, filters: dict) -> str:
    wanted = {key: _parse_filter(value) for key, value in filters.items()}
    matched = [
        index
        for index, record in enumerate(_task_records(cfg, dataset_selection))
        if all(
            allowed is None or record[key] in allowed
            for key, allowed in wanted.items()
        )
    ]
    return ",".join(str(index) for index in matched) + "\n"


def run(argv=None) -> str:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/benchmark.toml", type=Path)
    parser.add_argument(
        "--dataset",
        help="dataset key, comma-separated keys, or full (default: config default)",
    )
    parser.add_argument(
        "mode",
        choices=(
            "globals",
            "datasets",
            "dataset-env",
            "callers",
            "caller-env",
            "adapter",
            "labels",
            "tasks",
            "count",
            "indices",
        ),
    )
    parser.add_argument("name", nargs="?", help="dataset or caller name")
    parser.add_argument("--caller")
    parser.add_argument("--coverage")
    parser.add_argument("--divergence")
    parser.add_argument("--replicate")
    parser.add_argument("--sample")
    args = parser.parse_args(argv)
    cfg = load_config(args.config)

    if args.mode == "globals":
        return _globals(cfg)
    if args.mode == "datasets":
        return _datasets(cfg, args.dataset)
    if args.mode == "dataset-env":
        if not args.name:
            parser.error("dataset-env requires a dataset name")
        return _dataset_env(cfg, args.name)
    if args.mode == "callers":
        return _callers(cfg)
    if args.mode == "caller-env":
        if not args.name:
            parser.error("caller-env requires a caller name")
        return _caller_env(cfg, args.name)
    if args.mode == "adapter":
        if not args.name:
            parser.error("adapter requires a caller name")
        return _adapter(cfg, args.name)
    if args.mode == "labels":
        return _labels(cfg)
    if args.mode == "tasks":
        return _tasks(cfg, args.dataset)
    if args.mode == "count":
        return f"{len(_task_records(cfg, args.dataset))}\n"
    if args.mode == "indices":
        return _indices(
            cfg,
            args.dataset,
            {
                "caller": args.caller,
                "coverage": args.coverage,
                "divergence": args.divergence,
                "replicate": args.replicate,
                "sample": args.sample,
            },
        )
    raise AssertionError("unreachable")


def main(argv=None) -> int:
    try:
        sys.stdout.write(run(argv))
    except (KeyError, ValueError, FileNotFoundError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
