#!/usr/bin/env python3
"""Real-rice sample sheet and calling parameters for the Nextflow pipeline.

Reads the same config as ``run_relocate3.sh`` (``config.toml`` next to this
script, else ``config.example.toml``), so the pipeline run and the run-all run
cannot drift apart.

Modes:
  -o SHEET          write ``sample,fastq_1,fastq_2,genotype_bam`` (absolute
                    paths; ``-`` = stdout). Every FASTQ, genotype alignment and
                    its index (.crai for CRAM; .bai or .csi for BAM, as the
                    pipeline looks for them) must exist; all missing files are
                    reported together. An existing SHEET is kept when
                    identical and refused when different (``--force`` replaces
                    it).
  --print-params    print ``<pipeline param>\\t<value>`` for every calling
                    setting run_relocate3.sh passes to run-all.
  --print-samples   print the sample names, one per line.

Relative input paths in [paths] resolve against the config's directory, as in
run_relocate3.sh, which runs from there. ``RELOCATE3_VALIDATION_ROOT``
overrides ``paths.validation_root`` (_config.py).
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _config import load_config, load_samples  # noqa: E402

HEADER = ["sample", "fastq_1", "fastq_2", "genotype_bam"]

# Input keys under [paths] resolved against the config directory.
_INPUT_KEYS = ("genome", "te_library", "repeatmasker", "fq_dir", "sample_csv", "genome_aln_dir")

# (pipeline param, [relocate3] key): every setting run_relocate3.sh passes to
# `relocaTE3 run-all`, under the pipeline's name for it.
_PARAMS = (
    ("te_aligner", "aligner"),
    ("genome_aligner", "genome_aligner"),
    ("mismatch", "mismatch"),
    ("min_match", "min_match"),
    ("min_trimmed", "min_trimmed"),
    ("min_mapq", "min_mapq"),
    ("size", "insert_size"),
    ("te_name", "te_name"),
    ("target", "target"),
    ("tsd", "tsd"),
)


def default_config(harness_dir: Path) -> Path:
    """config.toml if present (a local edit), else the tracked example."""
    local = harness_dir / "config.toml"
    return local if local.is_file() else harness_dir / "config.example.toml"


def load(path: Path) -> dict[str, Any]:
    cfg = load_config(path)
    base = Path(cfg["_config_dir"])
    paths = cfg["paths"]
    for key in _INPUT_KEYS:
        if paths.get(key):
            paths[key] = os.path.normpath(base / paths[key])  # absolute paths unchanged
    return cfg


def _genotype_enabled(cfg: dict[str, Any]) -> bool:
    r3 = cfg["relocate3"]
    if not bool(r3.get("characterize", False)):
        return False
    source = r3.get("characterize_source", "external")
    if source != "external":
        raise SystemExit(
            f"ERROR: [relocate3].characterize_source = {source!r}; only 'external' "
            "(paths.genome_aln_dir / genome_aln_pattern) is supported, as in run_relocate3.sh"
        )
    return True


def pipeline_params(cfg: dict[str, Any]) -> list[tuple[str, str]]:
    r3 = cfg["relocate3"]
    out = [(name, str(r3[key])) for name, key in _PARAMS]
    out.append(("genotype", "true" if _genotype_enabled(cfg) else "false"))
    out += [(k, cfg["paths"][k]) for k in ("genome", "te_library", "repeatmasker")]
    return out


def _index_candidates(aln: str) -> list[str]:
    # As workflows/relocate3.nf genotypeBamIndex().
    return [aln + ".crai"] if aln.endswith(".cram") else [aln + ".bai", aln + ".csi"]


def sheet_rows(cfg: dict[str, Any]) -> list[list[str]]:
    """Rows for every sample; raises SystemExit listing every missing file."""
    paths, reads = cfg["paths"], cfg["reads"]
    genotype = _genotype_enabled(cfg)
    rows, problems = [], []
    for s in load_samples(paths["sample_csv"]):
        r1 = os.path.join(paths["fq_dir"], s + reads["r1_suffix"])
        r2 = os.path.join(paths["fq_dir"], s + reads["r2_suffix"])
        for fq in (r1, r2):
            if not os.path.isfile(fq):
                problems.append(f"{s}: FASTQ not found: {fq}")
        aln = ""
        if genotype:
            aln = os.path.join(paths["genome_aln_dir"], paths["genome_aln_pattern"].replace("${sample}", s))
            if not os.path.isfile(aln):
                problems.append(f"{s}: genotype alignment not found: {aln}")
            elif not any(os.path.isfile(i) for i in _index_candidates(aln)):
                problems.append(
                    f"{s}: genotype alignment has no index; expected "
                    f"{' or '.join(_index_candidates(aln))} (create it with `samtools index`)"
                )
        rows.append([s, r1, r2, aln])
    if problems:
        raise SystemExit("ERROR: sample sheet inputs missing:\n  " + "\n  ".join(problems))
    return rows


def render(rows: list[list[str]]) -> str:
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(HEADER)
    w.writerows(rows)
    return buf.getvalue()


def write_sheet(text: str, out: Path, force: bool) -> str:
    if out.exists():
        if out.read_text() == text:
            return f"sample sheet unchanged: {out}"
        if not force:
            raise SystemExit(
                f"ERROR: {out} exists and differs from the sheet this config gives; "
                "the inputs changed. Use a new output directory, or --force to replace it."
            )
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=out.parent, prefix=f".{out.name}.", delete=False) as tmp:
        tmp.write(text)
    os.replace(tmp.name, out)
    return f"sample sheet written: {out}"


def main(argv: list[str] | None = None) -> int:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--config", type=Path, default=default_config(here),
                    help="real-rice config (default: config.toml, else config.example.toml, next to this script)")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("-o", "--output", help="sample sheet to write ('-' = stdout)")
    mode.add_argument("--print-params", action="store_true", help="print pipeline params from the config")
    mode.add_argument("--print-samples", action="store_true", help="print sample names")
    ap.add_argument("--force", action="store_true", help="replace an existing, different sheet")
    args = ap.parse_args(argv)

    cfg = load(args.config)
    if args.print_params:
        for name, value in pipeline_params(cfg):
            print(f"{name}\t{value}")
    elif args.print_samples:
        print("\n".join(load_samples(cfg["paths"]["sample_csv"])))
    else:
        text = render(sheet_rows(cfg))
        if args.output == "-":
            sys.stdout.write(text)
        else:
            print(write_sheet(text, Path(args.output), args.force), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
