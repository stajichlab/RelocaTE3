"""Normalize legacy RelocaTE2 GFFs into a single TSV.

Walks ``paths.relocate2_results`` looking for
``<sample>/repeat/results/ALL.all_nonref_insert.gff`` files and produces a
flat TSV used by ``compare_nonref.py``.

Output columns:
    sample, chrom, start, end, strand, te_name, tsd, left_junction_reads,
    right_junction_reads, left_support_reads, right_support_reads, source_file
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import NONREF_COLUMNS as COLUMNS  # noqa: E402
from _config import load_config, load_samples  # noqa: E402

ATTR_RE = re.compile(r"([A-Za-z_]+)=([^;]*)")
LEGACY_GFF_RELPATH = Path("repeat/results/ALL.all_nonref_insert.gff")


def _attrs(field: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(field))


def parse_gff(path: Path, sample: str):
    """Yield normalized records from one ALL.all_nonref_insert.gff."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = _attrs(f[8])
            yield {
                "sample": sample,
                "chrom": f[0],
                "start": int(f[3]),
                "end": int(f[4]),
                "strand": f[6],
                "te_name": attrs.get("Name", "NA"),
                "tsd": attrs.get("TSD", "NA"),
                "left_junction_reads": int(attrs.get("Left_junction_reads", "0")),
                "right_junction_reads": int(attrs.get("Right_junction_reads", "0")),
                "left_support_reads": int(attrs.get("Left_support_reads", "0")),
                "right_support_reads": int(attrs.get("Right_support_reads", "0")),
                "source_file": str(path),
            }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.toml", help="Validation TOML config")
    ap.add_argument(
        "--out",
        help="Output TSV (default: <report_dir>/nonref/relocate2_calls.tsv)",
    )
    ap.add_argument(
        "--samples",
        action="append",
        default=None,
        help="Restrict to these sample names (repeatable). Default: read sample CSV.",
    )
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    samples = args.samples if args.samples else load_samples(cfg["paths"]["sample_csv"])
    r2_root = Path(cfg["paths"]["relocate2_results"])
    report_dir = Path(cfg["paths"]["report_dir"]) / "nonref"
    report_dir.mkdir(parents=True, exist_ok=True)
    out_path = Path(args.out) if args.out else report_dir / "relocate2_calls.tsv"

    n_total = 0
    n_missing = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for sample in samples:
            gff = r2_root / sample / LEGACY_GFF_RELPATH
            if not gff.is_file():
                print(f"WARN: missing legacy GFF for {sample}: {gff}", file=sys.stderr)
                n_missing += 1
                continue
            count = 0
            for rec in parse_gff(gff, sample):
                writer.writerow(rec)
                count += 1
            n_total += count
            print(f"  {sample}: {count} RelocaTE2 insertions", file=sys.stderr)

    print(
        f"Wrote {n_total} calls across {len(samples) - n_missing} samples -> {out_path}",
        file=sys.stderr,
    )
    return 0 if n_missing == 0 else 0  # missing samples warn but don't fail


if __name__ == "__main__":
    raise SystemExit(main())
