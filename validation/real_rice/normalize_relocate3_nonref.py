"""Normalize RelocaTE3 per-sample non-reference insertion tables into a TSV.

RelocaTE3's ``find-insertions`` writes one tab-delimited file per sample at
``<outdir>/results/<target>.<te_name>.all_nonref_insert.txt`` (see
src/RelocaTE3/insertions.py). Columns:

    family  tsd  sample  chrom  coor_start..coor_end  strand  T:N  R:N  L:N  ST:N  SR:N  SL:N

This script flattens those files into the same schema used by
``normalize_relocate2_nonref.py`` so ``compare_nonref.py`` can diff them directly.
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

_COORD_RE = re.compile(r"^(\d+)\.\.(\d+)$")


def _parse_count(field: str) -> int:
    """Parse a 'PREFIX:N' field into the integer N (returns 0 on malformed input)."""
    if ":" not in field:
        return 0
    try:
        return int(field.split(":", 1)[1])
    except ValueError:
        return 0


def parse_nonref_txt(path: Path, sample: str):
    """Yield normalized records from a single all_nonref_insert.txt file."""
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 12:
                continue
            family, tsd, _sample, chrom, coord, strand = f[0], f[1], f[2], f[3], f[4], f[5]
            m = _COORD_RE.match(coord)
            if not m:
                continue
            start, end = int(m.group(1)), int(m.group(2))
            total = _parse_count(f[6])     # T:N — total junction reads
            right = _parse_count(f[7])     # R:N — right-junction reads
            left = _parse_count(f[8])      # L:N — left-junction reads
            sup_right = _parse_count(f[10])  # SR:N — right-support reads
            sup_left = _parse_count(f[11])   # SL:N — left-support reads
            yield {
                "sample": sample,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "te_name": family,
                "tsd": tsd,
                "left_junction_reads": left,
                "right_junction_reads": right,
                "left_support_reads": sup_left,
                "right_support_reads": sup_right,
                "source_file": str(path),
            }
            _ = total  # unused; kept for clarity / future use


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.toml", help="Validation TOML config")
    ap.add_argument(
        "--out",
        help="Output TSV (default: <report_dir>/nonref/relocate3_calls.tsv)",
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
    r3_root = Path(cfg["paths"]["relocate3_outdir"])
    report_dir = Path(cfg["paths"]["report_dir"]) / "nonref"
    report_dir.mkdir(parents=True, exist_ok=True)
    out_path = Path(args.out) if args.out else report_dir / "relocate3_calls.tsv"

    te_name = cfg["relocate3"]["te_name"]
    target = cfg["relocate3"]["target"]
    fname = f"{target}.{te_name}.all_nonref_insert.txt"

    n_total = 0
    n_missing = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for sample in samples:
            txt = r3_root / sample / "results" / fname
            if not txt.is_file():
                print(
                    f"WARN: RelocaTE3 table not found for {sample}: {txt}",
                    file=sys.stderr,
                )
                n_missing += 1
                continue
            count = 0
            for rec in parse_nonref_txt(txt, sample):
                writer.writerow(rec)
                count += 1
            n_total += count
            print(f"  {sample}: {count} RelocaTE3 insertions", file=sys.stderr)

    print(
        f"Wrote {n_total} calls across {len(samples) - n_missing} samples -> {out_path}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
