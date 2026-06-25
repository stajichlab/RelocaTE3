"""Normalize legacy RelocaTE2 characterized insertion tables into a single TSV.

Walks ``paths.relocate2_results`` looking for
``<sample>/repeat/results/ALL.all_nonref_insert.characTErized.txt`` files and
flattens them into the schema consumed by ``compare_char.py``.

Input columns (RelocaTE2 + RelocaTE3 share an identical 8-col header):

    strain  TE  TSD  chromosome.pos  strand  avg_flankers  spanners  status

Output columns (see ``_common.CHAR_COLUMNS``):

    sample, chrom, start, end, strand, te_name, tsd, avg_flankers, spanners,
    status, source_file
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import CHAR_COLUMNS as COLUMNS  # noqa: E402
from _config import load_config, load_samples  # noqa: E402

LEGACY_CHAR_RELPATH = Path("repeat/results/ALL.all_nonref_insert.characTErized.txt")
_COORD_RE = re.compile(r"^([^:]+):(\d+)\.\.(\d+)$")


def parse_char_txt(path: Path, sample: str):
    """Yield normalized records from one characTErized.txt file."""
    with open(path) as fh:
        header = fh.readline()  # skip header
        if not header.lower().startswith("strain"):
            # File starts with data; rewind by reseeking.
            fh.seek(0)
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 8:
                continue
            strain, te_name, tsd, chrom_pos, strand, flank, span, status = f[:8]
            m = _COORD_RE.match(chrom_pos)
            if not m:
                continue
            chrom = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            yield {
                "sample": sample,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "te_name": te_name,
                "tsd": tsd,
                "avg_flankers": flank,
                "spanners": span,
                "status": status,
                "source_file": str(path),
            }
            _ = strain  # already supplied as ``sample``


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.toml", help="Validation TOML config")
    ap.add_argument(
        "--out",
        help="Output TSV (default: <report_dir>/characterized/relocate2_calls.tsv)",
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
    report_dir = Path(cfg["paths"]["report_dir"]) / "characterized"
    report_dir.mkdir(parents=True, exist_ok=True)
    out_path = Path(args.out) if args.out else report_dir / "relocate2_calls.tsv"

    n_total = 0
    n_missing = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for sample in samples:
            txt = r2_root / sample / LEGACY_CHAR_RELPATH
            if not txt.is_file():
                print(
                    f"WARN: missing legacy characterized TXT for {sample}: {txt}",
                    file=sys.stderr,
                )
                n_missing += 1
                continue
            count = 0
            for rec in parse_char_txt(txt, sample):
                writer.writerow(rec)
                count += 1
            n_total += count
            print(f"  {sample}: {count} RelocaTE2 characterized insertions", file=sys.stderr)

    print(
        f"Wrote {n_total} calls across {len(samples) - n_missing} samples -> {out_path}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
