"""Normalize RelocaTE3 per-sample characterized insertion tables into a TSV.

RelocaTE3's ``characterize`` step writes one tab-delimited file per sample at
``<outdir>/<sample>/results/<target>.<te_name>.all_nonref_insert.characTErized.txt``
(see src/RelocaTE3/characterize.py). Schema matches RelocaTE2's
characterized TXT exactly, so this script reuses the same parser as
``normalize_relocate2_char.py``.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import CHAR_COLUMNS as COLUMNS  # noqa: E402
from _config import load_config, load_samples  # noqa: E402
from normalize_relocate2_char import parse_char_txt  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.toml", help="Validation TOML config")
    ap.add_argument(
        "--out",
        help="Output TSV (default: <report_dir>/characterized/relocate3_calls.tsv)",
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
    report_dir = Path(cfg["paths"]["report_dir"]) / "characterized"
    report_dir.mkdir(parents=True, exist_ok=True)
    out_path = Path(args.out) if args.out else report_dir / "relocate3_calls.tsv"

    te_name = cfg["relocate3"]["te_name"]
    target = cfg["relocate3"]["target"]
    fname = f"{target}.{te_name}.all_nonref_insert.characTErized.txt"

    n_total = 0
    n_missing = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for sample in samples:
            txt = r3_root / sample / "results" / fname
            if not txt.is_file():
                print(
                    f"WARN: RelocaTE3 characterized TXT not found for {sample}: {txt}",
                    file=sys.stderr,
                )
                n_missing += 1
                continue
            count = 0
            for rec in parse_char_txt(txt, sample):
                writer.writerow(rec)
                count += 1
            n_total += count
            print(f"  {sample}: {count} RelocaTE3 characterized insertions", file=sys.stderr)

    print(
        f"Wrote {n_total} calls across {len(samples) - n_missing} samples -> {out_path}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
