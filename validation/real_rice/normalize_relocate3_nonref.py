"""Normalize RelocaTE3 per-sample non-reference calls into a TSV.

Reads ``<outdir>/<sample>/results/<target>.<te_name>.all_nonref_insert.gff``,
RelocaTE3's copy of RelocaTE2's headline tier (clean_false_positive.py: minus
calls at reference-TE boundaries and the singleton / insufficient_data /
supporting_reads classes). ``normalize_relocate2_nonref.py`` reads RelocaTE2's
``ALL.all_nonref_insert.gff``, the same tier, so ``compare_nonref.py`` compares
like with like. The unfiltered ``.txt`` is not comparable to it.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import NONREF_COLUMNS as COLUMNS  # noqa: E402
from _common import parse_nonref_gff  # noqa: E402
from _config import load_config, load_samples  # noqa: E402

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
    fname = f"{target}.{te_name}.all_nonref_insert.gff"

    n_total = 0
    n_missing = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for sample in samples:
            gff = r3_root / sample / "results" / fname
            if not gff.is_file():
                print(
                    f"WARN: RelocaTE3 GFF not found for {sample}: {gff}",
                    file=sys.stderr,
                )
                n_missing += 1
                continue
            count = 0
            for rec in parse_nonref_gff(gff, sample):
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
