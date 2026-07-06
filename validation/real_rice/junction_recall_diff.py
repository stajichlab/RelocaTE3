"""Per-site junction-read recall diff (R3 minimap2 vs R2 bwa-mem).

Joins R2 and R3 ``all_nonref_insert.txt`` files by (sample, chrom, +-window bp)
and prints sites where R2 has substantially more junction reads. Output:

  sample  chrom  pos  r2_T  r3_T  ratio

Used by ``plans/2026-06-26-genotype-status-parity.md`` to measure whether a
genome-align tuning change closes the recall gap that drives the residual
hom/excision -> het and het -> somatic status-confusion mismatches.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _config import load_config, load_samples  # noqa: E402

_COUNT = re.compile(r"T:(\d+)")


def parse_nonref(path: Path):
    """Yield (chrom, start_pos, total) tuples from an all_nonref_insert.txt."""
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or not cols[3].startswith("Chr"):
                continue
            chrom = cols[3]
            m_coord = re.match(r"(\d+)\.\.(\d+)", cols[4])
            if not m_coord:
                continue
            pos = int(m_coord.group(1))
            m_total = _COUNT.search(cols[6])
            if not m_total:
                continue
            total = int(m_total.group(1))
            yield chrom, pos, total


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.toml", help="Validation TOML config")
    ap.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help="Restrict to these samples (default: config sample CSV)",
    )
    ap.add_argument(
        "--window",
        type=int,
        default=5,
        help="+-bp window for site matching (default 5)",
    )
    ap.add_argument(
        "--min-ratio",
        type=float,
        default=2.0,
        help="Report sites where r2_T/r3_T >= this ratio (default 2.0)",
    )
    ap.add_argument("--out", default="-", help="Output TSV path or - for stdout")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    samples = args.samples or load_samples(cfg["paths"]["sample_csv"])
    r2_root = Path(cfg["paths"]["relocate2_results"])
    r3_root = Path(cfg["paths"]["relocate3_outdir"])
    te_name = cfg["relocate3"]["te_name"]
    target = cfg["relocate3"]["target"]
    r3_fname = f"{target}.{te_name}.all_nonref_insert.txt"

    rows: list[dict] = []
    n_r2_sites = n_r3_sites = 0
    for sample in samples:
        r2 = r2_root / sample / "repeat" / "results" / "ALL.all_nonref_insert.txt"
        r3 = r3_root / sample / "results" / r3_fname
        if not r2.exists() or not r3.exists():
            print(
                f"WARN: missing for {sample}: r2={r2.exists()} r3={r3.exists()}",
                file=sys.stderr,
            )
            continue
        r2_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for chrom, pos, total in parse_nonref(r2):
            r2_by_chrom[chrom].append((pos, total))
            n_r2_sites += 1
        r3_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for chrom, pos, total in parse_nonref(r3):
            r3_by_chrom[chrom].append((pos, total))
            n_r3_sites += 1
        for chrom, items in r2_by_chrom.items():
            r3_index = sorted(r3_by_chrom.get(chrom, []))
            for r2_pos, r2_total in items:
                r3_match: tuple[int, int] | None = None
                for r3_pos, r3_total in r3_index:
                    if abs(r3_pos - r2_pos) <= args.window:
                        if r3_match is None or abs(r3_pos - r2_pos) < abs(
                            r3_match[0] - r2_pos
                        ):
                            r3_match = (r3_pos, r3_total)
                if r3_match is None:
                    continue
                r3_total = r3_match[1]
                if r3_total == 0:
                    ratio = float("inf")
                else:
                    ratio = r2_total / r3_total
                if ratio >= args.min_ratio:
                    rows.append(
                        {
                            "sample": sample,
                            "chrom": chrom,
                            "pos": r2_pos,
                            "r2_T": r2_total,
                            "r3_T": r3_total,
                            "ratio": f"{ratio:.2f}" if ratio != float("inf") else "inf",
                        }
                    )

    out = sys.stdout if args.out == "-" else open(args.out, "w", newline="")
    w = csv.DictWriter(
        out,
        fieldnames=["sample", "chrom", "pos", "r2_T", "r3_T", "ratio"],
        delimiter="\t",
    )
    w.writeheader()
    w.writerows(rows)
    if args.out != "-":
        out.close()
    print(
        f"R2 sites: {n_r2_sites}  R3 sites: {n_r3_sites}  "
        f"low-recall sites (r2_T/r3_T >= {args.min_ratio}): {len(rows)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    raise SystemExit(main())
