"""Compare RelocaTE3 vs legacy RelocaTE2 non-reference insertion calls.

Inputs:
    Two normalized TSVs produced by normalize_relocate2_nonref.py and
    normalize_relocate3_nonref.py (same schema).

Matching rule (per-sample):
    Two calls match if they are on the same chromosome, name the same TE
    family, and their start coordinates are within ``compare.position_window``
    bp. RelocaTE2 calls are matched greedily to the nearest unmatched
    RelocaTE3 call within the window.

Outputs (under <report_dir>/nonref/):
    summary.tsv          per-sample shared / R2-only / R3-only counts + recall + Jaccard
    summary.txt          plain-text overall stats (human readable)
    matched_calls.tsv    every R2 call with its matched R3 call (or blank if unmatched)
    relocate2_only.tsv   R2 calls with no matching R3 call within the window
    relocate3_only.tsv   R3 calls with no matching R2 call within the window
    venn_total.png       overall Venn across every sample compared
    venn/venn_<sample>.png  one Venn per sample (only if matplotlib_venn imports)
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (  # noqa: E402
    draw_one_venn,
    filter_nonref,
    load_venn,
    match_sample,
    read_tsv,
    write_nonref_pairs,
    write_rows,
)
from _config import load_config  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.toml", help="Validation TOML config")
    ap.add_argument("--r2-tsv", help="RelocaTE2 calls TSV (default: <report>/nonref/relocate2_calls.tsv)")
    ap.add_argument("--r3-tsv", help="RelocaTE3 calls TSV (default: <report>/nonref/relocate3_calls.tsv)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    report_dir = Path(cfg["paths"]["report_dir"]) / "nonref"
    report_dir.mkdir(parents=True, exist_ok=True)

    r2_tsv = Path(args.r2_tsv) if args.r2_tsv else report_dir / "relocate2_calls.tsv"
    r3_tsv = Path(args.r3_tsv) if args.r3_tsv else report_dir / "relocate3_calls.tsv"

    if not r2_tsv.is_file():
        print(
            f"ERROR: missing {r2_tsv} (run normalize_relocate2_nonref.py first)",
            file=sys.stderr,
        )
        return 1
    if not r3_tsv.is_file():
        print(
            f"ERROR: missing {r3_tsv} (run normalize_relocate3_nonref.py first)",
            file=sys.stderr,
        )
        return 1

    window = int(cfg["compare"]["position_window"])
    te_family = cfg["compare"].get("te_family", "")
    min_jr = int(cfg["compare"].get("min_junction_reads", 1))

    r2_all = filter_nonref(read_tsv(r2_tsv), te_family, min_jr)
    r3_all = filter_nonref(read_tsv(r3_tsv), te_family, min_jr)

    by_sample_r2: dict[str, list[dict]] = defaultdict(list)
    by_sample_r3: dict[str, list[dict]] = defaultdict(list)
    for rec in r2_all:
        by_sample_r2[rec["sample"]].append(rec)
    for rec in r3_all:
        by_sample_r3[rec["sample"]].append(rec)

    samples = sorted(set(by_sample_r2) | set(by_sample_r3))

    all_matched: list[tuple[dict, dict]] = []
    all_r2_only: list[dict] = []
    all_r3_only: list[dict] = []
    summary_rows: list[dict] = []
    per_sample_counts: list[tuple[str, int, int, int]] = []
    total_shared = total_r2_only = total_r3_only = 0

    plt, venn2 = load_venn()
    venn_dir = report_dir / "venn"
    if plt is not None:
        if venn_dir.is_dir():
            for stale in venn_dir.glob("venn_*.png"):
                stale.unlink()
        venn_dir.mkdir(parents=True, exist_ok=True)

    for s in samples:
        matched, r2_only, r3_only = match_sample(
            by_sample_r2.get(s, []), by_sample_r3.get(s, []), window
        )
        all_matched.extend(matched)
        all_r2_only.extend(r2_only)
        all_r3_only.extend(r3_only)
        per_sample_counts.append((s, len(matched), len(r2_only), len(r3_only)))

        r2_n = len(by_sample_r2.get(s, []))
        r3_n = len(by_sample_r3.get(s, []))
        shared = len(matched)
        union = shared + len(r2_only) + len(r3_only)
        recall = shared / r2_n if r2_n else 0.0
        precision = shared / r3_n if r3_n else 0.0
        jaccard = shared / union if union else 0.0

        total_shared += shared
        total_r2_only += len(r2_only)
        total_r3_only += len(r3_only)

        summary_rows.append(
            {
                "sample": s,
                "relocate2_total": r2_n,
                "relocate3_total": r3_n,
                "shared": shared,
                "relocate2_only": len(r2_only),
                "relocate3_only": len(r3_only),
                "recall_vs_r2": f"{recall:.4f}",
                "precision_vs_r2": f"{precision:.4f}",
                "jaccard": f"{jaccard:.4f}",
            }
        )

    summary_tsv = report_dir / "summary.tsv"
    with open(summary_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    write_nonref_pairs(report_dir / "matched_calls.tsv", all_matched, all_r2_only)
    write_rows(report_dir / "relocate2_only.tsv", all_r2_only)
    write_rows(report_dir / "relocate3_only.tsv", all_r3_only)

    r2_total = sum(len(by_sample_r2.get(s, [])) for s in samples)
    r3_total = sum(len(by_sample_r3.get(s, [])) for s in samples)
    union = total_shared + total_r2_only + total_r3_only
    recall = total_shared / r2_total if r2_total else 0.0
    precision = total_shared / r3_total if r3_total else 0.0
    jaccard = total_shared / union if union else 0.0

    lines = [
        "RelocaTE3 vs RelocaTE2 comparison (non-reference / raw calls)",
        "=============================================================",
        f"Config            : {cfg['_config_path']}",
        f"TE family filter  : {te_family or '<all>'}",
        f"Position window   : {window} bp",
        f"Min junction reads: {min_jr}",
        f"Samples compared  : {len(samples)}",
        "",
        f"RelocaTE2 calls (filtered) : {r2_total}",
        f"RelocaTE3 calls (filtered) : {r3_total}",
        f"Shared                     : {total_shared}",
        f"RelocaTE2-only             : {total_r2_only}",
        f"RelocaTE3-only             : {total_r3_only}",
        "",
        f"Recall  (shared / R2 total): {recall:.4f}",
        f"Precision (shared / R3 tot): {precision:.4f}",
        f"Jaccard (shared / union)   : {jaccard:.4f}",
    ]
    summary_txt = report_dir / "summary.txt"
    summary_txt.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    if plt is not None:
        for s, shared, r2o, r3o in per_sample_counts:
            draw_one_venn(
                plt,
                venn2,
                shared,
                r2o,
                r3o,
                f"{s}: non-reference insertions",
                venn_dir / f"venn_{s}.png",
            )
        draw_one_venn(
            plt,
            venn2,
            total_shared,
            total_r2_only,
            total_r3_only,
            f"All samples (n={len(samples)}): non-reference insertions",
            report_dir / "venn_total.png",
        )

    print(f"\nReport directory: {report_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
