"""Compare RelocaTE3 vs legacy RelocaTE2 *characterized* (genotyped) calls.

Pairs are matched on position alone (greedy nearest-neighbor within
``compare_char.position_window`` bp on the same chrom + TE family). TSD and
genotype-status agreement are then reported as separate metrics over the
matched pairs.

Inputs:
    Two normalized TSVs produced by normalize_relocate2_char.py and
    normalize_relocate3_char.py (same schema; see _common.CHAR_COLUMNS).

Outputs (under <report_dir>/characterized/):
    summary.tsv           per-sample shared / only / agreement metrics
    summary.txt           plain-text overall stats
    matched_calls.tsv     paired calls with TSD/status/spanner detail
    relocate2_only.tsv    R2 calls with no R3 match within the window
    relocate3_only.tsv    R3 calls with no R2 match within the window
    tsd_confusion.tsv     R2 TSD x R3 TSD contingency over matched pairs
    status_confusion.tsv  R2 status x R3 status contingency
    venn_total.png        shared vs unique by position
    venn/venn_<sample>.png
    tsd_agreement.png     per-sample TSD-match-rate bar chart
    status_agreement.png  per-sample stacked status concordance
    status_confusion.png  heatmap of status_confusion.tsv
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
    filter_te_family,
    load_pyplot,
    load_venn,
    match_sample,
    read_tsv,
    write_rows,
)
from _config import load_config  # noqa: E402


def _norm_status(status: str, collapse_excision: bool) -> str:
    """Optionally drop ``/excision_with_footprint`` suffix for collapsed scoring."""
    if collapse_excision and "/excision_with_footprint" in status:
        return status.replace("/excision_with_footprint", "")
    return status


def _write_pairs(
    path: Path,
    matched: list[tuple[dict, dict]],
    r2_only: list[dict],
    collapse_excision: bool,
) -> None:
    cols = [
        "sample",
        "chrom",
        "te_name",
        "r2_start",
        "r2_end",
        "r2_strand",
        "r3_start",
        "r3_end",
        "r3_strand",
        "distance_bp",
        "r2_tsd",
        "r3_tsd",
        "tsd_match",
        "r2_status",
        "r3_status",
        "status_match",
        "status_match_collapsed",
        "r2_avg_flankers",
        "r3_avg_flankers",
        "r2_spanners",
        "r3_spanners",
        "matched",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for rec2, rec3 in matched:
            r2_status = rec2["status"]
            r3_status = rec3["status"]
            tsd_match = int(rec2["tsd"] == rec3["tsd"])
            status_match = int(r2_status == r3_status)
            collapsed = int(
                _norm_status(r2_status, collapse_excision)
                == _norm_status(r3_status, collapse_excision)
            )
            w.writerow(
                {
                    "sample": rec2["sample"],
                    "chrom": rec2["chrom"],
                    "te_name": rec2["te_name"],
                    "r2_start": rec2["start"],
                    "r2_end": rec2["end"],
                    "r2_strand": rec2["strand"],
                    "r3_start": rec3["start"],
                    "r3_end": rec3["end"],
                    "r3_strand": rec3["strand"],
                    "distance_bp": abs(int(rec3["start"]) - int(rec2["start"])),
                    "r2_tsd": rec2["tsd"],
                    "r3_tsd": rec3["tsd"],
                    "tsd_match": tsd_match,
                    "r2_status": r2_status,
                    "r3_status": r3_status,
                    "status_match": status_match,
                    "status_match_collapsed": collapsed,
                    "r2_avg_flankers": rec2["avg_flankers"],
                    "r3_avg_flankers": rec3["avg_flankers"],
                    "r2_spanners": rec2["spanners"],
                    "r3_spanners": rec3["spanners"],
                    "matched": 1,
                }
            )
        for rec2 in r2_only:
            w.writerow(
                {
                    "sample": rec2["sample"],
                    "chrom": rec2["chrom"],
                    "te_name": rec2["te_name"],
                    "r2_start": rec2["start"],
                    "r2_end": rec2["end"],
                    "r2_strand": rec2["strand"],
                    "r3_start": "",
                    "r3_end": "",
                    "r3_strand": "",
                    "distance_bp": "",
                    "r2_tsd": rec2["tsd"],
                    "r3_tsd": "",
                    "tsd_match": "",
                    "r2_status": rec2["status"],
                    "r3_status": "",
                    "status_match": "",
                    "status_match_collapsed": "",
                    "r2_avg_flankers": rec2["avg_flankers"],
                    "r3_avg_flankers": "",
                    "r2_spanners": rec2["spanners"],
                    "r3_spanners": "",
                    "matched": 0,
                }
            )


def _write_confusion(path: Path, confusion: dict[tuple[str, str], int]) -> None:
    """Write a (R2, R3) -> count contingency as a long-form TSV."""
    rows = [
        {"r2": k[0], "r3": k[1], "count": v}
        for k, v in sorted(confusion.items(), key=lambda kv: (-kv[1], kv[0]))
    ]
    if not rows:
        path.write_text("r2\tr3\tcount\n")
        return
    write_rows(path, rows)


def _plot_tsd_agreement(plt, summary_rows: list[dict], out_png: Path) -> None:
    samples = [r["sample"] for r in summary_rows]
    rates = [
        (r["tsd_match_n"] / r["shared"]) if r["shared"] else 0.0
        for r in summary_rows
    ]
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(samples) + 2), 4))
    ax.bar(samples, rates, color="#4C72B0")
    ax.set_ylim(0, 1)
    ax.set_ylabel("TSD match rate (matched pairs)")
    ax.set_title("TSD agreement: RelocaTE3 vs RelocaTE2")
    ax.axhline(1.0, color="grey", linewidth=0.5, linestyle="--")
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        tick.set_ha("right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def _plot_status_agreement(plt, summary_rows: list[dict], out_png: Path) -> None:
    samples = [r["sample"] for r in summary_rows]
    matches = [r["status_match_n"] for r in summary_rows]
    mismatches = [max(r["shared"] - r["status_match_n"], 0) for r in summary_rows]
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(samples) + 2), 4))
    ax.bar(samples, matches, label="status match", color="#55A868")
    ax.bar(samples, mismatches, bottom=matches, label="disagreement", color="#C44E52")
    ax.set_ylabel("Matched-pair count")
    ax.set_title("Status concordance per sample")
    ax.legend()
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        tick.set_ha("right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def _plot_status_confusion(
    plt, confusion: dict[tuple[str, str], int], out_png: Path
) -> None:
    if not confusion:
        return
    r2_labels = sorted({k[0] for k in confusion})
    r3_labels = sorted({k[1] for k in confusion})
    grid = [
        [confusion.get((r2, r3), 0) for r3 in r3_labels] for r2 in r2_labels
    ]
    fig, ax = plt.subplots(
        figsize=(max(4, 0.7 * len(r3_labels) + 2), max(3, 0.5 * len(r2_labels) + 2))
    )
    im = ax.imshow(grid, cmap="Blues", aspect="auto")
    ax.set_xticks(range(len(r3_labels)))
    ax.set_xticklabels(r3_labels, rotation=45, ha="right")
    ax.set_yticks(range(len(r2_labels)))
    ax.set_yticklabels(r2_labels)
    ax.set_xlabel("RelocaTE3 status")
    ax.set_ylabel("RelocaTE2 status")
    ax.set_title("Status confusion (matched pairs)")
    grid_max = max((max(row) for row in grid), default=0)
    for i, r2 in enumerate(r2_labels):
        for j, r3 in enumerate(r3_labels):
            ax.text(
                j,
                i,
                grid[i][j],
                ha="center",
                va="center",
                color="black" if grid[i][j] < grid_max / 2 else "white",
                fontsize=8,
            )
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.toml", help="Validation TOML config")
    ap.add_argument("--r2-tsv", help="RelocaTE2 char TSV (default: <report>/characterized/relocate2_calls.tsv)")
    ap.add_argument("--r3-tsv", help="RelocaTE3 char TSV (default: <report>/characterized/relocate3_calls.tsv)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    report_dir = Path(cfg["paths"]["report_dir"]) / "characterized"
    report_dir.mkdir(parents=True, exist_ok=True)

    r2_tsv = Path(args.r2_tsv) if args.r2_tsv else report_dir / "relocate2_calls.tsv"
    r3_tsv = Path(args.r3_tsv) if args.r3_tsv else report_dir / "relocate3_calls.tsv"

    if not r2_tsv.is_file():
        print(f"ERROR: missing {r2_tsv} (run normalize_relocate2_char.py first)", file=sys.stderr)
        return 1
    if not r3_tsv.is_file():
        print(f"ERROR: missing {r3_tsv} (run normalize_relocate3_char.py first)", file=sys.stderr)
        return 1

    char_cfg = cfg.get("compare_char", {})
    window = int(char_cfg.get("position_window", cfg["compare"].get("position_window", 100)))
    te_family = char_cfg.get("te_family", cfg["compare"].get("te_family", ""))
    collapse_excision = bool(char_cfg.get("collapse_excision_suffix", False))

    r2_all = filter_te_family(read_tsv(r2_tsv), te_family)
    r3_all = filter_te_family(read_tsv(r3_tsv), te_family)

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
    tsd_confusion: dict[tuple[str, str], int] = defaultdict(int)
    status_confusion: dict[tuple[str, str], int] = defaultdict(int)
    summary_rows: list[dict] = []
    per_sample_counts: list[tuple[str, int, int, int]] = []
    total_shared = total_r2_only = total_r3_only = 0
    total_tsd_match = total_status_match = total_status_match_collapsed = 0

    plt_v, venn2 = load_venn()
    venn_dir = report_dir / "venn"
    if plt_v is not None:
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

        tsd_match_n = 0
        status_match_n = 0
        status_match_n_collapsed = 0
        for rec2, rec3 in matched:
            if rec2["tsd"] == rec3["tsd"]:
                tsd_match_n += 1
            if rec2["status"] == rec3["status"]:
                status_match_n += 1
            if _norm_status(rec2["status"], collapse_excision) == _norm_status(
                rec3["status"], collapse_excision
            ):
                status_match_n_collapsed += 1
            tsd_confusion[(rec2["tsd"], rec3["tsd"])] += 1
            status_confusion[(rec2["status"], rec3["status"])] += 1

        total_shared += shared
        total_r2_only += len(r2_only)
        total_r3_only += len(r3_only)
        total_tsd_match += tsd_match_n
        total_status_match += status_match_n
        total_status_match_collapsed += status_match_n_collapsed

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
                "tsd_match_n": tsd_match_n,
                "tsd_match_rate": f"{(tsd_match_n / shared) if shared else 0.0:.4f}",
                "status_match_n": status_match_n,
                "status_match_rate": f"{(status_match_n / shared) if shared else 0.0:.4f}",
                "status_match_n_collapsed": status_match_n_collapsed,
            }
        )

    summary_tsv = report_dir / "summary.tsv"
    with open(summary_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    _write_pairs(report_dir / "matched_calls.tsv", all_matched, all_r2_only, collapse_excision)
    write_rows(report_dir / "relocate2_only.tsv", all_r2_only)
    write_rows(report_dir / "relocate3_only.tsv", all_r3_only)
    _write_confusion(report_dir / "tsd_confusion.tsv", tsd_confusion)
    _write_confusion(report_dir / "status_confusion.tsv", status_confusion)

    r2_total = sum(len(by_sample_r2.get(s, [])) for s in samples)
    r3_total = sum(len(by_sample_r3.get(s, [])) for s in samples)
    union = total_shared + total_r2_only + total_r3_only
    recall = total_shared / r2_total if r2_total else 0.0
    precision = total_shared / r3_total if r3_total else 0.0
    jaccard = total_shared / union if union else 0.0
    tsd_rate = total_tsd_match / total_shared if total_shared else 0.0
    status_rate = total_status_match / total_shared if total_shared else 0.0
    status_rate_coll = (
        total_status_match_collapsed / total_shared if total_shared else 0.0
    )

    lines = [
        "RelocaTE3 vs RelocaTE2 comparison (characterized / genotyped calls)",
        "===================================================================",
        f"Config            : {cfg['_config_path']}",
        f"TE family filter  : {te_family or '<all>'}",
        f"Position window   : {window} bp",
        f"Collapse excision : {collapse_excision}",
        f"Samples compared  : {len(samples)}",
        "",
        f"RelocaTE2 calls (filtered) : {r2_total}",
        f"RelocaTE3 calls (filtered) : {r3_total}",
        f"Shared (position-matched)  : {total_shared}",
        f"RelocaTE2-only             : {total_r2_only}",
        f"RelocaTE3-only             : {total_r3_only}",
        "",
        f"Recall  (shared / R2 total): {recall:.4f}",
        f"Precision (shared / R3 tot): {precision:.4f}",
        f"Jaccard (shared / union)   : {jaccard:.4f}",
        "",
        f"TSD agreement     : {total_tsd_match}/{total_shared} ({tsd_rate:.4f})",
        f"Status agreement  : {total_status_match}/{total_shared} ({status_rate:.4f})",
        f"Status agreement (collapsed): {total_status_match_collapsed}/{total_shared} ({status_rate_coll:.4f})",
    ]
    summary_txt = report_dir / "summary.txt"
    summary_txt.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    if plt_v is not None:
        for s, shared, r2o, r3o in per_sample_counts:
            draw_one_venn(
                plt_v,
                venn2,
                shared,
                r2o,
                r3o,
                f"{s}: characterized insertions",
                venn_dir / f"venn_{s}.png",
            )
        draw_one_venn(
            plt_v,
            venn2,
            total_shared,
            total_r2_only,
            total_r3_only,
            f"All samples (n={len(samples)}): characterized insertions",
            report_dir / "venn_total.png",
        )

    plt = load_pyplot()
    if plt is not None and summary_rows:
        # Convert string rates back where helpful (already ints for counts).
        rows_for_plot = [
            {**r, "shared": int(r["shared"]),
             "tsd_match_n": int(r["tsd_match_n"]),
             "status_match_n": int(r["status_match_n"])}
            for r in summary_rows
        ]
        _plot_tsd_agreement(plt, rows_for_plot, report_dir / "tsd_agreement.png")
        _plot_status_agreement(plt, rows_for_plot, report_dir / "status_agreement.png")
        _plot_status_confusion(plt, status_confusion, report_dir / "status_confusion.png")

    print(f"\nReport directory: {report_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
