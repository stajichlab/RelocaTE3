"""Tests for ``relocaTE3 sample-stats``: per-sample counts for the MultiQC report.

The counts come from files a sample's run already wrote (run-all layout):

- ``te_containing/<id>.read_repeat_name.txt``: one line per TE-containing read;
- ``flanking/<id>.*.flankingReads.fq``: FASTQ, 4 lines per read;
- ``results/<target>.<te_name>.all_nonref_insert.{raw,all,high_conf}.txt``:
  no header, one row per call;
- ``results/<target>.<te_name>.all_nonref_insert.gff``: the headline tier
  (RelocaTE2's set; the plain ``.txt`` is the unfiltered table), no header;
- ``results/<target>.<te_name>.all_nonref_insert.characTErized.txt``: one
  header line (optional file: only when genotyped);
- ``results/<id>.all_ref_insert.txt``: one header line (optional file: only
  with a RepeatMasker annotation).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from RelocaTE3 import cli
from RelocaTE3.cli import main
from RelocaTE3.sample_stats import (
    COLUMNS,
    SampleStatsError,
    sample_stats,
)

REAL_SAMPLE = Path(
    "/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/"
    "RelocaTE3/results/nextflow-m2/nextflow-real-rice/dev/G_1"
)


def _fq(n: int) -> str:
    return "".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(n))


def _rows(n: int, prefix: str = "mPing\tTAA\tS\tChr1") -> str:
    return "".join(f"{prefix}\t{i}..{i + 2}\t+\n" for i in range(n))


def _sample(tmp_path: Path, name: str = "S", target: str = "ALL", te: str = "repeat",
            characterized: bool = True, reference: bool = True) -> Path:  # fmt: skip
    d = tmp_path / name
    for sub in ("te_containing", "flanking", "results"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    (d / "te_containing" / f"{name}.read_repeat_name.txt").write_text(
        "".join(f"r{i}/1\tmPing\t+\n" for i in range(7))
    )
    (d / "te_containing" / f"{name}.te_hit_names.txt").write_text("x\n" * 99)
    (d / "flanking" / f"{name}.left.flankingReads.fq").write_text(_fq(5))
    (d / "flanking" / f"{name}.right.flankingReads.fq").write_text(_fq(4))
    stem = d / "results" / f"{target}.{te}.all_nonref_insert"
    Path(f"{stem}.raw.txt").write_text(_rows(6))
    Path(f"{stem}.all.txt").write_text(_rows(5))
    Path(f"{stem}.gff").write_text(_rows(4))
    Path(f"{stem}.high_conf.txt").write_text(_rows(3))
    Path(f"{stem}.txt").write_text(_rows(6))
    if characterized:
        Path(f"{stem}.characTErized.txt").write_text("strain\tTE\n" + "S\tmPing\n" * 2)
    if reference:
        (d / "results" / f"{name}.all_ref_insert.txt").write_text(
            "chrom\tstart\n" + "Chr1\t1\n" * 8
        )
    return d


# ---------------------------------------------------------------- counting


def test_counts_each_file_and_drops_headers(tmp_path):
    d = _sample(tmp_path)
    stats = sample_stats(d, "S", target="ALL", te_name="repeat")
    assert stats == {
        "te_containing_reads": 7,
        "flanking_reads": 9,
        "nonref_raw": 6,
        "nonref_all": 5,
        "nonref_headline": 4,
        "nonref_high_conf": 3,
        "characterized": 2,
        "reference": 8,
    }
    assert list(stats) == list(COLUMNS)


def test_optional_files_are_left_out_not_zero(tmp_path):
    d = _sample(tmp_path, characterized=False, reference=False)
    stats = sample_stats(d, "S", target="ALL", te_name="repeat")
    assert "characterized" not in stats
    assert "reference" not in stats
    assert stats["nonref_raw"] == 6


def test_header_only_characterized_table_counts_zero(tmp_path):
    d = _sample(tmp_path)
    (d / "results" / "ALL.repeat.all_nonref_insert.characTErized.txt").write_text("strain\tTE\n")
    assert sample_stats(d, "S", target="ALL", te_name="repeat")["characterized"] == 0


def test_empty_tier_counts_zero_and_blank_lines_are_ignored(tmp_path):
    d = _sample(tmp_path)
    (d / "results" / "ALL.repeat.all_nonref_insert.high_conf.txt").write_text("")
    (d / "results" / "ALL.repeat.all_nonref_insert.raw.txt").write_text(_rows(2) + "\n")
    stats = sample_stats(d, "S", target="ALL", te_name="repeat")
    assert stats["nonref_high_conf"] == 0
    assert stats["nonref_raw"] == 2


def test_target_and_te_name_pick_the_tables(tmp_path):
    d = _sample(tmp_path, target="Chr3", te="mping")
    stats = sample_stats(d, "S", target="Chr3", te_name="mping")
    assert stats["nonref_raw"] == 6
    with pytest.raises(SampleStatsError, match="ALL.repeat.all_nonref_insert.raw.txt"):
        sample_stats(d, "S", target="ALL", te_name="repeat")


def test_single_end_has_one_flanking_file(tmp_path):
    d = _sample(tmp_path)
    (d / "flanking" / "S.right.flankingReads.fq").unlink()
    assert sample_stats(d, "S", target="ALL", te_name="repeat")["flanking_reads"] == 5


def test_missing_required_file_fails(tmp_path):
    d = _sample(tmp_path)
    (d / "te_containing" / "S.read_repeat_name.txt").unlink()
    with pytest.raises(SampleStatsError, match="read_repeat_name"):
        sample_stats(d, "S", target="ALL", te_name="repeat")


def test_no_flanking_fastq_fails(tmp_path):
    d = _sample(tmp_path)
    for f in (d / "flanking").iterdir():
        f.unlink()
    with pytest.raises(SampleStatsError, match="flankingReads"):
        sample_stats(d, "S", target="ALL", te_name="repeat")


def test_truncated_fastq_fails(tmp_path):
    d = _sample(tmp_path)
    (d / "flanking" / "S.left.flankingReads.fq").write_text(_fq(2) + "@r9\nACGT\n")
    with pytest.raises(SampleStatsError, match="multiple of 4"):
        sample_stats(d, "S", target="ALL", te_name="repeat")


# ---------------------------------------------------------------- CLI + file


def test_cli_writes_a_multiqc_custom_content_file(tmp_path):
    d = _sample(tmp_path)
    out = tmp_path / "S.stats_mqc.tsv"
    assert main(["sample-stats", "-d", str(d), "-n", "S", "-o", str(out)]) == 0
    lines = out.read_text().splitlines()
    comments = [line for line in lines if line.startswith("#")]
    assert lines[: len(comments)] == comments
    text = "\n".join(c[2:] for c in comments)
    assert "id: 'relocate3_sample_stats'" in text
    assert "plot_type: 'table'" in text
    assert "format: 'tsv'" in text
    data = [line.split("\t") for line in lines[len(comments):]]
    assert data == [
        ["Sample", *COLUMNS],
        ["S", "7", "9", "6", "5", "4", "3", "2", "8"],
    ]


def test_cli_header_yaml_parses_and_describes_every_column(tmp_path):
    yaml = pytest.importorskip("yaml")
    d = _sample(tmp_path)
    out = tmp_path / "S.stats_mqc.tsv"
    assert main(["sample-stats", "-d", str(d), "-n", "S", "-o", str(out)]) == 0
    comments = [l[2:] for l in out.read_text().splitlines() if l.startswith("# ")]
    conf = yaml.safe_load("\n".join(comments))
    assert conf["id"] == "relocate3_sample_stats"
    assert conf["pconfig"]["id"] == "relocate3_sample_stats_table"
    assert set(conf["headers"]) == set(COLUMNS)
    assert all(h["title"] and h["description"] for h in conf["headers"].values())


def test_cli_passes_target_and_te_name(tmp_path, monkeypatch):
    seen = {}
    monkeypatch.setattr(cli, "cmd_sample_stats", lambda a: (seen.setdefault("a", a), 0)[1])
    rc = main(["sample-stats", "-d", "x", "-n", "S", "--target", "Chr3", "--te-name", "mping",
               "-o", "o.tsv"])  # fmt: skip
    assert rc == 0
    a = seen["a"]
    assert (a.sample_dir, a.name, a.target, a.te_name, a.output) == ("x", "S", "Chr3", "mping", "o.tsv")


def test_cli_defaults_match_run_all(tmp_path):
    a = cli.build_parser().parse_args(["sample-stats", "-d", "x", "-n", "S", "-o", "o"])
    r = cli.build_parser().parse_args(["run-all", "-1", "r1.fq", "-T", "te.fa", "-g", "g.fa", "-n", "s"])
    assert (a.target, a.te_name) == (r.target, r.te_name)


# ---------------------------------------------------------------- real data


def _wc_l(path: Path) -> int:
    with open(path, "rb") as fh:
        return sum(1 for _ in fh)


@pytest.mark.skipif(not REAL_SAMPLE.is_dir(), reason="real-rice M2 output G_1 not present")
def test_real_g1_counts_match_wc_l():
    d = REAL_SAMPLE
    stats = sample_stats(d, "G_1", target="ALL", te_name="mping")
    stem = d / "results" / "ALL.mping.all_nonref_insert"
    flank = sorted((d / "flanking").glob("G_1.*.flankingReads.fq"))
    assert len(flank) == 2
    assert stats == {
        "te_containing_reads": _wc_l(d / "te_containing" / "G_1.read_repeat_name.txt"),
        "flanking_reads": sum(_wc_l(f) for f in flank) // 4,
        "nonref_raw": _wc_l(Path(f"{stem}.raw.txt")),
        "nonref_all": _wc_l(Path(f"{stem}.all.txt")),
        "nonref_headline": _wc_l(Path(f"{stem}.gff")),
        "nonref_high_conf": _wc_l(Path(f"{stem}.high_conf.txt")),
        "characterized": _wc_l(Path(f"{stem}.characTErized.txt")) - 1,
        "reference": _wc_l(d / "results" / "G_1.all_ref_insert.txt") - 1,
    }
    # Values observed when this test was written (M2 real-rice run).
    assert stats["te_containing_reads"] == 60604
    assert (stats["nonref_raw"], stats["nonref_high_conf"], stats["characterized"]) == (623, 621, 623)
    assert stats["reference"] == 13


def test_fastq_trailing_blank_lines_are_tolerated(tmp_path):
    d = _sample(tmp_path)
    (d / "flanking" / "S.left.flankingReads.fq").write_text(_fq(5) + "\n\n")
    assert sample_stats(d, "S", target="ALL", te_name="repeat")["flanking_reads"] == 9


def test_flanking_reads_description_says_given_to_the_aligner():
    from RelocaTE3.sample_stats import _HEADERS

    desc = _HEADERS["flanking_reads"][1]
    assert "given to the genome aligner" in desc
    assert "aligned to" not in desc
