"""Tests for the small-panel checker."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation" / "longread_small"))

import check_panel as cp  # noqa: E402


def _fastq(path: Path, lengths: list[int]) -> None:
    with gzip.open(path, "wt") as handle:
        for index, length in enumerate(lengths):
            handle.write(f"@r{index}\n{'A' * length}\n+\n{'I' * length}\n")


def test_fastq_stats(tmp_path):
    path = tmp_path / "r.fastq.gz"
    _fastq(path, [100, 300, 600])
    stats = cp.fastq_stats(path)
    assert stats == {"reads": 3, "bases": 1000, "n50": 600, "max_length": 600}


def test_fastq_stats_accepts_header_comments(tmp_path):
    path = tmp_path / "r.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@m1/1/ccs np=10 rq=0.99\nACGT\n+\nIIII\n")
    assert cp.fastq_stats(path)["bases"] == 4


def test_fastq_stats_rejects_truncated(tmp_path):
    path = tmp_path / "r.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@r0\nACGT\n+\n")
    with pytest.raises(ValueError, match="truncated"):
        cp.fastq_stats(path)


def test_fastq_stats_rejects_quality_length_mismatch(tmp_path):
    path = tmp_path / "r.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@r0\nACGT\n+\nIII\n")
    with pytest.raises(ValueError, match="length") as caught:
        cp.fastq_stats(path)
    assert str(path) in str(caught.value)


def test_fastq_stats_rejects_empty(tmp_path):
    path = tmp_path / "r.fastq.gz"
    _fastq(path, [])
    with pytest.raises(ValueError, match="no reads"):
        cp.fastq_stats(path)


def test_fastq_stats_rejects_cut_gzip(tmp_path):
    path = tmp_path / "r.fastq.gz"
    _fastq(path, [500] * 50)
    data = path.read_bytes()
    path.write_bytes(data[: len(data) // 2])
    with pytest.raises((EOFError, OSError, ValueError)):
        cp.fastq_stats(path)


def test_count_clipped_reads():
    sam = [
        "@SQ\tSN:Chr3\tLN:10",
        "r1\t0\tChr3\t1\t60\t100M\t*\t0\t0\t*\t*",
        "r2\t0\tChr3\t1\t60\t30S70M\t*\t0\t0\t*\t*",
        "r3\t0\tChr3\t1\t60\t80M20S\t*\t0\t0\t*\t*",
        "r4\t0\tChr3\t1\t60\t90M10S\t*\t0\t0\t*\t*",
        "r5\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*",
    ]
    assert cp.count_clipped_reads(sam, min_clip=20) == 2


def test_count_clipped_reads_skips_secondary_and_supplementary():
    sam = [
        "r1\t256\tChr3\t1\t0\t30S70M\t*\t0\t0\t*\t*",
        "r2\t2048\tChr3\t1\t60\t30S70M\t*\t0\t0\t*\t*",
    ]
    assert cp.count_clipped_reads(sam, min_clip=20) == 0


def test_coverage_within_tolerance():
    assert (
        cp.coverage_ok(bases=38_000_000, genome_length=2_000_000, requested=20) is True
    )
    assert (
        cp.coverage_ok(bases=20_000_000, genome_length=2_000_000, requested=20) is False
    )


def test_clipped_short_reads_reports_minimap2_stderr(tmp_path):
    reads = tmp_path / "r.fq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    missing = tmp_path / "no_such_genome.fa"
    with pytest.raises(RuntimeError) as caught:
        cp.clipped_short_reads(missing, reads, reads, threads=1)
    message = str(caught.value)
    assert str(missing) in message
    assert "minimap2" in message
    assert "stderr" in message


# --- main on a tiny synthetic panel (minimap2 mocked out) -------------------

GENOME_LENGTH = 1000


def _panel(tmp_path: Path, read_bases: int = 5000) -> dict[str, Path]:
    """Build a one-row panel: 5x requested on a 1000 bp reference."""
    panel = tmp_path / "panel"
    reference = tmp_path / "ref.fa"
    reference.write_text(">Chr3\n" + "A" * GENOME_LENGTH + "\n")
    genomes = panel / "truth" / "genomes"
    genomes.mkdir(parents=True)
    (genomes / "sample_h1.fa").write_text(">Sample#1#Chr3\n" + "A" * 1200 + "\n")
    (panel / "truth_events.tsv").write_text(
        "event_id\n" + "".join(f"e{i}\n" for i in range(200))
    )
    rows = ["sample\tplatform\tcoverage\treplicate\treads\tcontrol_reads"]
    reads = "reads/ont-hq/cov5x_rep1/reads.fastq.gz"
    control = "reads/ont-hq/cov5x_rep1_reference_control/reads.fastq.gz"
    rows.append(f"cov5x_rep1\tont-hq\t5\t1\t{reads}\t{control}")
    for relative in (reads, control):
        path = panel / relative
        path.parent.mkdir(parents=True)
        _fastq(path, [read_bases // 10] * 10)
        (path.parent / ".complete").touch()
    (panel / "panel_manifest.tsv").write_text("\n".join(rows) + "\n")
    short = tmp_path / "short.fq.gz"
    _fastq(short, [100])
    return {"panel": panel, "reference": reference, "short": short}


def _argv(paths: dict[str, Path]) -> list[str]:
    return [
        "--panel",
        str(paths["panel"]),
        "--reference",
        str(paths["reference"]),
        "--short-r1",
        str(paths["short"]),
        "--short-r2",
        str(paths["short"]),
    ]


@pytest.fixture
def fake_minimap2(monkeypatch):
    """Reference leaves 1000 clipped reads, the sample genome 10."""

    def fake(genome, r1, r2, threads):
        return 10 if genome.name == "sample_h1.fa" else 1000

    monkeypatch.setattr(cp, "clipped_short_reads", fake)


def test_main_pass_writes_report_and_marker(tmp_path, fake_minimap2, capsys):
    paths = _panel(tmp_path)
    assert cp.main(_argv(paths)) == 0
    panel = paths["panel"]
    assert (panel / ".complete").is_file()
    lines = (panel / "qc" / "panel_check.tsv").read_text().splitlines()
    assert lines[0].split("\t")[0] == "platform"
    assert len([line for line in lines[1:] if not line.startswith("#")]) == 2
    assert all(line.endswith("\tok") for line in lines[1:3])
    assert "# clipped_short_reads_reference\t1000" in lines
    assert "# clipped_short_reads_sample\t10" in lines
    assert "RESULT: PASS" in capsys.readouterr().out


def test_main_coverage_off_fails_without_marker(tmp_path, fake_minimap2, capsys):
    paths = _panel(tmp_path, read_bases=2000)  # 2x delivered, 5x requested
    assert cp.main(_argv(paths)) == 1
    assert not (paths["panel"] / ".complete").exists()
    assert "coverage off target" in capsys.readouterr().err
    lines = (paths["panel"] / "qc" / "panel_check.tsv").read_text().splitlines()
    assert lines[1].endswith("\tFAIL")


def test_main_clip_ratio_fails(tmp_path, monkeypatch, capsys):
    paths = _panel(tmp_path)
    monkeypatch.setattr(cp, "clipped_short_reads", lambda genome, *a: 100)
    assert cp.main(_argv(paths)) == 1
    assert "does not explain" in capsys.readouterr().err
    assert not (paths["panel"] / ".complete").exists()


def test_main_refuses_already_checked_panel(tmp_path, fake_minimap2, capsys):
    paths = _panel(tmp_path)
    (paths["panel"] / ".complete").touch()
    assert cp.main(_argv(paths)) == 1
    err = capsys.readouterr().err
    assert "already checked" in err
    assert str(paths["panel"] / ".complete") in err
    assert not (paths["panel"] / "qc").exists()


@pytest.mark.parametrize(
    "remove, expected",
    [
        ("panel_manifest.tsv", "panel_manifest.tsv"),
        ("truth_events.tsv", "truth_events.tsv"),
        ("truth/genomes/sample_h1.fa", "sample_h1.fa"),
        (
            "reads/ont-hq/cov5x_rep1/.complete",
            "no completion marker",
        ),
        (
            "reads/ont-hq/cov5x_rep1/reads.fastq.gz",
            "reads.fastq.gz",
        ),
    ],
)
def test_main_missing_inputs_fail_cleanly(
    tmp_path, fake_minimap2, capsys, remove, expected
):
    paths = _panel(tmp_path)
    (paths["panel"] / remove).unlink()
    assert cp.main(_argv(paths)) == 1
    err = capsys.readouterr().err
    assert "FAIL:" in err
    assert expected in err
    assert not (paths["panel"] / ".complete").exists()


def test_main_truth_wrong_count_fails(tmp_path, fake_minimap2, capsys):
    paths = _panel(tmp_path)
    (paths["panel"] / "truth_events.tsv").write_text("event_id\ne0\n")
    assert cp.main(_argv(paths)) == 1
    assert "expected 200" in capsys.readouterr().err


def test_main_empty_manifest_fails(tmp_path, fake_minimap2, capsys):
    paths = _panel(tmp_path)
    manifest = paths["panel"] / "panel_manifest.tsv"
    manifest.write_text(manifest.read_text().splitlines()[0] + "\n")
    assert cp.main(_argv(paths)) == 1
    assert "no rows" in capsys.readouterr().err


def test_main_empty_fastq_fails(tmp_path, fake_minimap2, capsys):
    paths = _panel(tmp_path)
    _fastq(paths["panel"] / "reads/ont-hq/cov5x_rep1/reads.fastq.gz", [])
    assert cp.main(_argv(paths)) == 1
    assert "no reads" in capsys.readouterr().err
