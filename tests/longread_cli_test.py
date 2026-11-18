"""The lr-align, lr-signatures and lr-cluster subcommands.

Parser tests need no tools; the end-to-end tests align the synthetic genome
and reads of ``tests/longread_align_test.py`` (a 1 kb insert at 10 kb) with
the real minimap2 and samtools and are skipped without them.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pysam
import pytest

from RelocaTE3.cli import main
from RelocaTE3.longread.model import ClipSide, Locus, ReadType, Signature, SignatureKind
from RelocaTE3.longread.tables import read_clusters, read_signatures, write_signatures
from tests import longread_align_test
from tests.longread_align_test import HAVE_TOOLS, INS_AT

# The synthetic genome (1 kb insert at INS_AT) and three reads over it.
lr_synthetic = longread_align_test.synthetic

needs_tools = pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH")


def _timing(outdir: Path) -> list[dict[str, str]]:
    with (outdir / "timing.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- parser -----------------------------------------------------------------


@pytest.mark.parametrize("command", ["lr-align", "lr-signatures", "lr-cluster"])
def test_help_exits_zero(command: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert main([command, "--help"]) == 0
    out = capsys.readouterr().out
    assert f"usage: relocaTE3 {command}" in out
    assert "--sample" in out and "--outdir" in out


def test_top_level_help_lists_long_read_commands(
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert main(["--help"]) == 0
    out = capsys.readouterr().out
    for command in ("lr-align", "lr-signatures", "lr-cluster", "run-all"):
        assert command in out


def test_lr_align_rejects_reads_and_bam_together(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rc = main(["lr-align", "--sample", "S", "-o", str(tmp_path / "o"),
               "--genome", "g.fa", "--read-type", "hifi",
               "--reads", "r.fq", "--bam", "a.bam"])  # fmt: skip
    assert rc != 0
    assert "not allowed with argument" in capsys.readouterr().err
    assert not (tmp_path / "o").exists()


def test_lr_align_requires_reads_or_bam(tmp_path: Path) -> None:
    rc = main(["lr-align", "--sample", "S", "-o", str(tmp_path),
               "--genome", "g.fa", "--read-type", "hifi"])  # fmt: skip
    assert rc != 0


@pytest.mark.parametrize("command", ["lr-align", "lr-signatures"])
def test_bad_read_type_rejected(
    command: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    extra = (["--genome", "g.fa", "--reads", "r.fq"] if command == "lr-align"
             else ["--bam", "a.bam"])  # fmt: skip
    rc = main([command, "--sample", "S", "-o", str(tmp_path),
               "--read-type", "illumina", *extra])  # fmt: skip
    assert rc != 0
    assert "invalid choice: 'illumina'" in capsys.readouterr().err


def test_lr_signatures_missing_bam_is_an_error(tmp_path: Path) -> None:
    rc = main(["lr-signatures", "--sample", "S", "-o", str(tmp_path),
               "--bam", str(tmp_path / "nope.bam"), "--read-type", "hifi"])  # fmt: skip
    assert rc == 1
    assert not (tmp_path / "S.signatures.tsv.gz").exists()


def test_lr_cluster_missing_signatures_is_an_error(tmp_path: Path) -> None:
    rc = main(["lr-cluster", "--sample", "S", "-o", str(tmp_path),
               "--signatures", str(tmp_path / "nope.tsv.gz")])  # fmt: skip
    assert rc == 1


# --- lr-cluster window default ------------------------------------------------


def _ins(sig_id: str, read_id: str, pos: int) -> Signature:
    return Signature(
        signature_id=sig_id, kind=SignatureKind.INS,
        locus=Locus("ref", "chr1", pos, pos, "+"), insert_length=500,
        clip_side=ClipSide.NONE, read_id=read_id, query_start=1000,
        query_end=1500, mapq=60, sequence="A" * 500, left_flank="C",
        right_flank="G",
    )  # fmt: skip


def _two_signatures_60bp_apart(path: Path, read_type: ReadType) -> Path:
    sigs = [_ins("s1", "r1", 5000), _ins("s2", "r2", 5060)]
    write_signatures(path, sigs, read_type=read_type, reference_id="ref")
    return path


def _n_clusters(outdir: Path, sample: str) -> int:
    _, clusters = read_clusters(outdir / f"{sample}.clusters.tsv.gz")
    return len(list(clusters))


def test_lr_cluster_window_defaults_from_read_type(tmp_path: Path) -> None:
    sigs = _two_signatures_60bp_apart(tmp_path / "hifi.sigs.tsv.gz", ReadType.HIFI)
    assert main(["lr-cluster", "--sample", "A", "-o", str(tmp_path / "a"),
                 "--signatures", str(sigs)]) == 0  # fmt: skip
    assert _n_clusters(tmp_path / "a", "A") == 2  # hifi window 50 < 60

    assert main(["lr-cluster", "--sample", "B", "-o", str(tmp_path / "b"),
                 "--signatures", str(sigs), "--cluster-window", "100"]) == 0  # fmt: skip
    assert _n_clusters(tmp_path / "b", "B") == 1

    ont = _two_signatures_60bp_apart(tmp_path / "ont.sigs.tsv.gz", ReadType.ONT)
    assert main(["lr-cluster", "--sample", "C", "-o", str(tmp_path / "c"),
                 "--signatures", str(ont)]) == 0  # fmt: skip
    assert _n_clusters(tmp_path / "c", "C") == 1  # ont window 100 > 60
    meta, _ = read_clusters(tmp_path / "c" / "C.clusters.tsv.gz")
    assert meta["read_type"] == "ont" and meta["reference_id"] == "ref"


def test_lr_cluster_refuses_existing_output(tmp_path: Path) -> None:
    sigs = _two_signatures_60bp_apart(tmp_path / "s.tsv.gz", ReadType.HIFI)
    args = ["lr-cluster", "--sample", "A", "-o", str(tmp_path / "a"),
            "--signatures", str(sigs)]  # fmt: skip
    assert main(args) == 0
    assert main(args) == 1
    rows = _timing(tmp_path / "a")
    assert [r["status"] for r in rows] == ["ok", "FileExistsError"]


# --- end to end ---------------------------------------------------------------


@needs_tools
def test_end_to_end_align_signatures_cluster(
    lr_synthetic: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    fa, fq = lr_synthetic
    out = tmp_path / "out"
    assert main(["lr-align", "--sample", "S", "-o", str(out), "--genome", str(fa),
                 "--read-type", "hifi", "--reads", str(fq),
                 "--threads", "2"]) == 0  # fmt: skip
    bam = out / "S.lr.bam"
    assert bam.is_file() and Path(str(bam) + ".bai").is_file()

    assert main(["lr-signatures", "--sample", "S", "-o", str(out), "--bam", str(bam),
                 "--read-type", "hifi"]) == 0  # fmt: skip
    sig_path = out / "S.signatures.tsv.gz"
    meta, sigs = read_signatures(sig_path)
    n_sigs = len(list(sigs))
    assert meta["read_type"] == "hifi" and meta["reference_id"] == "ref"

    assert main(["lr-cluster", "--sample", "S", "-o", str(out),
                 "--signatures", str(sig_path)]) == 0  # fmt: skip
    meta, clusters = read_clusters(out / "S.clusters.tsv.gz")
    clusters = list(clusters)
    assert meta["read_type"] == "hifi" and meta["reference_id"] == "ref"
    near = [c for c in clusters if abs(c.locus.start - INS_AT) <= 20]
    assert len(near) == 1, clusters
    assert near[0].n_reads == 3

    rows = _timing(out)
    assert [r["stage"] for r in rows] == ["lr-align", "lr-signatures", "lr-cluster"]
    assert all(r["status"] == "ok" and r["sample"] == "S" for r in rows)
    with pysam.AlignmentFile(str(bam)) as fh:
        mapped = fh.mapped
    assert int(rows[0]["records_out"]) == mapped > 0
    assert int(rows[1]["records_out"]) == n_sigs
    assert int(rows[1]["records_in"]) > 0
    assert int(rows[2]["records_in"]) == n_sigs
    assert int(rows[2]["records_out"]) == len(clusters)

    # lr-align --bam: links the BAM and its index after checking them.
    linked = tmp_path / "linked"
    args = ["lr-align", "--sample", "L", "-o", str(linked), "--genome", str(fa),
            "--read-type", "hifi", "--bam", str(bam)]  # fmt: skip
    assert main(args) == 0
    link = linked / "L.lr.bam"
    link_bai = linked / "L.lr.bam.bai"
    assert link.is_symlink() and link_bai.is_symlink()
    assert Path(link.readlink()).is_absolute()
    assert link.resolve() == bam.resolve()
    assert link_bai.resolve() == Path(str(bam) + ".bai").resolve()
    assert _timing(linked)[0]["records_out"] == str(mapped)

    assert main(args) == 1  # the links exist now
    assert [r["status"] for r in _timing(linked)] == ["ok", "FileExistsError"]


@needs_tools
def test_lr_align_bam_rejects_wrong_reference(
    lr_synthetic: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    fa, fq = lr_synthetic
    out = tmp_path / "out"
    assert main(["lr-align", "--sample", "S", "-o", str(out), "--genome", str(fa),
                 "--read-type", "hifi", "--reads", str(fq)]) == 0  # fmt: skip
    other = tmp_path / "other.fa"
    other.write_text(">chrX\n" + "A" * 1000 + "\n")
    rc = main(["lr-align", "--sample", "L", "-o", str(tmp_path / "l"),
               "--genome", str(other), "--read-type", "hifi",
               "--bam", str(out / "S.lr.bam")])  # fmt: skip
    assert rc == 1
    assert not (tmp_path / "l" / "L.lr.bam").exists()
