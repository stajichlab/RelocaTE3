"""The lr-align, lr-signatures and lr-cluster subcommands.

Parser tests need no tools; the end-to-end tests align the synthetic genome
and reads of ``tests/longread_align_test.py`` (a 1 kb insert at 10 kb) with
the real minimap2 and samtools and are skipped without them.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import pysam
import pytest

from RelocaTE3.cli import main
from RelocaTE3.longread.model import (
    NO_HIT, CallTier, ClipSide, Cluster, FilterReason, Locus, ReadType, Signature,
    SignatureKind, TeHit,
)  # fmt: skip
from RelocaTE3.longread.tables import (
    read_calls, read_clusters, read_signatures, read_te_hits, write_clusters,
    write_signatures, write_te_hits,
)  # fmt: skip
from tests import longread_align_test
from tests.longread_align_test import HAVE_TOOLS, INS_AT, _bam

# The synthetic genome (1 kb insert at INS_AT) and three reads over it.
lr_synthetic = longread_align_test.synthetic

needs_tools = pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH")


def _timing(outdir: Path) -> list[dict[str, str]]:
    with (outdir / "timing.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_meta_file(path: Path) -> dict[str, str]:
    with path.open() as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert rows[0] == ["key", "value"]
    return dict(rows[1:])


# --- parser -----------------------------------------------------------------


LR_COMMANDS = ["lr-align", "lr-signatures", "lr-cluster", "lr-annotate", "lr-call"]


@pytest.mark.parametrize("command", LR_COMMANDS)
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
    for command in (*LR_COMMANDS, "run-all"):
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


# --- dangling links, missing files, CRAM without --genome -----------------------


def test_lr_align_reads_refuses_dangling_output_links(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    fq = tmp_path / "r.fq"
    fq.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "out"
    out.mkdir()
    link = out / "S.lr.bam"
    link.symlink_to(tmp_path / "gone.bam")
    rc = main(["lr-align", "--sample", "S", "-o", str(out), "--genome", str(fa),
               "--read-type", "hifi", "--reads", str(fq)])  # fmt: skip
    assert rc == 1
    assert link.is_symlink() and link.readlink() == tmp_path / "gone.bam"
    assert _timing(out)[0]["status"] == "FileExistsError"

    link.unlink()
    index_link = out / "S.lr.bam.bai"
    index_link.symlink_to(tmp_path / "gone.bam.bai")
    rc = main(["lr-align", "--sample", "S", "-o", str(out), "--genome", str(fa),
               "--read-type", "hifi", "--reads", str(fq)])  # fmt: skip
    assert rc == 1
    assert index_link.is_symlink() and not link.exists()


def test_missing_input_dangling_symlink_is_named(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    sigs = tmp_path / "s.tsv.gz"
    sigs.symlink_to(tmp_path / "gone.tsv.gz")
    rc = main(["lr-cluster", "--sample", "S", "-o", str(tmp_path / "o"),
               "--signatures", str(sigs)])  # fmt: skip
    assert rc == 1
    assert (
        f"dangling symlink (target missing): {sigs} -> {tmp_path / 'gone.tsv.gz'}"
        in caplog.text
    )


def test_cram_without_genome_fails_early(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    cram = tmp_path / "a.cram"
    cram.write_bytes(b"not read")
    rc = main(["lr-signatures", "--sample", "S", "-o", str(tmp_path / "o"),
               "--bam", str(cram), "--read-type", "hifi"])  # fmt: skip
    assert rc == 1
    assert "CRAM input needs --genome (the reference FASTA)" in caplog.text
    assert not (tmp_path / "o" / "S.signatures.tsv.gz").exists()


# --- reference-id sidecar between lr-align and lr-signatures --------------------


def _user_bam(tmp_path: Path, name: str = "x.lr.bam") -> Path:
    return _bam(tmp_path / name, [("chr1", 1000)])


def _write_sidecar(bam: Path, reference_id: str, read_type: str) -> Path:
    meta = bam.with_suffix(".meta.tsv")
    meta.write_text(f"key\tvalue\nreference_id\t{reference_id}\n"
                    f"read_type\t{read_type}\ngenome\t/g.fa\n")  # fmt: skip
    return meta


def test_lr_signatures_rejects_sidecar_reference_mismatch(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    bam = _user_bam(tmp_path)
    meta = _write_sidecar(bam, "other", "hifi")
    assert meta.name == "x.lr.meta.tsv"
    rc = main(["lr-signatures", "--sample", "S", "-o", str(tmp_path / "o"),
               "--bam", str(bam), "--read-type", "hifi",
               "--reference-id", "ref"])  # fmt: skip
    assert rc == 1
    assert "'other'" in caplog.text and "'ref'" in caplog.text
    assert str(meta) in caplog.text
    assert not (tmp_path / "o" / "S.signatures.tsv.gz").exists()


def test_lr_signatures_rejects_sidecar_read_type_mismatch(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    bam = _user_bam(tmp_path)
    _write_sidecar(bam, "ref", "ont")
    rc = main(["lr-signatures", "--sample", "S", "-o", str(tmp_path / "o"),
               "--bam", str(bam), "--read-type", "hifi"])  # fmt: skip
    assert rc == 1
    assert "'ont'" in caplog.text and "'hifi'" in caplog.text


def test_lr_signatures_matching_sidecar_runs(tmp_path: Path) -> None:
    bam = _user_bam(tmp_path)
    _write_sidecar(bam, "ref", "hifi")
    assert main(["lr-signatures", "--sample", "S", "-o", str(tmp_path / "o"),
                 "--bam", str(bam), "--read-type", "hifi"]) == 0  # fmt: skip


def test_lr_signatures_without_sidecar_runs(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level("INFO")
    bam = _user_bam(tmp_path)
    assert main(["lr-signatures", "--sample", "S", "-o", str(tmp_path / "o"),
                 "--bam", str(bam), "--read-type", "hifi"]) == 0  # fmt: skip
    assert (tmp_path / "o" / "S.signatures.tsv.gz").is_file()
    assert "no sidecar" in caplog.text


@needs_tools
def test_lr_align_writes_sidecar_and_refuses_second_run(
    lr_synthetic: tuple[Path, Path], tmp_path: Path
) -> None:
    fa, fq = lr_synthetic
    out = tmp_path / "out"
    args = ["lr-align", "--sample", "S", "-o", str(out), "--genome", str(fa),
            "--read-type", "ont", "--reads", str(fq),
            "--reference-id", "IRGSP"]  # fmt: skip
    assert main(args) == 0
    meta_path = out / "S.lr.meta.tsv"
    assert _read_meta_file(meta_path) == {
        "reference_id": "IRGSP", "read_type": "ont", "genome": str(fa.resolve()),
    }  # fmt: skip
    assert not (out / "S.lr.meta.tsv.partial").exists()
    before = meta_path.read_bytes()
    stat = meta_path.stat()

    # Only the sidecar is left: the BAM is moved away, the sidecar still blocks.
    (out / "S.lr.bam").rename(tmp_path / "moved.bam")
    (out / "S.lr.bam.bai").rename(tmp_path / "moved.bam.bai")
    assert main(args) == 1
    assert not (out / "S.lr.bam").exists()
    assert meta_path.read_bytes() == before
    assert meta_path.stat().st_mtime_ns == stat.st_mtime_ns

    # --bam branch writes the sidecar too, and checks it before linking.
    linked = tmp_path / "linked"
    largs = ["lr-align", "--sample", "L", "-o", str(linked), "--genome", str(fa),
             "--read-type", "ont", "--bam", str(tmp_path / "moved.bam"),
             "--reference-id", "IRGSP"]  # fmt: skip
    assert main(largs) == 0
    assert _read_meta_file(linked / "L.lr.meta.tsv")["reference_id"] == "IRGSP"
    # lr-signatures reads the sidecar beside the linked BAM.
    sig_args = ["lr-signatures", "--sample", "L", "-o", str(linked),
                "--bam", str(linked / "L.lr.bam"), "--read-type", "ont"]  # fmt: skip
    assert main([*sig_args, "--reference-id", "ref"]) == 1
    assert main([*sig_args, "--reference-id", "IRGSP"]) == 0

    (linked / "L.lr.bam").unlink()
    (linked / "L.lr.bam.bai").unlink()
    assert main(largs) == 1
    assert not os.path.lexists(linked / "L.lr.bam")


# --- lr-annotate and lr-call ----------------------------------------------------


def test_lr_annotate_missing_inputs_are_errors(tmp_path: Path) -> None:
    lib = tmp_path / "lib.fa"
    lib.write_text(">synTE\nACGT\n")
    sigs = _two_signatures_60bp_apart(tmp_path / "s.tsv.gz", ReadType.HIFI)
    base = ["lr-annotate", "--sample", "S", "-o", str(tmp_path / "o")]
    assert main([*base, "--signatures", str(tmp_path / "nope.tsv.gz"),
                 "--te-library", str(lib)]) == 1  # fmt: skip
    assert main([*base, "--signatures", str(sigs),
                 "--te-library", str(tmp_path / "nope.fa")]) == 1  # fmt: skip
    assert not (tmp_path / "o" / "S.te_hits.tsv.gz").exists()
    assert [r["stage"] for r in _timing(tmp_path / "o")] == ["lr-annotate"] * 2


def test_lr_annotate_refuses_existing_output(tmp_path: Path) -> None:
    lib = tmp_path / "lib.fa"
    lib.write_text(">synTE\nACGT\n")
    sigs = _two_signatures_60bp_apart(tmp_path / "s.tsv.gz", ReadType.HIFI)
    out = tmp_path / "o"
    out.mkdir()
    (out / "S.te_hits.tsv.gz").write_bytes(b"old")
    assert main(["lr-annotate", "--sample", "S", "-o", str(out),
                 "--signatures", str(sigs), "--te-library", str(lib)]) == 1  # fmt: skip
    assert (out / "S.te_hits.tsv.gz").read_bytes() == b"old"
    assert _timing(out)[0]["status"] == "FileExistsError"


def _clip(sig_id: str, read_id: str, pos: int) -> Signature:
    return Signature(
        signature_id=sig_id, kind=SignatureKind.CLIP,
        locus=Locus("ref", "chr1", pos, pos, "+"), insert_length=300,
        clip_side=ClipSide.RIGHT, read_id=read_id, query_start=0,
        query_end=300, mapq=60, sequence="T" * 300, left_flank="C",
        right_flank="",
    )  # fmt: skip


def _hit(sig: Signature, family: str = "synTE", identity: float = 0.95) -> TeHit:
    n = sig.insert_length
    if family == NO_HIT:
        return TeHit(sig.signature_id, sig.read_id, NO_HIT, ".", 0, 0, 0, 0, 0, n, 0.0)
    return TeHit(sig.signature_id, sig.read_id, family, "+", 0, n, n, 0, n, n,
                 identity)  # fmt: skip


def _cluster(cid: str, pos: int, sigs: list[Signature]) -> Cluster:
    kinds = [s.kind for s in sigs]
    return Cluster(
        cluster_id=cid, locus=Locus("ref", "chr1", pos, pos, "+"),
        n_reads=len({s.read_id for s in sigs}),
        n_ins=kinds.count(SignatureKind.INS),
        n_split_ins=kinds.count(SignatureKind.SPLIT_INS),
        n_clip_left=sum(s.clip_side is ClipSide.LEFT for s in sigs),
        n_clip_right=sum(s.clip_side is ClipSide.RIGHT for s in sigs),
        insert_lengths=tuple(sorted(s.insert_length for s in sigs
                                    if s.kind is not SignatureKind.CLIP)),
        signature_ids=tuple(s.signature_id for s in sigs),
    )  # fmt: skip


def _call_tables(tmp_path: Path, *, cluster_ref: str = "ref",
                 hits_read_type: ReadType = ReadType.HIFI,
                 drop_hit: str | None = None) -> list[str]:  # fmt: skip
    """Two clusters: 3 reads with a synTE INS at 5000; 1 unhit clip at 20000."""
    spanning = [_ins(f"s{i}", f"r{i}", 5000 + i) for i in range(3)]
    clip = _clip("c1", "r9", 20_000)
    sigs = [*spanning, clip]
    hits = [*(_hit(s) for s in spanning), _hit(clip, NO_HIT)]
    hits = [h for h in hits if h.signature_id != drop_hit]
    paths = {k: tmp_path / f"in.{k}.tsv.gz" for k in ("clusters", "sigs", "hits")}
    write_signatures(paths["sigs"], sigs, read_type=ReadType.HIFI,
                     reference_id="ref")  # fmt: skip
    write_te_hits(paths["hits"], hits, read_type=hits_read_type,
                  reference_id="ref")  # fmt: skip
    write_clusters(paths["clusters"],
                   [_cluster("k1", 5000, spanning), _cluster("k2", 20_000, [clip])],
                   read_type=ReadType.HIFI, reference_id=cluster_ref)  # fmt: skip
    return ["--clusters", str(paths["clusters"]), "--signatures", str(paths["sigs"]),
            "--te-hits", str(paths["hits"])]  # fmt: skip


def test_lr_call_on_hand_written_tables(tmp_path: Path) -> None:
    out = tmp_path / "o"
    assert main(["lr-call", "--sample", "S", "-o", str(out),
                 *_call_tables(tmp_path)]) == 0  # fmt: skip
    meta, calls = read_calls(out / "S.calls.tsv.gz")
    calls = list(calls)
    assert meta["read_type"] == "hifi" and meta["reference_id"] == "ref"
    assert [(c.locus.start, c.tier, c.family, c.filters) for c in calls] == [
        (5001, CallTier.SPANNING, "synTE", ()),
        (20_000, CallTier.CLIP_ONLY, NO_HIT,
         (FilterReason.NO_TE_HIT, FilterReason.LOW_SUPPORT)),
    ]  # fmt: skip
    assert calls[0].n_reads == 3 and calls[0].median_identity == 0.95
    row = _timing(out)[0]
    assert row["stage"] == "lr-call" and row["status"] == "ok"
    assert row["records_in"] == "2" and row["records_out"] == "2"


def test_lr_call_thresholds_from_options(tmp_path: Path) -> None:
    out = tmp_path / "o"
    assert main(["lr-call", "--sample", "S", "-o", str(out), *_call_tables(tmp_path),
                 "--min-support", "4", "--min-te-identity", "0.96"]) == 0  # fmt: skip
    _, calls = read_calls(out / "S.calls.tsv.gz")
    first = next(iter(calls))
    assert first.family == NO_HIT  # every hit is below 0.96 identity
    assert FilterReason.LOW_SUPPORT in first.filters


def test_lr_call_refuses_existing_output(tmp_path: Path) -> None:
    args = ["lr-call", "--sample", "S", "-o", str(tmp_path / "o"),
            *_call_tables(tmp_path)]  # fmt: skip
    assert main(args) == 0
    assert main(args) == 1
    assert [r["status"] for r in _timing(tmp_path / "o")] == ["ok", "FileExistsError"]


def test_lr_call_rejects_mismatched_reference_id(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    args = _call_tables(tmp_path, cluster_ref="other")
    rc = main(["lr-call", "--sample", "S", "-o", str(tmp_path / "o"), *args])
    assert rc == 1
    assert args[1] in caplog.text and args[3] in caplog.text
    assert not (tmp_path / "o" / "S.calls.tsv.gz").exists()


def test_lr_call_rejects_mismatched_read_type(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    args = _call_tables(tmp_path, hits_read_type=ReadType.ONT)
    rc = main(["lr-call", "--sample", "S", "-o", str(tmp_path / "o"), *args])
    assert rc == 1
    assert args[5] in caplog.text


def test_lr_call_rejects_unannotated_signature(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    args = _call_tables(tmp_path, drop_hit="s1")
    rc = main(["lr-call", "--sample", "S", "-o", str(tmp_path / "o"), *args])
    assert rc == 1
    assert "not annotated" in caplog.text and "'s1'" in caplog.text


def test_lr_call_missing_input_is_an_error(tmp_path: Path) -> None:
    args = _call_tables(tmp_path)
    args[5] = str(tmp_path / "nope.tsv.gz")
    assert main(["lr-call", "--sample", "S", "-o", str(tmp_path / "o"), *args]) == 1


@needs_tools
def test_end_to_end_five_stages(
    lr_synthetic: tuple[Path, Path], tmp_path: Path
) -> None:
    fa, fq = lr_synthetic
    # The insert itself is the TE: bases 3500-4500 of the "fwd" read
    # (sample[6_500:14_500], insert at sample[10_000:11_000]).
    reads = fq.read_text().split("\n")
    fwd = reads[reads.index("@fwd") + 1]
    lib = tmp_path / "lib.fa"
    lib.write_text(f">synTE test insert\n{fwd[3_500:4_500]}\n")
    out = tmp_path / "out"
    so = ["--sample", "S", "-o", str(out)]
    assert main(["lr-align", *so, "--genome", str(fa), "--read-type", "hifi",
                 "--reads", str(fq)]) == 0  # fmt: skip
    assert main(["lr-signatures", *so, "--bam", str(out / "S.lr.bam"),
                 "--read-type", "hifi"]) == 0  # fmt: skip
    sigs = out / "S.signatures.tsv.gz"
    assert main(["lr-cluster", *so, "--signatures", str(sigs)]) == 0
    assert main(["lr-annotate", *so, "--signatures", str(sigs),
                 "--te-library", str(lib), "--threads", "2"]) == 0  # fmt: skip
    hits_path = out / "S.te_hits.tsv.gz"
    meta, hits = read_te_hits(hits_path)
    hits = list(hits)
    assert meta["read_type"] == "hifi" and meta["reference_id"] == "ref"
    assert any(h.family == "synTE" for h in hits)
    # hifi default min support is 2; the synthetic set has 3 reads.
    assert main(["lr-call", *so, "--clusters", str(out / "S.clusters.tsv.gz"),
                 "--signatures", str(sigs), "--te-hits", str(hits_path)]) == 0  # fmt: skip
    _, calls = read_calls(out / "S.calls.tsv.gz")
    passed = [c for c in calls if c.passed]
    assert len(passed) == 1, passed
    call = passed[0]
    assert call.tier is CallTier.SPANNING and call.family == "synTE"
    assert abs(call.locus.start - INS_AT) <= 50
    assert call.n_reads == 3

    rows = _timing(out)
    assert [r["stage"] for r in rows] == ["lr-align", "lr-signatures", "lr-cluster",
                                          "lr-annotate", "lr-call"]  # fmt: skip
    assert all(r["status"] == "ok" for r in rows)
    n_sigs = len(list(read_signatures(sigs)[1]))
    assert int(rows[3]["records_in"]) == n_sigs
    assert int(rows[3]["records_out"]) == len(hits)
    assert int(rows[4]["records_in"]) == int(rows[2]["records_out"])
