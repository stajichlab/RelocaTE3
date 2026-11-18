"""Long-read alignment: minimap2 command, piped align, input BAM checks.

Unit tests build BAMs with pysam; the integration test aligns synthetic
reads (a known 1 kb insert) with the real minimap2 and samtools.
"""

from __future__ import annotations

import logging
import os
import random
import shutil
import subprocess
from pathlib import Path

import pysam
import pytest

from RelocaTE3.longread import align
from RelocaTE3.longread.align import (
    align_reads, check_alignment_input, minimap2_command,
)  # fmt: skip
from RelocaTE3.longread.model import Reference, ReadType, SignatureKind
from RelocaTE3.longread.placements import iter_placements
from RelocaTE3.longread.reference import LinearFastaReferenceSpace
from RelocaTE3.longread.signatures import SignatureParams, collect_signatures

HAVE_TOOLS = shutil.which("minimap2") is not None and shutil.which("samtools")
_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


# --- minimap2_command -------------------------------------------------------


@pytest.mark.parametrize(
    ("read_type", "preset"),
    [(ReadType.ONT, "map-ont"), (ReadType.ONT_HQ, "lr:hq"),
     (ReadType.HIFI, "map-hifi"), (ReadType.CLR, "map-pb")],
)  # fmt: skip
def test_minimap2_command_exact(read_type: ReadType, preset: str) -> None:
    cmd = minimap2_command(Path("g.fa"), [Path("a.fq"), Path("b.fq.gz")], read_type, 4)
    assert cmd == [
        "minimap2", "-a", "-x", preset, "-Y", "--MD", "--secondary=no",
        "-t", "4", "g.fa", "a.fq", "b.fq.gz",
    ]  # fmt: skip


def test_minimap2_command_rejects_zero_threads() -> None:
    with pytest.raises(ValueError, match="threads"):
        minimap2_command(Path("g.fa"), [Path("a.fq")], ReadType.ONT, 0)


# --- check_alignment_input --------------------------------------------------


def _space(tmp_path: Path) -> LinearFastaReferenceSpace:
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "A" * 1000 + "\n>chr2\n" + "C" * 500 + "\n")
    return LinearFastaReferenceSpace(Reference("ref", str(fa)))


def _bam(
    path: Path,
    sq: list[tuple[str, int]],
    *,
    sorted_: bool = True,
    index: bool = True,
    fasta: Path | None = None,
) -> Path:
    hd = {"VN": "1.6", "SO": "coordinate" if sorted_ else "unsorted"}
    header = {"HD": hd, "SQ": [{"SN": n, "LN": ln} for n, ln in sq]}
    if fasta is None:
        out = pysam.AlignmentFile(str(path), "wb", header=header)
    else:
        out = pysam.AlignmentFile(
            str(path), "wc", header=header, reference_filename=str(fasta)
        )
    with out:
        seg = pysam.AlignedSegment(out.header)
        seg.query_name = "r1"
        seg.query_sequence = "A" * 50
        seg.reference_id = 0
        seg.reference_start = 100
        seg.mapping_quality = 60
        seg.cigartuples = [(0, 50)]
        seg.query_qualities = pysam.qualitystring_to_array("I" * 50)
        out.write(seg)
    if index:
        pysam.index(str(path))
    return path


GOOD_SQ = [("chr1", 1000), ("chr2", 500)]


def test_check_accepts_good_bam_any_sq_order(tmp_path: Path) -> None:
    space = _space(tmp_path)
    check_alignment_input(_bam(tmp_path / "a.bam", GOOD_SQ), space)
    check_alignment_input(_bam(tmp_path / "b.bam", GOOD_SQ[::-1]), space)


def test_check_rejects_unsorted(tmp_path: Path) -> None:
    bam = _bam(tmp_path / "u.bam", GOOD_SQ, sorted_=False, index=False)
    with pytest.raises(ValueError, match="samtools sort"):
        check_alignment_input(bam, _space(tmp_path))


def test_check_rejects_unindexed(tmp_path: Path) -> None:
    bam = _bam(tmp_path / "n.bam", GOOD_SQ, index=False)
    with pytest.raises(ValueError, match="samtools index"):
        check_alignment_input(bam, _space(tmp_path))


def test_check_rejects_wrong_sq_name(tmp_path: Path) -> None:
    bam = _bam(tmp_path / "w.bam", [("chr1", 1000), ("chrX", 500)])
    with pytest.raises(ValueError, match="samtools") as err:
        check_alignment_input(bam, _space(tmp_path))
    msg = str(err.value)
    assert "chr2" in msg and "chrX" in msg
    assert "missing" in msg and "extra" in msg


def test_check_rejects_wrong_sq_length(tmp_path: Path) -> None:
    bam = _bam(tmp_path / "l.bam", [("chr1", 1000), ("chr2", 499)])
    with pytest.raises(ValueError, match="samtools") as err:
        check_alignment_input(bam, _space(tmp_path))
    msg = str(err.value)
    assert "chr2" in msg and "499" in msg and "500" in msg


def test_check_bam_needs_bai_or_csi_not_crai(tmp_path: Path) -> None:
    bam = _bam(tmp_path / "c.bam", GOOD_SQ, index=False)
    Path(str(bam) + ".crai").write_bytes(b"")
    with pytest.raises(ValueError, match="samtools index"):
        check_alignment_input(bam, _space(tmp_path))


def test_check_bam_accepts_csi(tmp_path: Path) -> None:
    bam = _bam(tmp_path / "c.bam", GOOD_SQ, index=False)
    pysam.index("-c", str(bam))
    assert Path(str(bam) + ".csi").is_file()
    check_alignment_input(bam, _space(tmp_path))


def test_check_accepts_indexed_cram(tmp_path: Path) -> None:
    space = _space(tmp_path)
    fa = tmp_path / "ref.fa"
    cram = _bam(tmp_path / "g.cram", GOOD_SQ, fasta=fa)
    assert Path(str(cram) + ".crai").is_file()
    check_alignment_input(cram, space, reference_fasta=fa)


def test_check_cram_needs_crai_not_bai(tmp_path: Path) -> None:
    space = _space(tmp_path)
    fa = tmp_path / "ref.fa"
    cram = _bam(tmp_path / "g.cram", GOOD_SQ, fasta=fa, index=False)
    Path(str(cram) + ".bai").write_bytes(b"")
    with pytest.raises(ValueError, match="samtools index"):
        check_alignment_input(cram, space, reference_fasta=fa)


def test_check_cram_needs_reference(tmp_path: Path) -> None:
    cram = tmp_path / "x.cram"
    cram.write_bytes(b"")
    with pytest.raises(ValueError, match="reference_fasta"):
        check_alignment_input(cram, _space(tmp_path))


# --- align_reads guards -----------------------------------------------------


def test_align_refuses_existing_output(tmp_path: Path) -> None:
    out = tmp_path / "x.bam"
    out.write_bytes(b"")
    with pytest.raises(FileExistsError):
        align_reads(tmp_path / "g.fa", [tmp_path / "r.fq"], ReadType.ONT, out, 1)


@pytest.mark.parametrize("missing", ["minimap2", "samtools"])
def test_align_missing_tool(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, missing: str
) -> None:
    real = shutil.which
    monkeypatch.setattr(
        align.shutil, "which", lambda t: None if t == missing else real(t) or t
    )
    with pytest.raises(RuntimeError, match=missing):
        align_reads(
            tmp_path / "g.fa", [tmp_path / "r.fq"], ReadType.ONT,
            tmp_path / "x.bam", 1,
        )  # fmt: skip


def test_align_rejects_zero_threads(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="threads"):
        align_reads(
            tmp_path / "g.fa", [tmp_path / "r.fq"], ReadType.ONT,
            tmp_path / "x.bam", 0,
        )  # fmt: skip


@pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH")
def test_align_failing_minimap2_cleans_up(tmp_path: Path) -> None:
    reads = tmp_path / "r.fq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "x.bam"
    with pytest.raises(RuntimeError, match="minimap2") as err:
        align_reads(tmp_path / "missing.fa", [reads], ReadType.ONT, out, 1)
    assert "failed to open file" in str(err.value)  # minimap2 stderr
    assert sorted(p.name for p in tmp_path.iterdir()) == ["r.fq"]


# --- integration ------------------------------------------------------------

INS_AT = 10_000
INS_LEN = 1_000


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    d = tmp_path_factory.mktemp("lr_align")
    rng = random.Random(20261117)
    genome = _seq(20_000, rng)
    insert = _seq(INS_LEN, rng)
    sample = genome[:INS_AT] + insert + genome[INS_AT:]
    fa = d / "genome.fa"
    fa.write_text(
        ">chr1\n"
        + "\n".join(genome[i : i + 80] for i in range(0, len(genome), 80))
        + "\n"
    )
    reads = {
        "fwd": sample[6_500:14_500],
        "rev": _rc(sample[6_500:14_500]),
        "off": sample[5_200:13_200],
    }
    fq = d / "reads.fq"
    fq.write_text("".join(f"@{n}\n{s}\n+\n{'I' * len(s)}\n" for n, s in reads.items()))
    return fa, fq


@pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH")
@pytest.mark.parametrize("read_type", list(ReadType))
def test_align_integration(
    synthetic: tuple[Path, Path], tmp_path: Path, read_type: ReadType
) -> None:
    fa, fq = synthetic
    out = tmp_path / f"{read_type.value}.bam"
    got = align_reads(fa, [fq], read_type, out, 2)
    assert got == out and out.is_file()
    assert Path(str(out) + ".bai").is_file()
    assert sorted(p.name for p in tmp_path.iterdir()) == [out.name, out.name + ".bai"]
    space = LinearFastaReferenceSpace(Reference("ref", str(fa)))
    check_alignment_input(out, space)

    sigs = list(
        collect_signatures(iter_placements(out, "ref", min_mapq=20), SignatureParams())
    )
    hits = [
        s
        for s in sigs
        if s.kind in (SignatureKind.INS, SignatureKind.SPLIT_INS)
        and abs(s.locus.start - INS_AT) <= 20
        and abs(s.insert_length - INS_LEN) <= 0.1 * INS_LEN
    ]
    print(
        read_type.value,
        sorted({(str(s.kind), s.locus.start, s.insert_length) for s in sigs}),
    )
    assert hits, f"{read_type}: no INS/SPLIT_INS near {INS_AT}: {sigs}"


# --- align_reads failure handling (real minimap2) --------------------------


@pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH")
def test_align_blames_samtools_when_sort_dies_first(
    synthetic: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Sort dies first; minimap2 then gets SIGPIPE. Blame samtools."""
    fa, fq = synthetic
    real = shutil.which("samtools")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    fake = bindir / "samtools"
    fake.write_text(
        "#!/usr/bin/env bash\n"
        'if [ "$1" = sort ]; then\n'
        "  head -c 100 >/dev/null\n"  # read stdin briefly
        "  exec 0<&-\n"  # close it: minimap2's next write gets SIGPIPE
        "  sleep 0.3\n"
        "  echo 'samtools sort: No space left on device' >&2\n"
        "  exit 1\n"
        "fi\n"
        f'exec {real} "$@"\n'
    )
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    with pytest.raises(RuntimeError) as err:
        align_reads(fa, [fq], ReadType.HIFI, out_dir / "x.bam", 1)
    msg = str(err.value)
    assert msg.startswith("samtools sort failed"), msg
    assert "No space left on device" in msg
    assert list(out_dir.iterdir()) == []


@pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH")
def test_align_interrupt_kills_children_and_cleans_up(
    synthetic: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fa, fq = synthetic
    procs: list[subprocess.Popen] = []
    raised = [False]

    class InterruptedPopen(subprocess.Popen):  # type: ignore[type-arg]
        def __init__(self, *args, **kwargs) -> None:  # noqa: ANN002, ANN003
            super().__init__(*args, **kwargs)
            procs.append(self)

        def wait(self, timeout: float | None = None) -> int:
            if not raised[0]:
                raised[0] = True
                raise KeyboardInterrupt
            return super().wait(timeout)

    monkeypatch.setattr(align.subprocess, "Popen", InterruptedPopen)
    with pytest.raises(KeyboardInterrupt):
        align_reads(fa, [fq], ReadType.HIFI, tmp_path / "x.bam", 1)
    assert len(procs) == 2
    assert all(p.poll() is not None for p in procs)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH")
def test_align_removes_stale_sort_temp_files(
    synthetic: tuple[Path, Path], tmp_path: Path
) -> None:
    fa, fq = synthetic
    stale = tmp_path / "x.sorttmp.0000.bam"
    stale.write_bytes(b"stale")
    align_reads(fa, [fq], ReadType.HIFI, tmp_path / "x.bam", 1)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["x.bam", "x.bam.bai"]


@pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH")
def test_align_sort_failure_removes_sort_temp_files(
    synthetic: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Temp files written under the -T prefix are removed on failure."""
    fa, fq = synthetic
    real = shutil.which("samtools")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    fake = bindir / "samtools"
    # Write a temp file at the -T prefix (argument after -T), then fail.
    fake.write_text(
        "#!/usr/bin/env bash\n"
        'if [ "$1" = sort ]; then\n'
        '  prev=""; for a in "$@"; do\n'
        '    if [ "$prev" = -T ]; then echo x > "$a.0000.bam"; fi; prev="$a"\n'
        "  done\n"
        "  cat >/dev/null; echo 'sort: boom' >&2; exit 1\n"
        "fi\n"
        f'exec {real} "$@"\n'
    )
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    with pytest.raises(RuntimeError, match="boom"):
        align_reads(fa, [fq], ReadType.HIFI, out_dir / "x.bam", 1)
    assert list(out_dir.iterdir()) == []


@pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH")
def test_align_removes_orphan_index(
    synthetic: tuple[Path, Path], tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """An index without its BAM is stale (died between renames): remove it."""
    fa, fq = synthetic
    orphan = tmp_path / "x.bam.bai"
    orphan.write_bytes(b"stale")
    with caplog.at_level(logging.WARNING, logger="RelocaTE3.longread.align"):
        align_reads(fa, [fq], ReadType.HIFI, tmp_path / "x.bam", 1)
    assert "x.bam.bai" in caplog.text and "stale" in caplog.text
    assert orphan.read_bytes() != b"stale"  # replaced by the real index
    pysam.AlignmentFile(str(tmp_path / "x.bam")).check_index()


@pytest.mark.skipif(not HAVE_TOOLS, reason="minimap2/samtools not on PATH")
def test_align_interrupt_during_index_cleans_up(
    synthetic: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fa, fq = synthetic

    def interrupted_run(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        raise KeyboardInterrupt

    monkeypatch.setattr(align.subprocess, "run", interrupted_run)
    with pytest.raises(KeyboardInterrupt):
        align_reads(fa, [fq], ReadType.HIFI, tmp_path / "x.bam", 1)
    assert list(tmp_path.iterdir()) == []
