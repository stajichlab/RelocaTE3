"""The BLAT TE search as three separately callable steps (Nextflow M4).

``relocaTE3 run`` searches each read side in 200k-read BLAT chunks inside one
process. To spread that search over many jobs the work is cut into three steps
that can run in different processes:

1. ``split_te_query``: one read side -> chunk FASTAs (+ the read count of each);
2. ``search_te_chunk``: one chunk FASTA -> one sorted, indexed chunk BAM;
3. ``merge_trim``: every side's chunk BAMs, in order -> the side BAMs and the
   trim outputs (``te_containing/``, ``flanking/``, ``te_portions/``).

The core rule is that the outputs are identical to ``identify_TE_reads``: text
byte for byte, BAMs record for record. A manifest (TSV) records the chunks so
the merge step can refuse a missing, extra or duplicated chunk.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pysam
import pytest

from RelocaTE3 import te_scatter
from RelocaTE3.aligners import (
    BLAT_CHUNK_SEQS,
    BlatBackend,
    TEChunk,
    TESide,
    search_te_chunk,
    split_te_query,
)
from RelocaTE3.librelocate import RelocaTE, merge_trim
from RelocaTE3.ReadLibrary import ReadLibrary

DATA = Path(__file__).parent / "data"
R1 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
RICETE = DATA / "TE_lib" / "RiceTE.fa"
READS_PER_SIDE = 75_910  # records in each Chr3 fixture FASTQ

needs_blat = pytest.mark.skipif(
    shutil.which("blat") is None or shutil.which("seqtk") is None,
    reason="blat/seqtk not on PATH (pixi env 'blat' or 'pipeline')",
)


def outputs(outdir: Path) -> dict[str, object]:
    """Every file under ``outdir``: BAMs as (@SQ lines, records), others as bytes.

    Indexes are skipped (they are derived from the BAM); @PG lines are not
    compared, as in validation/blat_memory/compare_outputs.sh.
    """
    found: dict[str, object] = {}
    for path in sorted(outdir.rglob("*")):
        if not path.is_file() or path.suffix == ".bai":
            continue
        rel = str(path.relative_to(outdir))
        if path.suffix == ".bam":
            with pysam.AlignmentFile(str(path), "rb") as fh:
                sq = [line for line in str(fh.header).splitlines() if line.startswith("@SQ")]
                found[rel] = (sq, [r.to_string() for r in fh.fetch(until_eof=True)])
        else:
            found[rel] = path.read_bytes()
    return found


def assert_same_outputs(got: Path, want: Path) -> None:
    a, b = outputs(got), outputs(want)
    assert sorted(a) == sorted(b)
    for rel in b:
        assert a[rel] == b[rel], f"{rel} differs"


@pytest.fixture(scope="module")
def reference(tmp_path_factory) -> Path:
    """``identify_TE_reads`` (what ``relocaTE3 run`` calls) at the default chunk size."""
    out = tmp_path_factory.mktemp("reference")
    reads = ReadLibrary([str(R1), str(R2)], "S")
    RelocaTE(TElib=str(RICETE), threads=4).identify_TE_reads(reads, out, te_aligner="blat")
    return out


def three_steps(work: Path, chunk_reads: int) -> list[TESide]:
    """Steps 1 and 2 for both sides, each chunk searched on its own."""
    sides = []
    for index, (direction, fq) in enumerate((("left", R1), ("right", R2))):
        chunks = split_te_query(fq, work / "split", f"S.{direction}", chunk_reads=chunk_reads)
        parts = [
            search_te_chunk(
                RICETE, fa, work / "bams" / f"{Path(fa).stem}.bam", mate=str(index + 1)
            )
            for fa, _ in chunks
        ]
        assert [p.n_reads for p in parts] == [n for _, n in chunks]
        assert sum(n for _, n in chunks) == READS_PER_SIDE
        sides.append(TESide(direction=direction, read_file=Path(fq), chunks=tuple(parts)))
    return sides


@pytest.fixture(scope="module")
def small_chunks(tmp_path_factory) -> list[TESide]:
    """Chunk BAMs of both sides at 5000 reads per chunk (16 chunks per side)."""
    return three_steps(tmp_path_factory.mktemp("small_chunks"), 5000)


@needs_blat
def test_three_steps_match_identify_te_reads(reference, small_chunks, tmp_path):
    assert [len(s.chunks) for s in small_chunks] == [16, 16]
    out = tmp_path / "out"
    merge_trim("S", small_chunks, out, RICETE, threads=2)
    assert_same_outputs(out, reference)
    # the merge step does not consume (delete) its inputs
    assert all(Path(c.bam).exists() for s in small_chunks for c in s.chunks)


@needs_blat
@pytest.mark.parametrize("chunk_reads", [30_000, BLAT_CHUNK_SEQS])
def test_chunk_size_does_not_change_outputs(reference, chunk_reads, tmp_path):
    sides = three_steps(tmp_path / "work", chunk_reads)
    expected = -(-READS_PER_SIDE // chunk_reads)
    assert [len(s.chunks) for s in sides] == [expected, expected]
    merge_trim("S", sides, tmp_path / "out", RICETE, threads=2)
    assert_same_outputs(tmp_path / "out", reference)


@needs_blat
def test_merge_trim_refuses_reads_the_chunks_do_not_cover(small_chunks, tmp_path):
    short = [TESide(s.direction, s.read_file, s.chunks[:-1]) for s in small_chunks]
    with pytest.raises(ValueError, match="more records than its TE chunks cover"):
        merge_trim("S", short, tmp_path / "out", RICETE)


def test_default_chunk_size_is_relocate2s():
    """200k reads per chunk, RelocaTE2's ``fastq_split -s`` (relocaTE2.py:91)."""
    assert BLAT_CHUNK_SEQS == BlatBackend.BLAT_CHUNK_SEQS == 200_000


def _fastq(path: Path, n: int) -> Path:
    path.write_text("".join(f"@r{i}/1\nACGTACGTAC\n+\nIIIIIIIIII\n" for i in range(n)))
    return path


@pytest.mark.parametrize("seqtk", [True, False], ids=["seqtk", "python"])
def test_split_te_query_chunks_and_counts(tmp_path, monkeypatch, seqtk):
    if seqtk and shutil.which("seqtk") is None:
        pytest.skip("seqtk not on PATH")
    if not seqtk:
        monkeypatch.setattr("RelocaTE3.aligners.shutil.which", lambda tool: None)
    fq = _fastq(tmp_path / "reads.fq", 7)
    chunks = split_te_query(fq, tmp_path / "out", "S.left", chunk_reads=3)
    assert [Path(p).name for p, _ in chunks] == [
        "chunk.S.left.00000.fa", "chunk.S.left.00001.fa", "chunk.S.left.00002.fa"
    ]
    assert [n for _, n in chunks] == [3, 3, 1]
    assert Path(chunks[2][0]).read_text() == ">r6/1\nACGTACGTAC\n"
    # only the chunks are left behind (no whole-side FASTA)
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == [
        Path(p).name for p, _ in chunks
    ]


def test_split_te_query_refuses_existing_chunks(tmp_path):
    fq = _fastq(tmp_path / "reads.fq", 2)
    split_te_query(fq, tmp_path / "out", "S.left", chunk_reads=1)
    with pytest.raises(FileExistsError, match="chunk.S.left"):
        split_te_query(fq, tmp_path / "out", "S.left", chunk_reads=1)


def test_split_te_query_empty_input_gives_no_chunks(tmp_path):
    fq = _fastq(tmp_path / "reads.fq", 0)
    assert split_te_query(fq, tmp_path / "out", "S.left") == []


# ---------------------------------------------------------------------------
# manifest
# ---------------------------------------------------------------------------


def _rows(side, n, mate, read_file="R1.fq.gz", reads=10):
    return [
        te_scatter.ManifestRow(side, i, f"chunk.S.{side}.{i:05d}.fa", reads, mate, read_file)
        for i in range(n)
    ]


def _bams(tmp_path, rows):
    return [tmp_path / te_scatter.chunk_bam_name(r.chunk) for r in rows]


def test_manifest_round_trip(tmp_path):
    rows = _rows("left", 2, "1") + _rows("right", 2, "2", "R2.fq.gz")
    path = tmp_path / "manifest.tsv"
    te_scatter.write_manifest(path, rows)
    assert path.read_text().splitlines()[0] == "side\tindex\tchunk\tn_reads\tmate\tread_file"
    assert path.read_text().splitlines()[1] == "left\t0\tchunk.S.left.00000.fa\t10\t1\tR1.fq.gz"
    assert te_scatter.read_manifest(path) == rows


def test_manifest_single_end_mate_is_dot(tmp_path):
    rows = _rows("left", 1, None)
    path = tmp_path / "manifest.tsv"
    te_scatter.write_manifest(path, rows)
    assert path.read_text().splitlines()[1].split("\t")[4] == "."
    assert te_scatter.read_manifest(path) == rows


def test_manifest_rejects_bad_header(tmp_path):
    path = tmp_path / "manifest.tsv"
    path.write_text("side\tindex\n")
    with pytest.raises(ValueError, match="not a TE chunk manifest"):
        te_scatter.read_manifest(path)


def test_sides_from_manifest_pairs_bams_in_index_order(tmp_path):
    rows = _rows("left", 3, "1") + _rows("right", 2, "2", "R2.fq.gz", reads=4)
    bams = _bams(tmp_path, rows)
    sides = te_scatter.sides_from_manifest(rows, list(reversed(bams)), ["R1.fq.gz", "R2.fq.gz"])
    assert [s.direction for s in sides] == ["left", "right"]
    assert [s.read_file for s in sides] == [Path("R1.fq.gz"), Path("R2.fq.gz")]
    assert sides[0].chunks == tuple(TEChunk(bam=b, n_reads=10) for b in bams[:3])
    assert sides[1].chunks == tuple(TEChunk(bam=b, n_reads=4) for b in bams[3:])


def test_sides_from_manifest_side_with_no_chunks_is_an_empty_side(tmp_path):
    rows = _rows("left", 1, "1")
    sides = te_scatter.sides_from_manifest(rows, _bams(tmp_path, rows), ["a.fq", "b.fq"])
    assert [(s.direction, len(s.chunks)) for s in sides] == [("left", 1), ("right", 0)]


@pytest.mark.parametrize(
    "case, message",
    [
        ("missing", "missing chunk BAM.*chunk.S.left.00001.bam"),
        ("extra", "not in the manifest.*other.bam"),
        ("duplicate_bam", "given more than once.*chunk.S.left.00000.bam"),
        ("gap", "indices for side 'left' must be 0..1"),
        ("duplicate_row", "indices for side 'left' must be 0..2"),
        ("unknown_side", "side 'right'.*one read file"),
        ("mate", "mate '2' for side 'left'.*expected '1'"),
    ],
)
def test_sides_from_manifest_refuses_inconsistent_chunks(tmp_path, case, message):
    rows = _rows("left", 2, "1")
    reads = ["R1.fq.gz", "R2.fq.gz"]
    bams = _bams(tmp_path, rows)
    if case == "missing":
        bams = bams[:1]
    elif case == "extra":
        bams.append(tmp_path / "other.bam")
    elif case == "duplicate_bam":
        bams.append(tmp_path / "elsewhere" / bams[0].name)
    elif case == "gap":
        rows[1] = te_scatter.ManifestRow("left", 2, rows[1].chunk, 10, "1", "R1.fq.gz")
    elif case == "duplicate_row":
        rows.append(rows[1])
    elif case == "unknown_side":
        rows = _rows("right", 1, None)
        bams = _bams(tmp_path, rows)
        reads = ["R1.fq.gz"]
    elif case == "mate":
        rows = _rows("left", 2, "2")
    with pytest.raises(ValueError, match=message):
        te_scatter.sides_from_manifest(rows, bams, reads)


def test_chunk_bam_name():
    assert te_scatter.chunk_bam_name("chunk.S.left.00003.fa") == "chunk.S.left.00003.bam"
    assert te_scatter.chunk_bam_name("dir/chunk.S.left.00003.fa") == "chunk.S.left.00003.bam"
