"""Per-chunk TE trimming must equal whole-file trimming, byte for byte.

RelocaTE2 runs its TE search and trim per 200k-read FASTQ chunk
(relocaTE2.py ``--split`` + fastq_split.pl; relocaTE_trim.py per chunk), so it
never holds more than one chunk's TE hits. RelocaTE3 ran BLAT per chunk but then
merged everything and trimmed from one side-wide dictionary that held every
TE-hit read with its sequence and qualities (~16 of 19 GiB on riceTElib 30x).

The chunked path keeps RelocaTE3's existing output order (first TE hit in
coordinate order), because ``bwa aln``/``samse`` picks among equal hits with a
sequential RNG, so flanking-read order is not cosmetic.
"""

from __future__ import annotations

import gc
import gzip
import random
import shutil
import weakref
from itertools import islice
from pathlib import Path

import pysam
import pytest

from RelocaTE3.aligners import BlatBackend, ChunkedTEAlignment, TEChunk
from RelocaTE3.librelocate import RelocaTE, reverse_complement
from RelocaTE3.ReadLibrary import ReadLibrary

DATA = Path(__file__).parent / "data"
R1 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_1.fq.gz"
R2 = DATA / "sim_reads" / "MSU7.Chr3_2M.ALL_reads_6X_100_500_2.fq.gz"
RICETE = DATA / "TE_lib" / "RiceTE.fa"

TE_LEN = 100
READ_LEN = 40
# (cigar, reference_start) templates that exercise every trim branch:
# 5' junctions, 3' junctions, fully internal (:middle), and unclassified hits.
TEMPLATES = [
    ([(4, 15), (0, 25)], 0),  # 15S25M at TE start -> 5' junction
    ([(0, 25), (4, 15)], 0),  # 25M15S at TE start
    ([(0, 25), (4, 15)], 75),  # 25M15S ending at TE end -> 3' junction
    ([(4, 15), (0, 25)], 75),  # 15S25M at TE end
    ([(0, 40)], 30),  # fully internal -> :middle
    ([(4, 5), (0, 30), (4, 5)], 40),  # internal, read ends clipped -> unclassified
    ([(4, 12), (0, 28)], 1),  # near TE start, shorter flank
]


def _header():
    return {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "teA", "LN": TE_LEN}, {"SN": "teB", "LN": TE_LEN}],
    }


def _segment(header, name, seq, te, template, reverse, nm):
    cigar, start = template
    seg = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header))
    seg.query_name = name
    seg.flag = 16 if reverse else 0
    seg.reference_id = te
    seg.reference_start = start
    seg.mapping_quality = 255
    seg.cigartuples = cigar
    seg.query_sequence = reverse_complement(seq) if reverse else seq
    seg.set_tag("NM", nm)
    return seg


def _write_sorted_bam(header, segments, path: Path) -> Path:
    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for seg in segments:
            out.write(seg)
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))
    unsorted.unlink()
    return path


def _synthetic_side(tmp: Path, side: str, mate: str, n_reads=61, chunk=7, seed=7):
    """One mate side: a FASTQ, per-chunk TE BAMs, and the whole-side BAM.

    Hits are drawn from a few templates on two TEs so equal coordinates (and
    equal-scoring family ties) recur both within and across chunks. FASTQ names
    use Illumina's separate mate field, so restoration must canonicalize.
    """
    rng = random.Random(seed)
    header = _header()
    reads = []
    for i in range(n_reads):
        seq = "".join(rng.choice("ACGT") for _ in range(READ_LEN))
        qual = "".join(rng.choice("#5?FI") for _ in range(READ_LEN))
        reads.append((f"r{i:03d}", seq, qual))

    fastq = tmp / f"{side}.fq"
    with open(fastq, "w") as fh:
        for name, seq, qual in reads:
            fh.write(f"@{name} {mate}:N:0:ACGT\n{seq}\n+\n{qual}\n")

    per_chunk = []
    for c in range(0, n_reads, chunk):
        segments = []
        for name, seq, _ in reads[c : c + chunk]:
            if rng.random() < 0.15:
                continue  # no TE hit
            for _ in range(rng.randint(1, 3)):
                segments.append(
                    _segment(
                        header,
                        f"{name}/{mate}",
                        seq,
                        rng.randint(0, 1),
                        rng.choice(TEMPLATES),
                        rng.random() < 0.5,
                        rng.choice([0, 0, 1, 3]),
                    )
                )
        rng.shuffle(segments)  # aligner output order is arbitrary
        per_chunk.append((segments, len(reads[c : c + chunk])))

    chunk_dir = tmp / f"{side}.te_chunks"
    chunk_dir.mkdir()
    chunks = tuple(
        TEChunk(
            bam=_write_sorted_bam(header, segs, chunk_dir / f"chunk.{i:05d}.bam"),
            n_reads=n,
        )
        for i, (segs, n) in enumerate(per_chunk)
    )
    whole = _write_sorted_bam(
        header, [s for segs, _ in per_chunk for s in segs], tmp / f"S.{side}.bam"
    )
    return fastq, ChunkedTEAlignment(bam=whole, chunks=chunks, chunk_dir=chunk_dir)


def _outputs(outdir: Path) -> dict[str, bytes]:
    return {
        str(p.relative_to(outdir)): p.read_bytes()
        for p in sorted(outdir.rglob("*"))
        if p.is_file()
    }


def _paired_fixture(tmp_path):
    left_fq, left = _synthetic_side(tmp_path, "left", "1", seed=7)
    right_fq, right = _synthetic_side(tmp_path, "right", "2", seed=11)
    return [left_fq, right_fq], [("left", left), ("right", right)]


def test_chunked_trim_is_byte_identical_to_whole_file_trim(tmp_path):
    fastqs, sides = _paired_fixture(tmp_path)

    whole_out = tmp_path / "whole"
    n_whole = RelocaTE().write_trimmed_reads(
        "S",
        [(direction, aln.bam) for direction, aln in sides],
        whole_out,
        source_fastqs=fastqs,
    )
    chunked_out = tmp_path / "chunked"
    n_chunked = RelocaTE().write_trimmed_reads_chunked(
        "S", sides, chunked_out, source_fastqs=fastqs
    )

    whole, chunked = _outputs(whole_out), _outputs(chunked_out)
    assert n_whole > 0
    assert n_chunked == n_whole
    assert sorted(chunked) == sorted(whole)
    for path in whole:
        assert chunked[path] == whole[path], f"{path} differs"
    # The fixture must actually exercise the classified branches.
    assert whole["te_containing/S.read_repeat_name.txt"].count(b"\n") > 10
    assert b"TE_best_families" in whole["te_containing/S.read_repeat_name.txt"]


def test_chunked_trim_holds_one_chunk_of_records_at_a_time(tmp_path, monkeypatch):
    fastqs, sides = _paired_fixture(tmp_path)
    chunk_bams = {str(c.bam) for _, aln in sides for c in aln.chunks}
    whole_bams = {str(aln.bam) for _, aln in sides}

    class Records(dict):  # defaultdict cannot be weakly referenced
        pass

    rt = RelocaTE()
    original = rt._parse_te_bam
    seen = []
    previous = []

    def parse(bam, *args, **kwargs):
        gc.collect()
        assert all(ref() is None for ref in previous), "previous chunk still held"
        coord = Records(original(bam, *args, **kwargs))
        seen.append((str(bam), len(coord)))
        previous.append(weakref.ref(coord))
        return coord

    monkeypatch.setattr(rt, "_parse_te_bam", parse)
    rt.write_trimmed_reads_chunked("S", sides, tmp_path / "out", source_fastqs=fastqs)

    parsed = [bam for bam, _ in seen]
    assert set(parsed) == chunk_bams and len(parsed) == len(chunk_bams)
    assert not whole_bams & set(parsed), "the whole-side BAM must not be parsed"
    assert max(n for _, n in seen) <= 7


def test_chunked_trim_rejects_fastq_that_does_not_match_chunks(tmp_path):
    fastqs, sides = _paired_fixture(tmp_path)
    direction, aln = sides[0]
    short = ChunkedTEAlignment(
        bam=aln.bam,
        chunks=(TEChunk(bam=aln.chunks[0].bam, n_reads=aln.chunks[0].n_reads),),
        chunk_dir=aln.chunk_dir,
    )
    with pytest.raises(ValueError, match="chunk"):
        RelocaTE().write_trimmed_reads_chunked(
            "S", [(direction, short)], tmp_path / "out", source_fastqs=fastqs[:1]
        )


def _subset_fastq(src: Path, dest: Path, n: int) -> Path:
    with gzip.open(src, "rt") as fh, open(dest, "w") as out:
        out.writelines(islice(fh, 4 * n))
    return dest


def _records(bam: Path) -> list[str]:
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        return [rec.to_string() for rec in fh.fetch(until_eof=True)]


@pytest.mark.skipif(
    shutil.which("blat") is None or shutil.which("seqtk") is None,
    reason="blat/seqtk not on PATH (pixi env 'blat')",
)
def test_blat_chunked_run_matches_single_chunk_run(tmp_path, monkeypatch):
    """Real BLAT: many small chunks give the same BAM and trim outputs as one."""
    left = _subset_fastq(R1, tmp_path / "r1.fq", 1500)
    right = _subset_fastq(R2, tmp_path / "r2.fq", 1500)
    reads = ReadLibrary([str(left), str(right)], "S")

    # Reference: one BLAT chunk per side, whole-file trim.
    single = tmp_path / "single"
    single.mkdir()
    bams = BlatBackend(threads=2).map_te_library(reads, str(RICETE), str(single))
    RelocaTE().write_trimmed_reads(
        "S",
        list(zip(["left", "right"], bams)),
        single,
        source_fastqs=reads.file_set,
    )

    # Chunked: 4 chunks per side, per-chunk trim through the pipeline entry.
    monkeypatch.setattr(BlatBackend, "BLAT_CHUNK_SEQS", 400)
    real = BlatBackend.map_te_library_chunked
    chunk_counts = []

    def spy(self, *args, **kwargs):
        sides = real(self, *args, **kwargs)
        chunk_counts.extend(len(s.chunks) for s in sides)
        assert all(sum(c.n_reads for c in s.chunks) == 1500 for s in sides)
        return sides

    monkeypatch.setattr(BlatBackend, "map_te_library_chunked", spy)
    chunked = tmp_path / "chunked"
    chunked.mkdir()
    RelocaTE(TElib=str(RICETE), threads=2).identify_TE_reads(
        reads, chunked, te_aligner="blat"
    )

    assert chunk_counts == [4, 4]
    for side in ("left", "right"):
        assert _records(chunked / f"S.{side}.bam") == _records(single / f"S.{side}.bam")
        assert not (chunked / f"S.{side}.te_chunks").exists(), "chunk BAMs left behind"
    trim_dirs = ("flanking", "te_containing", "te_portions")
    want = {k: v for k, v in _outputs(single).items() if k.split("/")[0] in trim_dirs}
    got = {k: v for k, v in _outputs(chunked).items() if k.split("/")[0] in trim_dirs}
    assert sorted(got) == sorted(want)
    for path in want:
        assert got[path] == want[path], f"{path} differs"
    assert want["te_containing/S.read_repeat_name.txt"]
