"""Pin record lifetimes and streaming without measuring allocator-dependent RSS."""

from io import StringIO
from pathlib import Path
import subprocess
import weakref

import pytest

from RelocaTE3.aligners import BlatBackend
from RelocaTE3.librelocate import RelocaTE


def test_completed_mate_records_released_before_parsing_next(tmp_path, monkeypatch):
    class Records(dict):
        pass

    previous = []
    written = []
    rt = RelocaTE()

    def parse(bam, mismatch_allowance):
        assert all(ref() is None for ref in previous), "Previous mate still retained"
        records = Records({bam.stem: {"seq": "ACGT", "qual": "IIII"}})
        previous.append(weakref.ref(records))
        return records

    def write(records, *args):
        written.extend(records)
        return len(records)

    monkeypatch.setattr(rt, "_parse_te_bam", parse)
    monkeypatch.setattr(rt, "_write_direction", write)
    assert (
        rt.write_trimmed_reads(
            "sample",
            [("left", Path("left.bam")), ("right", Path("right.bam"))],
            tmp_path,
        )
        == 2
    )
    assert written == ["left", "right"]
    assert all(ref() is None for ref in previous)


class IterationOnly(StringIO):
    def read(self, *args):
        raise AssertionError("FASTA must be iterated, not read into a whole string")

    def readlines(self, *args):
        raise AssertionError("FASTA must not become a line list")


def _psl(tmp_path, names):
    path = tmp_path / "hits.psl"
    path.write_text(
        "".join("\t".join(["0"] * 9 + [name] + ["0"] * 11) + "\n" for name in names)
    )
    return path


def test_query_fallback_streams_multiline_selected_sequences(tmp_path, monkeypatch):
    query = tmp_path / "query.fa"
    psl = _psl(tmp_path, ["keep/1", "keep/1", "missing"])
    query.write_text(">skip\nTTTT\n>keep/1 comment\nAC\nGT\n>skip2\nCCCC\n")
    real_open = open

    def guarded_open(path, *args, **kwargs):
        if Path(path) == query:
            return IterationOnly(query.read_text())
        return real_open(path, *args, **kwargs)

    monkeypatch.setattr("RelocaTE3.aligners.shutil.which", lambda tool: None)
    monkeypatch.setattr("RelocaTE3.aligners.open", guarded_open, raising=False)
    assert BlatBackend._query_sequences(query, psl, tmp_path) == {"keep/1": "ACGT"}


def test_query_seqtk_output_goes_to_file_and_is_streamed(tmp_path, monkeypatch):
    psl = _psl(tmp_path, ["keep/2"])
    monkeypatch.setattr("RelocaTE3.aligners.shutil.which", lambda tool: "seqtk")
    handles = []

    def fake_run(command, **kwargs):
        assert command[:2] == ["seqtk", "subseq"]
        assert not kwargs.get("capture_output")
        assert kwargs["check"] is True
        assert kwargs["stdout"].fileno() >= 0
        handles.append(kwargs["stdout"])
        assert set(Path(command[-1]).read_text().splitlines()) == {"keep/2"}
        kwargs["stdout"].write(">keep/2 comment\nAC\nGT\n")

    monkeypatch.setattr("RelocaTE3.aligners.subprocess.run", fake_run)
    assert BlatBackend._query_sequences(tmp_path / "query.fa", psl, tmp_path) == {
        "keep/2": "ACGT"
    }
    assert handles[0].closed


def test_query_seqtk_failure_propagates_and_cleans_partial_output(
    tmp_path, monkeypatch
):
    psl = _psl(tmp_path, ["keep"])
    monkeypatch.setattr("RelocaTE3.aligners.shutil.which", lambda tool: "seqtk")
    handles = []

    def fail(command, **kwargs):
        handles.append(kwargs["stdout"])
        kwargs["stdout"].write(">keep\nPARTIAL\n")
        raise subprocess.CalledProcessError(9, command)

    monkeypatch.setattr("RelocaTE3.aligners.subprocess.run", fail)
    with pytest.raises(subprocess.CalledProcessError, match="9"):
        BlatBackend._query_sequences(tmp_path / "query.fa", psl, tmp_path)
    assert handles[0].closed


def test_query_no_hits_never_opens_query_or_starts_seqtk(tmp_path, monkeypatch):
    psl = _psl(tmp_path, [])

    def unexpected(*args, **kwargs):
        raise AssertionError("No subprocess needed for empty PSL")

    monkeypatch.setattr("RelocaTE3.aligners.subprocess.run", unexpected)
    assert BlatBackend._query_sequences(tmp_path / "absent.fa", psl, tmp_path) == {}


def test_blat_loads_one_chunk_of_sequences_and_releases_before_sort(
    tmp_path, monkeypatch
):
    """SEQ is filled per chunk from that chunk's FASTA, never the whole query.

    RelocaTE3 used to merge every chunk's PSL and load the matched sequences of
    the whole side into one dictionary before conversion. Each chunk now loads
    only its own matched sequences, and drops them before its BAM sort.
    """
    import pysam

    class Sequences(dict):
        pass

    backend = BlatBackend()
    backend.BLAT_CHUNK_SEQS = 2  # 5 reads -> 3 chunks
    library = tmp_path / "te.fa"
    library.write_text(">TE\nACGTACGTAC\n")
    pysam.faidx(str(library))
    reads = tmp_path / "reads.fa"
    reads.write_text("".join(f">r{i}\nACGTACGTAC\n" for i in range(5)))

    loaded = []

    def sequences(query_fa, psl, tmpdir):
        names = [line.split("\t")[9] for line in open(psl)]
        value = Sequences({name: "ACGTACGTAC" for name in names})
        loaded.append((Path(query_fa).name, set(value), weakref.ref(value)))
        return value

    real_run = subprocess.run

    def run(cmd, **kwargs):
        if cmd[0] == "blat":
            # one full-length hit per read in this chunk
            names = [ln[1:].strip() for ln in open(cmd[2]) if ln.startswith(">")]
            Path(cmd[-1]).write_text(
                "".join(
                    "\t".join(
                        ["10", "0", "0", "0", "0", "0", "0", "0", "+", n, "10", "0",
                         "10", "TE", "10", "0", "10", "1", "10,", "0,", "0,"]
                    ) + "\n"
                    for n in names
                )
            )
            return None
        if cmd[:2] == ["samtools", "sort"]:
            assert all(ref() is None for _, _, ref in loaded), "sequences retained"
        return real_run(cmd, **kwargs)

    monkeypatch.setattr("RelocaTE3.aligners.shutil.which", lambda tool: None)
    monkeypatch.setattr("RelocaTE3.aligners.subprocess.run", run)
    monkeypatch.setattr(backend, "_query_sequences", sequences)
    result = backend._blat_side(
        str(library), reads, tmp_path / "out.left.bam", 1, str(tmp_path)
    )

    assert [names for _, names, _ in loaded] == [
        {"r0", "r1"}, {"r2", "r3"}, {"r4"}
    ], "each chunk loads only its own matched reads"
    assert all(name.startswith("chunk.") for name, _, _ in loaded)
    assert [c.n_reads for c in result.chunks] == [2, 2, 1]
    with pysam.AlignmentFile(str(result.bam)) as fh:
        assert sorted(r.query_name for r in fh) == ["r0", "r1", "r2", "r3", "r4"]
