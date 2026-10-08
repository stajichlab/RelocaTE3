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


def test_blat_sequences_released_before_bam_sort(tmp_path, monkeypatch):
    import pysam

    class Sequences(dict):
        pass

    references = []
    sam_output = StringIO()
    real_open = open

    class Writer:
        def __enter__(self):
            return sam_output

        def __exit__(self, *args):
            pass

    def capture_sam(path, *args, **kwargs):
        if Path(path) == tmp_path / "aln.sam":
            return Writer()
        return real_open(path, *args, **kwargs)

    backend = BlatBackend()
    library = tmp_path / "te.fa"
    library.write_text(">TE\nACGT\n")
    pysam.faidx(str(library))
    monkeypatch.setattr(backend, "_write_query_chunks", lambda *args: [])

    def sequences(*args):
        value = Sequences({"read": "ACGT"})
        references.append(weakref.ref(value))
        return value

    def convert(lines, sequences):
        assert sequences == {"read": "ACGT"}
        record = "read\t0\tTE\t1\t60\t4M\t*\t0\t0\tACGT\tIIII"
        yield record
        assert sam_output.getvalue().endswith(record + "\n"), (
            "Writer buffered the iterator"
        )

    def sort(sam, out_bam, threads):
        assert all(ref() is None for ref in references)
        assert "ACGT\tIIII" in sam_output.getvalue()
        return out_bam

    monkeypatch.setattr(backend, "_query_sequences", sequences)
    monkeypatch.setattr("RelocaTE3.aligners.open", capture_sam, raising=False)
    monkeypatch.setattr("RelocaTE3.aligners._iter_psl_to_sam", convert)
    monkeypatch.setattr("RelocaTE3.aligners._sam_to_sorted_mapped_bam", sort)
    assert backend._blat_side(library, "reads.fq", "out.bam", 1, tmp_path) == "out.bam"
