"""Controls for the experimental intact-read veto; production is unchanged."""

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest
from RelocaTE3 import insertions as ins
from RelocaTE3.models import Insertion

ROOT = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location(
    "intact_replay", ROOT / "scripts/replay_intact_fullreads.py"
)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


@pytest.mark.parametrize(
    "cigar,expected",
    [
        ([(0, 150)], True),
        ([(4, 10), (0, 140)], True),
        ([(4, 11), (0, 139)], False),
        ([(0, 130), (1, 10), (4, 10)], True),
        ([(7, 135), (8, 5), (4, 10)], True),
        ([(0, 130), (2, 10), (4, 20)], False),
    ],
)
def test_query_coverage_not_reference_span_or_proper_pair(cigar, expected):
    record = SimpleNamespace(
        is_unmapped=False, query_sequence="A" * 150, cigartuples=cigar
    )
    assert audit.is_intact(record) == expected
    assert ins._is_intact_fullread(record) == expected
    record.is_unmapped = True
    assert not audit.is_intact(record)


def test_mate_identity_and_empty_sides():
    call = Insertion(
        chrom="Chr1",
        start=100,
        end=102,
        te_name="TE",
        strand="+",
        tsd="UNK",
        left_junction_reads=1,
        read_names=["r/1:end:5"],
    )

    def check(keys):
        return audit.veto(call, keys, ins._junction_fullread_key, ins._fullread_key)[0]

    assert not check({"r/2"})
    assert check({"r/1"})
    assert check({"r"})  # historical unpaired sidecar
    call.left_junction_reads = 0
    call.read_names = []
    assert not check({"r/1"})


def test_real_displaced_fullread_explains_reference_sine_call():
    fixture = json.loads((ROOT / "tests/data/intact_fullreads/Os3912.json").read_text())

    def records(label):
        data = fixture["alignments"][label]
        header = pysam.AlignmentHeader.from_dict(data["header"])
        return [pysam.AlignedSegment.fromstring(s, header) for s in data["sam"]]

    junctions = records("R3_junctions")
    left, right = [], []
    for r in junctions:
        info = ins._junction_info(
            r.query_name,
            "-" if r.is_reverse else "+",
            r.reference_start + 1,
            r.reference_end,
        )
        (left if info[0] == "left" else right).append(r.query_name)
    assert (len(left), len(right)) == (19, 1)
    call = Insertion(
        chrom="Chr4",
        start=19264325,
        end=19264335,
        te_name="Os3912",
        strand="-",
        tsd="CTTAGTCGATA",
        left_junction_reads=19,
        right_junction_reads=1,
        read_names=left + right,
    )
    full = records("R3_fullreads")

    class Bam:
        def fetch(self, *args):
            return iter(full)

    assert ins._fullread_false_junction(Bam(), call)
    assert audit.veto(
        call,
        audit.local_intact(Bam(), call, ins),
        ins._junction_fullread_key,
        ins._fullread_key,
    ) == (True, 19, 1)
    intact = {ins._fullread_record_key(r) for r in full if audit.is_intact(r)}
    assert audit.veto(call, intact, ins._junction_fullread_key, ins._fullread_key) == (
        True,
        19,
        1,
    )
    displaced = next(
        r
        for r in full
        if ins._fullread_record_key(r) == ins._junction_fullread_key(right[0])
    )
    assert (
        displaced.reference_start + 1,
        displaced.reference_end,
        displaced.cigarstring,
    ) == (19264493, 19264642, "150M")
    # Same read end is intact in the archived R2 full-read alignment.
    assert any(
        ins._fullread_record_key(r) == ins._fullread_record_key(displaced)
        and audit.is_intact(r)
        for r in records("R2_fullreads")
    )


def test_local_policy_excludes_remote_and_other_contig_alignments(tmp_path):
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"SO": "coordinate"},
            "SQ": [{"SN": "Chr1", "LN": 20000}, {"SN": "Chr2", "LN": 20000}],
        }
    )
    path = tmp_path / "full.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, chrom, start, cigar in (
            ("near", 0, 1100, "150M"),
            ("clipped", 0, 1100, "100M50S"),
            ("remote", 0, 10000, "150M"),
            ("other", 1, 1100, "150M"),
        ):
            record = pysam.AlignedSegment(header)
            record.query_name = name
            record.query_sequence = "A" * 150
            record.reference_id = chrom
            record.reference_start = start
            record.flag = 65
            record.cigarstring = cigar
            out.write(record)
    pysam.index(str(path))
    call = Insertion(
        chrom="Chr1", start=1000, end=1002, te_name="TE", strand="+", tsd="UNK"
    )
    with pysam.AlignmentFile(str(path), "rb") as bam:
        assert audit.local_intact(bam, call, ins) == {"near/1"}
        call.left_junction_reads = 1
        for name, expected in (
            ("near/1", True),
            ("clipped/1", False),
            ("remote/1", False),
            ("other/1", False),
            ("near/2", False),
        ):
            call.read_names = [name + ":end:5"]
            assert ins._fullread_false_junction(bam, call) == expected
    assert audit.local_intact(None, call, ins) == set()


def test_spanning_and_intact_thresholds_are_not_mixed_across_sides():
    call = Insertion(
        chrom="Chr1",
        start=1000,
        end=1002,
        te_name="TE",
        strand="+",
        tsd="UNK",
        left_junction_reads=1,
        right_junction_reads=1,
        read_names=["left:end:5", "right:start:3"],
    )

    def record(name, start, cigar):
        r = pysam.AlignedSegment()
        r.query_name = name
        r.query_sequence = "A" * 150
        r.reference_start = start
        r.cigarstring = cigar
        return r

    # Left spans but is clipped; right is intact nearby but does not span.
    records = [record("left", 950, "100M50S"), record("right", 1200, "150M")]

    class Bam:
        def fetch(self, *args):
            return iter(records)

    assert not ins._fullread_false_junction(Bam(), call)
