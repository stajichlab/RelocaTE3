"""Exact SAM output and incremental consumption for the production converter."""

import pytest

from RelocaTE3.aligners import _iter_psl_to_sam, psl_to_sam


def psl(**changes):
    fields = [
        "4",
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        "+",
        "read",
        "8",
        "2",
        "6",
        "TE",
        "500",
        "100",
        "104",
        "1",
        "4,",
        "2,",
        "100,",
    ]
    for index, value in changes.items():
        fields[int(index[1:])] = str(value)
    return "\t".join(fields)


@pytest.mark.parametrize(
    "line,sequences,expected",
    [
        (
            psl(),
            {"read": "ACGTACGT"},
            "read\t0\tTE\t101\t255\t2S4M2S\t*\t0\t0\tACGTACGT\t*\tNM:i:0",
        ),
        (
            psl(c8="-", c11=1, c12=5, c19="1,", c1=1),
            {"read": "AACCGGTT"},
            "read\t16\tTE\t101\t255\t3S4M1S\t*\t0\t0\tAACCGGTT\t*\tNM:i:1",
        ),
        (
            psl(c8="-"),
            {"read": "AAAACCGT"},
            "read\t16\tTE\t101\t255\t2S4M2S\t*\t0\t0\tACGGTTTT\t*\tNM:i:0",
        ),
        (
            psl(
                c4=1,
                c5=1,
                c6=1,
                c7=2,
                c11=1,
                c12=6,
                c16=106,
                c17=2,
                c18="2,2,",
                c19="1,4,",
                c20="100,104,",
                c1=1,
            ),
            None,
            "read\t0\tTE\t101\t255\t1S2M1I2D2M2S\t*\t0\t0\t*\t*\tNM:i:4",
        ),
        (psl(), {}, "read\t0\tTE\t101\t255\t2S4M2S\t*\t0\t0\t*\t*\tNM:i:0"),
    ],
)
def test_exact_records_and_list_compatibility(line, sequences, expected):
    assert list(_iter_psl_to_sam([line], sequences)) == [expected]
    result = psl_to_sam([line], sequences)
    assert isinstance(result, list)
    assert result == [expected]


def test_iterator_is_lazy_and_writes_input_order_without_deduplication():
    consumed = []

    def source():
        for name in ("b", "a", "b"):
            consumed.append(name)
            yield psl(c9=name)

    result = _iter_psl_to_sam(source())
    assert consumed == []
    assert next(result).startswith("b\t")
    assert consumed == ["b"]
    assert [line.split("\t")[0] for line in result] == ["a", "b"]


@pytest.mark.parametrize(
    "changes", [{"c4": 2}, {"c5": 4}, {"c6": 2}, {"c7": 4}, {"c17": 3}]
)
def test_rejected_alignments_headers_and_empty_input(changes):
    lines = ["psLayout version 3", "", "match\tmis", psl(**changes), psl()]
    assert list(_iter_psl_to_sam(lines)) == psl_to_sam([psl()])
    assert list(_iter_psl_to_sam([])) == []
    assert psl_to_sam([]) == []


def test_late_parse_failure_is_not_silently_ignored():
    lines = [psl(), psl(c1="invalid")]
    iterator = _iter_psl_to_sam(lines)
    assert next(iterator) == psl_to_sam([psl()])[0]
    with pytest.raises(ValueError):
        next(iterator)
    with pytest.raises(ValueError):
        psl_to_sam(lines)
