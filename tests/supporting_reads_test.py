"""Support-only insertions: RelocaTE2's ``all_nonref_supporting`` output.

RelocaTE2 calls insertions from paired-end mates alone when no junction read
mapped, and writes them to a *separate* file
(``<target>.<TE>.all_nonref_supporting.{txt,gff}``, relocaTE_insertionFinder.py:
151, 431-459). They never enter ``all_nonref_insert`` -- they carry T:0 R:0 L:0
and are much weaker evidence, so RelocaTE2 segregates them.

RelocaTE3 had ``_call_support_only`` but it was reachable only from the
module-level ``find_insertions`` used by ``pipeline.run_sample`` -- *not* from
``InsertionFinder``, which is what the CLI and the benchmark run. So the shipped
pipeline emitted no support-only calls at all, and no such file.

RelocaTE2's three cases (relocaTE_insertionFinder.py:431-459):

    both strands   ins_start = boundary(+ reads), ins_end = boundary(- reads)
    + strand only  ins_end   = ins_start + lib_size * 1.2
    - strand only  ins_start = ins_end   - lib_size * 1.2

where lib_size is ``--size`` (default 500), and the 1.2 is "insert size * (1 +
sd of library)" in RelocaTE2's own comment.
"""

from __future__ import annotations

from RelocaTE3.insertions import _Cluster, call_support_only

INSERT = 500
SPAN = int(INSERT * 1.2)  # 600


def _cluster(reads: list[tuple[int, int, str]], chrom: str = "Chr3") -> _Cluster:
    """Cluster carrying only supporting reads: (gstart, gend, strand)."""
    cluster = _Cluster(chrom)
    for index, (gstart, gend, strand) in enumerate(reads):
        cluster.support.append((f"read{index}", gstart, gend, strand, "ACGT"))
        cluster.extend(gstart, gend)
    return cluster


def test_both_strands_bracket_the_insertion():
    """The site lies in the gap between the innermost mates."""
    ins = call_support_only(
        _cluster([(100, 200, "+"), (120, 220, "+"), (900, 1000, "-")]),
        insert_size=INSERT,
    )
    assert ins is not None
    assert (ins.start, ins.end) == (220, 900), "rightmost + end .. leftmost - start"
    assert ins.left_support_reads == 2
    assert ins.right_support_reads == 1
    assert ins.tsd == "supporting_reads"
    assert ins.left_junction_reads == 0 and ins.right_junction_reads == 0


def test_plus_strand_only_extends_by_the_library_insert_size():
    """relocaTE_insertionFinder.py:446 -- ins_end = ins_start + lib_size * 1.2."""
    ins = call_support_only(_cluster([(100, 200, "+")]), insert_size=INSERT)
    assert ins is not None
    assert (ins.start, ins.end) == (200, 200 + SPAN)
    assert ins.left_support_reads == 1 and ins.right_support_reads == 0


def test_minus_strand_only_extends_backwards():
    """relocaTE_insertionFinder.py:455 -- ins_start = ins_end - lib_size * 1.2."""
    ins = call_support_only(_cluster([(900, 1000, "-")]), insert_size=INSERT)
    assert ins is not None
    assert (ins.start, ins.end) == (900 - SPAN, 900)
    assert ins.right_support_reads == 1 and ins.left_support_reads == 0


def test_insert_size_is_honoured():
    """--size is RelocaTE2's -s/--size, default 500."""
    ins = call_support_only(_cluster([(100, 200, "+")]), insert_size=1000)
    assert ins is not None
    assert ins.end == 200 + int(1000 * 1.2)


def test_overlapping_brackets_are_rejected():
    """relocaTE_insertionFinder.py:440 -- `if ins_start > ins_end: continue`.

    Mates that overlap cannot bracket a gap, so the site is ambiguous.
    """
    assert call_support_only(
        _cluster([(900, 1000, "+"), (100, 200, "-")]), insert_size=INSERT
    ) is None


def test_cluster_without_supporting_reads_yields_nothing():
    assert call_support_only(_cluster([]), insert_size=INSERT) is None


def test_negative_coordinates_are_clamped():
    """A - only cluster near the contig start must not produce a negative start."""
    ins = call_support_only(_cluster([(10, 50, "-")]), insert_size=INSERT)
    assert ins is not None
    assert ins.start >= 1, "GFF coordinates are 1-based and must stay positive"


def test_support_only_call_retains_family_evidence_without_assigning_primary():
    cluster = _cluster([(100, 200, "+"), (900, 1000, "-")])
    ins = call_support_only(
        cluster,
        insert_size=INSERT,
        read_repeat={
            "read0/1": ("mPing", "+"),
            "read1/2": ("mPing", "+"),
        },
    )

    assert ins is not None
    assert ins.te_name == "NA", "indirect mate evidence must not become primary"
    assert ins.te_supporting_family_support == {"mPing": 2}
    assert ins.te_supporting_family_confidence == 1.0
    assert ins.te_supporting_family_status == "unique"
    assert ins.te_family_concordance == "supporting_only"
