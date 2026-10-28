"""Step 4 loads only the TE-hit names it can ever query.

``te_hit_names.txt`` lists every read with a selected TE hit: ~34M names on
riceTElib 30x with BLAT (1.7 GB), several GiB as a Python set. The step-4
planner only asks about flanking reads and their two mates, so restricting the
set to names whose pair base appears in the flanking FASTQs must give identical
plans and FASTQs. RelocaTE2 likewise only resolves mates of the reads it trimmed
(clean_pairs_memory.py).
"""

from __future__ import annotations

from pathlib import Path

import pytest

import RelocaTE3.genome_align as genome_align
from RelocaTE3.genome_align import (
    flanking_read_bases,
    plan_genome_alignment_inputs,
    read_te_hit_names,
)
from RelocaTE3.ReadLibrary import ReadLibrary


def _write_fastq(path: Path, records: list[tuple[str, str]]) -> Path:
    with open(path, "w") as out:
        for name, sequence in records:
            out.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")
    return path


BASES = ("jj", "jm", "ju", "jn", "mm", "mu", "mn")
TE_HITS = [
    "jj/1:end:5", "jj/2:start:3", "jm/1:end:5", "jm/2:middle", "ju/1:end:5",
    "ju/2", "jn/2:start:3", "mm/1:middle", "mm/2:middle", "mu/1:middle", "mu/2",
    "mn/1:middle",
    # TE hits of reads with no flank on either end: never queried by step 4.
    "noise/1", "noise/2:middle", "other/2", "jjx/1", "single",
]


def _fixture(tmp_path: Path):
    r1 = _write_fastq(
        tmp_path / "reads_R1.fq",
        [(f"{b}/1", "ACAC" if b == "jn" else "AAAA") for b in BASES + ("noise",)],
    )
    r2 = _write_fastq(
        tmp_path / "reads_R2.fq",
        [(f"{b}/2", "GTGT" if b == "mn" else "TTTT") for b in BASES + ("noise",)],
    )
    left = _write_fastq(
        tmp_path / "S.left.flankingReads.fq",
        [
            ("jj/1:end:5", "JJJ1"),
            ("jm/1:end:5", "JMJ1"),
            ("ju/1:end:5", "JUJ1"),
            ("mm/1:middle", "MMM1"),
            ("mu/1:middle", "MUM1"),
            ("mn/1:middle", "MNM1"),
        ],
    )
    right = _write_fastq(
        tmp_path / "S.right.flankingReads.fq",
        [
            ("jj/2:start:3", "JJJ2"),
            ("jm/2:middle", "JMM2"),
            ("jn/2:start:3", "JNJ2"),
            ("mm/2:middle", "MMM2"),
        ],
    )
    names = tmp_path / "S.te_hit_names.txt"
    names.write_text("".join(f"{n}\n" for n in TE_HITS))
    return ReadLibrary([str(r1), str(r2)], "S"), [str(left), str(right)], names


def test_flanking_bases_and_filtered_loader(tmp_path):
    _, flanks, names = _fixture(tmp_path)
    bases = flanking_read_bases(flanks)
    assert bases == {"jj", "jm", "ju", "jn", "mm", "mu", "mn"}
    assert read_te_hit_names(names, bases=bases) == read_te_hit_names(names) - {
        "noise/1", "noise/2", "other/2", "jjx/1", "single"
    }


def test_single_end_names_without_mate_suffix_are_kept(tmp_path):
    flank = _write_fastq(tmp_path / "S.left.flankingReads.fq", [("single:end:5", "AA")])
    names = tmp_path / "S.te_hit_names.txt"
    names.write_text("single\nsingle/1\nother\n")
    bases = flanking_read_bases([str(flank)])
    assert read_te_hit_names(names, bases=bases) == {"single", "single/1"}


def test_planner_output_identical_with_filtered_te_hits(tmp_path):
    reads, flanks, names = _fixture(tmp_path)

    def plan(all_te_hits, tag):
        outs = [tmp_path / f"{tag}.{n}.fq" for n in ("p1", "p2", "junction", "support")]
        result = plan_genome_alignment_inputs(flanks, all_te_hits, reads, *outs)
        return result, [p.read_bytes() for p in outs]

    full = plan(read_te_hit_names(names), "full")
    filtered = plan(read_te_hit_names(names, bases=flanking_read_bases(flanks)), "filt")
    assert filtered == full
    assert (full[0].paired, full[0].single_junctions, full[0].support_reads) == (2, 2, 1)


def test_align_to_genome_hands_planner_the_filtered_set(tmp_path, monkeypatch):
    reads, flanks, names = _fixture(tmp_path)
    outdir = tmp_path / "out"
    (outdir / "flanking").mkdir(parents=True)
    (outdir / "te_containing").mkdir()
    for flank in flanks:
        Path(flank).rename(outdir / "flanking" / Path(flank).name)
    names.rename(outdir / "te_containing" / names.name)

    seen = {}

    def capture(backend, genome, flanking_files, all_te_hits, *args, **kwargs):
        seen["hits"] = set(all_te_hits)
        raise RuntimeError("stop after planning inputs")

    monkeypatch.setattr(genome_align, "align_flanks_anchored", capture)
    with pytest.raises(RuntimeError, match="stop after"):
        genome_align.align_to_genome(reads, "genome.fa", outdir, genome_aligner="bwaaln")
    assert "jn/2" in seen["hits"]
    assert not {"noise/1", "other/2", "single"} & seen["hits"]
