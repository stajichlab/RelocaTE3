"""align-reads builds the genotyping alignment, and run-all delegates to it.

The genotyping BAM (all original reads -> genome) used to be built inline in
run-all's characterize step. In the Nextflow pipeline it runs as its own
process, in parallel with the TE search, so it has to be a subcommand that
writes the same file run-all always wrote.
"""

from __future__ import annotations

import random
from pathlib import Path

import pysam
import pytest

from RelocaTE3 import cli
from RelocaTE3.cli import build_parser, main

from tests.main_test import _fake_stage_outputs

_COMP = str.maketrans("ACGT", "TGCA")


def _genome(tmp_path) -> Path:
    """A small random two-contig genome FASTA."""
    rng = random.Random(1)
    fa = tmp_path / "g.fa"
    fa.write_text(
        "".join(
            f">c{i}\n{''.join(rng.choice('ACGT') for _ in range(3000))}\n"
            for i in range(2)
        )
    )
    return fa


def _simulate_pairs(fa: Path, tmp_path, n: int = 20) -> tuple[Path, Path]:
    """n read pairs: 100 bp reads, 300 bp insert, R2 reverse-complemented."""
    seqs = {}
    name = None
    for line in fa.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        else:
            seqs[name] += line
    rng = random.Random(2)
    r1, r2 = [], []
    contigs = sorted(seqs)
    for i in range(n):
        chrom = seqs[contigs[i % len(contigs)]]
        start = rng.randrange(0, len(chrom) - 300)
        frag = chrom[start : start + 300]
        q = "I" * 100
        r1.append(f"@p{i}/1\n{frag[:100]}\n+\n{q}\n")
        r2.append(f"@p{i}/2\n{frag[-100:].translate(_COMP)[::-1]}\n+\n{q}\n")
    f1, f2 = tmp_path / "r_1.fq", tmp_path / "r_2.fq"
    f1.write_text("".join(r1))
    f2.write_text("".join(r2))
    return f1, f2


def test_align_reads_parser_required_args_and_threads_default():
    ns = build_parser().parse_args(
        ["align-reads", "-1", "a.fq", "-g", "g.fa", "-n", "S", "-o", "out"]
    )
    assert ns.left == "a.fq"
    assert ns.right is None
    assert ns.genome_fasta == "g.fa"
    assert ns.name == "S"
    assert ns.outdir == "out"
    assert ns.threads == 1
    assert ns.func is cli.cmd_align_reads


@pytest.mark.parametrize("missing", ["-1", "-g", "-n"])
def test_align_reads_parser_requires(missing):
    argv = {"-1": "a.fq", "-g": "g.fa", "-n": "S"}
    del argv[missing]
    flat = [x for kv in argv.items() for x in kv]
    with pytest.raises(SystemExit):
        build_parser().parse_args(["align-reads", *flat])


def test_align_reads_writes_sorted_indexed_bam(tmp_path):
    fa = _genome(tmp_path)
    r1, r2 = _simulate_pairs(fa, tmp_path)
    out = tmp_path / "out"
    rc = main(
        [
            "align-reads",
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-g",
            str(fa),
            "-n",
            "S",
            "-o",
            str(out),
            "--threads",
            "2",
        ]
    )
    assert rc == 0
    bam = out / "genome_aln" / "S.reads.genome.bam"
    assert bam.is_file()
    assert Path(f"{bam}.bai").is_file()
    with pysam.AlignmentFile(str(bam)) as fh:
        assert fh.header.to_dict()["HD"]["SO"] == "coordinate"
        mapped = sum(1 for r in fh.fetch() if not r.is_unmapped)
    assert mapped == 40


def _record_run_all(monkeypatch):
    calls: list[tuple[str, object]] = []

    def rec(tag):
        def _inner(args):
            calls.append((tag, args))
            return 0

        return _inner

    for name, tag in [
        ("cmd_index_genome", "index-genome"),
        ("cmd_run", "run"),
        ("cmd_align_genome", "align-genome"),
        ("cmd_find_insertions", "find-insertions"),
        ("cmd_find_reference", "find-reference"),
        ("cmd_align_reads", "align-reads"),
        ("cmd_characterize", "characterize"),
    ]:
        monkeypatch.setattr(cli, name, rec(tag))
    return calls


def _run_all_argv(tmp_path, *extra):
    return [
        "run-all",
        "-1",
        "r1.fq",
        "-2",
        "r2.fq",
        "-T",
        "te.fa",
        "-g",
        "ref.fa",
        "-n",
        "HEG4",
        "-o",
        str(tmp_path),
        "--threads",
        "3",
        "--genotype",
        *extra,
    ]


def _fake_sites(tmp_path):
    _fake_stage_outputs(tmp_path, "HEG4")
    (tmp_path / "results" / "ALL.repeat.all_nonref_insert.txt").write_text("x\n")


def test_run_all_genotype_delegates_to_align_reads(monkeypatch, tmp_path):
    calls = _record_run_all(monkeypatch)
    _fake_sites(tmp_path)

    assert main(_run_all_argv(tmp_path)) == 0
    tags = [t for t, _ in calls]
    assert "align-reads" in tags
    assert tags.index("align-reads") < tags.index("characterize")
    assert tags.index("align-reads") > tags.index("find-insertions")

    ar = dict(calls)["align-reads"]
    assert ar.left == "r1.fq" and ar.right == "r2.fq"
    assert ar.genome_fasta == "ref.fa"
    assert ar.name == "HEG4"
    assert ar.outdir == str(tmp_path)
    assert ar.threads == 3

    ch = dict(calls)["characterize"]
    assert ch.bam == [str(tmp_path / "genome_aln" / "HEG4.reads.genome.bam")]


def test_run_all_genotype_bam_skips_align_reads(monkeypatch, tmp_path):
    calls = _record_run_all(monkeypatch)
    _fake_sites(tmp_path)

    assert main(_run_all_argv(tmp_path, "--genotype-bam", "x.cram")) == 0
    tags = [t for t, _ in calls]
    assert "align-reads" not in tags
    assert dict(calls)["characterize"].bam == ["x.cram"]
