"""``te-split`` / ``te-search-chunk`` / ``te-merge-trim``: the BLAT TE search across jobs.

The three subcommands run the steps ``relocaTE3 run`` runs in one process
(:mod:`RelocaTE3.te_scatter`), so the Nextflow pipeline can search the chunks
of one sample in many jobs. Run one after another, in separate directories,
they must write exactly the files ``relocaTE3 run`` writes.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pytest

from RelocaTE3 import cli, te_scatter
from RelocaTE3.aligners import TEChunk

from tests.te_scatter_test import R1, R2, RICETE, needs_blat, outputs

SCATTER = ("te-split", "te-search-chunk", "te-merge-trim")


def _subparser(name: str) -> argparse.ArgumentParser:
    for action in cli.build_parser()._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action.choices[name]
    raise AssertionError("build_parser() has no subcommands")


def _actions(name: str) -> dict[str, argparse.Action]:
    return {a.dest: a for a in _subparser(name)._actions}


def _parse(*argv: str) -> argparse.Namespace:
    return cli.build_parser().parse_args(list(argv))


# ---------------------------------------------------------------------------
# parsers
# ---------------------------------------------------------------------------


def test_te_split_parser():
    args = _parse("te-split", "-1", "a.fq", "-2", "b.fq", "-n", "S", "-o", "out")
    assert (args.left, args.right, args.name, args.outdir) == ("a.fq", "b.fq", "S", "out")
    assert args.chunk_reads == 200_000
    assert args.func is cli.cmd_te_split
    # run's spelling works too
    args = _parse("te-split", "--left", "a.fq", "-n", "S", "--chunk-reads", "5000")
    assert (args.left, args.right, args.outdir, args.chunk_reads) == ("a.fq", None, ".", 5000)


def test_te_search_chunk_parser():
    args = _parse(
        "te-search-chunk", "-T", "te.fa", "--manifest", "m.tsv",
        "--chunk", "c0.fa", "c1.fa", "--chunk", "c2.fa", "--te-opts=-minIdentity=80",
        "--threads", "4",
    )
    assert args.chunks == ["c0.fa", "c1.fa", "c2.fa"]
    assert (args.te_library, args.manifest, args.te_opts, args.threads, args.outdir) == (
        "te.fa", "m.tsv", "-minIdentity=80", 4, "."
    )
    assert args.func is cli.cmd_te_search_chunk


def test_te_merge_trim_parser():
    args = _parse(
        "te-merge-trim", "-1", "a.fq", "-2", "b.fq", "-T", "te.fa", "-n", "S",
        "--manifest", "m.tsv", "--chunk-bams", "x.bam", "y.bam",
    )
    assert args.chunk_bams == ["x.bam", "y.bam"]
    assert (args.minimum_match_length, args.minimum_trimmed_length, args.mismatch_allowance) == (
        10, 10, 2
    )
    assert args.func is cli.cmd_te_merge_trim


@pytest.mark.parametrize("missing", ["--manifest", "--chunk"])
def test_te_search_chunk_requires(missing):
    argv = ["te-search-chunk", "-T", "te.fa", "--manifest", "m.tsv", "--chunk", "c.fa"]
    i = argv.index(missing)
    del argv[i : i + 2]
    with pytest.raises(SystemExit):
        _parse(*argv)


# ---------------------------------------------------------------------------
# guard: every `run` option that shapes the BLAT search or the trim must be
# accepted, with the same spelling, by the scatter subcommand that does it
# ---------------------------------------------------------------------------

#: Every `relocaTE3 run` dest -> scatter subcommands that must accept it.
#: A new `run` option fails test_every_run_option_is_classified until it is
#: added here (and, if it shapes the search or trim, to the subcommands).
RUN_DEST_TO_SCATTER = {
    "left": ("te-split", "te-merge-trim"),
    "right": ("te-split", "te-merge-trim"),
    "name": ("te-split", "te-merge-trim"),
    "outdir": SCATTER,
    "threads": ("te-search-chunk", "te-merge-trim"),
    "te_library": ("te-search-chunk", "te-merge-trim"),
    "te_opts": ("te-search-chunk",),  # BLAT command line
    "minimum_match_length": ("te-merge-trim",),
    "minimum_trimmed_length": ("te-merge-trim",),
    "mismatch_allowance": ("te-merge-trim",),
    "verbose": SCATTER,
    # Scatter is the BLAT path only (the pipeline uses it when te_aligner is blat).
    "te_aligner": (),
    # argparse internals
    "help": SCATTER,
    "version": SCATTER,
}


def test_every_run_option_is_classified():
    run_dests = set(_actions("run")) - {"func"}
    assert run_dests == set(RUN_DEST_TO_SCATTER), (
        "relocaTE3 run's options changed: classify the new ones in "
        "RUN_DEST_TO_SCATTER (does the option change the BLAT search or the trim?)"
    )


@pytest.mark.parametrize(
    "dest, sub",
    [(d, s) for d, subs in sorted(RUN_DEST_TO_SCATTER.items()) for s in subs],
)
def test_scatter_subcommand_accepts_run_option(dest, sub):
    run = _actions("run")[dest]
    actions = _actions(sub)
    assert dest in actions, f"{sub} lacks run's {'/'.join(run.option_strings)}"
    mine = actions[dest]
    # run's spellings all work (the pipeline can pass the same flags) ...
    assert set(run.option_strings) <= set(mine.option_strings)
    # ... with the same meaning
    assert (mine.type, mine.default) == (run.type, run.default), dest


# ---------------------------------------------------------------------------
# behaviour
# ---------------------------------------------------------------------------


def test_te_search_chunk_takes_mate_from_manifest_and_passes_te_opts(tmp_path, monkeypatch):
    rows = [
        te_scatter.ManifestRow("left", 0, "chunk.S.left.00000.fa", 3, "1", "a.fq"),
        te_scatter.ManifestRow("right", 0, "chunk.S.right.00000.fa", 3, "2", "b.fq"),
    ]
    manifest = te_scatter.write_manifest(tmp_path / "manifest.tsv", rows)
    calls = []

    def fake_search(te_library, chunk_fa, chunk_bam, mate=None, te_opts=(), **kw):
        calls.append((Path(chunk_fa).name, Path(chunk_bam), mate, list(te_opts)))
        return TEChunk(bam=Path(chunk_bam), n_reads=3)

    monkeypatch.setattr(te_scatter, "search_te_chunk", fake_search)
    monkeypatch.setattr(te_scatter, "te_sam_header", lambda lib: "@HD\n")
    out = tmp_path / "out"
    rc = cli.main([
        "te-search-chunk", "-T", "te.fa", "--manifest", str(manifest), "-o", str(out),
        "--chunk", "split/chunk.S.right.00000.fa", "split/chunk.S.left.00000.fa",
        "--te-opts", "-minIdentity=80 -fine",
    ])
    assert rc == 0
    assert sorted(calls) == [
        ("chunk.S.left.00000.fa", out / "chunk.S.left.00000.bam", "1",
         ["-minIdentity=80", "-fine"]),
        ("chunk.S.right.00000.fa", out / "chunk.S.right.00000.bam", "2",
         ["-minIdentity=80", "-fine"]),
    ]


def test_te_search_chunk_refuses_chunk_not_in_manifest(tmp_path):
    rows = [te_scatter.ManifestRow("left", 0, "chunk.S.left.00000.fa", 3, None, "a.fq")]
    manifest = te_scatter.write_manifest(tmp_path / "manifest.tsv", rows)
    args = _parse("te-search-chunk", "-T", "te.fa", "--manifest", str(manifest),
                  "--chunk", "chunk.S.left.00001.fa")
    with pytest.raises(ValueError, match="not in the manifest.*chunk.S.left.00001.fa"):
        cli.cmd_te_search_chunk(args)


def test_te_split_refuses_existing_chunk_dir(tmp_path):
    fq = tmp_path / "a.fq"
    fq.write_text("@r/1\nACGT\n+\nIIII\n")
    (tmp_path / "out" / "S.te_chunks").mkdir(parents=True)
    args = _parse("te-split", "-1", str(fq), "-n", "S", "-o", str(tmp_path / "out"))
    with pytest.raises(FileExistsError, match="S.te_chunks"):
        cli.cmd_te_split(args)


@pytest.fixture(scope="module")
def run_outdir(tmp_path_factory) -> Path:
    """``relocaTE3 run`` on the Chr3 fixture (BLAT, defaults)."""
    out = tmp_path_factory.mktemp("run")
    assert cli.main([
        "run", "-l", str(R1), "-r", str(R2), "-T", str(RICETE), "-n", "S",
        "-o", str(out), "--threads", "4",
    ]) == 0
    return out


@pytest.fixture(scope="module")
def scattered(tmp_path_factory):
    """te-split at 5000 reads per chunk; te-search-chunk as two separate 'jobs'."""
    root = tmp_path_factory.mktemp("scatter")
    assert cli.main([
        "te-split", "-1", str(R1), "-2", str(R2), "-n", "S", "-o", str(root / "split"),
        "--chunk-reads", "5000",
    ]) == 0
    chunk_dir = root / "split" / "S.te_chunks"
    manifest = chunk_dir / "manifest.tsv"
    fastas = sorted(chunk_dir.glob("chunk.*.fa"))
    jobs = [fastas[0::2], fastas[1::2]]  # each job gets chunks of both sides
    for number, group in enumerate(jobs):
        assert cli.main([
            "te-search-chunk", "-T", str(RICETE), "--manifest", str(manifest),
            "-o", str(root / f"search{number}"), "--threads", "3",
            "--chunk", *map(str, group),
        ]) == 0
    bams = sorted((root / "search0").glob("*.bam")) + sorted((root / "search1").glob("*.bam"))
    return root, manifest, fastas, bams


@needs_blat
def test_te_split_writes_chunks_and_manifest(scattered):
    root, manifest, fastas, _ = scattered
    rows = te_scatter.read_manifest(manifest)
    assert [(r.side, r.index) for r in rows] == (
        [("left", i) for i in range(16)] + [("right", i) for i in range(16)]
    )
    assert {r.mate for r in rows if r.side == "left"} == {"1"}
    assert sum(r.n_reads for r in rows if r.side == "right") == 75_910
    assert sorted(r.chunk for r in rows) == [f.name for f in fastas]
    assert sorted(p.name for p in (root / "split").iterdir()) == ["S.te_chunks"]


@needs_blat
def test_te_search_chunk_writes_one_indexed_bam_per_chunk(scattered):
    root, _, fastas, bams = scattered
    assert sorted(b.name for b in bams) == sorted(te_scatter.chunk_bam_name(f) for f in fastas)
    for d in ("search0", "search1"):
        names = sorted(p.name for p in (root / d).iterdir())
        assert all(n.endswith((".bam", ".bam.bai")) for n in names), names
        assert len(names) == 2 * len(list((root / d).glob("*.bam")))


@needs_blat
def test_scatter_subcommands_match_run(run_outdir, scattered, tmp_path):
    _, manifest, _, bams = scattered
    out = tmp_path / "merged"
    assert cli.main([
        "te-merge-trim", "-1", str(R1), "-2", str(R2), "-T", str(RICETE), "-n", "S",
        "-o", str(out), "--manifest", str(manifest), "--threads", "2",
        "--chunk-bams", *map(str, bams),
    ]) == 0
    # exactly run's files, .bai included; no chunk or scratch leftovers
    files = lambda d: sorted(str(p.relative_to(d)) for p in d.rglob("*") if p.is_file())
    assert files(out) == files(run_outdir)
    assert len(files(out)) == 12
    want, got = outputs(run_outdir), outputs(out)
    for rel in want:
        assert got[rel] == want[rel], f"{rel} differs"


@needs_blat
def test_te_merge_trim_refuses_missing_chunk_bam(scattered, tmp_path):
    _, manifest, _, bams = scattered
    dropped = bams[3]
    args = _parse(
        "te-merge-trim", "-1", str(R1), "-2", str(R2), "-T", str(RICETE), "-n", "S",
        "-o", str(tmp_path / "merged"), "--manifest", str(manifest),
        "--chunk-bams", *[str(b) for b in bams if b != dropped],
    )
    with pytest.raises(ValueError, match=f"missing chunk BAM.*{dropped.name}"):
        cli.cmd_te_merge_trim(args)
    assert not (tmp_path / "merged").exists() or not any((tmp_path / "merged").iterdir())
    assert cli.main([
        "te-merge-trim", "-1", str(R1), "-2", str(R2), "-T", str(RICETE), "-n", "S",
        "-o", str(tmp_path / "merged2"), "--manifest", str(manifest),
        "--chunk-bams", *[str(b) for b in bams if b != dropped],
    ]) == 1
