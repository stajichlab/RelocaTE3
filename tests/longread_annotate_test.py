"""Long-read TE annotation: signature sequences aligned to the TE library.

The library is built from random sequence with a fixed seed, so the expected
hits follow from how each query was made (copy, reverse complement, nesting),
not from the code under test. Random "TEs" are not real TEs; they only need
to be distinct from each other and from random background.
"""

from __future__ import annotations

import random
import shutil
from pathlib import Path

import pytest

from RelocaTE3.longread.annotate import (
    annotate_signatures, parse_paf_line, te_minimap2_command,
)  # fmt: skip
from RelocaTE3.longread.model import (
    NO_HIT, ClipSide, Locus, Signature, SignatureKind, TeHit,
)  # fmt: skip

needs_minimap2 = pytest.mark.skipif(
    shutil.which("minimap2") is None, reason="minimap2 not on PATH"
)

BASES = "ACGT"
_COMP = str.maketrans("ACGT", "TGCA")


def _rand(n: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _mutate(seq: str, rate: float, rng: random.Random) -> str:
    """Apply substitutions, insertions and deletions, each at ``rate / 3``."""
    out = []
    for base in seq:
        r = rng.random()
        if r < rate / 3:
            out.append(rng.choice([b for b in BASES if b != base]))
        elif r < 2 * rate / 3:
            out.append(base + rng.choice(BASES))
        elif r < rate:
            continue
        else:
            out.append(base)
    return "".join(out)


@pytest.fixture(scope="module")
def tes() -> dict[str, str]:
    rng = random.Random(7)
    return {"mGing": _rand(145, rng), "mPing": _rand(430, rng),
            "Retro1": _rand(1079, rng)}  # fmt: skip


@pytest.fixture
def library(tmp_path: Path, tes: dict[str, str]) -> Path:
    path = tmp_path / "lib.fa"
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in tes.items()))
    return path


def _sig(
    sig_id: str,
    read_id: str,
    seq: str,
    kind: SignatureKind = SignatureKind.INS,
    side: ClipSide = ClipSide.NONE,
) -> Signature:
    return Signature(
        signature_id=sig_id,
        kind=kind,
        locus=Locus("ref", "Chr3", 1000, 1000, "+"),
        insert_length=len(seq),
        clip_side=side,
        read_id=read_id,
        query_start=0,
        query_end=len(seq),
        mapq=60,
        sequence=seq,
        left_flank="",
        right_flank="",
    )


def _run(
    sigs: list[Signature], library: Path, tmp_path: Path
) -> tuple[list[TeHit], dict[str, int]]:
    stats: dict[str, int] = {}
    hits = annotate_signatures(sigs, library, threads=2, workdir=tmp_path, stats=stats)
    return hits, stats


# --- unit ----------------------------------------------------------------


def test_command_uses_measured_flags() -> None:
    assert te_minimap2_command(Path("lib.fa"), Path("q.fa"), 4) == [
        "minimap2", "-c", "-x", "map-ont", "-k", "11", "-w", "5", "-r",
        "100,100", "--secondary=no", "-t", "4", "lib.fa", "q.fa",
    ]  # fmt: skip


def test_command_rejects_zero_threads() -> None:
    with pytest.raises(ValueError, match="threads"):
        te_minimap2_command(Path("lib.fa"), Path("q.fa"), 0)


PAF = "s1\t440\t5\t435\t-\tmPing\t430\t0\t430\t425\t431\t60\ttp:A:P\n"


def test_parse_primary_line() -> None:
    hit = parse_paf_line(PAF, {"s1": "r1"})
    assert hit == TeHit(
        signature_id="s1", read_id="r1", family="mPing", te_strand="-",
        te_start=0, te_end=430, te_length=430, query_start=5, query_end=435,
        query_length=440, identity=round(425 / 431, 4),
    )  # fmt: skip


def test_parse_secondary_line_is_none() -> None:
    assert parse_paf_line(PAF.replace("tp:A:P", "tp:A:S"), {"s1": "r1"}) is None


# --- integration (real minimap2) -------------------------------------------


@needs_minimap2
def test_mutated_short_te_hits_forward(
    library: Path, tes: dict[str, str], tmp_path: Path
) -> None:
    rng = random.Random(11)
    seq = _mutate(tes["mGing"], 0.05, rng) + _rand(5, rng)
    hits, _ = _run([_sig("s1", "r1", seq)], library, tmp_path)
    assert len(hits) == 1
    assert (hits[0].family, hits[0].te_strand) == ("mGing", "+")
    assert hits[0].query_length == len(seq)


@needs_minimap2
def test_reverse_complement_hits_minus(
    library: Path, tes: dict[str, str], tmp_path: Path
) -> None:
    hits, _ = _run([_sig("s1", "r1", _revcomp(tes["mPing"]))], library, tmp_path)
    assert [(h.family, h.te_strand) for h in hits] == [("mPing", "-")]


@needs_minimap2
def test_random_sequence_is_no_hit(library: Path, tmp_path: Path) -> None:
    seq = _rand(500, random.Random(3))
    hits, stats = _run([_sig("s1", "r1", seq)], library, tmp_path)
    assert len(hits) == 1
    h = hits[0]
    assert (h.family, h.te_strand, h.query_length) == (NO_HIT, ".", 500)
    assert (stats["with_hit"], stats["without_hit"], stats["hits"]) == (0, 1, 0)


@needs_minimap2
def test_nested_insert_splits(
    library: Path, tes: dict[str, str], tmp_path: Path
) -> None:
    retro, mping = tes["Retro1"], tes["mPing"]
    seq = retro[:500] + mping + retro[500:]
    hits, _ = _run([_sig("s1", "r1", seq)], library, tmp_path)
    families = {h.family for h in hits}
    assert families == {"Retro1", "mPing"}
    (m,) = [h for h in hits if h.family == "mPing"]
    assert abs(m.query_start - 500) <= 10
    assert abs(m.query_end - 930) <= 10


@needs_minimap2
def test_clip_signature_is_annotated(
    library: Path, tes: dict[str, str], tmp_path: Path
) -> None:
    sig = _sig("c1", "r1", tes["Retro1"][:300], kind=SignatureKind.CLIP,
               side=ClipSide.RIGHT)  # fmt: skip
    hits, _ = _run([sig], library, tmp_path)
    assert [h.family for h in hits] == ["Retro1"]


@needs_minimap2
def test_empty_sequence_is_no_hit_and_not_queried(
    library: Path, tes: dict[str, str], tmp_path: Path
) -> None:
    sigs = [_sig("e1", "r1", ""), _sig("s2", "r2", tes["mPing"])]
    hits, stats = _run(sigs, library, tmp_path)
    assert hits[0] == TeHit("e1", "r1", NO_HIT, ".", 0, 0, 0, 0, 0, 0, 0.0)
    assert hits[1].family == "mPing"
    assert stats["queries"] == 1
    assert stats["signatures"] == 2


@needs_minimap2
def test_order_follows_input_then_query_position(
    library: Path, tes: dict[str, str], tmp_path: Path
) -> None:
    retro, mping = tes["Retro1"], tes["mPing"]
    nested = retro[:500] + mping + retro[500:]
    sigs = [
        _sig("z", "r1", tes["mGing"]),
        _sig("a", "r2", _rand(500, random.Random(5))),
        _sig("m", "r3", nested),
    ]
    hits, _ = _run(sigs, library, tmp_path)
    ids = [h.signature_id for h in hits]
    first = {sid: ids.index(sid) for sid in ("z", "a", "m")}
    assert first["z"] < first["a"] < first["m"]
    assert ids == sorted(ids, key=lambda s: first[s])
    m_hits = [h for h in hits if h.signature_id == "m"]
    keys = [(h.query_start, h.query_end, h.family) for h in m_hits]
    assert keys == sorted(keys)
    assert len(m_hits) >= 2


def test_duplicate_signature_id_raises(library: Path, tmp_path: Path) -> None:
    sigs = [_sig("s1", "r1", "ACGT"), _sig("s1", "r2", "ACGT")]
    with pytest.raises(ValueError, match="s1"):
        annotate_signatures(sigs, library, threads=1, workdir=tmp_path)


def test_zero_signatures_skips_minimap2(
    library: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("PATH", "")  # minimap2 would not be found
    stats: dict[str, int] = {}
    assert annotate_signatures([], library, threads=1, workdir=tmp_path,
                               stats=stats) == []  # fmt: skip
    assert stats == {"signatures": 0, "queries": 0, "with_hit": 0,
                     "without_hit": 0, "hits": 0}  # fmt: skip


def test_missing_library_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        annotate_signatures([_sig("s1", "r1", "ACGT")], tmp_path / "nope.fa",
                            threads=1, workdir=tmp_path)  # fmt: skip


def test_missing_minimap2_raises(
    library: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="minimap2"):
        annotate_signatures([_sig("s1", "r1", "ACGT")], library, threads=1,
                            workdir=tmp_path)  # fmt: skip


@needs_minimap2
def test_stats_and_no_leftovers(
    library: Path, tes: dict[str, str], tmp_path: Path
) -> None:
    retro, mping = tes["Retro1"], tes["mPing"]
    sigs = [
        _sig("s1", "r1", tes["mPing"]),
        _sig("s2", "r2", _rand(500, random.Random(9))),
        _sig("s3", "r3", ""),
        _sig("s4", "r4", retro[:500] + mping + retro[500:]),
    ]
    hits, stats = _run(sigs, library, tmp_path)
    assert set(stats) == {"signatures", "queries", "with_hit", "without_hit",
                          "hits"}  # fmt: skip
    real = [h for h in hits if h.family != NO_HIT]
    assert stats["signatures"] == 4
    assert stats["queries"] == 3
    assert stats["with_hit"] == 2
    assert stats["without_hit"] == 2
    assert stats["hits"] == len(real)
    assert not list(tmp_path.glob(".lr-annotate.*"))


def test_minimap2_failure_names_command(
    library: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = tmp_path / "bin"
    fake.mkdir()
    script = fake / "minimap2"
    script.write_text("#!/bin/sh\necho 'boom: bad index' >&2\nexit 3\n")
    script.chmod(0o755)
    monkeypatch.setenv("PATH", str(fake))
    work = tmp_path / "work"
    work.mkdir()
    with pytest.raises(RuntimeError, match=r"(?s)minimap2.*exit 3.*boom"):
        annotate_signatures([_sig("s1", "r1", "ACGT")], library, threads=1,
                            workdir=work)  # fmt: skip
    assert not list(work.glob(".lr-annotate.*"))
