"""Tests for ``relocaTE3 cohort``: cohort call table, genotype matrix, summary."""

from __future__ import annotations

import collections
import random
from pathlib import Path

import pytest

from RelocaTE3 import cli
from RelocaTE3.cli import main
from RelocaTE3.cohort import (
    CohortError,
    build_cohort,
    natural_key,
    status_class,
)

HEADER = (
    "strain\tTE\tTSD\tchromosome.pos\tstrand\tavg_flankers\tspanners\t"
    "status\tTE_family_support\tTE_family_confidence\tTE_family_status\t"
    "TE_supporting_family_support\tTE_supporting_family_confidence\t"
    "TE_supporting_family_status\tTE_family_concordance\t"
    "TE_family_candidate_support\tTE_family_ambiguous_reads\tTE_family_resolution"
)

REAL_DIR = Path(
    "/bigdata/stajichlab/nmath020/github/github_tools/RelocaTE/RelocaTE3_jason/"
    "RelocaTE3/results/nextflow-m2/nextflow-real-rice/dev"
)
REAL_SAMPLES = [
    "B_10", "B_120", "B_185", "B_42", "B_68",
    "G_1", "G_104", "G_207", "G_260", "G_52",
]  # fmt: skip


def _row(sample, chrom, start, end, strand="+", status="homozygous", te="mPing", flank="10"):
    return (
        f"{sample}\t{te}\tTAA\t{chrom}:{start}..{end}\t{strand}\t{flank}\t0\t{status}\t"
        f"{te}=5\t1.000000\tunique\t{te}=2\t1.000000\tunique\tconcordant\t\t0\tselected_votes"
    )


def _write(path: Path, rows, header=HEADER) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(line + "\n" for line in [header, *rows]))
    return path


def _run(tmp_path, specs, outdir=None):
    outdir = outdir or tmp_path / "cohort"
    rc = main(["cohort", "-i", *[str(s) for s in specs], "-o", str(outdir)])
    return rc, outdir


def _lines(path: Path) -> list[str]:
    return path.read_text().splitlines()


def _data(path: Path) -> list[str]:
    """Non-comment lines of a TSV."""
    return [line for line in _lines(path) if not line.startswith("#")]


# ---------------------------------------------------------------- CLI wiring


def test_cohort_is_registered_and_dispatches(monkeypatch, tmp_path):
    seen = {}
    monkeypatch.setattr(cli, "cmd_cohort", lambda a: (seen.setdefault("a", a), 0)[1])
    rc = main(["cohort", "-i", "a.txt", "S=b.txt", "-i", "c.txt", "-o", str(tmp_path)])
    assert rc == 0
    assert seen["a"].inputs == ["a.txt", "S=b.txt", "c.txt"]
    assert seen["a"].outdir == str(tmp_path)


def test_cohort_help_mentions_inputs(capsys):
    with pytest.raises(SystemExit) as err:
        cli.build_parser().parse_args(["cohort", "--help"])
    assert err.value.code == 0
    out = capsys.readouterr().out
    assert "-i" in out and "--outdir" in out


# ------------------------------------------------------- cohort call table


def test_cohort_table_orders_by_sample_then_original_row_order(tmp_path):
    # original row order inside a sample is NOT coordinate order: it must be kept
    b = _write(tmp_path / "b.txt", [_row("B", "Chr2", 50, 52), _row("B", "Chr1", 9, 11)])
    a = _write(tmp_path / "a.txt", [_row("A", "Chr10", 5, 7), _row("A", "Chr1", 1, 3)])
    rc, out = _run(tmp_path, [b, a])
    assert rc == 0
    lines = _lines(out / "cohort.characTErized.txt")
    assert lines[0] == HEADER
    assert lines[1:] == _lines(a)[1:] + _lines(b)[1:]


def test_cohort_output_is_independent_of_input_order(tmp_path):
    files = [
        _write(tmp_path / f"{s}.txt", [_row(s, "Chr1", i * 10, i * 10 + 2) for i in range(3)])
        for s in ("S1", "S3", "S2", "S10")
    ]
    _, out1 = _run(tmp_path, files, tmp_path / "o1")
    shuffled = files[:]
    random.Random(7).shuffle(shuffled)
    _, out2 = _run(tmp_path, shuffled, tmp_path / "o2")
    for name in (
        "cohort.characTErized.txt",
        "cohort.genotype_matrix.tsv",
        "cohort.summary_mqc.tsv",
    ):
        assert (out1 / name).read_bytes() == (out2 / name).read_bytes(), name


def test_header_mismatch_fails(tmp_path):
    a = _write(tmp_path / "a.txt", [_row("A", "Chr1", 1, 3)])
    other = HEADER.replace("TE_family_resolution", "something_else")
    b = _write(tmp_path / "b.txt", [_row("B", "Chr1", 1, 3)], header=other)
    rc, out = _run(tmp_path, [a, b])
    assert rc != 0
    assert not (out / "cohort.characTErized.txt").exists()
    with pytest.raises(CohortError, match="header"):
        build_cohort([str(a), str(b)])


def test_header_missing_required_column_fails(tmp_path):
    bad = HEADER.replace("status", "state", 1)
    a = _write(tmp_path / "a.txt", [_row("A", "Chr1", 1, 3)], header=bad)
    with pytest.raises(CohortError, match="status"):
        build_cohort([str(a)])


def test_duplicate_sample_fails(tmp_path):
    a1 = _write(tmp_path / "x" / "a.txt", [_row("A", "Chr1", 1, 3)])
    a2 = _write(tmp_path / "y" / "a.txt", [_row("A", "Chr2", 1, 3)])
    with pytest.raises(CohortError, match="'A'.*more than one input"):
        build_cohort([str(a1), str(a2)])
    rc, _ = _run(tmp_path, [a1, a2])
    assert rc != 0


def test_mixed_strains_in_one_file_fail(tmp_path):
    a = _write(tmp_path / "a.txt", [_row("A", "Chr1", 1, 3), _row("B", "Chr1", 5, 7)])
    with pytest.raises(CohortError, match="more than one strain"):
        build_cohort([str(a)])


def test_header_only_plain_path_fails_but_named_input_is_allowed(tmp_path):
    empty = _write(tmp_path / "e.txt", [])
    a = _write(tmp_path / "a.txt", [_row("A", "Chr1", 1, 3)])
    with pytest.raises(CohortError, match="no data rows.*NAME=PATH"):
        build_cohort([str(a), str(empty)])

    rc, out = _run(tmp_path, [a, f"E={empty}"])
    assert rc == 0
    assert _lines(out / "cohort.characTErized.txt") == [HEADER, _row("A", "Chr1", 1, 3)]
    matrix = _data(out / "cohort.genotype_matrix.tsv")
    assert matrix[0].split("\t")[5:] == ["A", "E"]
    assert matrix[1].split("\t")[5:] == ["homozygous", "."]
    summary = {r.split("\t")[0]: r.split("\t") for r in _data(out / "cohort.summary_mqc.tsv")}
    assert summary["E"][1] == "0"


def test_named_input_must_match_strain_column(tmp_path):
    a = _write(tmp_path / "a.txt", [_row("A", "Chr1", 1, 3)])
    with pytest.raises(CohortError, match="strain"):
        build_cohort([f"B={a}"])
    # the matching name is accepted
    assert build_cohort([f"A={a}"]).samples == ["A"]


def test_empty_file_fails(tmp_path):
    p = tmp_path / "zero.txt"
    p.write_text("")
    with pytest.raises(CohortError, match="empty"):
        build_cohort([f"Z={p}"])


def test_missing_file_fails(tmp_path):
    with pytest.raises(CohortError, match="not found"):
        build_cohort([str(tmp_path / "nope.txt")])


def test_no_inputs_fails():
    with pytest.raises(CohortError):
        build_cohort([])


# ---------------------------------------------------------- genotype matrix


def test_matrix_dot_cells_and_exact_site_keys(tmp_path):
    a = _write(
        tmp_path / "a.txt",
        [
            _row("A", "Chr1", 100, 102, "+", "homozygous"),
            _row("A", "Chr1", 200, 202, "-", "heterozygous"),
        ],
    )
    b = _write(
        tmp_path / "b.txt",
        [
            _row("B", "Chr1", 100, 102, "+", "somatic_insertion"),
            # 1 bp away from A's site: nearby sites are NOT merged
            _row("B", "Chr1", 201, 203, "-", "homozygous/excision_no_footprint"),
            # same position as A's, other strand: a different key
            _row("B", "Chr1", 200, 202, "+", "heterozygous"),
        ],
    )
    rc, out = _run(tmp_path, [a, b])
    assert rc == 0
    rows = [r.split("\t") for r in _data(out / "cohort.genotype_matrix.tsv")]
    assert rows[0] == ["chromosome", "start", "end", "strand", "TE", "A", "B"]
    assert rows[1:] == [
        ["Chr1", "100", "102", "+", "mPing", "homozygous", "somatic_insertion"],
        ["Chr1", "200", "202", "+", "mPing", ".", "heterozygous"],
        ["Chr1", "200", "202", "-", "mPing", "heterozygous", "."],
        ["Chr1", "201", "203", "-", "mPing", ".", "homozygous/excision_no_footprint"],
    ]


def test_matrix_is_a_plain_tsv_without_comment_lines(tmp_path):
    a = _write(tmp_path / "a.txt", [_row("A", "Chr1", 1, 3)])
    _, out = _run(tmp_path, [a])
    lines = _lines(out / "cohort.genotype_matrix.tsv")
    assert not any(line.startswith("#") for line in lines)
    assert len({len(line.split("\t")) for line in lines}) == 1


def test_te_family_is_part_of_the_site_key(tmp_path):
    a = _write(
        tmp_path / "a.txt",
        [_row("A", "Chr1", 1, 3, te="mPing"), _row("A", "Chr1", 1, 3, te="Ping")],
    )
    rc, out = _run(tmp_path, [a])
    assert rc == 0
    rows = [r.split("\t") for r in _data(out / "cohort.genotype_matrix.tsv")[1:]]
    assert [r[4] for r in rows] == ["Ping", "mPing"]


def test_duplicate_site_within_one_sample_fails_loudly(tmp_path):
    a = _write(
        tmp_path / "a.txt",
        [_row("A", "Chr1", 1, 3, status="homozygous"), _row("A", "Chr1", 1, 3, status="heterozygous")],
    )
    with pytest.raises(CohortError, match=r"duplicate site.*'A'.*Chr1:1\.\.3"):
        build_cohort([str(a)])
    rc, out = _run(tmp_path, [a])
    assert rc != 0
    assert not out.exists() or not any(out.iterdir())


def test_unparsable_coordinates_fail(tmp_path):
    a = _write(tmp_path / "a.txt", [_row("A", "Chr1", 1, 3).replace("Chr1:1..3", "Chr1-1-3")])
    with pytest.raises(CohortError, match="chromosome.pos"):
        build_cohort([str(a)])


def test_natural_chromosome_order(tmp_path):
    chroms = ["Chr10", "Chr2", "ChrUn", "Chr1", "scaffold_10", "scaffold_9", "ChrSy", "Chr11"]
    a = _write(tmp_path / "a.txt", [_row("A", c, 5, 7) for c in chroms])
    _, out = _run(tmp_path, [a])
    got = [r.split("\t")[0] for r in _data(out / "cohort.genotype_matrix.tsv")[1:]]
    assert got == ["Chr1", "Chr2", "Chr10", "Chr11", "ChrSy", "ChrUn", "scaffold_9", "scaffold_10"]


def test_natural_key_ties_and_colons():
    assert sorted(["Chr01", "Chr1"], key=natural_key) == ["Chr01", "Chr1"]
    assert natural_key("Chr2") < natural_key("Chr10")


def test_chromosome_names_with_colons_parse(tmp_path):
    a = _write(tmp_path / "a.txt", [_row("A", "HLA-A*01:01", 5, 7)])
    _, out = _run(tmp_path, [a])
    row = _data(out / "cohort.genotype_matrix.tsv")[1].split("\t")
    assert row[:3] == ["HLA-A*01:01", "5", "7"]


def test_matrix_sorts_numerically_by_start_then_end(tmp_path):
    a = _write(tmp_path / "a.txt", [_row("A", "Chr1", 1000, 1002), _row("A", "Chr1", 999, 1001), _row("A", "Chr1", 999, 1000)])
    _, out = _run(tmp_path, [a])
    got = [tuple(r.split("\t")[1:3]) for r in _data(out / "cohort.genotype_matrix.tsv")[1:]]
    assert got == [("999", "1000"), ("999", "1001"), ("1000", "1002")]


# ------------------------------------------------------------ MultiQC summary


def test_status_class_rule():
    assert status_class("homozygous") == ("homozygous", False)
    assert status_class("homozygous/excision_no_footprint") == ("homozygous", True)
    assert status_class("heterozygous/excision_with_footprint") == ("heterozygous", True)
    assert status_class("somatic_insertion") == ("somatic_insertion", False)
    assert status_class("other") == ("other", False)


def test_summary_counts_and_multiqc_header(tmp_path):
    a = _write(
        tmp_path / "a.txt",
        [
            _row("A", "Chr1", 1, 3, status="homozygous"),
            _row("A", "Chr1", 5, 7, status="homozygous/excision_no_footprint"),
            _row("A", "Chr1", 9, 11, status="heterozygous"),
            _row("A", "Chr1", 13, 15, status="somatic_insertion"),
            _row("A", "Chr1", 17, 19, status="heterozygous/excision_with_footprint"),
            _row("A", "Chr1", 21, 23, status="other"),
        ],
    )
    b = _write(tmp_path / "b.txt", [_row("B", "Chr1", 1, 3, status="heterozygous")])
    _, out = _run(tmp_path, [b, a])
    lines = _lines(out / "cohort.summary_mqc.tsv")
    comments = [line for line in lines if line.startswith("#")]
    # MultiQC custom-content header: YAML in '# ' comment lines before the data
    assert lines[: len(comments)] == comments
    text = "\n".join(c[2:] for c in comments)
    assert "id: 'relocate3_cohort_summary'" in text
    assert "plot_type: 'table'" in text
    assert "format: 'tsv'" in text
    assert "section_name:" in text and "description:" in text

    data = [r.split("\t") for r in _data(out / "cohort.summary_mqc.tsv")]
    assert data[0] == [
        "Sample",
        "total_calls",
        "homozygous",
        "heterozygous",
        "somatic_insertion",
        "other",
        "excision",
    ]
    assert data[1] == ["A", "6", "2", "2", "1", "1", "2"]
    assert data[2] == ["B", "1", "0", "1", "0", "0", "0"]


def test_summary_yaml_header_parses(tmp_path):
    yaml = pytest.importorskip("yaml")
    a = _write(tmp_path / "a.txt", [_row("A", "Chr1", 1, 3)])
    _, out = _run(tmp_path, [a])
    comments = [l[1:] for l in _lines(out / "cohort.summary_mqc.tsv") if l.startswith("#")]
    conf = yaml.safe_load("\n".join(c[1:] if c.startswith(" ") else c for c in comments))
    assert conf["id"] == "relocate3_cohort_summary"
    assert conf["plot_type"] == "table"
    assert conf["pconfig"]["id"]


def test_summary_unknown_genotype_class_gets_its_own_column(tmp_path):
    a = _write(tmp_path / "a.txt", [_row("A", "Chr1", 1, 3, status="weird/excision_x")])
    _, out = _run(tmp_path, [a])
    data = [r.split("\t") for r in _data(out / "cohort.summary_mqc.tsv")]
    assert data[0][-2:] == ["excision", "weird"]
    assert data[1] == ["A", "1", "0", "0", "0", "0", "1", "1"]


# ------------------------------------------------------------- real rice


@pytest.mark.skipif(not REAL_DIR.is_dir(), reason="real-rice M2 outputs not present")
def test_real_rice_cohort_matches_per_sample_files(tmp_path):
    files = {
        s: REAL_DIR / s / "results" / "ALL.mping.all_nonref_insert.characTErized.txt"
        for s in REAL_SAMPLES
    }
    rc, out = _run(tmp_path, list(files.values()))
    assert rc == 0

    cohort = _lines(out / "cohort.characTErized.txt")
    per_sample = {s: _lines(p) for s, p in files.items()}
    assert len(cohort) == 1 + sum(len(v) - 1 for v in per_sample.values())
    assert cohort[0] == per_sample["G_1"][0]
    # each sample's block equals its file's rows, exactly and in order
    by_sample = collections.defaultdict(list)
    for line in cohort[1:]:
        by_sample[line.split("\t", 1)[0]].append(line)
    assert sorted(by_sample) == sorted(REAL_SAMPLES)
    for s, lines in per_sample.items():
        assert by_sample[s] == lines[1:], s
    # and as a multiset over the whole cohort
    expected = collections.Counter(line for v in per_sample.values() for line in v[1:])
    assert collections.Counter(cohort[1:]) == expected

    matrix = [r.split("\t") for r in _data(out / "cohort.genotype_matrix.tsv")]
    assert matrix[0][5:] == sorted(REAL_SAMPLES)
    n_called = [sum(c != "." for c in r[5:]) for r in matrix[1:]]
    # every per-sample call lands in exactly one matrix cell
    assert sum(n_called) == len(cohort) - 1
    shared = sum(n >= 2 for n in n_called)
    print(f"real-rice matrix: {len(matrix) - 1} sites x {len(matrix[0]) - 5} samples; "
          f"{shared} sites called in >=2 samples")

    summary = {r.split("\t")[0]: r.split("\t") for r in _data(out / "cohort.summary_mqc.tsv")[1:]}
    for s, lines in per_sample.items():
        assert int(summary[s][1]) == len(lines) - 1
