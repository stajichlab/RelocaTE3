"""Extraction checks require no benchmark files or large scans."""
import hashlib
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "trace_line", Path(__file__).parents[1] / "scripts/trace_line_families.py")
trace = importlib.util.module_from_spec(spec)
spec.loader.exec_module(trace)


def test_exact_names_duplicate_hits_and_hash(tmp_path):
    path = tmp_path / "mapping.tsv"
    data = b"read/1\tA\nread/10\tB\nread/1\tC\n"
    path.write_bytes(data)
    hashes = {}
    rows = trace.select_lines(path, {"read/1"}, 0, hashes)
    assert [r["line_number"] for r in rows] == [1, 3]
    assert hashes[str(path)] == hashlib.sha256(data).hexdigest()


def test_psl_column_and_short_header(tmp_path):
    path = tmp_path / "hits.psl"
    path.write_text("header\n" + "\t".join(["0"] * 9 + ["wanted", "100"]) + "\n")
    rows = trace.select_lines(path, {"wanted"}, 9, {})
    assert len(rows) == 1 and rows[0]["line_number"] == 2
