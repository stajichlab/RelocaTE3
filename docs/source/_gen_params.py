#!/usr/bin/env python3
"""Regenerate the parameter tables of docs/source/pipeline.rst.

The tables come from nextflow_schema.json and sit between the marker comments
``.. BEGIN GENERATED PARAMETERS`` and ``.. END GENERATED PARAMETERS``. The
generated text is committed, so the page reads the same on GitHub (which does
not process ``.. include::``) and in Sphinx. tests/pipeline_docs_test.py fails
when the committed tables differ from a fresh generation.

    python docs/source/_gen_params.py            # rewrite pipeline.rst in place
    python docs/source/_gen_params.py --check    # exit 1 if it is stale

Standard library only; run from any directory.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = ROOT / "nextflow_schema.json"
GUIDE = ROOT / "docs" / "source" / "pipeline.rst"

BEGIN = ".. BEGIN GENERATED PARAMETERS (docs/source/_gen_params.py from nextflow_schema.json; do not edit by hand)"
END = ".. END GENERATED PARAMETERS"


def _escape(text: str) -> str:
    """Make schema prose safe as reStructuredText inline text."""
    text = text.replace("\\", "\\\\")
    for ch in "*`|":
        text = text.replace(ch, "\\" + ch)
    # word_ followed by a space or end would be a hyperlink reference
    return re.sub(r"(\w)_(?=\s|$|[.,;:)])", r"\1\\_", text)


def _type(prop: dict) -> str:
    t = prop.get("type", "")
    return " or ".join(t) if isinstance(t, list) else t


def _default(prop: dict) -> str:
    if "default" not in prop:
        return "(none)"
    d = prop["default"]
    if d == "":
        return "``''``"
    if isinstance(d, bool):
        return f"``{str(d).lower()}``"
    return f"``{d}``"


def render(schema: dict) -> str:
    """The generated block, markers included."""
    out = [BEGIN, ""]
    required = set()
    for group in schema["$defs"].values():
        required.update(group.get("required", []))
    for group in schema["$defs"].values():
        title = group["title"]
        out += [title, "~" * len(title), ""]
        out += [_escape(group["description"]), ""]
        if "help_text" in group:
            out += [_escape(group["help_text"]), ""]
        out += [
            ".. list-table::",
            "   :header-rows: 1",
            "   :widths: 22 10 16 52",
            "",
            "   * - Parameter",
            "     - Type",
            "     - Default",
            "     - Description",
        ]
        for name, prop in group["properties"].items():
            desc = _escape(prop.get("description", ""))
            if "enum" in prop:
                desc += " One of: " + ", ".join(f"``{v}``" for v in prop["enum"]) + "."
            notes = []
            if name in required:
                notes.append("required")
            if prop.get("hidden"):
                notes.append("hidden from ``--help``")
            if notes:
                desc += " (" + "; ".join(notes) + ")"
            out += [
                f"   * - ``--{name}``",
                f"     - {_type(prop)}",
                f"     - {_default(prop)}",
                f"     - {desc}",
            ]
        out.append("")
    out.append(END)
    return "\n".join(out)


def _span(text: str) -> tuple[int, int]:
    start = text.find(BEGIN)
    end = text.find(END)
    if start < 0 or end < start:
        raise ValueError(f"markers not found: {BEGIN!r} ... {END!r}")
    return start, end + len(END)


def generated_region(text: str) -> str:
    start, end = _span(text)
    return text[start:end]


def update_text(text: str, schema: dict) -> str:
    start, end = _span(text)
    return text[:start] + render(schema) + text[end:]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="exit 1 if the page is stale; write nothing")
    ap.add_argument("guide", nargs="?", type=Path, default=GUIDE)
    ap.add_argument("--schema", type=Path, default=SCHEMA)
    args = ap.parse_args(argv)
    text = args.guide.read_text()
    new = update_text(text, json.loads(args.schema.read_text()))
    if args.check:
        if new != text:
            print(f"{args.guide}: parameter tables are stale; run {Path(__file__).name}", file=sys.stderr)
            return 1
        return 0
    if new != text:
        args.guide.write_text(new)
        print(f"updated {args.guide}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
