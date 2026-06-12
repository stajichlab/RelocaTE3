"""Sphinx configuration for RelocaTE3 documentation."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

project = "RelocaTE3"
author = "Jason Stajich, Nathan Mathieu"
copyright = "2024, Jason Stajich, Nathan Mathieu"

try:
    from importlib.metadata import metadata as _meta

    release = _meta("RelocaTE3")["Version"]
except Exception:
    release = "unknown"

version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "alabaster"
html_static_path = ["_static"]
html_theme_options = {
    "description": "TE insertion polymorphism detection from resequencing data",
    "github_user": "stajichlab",
    "github_repo": "RelocaTE3",
    "fixed_sidebar": True,
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pysam": ("https://pysam.readthedocs.io/en/latest/", None),
    "biopython": ("https://biopython.org/docs/latest/api/", None),
}

autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = False
