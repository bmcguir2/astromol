from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import astromol  # noqa: E402


project = "astromol"
author = "Brett A. McGuire"
copyright = "2026, Brett A. McGuire"
release = astromol.__version__
version = astromol.__version__

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
]

source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}
master_doc = "index"

autosummary_generate = True
autodoc_member_order = "bysource"
napoleon_google_docstring = True
napoleon_numpy_docstring = True

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

html_theme = "sphinx_rtd_theme"
html_title = "astromol"
html_static_path: list[str] = []
exclude_patterns = [
    "_build",
    "**/.ipynb_checkpoints",
    "**/.ipynb_checkpoints/*",
    "Thumbs.db",
    ".DS_Store",
]
