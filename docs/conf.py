from __future__ import annotations

import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath("../src"))

project = 'SenID'
author = 'Axel A. Almet'
release = '0.1.1'

try:
    from importlib.metadata import version as _version
    release = _version("senid")
except Exception:
    release = "0.1.1"

root_doc = "index"

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "myst_parser",
    "autoapi.extension",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "build", "Thumbs.db", ".DS_Store"]

html_theme = "shibuya"

autoapi_type = "python"
autoapi_dirs = ["../src"]               # point to src/, not src/senid/
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]
autoapi_ignore = [
    "*/_version*",
    "*/tests/*",
    "*-checkpoint*",
]
autoapi_keep_files = True
autoapi_add_toctree_entry = False       # we manage the toctree manually
autoapi_python_class_content = "both"

# -- Autodoc mock imports ----------------------------------------------------
autodoc_mock_imports = [
    "matplotlib",
    "matplotlib.pyplot",
    "scanpy",
    "anndata",
    "squidpy",
    "commot",
    "cnmf",
    "monod",
    "pingouin",
    "joblib",
    "tqdm",
    "numpy",
    "pandas",
    "scipy",
]

autodoc_typehints = "description"

# -- Napoleon (NumPy docstrings) ---------------------------------------------
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# -- MyST (Markdown) settings ------------------------------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "tasklist",
]
myst_heading_anchors = 3

# -- Intersphinx links -------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable", None),
    "anndata": ("https://anndata.readthedocs.io/en/stable", None),
}

# -- HTML output -------------------------------------------------------------
html_theme = "shibuya"

# -- AutoAPI skip handler ----------------------------------------------------
def _autoapi_skip_member(app, what, name, obj, skip, options):
    short_name = name.split(".")[-1]
    if what in ("function", "method", "attribute", "property", "class"):
        if short_name.startswith("_") and not short_name.startswith("__"):
            return True
    if short_name.startswith("__") and short_name not in ("__init__",):
        return True
    return skip

def setup(app):
    app.connect("autoapi-skip-member", _autoapi_skip_member)
