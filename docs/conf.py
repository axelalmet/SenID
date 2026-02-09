from __future__ import annotations

import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath("../src"))

project = 'SenID'
author = 'Axel A. Almet'
release = '0.1.1'

root_doc = "index"


extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
]

html_theme = "sphinx_rtd_theme"

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "tasklist",
]

myst_heading_anchors = 3  # nice for deep-linking

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
napoleon_google_docstring = True
napoleon_numpy_docstring = True

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