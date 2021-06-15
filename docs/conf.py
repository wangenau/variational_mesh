#!/usr/bin/env python3
# Configuration file for the Sphinx documentation builder.
# For a full list of options see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sphinx_rtd_theme
import var_mesh

# -- Project information -----------------------------------------------------

project = 'Variational Mesh'
author = 'Wanja Schulze'
copyright = '2020-2021, Wanja Schulze'
version = var_mesh.__version__
release = var_mesh.__version__

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme'
]
templates_path = []
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = []
html_show_sphinx = False
