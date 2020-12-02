#!/usr/bin/env python3
# Configuration file for the Sphinx documentation builder.
# For a full list of options see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Module imports ----------------------------------------------------------
import os
import sys
import sphinx_rtd_theme
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'Variational Mesh'
author = 'Wanja Schulze'
copyright = '2020, Wanja Schulze'
version = '1.0'
release = '1.0.0'


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme'
]
exclude_patterns = []
templates_path = ['_templates']
pygments_style = 'sphinx'


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
