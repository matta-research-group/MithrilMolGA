# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/stable/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# Incase the project was not installed
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import MithrilMolGA


# -- Project information -----------------------------------------------------

project = 'MithrilMolGA'
author = 'Tristan Stephens-Jones, Micaela Matta'
release = '0.1.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.extlinks',
]

autosummary_generate = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True


templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'imported-members': False,
    'private-members': False,
    'show-inheritance': True,
}
autodoc_member_order = 'bysource'

# Mock imports for external dependencies
autodoc_mock_imports = [
    'rdkit',
    'numpy',
    'scipy',
    'pandas',
    'matplotlib',
    'psi4',
]

templates_path = ['_templates']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'default'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


htmlhelp_basename = 'MithrilMolGAdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'MithrilMolGA.tex', 'MithrilMolGA Documentation',
     'MithrilMolGA', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'MithrilMolGA', 'MithrilMolGA Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'MithrilMolGA', 'MithrilMolGA Documentation',
     author, 'MithrilMolGA', 'A Genetic Algorithm (GA) package for discovery of OMIECs',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------
