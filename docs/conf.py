import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

project = 'Rheology Explorer'
copyright = '2026, Tilman May'
author = 'Tilman May'
release = '0.1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
