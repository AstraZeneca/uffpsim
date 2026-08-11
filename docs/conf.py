# Configuration file for the Sphinx documentation builder.
import os
import sys

sys.path.insert(0, os.path.abspath('..'))

project = 'uffpsim'
copyright = '2026, Rajendra Kumar, Gian Marco Ghiandoni, Young Mi Park, Prakash Chandra Rathi'
author = 'Rajendra Kumar, Gian Marco Ghiandoni'

try:
    from uffpsim import __version__
    release = __version__
except Exception:
    release = '0.0.1'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_inline_tabs',
]

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
source_encoding = 'utf-8'

autodoc_mock_imports = [
    'uffpsimLib',
    'rdkit',
    'rdkit.Chem',
    'rdkit.Chem.Draw',
    'rdkit.Chem.rdMolDescriptors',
    'rdkit.Chem.rdFingerprintGenerator',
    'rdkit.Avalon',
    'flask',
]

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'style_nav_header_background': '#2980B9',
}
html_domain_indices = False
master_doc = 'index'
language = 'en'
html_static_path = ['_static']
autodoc_typehints = 'description'
