# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PredSim'
copyright = "%Y, KU Leuven Neuromechanics"
author = "Lars D'Hondt"
github_username = 'KULeuvenNeuromechanics'
github_repository = 'https://github.com/KULeuvenNeuromechanics/PredSim'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser']
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_theme_options = {
    'fixed_sidebar': 'true', 
    'globaltoc_maxdepth': '3',
    'sidebar_collapse': 'false',
    'sidebar_width': '250px',
    'page_width': '90%',
    'github_user': 'KULeuvenNeuromechanics',
    'github_repo': 'PredSim',
    #'github_banner': 'true',
    #'github_button': 'true',
}

html_sidebars = {
    #'**': ['navigation.html', 'searchfield.html'],
    #'**': ['searchbox.html','globaltoc.html'],
    '**': ['searchfield.html','globaltoc.html'],
}

html_static_path = ['_static']
#html_logo = '_static/PredSimlogo_notext.png'
html_logo = '_static/PredSimlogo_text.png'
html_favicon = '_static/PredSimlogo_notext (small).png'
