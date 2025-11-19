# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PredSim'
copyright = "%Y, KU Leuven Neuromechanics"

author = "Lars D'Hondt" # add you name when you contribute

github_username = 'KULeuvenNeuromechanics'
github_repository = 'https://github.com/KULeuvenNeuromechanics/PredSim'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_parser', # to use markdown source files
    'sphinx_design', # extra formatting
    ]

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}


# directory in source with images
html_static_path = ['_static']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# template for layout
html_theme = 'alabaster'

# configure template
html_theme_options = {
    'fixed_sidebar': 'true', # sidebar stays on screen when scrolling down page
    'globaltoc_maxdepth': '3', # table of contents in sidebar shows max 3 levels
    'logo': 'PredSimlogo_text.png',
    'touch_icon': 'PredSimlogo_notext (small).png',
    'github_user': 'KULeuvenNeuromechanics',
    'github_repo': 'PredSim',
    'show_powered_by': 'false',
}

# contents of sidebar (in this order)
html_sidebars = {
    '**': [
        'about.html', # logo, github button
        'searchfield.html', # search
        'globaltoc.html', # table of contents
    ]
}

# icon shown in browser tab title and favourites
html_favicon = '_static/PredSimlogo_notext (small).png'
