# Configuration file for the Sphinx documentation builder.

from sitesyncro import __version__

project = 'SiteSyncro'
copyright = '2024, Peter Demján'
author = 'Peter Demján'
release = __version__

extensions = ['sphinx.ext.autodoc', 'myst_parser']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

source_suffix = {
	'.rst': 'restructuredtext',
	'.md': 'myst_parser',
}

html_theme = 'nature'
html_static_path = ['_static']

# Add a link to your GitHub repository
html_context = {
	'display_github': True,
	'github_user': 'demjanp',
	'github_repo': 'SiteSyncro',
	'github_version': 'master/',
}

html_sidebars = {
	'**': ['localtoc.html', 'project_home.html', 'searchbox.html'],
}