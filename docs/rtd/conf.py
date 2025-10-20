import os
from datetime import date

project = "Magrathea"
author = "Magrathea Developers"
copyright = f"{date.today().year}, {author}"

extensions = [
    "myst_parser",
    "breathe",
    "exhale",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "_doxygen"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# Breathe configuration
breathe_projects = { "Magrathea": "./_doxygen/xml" }
breathe_default_project = "Magrathea"

# Exhale configuration
exhale_args = {
    "containmentFolder":     "./docs_api",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "C++ API",
    "doxygenStripFromPath":  "../../src",
    "createTreeView":        True,
}

myst_enable_extensions = ["colon_fence", "dollarmath", "linkify"]

# Ensure directories exist for local builds
os.makedirs("./_doxygen/xml", exist_ok=True)
os.makedirs("./api", exist_ok=True)
