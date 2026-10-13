# SPDX-License-Identifier: MIT
"""Check citation catalogs bundled with the installed Python package."""

from pyalps._resources import runtime_directory


def test_installed_citation_catalog():
    catalog = runtime_directory() / "share/alps"
    for filename in ("CITATION.cff", "CITATIONS.yaml", "CITATION.md"):
        assert (catalog / filename).is_file(), filename
