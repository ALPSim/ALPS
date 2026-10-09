# SPDX-License-Identifier: MIT
"""Check citation catalogs bundled with the installed Python package."""

from pathlib import Path

import pyalps


def test_installed_citation_catalog():
    catalog = Path(pyalps.__file__).resolve().parent / "share/alps"
    for filename in ("CITATION.cff", "CITATIONS.yaml", "CITATION.md"):
        assert (catalog / filename).is_file(), filename
