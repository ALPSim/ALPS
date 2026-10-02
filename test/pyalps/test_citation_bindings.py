# SPDX-License-Identifier: MIT
"""Exercise saved citation metadata through the installed Python bindings."""

from pathlib import Path

import pytest

import pyalps
from pyalps.hdf5 import archive
from pyalps.ngs import params, results, saveResults


def test_save_results_preserves_readable_citations(tmp_path):
    source = tmp_path / "results.h5"
    parameters = params()
    parameters["SEED"] = 17
    with archive(str(source), "w") as output:
        # Empty result sets still record the build that produced the file.
        saveResults(results(), parameters, output, "/simulation/results")
        first = pyalps.read_citations(output)
        saveResults(results(), parameters, output, "/simulation/results")
        assert pyalps.read_citations(output) == first

    assert len(first) == 1
    record = first[0]
    assert record["component"] == "framework"
    assert record["activity"] == "unspecified"
    assert record["framework"]
    assert record["algorithm"] == []

    before = source.read_bytes()
    assert pyalps.read_citations(source) == first
    destination = tmp_path / "CITATION.cff"
    assert pyalps.export_citations(source, destination) == destination
    assert destination.read_text(encoding="utf-8") == record["bibliography_cff"]
    assert source.read_bytes() == before
    with pytest.raises(FileExistsError):
        pyalps.export_citations(source, destination)

    with archive(str(source), "a") as output:
        path = f"/provenance/alps/citations/records/{record['id']}/notice"
        output.delete_data(path)
        output[path] = "changed"
    with pytest.raises(ValueError, match="fingerprint mismatch"):
        pyalps.read_citations(source)


def test_legacy_file_has_no_invented_citations(tmp_path):
    source = tmp_path / "legacy.h5"
    with archive(str(source), "w") as output:
        output["/parameters/SEED"] = 17
    before = source.read_bytes()
    assert pyalps.read_citations(source) == []
    assert source.read_bytes() == before


def test_installed_citation_catalog():
    catalog = Path(pyalps.__file__).resolve().parent / "share/alps"
    for filename in ("CITATION.cff", "CITATIONS.yaml", "CITATION.md"):
        assert (catalog / filename).is_file(), filename
