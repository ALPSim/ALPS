"""Check the page listing in .github/scripts/website_pages.py."""

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def load_helper():
    spec = importlib.util.spec_from_file_location("website_pages", ROOT / ".github/scripts/website_pages.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_stale_patterns_are_named():
    helper = load_helper()
    page_map = {"entries": [
        {"name": "a", "sources": ["src/a/**", "src/a.h"], "pages": ["a.md"]},
        {"name": "b", "sources": ["src/gone/**"], "pages": ["b.md"]},
    ]}
    gone = helper.stale(page_map, ["src/a/x.cpp", "src/a.h"])
    assert gone == [("b", "src/gone/**")]
    assert "- b: `src/gone/**`" in helper.report([], gone)


def test_affected_entries_list_the_changed_sources():
    helper = load_helper()
    page_map = {"ignored": ["**/tests/**"], "entries": [
        {"name": "fulldiag", "sources": ["src/apps/diag/fulldiag/**", "src/apps/diag/diag.h"],
         "pages": ["ed.md", "fulldiag.md"]},
        {"name": "sparsediag", "sources": ["src/apps/diag/sparsediag/**", "src/apps/diag/diag.h"],
         "pages": ["ed.md", "sparsediag.md"]},
        {"name": "dmrg", "sources": ["src/apps/dmrg/**"], "pages": ["dmrg.md"]},
    ]}
    changed = ["src/apps/diag/fulldiag/main.C", "src/apps/diag/diag.h", "README.md",
               "src/apps/dmrg/tests/test.py"]
    found = helper.affected(page_map, changed)
    assert found == [
        ("fulldiag", ["src/apps/diag/diag.h", "src/apps/diag/fulldiag/main.C"], ["ed.md", "fulldiag.md"]),
        ("sparsediag", ["src/apps/diag/diag.h"], ["ed.md", "sparsediag.md"]),
    ]
    assert helper.affected(page_map, ["README.md", "src/apps/dmrg/tests/test.py"]) == []
    assert "No website pages" in helper.report([])
    text = helper.report(found)
    assert "3 page(s)" in text
    assert "| fulldiag | `src/apps/diag/diag.h`, `src/apps/diag/fulldiag/main.C` | [ed](" in text
