"""Check .github/website-pages.json and the page listing in website_pages.py."""

import importlib.util
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def load_helper():
    spec = importlib.util.spec_from_file_location("website_pages", ROOT / ".github/scripts/website_pages.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_map_matches_tracked_files():
    helper = load_helper()
    files = subprocess.run(["git", "ls-files"], cwd=ROOT, stdout=subprocess.PIPE, check=True,
                           text=True).stdout.split()
    assert helper.problems(helper.load(), files) == []


def test_problems_name_stale_patterns_and_bad_pages():
    helper = load_helper()
    page_map = {"entries": [
        {"name": "a", "sources": ["src/a/**"], "pages": ["a.md"]},
        {"name": "b", "sources": ["src/gone/**"], "pages": ["/b.md"]},
        {"name": "c", "sources": [], "pages": ["c.md"]},
    ], "unmapped": {"a.md": "index"}}
    assert helper.problems(page_map, ["src/a/x.cpp"]) == [
        "b: page '/b.md' should be a .md path under content/en",
        "b: 'src/gone/**' matches no tracked file",
        "c: needs both sources and pages",
        "a.md is both mapped and listed as unmapped",
    ]


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


def test_map_is_valid_json_with_the_documented_keys():
    page_map = json.loads((ROOT / ".github/website-pages.json").read_text(encoding="utf-8"))
    assert set(page_map) == {"description", "site", "watched", "ignored", "entries", "unmapped"}
    assert page_map["site"] == "ALPSim/alpsim.github.io"
