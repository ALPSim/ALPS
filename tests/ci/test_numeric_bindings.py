"""Keep the vendored subset closed under all conditional include branches."""
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]
INCLUDE = ROOT / "third_party/boost_numeric_bindings/include"
PREFIX = "boost/numeric/bindings/"
SOURCE_SUFFIXES = {".h", ".hpp", ".ipp", ".cpp", ".c", ".cc", ".cxx",
                   ".hh", ".in", ".tpp", ".tcc"}


def binding_includes(path):
    text = re.sub(r"/\*.*?\*/", "", path.read_text(errors="replace"), flags=re.S)
    # Inspect both sides of #if: the local compiler cannot cover every backend.
    return set(re.findall(
        r'^\s*#\s*include\s*[<"](boost/numeric/bindings/[^>"\n]+)[>"]',
        text, re.M))


def test_numeric_bindings_include_closure():
    # Include new sources during development, but exclude ignored build products.
    sources = subprocess.check_output([
        "git", "ls-files", "--cached", "--others", "--exclude-standard", "-z",
        "--", "src", "tests", "tutorials", "python", ".github",
    ], cwd=ROOT).decode().split("\0")
    pending = set()
    for name in sources:
        path = ROOT / name
        if path.is_file() and path.suffix.lower() in SOURCE_SUFFIXES:
            pending.update(binding_includes(path))
    assert pending, "No Numeric Bindings consumers found"
    visited = set()
    while pending:
        name = pending.pop()
        if name in visited:
            continue
        path = INCLUDE / name
        assert path.is_file(), f"Missing vendored include: {name}"
        visited.add(name)
        pending.update(binding_includes(path) - visited)

    headers = {path.relative_to(INCLUDE).as_posix()
               for path in (INCLUDE / PREFIX).rglob("*") if path.is_file()}
    assert headers == visited, f"Unreferenced vendored headers: {sorted(headers - visited)}"
