"""Coarse CI selection. Unknown paths and unavailable diffs get full coverage."""

import json
import os
from pathlib import Path
import subprocess


def changed_paths(event, payload):
    if event == "pull_request":
        base = payload.get("pull_request", {}).get("base", {}).get("sha")
    elif event == "merge_group":
        base = payload.get("merge_group", {}).get("base_sha")
    else:
        return None
    if not base:
        return None
    try:
        result = subprocess.run(
            ["git", "diff", "--name-only", "--no-renames", "-z", base, "HEAD", "--"],
            check=True, capture_output=True,
        )
    except subprocess.CalledProcessError:
        return None
    return result.stdout.decode("utf-8", errors="surrogateescape").strip("\0").split("\0")


def documentation(path):
    return (path.endswith((".md", ".rst"))
            or path in {"LICENSE.txt", "CITATION.md"}
            or path.startswith(".github/ISSUE_TEMPLATE/"))


def select(event, paths):
    if event not in {"pull_request", "merge_group"} or paths is None:
        return {"source": True, "packaging": "full"}
    paths = [path for path in paths if not documentation(path)]
    if not paths:
        return {"source": False, "packaging": "none"}
    # Build configuration and public headers can affect every binary platform.
    broad = any(
        Path(path).name == "CMakeLists.txt"
        or path.endswith((".h", ".hh", ".hpp", ".hxx", ".ipp", ".tpp", ".inl",
                          ".h.in", ".hpp.in", ".hxx.in", ".cmake"))
        or path.startswith((".github/", "cmake/", "third_party/",
                            "python/pyalps/_build_support/"))
        or path in {"CMakePresets.json",
                    "python/pyalps/pyproject.toml"}
        or not path.startswith(("src/", "tests/", "python/", "tutorials/"))
        for path in paths
    )
    python_only = all(path.startswith(("python/", "tests/pyalps/", "tests/packaging/"))
                      for path in paths)
    # Ordinary C++ edits get bindings and installed-SDK tests in the primary
    # source job. Python edits additionally exercise a real repaired wheel.
    python_changes = any(path.startswith(("python/", "tests/pyalps/", "tests/packaging/"))
                         for path in paths)
    return {"source": broad or not python_only,
            "packaging": "full" if broad else "quick" if python_changes else "none"}


def current():
    event = os.environ.get("GITHUB_EVENT_NAME", "")
    event_path = os.environ.get("GITHUB_EVENT_PATH")
    payload = json.loads(Path(event_path).read_text()) if event_path else {}
    return select(event, changed_paths(event, payload))


def emit(values):
    lines = []
    for key, value in values.items():
        serialized = json.dumps(value, separators=(",", ":")) if not isinstance(value, str) else value
        lines.append(f"{key}={serialized}")
    print("\n".join(lines))
    if output := os.environ.get("GITHUB_OUTPUT"):
        with Path(output).open("a", encoding="utf-8") as stream:
            stream.write("\n".join(lines) + "\n")


def packaging_matrices(tier):
    platforms = [
        {"os": "ubuntu-24.04", "target": "", "arch": "x86_64", "family": "manylinux", "build": "cp3{10,11,12,13,14}-manylinux*"},
        {"os": "ubuntu-24.04", "target": "", "arch": "x86_64", "family": "musllinux", "build": "cp3{10,11,12,13,14}-musllinux*"},
        {"os": "macos-15", "target": "15.0", "arch": "arm64", "family": "macos", "build": "cp3{10,11,12,13,14}-macosx*"},
    ]
    smoke_platforms = [
        {"os": "ubuntu-24.04", "architecture": "x64", "artifact": "cibw-wheels-manylinux"},
        {"os": "macos-15", "architecture": "arm64", "artifact": "cibw-wheels-macos"},
        {"os": "macos-26", "architecture": "arm64", "artifact": "cibw-wheels-macos"},
    ]
    return {
        "wheel_matrix": {"plat": platforms if tier == "full" else platforms[:1]},
        "smoke_matrix": {
            "plat": smoke_platforms if tier == "full" else smoke_platforms[:1],
            "python": ["3.10", "3.11", "3.12", "3.13", "3.14"] if tier == "full" else ["3.10", "3.14"],
        },
    }


if __name__ == "__main__":
    policy = current()
    emit(policy | packaging_matrices(policy["packaging"]))
