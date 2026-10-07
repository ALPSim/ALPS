"""Expand the complete source matrix without event- or path-based reductions."""

import json
import os
from pathlib import Path
import re


def select_matrix(manifest):
    defaults = {
        "boost": manifest["boost_default"], "standard": 17,
        "packages": "", "repository": "", "extras": False,
        "extensive": False, "python": "3.14", "mpi": "ON",
    }
    builds = []
    identifiers = set()
    for entry in manifest["builds"]:
        build = defaults | entry
        identifier = build["id"]
        if not re.fullmatch(r"[a-z0-9-]+", identifier) or identifier in identifiers:
            raise ValueError(f"Invalid or duplicate build id: {identifier}")
        identifiers.add(identifier)
        checksum = manifest["boost"][build["boost"]]
        if not re.fullmatch(r"[a-f0-9]{64}", checksum):
            raise ValueError(f"Invalid Boost checksum for {identifier}")
        build["boost_sha256"] = checksum
        build["dependency"] = f"boost-{build['boost']}"
        builds.append(build)
    if not builds:
        raise ValueError("Empty source matrix")
    return {"include": builds}


def main():
    manifest = json.loads((Path(__file__).resolve().parents[1] / "ci-matrix.json").read_text())
    matrix = select_matrix(manifest)
    output = "matrix=" + json.dumps(matrix, separators=(",", ":")) + "\n"
    print(output, end="")
    if path := os.environ.get("GITHUB_OUTPUT"):
        with Path(path).open("a", encoding="utf-8") as stream:
            stream.write(output)
    if summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        with Path(summary).open("a", encoding="utf-8") as stream:
            stream.write(f"Source matrix: all {len(matrix['include'])} builds.\n")


if __name__ == "__main__":
    main()
