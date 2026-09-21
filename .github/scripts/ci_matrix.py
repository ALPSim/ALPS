"""Select the routine or full source-build matrix from one manifest."""

import argparse
import json
import os
from pathlib import Path
import re


def select_matrix(manifest, tier):
    defaults = {
        "quick": False, "boost": "1.91.0", "standard": 17,
        "packages": "", "repository": "", "extras": False, "sanitizer": False,
        "extensive": False,
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
        build["mpi"] = "OFF" if build["sanitizer"] else "ON"
        if tier == "full" or build["quick"]:
            builds.append(build)
    if not builds:
        raise ValueError(f"Empty {tier} matrix")
    return {"include": builds}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tier", choices=("auto", "quick", "full"), default="auto")
    args = parser.parse_args()
    tier = args.tier
    if tier == "auto":
        periodic = os.environ.get("GITHUB_EVENT_NAME") == "schedule"
        release = os.environ.get("GITHUB_REF", "").startswith("refs/tags/v")
        tier = "full" if periodic or release else "quick"
    manifest = json.loads((Path(__file__).resolve().parents[1] / "ci-matrix.json").read_text())
    matrix = select_matrix(manifest, tier)
    serialized = json.dumps(matrix, separators=(",", ":"))
    print(serialized)
    if output := os.environ.get("GITHUB_OUTPUT"):
        with Path(output).open("a", encoding="utf-8") as stream:
            stream.write(f"matrix={serialized}\n")
    if summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        with Path(summary).open("a", encoding="utf-8") as stream:
            stream.write(f"Source matrix: **{tier}**, {len(matrix['include'])} builds.\n")


if __name__ == "__main__":
    main()
