"""Select the quick, routine, or full source-build matrix from one manifest."""

import argparse
import json
import os
from pathlib import Path
import re


def boost_package(version, runner):
    # GCC 11/libstdc++ and Open MPI 4 binaries also serve newer Ubuntu runners
    # and Clang with libstdc++; the full matrix checks this compatibility.
    platform = {
        "ubuntu-22.04": "linux-x64",
        "ubuntu-24.04": "linux-x64",
        "ubuntu-24.04-arm": "linux-arm64",
        "macos-15": "macos-arm64",
        "macos-26": "macos-arm64",
        "macos-15-intel": "macos-x64",
    }[runner]
    return f"boost-{version}-{platform}"


def select_matrix(manifest, tier, only=""):
    defaults = {
        "quick": False, "boost": manifest["boost_default"], "standard": 17,
        "packages": "", "repository": "", "extras": False, "sanitizer": False,
        "extensive": False, "python": False,
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
        build["dependency"] = boost_package(build["boost"], build["os"])
        build["mpi"] = "OFF" if build["sanitizer"] or build["python"] else "ON"
        if tier == "quick":
            build["extras"] = False
        if tier == "full" or (tier == "routine" and build["quick"]) or (tier == "quick" and build["python"]):
            builds.append(build)
    if only:
        builds = [build for build in builds if build["id"] == only]
    if not builds:
        raise ValueError(f"Empty {tier} matrix (requested build: {only or 'all'})")
    return {"include": builds}


def main():
    from ci_policy import current, emit

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tier", choices=("auto", "quick", "routine", "full"), default="auto")
    parser.add_argument("--only", default="", help="Select one build ID from the chosen tier")
    args = parser.parse_args()
    tier = args.tier
    if tier == "auto":
        periodic = os.environ.get("GITHUB_EVENT_NAME") == "schedule"
        release = os.environ.get("GITHUB_REF", "").startswith("refs/tags/v")
        tier = "full" if periodic or release else (
            "quick" if os.environ.get("GITHUB_EVENT_NAME") in {"pull_request", "merge_group"}
            else "routine")
    manifest = json.loads((Path(__file__).resolve().parents[1] / "ci-matrix.json").read_text())
    matrix = select_matrix(manifest, tier, args.only)
    policy = current()
    changed_build_setup = os.environ.get("GITHUB_EVENT_NAME") in {"pull_request", "merge_group"} and policy["packaging"] == "full"
    emit({"matrix": matrix, "source": policy["source"],
          "developer": os.environ.get("DEVELOPER_SETUP") == "true" or
                       (tier == "full" and not args.only) or changed_build_setup})
    if summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        with Path(summary).open("a", encoding="utf-8") as stream:
            stream.write(f"Source matrix: **{tier}**, {len(matrix['include'])} builds.\n")


if __name__ == "__main__":
    main()
