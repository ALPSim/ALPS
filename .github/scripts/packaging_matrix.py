"""Expand all supported wheel builds and Python runtime smoke tests."""

import json
import os
from pathlib import Path


def packaging_matrices():
    platforms = [
        {"os": "ubuntu-24.04", "target": "", "arch": "x86_64", "family": "manylinux", "build": "cp3{11,12,13,14}-manylinux*"},
        {"os": "ubuntu-24.04", "target": "", "arch": "x86_64", "family": "musllinux", "build": "cp3{11,12,13,14}-musllinux*"},
        {"os": "macos-15", "target": "15.0", "arch": "arm64", "family": "macos", "build": "cp3{11,12,13,14}-macosx*"},
    ]
    smoke_platforms = [
        {"os": "ubuntu-24.04", "architecture": "x64", "artifact": "cibw-wheels-manylinux"},
        {"os": "macos-15", "architecture": "arm64", "artifact": "cibw-wheels-macos"},
        {"os": "macos-26", "architecture": "arm64", "artifact": "cibw-wheels-macos"},
    ]
    return {
        "wheel_matrix": {"plat": platforms},
        "smoke_matrix": {
            "plat": smoke_platforms,
            "python": ["3.11", "3.12", "3.13", "3.14"],
        },
    }


if __name__ == "__main__":
    output = "".join(f"{key}={json.dumps(value, separators=(',', ':'))}\n"
                     for key, value in packaging_matrices().items())
    print(output, end="")
    if path := os.environ.get("GITHUB_OUTPUT"):
        with Path(path).open("a", encoding="utf-8") as stream:
            stream.write(output)
