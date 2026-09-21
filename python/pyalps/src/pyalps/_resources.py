# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
"""Locate installed runtime resources, including in an editable installation."""

from importlib.resources import files
from pathlib import Path


def runtime_directory():
    # Editable packages combine live Python sources and CMake-installed resources.
    # The manifest belongs to the latter and anchors all native runtime paths.
    manifest = files(__package__).joinpath("runtime.json")
    if not manifest.is_file():
        raise ImportError("pyalps runtime is missing; install the package with pip")
    return Path(str(manifest)).parent
