# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
"""Check pip install/uninstall ownership without compiling ALPS."""

import os
from pathlib import Path
import runpy
import subprocess
import sys
import venv

PROJECT = Path(__file__).resolve().parents[2] / "python/pyalps"
CLI = runpy.run_path(str(PROJECT / "src/pyalps_cli/__init__.py"))
COMMANDS = (*CLI["NATIVE_PROGRAMS"], *CLI["EXPORTERS"])


def _run(command, **kwargs):
    result = subprocess.run(command, text=True, capture_output=True, timeout=60, **kwargs)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def test_bindings_only_install_and_uninstall_preserve_sdk_commands(tmp_path):
    environment = tmp_path / "environment with spaces"
    venv.EnvBuilder(with_pip=False, symlinks=True).create(environment)
    python = environment / "bin/python"
    sdk = python.parent
    for name in COMMANDS:
        executable = sdk / name
        executable.write_text("#!/bin/sh\nexit 0\n")
        executable.chmod(0o755)
    original_commands = {name: (sdk / name).read_bytes() for name in COMMANDS}
    source = tmp_path / "source"
    source.mkdir()
    # Use the production metadata provider: any accidental entry points would
    # overwrite the SDK commands in this same prefix. No native build is needed.
    (source / "pyproject.toml").write_text(f'''
[build-system]
requires = ["scikit-build-core>=1.0"]
build-backend = "scikit_build_core.build"
[project]
name = "pyalps-cli-test"
dynamic = ["version", "scripts"]
[[tool.dynamic-metadata]]
provider = {{ path = "{PROJECT.as_posix()}/_build_support", module = "alps_version" }}
[tool.scikit-build]
wheel.packages = []
''')
    (source / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.22)
project(cli_test LANGUAGES NONE)
''')
    wheels = tmp_path / "wheels"
    wheels.mkdir()
    build_env = {**os.environ, "PYALPS_BUNDLE_APPLICATIONS": "OFF"}
    _run(
        [sys.executable, "-c",
         "from scikit_build_core.build import build_wheel, prepare_metadata_for_build_wheel; "
         "import sys; metadata = prepare_metadata_for_build_wheel('metadata'); "
         "build_wheel(sys.argv[1], metadata_directory='metadata/' + metadata)", str(wheels)],
        cwd=source, env=build_env,
    )
    _run([
        sys.executable, "-m", "pip", "--python", str(python), "install",
        "--no-deps", "--no-index", "--no-cache-dir", str(next(wheels.glob("*.whl"))),
    ])
    _run([str(python), "-I", "-c", "from importlib.metadata import distribution; "
          "assert not distribution('pyalps-cli-test').entry_points"])
    assert {name: (sdk / name).read_bytes() for name in COMMANDS} == original_commands

    # Pip must leave all SDK commands intact on uninstall, too.
    _run([
        sys.executable, "-m", "pip", "--python", str(python),
        "uninstall", "-y", "pyalps-cli-test",
    ])
    assert {name: (sdk / name).read_bytes() for name in COMMANDS} == original_commands
