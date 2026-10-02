# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
"""Exercise pip's real launchers with a tiny wheel, without compiling ALPS."""

import ast
import os
from pathlib import Path
import subprocess
import sys
import venv

import pytest

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


PROJECT = Path(__file__).resolve().parents[2] / "python/pyalps"


def _run(command, **kwargs):
    return subprocess.run(
        command, text=True, capture_output=True, timeout=60, **kwargs,
    )


def _sdk(directory, label):
    directory.mkdir(parents=True)
    executable = directory / "spinmc"
    executable.write_text(f'#!/bin/sh\nprintf "%s\\n" "{label}" "$ALPS_BIN_PATH" "$@"\n')
    executable.chmod(0o755)
    return directory


@pytest.fixture(scope="module")
def installed_cli(tmp_path_factory):
    root = tmp_path_factory.mktemp("cli-installed")
    sdk = _sdk(root / "SDK with spaces" / "bin", "configured SDK")
    source = root / "source"
    source.mkdir()
    # Use the production entry points, launcher source and generated config
    # template, with no native targets or scientific Python dependencies.
    project = tomllib.loads((PROJECT / "pyproject.toml").read_text())
    entries = "\n".join(f'{key} = "{value}"' for key, value in project["project"]["scripts"].items())
    (source / "pyproject.toml").write_text(
        '[build-system]\nrequires = ["scikit-build-core>=1.0"]\n'
        'build-backend = "scikit_build_core.build"\n'
        '[project]\nname = "pyalps-cli-test"\nversion = "0.0.0"\n'
        f'[project.scripts]\n{entries}\n'
        '[tool.scikit-build]\nwheel.packages = []\n'
    )
    (source / "CMakeLists.txt").write_text(f'''
cmake_minimum_required(VERSION 3.22)
project(cli_test LANGUAGES NONE)
set(PYALPS_ALPS_BIN_FALLBACK "{sdk.as_posix()}")
configure_file("{PROJECT.as_posix()}/src/pyalps/pyalps_config.py.in"
  "${{CMAKE_CURRENT_BINARY_DIR}}/pyalps_config.py" @ONLY)
install(FILES "${{CMAKE_CURRENT_BINARY_DIR}}/pyalps_config.py" DESTINATION pyalps)
install(FILES "${{CMAKE_CURRENT_SOURCE_DIR}}/__init__.py" DESTINATION pyalps)
install(DIRECTORY "{PROJECT.as_posix()}/src/pyalps_cli" DESTINATION .
  FILES_MATCHING PATTERN "*.py")
''')
    (source / "__init__.py").write_text('raise AssertionError("pyalps must not be imported")\n')
    wheels = root / "wheels"
    wheels.mkdir()
    built = _run(
        [sys.executable, "-c", "from scikit_build_core.build import build_wheel; "
         "import sys; build_wheel(sys.argv[1])", str(wheels)], cwd=source,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    environment = root / "environment with spaces"
    venv.EnvBuilder(with_pip=False, symlinks=True).create(environment)
    python = environment / "bin/python"
    installed = _run([
        sys.executable, "-m", "pip", "--python", str(python), "install",
        "--no-deps", "--no-index", "--no-cache-dir", str(next(wheels.glob("*.whl"))),
    ])
    assert installed.returncode == 0, installed.stdout + installed.stderr
    env = os.environ.copy()
    for key in ("ALPS_XML_PATH", "ALPS_BIN_PATH", "PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    env["PATH"] = os.pathsep.join((str(python.parent), str(sdk), os.defpath))
    return python, sdk, env


def test_installed_bindings_only_launcher_and_python_helper(installed_cli):
    python, sdk, env = installed_cli
    # The interpreter path itself contains spaces: pip must generate a safe
    # console script, and that script must not hide the configured SDK.
    result = _run(["spinmc", "argument with spaces"], env=env)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["configured SDK", str(sdk), "argument with spaces"]

    # Exercise the production Python command-discovery path without loading
    # native bindings or installing NumPy/SciPy for this packaging test.
    tree = ast.parse((PROJECT / "src/pyalps/tools.py").read_text())
    names = {"check_existence", "list2cmdline", "executeCommand", "runApplication"}
    helpers = "\n".join(
        ast.unparse(node) for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    )
    code = (
        "import os, platform, subprocess\n"
        "log = lambda message: None\n"
        + helpers + "\nassert runApplication('spinmc', 'job.in.xml')[0] == 0\n"
    )
    result = _run([str(python), "-c", code], env=env)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["configured SDK", str(sdk), "job.in.xml"]


def test_bindings_only_explicit_sdk_override(installed_cli, tmp_path):
    _, _, env = installed_cli
    sdk = _sdk(tmp_path / "other SDK" / "bin", "selected SDK")
    result = _run(["spinmc"], env={**env, "ALPS_BIN_PATH": str(sdk)})
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["selected SDK", str(sdk)]


def test_missing_sdk_does_not_search_path(installed_cli, tmp_path):
    _, _, env = installed_cli
    missing = tmp_path / "missing SDK"
    result = _run(["spinmc"], env={**env, "ALPS_BIN_PATH": str(missing)})
    assert result.returncode == 127
    assert str(missing / "spinmc") in result.stderr
    assert "Traceback" not in result.stderr


def test_sdk_cannot_delegate_to_its_own_launcher(installed_cli):
    python, _, env = installed_cli
    result = _run(["spinmc"], env={**env, "ALPS_BIN_PATH": str(python.parent)})
    assert result.returncode == 127
    assert "is this pyalps launcher" in result.stderr
    assert "Traceback" not in result.stderr
