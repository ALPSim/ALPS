# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
"""Exercise wheel ownership and command selection with pip, without compiling ALPS."""

import ast
import os
from pathlib import Path
import runpy
import subprocess
import sys
import venv

import pytest


PROJECT = Path(__file__).resolve().parents[2] / "python/pyalps"
CLI = runpy.run_path(str(PROJECT / "src/pyalps_cli/__init__.py"))
COMMANDS = (*CLI["NATIVE_PROGRAMS"], *CLI["EXPORTERS"])


def _run(command, **kwargs):
    return subprocess.run(command, text=True, capture_output=True, timeout=60, **kwargs)


def _sdk(directory, label):
    directory.mkdir(parents=True, exist_ok=True)
    for name in COMMANDS:
        executable = directory / name
        executable.write_text(f'#!/bin/sh\nprintf "%s\\n" "{label}" "$@"\n')
        executable.chmod(0o755)
    return directory


@pytest.fixture(scope="module", params=["bundled", "bindings-only", "bindings-same-prefix"])
def installed_cli(tmp_path_factory, request):
    root = tmp_path_factory.mktemp("cli-installed")
    bundled = request.param == "bundled"
    environment = root / "environment with spaces"
    venv.EnvBuilder(with_pip=False, symlinks=True).create(environment)
    python = environment / "bin/python"
    sdk = _sdk(
        python.parent if request.param == "bindings-same-prefix" else root / "SDK with spaces/bin",
        "configured SDK",
    )
    original_commands = {name: (sdk / name).read_bytes() for name in COMMANDS}
    source = root / "source"
    source.mkdir()
    # The real metadata provider decides whether pip gets any console scripts.
    # Fake native programs let us test installation without a full C++ build.
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
    payload = (
        f'install(DIRECTORY "{sdk.as_posix()}/" DESTINATION pyalps/bin USE_SOURCE_PERMISSIONS)'
        if bundled else ""
    )
    (source / "CMakeLists.txt").write_text(f'''
cmake_minimum_required(VERSION 3.22)
project(cli_test LANGUAGES NONE)
set(PYALPS_ALPS_BIN_FALLBACK "{'' if bundled else sdk.as_posix()}")
configure_file("{PROJECT.as_posix()}/src/pyalps/pyalps_config.py.in"
  "${{CMAKE_CURRENT_BINARY_DIR}}/pyalps_config.py" @ONLY)
install(FILES "${{CMAKE_CURRENT_BINARY_DIR}}/pyalps_config.py" DESTINATION pyalps)
install(FILES "${{CMAKE_CURRENT_SOURCE_DIR}}/__init__.py" DESTINATION pyalps)
install(DIRECTORY "{PROJECT.as_posix()}/src/pyalps_cli" DESTINATION .
  FILES_MATCHING PATTERN "*.py")
{payload}
''')
    (source / "__init__.py").write_text('raise AssertionError("pyalps must not be imported")\n')
    wheels = root / "wheels"
    wheels.mkdir()
    build_env = {**os.environ, "PYALPS_BUNDLE_APPLICATIONS": "ON" if bundled else "OFF"}
    built = _run(
        [sys.executable, "-c",
         "from scikit_build_core.build import build_wheel, prepare_metadata_for_build_wheel; "
         "import sys; metadata = prepare_metadata_for_build_wheel('metadata'); "
         "build_wheel(sys.argv[1], metadata_directory='metadata/' + metadata)", str(wheels)],
        cwd=source, env=build_env,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    installed = _run([
        sys.executable, "-m", "pip", "--python", str(python), "install",
        "--no-deps", "--no-index", "--no-cache-dir", str(next(wheels.glob("*.whl"))),
    ])
    assert installed.returncode == 0, installed.stdout + installed.stderr
    env = os.environ.copy()
    for key in ("ALPS_XML_PATH", "ALPS_BIN_PATH", "PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    env["PATH"] = os.pathsep.join((str(python.parent), str(sdk), os.defpath))
    yield python, sdk, env, bundled

    # Uninstall must not remove any source-installed commands, including when
    # bindings and the SDK share a prefix. Check every overlapping command.
    removed = _run([
        sys.executable, "-m", "pip", "--python", str(python),
        "uninstall", "-y", "pyalps-cli-test",
    ])
    assert removed.returncode == 0, removed.stdout + removed.stderr
    assert {name: (sdk / name).read_bytes() for name in COMMANDS} == original_commands
    if bundled:
        assert all(not (python.parent / name).exists() for name in COMMANDS)


def _python_helper(python, env, expression, **kwargs):
    # Load the production helper functions without the native/scientific stack.
    tree = ast.parse((PROJECT / "src/pyalps/tools.py").read_text())
    names = {
        "check_existence", "list2cmdline", "_execute_application", "runApplication",
        "make_list", "runDMFT", "evaluateLoop", "evaluateSpinMC",
    }
    helpers = "\n".join(
        ast.unparse(node) for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    )
    code = (
        "import os, subprocess\n"
        "from pyalps_cli import resolve_executable as _resolve_executable\n"
        "log = lambda message: None\n" + helpers + "\n" + expression + "\n"
    )
    return _run([str(python), "-c", code], env=env, **kwargs)


def test_only_bundled_wheels_install_launchers(installed_cli):
    python, sdk, env, bundled = installed_cli
    result = _run([str(python), "-c", "import importlib.metadata as m; "
                   "print(','.join(sorted(e.name for e in "
                   "m.distribution('pyalps-cli-test').entry_points)))"], env=env)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == (",".join(sorted(COMMANDS)) if bundled else "")
    # Existing SDK commands survive installation, as well as fixture uninstall.
    assert all("pyalps_cli" not in (sdk / name).read_text() for name in COMMANDS)
    result = _run(["spinmc", "argument with spaces"], env=env)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["configured SDK", "argument with spaces"]


@pytest.mark.parametrize("expression", [
    "runApplication('spinmc', 'job.in.xml')",
    "runDMFT(['job.in.xml'])",
    "evaluateLoop(['job.in.xml'])",
    "evaluateSpinMC(['job.in.xml'])",
])
def test_python_helpers_ignore_conflicting_path(installed_cli, tmp_path, expression):
    python, _, env, _ = installed_cli
    conflicting = _sdk(tmp_path / "conflicting SDK", "wrong SDK")
    env = {**env, "PATH": str(conflicting) + os.pathsep + env["PATH"]}
    result = _python_helper(python, env, expression)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines()[0] == "configured SDK"


def test_sdk_override_only_applies_to_bindings_only(installed_cli, tmp_path):
    python, _, env, bundled = installed_cli
    sdk = _sdk(tmp_path / "selected SDK", "selected SDK")
    env = {**env, "ALPS_BIN_PATH": str(sdk)}
    result = _python_helper(python, env, "runApplication('spinmc', 'job.in.xml')")
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["configured SDK" if bundled else "selected SDK", "job.in.xml"]


@pytest.mark.parametrize("relative", [False, True])
def test_explicit_executable_path_and_literal_arguments(installed_cli, tmp_path, relative):
    python, _, env, _ = installed_cli
    sdk = _sdk(tmp_path / "explicit SDK $(not-a-shell)", "explicit SDK")
    executable = sdk / "spinmc"
    if relative:
        executable = Path(".") / executable.relative_to(tmp_path)
    argument = "input $(not-a-shell); with spaces.in.xml"
    result = _python_helper(
        python, env, f"runApplication({str(executable)!r}, {argument!r})", cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["explicit SDK", argument]


def test_missing_selected_sdk_does_not_search_path(installed_cli, tmp_path):
    python, _, env, bundled = installed_cli
    if bundled:
        pytest.skip("SDK overrides do not apply to bundled wheels")
    missing = tmp_path / "missing SDK"
    result = _python_helper(
        python, {**env, "ALPS_BIN_PATH": str(missing)},
        "runApplication('spinmc', 'job.in.xml')",
    )
    assert result.returncode != 0
    assert str(missing / "spinmc") in result.stderr
    assert not result.stdout


def test_explicit_mpi_application_keeps_its_sdk_and_diagonalization_flags(installed_cli, tmp_path):
    python, _, env, _ = installed_cli
    sdk = _sdk(tmp_path / "MPI SDK", "MPI SDK")
    mpirun = tmp_path / "MPI launcher"
    mpirun.write_text('#!/bin/sh\nprintf "%s\\n" "$ALPS_BIN_PATH" "$@"\n')
    mpirun.chmod(0o755)
    result = _python_helper(
        python, {**env, "ALPS_BIN_PATH": "/conflicting/sdk"},
        f"runApplication({str(sdk / 'fulldiag')!r}, 'job.in.xml', MPI=2, mpirun={str(mpirun)!r})",
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [
        str(sdk), "-np", "2", str(sdk / "fulldiag"), "--mpi", "--Nmax", "1", "job.in.xml",
    ]
