# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
"""Exercise launchers without any ALPS, NumPy, or SciPy installation."""

import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys

import pytest


@pytest.fixture
def launcher(tmp_path):
    root = tmp_path / "install with spaces"
    root.mkdir()
    source = Path(__file__).resolve().parents[2] / "python/pyalps/src/pyalps_cli"
    shutil.copytree(source, root / "pyalps_cli", ignore=shutil.ignore_patterns("__pycache__"))
    package = root / "pyalps"
    (package / "bin").mkdir(parents=True)
    (package / "xml").mkdir()
    (package / "__init__.py").write_text('raise AssertionError("pyalps must not be imported")\n')
    executable = package / "bin/spinmc"
    env = os.environ.copy()
    for key in ("ALPS_XML_PATH", "ALPS_BIN_PATH", "PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    # Include a conflicting command; an absolute bundled path must win.
    other = tmp_path / "other"
    other.mkdir()
    (other / "spinmc").write_text("#!/bin/sh\nexit 99\n")
    (other / "spinmc").chmod(0o755)
    env["PATH"] = str(other)

    def run(body=None, args=(), **kwargs):
        if body is not None:
            executable.write_text(f"#!{sys.executable}\n{body}\n")
            executable.chmod(0o755)
        return subprocess.run(
            [sys.executable, "-S", "-c", "import pyalps_cli; raise SystemExit(pyalps_cli.spinmc())", *args],
            cwd=root, env=env, text=True, capture_output=True, timeout=10, **kwargs,
        )

    return run, package, env


def test_arguments_streams_exit_status_and_default_resources(launcher):
    run, package, _ = launcher
    result = run('''import json, os, sys
print(json.dumps([sys.argv[1:], sys.stdin.read(), os.environ["ALPS_XML_PATH"], os.environ["ALPS_BIN_PATH"]]))
print("native stderr", file=sys.stderr)
sys.exit(23)''', args=["a b", "", "$(not-a-shell)", "--flag"], input="input data\n")
    assert result.returncode == 23
    assert json.loads(result.stdout) == [
        ["a b", "", "$(not-a-shell)", "--flag"], "input data\n",
        str(package / "xml"), str(package / "bin"),
    ]
    assert result.stderr == "native stderr\n"


def test_explicit_resource_overrides(launcher):
    run, _, env = launcher
    env.update(ALPS_XML_PATH="/custom/xml", ALPS_BIN_PATH="/custom/bin")
    result = run('import os; print(os.environ["ALPS_XML_PATH"]); print(os.environ["ALPS_BIN_PATH"])')
    assert result.returncode == 0
    assert result.stdout.splitlines() == ["/custom/xml", "/custom/bin"]


def test_missing_binary_is_actionable_and_does_not_search_path(launcher):
    run, _, _ = launcher
    result = run()
    assert result.returncode == 127
    assert "PYALPS_BUNDLE_APPLICATIONS=ON" in result.stderr
    assert "Traceback" not in result.stderr


def test_missing_bundled_binary_never_uses_sdk_override(launcher):
    run, package, env = launcher
    (package / "pyalps_config.py").write_text('ALPS_BIN_INSTALL_DIR = ""\n')
    env["ALPS_BIN_PATH"] = env["PATH"]  # contains a working, conflicting spinmc
    result = run()
    assert result.returncode == 127
    assert "PYALPS_BUNDLE_APPLICATIONS=ON" in result.stderr


def test_nonexecutable_binary_is_reported(launcher):
    run, package, _ = launcher
    (package / "bin/spinmc").write_text("not executable")
    result = run()
    assert result.returncode == 126
    assert "cannot execute" in result.stderr


@pytest.mark.skipif(os.name != "posix", reason="POSIX exec signal semantics")
def test_native_signal_is_preserved(launcher):
    run, _, _ = launcher
    result = run('import os, signal; os.kill(os.getpid(), signal.SIGTERM)')
    assert result.returncode == -signal.SIGTERM


@pytest.mark.skipif(os.name != "posix", reason="POSIX exec signal semantics")
@pytest.mark.parametrize("signal_name", ["SIGPIPE", "SIGXFSZ"])
def test_python_ignored_signals_have_native_defaults(launcher, signal_name):
    native_signal = getattr(signal, signal_name, None)
    if native_signal is None:
        pytest.skip(f"{signal_name} is unavailable on this platform")
    run, package, _ = launcher
    executable = package / "bin/spinmc"
    # Use a native shell: another Python interpreter would ignore these
    # signals again. Disable core dumps before sending SIGXFSZ.
    executable.write_text(
        f"#!/bin/sh\nulimit -c 0\nkill -{native_signal} $$\nexit 99\n"
    )
    executable.chmod(0o755)
    direct = subprocess.run(
        [str(executable)], capture_output=True, text=True, timeout=10,
    )
    assert direct.returncode == -native_signal
    result = run()
    assert result.returncode == direct.returncode
