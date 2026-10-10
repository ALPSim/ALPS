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

    def run(body=None, args=(), program="spinmc", isolated=True, **kwargs):
        if body is not None:
            executable.write_text(f"#!{sys.executable}\n{body}\n")
            executable.chmod(0o755)
        return subprocess.run(
            [sys.executable, *(["-S"] if isolated else []), "-c",
             f"import sys, pyalps_cli; sys.argv[0] = {program!r}; "
             "raise SystemExit(pyalps_cli.main())", *args],
            cwd=root, env=env, text=True, capture_output=True, timeout=10, **kwargs,
        )

    return run, package, env


def test_exporter_without_pyalps_or_external_programs(launcher, tmp_path):
    run, package, _ = launcher
    stylesheets = Path(__file__).resolve().parents[2] / "src/alps/resources"
    shutil.copytree(stylesheets, package / "xml", dirs_exist_ok=True)
    source = tmp_path / "input with spaces.xml"
    source.write_text('<!DOCTYPE plot SYSTEM "https://invalid.example/unavailable.dtd">\n'
                      '<plot><set><point><x>1</x><y>2</y></point></set></plot>')
    result = run(program="plot2text", isolated=False, args=[str(source)])
    assert result.returncode == 0, result.stderr
    assert "1\t2" in result.stdout
    assert not result.stderr


def test_exporter_stdin_and_stylesheet_override(launcher, tmp_path):
    run, _, env = launcher
    custom = tmp_path / "custom stylesheets"
    custom.mkdir()
    (custom / "plot2text.xsl").write_text('''<xsl:stylesheet version="1.0"
      xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
      <xsl:output method="text"/><xsl:template match="/">override</xsl:template>
      </xsl:stylesheet>''')
    env["ALPS_XML_PATH"] = str(custom)
    result = run(program="plot2text", isolated=False, args=["-"], input="<plot/>")
    assert result.returncode == 0, result.stderr
    assert result.stdout == "override"


@pytest.mark.parametrize("xml", [None, "<broken>"])
def test_exporter_reports_bad_input(launcher, tmp_path, xml):
    run, _, _ = launcher
    source = tmp_path / "input.xml"
    if xml is not None:
        source.write_text(xml)
    result = run(program="plot2text", isolated=False, args=[str(source)])
    assert result.returncode == 1
    assert "plot2text:" in result.stderr
    assert "Traceback" not in result.stderr


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


def test_xml_override_does_not_select_another_binary_installation(launcher):
    run, package, env = launcher
    env.update(ALPS_XML_PATH="/custom/xml", ALPS_BIN_PATH="/custom/bin")
    result = run('import os; print(os.environ["ALPS_XML_PATH"]); print(os.environ["ALPS_BIN_PATH"])')
    assert result.returncode == 0
    assert result.stdout.splitlines() == ["/custom/xml", str(package / "bin")]


def test_missing_bundled_binary_does_not_use_path_or_sdk(launcher):
    run, package, env = launcher
    (package / "pyalps_config.py").write_text('ALPS_BIN_INSTALL_DIR = ""\n')
    env["ALPS_BIN_PATH"] = env["PATH"]  # contains a working, conflicting spinmc
    result = run()
    assert result.returncode == 127
    assert "PYALPS_BUNDLE_APPLICATIONS=ON" in result.stderr
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize("selection", ["configured", "override", "missing"])
def test_bindings_only_sdk_selection(launcher, tmp_path, selection):
    run, package, env = launcher
    sdk = tmp_path / "configured SDK"
    override = tmp_path / "selected SDK"
    for directory in (sdk, override):
        directory.mkdir()
        executable = directory / "spinmc"
        executable.write_text('#!/bin/sh\nprintf "%s\\n" "$0"\n')
        executable.chmod(0o755)
    (package / "pyalps_config.py").write_text(f"ALPS_BIN_INSTALL_DIR = {str(sdk)!r}\n")
    selected = sdk if selection == "configured" else override
    if selection == "missing":
        selected = tmp_path / "missing SDK"
    if selection != "configured":
        env["ALPS_BIN_PATH"] = str(selected)
    result = run()
    if selection == "missing":
        assert result.returncode == 127
        assert str(selected / "spinmc") in result.stderr
        assert not result.stdout
    else:
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == str(selected / "spinmc")


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
