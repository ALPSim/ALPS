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
    (package / "runtime.json").write_text(json.dumps({"sdk_bin": ""}))
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


RESOURCES = Path(__file__).resolve().parents[2] / "src/alps/resources"
PLOT = '<plot><set><point><x>1</x><y>2</y><dy>0.5</dy></point></set></plot>'
SPECIFICATION = '''<plot name="Energy"><legend show="false"/>
<xaxis label="Temperature" type="PARAMETER" name="T"/>
<yaxis label="Energy" type="SCALAR_AVERAGE" name="Energy"/></plot>'''


def simulation(x, y):
    return f'''<SIMULATION><PARAMETERS><PARAMETER name="T">{x}</PARAMETER></PARAMETERS>
<AVERAGES><SCALAR_AVERAGE name="Energy"><MEAN>{y}</MEAN>
<ERROR>0.1</ERROR></SCALAR_AVERAGE></AVERAGES></SIMULATION>'''


def test_exporter_without_pyalps_or_external_programs(launcher, tmp_path):
    run, package, _ = launcher
    shutil.copytree(RESOURCES, package / "xml", dirs_exist_ok=True)
    source = tmp_path / "input with spaces.xml"
    source.write_text('<!DOCTYPE plot SYSTEM "https://invalid.example/unavailable.dtd">\n'
                      '<plot><set><point><x>1</x><y>2</y></point></set></plot>')
    result = run(program="alps-xml", isolated=False, args=["plot", "text", str(source)])
    assert result.returncode == 0, result.stderr
    assert "1\t2" in result.stdout
    assert not result.stderr


@pytest.mark.parametrize("program, command", [
    ("plot2text", ["plot", "text"]), ("plot2gp", ["plot", "gnuplot"]),
    ("plot2xmgr", ["plot", "grace"]), ("convert2text", ["convert", "text"]),
])
def test_exporter_alias_matches_alps_xml(launcher, tmp_path, program, command):
    run, package, _ = launcher
    shutil.copytree(RESOURCES, package / "xml", dirs_exist_ok=True)
    source = tmp_path / "input.xml"
    source.write_text(PLOT if command[0] == "plot" else simulation(1, 3))
    expected = run(program="alps-xml", isolated=False, args=[*command, str(source)])
    assert expected.returncode == 0, expected.stderr
    result = run(program=program, isolated=False, args=[str(source)])
    assert result.returncode == 0, result.stderr
    assert result.stdout == expected.stdout
    assert not result.stderr


@pytest.mark.parametrize("format, marker", [
    ("text", "1\t2\t0.5"), ("html", "<html"), ("gnuplot", "plot '-'"),
    ("grace", "# Grace project file"), ("matplotlib", "errorbar("),
])
def test_alps_xml_plot_formats(launcher, tmp_path, format, marker):
    run, package, _ = launcher
    shutil.copytree(RESOURCES, package / "xml", dirs_exist_ok=True)
    source = tmp_path / "plot.xml"
    source.write_text(PLOT)
    output = tmp_path / "result 'quoted'; name"
    result = run(program="alps-xml", isolated=False,
                 args=["plot", format, str(source), "--output", str(output)])
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert marker in output.read_text()


@pytest.mark.parametrize("archive", [False, True])
def test_alps_xml_extract_simulations_and_archives(launcher, tmp_path, archive):
    run, package, _ = launcher
    shutil.copytree(RESOURCES, package / "xml", dirs_exist_ok=True)
    specification = tmp_path / "plot specification.xml"
    specification.write_text(SPECIFICATION)
    inputs = [tmp_path / "simulation one.xml", tmp_path / "simulation two.xml"]
    for path, x, y in zip(inputs, (2, 1), (4, 3)):
        content = simulation(x, y)
        path.write_text(f"<ARCHIVE>{content}</ARCHIVE>" if archive else content)
    result = run(program="alps-xml", isolated=False,
                 args=["extract", "text", str(specification), *map(str, inputs)])
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["1\t3\t0.1", "2\t4\t0.1"]


@pytest.mark.parametrize("arguments", [
    ["plot", "text", "missing.xml"],
    ["extract", "text", "missing spec.xml", "missing.xml"],
    ["convert", "pdf", "input.xml"],
])
def test_alps_xml_reports_errors(launcher, tmp_path, arguments):
    run, package, _ = launcher
    shutil.copytree(RESOURCES, package / "xml", dirs_exist_ok=True)
    output = tmp_path / "existing.txt"
    output.write_text("previous result")
    result = run(program="alps-xml", isolated=False, args=[*arguments, "-o", str(output)])
    assert result.returncode != 0
    assert "alps-xml" in result.stderr
    assert "Traceback" not in result.stderr
    assert output.read_text() == "previous result"


@pytest.mark.skipif(shutil.which("xsltproc") is None, reason="requires xsltproc")
@pytest.mark.parametrize("arguments", [
    ["plot", "text"], ["plot", "html"], ["plot", "gnuplot"], ["plot", "grace"],
    ["plot", "matplotlib"], ["convert", "text"], ["convert", "html"],
    ["extract", "text"], ["extract", "html"], ["extract", "gnuplot"], ["extract", "grace"],
])
def test_alps_xml_matches_sdk_command(launcher, tmp_path, arguments):
    run, package, _ = launcher
    shutil.copytree(RESOURCES, package / "xml", dirs_exist_ok=True)
    sdk = tmp_path / "sdk/bin/alps-xml"
    sdk.parent.mkdir(parents=True)
    template = (RESOURCES.parents[2] / "src/tools/xml/alps-xml.in").read_text()
    sdk.write_text(template.replace("@ALPS_XML_FROM_BINDIR@", os.path.relpath(RESOURCES, sdk.parent)))
    if arguments[0] == "extract":
        specification = tmp_path / "specification.xml"
        specification.write_text(SPECIFICATION)
        inputs = [tmp_path / "one.xml", tmp_path / "two.xml"]
        for path, x, y in zip(inputs, (2, 1), (4, 3)):
            path.write_text(simulation(x, y))
        files = [specification, *inputs]
    else:
        source = tmp_path / "input.xml"
        source.write_text(PLOT if arguments[0] == "plot" else simulation(1, 3))
        files = [source]
    expected = subprocess.run([sys.executable, str(sdk), *arguments, *map(str, files)],
                              capture_output=True, text=True, timeout=30)
    assert expected.returncode == 0, expected.stderr
    result = run(program="alps-xml", isolated=False, args=[*arguments, *map(str, files)])
    assert result.returncode == 0, result.stderr
    assert result.stdout == expected.stdout


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
    (package / "runtime.json").write_text(json.dumps({"sdk_bin": ""}))
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
    (package / "runtime.json").write_text(json.dumps({"sdk_bin": str(sdk)}))
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
