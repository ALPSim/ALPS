"""Exercise the installed XML CLI with real transformations after relocation.

Set ALPS_XML_BUILD to a configured Unix build with applications enabled.
Only the xml install component is needed; no compilation is required.
"""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import types

import pytest


pytestmark = pytest.mark.skipif(
    os.name == "nt" or not os.environ.get("ALPS_XML_BUILD"),
    reason="requires ALPS_XML_BUILD with the Unix XML tool enabled",
)

PLOT = '''<plot name="Energy"><legend show="false"/>
<xaxis label="Temperature"/><yaxis label="Energy"/>
<set><point><x>2</x><y>4</y><dy>0.2</dy></point>
<point><x>1</x><y>3</y><dy>0.1</dy></point></set></plot>'''
SPECIFICATION = '''<plot name="Energy"><legend show="false"/>
<xaxis label="Temperature" type="PARAMETER" name="T"/>
<yaxis label="Energy" type="SCALAR_AVERAGE" name="Energy"/></plot>'''


def simulation(x, y):
    return f'''<SIMULATION><PARAMETERS><PARAMETER name="T">{x}</PARAMETER></PARAMETERS>
<AVERAGES><SCALAR_AVERAGE name="Energy"><MEAN>{y}</MEAN>
<ERROR>0.1</ERROR></SCALAR_AVERAGE></AVERAGES></SIMULATION>'''


@pytest.fixture(scope="module")
def cli(tmp_path_factory):
    assert shutil.which("xsltproc"), "Install xsltproc to run the XML tool tests"
    directory = tmp_path_factory.mktemp("xml install")
    original = directory / "original prefix"
    subprocess.run([
        "cmake", "--install", os.environ["ALPS_XML_BUILD"], "--component", "xml",
        "--prefix", str(original),
    ], check=True, capture_output=True, text=True)
    relocated = directory / "relocated prefix with spaces"
    shutil.move(original, relocated)
    assert not original.exists()
    scripts = list(relocated.rglob("alps-xml"))
    assert len(scripts) == 1, "Configure with ALPS_BUILD_APPLICATIONS=ON"
    return scripts[0]


def run(cli, *args, **kwargs):
    return subprocess.run([str(cli), *map(str, args)], capture_output=True, text=True, **kwargs)


@pytest.fixture
def plot(tmp_path):
    source = tmp_path / "plot 'quoted'; data.xml"
    source.write_text(PLOT)
    return source


@pytest.mark.parametrize("format, marker", [
    ("text", "1\t3\t0.1"), ("html", "<html"),
    ("gnuplot", 'set title "Energy"'), ("grace", "@"),
    ("matplotlib", "errorbar(data[0],data[1]"),
])
def test_plot_formats_from_relocated_install(cli, plot, format, marker):
    result = run(cli, "plot", format, plot)
    assert result.returncode == 0, result.stderr
    assert marker in result.stdout


def test_generated_matplotlib_program_runs_on_python3(cli, plot, monkeypatch):
    result = run(cli, "plot", "matplotlib", plot)
    assert result.returncode == 0, result.stderr
    calls = []
    pylab = types.ModuleType("pylab")
    for name in ("subplot", "title", "xlabel", "ylabel", "legend"):
        setattr(pylab, name, lambda *args, **kwargs: None)
    pylab.errorbar = lambda *args, **kwargs: calls.append((args, kwargs))
    monkeypatch.setitem(sys.modules, "pylab", pylab)
    monkeypatch.setitem(sys.modules, "numpy", types.ModuleType("numpy"))
    exec(compile(result.stdout, "generated-plot.py", "exec"), {})
    assert calls == [(((1, 2), (3, 4)), {"yerr": (0.1, 0.2)})]


@pytest.mark.parametrize("format, marker", [("text", "Energy"), ("html", "<html")])
def test_convert_simulation(cli, tmp_path, format, marker):
    source = tmp_path / "simulation.xml"
    source.write_text(simulation(1, 3))
    result = run(cli, "convert", format, source)
    assert result.returncode == 0, result.stderr
    assert marker in result.stdout


@pytest.mark.parametrize("archive", [False, True])
def test_extract_simulations_and_archives(cli, tmp_path, archive):
    specification = tmp_path / "plot specification.xml"
    specification.write_text(SPECIFICATION)
    inputs = [tmp_path / "simulation one.xml", tmp_path / "simulation two.xml"]
    for path, x, y in zip(inputs, (2, 1), (4, 3)):
        content = simulation(x, y)
        path.write_text(f"<ARCHIVE>{content}</ARCHIVE>" if archive else content)
    output = tmp_path / "result 'quoted'; name.txt"
    result = run(cli, "extract", "text", specification, *inputs, "--output", output)
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert output.read_text().splitlines() == ["1\t3\t0.1", "2\t4\t0.1"]


@pytest.mark.parametrize("content", ["<broken", None])
def test_transform_failure_preserves_existing_output(cli, tmp_path, content):
    source = tmp_path / "bad input.xml"
    if content is not None:
        source.write_text(content)
    output = tmp_path / "existing.txt"
    output.write_text("previous result")
    result = run(cli, "plot", "text", source, "--output", output)
    assert result.returncode != 0
    assert result.stdout == ""
    assert "alps-xml:" in result.stderr
    assert output.read_text() == "previous result"


def test_failed_extraction_cleans_temporary_files(cli, tmp_path):
    specification = tmp_path / "bad spec.xml"
    specification.write_text("<broken")
    source = tmp_path / "simulation.xml"
    source.write_text(simulation(1, 3))
    temporary = tmp_path / "temporary"
    temporary.mkdir()
    result = run(cli, "extract", "text", specification, source,
                 env={**os.environ, "TMPDIR": str(temporary)})
    assert result.returncode != 0
    assert result.stdout == ""
    assert list(temporary.iterdir()) == []


def test_repository_spin_plot_definition(cli, tmp_path):
    specification = Path(__file__).with_name("fixtures") / "spin-energy.xml"
    source = tmp_path / "simulation.xml"
    source.write_text(simulation(1, 3))
    result = run(cli, "extract", "text", specification, source)
    assert result.returncode == 0, result.stderr
    assert "1\t3\t0.1" in result.stdout


def test_missing_transformer_is_an_error(cli, plot, tmp_path):
    # Invoke Python explicitly so PATH can hide xsltproc without hiding Python.
    result = subprocess.run([sys.executable, str(cli), "plot", "text", str(plot)],
                            env={**os.environ, "PATH": str(tmp_path)},
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "xsltproc" in result.stderr
    assert result.stdout == ""
