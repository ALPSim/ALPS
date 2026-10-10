"""Install legacy XML wrappers without requiring the compiled SDK dependencies."""
import os
from pathlib import Path
import shutil
import stat
import subprocess

import pytest

from test_xml import PLOT, simulation


SOURCE = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(os.name == "nt", reason="Unix shell wrappers")


@pytest.fixture(scope="module", params=[("bin", "share"),
                                       ("tools/bin", "resources with spaces")])
def legacy_install(request, tmp_path_factory):
    assert shutil.which("xsltproc"), "Install xsltproc to run the XML tool tests"
    directory = tmp_path_factory.mktemp("legacy xml install")
    project = directory / "project"
    project.mkdir()
    # Exercise the production wrapper configure/install rules in isolation.
    (project / "CMakeLists.txt").write_text(
        'cmake_minimum_required(VERSION 3.27)\n'
        'project(xml_wrappers LANGUAGES C)\n'
        'include(GNUInstallDirs)\n'
        f'add_subdirectory("{SOURCE / "src/tools/xml"}" xml)\n'
        f'install(DIRECTORY "{SOURCE / "src/alps/resources"}/"\n'
        '  DESTINATION "${CMAKE_INSTALL_DATADIR}/alps/xml" COMPONENT xml\n'
        '  FILES_MATCHING PATTERN "*.xml" PATTERN "*.xsl")\n'
    )
    bindir, datadir = request.param
    build = directory / "build"
    original = directory / "original prefix"
    subprocess.run([
        "cmake", "-S", str(project), "-B", str(build),
        f"-DCMAKE_INSTALL_PREFIX={original}",
        f"-DCMAKE_INSTALL_BINDIR={bindir}", f"-DCMAKE_INSTALL_DATADIR={datadir}",
    ], check=True, capture_output=True, text=True)
    for component in ("xml", "tools"):
        subprocess.run(["cmake", "--install", str(build), "--component", component],
                       check=True, capture_output=True, text=True)
    relocated = directory / "relocated prefix with spaces"
    shutil.move(original, relocated)
    return relocated / bindir, relocated / datadir / "alps/xml"


def run(installation, command, *args, **kwargs):
    return subprocess.run([str(installation[0] / command), *map(str, args)],
                          capture_output=True, text=True, **kwargs)


@pytest.mark.parametrize("command, marker", [
    ("plot2text", "1\t3\t0.1"), ("plot2html", "<html"),
    ("plot2gp", 'set title "Energy"'), ("plot2xmgr", "@"),
    ("convert2text", "Energy"), ("convert2html", "<html"),
])
def test_legacy_conversion_after_relocation(legacy_install, tmp_path, command, marker):
    source = tmp_path / "input 'quoted'; data.xml"
    # Legacy convert commands replace the two-line declaration/style header.
    source.write_text('<?xml version="1.0"?>\n'
                      '<?xml-stylesheet type="text/xsl" href="ALPS.xsl"?>\n'
                      + (PLOT if command.startswith("plot") else simulation(1, 3)))
    result = run(legacy_install, command, source)
    assert result.returncode == 0, result.stderr
    assert marker in result.stdout


def test_local_stylesheet_preserves_files_and_permissions(legacy_install, tmp_path):
    sources = [tmp_path / "input 'quoted'; data.xml", tmp_path / "-second.xml"]
    for source in sources:
        source.write_text(simulation(1, 3))
        source.chmod(0o640)
    result = run(legacy_install, "use_local_stylesheet",
                 *(source.name for source in sources), cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    for source in sources:
        assert 'href="ALPS.xsl"' in source.read_text()
        assert '<MEAN>3</MEAN>' in source.read_text()
        assert stat.S_IMODE(source.stat().st_mode) == 0o640
    assert (tmp_path / "ALPS.xsl").read_bytes() == (legacy_install[1] / "ALPS.xsl").read_bytes()
    assert set(tmp_path.iterdir()) == {*sources, tmp_path / "ALPS.xsl"}


@pytest.mark.parametrize("failure", ["malformed", "partial_output", "missing_stylesheet"])
def test_failed_local_transformation_keeps_input(legacy_install, tmp_path, failure):
    source = tmp_path / "result.xml"
    original = b"<broken" if failure == "malformed" else simulation(1, 3).encode()
    source.write_bytes(original)
    source.chmod(0o640)
    env = dict(os.environ)
    installation = legacy_install
    if failure == "partial_output":
        transformer = tmp_path / "failing transformer"
        transformer.write_text('#!/bin/sh\nif [ "$#" -eq 0 ]; then echo xsltproc; exit 0; fi\n'
                               'echo partial output\nexit 7\n')
        transformer.chmod(0o755)
        env["ALPS_XSLT_TRANSFORMER"] = str(transformer)
    elif failure == "missing_stylesheet":
        # Copy the install so the module fixture remains intact for later tests.
        isolated = tmp_path / "isolated"
        shutil.copytree(legacy_install[0], isolated / "bin")
        # The configured relative layout must be kept when removing a resource.
        relative = os.path.relpath(legacy_install[1], legacy_install[0])
        resources = (isolated / "bin" / relative).resolve()
        shutil.copytree(legacy_install[1], resources)
        (resources / "changestylesheet.xsl").unlink()
        installation = (isolated / "bin", resources)
    before = set(tmp_path.iterdir())
    result = run(installation, "use_local_stylesheet", source.name, cwd=tmp_path, env=env)
    assert result.returncode != 0
    assert source.read_bytes() == original
    assert stat.S_IMODE(source.stat().st_mode) == 0o640
    assert set(tmp_path.iterdir()) == before | {tmp_path / "ALPS.xsl"}


@pytest.mark.parametrize("command", ["convert2html", "convert2text"])
def test_conversion_failure_is_reported_and_cleaned(legacy_install, tmp_path, command):
    source = tmp_path / "broken.xml"
    source.write_text('<?xml version="1.0"?>\n<!-- header -->\n<broken')
    temporary = tmp_path / "temporary"
    temporary.mkdir()
    result = run(legacy_install, command, source,
                 env={**os.environ, "TMPDIR": str(temporary)})
    assert result.returncode != 0
    assert list(temporary.iterdir()) == []
