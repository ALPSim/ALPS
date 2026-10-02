# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
"""Run installed shell commands outside the checkout, as tutorial users do."""

import importlib.metadata
import importlib.util
import math
import os
from pathlib import Path
import subprocess
import sysconfig
import xml.etree.ElementTree as ET

import pytest


@pytest.fixture
def wheel_cli(tmp_path):
    # find_spec locates the package without importing it and setting ALPS_*.
    package = Path(importlib.util.find_spec("pyalps").origin).parent
    if not (package / "bin").is_dir():
        pytest.skip("this installation does not bundle ALPS programs")
    scripts = Path(sysconfig.get_path("scripts"))
    env = os.environ.copy()
    for key in ("ALPS_XML_PATH", "ALPS_BIN_PATH", "ALPS_ROOT", "PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    env["PATH"] = str(scripts) + os.pathsep + os.defpath

    def run(command, *args, **kwargs):
        assert (scripts / command).is_file(), f"pip did not install {command}"
        result = subprocess.run(
            [command, *args], cwd=tmp_path, env=env,
            text=True, capture_output=True, timeout=60, **kwargs,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return result

    return run, package, scripts


def test_entry_points_cover_bundled_programs(wheel_cli):
    _, package, scripts = wheel_cli
    entries = {
        entry.name: entry for entry in importlib.metadata.distribution("pyalps").entry_points
        if entry.group == "console_scripts"
    }
    programs = {p.name for p in (package / "bin").iterdir() if p.is_file()}
    assert set(entries) == programs
    for name, entry in entries.items():
        assert (scripts / name).is_file()
        assert callable(entry.load())
        assert entry.module == "pyalps_cli"


@pytest.mark.parametrize("via_python", [False, True])
def test_parameter2xml_then_spinmc(wheel_cli, tmp_path, via_python):
    run, _, _ = wheel_cli
    parameters = tmp_path / "simulation input"
    parameters.write_text('''LATTICE="square lattice"
MODEL="Ising"
L=4
J=1
T=2
UPDATE="cluster"
THERMALIZATION=8
SWEEPS=32
SEED=42
{}
''')
    run("parameter2xml", parameters.name)
    job = parameters.name + ".in.xml"
    assert ET.parse(tmp_path / job).getroot().tag == "JOB"
    if via_python:
        # runApplication now encounters pip's launcher on PATH as well.
        run("python", "-c",
            "import pyalps, sys; sys.exit(pyalps.runApplication("
            "'spinmc', sys.argv[1], Tmin=1, writexml=True)[0])", job)
    else:
        run("spinmc", "--Tmin", "1", "--write-xml", job)
    output = tmp_path / (parameters.name + ".task1.out.xml")
    root = ET.parse(output).getroot()
    means = root.findall(".//SCALAR_AVERAGE[@name='Energy']/MEAN")
    assert means
    assert all(math.isfinite(float(mean.text)) for mean in means)


def test_printgraph_uses_bundled_lattice_library(wheel_cli):
    run, _, _ = wheel_cli
    result = run("printgraph", input='LATTICE="chain lattice"\nL=4\n')
    graph = ET.fromstring(result.stdout)
    assert graph.tag == "GRAPH"
    assert len(graph.findall("VERTEX")) == 4
