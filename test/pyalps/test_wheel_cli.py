# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
"""Exercise installed commands and Python application helpers."""

import importlib.metadata
import math
import os
from pathlib import Path
import subprocess
import sys
import sysconfig
import xml.etree.ElementTree as ET

import pytest
import pyalps


@pytest.mark.parametrize("sdk_override", [None, "/explicit/sdk/bin"])
def test_import_preserves_binary_selection_environment(tmp_path, sdk_override):
    env = os.environ.copy()
    for key in ("ALPS_BIN_PATH", "PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    if sdk_override is not None:
        env["ALPS_BIN_PATH"] = sdk_override
    result = subprocess.run(
        [sys.executable, "-c", "import os, pyalps; "
         f"assert os.environ.get('ALPS_BIN_PATH') == {sdk_override!r}"],
        cwd=tmp_path, env=env, text=True, capture_output=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.fixture
def wheel_cli(tmp_path):
    package = Path(pyalps.__file__).resolve().parent
    if not (package / "bin").is_dir():
        pytest.skip("this installation does not bundle ALPS programs")
    scripts = Path(sysconfig.get_path("scripts"))
    env = os.environ.copy()
    for key in ("ALPS_XML_PATH", "ALPS_BIN_PATH", "ALPS_ROOT", "PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    # Exclude system xsltproc: exporters must work with pip dependencies alone.
    env["PATH"] = str(scripts)

    def run(command, *args, **kwargs):
        assert (scripts / command).is_file(), f"pip did not install {command}"
        result = subprocess.run(
            [command, *args], cwd=tmp_path, env=env,
            text=True, capture_output=True, timeout=60, **kwargs,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return result

    return run, package, scripts


@pytest.fixture
def sdk(tmp_path, monkeypatch):
    sdk = tmp_path / "SDK $(not-a-shell) with spaces"
    sdk.mkdir()
    for name in ("spinmc", "dmft", "loop", "spinmc_evaluate", "fulldiag"):
        program = sdk / name
        program.write_text('#!/bin/sh\nprintf "%s\\n" "$0" "$ALPS_BIN_PATH" "$@"\n')
        program.chmod(0o755)
    monkeypatch.delenv("ALPS_BIN_PATH", raising=False)
    monkeypatch.setenv("PATH", "")
    return sdk


@pytest.mark.parametrize("helper, args, command, flags", [
    (pyalps.runApplication, ("spinmc", "job.in.xml"), "spinmc", []),
    (pyalps.runDMFT, (["job.in.xml"],), "dmft", []),
    (pyalps.evaluateLoop, (["job.in.xml"],), "loop", ["--evaluate"]),
    (pyalps.evaluateSpinMC, (["job.in.xml"],), "spinmc_evaluate", []),
])
def test_python_helpers_use_resolved_executable(sdk, monkeypatch, capfd, helper, args, command, flags):
    # Selection is covered by launcher tests; here check how helpers use it.
    monkeypatch.setattr(pyalps.tools, "_resolve_executable", lambda name: str(sdk / name))
    helper(*args)
    expected = [str(sdk / command), str(sdk), *flags, "job.in.xml"]
    assert "\n".join(expected) + "\n" in capfd.readouterr().out
    assert "ALPS_BIN_PATH" not in os.environ


@pytest.mark.parametrize("relative", [False, True])
def test_explicit_path_and_literal_arguments(sdk, tmp_path, monkeypatch, capfd, relative):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("ALPS_BIN_PATH", "/conflicting/sdk")
    executable = sdk / "spinmc"
    argument = "input $(not-a-shell); with spaces.in.xml"
    selected = executable.relative_to(tmp_path) if relative else executable
    status, _ = pyalps.runApplication(str(selected), argument)
    assert status == 0
    assert f"{executable}\n{sdk}\n{argument}\n" in capfd.readouterr().out
    assert os.environ["ALPS_BIN_PATH"] == "/conflicting/sdk"


def test_explicit_mpi_application_keeps_sdk_and_diagonalization_flags(sdk, tmp_path, monkeypatch, capfd):
    monkeypatch.setenv("ALPS_BIN_PATH", "/conflicting/sdk")
    mpirun = tmp_path / "MPI launcher"
    mpirun.write_text('#!/bin/sh\nprintf "%s\\n" "$ALPS_BIN_PATH" "$@"\n')
    mpirun.chmod(0o755)
    status, _ = pyalps.runApplication(str(sdk / "fulldiag"), "job.in.xml", MPI=2, mpirun=str(mpirun))
    assert status == 0
    expected = [str(sdk), "-np", "2", str(sdk / "fulldiag"), "--mpi", "--Nmax", "1", "job.in.xml"]
    assert "\n".join(expected) + "\n" in capfd.readouterr().out
    assert os.environ["ALPS_BIN_PATH"] == "/conflicting/sdk"


def test_entry_points_cover_bundled_programs(wheel_cli):
    _, package, scripts = wheel_cli
    entries = {
        entry.name: entry for entry in importlib.metadata.distribution("pyalps").entry_points
        if entry.group == "console_scripts"
    }
    programs = {p.name for p in (package / "bin").iterdir() if p.is_file()}
    assert set(entries) == programs | {"convert2text", "plot2text", "plot2gp", "plot2xmgr"}
    for name, entry in entries.items():
        assert (scripts / name).is_file()
        assert callable(entry.load())


def test_parameter2xml_then_spinmc(wheel_cli, tmp_path):
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
    run("spinmc", "--Tmin", "1", "--write-xml", job)
    output = tmp_path / (parameters.name + ".task1.out.xml")
    root = ET.parse(output).getroot()
    means = root.findall(".//SCALAR_AVERAGE[@name='Energy']/MEAN")
    assert means
    assert all(math.isfinite(float(mean.text)) for mean in means)
    # convert2xml takes the legacy run file, not its HDF5 companion.
    checkpoints = sorted(tmp_path.glob("*.run[0-9]"))
    assert checkpoints
    run("convert2xml", *[str(path) for path in checkpoints])
    converted = list(tmp_path.glob("*.run*.xml"))
    assert converted
    assert "Energy" in run("convert2text", str(converted[0])).stdout


def test_printgraph_uses_bundled_lattice_library(wheel_cli):
    run, _, _ = wheel_cli
    result = run("printgraph", input='LATTICE="chain lattice"\nL=4\n')
    graph = ET.fromstring(result.stdout)
    assert graph.tag == "GRAPH"
    assert len(graph.findall("VERTEX")) == 4


def test_fulldiag_plot_exporters(wheel_cli, tmp_path):
    run, _, _ = wheel_cli
    parameters = tmp_path / "diagonalization"
    parameters.write_text('''LATTICE="chain lattice"
MODEL="spin"
local_S=0.5
J=1
{L=2}
''')
    run("parameter2xml", parameters.name)
    # Cover Python dispatch here; spinmc above covers the shell launcher.
    run("python", "-c", "import pyalps, sys; "
        "sys.exit(pyalps.runApplication('fulldiag', sys.argv[1])[0])",
        parameters.name + ".in.xml")
    run("fulldiag_evaluate", "--T_MIN", "0.5", "--T_MAX", "2", "--DELTA_T", "0.5",
        parameters.name + ".task1.out.xml")
    plots = sorted(tmp_path.glob("*.plot.xml"))
    assert plots
    source = str(plots[0])
    rows = run("plot2text", source).stdout.splitlines()
    # The stylesheet also emits a legend label before each data series.
    values = [[float(value) for value in row.split()] for row in rows if "\t" in row]
    assert len(values) == 4
    assert all(math.isfinite(value) for row in values for value in row)
    assert "set xlabel" in run("plot2gp", source).stdout
    assert "# Grace project file" in run("plot2xmgr", source).stdout


def test_simplemc_snapshot_to_vtk(wheel_cli, tmp_path):
    run, _, _ = wheel_cli
    parameters = tmp_path / "snapshot"
    parameters.write_text('''LATTICE="square lattice"
ALGORITHM="ising"
L=4
J=1
T=2
THERMALIZATION=0
SWEEPS=8
SNAPSHOT_INTERVAL=8
SEED=42
{}
''')
    run("parameter2xml", parameters.name)
    run("simplemc", "--Tmin", "1", parameters.name + ".in.xml")
    snapshots = sorted(tmp_path.glob("*.snap"))
    assert snapshots
    run("snap2vtk", str(snapshots[0]))
    vtk = snapshots[0].with_suffix(".vtk").read_text()
    assert "POINTS 16 float" in vtk
    assert "POINT_DATA 16" in vtk


def test_maxent_from_hdf5_parameters(wheel_cli, tmp_path):
    import numpy as np

    run, _, _ = wheel_cli
    parameters = {
        "BETA": 2.0, "NDAT": 6, "NFREQ": 20, "N_ALPHA": 2,
        "ALPHA_MIN": 0.1, "ALPHA_MAX": 1.0, "MAX_IT": 2,
        "OMEGA_MAX": 4.0, "FREQUENCY_GRID": "linear", "KERNEL": "fermionic",
        "DATASPACE": "time", "TEXT_OUTPUT": 0, "VERBOSE": 0,
        "PARTICLE_HOLE_SYMMETRY": 1, "NORM": 1.0, "MAX_TIME": 1,
        "BASENAME": str(tmp_path / "spectrum"),
    }
    for index in range(6):
        parameters[f"X_{index}"] = -0.5
        parameters[f"SIGMA_{index}"] = 0.01
    inputs = pyalps.writeInputH5Files(str(tmp_path / "maxent input"), [parameters])
    run("maxent", inputs[0])
    with pyalps.hdf5.archive(str(tmp_path / "spectrum.out.h5"), "r") as result:
        spectrum = result["/spectrum/maximum"]
        assert len(spectrum) == 20
        assert np.isfinite(spectrum).all()
