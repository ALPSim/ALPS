"""Exercise live Python sources with separately installed native resources."""

import json
import os
from pathlib import Path
import shutil
import site
import subprocess
import sys
import venv

import pytest

SOURCE = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(not os.environ.get("ALPS_DIR"), reason="requires an installed SDK")


def test_editable_core_install(tmp_path):
    for dependency in ("nanobind", "scikit_build_core", "numpy", "scipy"):
        pytest.importorskip(dependency)
    source = tmp_path / "source"
    shutil.copytree(SOURCE / "python/pyalps", source, ignore=shutil.ignore_patterns(
        "dist", "build", "_skbuild", "__pycache__", ".pytest_cache"))
    for name in ("ALPS_VERSION.txt", "LICENSE.txt"):
        shutil.copy2(SOURCE / name, source)
    environment = tmp_path / "venv"
    venv.EnvBuilder().create(environment)
    python = environment / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    purelib = Path(subprocess.check_output([
        str(python), "-c", "import sysconfig; print(sysconfig.get_path('purelib'))",
    ], text=True).strip())
    # Reuse test dependencies without network access or changing the caller's env.
    (purelib / "test_dependencies.pth").write_text(
        "\n".join(site.getsitepackages()) + "\n", encoding="utf-8")
    build_env = dict(os.environ)
    build_env.pop("PYTHONPATH", None)
    build_env["CMAKE_GENERATOR"] = "Ninja"
    build_env["CMAKE_ARGS"] = " ".join(
        '"' + argument + '"' for argument in
        json.loads(os.environ.get("ALPS_TEST_CMAKE_ARGS", "[]")))
    subprocess.run([
        sys.executable, "-m", "pip", "--python", str(python), "install",
        "--no-deps", "--no-build-isolation", "--editable", str(source),
        "--config-settings=build-dir=" + str(tmp_path / "native"),
        "--config-settings=cmake.define.PYALPS_BUILD_SOLVERS=OFF",
        "--config-settings=cmake.define.PYALPS_BUNDLE_APPLICATIONS=OFF",
    ], check=True, env=build_env, cwd=tmp_path)
    probe = """
import json
from pathlib import Path
import sys
import pyalps
from pyalps._resources import runtime_directory
source, purelib = map(Path, sys.argv[1:])
assert Path(pyalps.__file__) == source / 'src/pyalps/__init__.py'
runtime = runtime_directory()
assert runtime == purelib / 'pyalps'
assert Path(pyalps.get_cmake_dir()) == runtime / 'cmake'
assert (runtime / 'cmake/pyalpsConfig.cmake').is_file()
assert (runtime / 'include/pyalps/export_simulation.hpp').is_file()
assert (runtime / 'xml/ALPS.xsl').is_file()
assert not (runtime / '__init__.py').exists()
assert not (runtime / 'pyalps_config.py').exists()
assert not (runtime / 'bin/spinmc').exists()
assert not (runtime / 'bin/spinmc.exe').exists()
assert not hasattr(pyalps, 'cthyb')
manifest = json.loads((runtime / 'runtime.json').read_text())
assert any('alps' in entry['path'] for entry in manifest['libraries'])
for entry in manifest['libraries']:
    assert (runtime / entry['path']).is_file()
"""
    def check(script):
        subprocess.run([str(python), "-c", script, str(source), str(purelib)],
                       check=True, env=build_env, cwd=tmp_path)
    check(probe)
    with (source / "src/pyalps/__init__.py").open("a", encoding="utf-8") as output:
        output.write("\n_editable_test_marker = 'live source'\n")
    check(probe + "\nassert pyalps._editable_test_marker == 'live source'\n")
