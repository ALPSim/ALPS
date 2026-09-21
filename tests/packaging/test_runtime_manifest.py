"""Wheel metadata must survive relocation and identify repaired filenames exactly."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "runtime_manifest", ROOT / "python/pyalps/_build_support/runtime_manifest.py")
runtime = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runtime)


def test_manifest_records_repaired_paths_without_sdk_name_matching(tmp_path):
    package = tmp_path / "pyalps"
    (package / "lib").mkdir(parents=True)
    (package / "lib/libalps.so.3").touch()
    dependencies = tmp_path / "pyalps.libs"
    dependencies.mkdir()
    (dependencies / "libhdf5_serial-a1b2c3.so.103").touch()
    (dependencies / "libgfortran-1234abcd.so.5").touch()
    (package / "pyngshdf5_c.cpython-313-x86_64-linux-gnu.so").touch()
    (package / "lib/notes.txt").touch()
    runtime.write_manifest(tmp_path, "3.0.0", repaired=True)
    manifest = json.loads((package / "runtime.json").read_text())
    assert manifest == {
        "schema": 1, "alps_version": "3.0.0", "repaired": True,
        "libraries": [
            {"path": "lib/libalps.so.3"},
            {"path": "../pyalps.libs/libgfortran-1234abcd.so.5"},
            {"path": "../pyalps.libs/libhdf5_serial-a1b2c3.so.103"},
        ],
    }
    # Regeneration preserves SDK version but drops paths removed by repair.
    (dependencies / "libgfortran-1234abcd.so.5").unlink()
    runtime.write_manifest(tmp_path, repaired=True)
    assert len(json.loads((package / "runtime.json").read_text())["libraries"]) == 2


def test_windows_manifest_records_dlls_only(tmp_path):
    binary = tmp_path / "pyalps/bin"
    binary.mkdir(parents=True)
    for name in ("alps.dll", "hdf5.dll", "spinmc.exe", "alps.lib"):
        (binary / name).touch()
    runtime.write_manifest(tmp_path, "3.0.0")
    manifest = json.loads((binary.parent / "runtime.json").read_text())
    assert not manifest["repaired"]
    assert manifest["libraries"] == [{"path": "bin/alps.dll"}, {"path": "bin/hdf5.dll"}]
