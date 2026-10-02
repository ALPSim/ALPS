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
    (package / "lib/libalps_osiris-a1b2c3.so.3").touch()
    dependencies = tmp_path / "pyalps.libs"
    dependencies.mkdir()
    (dependencies / "libhdf5_serial-a1b2c3.so.103").touch()
    (dependencies / "libgfortran-1234abcd.so.5").touch()
    (package / "pyngshdf5_c.cpython-313-x86_64-linux-gnu.so").touch()
    (package / "lib/notes.txt").touch()
    runtime.write_manifest(tmp_path, "3.0.0", repaired=True)
    manifest = json.loads((package / "runtime.json").read_text())
    assert manifest == {
        "schema": 1, "nanobind": {}, "alps_version": "3.0.0", "repaired": True,
        "libraries": [
            {"path": "lib/libalps.so.3"},
            {"path": "lib/libalps_osiris-a1b2c3.so.3"},
            {"path": "../pyalps.libs/libgfortran-1234abcd.so.5"},
            {"path": "../pyalps.libs/libhdf5_serial-a1b2c3.so.103"},
        ],
    }
    # Regeneration preserves SDK version but drops paths removed by repair.
    (dependencies / "libgfortran-1234abcd.so.5").unlink()
    runtime.write_manifest(tmp_path, repaired=True)
    assert len(json.loads((package / "runtime.json").read_text())["libraries"]) == 3


def test_windows_manifest_records_dlls_only(tmp_path):
    binary = tmp_path / "pyalps/bin"
    binary.mkdir(parents=True)
    for name in ("alps.dll", "alps_osiris.dll", "hdf5.dll", "spinmc.exe", "alps.lib"):
        (binary / name).touch()
    runtime.write_manifest(tmp_path, "3.0.0")
    manifest = json.loads((binary.parent / "runtime.json").read_text())
    assert not manifest["repaired"]
    assert manifest["libraries"] == [
        {"path": "bin/alps.dll"}, {"path": "bin/alps_osiris.dll"}, {"path": "bin/hdf5.dll"}]


def test_repair_preserves_nanobind_abi(tmp_path):
    (tmp_path / "pyalps").mkdir()
    runtime.write_manifest(tmp_path, "3.0.0", nanobind_version="2.15.0", nanobind_abi="21")
    runtime.write_manifest(tmp_path, repaired=True)
    metadata = json.loads((tmp_path / "pyalps/runtime.json").read_text())
    assert metadata["nanobind"] == {"version": "2.15.0", "internals_abi": "21"}


def test_macos_component_dependencies_use_packaged_libraries(tmp_path, monkeypatch):
    package = tmp_path / "pyalps"
    (package / "lib").mkdir(parents=True)
    (package / "_ext").mkdir()
    core = package / "lib/libalps.3.dylib"
    utilities = package / "lib/libalps_utilities.3.dylib"
    hdf5 = package / "lib/libalps_hdf5.3.dylib"
    params = package / "lib/libalps_params.3.dylib"
    osiris = package / "lib/libalps_osiris.3.dylib"
    xml = package / "lib/libalps_xml.3.dylib"
    cli = package / "lib/libalps_cli.3.dylib"
    extension = package / "_ext/example.so"
    for binary in (core, params, hdf5, utilities, osiris, xml, cli, extension):
        binary.touch()
    monkeypatch.setattr(runtime.sys, "platform", "darwin")
    monkeypatch.setattr(runtime.subprocess, "check_output",
                        lambda command, **kwargs: f"{command[-1]}:\n@rpath/{Path(command[-1]).name}\n")
    commands = []
    monkeypatch.setattr(runtime.subprocess, "run", lambda command, **kwargs: commands.append(command))
    runtime.write_manifest(tmp_path, "3.0.0")
    rewrites = {command[-1]: command for command in commands if command[0] == "install_name_tool"}
    assert set(rewrites) == {str(core), str(params), str(hdf5), str(utilities), str(osiris), str(xml), str(cli), str(extension)}
    assert "@loader_path/libalps_utilities.3.dylib" in rewrites[str(core)]
    assert "@loader_path/../lib/libalps_utilities.3.dylib" in rewrites[str(extension)]
    assert "@loader_path/libalps_hdf5.3.dylib" in rewrites[str(core)]
    assert "@loader_path/libalps_utilities.3.dylib" in rewrites[str(hdf5)]
    assert "@loader_path/libalps_params.3.dylib" in rewrites[str(core)]
    assert "@loader_path/libalps_hdf5.3.dylib" in rewrites[str(params)]
    assert "@loader_path/../lib/libalps_params.3.dylib" in rewrites[str(extension)]
    manifest = json.loads((package / "runtime.json").read_text())
    for name in ("osiris", "xml", "cli"):
        filename = f"libalps_{name}.3.dylib"
        for binary, destination in (
            (core, f"@loader_path/{filename}"),
            (extension, f"@loader_path/../lib/{filename}"),
        ):
            command = rewrites[str(binary)]
            index = command.index(f"@rpath/{filename}")
            assert command[index - 1:index + 2] == [
                "-change", f"@rpath/{filename}", destination]
        assert {"path": f"lib/{filename}",
                "install_name": f"@rpath/{filename}"} in manifest["libraries"]
    signed = {command[-1] for command in commands if command[0] == "codesign"}
    assert signed == set(rewrites)
