"""Record wheel runtime paths after CMake installation or wheel repair.

The consumer never derives names from the SDK: auditwheel/delocate may rename
every dependency. Repacking with wheel also regenerates the required RECORD.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def write_manifest(tree: Path, version: str | None = None, *, repaired=False):
    package = tree / "pyalps"
    manifest = package / "runtime.json"
    if version is None:
        version = json.loads(manifest.read_text(encoding="utf-8"))["alps_version"]
    libraries = []
    for relative in ("lib", "bin", ".dylibs", "../pyalps.libs"):
        directory = package / relative
        if not directory.is_dir():
            continue
        for library in sorted(directory.iterdir()):
            if not library.is_file() or not (
                ".so" in library.suffixes or library.suffix.lower() in {".dll", ".dylib"}
            ):
                continue
            entry = {"path": f"{relative}/{library.name}"}
            if library.suffix == ".dylib":
                if repaired:
                    # New consumers must be able to link these libraries directly.
                    # Existing wheel references remain relative to their loaders.
                    subprocess.run([
                        "install_name_tool", "-id", f"@rpath/{library.name}", str(library)
                    ], check=True)
                    subprocess.run(["codesign", "--force", "--sign", "-", str(library)], check=True)
                output = subprocess.check_output(["otool", "-D", str(library)], text=True)
                names = {line.strip() for line in output.splitlines() if line and not line.endswith(":")}
                if len(names) != 1:
                    raise ValueError(f"Expected one install name for {library}: {names}")
                entry["install_name"] = names.pop()
            libraries.append(entry)
    manifest.write_text(json.dumps({
        "schema": 1,
        "alps_version": version,
        "repaired": repaired,
        "libraries": libraries,
    }, indent=2) + "\n", encoding="utf-8")


def finalize_wheel(wheel: Path, destination: Path):
    with tempfile.TemporaryDirectory(prefix="pyalps-runtime-") as temporary:
        subprocess.run([sys.executable, "-m", "wheel", "unpack", str(wheel), "-d", temporary], check=True)
        tree, = Path(temporary).iterdir()
        write_manifest(tree, repaired=True)
        destination.mkdir(parents=True, exist_ok=True)
        subprocess.run([sys.executable, "-m", "wheel", "pack", str(tree), "-d", str(destination)], check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--tree", type=Path)
    source.add_argument("--wheel", type=Path, nargs="+")
    parser.add_argument("--alps-version")
    parser.add_argument("--destination", type=Path)
    arguments = parser.parse_args()
    if arguments.tree:
        write_manifest(arguments.tree, arguments.alps_version)
    else:
        for wheel in arguments.wheel:
            finalize_wheel(wheel, arguments.destination or wheel.parent)
