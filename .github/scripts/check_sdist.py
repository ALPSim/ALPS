"""Build and import the actual release sdist using an already installed SDK."""

import argparse
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import venv


def check(archive, settings):
    # Windows can retain a DLL handle briefly after the import process exits.
    # Best-effort cleanup must not change the build/import validation result.
    with tempfile.TemporaryDirectory(prefix="alps-sdist-", ignore_cleanup_errors=sys.platform == "win32") as temporary:
        root = Path(temporary)
        with tarfile.open(archive) as source:
            source.extractall(root / "source", filter="data")
        projects = list((root / "source").glob("*/pyproject.toml"))
        if len(projects) != 1:
            raise ValueError("Expected one project in the source distribution")
        output = root / "wheels"
        command = [sys.executable, "-m", "build", "--wheel", "--no-isolation",
                   "--outdir", str(output), str(projects[0].parent)]
        for setting in settings:
            command.extend(("--config-setting", setting))
        subprocess.run(command, cwd=root, check=True)
        wheels = list(output.glob("*.whl"))
        if len(wheels) != 1:
            raise ValueError("Expected one rebuilt wheel")
        environment = root / "venv"
        # A separate interpreter proves imports come from this rebuilt wheel.
        # Runtime requirements use binary wheels and the normal pip cache.
        venv.EnvBuilder(with_pip=True).create(environment)
        python = environment / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
        subprocess.run([str(python), "-m", "pip", "install", "--only-binary=:all:",
                        str(wheels[0])], cwd=root, check=True)
        subprocess.run([str(python), "-c", (
            "from pathlib import Path; import sys, pyalps; "
            "assert Path(pyalps.__file__).is_relative_to(sys.prefix); "
            "from pyalps._ext import cthyb, ctint, maxent_c; "
            "print('Release sdist rebuilt and imported:', pyalps.__file__)"
        )], cwd=root, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--config-setting", action="append", default=[])
    args = parser.parse_args()
    check(args.archive.resolve(), args.config_setting)
