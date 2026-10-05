"""Exercise the information CLI contract without simulation files, serial or MPI."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("generate_citations", root / "script/generate_citations.py")
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)
_, policy, references, framework = generator.load_catalog(root)
binary, component, *launcher = sys.argv[1:]
expected = generator.detailed_notice(policy, references, framework, component)
for disabled in (False, True):
    env = os.environ.copy()
    env.pop("ALPS_NO_CITATIONS", None)
    if disabled:
        env["ALPS_NO_CITATIONS"] = "1"
    for argument in ("--help", "-h", "--license", "-l", "--citations"):
        for mpi_flag in ([], ["--mpi"]):
            with tempfile.TemporaryDirectory(prefix="alps-citation-cli-") as cwd:
                command = launcher + [binary, argument] + mpi_flag
                result = subprocess.run(command, cwd=cwd, env=env, input="", text=True, capture_output=True, timeout=45)
                assert not list(Path(cwd).iterdir()), f"Information query created files: {command}"
            assert result.returncode == 0, f"{command}\n{result.stdout}\n{result.stderr}"
            if argument == "--citations":
                assert result.stdout == expected, f"{command}: expected exactly the catalog notice\n{result.stdout}"
            else:
                assert "Recommended citations for " not in result.stdout, command
                if argument in ("--help", "-h"):
                    assert result.stdout.count("Information options (no simulation input required):") == 1, command
                else:
                    assert result.stdout.count("Licensed under the MIT License") == 1, command
# Ambiguous queries must fail before opening input or printing citations.
for arguments in (["--help", "--citations"], ["--citations", "missing.in.h5"]):
    with tempfile.TemporaryDirectory(prefix="alps-citation-cli-") as cwd:
        result = subprocess.run(launcher + [binary] + arguments, cwd=cwd, input="", text=True,
                                capture_output=True, timeout=45)
        assert not list(Path(cwd).iterdir())
    assert result.returncode != 0, arguments
    assert "Recommended citations for " not in result.stdout, arguments
    assert "without calculation arguments" in result.stderr, result.stderr
