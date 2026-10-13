"""Cover embedded scheduler startup, parallel stdin startup, and MPI ownership."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("generate_citations", root / ".github/scripts/generate_citations.py")
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)
_, policy, references, framework = generator.load_catalog(root)
binary, *launcher = sys.argv[1:]
for mode, component in (("single", "spinmc"), ("parapack", "looper"), ("owned-query", "spinmc"), ("single-query", "spinmc")):
    arguments = ["--citations"] if mode in ("owned-query", "single-query") else ["--mpi"] if launcher and mode == "parapack" else []
    render = generator.detailed_notice if mode in ("owned-query", "single-query") else generator.notice
    expected = render(policy, references, framework, component)
    for disabled in (False, True):
        env = os.environ.copy()
        env.pop("ALPS_NO_CITATIONS", None)
        if disabled:
            env["ALPS_NO_CITATIONS"] = "1"
        with tempfile.TemporaryDirectory(prefix="alps-citation-runtime-") as cwd:
            result = subprocess.run(launcher + [binary, mode] + arguments, cwd=cwd, env=env,
                                    input="SEED=17; {} {}", text=True, capture_output=True, timeout=45)
        assert result.returncode == 0, result.stdout + result.stderr
        notice_output = result.stderr if mode == "parapack" else result.stdout
        count = 0 if disabled and mode in ("single", "parapack") else 1
        assert notice_output.count(expected) == count, notice_output
        assert notice_output.count("Recommended citations for ") == count, notice_output
        if mode == "parapack":
            assert "Recommended citations for " not in result.stdout, result.stdout
            assert result.stdout.startswith("[input parameters]"), result.stdout
            assert result.stdout.count("[input parameters]") == 2, result.stdout
            if launcher:
                for rank in range(2):
                    assert result.stdout.count(f"[results {rank}]") == 2, result.stdout
            else:
                assert result.stdout.count("[results]") == 2, result.stdout
            if disabled:
                assert result.stdout == numerical_stdout, result.stdout
            else:
                numerical_stdout = result.stdout
        elif mode in ("owned-query", "single-query"):
            assert result.stdout == expected, result.stdout
        if mode in ("single", "parapack"):
            assert "copyright (c)" in notice_output, notice_output
            assert "Licensed under the MIT License." in notice_output, notice_output
