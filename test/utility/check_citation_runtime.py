"""Cover embedded scheduler startup, parallel stdin startup, and MPI ownership."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("generate_citations", root / "script/generate_citations.py")
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)
_, policy, references, framework = generator.load_catalog(root)
binary, *launcher = sys.argv[1:]
for mode, component in (("single", "spinmc"), ("parapack", "looper"), ("owned-query", "spinmc"), ("single-query", "spinmc")):
    arguments = ["--citations"] if mode in ("owned-query", "single-query") else ["--mpi"] if launcher and mode == "parapack" else []
    with tempfile.TemporaryDirectory(prefix="alps-citation-runtime-") as cwd:
        result = subprocess.run(launcher + [binary, mode] + arguments, cwd=cwd,
                                input="SEED=17; {} {}", text=True, capture_output=True, timeout=45)
    render = generator.detailed_notice if mode in ("owned-query", "single-query") else generator.notice
    expected = render(policy, references, framework, component)
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.count(expected) == 1, result.stdout
    assert result.stdout.count("Recommended citations for ") == 1, result.stdout
    if mode == "parapack":
        assert result.stdout.count("[input parameters]") == 2, result.stdout
        if launcher:
            for rank in range(2):
                assert result.stdout.count(f"[results {rank}]") == 2, result.stdout
        else:
            assert result.stdout.count("[results]") == 2, result.stdout
    elif mode in ("owned-query", "single-query"):
        assert result.stdout == expected, result.stdout
