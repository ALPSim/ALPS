"""Test the installed wheel and retain results across cibuildwheel interpreters."""

from pathlib import Path
import platform
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[2]


def main():
    import pyalps

    package = Path(pyalps.__file__).resolve()
    if package.is_relative_to(ROOT / "python"):
        raise RuntimeError(f"Wheel tests imported source-tree pyalps: {package}")
    reports = ROOT / "_build/wheel-reports"
    reports.mkdir(parents=True, exist_ok=True)
    identifier = f"{platform.python_implementation()}-{platform.python_version()}-{platform.machine()}"
    (reports / f"{identifier}-import.txt").write_text(f"Python: {sys.executable}\npyalps: {package}\n")
    # The C++ SDK is identical across wheel ABIs. Compile its consumers once
    # per platform; Python/downstream-extension tests still run on every ABI.
    suites = [str(ROOT / 'python/pyalps/tests')]
    if sys.version_info[:2] == (3, 12):
        suites.append(str(ROOT / 'tests/cmake'))
    elif sys.version_info[:2] > (3, 12):
        # Reused abi3 artifact: check interpreter-facing behavior, without
        # repeating its numerical regressions on every newer interpreter.
        suites = [str(ROOT / 'python/pyalps/tests' / name) for name in (
            'test_binding_surface.py', 'test_mapping_lifetimes.py', 'test_wheel_payload.py',
        )]
    with tempfile.TemporaryDirectory(prefix="alps-wheel-tests-") as work:
        return subprocess.call([
            sys.executable, "-m", "pytest", "-v", "-o", "faulthandler_timeout=300",
            *suites,
            f"--junitxml={reports / f'{identifier}.xml'}",
        ], cwd=work)


if __name__ == "__main__":
    raise SystemExit(main())
