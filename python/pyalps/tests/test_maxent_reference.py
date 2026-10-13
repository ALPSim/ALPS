"""Apply the same physical oracle to the Python binding and native executable."""
import importlib.util
import os
from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[3] / "src/apps/maxent/tests/reference.py"
if not SOURCE.is_file():
    pytest.skip("shared MaxEnt reference requires the ALPS checkout", allow_module_level=True)
SPEC = importlib.util.spec_from_file_location("maxent_reference", SOURCE)
reference = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(reference)


@pytest.fixture(scope="module")
def reference_runs(tmp_path_factory):
    maxent = pytest.importorskip("pyalps.maxent_c")
    return reference.run_cases(tmp_path_factory.mktemp("maxent-python"),
                               maxent.AnalyticContinuation)


def test_maxent_physical_reference_cases(reference_runs):
    # The shared runner validates every case, covariance and normalization pairs.
    assert reference_runs


def test_maxent_cli_matches_python(reference_runs, tmp_path):
    from pyalps_cli import resolve_executable

    executable = os.environ.get("ALPS_MAXENT_EXECUTABLE")
    if executable:
        executable = Path(executable).resolve(strict=True)
    else:
        try:
            executable = resolve_executable("maxent.exe" if os.name == "nt" else "maxent")
        except RuntimeError as error:
            pytest.skip(str(error))
    actual = reference.run_cases(tmp_path, lambda parameters: reference.run_cli(executable, parameters))
    for name, expected in reference_runs.items():
        reference.compare_results(actual[name], expected, relative=2e-6, absolute=1e-9)
