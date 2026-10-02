"""Physical reference cases shared with the native MaxEnt acceptance test."""

from contextlib import contextmanager
import json
import os
from pathlib import Path
import subprocess

import numpy as np
import pytest

from pyalps import hdf5, ngs


CASES_FILE = (Path(__file__).resolve().parents[2] /
              "src/apps/maxent/tests/reference_cases.json")


@contextmanager
def working_directory(directory):
    # MaxEnt currently writes deltaOmega.dat even with TEXT_OUTPUT disabled.
    previous = Path.cwd()
    os.chdir(directory)
    try:
        yield
    finally:
        os.chdir(previous)


def kernel(case, beta, points, omega):
    if case["dataspace"] == "time":
        return -np.exp(-np.outer(points, omega)) / (1 + np.exp(-beta * omega))
    frequency = (2 * points + (case["kernel"] == "fermionic")) * np.pi / beta
    denominator = frequency[:, None] ** 2 + omega[None, :] ** 2
    if case["kernel"] == "fermionic":
        return -frequency[:, None] / denominator
    return -omega[None, :] ** 2 / denominator


def prepare(defaults, case, directory):
    count = defaults["ndat"]
    points = np.arange(count, dtype=float)
    if case["dataspace"] == "time":
        points *= defaults["beta"] / (count - 1)
    errors = defaults["sigma"] * (1 + np.arange(count) / count)
    poles = np.asarray(case["poles"], dtype=float)
    values = kernel(case, defaults["beta"], points, poles[:, 0]) @ poles[:, 1]
    if case.get("omit_last", False):
        points, errors, values = points[:-1], errors[:-1], values[:-1]
    norm = case.get("norm", 1.0)
    parameters = {
        "BETA": defaults["beta"], "NDAT": len(values), "NFREQ": defaults["nfreq"],
        "N_ALPHA": defaults["n_alpha"], "ALPHA_MIN": defaults["alpha_min"],
        "ALPHA_MAX": defaults["alpha_max"], "MAX_IT": defaults["max_it"],
        "OMEGA_MIN": defaults["omega_min"], "OMEGA_MAX": defaults["omega_max"],
        "NORM": norm, "KERNEL": case["kernel"], "DATASPACE": case["dataspace"],
        "FREQUENCY_GRID": case["grid"], "DEFAULT_MODEL": "flat",
        "PARTICLE_HOLE_SYMMETRY": int(case["dataspace"] == "frequency"),
        "TEXT_OUTPUT": 0, "VERBOSE": 0, "MAX_TIME": 60,
        "DATA_IN_HDF5": 1, "DATA": str(directory / "input.h5"),
        "BASENAME": str(directory / "python"),
    }
    if case["dataspace"] == "time":
        parameters.update({f"TAU_{index}": float(point) for index, point in enumerate(points)})
    with hdf5.archive(parameters["DATA"], "w") as archive:
        archive["/Data"] = values * norm
        archive["/Error"] = errors * norm
        if case.get("covariance", False):
            diagonal = (errors * norm) ** 2
            if case.get("singular_last", False):
                diagonal[-1] = 0
            archive["/Covariance"] = np.diag(diagonal).ravel()
            parameters["COVARIANCE_MATRIX"] = "HDF5"
        archive["/parameters"] = ngs.params(parameters)
    return parameters, points, values, errors


def read_result(filename, case):
    paths = {key: "/spectrum/" + key for key in
             ("omega", "average", "maximum", "chi", "variance")}
    paths.update(alpha="/alpha/values", probability="/alpha/probability")
    if case["kernel"] == "bosonic":
        paths.update(bosonic_average="/spectrum/bosonic/average",
                     bosonic_maximum="/spectrum/bosonic/maximum")
    with hdf5.archive(str(filename), "r") as archive:
        return {key: np.asarray(archive[path]) for key, path in paths.items()}


def validate(defaults, case, result, points, values, errors):
    norm = case.get("norm", 1.0)
    knots = np.linspace(0., 1., defaults["nfreq"] + 1)
    if case["grid"] == "Lorentzian":
        transformed = np.tan(np.pi * (0.01 + 0.98 * knots - 0.5))
        knots = (transformed - transformed[0]) / (transformed[-1] - transformed[0])
    edges = defaults["omega_min"] + (defaults["omega_max"] - defaults["omega_min"]) * knots
    widths = np.diff(edges)
    expected_omega = (edges[1:] + edges[:-1]) / 2
    np.testing.assert_allclose(result["omega"], expected_omega, rtol=1e-12, atol=1e-13)
    for name, array in result.items():
        assert array.ndim == 1 and np.all(np.isfinite(array)), (case["name"], name)
    for name in ("average", "maximum", "chi", "variance"):
        assert result[name].shape == result["omega"].shape
        assert np.all(result[name] >= 0), (case["name"], name)
    weights = result["average"] * widths / norm
    assert abs(weights.sum() - 1) < 0.05, (case["name"], "spectral mass", weights.sum())
    reconstructed = kernel(case, defaults["beta"], points, result["omega"]) @ weights
    residual = (reconstructed - values) / errors
    if case.get("singular_last", False):
        residual = residual[:-1]
    rms = np.sqrt(np.mean(residual ** 2))
    assert rms < 5, (case["name"], "forward RMS in sigma units", rms)
    alpha, probability = result["alpha"], result["probability"]
    assert len(alpha) == defaults["n_alpha"] and probability.shape == alpha.shape
    assert np.all(np.diff(alpha) < 0) and np.all(probability >= 0)
    np.testing.assert_allclose(alpha, np.geomspace(defaults["alpha_max"], defaults["alpha_min"],
                                                defaults["n_alpha"]), rtol=1e-12)
    integral = np.sum((probability[:-1] + probability[1:]) * -np.diff(alpha) / 2)
    assert abs(integral - 1) < 1e-12
    if case["kernel"] == "bosonic":
        for name in ("average", "maximum"):
            np.testing.assert_allclose(result["bosonic_" + name],
                                       np.pi * result["omega"] * result[name], rtol=1e-12)


@pytest.fixture(scope="module")
def reference_runs(tmp_path_factory):
    maxent = pytest.importorskip("pyalps.maxent_c")
    definitions = json.loads(CASES_FILE.read_text())
    root = tmp_path_factory.mktemp("maxent-reference")
    runs = {}
    for case in definitions["cases"]:
        directory = root / case["name"]
        directory.mkdir()
        parameters, points, values, errors = prepare(definitions["defaults"], case, directory)
        with working_directory(directory):
            maxent.AnalyticContinuation(parameters)
        result = read_result(directory / "python.out.h5", case)
        runs[case["name"]] = (case, directory, result, points, values, errors)
    return definitions["defaults"], runs


def test_maxent_physical_reference_cases(reference_runs):
    defaults, runs = reference_runs
    for case, _, result, points, values, errors in runs.values():
        validate(defaults, case, result, points, values, errors)


def compare_results(actual, expected, scale=1., relative=2e-5, absolute=1e-8):
    for name in expected:
        factor = (scale ** 2 if name == "variance" else scale
                  if name in ("average", "maximum", "chi", "bosonic_average", "bosonic_maximum")
                  else 1.)
        normalized = actual[name] / factor
        assert normalized.shape == expected[name].shape
        assert np.all(np.isfinite(normalized)) and np.all(np.isfinite(expected[name]))
        # Match the native normwise comparison: tiny tails and variances should
        # not amplify harmless roundoff into large pointwise relative errors.
        magnitude = max(np.max(np.abs(normalized)), np.max(np.abs(expected[name])))
        np.testing.assert_allclose(normalized, expected[name], rtol=0,
                                   atol=absolute + relative * magnitude, err_msg=name)


def test_maxent_covariance_and_normalization(reference_runs):
    _, runs = reference_runs
    result = lambda name: runs[name][2]
    compare_results(result("time_covariance"), result("time_linear"))
    compare_results(result("time_covariance_singular"), result("time_omitted"))
    compare_results(result("time_scaled"), result("time_linear"), runs["time_scaled"][0]["norm"])


def test_maxent_cli_matches_python(reference_runs):
    # SDK integration jobs set ALPS_DIR. Wheel-only smoke jobs may lack the SDK.
    executable = os.environ.get("ALPS_MAXENT_EXECUTABLE")
    if executable:
        executable = Path(executable).resolve()
        assert executable.is_file(), f"ALPS_MAXENT_EXECUTABLE does not exist: {executable}"
    if not executable and os.environ.get("ALPS_DIR"):
        executable = (Path(os.environ["ALPS_DIR"]).resolve().parents[1] /
                      "bin" / ("maxent.exe" if os.name == "nt" else "maxent"))
    if not executable or not Path(executable).is_file():
        pytest.skip("requires the MaxEnt CLI from an SDK or ALPS_MAXENT_EXECUTABLE")
    defaults, runs = reference_runs
    for case, directory, expected, points, values, errors in runs.values():
        input_file = directory / "input.h5"
        with hdf5.archive(str(input_file), "a") as archive:
            archive["/parameters/BASENAME"] = str(directory / "cli")
        process = subprocess.run([str(executable), str(input_file)], cwd=directory,
                                 text=True, capture_output=True, timeout=60)
        assert process.returncode == 0, process.stdout + process.stderr
        actual = read_result(directory / "cli.out.h5", case)
        validate(defaults, case, actual, points, values, errors)
        compare_results(actual, expected, relative=2e-6, absolute=1e-9)
