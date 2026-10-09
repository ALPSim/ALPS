"""Shared physical reference checks for the MaxEnt executable and Python binding."""

import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import h5py


CASES_FILE = Path(__file__).with_name("reference_cases.json")


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
        "TEXT_OUTPUT": 0, "VERBOSE": 0, "MAX_TIME": 60, "CUT": 0.01,
        "DATA_IN_HDF5": 1, "DATA": str(directory / "input.h5"),
        "BASENAME": str(directory / "result"),
    }
    if case["dataspace"] == "time":
        parameters.update({f"TAU_{index}": float(point) for index, point in enumerate(points)})
    with h5py.File(parameters["DATA"], "w") as archive:
        archive["/Data"] = values * norm
        archive["/Error"] = errors * norm
        if case.get("covariance", False):
            diagonal = (errors * norm) ** 2
            if case.get("singular_last", False):
                diagonal[-1] = 0
            archive["/Covariance"] = np.diag(diagonal).ravel()
            parameters["COVARIANCE_MATRIX"] = "HDF5"
        for name, value in parameters.items():
            archive["/parameters/" + name] = value
    return parameters, points, values, errors


def read_result(filename, case):
    paths = {key: "/spectrum/" + key for key in
             ("omega", "average", "maximum", "chi", "variance")}
    paths.update(alpha="/alpha/values", probability="/alpha/probability")
    if case["kernel"] == "bosonic":
        paths.update(bosonic_average="/spectrum/bosonic/average",
                     bosonic_maximum="/spectrum/bosonic/maximum")
    with h5py.File(filename, "r") as archive:
        result = {key: np.asarray(archive[path]) for key, path in paths.items()}
    # Preserve the native test's check of the actual quadrature diagnostic.
    diagnostic = np.loadtxt(filename.parent / "deltaOmega.dat", ndmin=2)
    np.testing.assert_array_equal(diagnostic[:, 0], np.arange(len(result["omega"])))
    result["bin_widths"] = diagnostic[:, 1]
    return result


def validate(defaults, case, result, points, values, errors):
    norm = case.get("norm", 1.0)
    knots = np.linspace(0., 1., defaults["nfreq"] + 1)
    if case["grid"] == "Lorentzian":
        transformed = np.tan(np.pi * (0.01 + 0.98 * knots - 0.5))
        knots = (transformed - transformed[0]) / (transformed[-1] - transformed[0])
    edges = defaults["omega_min"] + (defaults["omega_max"] - defaults["omega_min"]) * knots
    widths = np.diff(edges)
    expected_omega = (edges[1:] + edges[:-1]) / 2
    assert result["omega"].shape == expected_omega.shape
    np.testing.assert_allclose(result["omega"], expected_omega, rtol=1e-13, atol=1e-13)
    for name, array in result.items():
        assert array.ndim == 1 and np.all(np.isfinite(array)), (case["name"], name)
    for name in ("average", "maximum", "chi", "variance"):
        assert result[name].shape == result["omega"].shape
        assert np.all(result[name] >= 0), (case["name"], name)
    np.testing.assert_allclose(result["bin_widths"], widths, rtol=5.1e-6, atol=1e-12)
    assert np.all(result["bin_widths"] > 0)
    for name in ("average", "maximum", "chi"):
        weights = result[name] * widths / norm
        assert abs(weights.sum() - 1) <= 0.05, (case["name"], name, "spectral mass", weights.sum())
        reconstructed = kernel(case, defaults["beta"], points, result["omega"]) @ weights
        residual = (reconstructed - values) / errors
        if case.get("singular_last", False):
            residual = residual[:-1]
        rms = np.sqrt(np.mean(residual ** 2))
        assert rms <= 5, (case["name"], name, "forward RMS in sigma units", rms)
    alpha, probability = result["alpha"], result["probability"]
    assert len(alpha) == defaults["n_alpha"] and probability.shape == alpha.shape
    assert np.all(np.diff(alpha) < 0) and np.all(probability >= 0)
    np.testing.assert_allclose(alpha, np.geomspace(defaults["alpha_max"], defaults["alpha_min"],
                                                defaults["n_alpha"]), rtol=1e-12)
    integral = np.sum((probability[:-1] + probability[1:]) * -np.diff(alpha) / 2)
    assert abs(integral - 1) < 1e-12
    if case["kernel"] == "bosonic":
        for name in ("average", "maximum"):
            assert result["bosonic_" + name].shape == result["omega"].shape
            np.testing.assert_allclose(result["bosonic_" + name],
                                       np.pi * result["omega"] * result[name], rtol=1e-12)


def run_cases(root, solve):
    definitions = json.loads(CASES_FILE.read_text())
    runs = {}
    for case in definitions["cases"]:
        directory = root / case["name"]
        directory.mkdir()
        parameters, points, values, errors = prepare(definitions["defaults"], case, directory)
        with working_directory(directory):
            solve(parameters)
        result = read_result(directory / "result.out.h5", case)
        validate(definitions["defaults"], case, result, points, values, errors)
        runs[case["name"]] = result
    compare_results(runs["time_covariance"], runs["time_linear"])
    compare_results(runs["time_covariance_singular"], runs["time_omitted"])
    scaled = next(case for case in definitions["cases"] if case["name"] == "time_scaled")
    compare_results(runs["time_scaled"], runs["time_linear"], scaled["norm"])
    return runs


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
                                   atol=(1e-13 + 1e-13 * magnitude if name in ("omega", "alpha")
                                         else absolute + relative * magnitude), err_msg=name)


def run_cli(executable, parameters):
    process = subprocess.run([str(executable), parameters["DATA"]],
                             text=True, capture_output=True, timeout=60)
    assert process.returncode == 0, process.stdout + process.stderr


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="alps-maxent-reference-") as temporary:
        runs = run_cases(Path(temporary), lambda parameters: run_cli(executable, parameters))
        print(f"MaxEnt reference checks passed: {len(runs)} cases")
