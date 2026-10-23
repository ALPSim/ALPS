"""Check MaxEnt text columns against its scientific HDF5 output, and bootstrap seeding."""

import argparse
from pathlib import Path
import os
import subprocess
import tempfile

import h5py
import numpy as np


BETA, COUNT, FREQUENCIES, ALPHAS, NORM = 4.0, 16, 48, 6, 2.5


def solve(executable, directory, returncode=0, covariance=None, **overrides):
    beta, count, frequencies, alphas, norm = BETA, COUNT, FREQUENCIES, ALPHAS, NORM
    tau = np.linspace(0.0, beta, count)
    poles = np.array([-1.0, 1.0])
    data = -np.exp(-np.outer(tau, poles)) / (1.0 + np.exp(-beta * poles))
    data = data @ np.array([0.6, 0.4])
    source = directory / "input.h5"
    parameters = {
        "BETA": beta, "NDAT": count, "NFREQ": frequencies,
        "N_ALPHA": alphas, "ALPHA_MIN": 0.1, "ALPHA_MAX": 20.0,
        "MAX_IT": 120, "MAX_TIME": 60, "NORM": norm,
        "OMEGA_MIN": -4.0, "OMEGA_MAX": 4.0,
        "DEFAULT_MODEL": "flat", "FREQUENCY_GRID": "linear",
        "DATASPACE": "time", "KERNEL": "fermionic",
        "PARTICLE_HOLE_SYMMETRY": 0, "TEXT_OUTPUT": 1, "SELF": 1,
        "DATA_IN_HDF5": 1, "DATA": str(source),
        "BASENAME": str(directory / "result"),
        **overrides,
    }
    with h5py.File(source, "w") as archive:
        archive["Data"] = norm * data
        archive["Error"] = np.full(count, norm * 0.001)
        if covariance is not None:
            archive["Covariance"] = covariance
        for name, value in parameters.items():
            archive["parameters/" + name] = np.int32(value) if isinstance(value, int) else value
    environment = {k: v for k, v in os.environ.items() if k != "ALPS_NO_CITATIONS"}
    completed = subprocess.run([str(executable), str(source)], cwd=directory, env=environment,
                               capture_output=True, text=True, timeout=60)
    assert completed.returncode == returncode, completed.stdout + completed.stderr
    return completed


def check_citations(executable, directory):
    """The solver prints its citations once per run, and not for rejected input."""
    completed = solve(executable, directory)
    assert completed.stdout.count("Recommended citations for ALPS maximum entropy") == 1, completed.stdout
    for path in directory.glob("result*"):
        path.unlink()
    rejected = solve(executable, directory, returncode=1, NDAT=0)
    assert "Recommended citations for" not in rejected.stdout, rejected.stdout
    assert not list(directory.glob("result*")), "Rejected input must not write output"


def check_rejected_input(executable, directory):
    """Input the solver cannot use is an error naming the cause, before any output."""
    model = directory / "narrow_model.dat"
    # Linear extrapolation of this table is negative beyond |omega| = 3.
    np.savetxt(model, [[-2.0, 0.1], [0.0, 0.4], [2.0, 0.1]])
    cases = {
        "all covariance eigenvalues": dict(
            covariance=np.eye(COUNT).ravel() * 1e-12, COVARIANCE_MATRIX="HDF5", CM_EIGENVALUE_CUTOFF=1e-10),
        "Default model point < 0": dict(DEFAULT_MODEL=str(model)),
        "N_ALPHA must be at least 2": dict(N_ALPHA=1),
    }
    for diagnostic, overrides in cases.items():
        run = directory / ("rejected-" + diagnostic.split()[0])
        run.mkdir()
        rejected = solve(executable, run, returncode=1, **overrides)
        assert diagnostic in rejected.stderr, rejected.stdout + rejected.stderr
        assert not list(run.glob("result*")), diagnostic
        assert not (run / "deltaOmega.dat").exists(), diagnostic


def check(executable, directory):
    count, frequencies, alphas, norm = COUNT, FREQUENCIES, ALPHAS, NORM
    solve(executable, directory)

    def table(name, rows, columns):
        values = np.loadtxt(directory / name, ndmin=2)
        assert values.shape == (rows, columns), (name, values.shape)
        assert np.all(np.isfinite(values)), name
        return values

    # The original diagnostic has three columns, including normalized bin mass.
    diagnostic = table("deltaOmega.dat", frequencies, 3)
    np.testing.assert_array_equal(diagnostic[:, 0], np.arange(frequencies))
    np.testing.assert_allclose(diagnostic[:, 1], 8.0 / frequencies, rtol=5.1e-6)
    np.testing.assert_allclose(diagnostic[:, 2], 1.0 / frequencies, rtol=5.1e-6)

    with h5py.File(directory / "result.out.h5", "r") as archive:
        omega = archive["spectrum/omega"][:]
        for label, key in (("avspec", "average"), ("maxspec", "maximum"), ("chispec", "chi")):
            spectrum = archive["spectrum/" + key][:]
            ordinary = table("result." + label + ".dat", frequencies, 3)
            np.testing.assert_allclose(ordinary[:, 0], omega, rtol=5.1e-6)
            np.testing.assert_allclose(ordinary[:, 1], spectrum, rtol=5.1e-6, atol=1e-12)
            np.testing.assert_allclose(ordinary[:, 2], norm / 8.0, rtol=5.1e-6)
            if label != "chispec":
                # Existing SELF outputs have two columns and Im Sigma = -pi*A.
                self_energy = table("result." + label + "_self.dat", frequencies, 2)
                np.testing.assert_allclose(self_energy[:, 0], omega, rtol=5.1e-6)
                np.testing.assert_allclose(self_energy[:, 1], -np.pi * spectrum,
                                           rtol=5.1e-6, atol=1e-12)
        probability = table("result.prob.dat", alphas, 2)
        np.testing.assert_allclose(probability[:, 0], archive["alpha/values"][:], rtol=5.1e-6)
        np.testing.assert_allclose(probability[:, 1], archive["alpha/probability"][:], rtol=5.1e-6)
    table("result.chi2.dat", alphas, 2)
    table("result.spex.dat", alphas * frequencies, 2)
    table("result.fits.dat", alphas * count, 4)


def check_bootstrap_seed(executable, directory):
    # GENERATE_ERR draws from SEED (default 42), so equal seeds reproduce the
    # bootstrap columns and a different seed changes them.
    errors = {}
    for label, seed in (("default", None), ("42", 42), ("7", 7)):
        run = directory / ("seed-" + label)
        run.mkdir()
        overrides = {"GENERATE_ERR": 1, "TEXT_OUTPUT": 0, "SELF": 0}
        if seed is not None:
            overrides["SEED"] = seed
        solve(executable, run, **overrides)
        errors[label] = np.loadtxt(run / "result.booterr.dat", ndmin=2)
        assert errors[label].shape == (FREQUENCIES, 4) and np.all(np.isfinite(errors[label]))
    np.testing.assert_array_equal(errors["default"], errors["42"])
    np.testing.assert_array_equal(errors["7"][:, :2], errors["42"][:, :2])
    assert np.any(errors["7"][:, 2:] != errors["42"][:, 2:])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    arguments = parser.parse_args()
    executable = arguments.executable.resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="alps-maxent-output-") as temporary:
        check(executable, Path(temporary))
        check_bootstrap_seed(executable, Path(temporary))
        check_rejected_input(executable, Path(temporary))
    with tempfile.TemporaryDirectory(prefix="alps-maxent-citations-") as temporary:
        check_citations(executable, Path(temporary))
    print("MaxEnt text and self-energy output checks passed")
