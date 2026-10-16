"""Check MaxEnt text columns against its scientific HDF5 output, and bootstrap seeding."""

import argparse
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


BETA, COUNT, FREQUENCIES, ALPHAS, NORM = 4.0, 16, 48, 6, 2.5


def solve(executable, directory, **overrides):
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
        for name, value in parameters.items():
            archive["parameters/" + name] = np.int32(value) if isinstance(value, int) else value
    completed = subprocess.run([str(executable), str(source)], cwd=directory,
                               capture_output=True, text=True, timeout=60)
    assert completed.returncode == 0, completed.stdout + completed.stderr


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
    print("MaxEnt text and self-energy output checks passed")
