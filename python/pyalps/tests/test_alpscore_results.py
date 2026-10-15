"""Read results written by ALPSCore's accumulators, and legacy binning analyses."""
import gzip
import json
from pathlib import Path

import numpy as np
import pytest
import pyalps

FIXTURES = Path(__file__).with_name("fixtures")
ALPSCORE = FIXTURES / "alpscore_results" / "results_7146b9e1.h5"
EXPECTED = json.loads((FIXTURES / "alpscore_results" / "expected.json").read_text())


def by_observable(datasets):
    return {dataset.props["observable"]: dataset for dataset in pyalps.flatten(datasets)}


def test_alpscore_parameters_and_measurements():
    data = by_observable(pyalps.loadMeasurements([str(ALPSCORE)]))
    assert set(data) == set(EXPECTED)
    for name, dataset in data.items():
        # ALPSCore's /parameters group: per-key datasets plus attributes pyalps ignores.
        assert dataset.props["L"] == 4
        assert dataset.props["T"] == 1.5
        assert dataset.props["MODEL"] == "ising"
        assert dataset.props["MEASURE"]
        expected = EXPECTED[name]
        if "error" not in expected:
            np.testing.assert_array_equal(dataset.y, [expected["mean"]])
            continue
        if np.ndim(expected["mean"]) == 0:
            values = list(dataset.y)
        elif "error_bins" in expected and name.endswith("FullBinning"):
            # Full binning stores a time series, so pyalps reads a vector result.
            values = [dataset.y[i] for i in range(len(expected["mean"]))]
        else:
            # Without a time series pyalps reads vector means only, as for legacy files.
            np.testing.assert_array_equal(dataset.y, expected["mean"])
            continue
        np.testing.assert_array_equal([value.mean for value in values], np.atleast_1d(expected["mean"]))
        np.testing.assert_array_equal([value.error for value in values], np.atleast_1d(expected["error"]))
        assert all(value.count == expected["count"] for value in values)


def test_alpscore_binning_analysis():
    data = by_observable(pyalps.loadBinningAnalysis([str(ALPSCORE)]))
    with_bins = {name: value for name, value in EXPECTED.items() if "error_bins" in value}
    assert set(data) == {f"binning analysis of {name}" for name in with_bins}
    for name, expected in with_bins.items():
        dataset = data[f"binning analysis of {name}"]
        np.testing.assert_array_equal(dataset.y, expected["error_bins"])
        np.testing.assert_array_equal(dataset.x, np.arange(len(expected["error_bins"])))
        assert dataset.props["hdf5_path"] == f"/simulation/results/{name}"
        # The last level is the error ALPSCore reports.
        np.testing.assert_array_equal(dataset.y[-1], expected["error"])


def test_legacy_binning_analysis(tmp_path):
    archive = tmp_path / "spinmc.h5"
    archive.write_bytes(gzip.decompress(
        (FIXTURES / "application_results" / "spinmc.h5.gz").read_bytes()))
    data = by_observable(pyalps.loadBinningAnalysis([str(archive)], what="Energy"))
    dataset = data["binning analysis of Energy"]
    # pyalps's estimate from the frozen file's timeseries/logbinning datasets,
    # pinned when the reader learned ALPSCore's layout.
    np.testing.assert_allclose(dataset.y, [
        0.0037538385945939637, 0.004207532555667101, 0.004544246332006187,
        0.004714680648957985, 0.004823537789552665, 0.004867844388314303,
        0.004923664488385224, 0.004938113486348119, 0.00496460135933891,
        0.005086211181013936, 0.005037883814616264, 0.005259472282441837], rtol=1e-14)
    assert dataset.props["hdf5_path"] == "/simulation/results/Energy"
