# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
"""The public no-linear-drift diagnostic and its confidence-level contract."""

import numpy as np
import pytest

import pyalps


def check_series(monkeypatch, values, **kwargs):
    monkeypatch.setattr(pyalps, "loadTimeSeries", lambda *_: np.asarray(values))
    return pyalps.checkSteadyState(outfile="sample.h5", observable="Energy", **kwargs)


@pytest.mark.parametrize("confidence,cutoff", [
    (0.6827, 1.0000217133229992),  #two-sided standard normal quantiles, (confidence,num_std)
    (0.9, 1.6448536269514722),
    (0.95, 1.959963984540054),
])
def test_central_confidence_threshold(monkeypatch, confidence, cutoff):
    result = check_series(monkeypatch, [1., 0., 0., 1.],
                          confidenceInterval=confidence, includeLog=True)
    assert result["statistics"]["z0"] == pytest.approx(cutoff)
    assert result["statistics"]["confidenceInterval"] == confidence


def test_higher_confidence_widens_no_drift_acceptance(monkeypatch):
    # For [0, 0, 1, 0], slope=0.1 and the legacy slope std=sqrt(0.05).
    # z=sqrt(0.2) lies between the 20% and 68.27% central cutoffs.
    low = check_series(monkeypatch, [0., 0., 1., 0.], confidenceInterval=0.2)
    default = check_series(monkeypatch, [0., 0., 1., 0.])
    assert not low["value"]
    assert default["value"]


def test_preserves_slope_statistic_and_log_schema(monkeypatch):
    result = check_series(monkeypatch, [0., 0., 1., 0.], includeLog=True)
    assert set(result) == {"value", "props", "statistics"}
    assert result["props"] == {"outfile": "sample.h5", "observable": "Energy"}
    stats = result["statistics"]
    assert set(stats) == {"beta1", "confidenceInterval", "z", "z0"}
    assert stats["beta1"]["value"] == pytest.approx(0.1)
    assert stats["beta1"]["std"] == pytest.approx(np.sqrt(0.05))
    assert stats["z"] == pytest.approx(np.sqrt(0.2))


@pytest.mark.parametrize("value", [0., 3., -2.5])
def test_constant_series_has_no_linear_drift(monkeypatch, value):
    result = check_series(monkeypatch, np.full(8, value), includeLog=True)
    assert result["value"]
    assert result["statistics"]["beta1"] == {"value": 0., "std": 0.}
    assert result["statistics"]["z"] == 0.


@pytest.mark.parametrize("sign", [-1., 1.])
def test_clear_linear_drift_is_rejected(monkeypatch, sign):
    result = check_series(monkeypatch, sign * np.arange(100.), confidenceInterval=0.95)
    assert result == {"value": False}


@pytest.mark.parametrize("values", [[], [1.], [1., np.nan, 2.], [1., np.inf, 2.],
                                    [[1., 2.], [3., 4.]], [1j, 2j], ["x", "y"]])
def test_invalid_series_is_rejected(monkeypatch, values):
    with pytest.raises(ValueError, match="time series"):
        check_series(monkeypatch, values)


@pytest.mark.parametrize("confidence", [0., 1., -0.1, 1.1, np.nan, np.inf, None, [0.9]])
def test_invalid_confidence_is_rejected_before_loading(monkeypatch, confidence):
    def unexpected_load(*_):
        pytest.fail("invalid confidence should be rejected before reading a file")
    monkeypatch.setattr(pyalps, "loadTimeSeries", unexpected_load)
    with pytest.raises(ValueError, match="confidenceInterval"):
        pyalps.checkSteadyState(outfile="unused", observable="Energy",
                               confidenceInterval=confidence)


def test_two_sample_series_remains_supported(monkeypatch):
    result = check_series(monkeypatch, [0., 1.], confidenceInterval=0.9, includeLog=True)
    assert result["statistics"]["z"] == pytest.approx(1.)
    assert result["value"]


def test_nested_datasets_are_annotated_in_place(monkeypatch):
    first, second = pyalps.DataSet(), pyalps.DataSet()
    first.props.update(filename="first.h5", observable="Energy")
    second.props.update(filename="second.h5", observable="Magnetization")
    samples = {"first.h5": np.zeros(4), "second.h5": np.arange(100.)}
    monkeypatch.setattr(pyalps, "loadTimeSeries", lambda filename, _: samples[filename])
    result = pyalps.checkSteadyState([[first], [second]], confidenceInterval=0.9)
    assert len(result) == 2 and result[0] is first and result[1] is second
    for dataset, expected in [(first, True), (second, False)]:
        log = dataset.props["checkSteadyState"]
        assert bool(log["value"]) == expected
        assert log["props"]["observable"] == dataset.props["observable"]
        assert log["statistics"]["confidenceInterval"] == 0.9


@pytest.mark.parametrize("suffix", [".h5", ".xml"])
def test_reads_real_hdf5_timeseries(tmp_path, suffix):
    from pyalps.hdf5 import archive
    filename = str(tmp_path / "measurement.h5")
    with archive(filename, "w") as output:
        output["/simulation/results/Energy/timeseries/data"] = np.array([0., 0., 1., 0.])
    result = pyalps.checkSteadyState(
        outfile=str(tmp_path / ("measurement" + suffix)), observable="Energy", includeLog=True)
    assert result["value"]
    assert result["statistics"]["z"] == pytest.approx(np.sqrt(0.2))


def test_threshold_near_one_stays_finite(monkeypatch):
    result = check_series(monkeypatch, [1., 0., 0., 1.],
                          confidenceInterval=np.nextafter(1., 0.), includeLog=True)
    assert np.isfinite(result["statistics"]["z0"])
