"""pyalps.alea after the C++ ALEA bindings: exact numbers for the NumPy time
series functions, and how errors combine."""

import numpy as np
import pytest

from pyalps import alea


def test_autocorrelation_of_an_exact_series():
    x = np.array([1., -1., 1., -1., 1., -1.])
    series = alea.MCScalarTimeseries(x)
    # mean 0, variance 6/5; lag i sums (6 - i) products of (-1)^i
    expected = [(-1) ** i * (6 - i) / (1.2 * (6 - i)) for i in (1, 2, 3)]
    np.testing.assert_allclose(alea.autocorrelation(series, _distance=3).timeseries(), expected)
    # stops at the first value below 0.5 times the first (-5/6): lag 2 is 5/6
    np.testing.assert_allclose(alea.autocorrelation(series, _limit=0.5).timeseries(), expected[:1])

    vector = alea.MCVectorTimeseries(np.column_stack([x, 2 * x]))
    np.testing.assert_allclose(alea.autocorrelation(vector, _distance=2).timeseries(),
                               np.column_stack([expected[:2], expected[:2]]))


def test_cuts_share_the_measurements():
    series = alea.MCScalarTimeseries(np.array([1., .8, .5, .3, .1]))
    head = alea.cut_head(series, _distance=2)
    assert isinstance(head, alea.MCScalarTimeseriesView)
    np.testing.assert_array_equal(head.timeseries(), [.5, .3, .1])
    np.testing.assert_array_equal(alea.cut_tail(head, _distance=-1).timeseries(), [.5])
    # up to and including the first value <= 0.5 times the first
    np.testing.assert_array_equal(alea.cut_head(series, _limit=0.5).timeseries(), [.3, .1])
    np.testing.assert_array_equal(alea.cut_tail(series, _limit=0.5).timeseries(), [1., .8, .5])


def test_exponential_fit_recovers_an_exact_exponential():
    lags = np.arange(1, 11)
    correlation = alea.MCScalarTimeseries(0.9 * np.exp(-lags / 3.))
    fit = alea.exponential_autocorrelation_time(correlation, _from=2, _to=8)
    assert isinstance(fit, alea.StdPairDouble)
    assert fit.first == pytest.approx(0.9)
    assert fit.second == pytest.approx(-1. / 3.)
    ranged = alea.exponential_autocorrelation_time(correlation, _max=0.8, _min=0.1)
    assert (ranged.first, ranged.second) == pytest.approx((0.9, -1. / 3.))
    tail = (0.9 * 3.) * np.exp(-10.5 / 3.)
    assert alea.integrated_autocorrelation_time(correlation, fit) == \
        pytest.approx(correlation.timeseries().sum() + tail)


def test_binning_error_is_that_of_a_real_observable():
    samples = np.random.default_rng(3).normal(size=1024)
    observable = alea.RealObservable("x")
    for sample in samples:
        observable << sample
    assert alea.error(samples, alea.binning) == observable.error
    assert observable.mean == pytest.approx(samples.mean(), rel=1e-12)


def test_numbers_with_errors_remember_what_they_depend_on():
    a = alea.MCScalarData(2.0, 0.1)
    b = alea.MCScalarData(3.0, 0.2)
    assert (a + a).error == pytest.approx(0.2)
    assert (a - a).error == 0.0
    assert (a * b / a).error == pytest.approx(0.2)
    assert (a + b).error == pytest.approx(np.hypot(0.1, 0.2))
    vector = alea.MCVectorData([1.0, 2.0], [0.1, 0.1])
    assert (vector[0] - vector[0]).error == 0.0
    assert (vector[0] + vector[1]).error == pytest.approx(np.hypot(0.1, 0.1))
    assert (vector * a - a * vector).error == pytest.approx([0.0, 0.0])


def test_measured_results_combine_by_jackknife():
    rng = np.random.default_rng(5)
    observable = alea.RealVectorObservable("x")
    for sample in rng.normal(size=(4096, 2)):
        observable << sample
    x = observable.result()
    assert x.has_bins and x.count == 4096
    assert (x - x).error == pytest.approx([0.0, 0.0], abs=1e-15)
    np.testing.assert_allclose((2 * x + 1).error, 2 * x.error, rtol=1e-12)
    # A number with an error is independent of the measurements.
    c = alea.MCScalarData(1.0, 0.5)
    np.testing.assert_allclose((x[0] + c).error, np.hypot(x.error[0], 0.5), rtol=1e-12)


def test_results_without_bins_only_scale_and_shift(tmp_path):
    path = str(tmp_path / "stored.h5")
    alea.MCScalarData(2.0, 0.1).save(path, "/simulation/results/a")
    a = alea.MCScalarData()
    a.load(path, "/simulation/results/a")
    assert not a.has_bins
    assert (3 * a - 1).error == pytest.approx(0.3)
    assert (a / alea.MCScalarData(4.0, 0.2)).error == pytest.approx(np.hypot(0.1 / 4, 2 * 0.2 / 16))
    with pytest.raises(ValueError, match="no bins"):
        a.sqrt()
    with pytest.raises(ValueError, match="no bins"):
        a + a
