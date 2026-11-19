# ALPS Project: Algorithms and Libraries for Physics Simulations
# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
"""Time series analysis of pyalps.alea, in NumPy.

A time series is a sequence of measurements: an ``MCScalarTimeseries`` or
``MCVectorTimeseries``, a view of one that leaves out its head or tail, a
NumPy array, or the bins of an ``MCScalarData`` or ``MCVectorData``.
"""

import math as _math
import warnings as _warnings

import numpy as _np

from ._mcdata import _MCData

__all__ = [
    "MCScalarTimeseries", "MCScalarTimeseriesView", "MCVectorTimeseries", "MCVectorTimeseriesView",
    "ValueWithError", "StdPairDouble",
    "size", "mean", "variance", "uncorrelated_error", "binning_error",
    "running_mean", "reverse_running_mean",
    "autocorrelation_distance", "autocorrelation_limit",
    "cut_head_distance", "cut_tail_distance", "cut_head_limit", "cut_tail_limit",
    "exponential_autocorrelation_time_distance", "exponential_autocorrelation_time_limit",
    "integrated_autocorrelation_time",
]


class NotEnoughMeasurementsError(ValueError):
    def __init__(self):
        ValueError.__init__(self, "Not enough measurements available.")


class _Timeseries:
    _ndim = 1

    def __init__(self, data=None):
        if data is None:
            values = _np.zeros((0,) * self._ndim)
        elif isinstance(data, _MCData):
            values = data.bins
        elif isinstance(data, _Timeseries):
            values = data._values
        else:
            values = data
        values = _np.array(values, dtype=float)
        if values.ndim != self._ndim:
            raise ValueError("a %s takes a %d-dimensional array" % (type(self).__name__, self._ndim))
        self._values = values

    def timeseries(self):
        """The measurements, as a new array."""
        return self._values.copy()

    @property
    def size(self):
        return len(self._values)

    def __len__(self):
        return len(self._values)

    def __repr__(self):
        return "%s(%s)" % (type(self).__name__, _np.array2string(self._values, separator=", "))


class MCScalarTimeseries(_Timeseries):
    """A time series of numbers."""


class MCVectorTimeseries(_Timeseries):
    """A time series of vectors, one row per measurement."""
    _ndim = 2


class MCScalarTimeseriesView(MCScalarTimeseries):
    """A part of a time series; it shares the measurements, not a copy."""

    def __init__(self, data):
        self._values = _values(data)


class MCVectorTimeseriesView(MCVectorTimeseries):
    """A part of a vector time series; it shares the measurements."""

    def __init__(self, data):
        self._values = _values(data)


class ValueWithError:
    """A mean and its error."""
    __slots__ = ("mean", "error")

    def __init__(self, mean, error):
        self.mean = float(mean)
        self.error = float(error)

    def __repr__(self):
        return "%g +/- %g" % (self.mean, self.error)


class StdPairDouble:
    """The parameters of an exponential fit, ``first * exp(second * t)``."""
    __slots__ = ("first", "second")

    def __init__(self, first=0.0, second=0.0):
        self.first = float(first)
        self.second = float(second)

    def __iter__(self):
        return iter((self.first, self.second))

    def __eq__(self, other):
        return tuple(self) == tuple(other)

    def __repr__(self):
        return "(%r, %r)" % (self.first, self.second)


# ---- helpers ---------------------------------------------------------

def _values(x):
    if isinstance(x, _Timeseries):
        return x._values
    if isinstance(x, _MCData):
        return x.bins
    return _np.asarray(x, dtype=float)


def _is_vector(x):
    return _values(x).ndim == 2


def _series(values, view=False):
    if values.ndim == 2:
        kind = MCVectorTimeseriesView if view else MCVectorTimeseries
    else:
        kind = MCScalarTimeseriesView if view else MCScalarTimeseries
    out = kind.__new__(kind)
    out._values = values
    return out


def _number(values):
    return float(values) if _np.ndim(values) == 0 else values


def _scalar_only(x, name):
    values = _values(x)
    if values.ndim != 1:
        raise TypeError("%s takes a time series of numbers" % name)
    return values


# ---- statistics ------------------------------------------------------

def size(x):
    """Number of measurements; the number of bins of an MCScalarData or
    MCVectorData."""
    return len(_values(x))


def mean(x):
    values = _values(x)
    if len(values) == 0:
        raise NotEnoughMeasurementsError()
    return _number(values.mean(axis=0))


def variance(x):
    """Variance of a single measurement, normalised by N - 1."""
    values = _values(x)
    if len(values) < 2:
        raise NotEnoughMeasurementsError()
    return _number(values.var(axis=0, ddof=1))


def uncorrelated_error(x):
    """Error of the mean if the measurements were independent."""
    return _number(_np.sqrt(_np.asarray(variance(x)) / len(_values(x))))


def binning_error(x):
    """Error of the mean from the binning analysis of the observable layer,
    as a RealObservable fed with these measurements would give."""
    from .cxx import pyngsobservable_c, pyngsresult_c
    values = _values(x)
    if values.ndim == 1:
        observable = pyngsobservable_c.createRealObservable("timeseries")
        for value in values:
            observable << float(value)
    else:
        observable = pyngsobservable_c.createRealVectorObservable("timeseries")
        for row in values:
            observable << _np.ascontiguousarray(row)
    result = pyngsresult_c.observable2result(observable)
    if any(flag != 0 for flag in result.converged_vector):
        _warnings.warn("The binning analysis of %d measurements has not converged; "
                       "the error may be too small" % len(values))
    return _number(_np.asarray(result.error_vector)[0] if values.ndim == 1 else result.error_vector)


def running_mean(x):
    """Means of the first 1, 2, ... measurements."""
    values = _scalar_only(x, "running_mean")
    return _series(_np.cumsum(values) / _np.arange(1, len(values) + 1))


def reverse_running_mean(x):
    """Means of the last N, N - 1, ... measurements."""
    values = _scalar_only(x, "reverse_running_mean")
    sums = _np.cumsum(values[::-1])[::-1]
    return _series(sums / _np.arange(len(values), 0, -1))


# ---- autocorrelation -------------------------------------------------

def _autocorrelation(values, lag, m, var):
    n = len(values)
    return ((values[:n - lag] - m) * (values[lag:] - m)).sum(axis=0) / (var * (n - lag))


def autocorrelation_distance(x, up_to):
    """Autocorrelation at the distances 1, 2, ..., up_to."""
    values = _values(x)
    n = len(values)
    if n < 2:
        raise NotEnoughMeasurementsError()
    if up_to < 0:
        up_to += n
    m, var = values.mean(axis=0), values.var(axis=0, ddof=1)
    out = []
    for lag in range(1, up_to + 1):
        if lag == n:
            _warnings.warn("Autocorrelation fully calculated with a size of %d" % (lag - 1))
            break
        out.append(_autocorrelation(values, lag, m, var))
    return _series(_np.array(out).reshape((len(out),) + values.shape[1:]))


def autocorrelation_limit(x, limit):
    """Autocorrelation from distance 1 until it, or one of its components,
    falls below `limit` times its value at distance 1."""
    values = _values(x)
    n = len(values)
    if n < 2:
        raise NotEnoughMeasurementsError()
    m, var = values.mean(axis=0), values.var(axis=0, ddof=1)
    out = []
    for lag in range(1, n):
        out.append(_autocorrelation(values, lag, m, var))
        if _np.any(out[-1] < limit * out[0]):
            break
    else:
        _warnings.warn("Autocorrelation fully calculated with a size of %d" % (n - 1))
    return _series(_np.array(out).reshape((len(out),) + values.shape[1:]))


# ---- views -----------------------------------------------------------

def cut_head_distance(x, cutoff):
    """The time series without its first `cutoff` measurements."""
    values = _values(x)
    if cutoff < 0:
        cutoff += len(values)
    return _series(values[cutoff:], view=True)


def cut_tail_distance(x, cutoff):
    """The time series without its last `cutoff` measurements."""
    values = _values(x)
    if cutoff < 0:
        cutoff += len(values)
    return _series(values[:len(values) - cutoff], view=True)


def _position_at_most(values, limit):
    """1-based position of the first value <= limit; the size if none is."""
    below = _np.flatnonzero(values <= limit)
    return int(below[0]) + 1 if len(below) else len(values)


def cut_head_limit(x, limit):
    """The time series from after the first measurement that is at most
    `limit` times the first one."""
    values = _scalar_only(x, "cut_head_limit")
    return _series(values[_position_at_most(values, limit * values[0]):], view=True)


def cut_tail_limit(x, limit):
    """The time series up to and including the first measurement that is at
    most `limit` times the first one."""
    values = _scalar_only(x, "cut_tail_limit")
    return _series(values[:_position_at_most(values, limit * values[0])], view=True)


# ---- exponential fits ------------------------------------------------

def exponential_autocorrelation_time_distance(x, start, stop):
    """Least-squares fit of ``first * exp(second * t)`` to the autocorrelation
    at the distances `start` to `stop`; the autocorrelation starts at
    distance 1."""
    values = _scalar_only(x, "exponential_autocorrelation_time")
    n = len(values)
    if start < 0:
        start += n
    if stop < 0:
        stop += n
    fitted = values[start - 1:stop]
    if len(fitted) < 2:
        raise ValueError("Regression Error: System is not determined!")
    second, intercept = _np.polyfit(_np.arange(start, start + len(fitted), dtype=float),
                                    _np.log(fitted), 1)
    return StdPairDouble(_math.exp(intercept), second)


def exponential_autocorrelation_time_limit(x, max, min):
    """The fit between the first distances where the autocorrelation is at
    most `max` and `min` times its value at distance 1."""
    values = _scalar_only(x, "exponential_autocorrelation_time")
    start = _position_at_most(values, max * values[0])
    stop = _position_at_most(values, min * values[0]) - 1
    if stop < start:
        _warnings.warn("Invalid Range! If you want to fit a positive exponential, exchange min and max.")
    return exponential_autocorrelation_time_distance(values, start, stop)


def integrated_autocorrelation_time(x, fit):
    """Sum of the autocorrelation, plus the tail of the exponential `fit`
    beyond it."""
    values = _scalar_only(x, "integrated_autocorrelation_time")
    first, second = fit
    return float(values.sum() - (first / second) * _math.exp(second * (len(values) + 0.5)))
