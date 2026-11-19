# ALPS Project: Algorithms and Libraries for Physics Simulations
# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
"""RealObservable and friends: observables of <alps/alea/observables.hpp>
under the names of ALEA, for Python simulations such as the code-01 tutorials.
"""

import numpy as _np

from .cxx import pyngsobservable_c as _observables
from .cxx import pyngsresult_c as _results
from .cxx.pyngshdf5_c import hdf5_archive_impl as _archive_type

__all__ = ["RealObservable", "RealVectorObservable",
           "RealTimeSeriesObservable", "RealVectorTimeSeriesObservable"]


class RealObservable:
    """Measurements of a number, with binning analysis and bins for functions
    of results. ``observable << value`` adds a measurement; ``save(filename)``
    writes the result to ``/simulation/results/<name>``."""
    _vector = False

    def __init__(self, name, bins=0):
        self._name = name
        self._bins = bins
        create = (_observables.createRealVectorObservable if self._vector
                  else _observables.createRealObservable)
        self._observable = create(name, bins)

    @property
    def name(self):
        return self._name

    def __lshift__(self, value):
        if self._vector:
            value = _np.ascontiguousarray(value, dtype=float)
        elif not isinstance(value, (int, float)):
            value = float(value)
        self._observable << value
        return self

    def __deepcopy__(self, memo):
        out = type(self)(self._name, self._bins)
        out._observable.merge(self._observable)
        return out

    def __repr__(self):
        return self._name

    def _result(self):
        return _results.observable2result(self._observable)

    def _numbers(self, name):
        values = _np.asarray(getattr(self._result(), name))
        return values if self._vector else values[0]

    @property
    def count(self):
        return self._observable.count

    @property
    def mean(self):
        return self._numbers("mean_vector")

    @property
    def error(self):
        return self._numbers("error_vector")

    @property
    def tau(self):
        return self._numbers("tau_vector")

    @property
    def variance(self):
        return self._numbers("variance_vector")

    @property
    def converged_errors(self):
        """0 converged, 1 maybe converged, 2 not converged."""
        return _np.asarray(self._result().converged_vector, dtype=_np.int32) if self._vector \
            else int(self._result().converged_vector[0])

    def result(self):
        """The result, as an MCScalarData or MCVectorData."""
        from ._mcdata import MCScalarData, MCVectorData
        return (MCVectorData if self._vector else MCScalarData)._wrap(self._result())

    def save(self, target):
        """Write the result to ``/simulation/results/<name>`` in the file
        ``target``, or at the context of the archive ``target``."""
        from .hdf5 import archive
        if isinstance(target, _archive_type):
            self._result().save(target)
            return
        with archive(target, "a") as ar:
            ar["/simulation/results/" + self._name] = self._result()


class RealVectorObservable(RealObservable):
    """Measurements of a vector, all of the size of the first."""
    _vector = True


class RealTimeSeriesObservable(RealObservable):
    """A RealObservable. The time series is kept as its bins."""


class RealVectorTimeSeriesObservable(RealVectorObservable):
    """A RealVectorObservable. The time series is kept as its bins."""
