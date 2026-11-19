# ALPS Project: Algorithms and Libraries for Physics Simulations
# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
"""MCScalarData and MCVectorData: Monte Carlo results with errors.

A value is the sum of two independent parts:

- a result of ``<alps/alea/observables.hpp>`` (a ``pyalps.ngs.result``), for
  values measured or loaded from a file. Functions of results with bins
  propagate errors by the jackknife over the bins, so correlations between
  observables of one simulation are kept. A result without bins, such as one
  read from XML, only scales and shifts by numbers; other functions of it are
  refused.
- a linear dependence on independent numbers with errors, made with
  ``MCScalarData(mean, error)`` or ``MCVectorData(mean, error)``. Functions of
  them propagate errors to first order and remember which numbers a value
  depends on, so that ``x + x`` has twice the error of ``x`` and ``x - x``
  none, as in the ``uncertainties`` package.
"""

import copy as _copy
import numbers as _numbers

import numpy as _np

from .cxx import pyngsresult_c as _results
from .cxx.pyngshdf5_c import hdf5_archive_impl as _archive_type

__all__ = ["MCScalarData", "MCVectorData"]


class _Source:
    """An independent number (or vector of numbers) with these errors."""
    __slots__ = ("sigma",)

    def __init__(self, sigma):
        self.sigma = sigma


def _vector(values):
    out = _np.array(values, dtype=float, ndmin=1)
    if out.ndim != 1:
        raise ValueError("means and errors must be numbers or one-dimensional")
    return out


class _MCData:
    """Common implementation; see the module documentation."""

    # numpy defers to our reflected operators instead of looping over us.
    __array_ufunc__ = None

    def _init(self, core=None, mean=None, grad=None, count=0):
        self._core = core
        self._mean = mean
        self._grad = {} if grad is None else grad
        self._count = count

    def _set(self, other):
        self._init(other._core, other._mean, other._grad, other._count)
        return self

    # ---- construction -------------------------------------------------

    @classmethod
    def _make(cls, core=None, mean=None, grad=None, count=0):
        out = cls.__new__(cls)
        out._init(core, mean, grad, count)
        return out

    @classmethod
    def _wrap(cls, core):
        """A value of the measured or loaded result `core`."""
        if core.count == 0:
            return cls._make(mean=_np.zeros(0 if cls is MCVectorData else 1))
        if cls is MCScalarData and core.size != 1:
            raise TypeError("a result with %d components is not a MCScalarData" % core.size)
        return cls._make(core=core)

    # ---- numbers ------------------------------------------------------

    def _means(self):
        if self._core is not None:
            return _np.asarray(self._core.mean_vector, dtype=float)
        return self._mean

    def _size(self):
        return len(self._means())

    def _errors(self):
        n = self._size()
        square = _np.zeros(n)
        if self._core is not None:
            square += _np.asarray(self._core.error_vector) ** 2
        for source, g in self._grad.items():
            square += _np.broadcast_to(((g * source.sigma) ** 2).sum(axis=1), (n,))
        return _np.sqrt(square)

    def _measured(self, name):
        if self._core is None:
            return _np.full(self._size(), _np.nan)
        return _np.asarray(getattr(self._core, name), dtype=float)

    @property
    def count(self):
        """Number of samples of a measured result; 1 for a number with an error."""
        return self._core.count if self._core is not None else self._count

    @property
    def bins(self):
        """Bin means of the measured part, one row per bin for a vector."""
        if self._core is None or not self._core.has_bins:
            return _np.zeros((0,) if self._scalar else (0, self._size()))
        matrix = _np.asarray(self._core.bin_matrix)
        return matrix[:, 0].copy() if self._scalar else matrix

    @property
    def has_bins(self):
        """Can the measured part enter functions of results other than scaling?"""
        return self._core is None or bool(self._core.has_bins)

    # ---- arithmetic ---------------------------------------------------

    @staticmethod
    def _operand(x):
        if isinstance(x, _MCData):
            return x
        return MCVectorData._make(mean=_vector(x))

    def _kind(self, other):
        scalar = self._scalar and (other._scalar if isinstance(other, _MCData)
                                   else _np.ndim(other) == 0)
        return MCScalarData if scalar else MCVectorData

    @staticmethod
    def _combine(terms, n):
        """Sum of d_k * G_k over the operands' gradients, d_k of shape (n,)."""
        grad = {}
        for derivative, g in terms:
            for source, matrix in g.items():
                part = _np.asarray(derivative)[:, None] * matrix
                part = _np.broadcast_to(part, (n, matrix.shape[1]))
                grad[source] = grad[source] + part if source in grad else part
        return {s: m for s, m in grad.items() if _np.any(m)}

    def _binary(self, other, op, reflected=False):
        kind = self._kind(other)
        a, b = self, self._operand(other)
        if reflected:
            a, b = b, a
        f, da, db = _BINARY[op]
        ma, mb = a._means(), b._means()
        if len(ma) != len(mb) and 1 not in (len(ma), len(mb)):
            raise ValueError("operands of %d and %d components do not combine" % (len(ma), len(mb)))
        mean = _np.asarray(f(ma, mb), dtype=float)
        n = len(mean)
        grad = self._combine([(_np.broadcast_to(da(ma, mb), (n,)), a._grad),
                              (_np.broadcast_to(db(ma, mb), (n,)), b._grad)], n)
        if a._core is None and b._core is None:
            return kind._make(mean=mean, grad=grad, count=max(a._count, b._count))
        if a._core is not None and b._core is not None:
            core = _LAYER[op](a._core, b._core)
        elif a._core is not None and len(mb) == 1:
            core = _LAYER[op](a._core, float(mb[0]))
        elif b._core is not None and len(ma) == 1:
            core = _LAYER[op](float(ma[0]), b._core)
        elif a._core is not None:
            core = _results.transform(lambda x: f(x, mb), [a._core])
        else:
            core = _results.transform(lambda x: f(ma, x), [b._core])
        return kind._make(core=core, grad=grad)

    def _unary(self, op, *args):
        f, df = _UNARY[op]
        m = self._means()
        grad = self._combine([(df(m, *args), self._grad)], len(m))
        if self._core is None:
            return type(self)._make(mean=_np.asarray(f(m, *args), dtype=float),
                                    grad=grad, count=self._count)
        if op == "neg":
            core = -self._core
        else:
            core = _results.transform(lambda x: f(x, *args), [self._core])
        return type(self)._make(core=core, grad=grad)

    def __add__(self, other): return self._binary(other, "add")
    def __radd__(self, other): return self._binary(other, "add", True)
    def __sub__(self, other): return self._binary(other, "sub")
    def __rsub__(self, other): return self._binary(other, "sub", True)
    def __mul__(self, other): return self._binary(other, "mul")
    def __rmul__(self, other): return self._binary(other, "mul", True)
    def __truediv__(self, other): return self._binary(other, "div")
    def __rtruediv__(self, other): return self._binary(other, "div", True)

    # In place, as the C++ types were: other references see the new value.
    def __iadd__(self, other): return self._set(self + other)
    def __isub__(self, other): return self._set(self - other)
    def __imul__(self, other): return self._set(self * other)
    def __itruediv__(self, other): return self._set(self / other)

    def __pos__(self): return type(self)._make(self._core, self._mean, dict(self._grad), self._count)
    def __neg__(self): return self._unary("neg")
    def __abs__(self): return self._unary("abs")

    def __pow__(self, exponent):
        if not isinstance(exponent, _numbers.Real):
            return NotImplemented
        return self._unary("pow", float(exponent))

    def sq(self): return self._unary("sq")
    def cb(self): return self._unary("cb")
    def sqrt(self): return self._unary("sqrt")
    def cbrt(self): return self._unary("cbrt")
    def exp(self): return self._unary("exp")
    def log(self): return self._unary("log")
    def sin(self): return self._unary("sin")
    def cos(self): return self._unary("cos")
    def tan(self): return self._unary("tan")
    def asin(self): return self._unary("asin")
    def acos(self): return self._unary("acos")
    def atan(self): return self._unary("atan")
    def sinh(self): return self._unary("sinh")
    def cosh(self): return self._unary("cosh")
    def tanh(self): return self._unary("tanh")
    def asinh(self): return self._unary("asinh")
    def acosh(self): return self._unary("acosh")
    def atanh(self): return self._unary("atanh")

    # ---- copies, merging, files ---------------------------------------

    def __copy__(self):
        return +self

    def __deepcopy__(self, memo):
        # The copy depends on the same numbers as the original: x - copy(x) is 0.
        return +self

    def merge(self, other):
        """Add the measurements of an independent run, in place."""
        if self._core is not None and other._core is not None and not self._grad and not other._grad:
            self._core = _results.join([self._core, other._core])
            return self
        if self._core is None and other._core is None:
            total = self._count + other._count
            if total == 0:
                return self
            merged = (self * self._count + other * other._count) / total
            merged._count = total
            return self._set(merged)
        raise ValueError("measurements and numbers with errors cannot be merged")

    def _result(self):
        """The value as a result, to write to a file."""
        if self._core is not None and not self._grad:
            return self._core
        if self._core is not None:
            raise ValueError("a function of measured results and numbers with errors cannot be "
                             "written; write its parts")
        return _results.stored(list(self._mean), list(self._errors()), max(self._count, 1),
                               self._scalar)

    def save(self, target, path=None):
        """Write to ``path`` in the file ``target``, or at the context of the
        archive ``target``."""
        if isinstance(target, _archive_type):
            if path is not None:
                target[path] = self._result()
            else:
                self._result().save(target)
            return
        from .hdf5 import archive
        with archive(target, "a") as ar:
            ar[path] = self._result()

    def load(self, source, path=None):
        """Read from ``path`` in the file ``source``, or from the context of
        the archive ``source``."""
        core = _results.result()
        if isinstance(source, _archive_type):
            if path is None:
                core.load(source)
            else:
                context = source.context
                full = path if path.startswith("/") else context.rstrip("/") + "/" + path
                source.set_context(full)
                try:
                    core.load(source)
                finally:
                    source.set_context(context)
        else:
            from .hdf5 import archive
            with archive(source, "r") as ar:
                ar.set_context(path)
                core.load(ar)
        return self._set(type(self)._wrap(core))

    # ---- printing -----------------------------------------------------

    def _lines(self, spec):
        return ["%s +/- %s" % (format(m, spec), format(e, spec))
                for m, e in zip(self._means(), self._errors())]


class MCScalarData(_MCData):
    """A number with an error: ``MCScalarData(mean, error)``, or the result
    of an observable from a file (see ``load``)."""
    _scalar = True

    def __init__(self, mean=None, error=0.0):
        if mean is None:
            self._init(mean=_np.zeros(1))
            return
        mean, error = float(mean), float(error)
        grad = {_Source(_np.array([error])): _np.ones((1, 1))} if error else {}
        self._init(mean=_np.array([mean]), grad=grad, count=1)

    @property
    def mean(self): return float(self._means()[0])

    @property
    def error(self): return float(self._errors()[0])

    @property
    def tau(self): return float(self._measured("tau_vector")[0])

    @property
    def variance(self): return float(self._measured("variance_vector")[0])

    def __repr__(self):
        return self._lines("g")[0]

    def __format__(self, spec):
        return self._lines(spec)[0]


class MCVectorData(_MCData):
    """Numbers with errors, ``MCVectorData(means, errors)``, independent of
    each other, or the result of a vector observable from a file."""
    _scalar = False

    def __init__(self, mean=None, error=None):
        if mean is None:
            self._init(mean=_np.zeros(0))
            return
        mean = _vector(mean)
        error = _np.zeros(len(mean)) if error is None else _vector(error)
        if len(error) != len(mean):
            raise ValueError("%d means but %d errors" % (len(mean), len(error)))
        grad = {_Source(error): _np.eye(len(mean))} if _np.any(error) else {}
        self._init(mean=mean, grad=grad, count=1)

    @property
    def mean(self): return self._means().copy()

    @property
    def error(self): return self._errors()

    @property
    def tau(self): return self._measured("tau_vector")

    @property
    def variance(self): return self._measured("variance_vector")

    def __len__(self):
        return self._size()

    def _select(self, indices, kind):
        indices = list(indices)
        grad = {s: g[indices] if g.shape[0] > 1 else _np.repeat(g, len(indices), axis=0)
                for s, g in self._grad.items()}
        if self._core is None:
            return kind._make(mean=self._mean[indices], grad=grad, count=self._count)
        if not indices:
            return kind._make(mean=_np.zeros(0))
        core = _results.select(self._core, indices, kind is MCScalarData)
        return kind._make(core=core, grad=grad)

    def __getitem__(self, key):
        n = self._size()
        if isinstance(key, slice):
            start, stop, step = key.indices(n)
            if step != 1:
                raise IndexError("slice step size not supported.")
            return self._select(range(start, max(start, stop)), MCVectorData)
        index = int(key)
        if index < 0:
            index += n
        if not 0 <= index < n:
            raise IndexError("Index out of range")
        return self._select([index], MCScalarData)

    def __iter__(self):
        return (self[i] for i in range(len(self)))

    def __eq__(self, other):
        if not isinstance(other, MCVectorData):
            return NotImplemented
        return (_np.array_equal(self._means(), other._means())
                and _np.array_equal(self._errors(), other._errors()))

    __hash__ = None

    def __contains__(self, value):
        if not isinstance(value, MCScalarData):
            return False
        return any(m == value.mean and e == value.error
                   for m, e in zip(self._means(), self._errors()))

    def __repr__(self):
        return "\n".join(self._lines("g"))

    def __format__(self, spec):
        return "\n".join(self._lines(spec))


def _cbrt_derivative(x):
    return 1.0 / (3.0 * _np.cbrt(x) ** 2)


_BINARY = {
    "add": (lambda a, b: a + b, lambda a, b: 1.0, lambda a, b: 1.0),
    "sub": (lambda a, b: a - b, lambda a, b: 1.0, lambda a, b: -1.0),
    "mul": (lambda a, b: a * b, lambda a, b: b, lambda a, b: a),
    "div": (lambda a, b: a / b, lambda a, b: 1.0 / b, lambda a, b: -a / (b * b)),
}

_LAYER = {
    "add": lambda a, b: a + b,
    "sub": lambda a, b: a - b,
    "mul": lambda a, b: a * b,
    "div": lambda a, b: a / b,
}

_UNARY = {
    "neg": (lambda x: -x, lambda x: -_np.ones_like(x)),
    "abs": (_np.abs, _np.sign),
    "pow": (lambda x, e: x ** e, lambda x, e: e * x ** (e - 1.0)),
    "sq": (lambda x: x * x, lambda x: 2.0 * x),
    "cb": (lambda x: x * x * x, lambda x: 3.0 * x * x),
    "sqrt": (_np.sqrt, lambda x: 0.5 / _np.sqrt(x)),
    "cbrt": (_np.cbrt, _cbrt_derivative),
    "exp": (_np.exp, _np.exp),
    "log": (_np.log, lambda x: 1.0 / x),
    "sin": (_np.sin, _np.cos),
    "cos": (_np.cos, lambda x: -_np.sin(x)),
    "tan": (_np.tan, lambda x: 1.0 / _np.cos(x) ** 2),
    "asin": (_np.arcsin, lambda x: 1.0 / _np.sqrt(1.0 - x * x)),
    "acos": (_np.arccos, lambda x: -1.0 / _np.sqrt(1.0 - x * x)),
    "atan": (_np.arctan, lambda x: 1.0 / (1.0 + x * x)),
    "sinh": (_np.sinh, _np.cosh),
    "cosh": (_np.cosh, _np.sinh),
    "tanh": (_np.tanh, lambda x: 1.0 - _np.tanh(x) ** 2),
    "asinh": (_np.arcsinh, lambda x: 1.0 / _np.sqrt(x * x + 1.0)),
    "acosh": (_np.arccosh, lambda x: 1.0 / _np.sqrt(x * x - 1.0)),
    "atanh": (_np.arctanh, lambda x: 1.0 / (1.0 - x * x)),
}
