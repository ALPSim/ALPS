"""Regression checks found by comparing the Boost.Python conversion contracts."""

import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from pyalps import alea, ngs
from pyalps.cxx.pyngsaccumulator_c import count_accumulator


def test_count_accumulator_accepts_non_scalar_samples():
    accumulator = count_accumulator()
    for sample in (np.ones((2, 3)), 1 + 2j, [1, 2], {"x": 1}, None):
        accumulator(sample)
    assert accumulator.count() == 5
    assert accumulator.result().count() == 5


@pytest.mark.parametrize("vector", [False, True])
def test_timeseries_from_mcdata(tmp_path, vector):
    observable_type = alea.RealVectorTimeSeriesObservable if vector else alea.RealTimeSeriesObservable
    data_type = alea.MCVectorData if vector else alea.MCScalarData
    series_type = alea.MCVectorTimeseries if vector else alea.MCScalarTimeseries
    observable = observable_type("samples")
    for i in range(32):
        observable << (np.array([float(i), 2.0 * i]) if vector else float(i))
    filename = str(tmp_path / "samples.h5")
    observable.save(filename)
    data = data_type()
    data.load(filename, "/simulation/results/samples")
    np.testing.assert_array_equal(series_type(data).timeseries(), data.bins)


@pytest.mark.parametrize("observable_type", [alea.RealVectorObservable, alea.RealVectorTimeSeriesObservable])
def test_vector_convergence_is_an_integer_array(observable_type):
    observable = observable_type("samples")
    for i in range(128):
        observable << np.array([float(i), 2.0 * i])
    convergence = observable.converged_errors
    assert isinstance(convergence, np.ndarray)
    assert convergence.shape == (2,)
    assert convergence.dtype.kind == "i"
    assert np.isin(convergence, [0, 1, 2]).all()


@pytest.mark.parametrize("layout", ["readonly", "strided", "fortran"])
def test_array_consumers_do_not_require_writable_samples(layout):
    values = np.arange(24.0).reshape(4, 6)
    if layout == "readonly":
        values.flags.writeable = False
    elif layout == "strided":
        values = values[:, ::2]
    else:
        values = np.asfortranarray(values)
    np.testing.assert_array_equal(alea.MCScalarTimeseries(values[0]).timeseries(), values[0])
    np.testing.assert_array_equal(alea.MCVectorTimeseries(values).timeseries(), values)
    observable = alea.RealVectorObservable("samples")
    ngs_observable = ngs.createRealVectorObservable("samples")
    for row in values:
        observable << row
        ngs_observable << row
    np.testing.assert_allclose(observable.mean, values.mean(axis=0))
    np.testing.assert_allclose(ngs.observable2result(ngs_observable).mean, values.mean(axis=0))


def test_mcvector_constructor_sizes_errors_before_indexing():
    # Printing/indexing a mean-only vector previously read an empty error
    # vector and crashed the interpreter. Keep native crash checks isolated.
    code = """
        import numpy as np
        from pyalps.alea import MCVectorData
        data = MCVectorData([1.0, 2.0, 3.0])
        np.testing.assert_array_equal(data.error, [0.0, 0.0, 0.0])
        assert data[1].mean == 2.0
        assert data[1].error == 0.0
        assert '2' in repr(data)
        assert '2.00' in format(data, '.2f')
        try:
            MCVectorData([1.0, 2.0], [0.1])
        except ValueError:
            pass
        else:
            raise AssertionError('mismatched mean/error lengths must be rejected')
    """
    result = subprocess.run(
        [sys.executable, "-X", "faulthandler", "-c", textwrap.dedent(code)],
        capture_output=True, text=True, timeout=30,
        env={**os.environ, "MallocScribble": "1"},
    )
    assert result.returncode == 0, result.stdout + result.stderr
