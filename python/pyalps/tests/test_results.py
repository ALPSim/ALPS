"""Arithmetic and functions of pyalps.ngs results, which wrap ALPSCore's.

These checks replace those of the retired pyalps.pyngsaccumulator_c on the NGS
accumulators. Constant samples give an independent analytic value.
"""

import numpy as np
import pytest

from pyalps import ngs


def constant(value, count=64):
    observable = (ngs.createRealVectorObservable("x") if np.ndim(value)
                  else ngs.createRealObservable("x"))
    for _ in range(count):
        observable << value
    return ngs.observable2result(observable)


def samples(values):
    observable = ngs.createRealObservable("x")
    for value in values:
        observable << float(value)
    return ngs.observable2result(observable)


FUNCTIONS = {
    "sin": np.sin, "cos": np.cos, "tan": np.tan, "sinh": np.sinh,
    "cosh": np.cosh, "tanh": np.tanh, "sqrt": np.sqrt, "cbrt": np.cbrt,
    "log": np.log, "exp": np.exp, "sq": np.square,
    "cb": lambda value: np.power(value, 3),
}


@pytest.mark.parametrize("value", [0.25, np.array([0.25, 0.5])])
@pytest.mark.parametrize("name", sorted(FUNCTIONS))
def test_functions_match_numpy(value, name):
    result = constant(value)
    np.testing.assert_allclose(getattr(result, name)().mean, FUNCTIONS[name](value), rtol=1e-12)
    # The function returns a new result.
    np.testing.assert_allclose(result.mean, value, rtol=1e-12)


@pytest.mark.parametrize("value", [-0.25, np.array([-0.25, 0.5])])
def test_abs_and_power(value):
    result = constant(value)
    np.testing.assert_allclose(abs(result).mean, np.abs(value), rtol=1e-12)
    np.testing.assert_allclose((result ** 2).mean, np.square(value), rtol=1e-12)


def test_arithmetic_against_scalars_and_results():
    left, right = samples([8, 12] * 32), samples([3, 7] * 32)
    assert left.mean == pytest.approx(10.0)
    assert (left * 2).mean == pytest.approx(20.0)
    assert (left + 2).mean == pytest.approx(12.0)
    assert (2 + left).mean == pytest.approx(12.0)
    assert (left - 2).mean == pytest.approx(8.0)
    assert (2 - left).mean == pytest.approx(-8.0)
    assert (left / 2).mean == pytest.approx(5.0)
    assert (-left).mean == pytest.approx(-10.0)
    assert (left + right).mean == pytest.approx(15.0)
    assert (left - right).mean == pytest.approx(5.0)
    assert left.mean == pytest.approx(10.0)


def test_inplace_operators_modify_and_return_the_same_result():
    result = samples([8, 12] * 32)
    alias = result
    result *= 2
    result += 5
    result -= 5
    result /= 2
    assert result is alias
    assert result.mean == pytest.approx(10.0)


def test_binning_analysis_reports_an_error_and_a_count():
    result = samples(np.linspace(0.0, 1.0, 1024))
    assert result.count == 1024
    assert result.mean == pytest.approx(0.5, abs=1e-9)
    assert result.error > 0.0
    assert isinstance(str(result), str) and str(result)
