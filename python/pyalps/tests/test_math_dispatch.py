"""Every pyalps.math dispatcher reaches a function that exists for its input.

numpy>=1.26 is supported, so the array branch must not rely on NumPy 2 aliases.
"""

import math

import numpy as np
import pytest

from pyalps import math as pam
from pyalps.alea import MCScalarData

EXPECTED = {
    "sq": lambda x: x * x,
    "sqrt": math.sqrt,
    "cb": lambda x: x * x * x,
    "cbrt": lambda x: math.copysign(abs(x) ** (1.0 / 3.0), x),
    "exp": math.exp,
    "log": math.log,
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "asin": math.asin,
    "acos": math.acos,
    "atan": math.atan,
    "sinh": math.sinh,
    "cosh": math.cosh,
    "tanh": math.tanh,
    "asinh": math.asinh,
    "acosh": math.acosh,
    "atanh": math.atanh,
}

# acosh needs x >= 1, the others accept a value in (0, 1).
ARGUMENT = {"acosh": 1.5}


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_scalar_and_array_dispatch(name):
    x = ARGUMENT.get(name, 0.5)
    function = getattr(pam, name)
    assert function(x) == pytest.approx(EXPECTED[name](x))
    np.testing.assert_allclose(function(np.array([x, x])), [EXPECTED[name](x)] * 2)


@pytest.mark.parametrize("name", ["sq", "cb"])
def test_integer_dispatch(name):
    assert getattr(pam, name)(2) == EXPECTED[name](2)


def test_cbrt_of_negative_scalar():
    assert pam.cbrt(-8.0) == pytest.approx(-2.0)


@pytest.mark.parametrize("name", ["sq", "cb", "sqrt"])
def test_observable_dispatch_uses_method(name):
    value = MCScalarData(4.0, 0.1)
    assert getattr(pam, name)(value).mean == pytest.approx(getattr(value, name)().mean)
