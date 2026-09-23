# ****************************************************************************
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# ALPS Libraries
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
#
# ****************************************************************************

# Formerly accumulators.py: a print-only script importing the extension by its
# flat module name, which only resolved when tests ran with PYTHONPATH pointed
# at the old in-tree build directory. pytest never collected it (the name
# matched neither test_*.py nor *_test.py), so it had gone unrun long enough
# for str() on every accumulator to regress unnoticed.

import numpy as np
import pytest

from pyalps.cxx import pyngsaccumulator_c as accumulator


ACCUMULATORS = (
    "count_accumulator",
    "mean_accumulator",
    "error_accumulator",
    "binning_analysis_accumulator",
    "max_num_binning_accumulator",
)


def feed(name, samples):
    """Build an accumulator of the named kind and feed it `samples`."""
    accum = getattr(accumulator, name)()
    for sample in samples:
        accum(float(sample))
    return accum


@pytest.mark.parametrize("name", ACCUMULATORS)
def test_accumulator_and_result_are_printable(name):
    # __str__ is bound as a std::string-returning function; without
    # nanobind/stl/string.h in that translation unit it raised TypeError
    # instead of printing.
    accum = feed(name, range(64))
    assert isinstance(str(accum), str)
    assert str(accum) != ""
    assert isinstance(str(accum.result()), str)
    assert str(accum.result()) != ""


@pytest.mark.parametrize("name", ACCUMULATORS)
def test_count_survives_the_accumulator_to_result_transition(name):
    accum = feed(name, range(37))
    assert accum.count() == 37
    assert accum.result().count() == 37


@pytest.mark.parametrize(
    "name", ("mean_accumulator", "error_accumulator",
             "binning_analysis_accumulator", "max_num_binning_accumulator"))
def test_mean_matches_numpy(name):
    samples = [1.5, 2.5, 3.5, 4.5]
    accum = feed(name, samples)
    assert accum.mean() == pytest.approx(np.mean(samples))
    assert accum.result().mean() == pytest.approx(np.mean(samples))


def test_error_of_two_samples():
    # mean 10, sample standard error of {8, 12} is 2.
    result = feed("error_accumulator", (8, 12)).result()
    assert result.mean() == pytest.approx(10.0)
    assert result.error() == pytest.approx(2.0)


def test_binning_analysis_reports_a_positive_error():
    for name in ("binning_analysis_accumulator", "max_num_binning_accumulator"):
        result = feed(name, np.linspace(0.0, 1.0, 1024)).result()
        assert result.mean() == pytest.approx(0.5, abs=1e-9)
        assert result.error() > 0.0


def test_result_arithmetic_against_scalars():
    result = feed("error_accumulator", (8, 12)).result()

    assert (result * 2).mean() == pytest.approx(20.0)
    assert (result + 2).mean() == pytest.approx(12.0)
    assert (2 + result).mean() == pytest.approx(12.0)
    assert (result - 2).mean() == pytest.approx(8.0)
    assert (2 - result).mean() == pytest.approx(-8.0)
    assert (result / 2).mean() == pytest.approx(5.0)
    assert (-result).mean() == pytest.approx(-10.0)

    # Non-mutating: the operators above take their operand by value.
    assert result.mean() == pytest.approx(10.0)


def test_result_arithmetic_against_results():
    left = feed("error_accumulator", (8, 12)).result()
    right = feed("error_accumulator", (3, 7)).result()

    assert (left + right).mean() == pytest.approx(15.0)
    assert (left - right).mean() == pytest.approx(5.0)
    assert (left * right).mean() == pytest.approx(50.0)
    assert (left / right).mean() == pytest.approx(2.0)


def test_inplace_operators_mutate_and_return_the_same_object():
    result = feed("error_accumulator", (8, 12)).result()
    identity = id(result)

    result *= 2
    assert result.mean() == pytest.approx(20.0)
    result += 5
    assert result.mean() == pytest.approx(25.0)
    result -= 5
    assert result.mean() == pytest.approx(20.0)
    result /= 2
    assert result.mean() == pytest.approx(10.0)

    # rv_policy::none must hand back the original wrapper, not a copy.
    assert id(result) == identity


def test_transcendental_functions_match_numpy():
    result = feed("mean_accumulator", (0.25,)).result()

    for name in ("sin", "cos", "tan", "sinh", "cosh", "tanh", "sqrt", "log", "abs"):
        expected = getattr(np, name if name != "abs" else "fabs")(0.25)
        assert getattr(result, name)().mean() == pytest.approx(expected), name

    # numpy ufuncs route through the bound methods for a scalar-like object.
    assert np.sin(result).mean() == pytest.approx(np.sin(0.25))


def test_reset_clears_the_samples():
    accum = feed("mean_accumulator", range(10))
    assert accum.count() == 10
    accum.reset()
    assert accum.count() == 0
