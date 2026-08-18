# ****************************************************************************
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# ALPS Libraries
#
# Copyright (C) 2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch> ,
#                       Lukas Gamper <gamperl@gmail.com>
#                       Matthias Troyer <troyer@itp.phys.ethz.ch>
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
#
# ****************************************************************************

# Assertion-based MCScalarData / MCVectorData arithmetic test. The
# expected values are the ones recorded in the historic mcdata.output
# fixture (error propagation without covariance).

from pyalps.alea import *
import numpy as np


def assert_scalar(value, mean, error):
    assert np.isclose(value.mean, mean, rtol=1e-9), (value.mean, mean)
    assert np.isclose(value.error, error, rtol=1e-9), (value.error, error)


def assert_vector(value, means, errors):
    np.testing.assert_allclose(value.mean, means, rtol=1e-9)
    np.testing.assert_allclose(value.error, errors, rtol=1e-9)


def test_mcdata_scalar():
    b = MCScalarData(1.21, 0.15)
    c = MCScalarData(-1.5, 0.2)

    a = MCScalarData(0.81, 0.1)
    a += b
    assert_scalar(a, 2.02, 0.180277563773)

    a = MCScalarData(1.2, 0.1)
    a -= b
    assert_scalar(a, -0.01, 0.180277563773)

    a = MCScalarData(1.2, 0.1)
    a *= b
    assert_scalar(a, 1.452, 0.216889372723)

    a = MCScalarData(1.2, 0.1)
    a /= b
    assert_scalar(a, 0.991735537190, 0.148138359895)

    a = MCScalarData(1.2, 0.1)
    a += 2.0
    assert_scalar(a, 3.2, 0.1)

    a = MCScalarData(1.2, 0.1)
    a -= 2.0
    assert_scalar(a, -0.8, 0.1)

    a = MCScalarData(1.2, 0.1)
    a *= 2.0
    assert_scalar(a, 2.4, 0.2)

    a = MCScalarData(1.2, 0.1)
    a /= 2.0
    assert_scalar(a, 0.6, 0.05)

    a = MCScalarData(1.2, 0.1)
    assert_scalar(a + b, 2.41, 0.180277563773)
    assert_scalar(a - b, -0.01, 0.180277563773)
    assert_scalar(a * b, 1.452, 0.216889372723)
    assert_scalar(a / b, 0.991735537190, 0.148138359895)
    assert_scalar(a + 2.0, 3.2, 0.1)
    assert_scalar(a - 2.0, -0.8, 0.1)
    assert_scalar(a * 2.0, 2.4, 0.2)
    assert_scalar(a / 2.0, 0.6, 0.05)
    assert_scalar(2.0 / a, 1.666666666667, 0.138888888889)

    # NOTE: documents a long-standing libalps bug, present in the old
    # Boost.Python build too (the historic fixture also shows +1.2):
    # mcdata<T>::operator-() (src/alps/alea/mcdata.hpp) negates a copy
    # and returns *this unchanged, so unary minus is a no-op. When the
    # C++ operator is fixed, flip these expectations to -1.2 / negated
    # means.
    assert_scalar(-a, 1.2, 0.1)
    assert_scalar(abs(c), 1.5, 0.2)

    assert_scalar(pow(a, 2.71), 1.639008390308, 0.370142728145)
    assert_scalar(a.sq(), 1.44, 0.24)
    assert_scalar(a.sqrt(), 1.095445115010, 0.045643546459)
    assert_scalar(a.cb(), 1.728, 0.432)
    assert_scalar(a.cbrt(), 1.062658569183, 0.029518293588)
    assert_scalar(a.exp(), 3.320116922737, 0.332011692274)
    assert_scalar(a.log(), 0.182321556794, 0.083333333333)
    assert_scalar(a.sin(), 0.932039085967, 0.036235775448)
    assert_scalar(a.cos(), 0.362357754477, 0.093203908597)
    assert_scalar(a.tan(), 2.572151622126, 0.761596396721)
    assert_scalar(a.tanh(), 0.833654607012, 0.030501999621)


def test_mcdata_vector():
    X = MCVectorData(np.array([2.3, 1.2, 0.7]), np.array([0.01, 0.01, 0.01]))
    Y = X + 1.0

    assert_vector(X, [2.3, 1.2, 0.7], [0.01] * 3)
    assert_vector(Y, [3.3, 2.2, 1.7], [0.01] * 3)

    assert_vector(X + Y, [5.6, 3.4, 2.4], [0.014142135624] * 3)
    assert_vector(X + 2.0, [4.3, 3.2, 2.7], [0.01] * 3)
    assert_vector(2.0 + X, [4.3, 3.2, 2.7], [0.01] * 3)

    assert_vector(X / Y,
                  [0.696969696970, 0.545454545455, 0.411764705882],
                  [0.003693697954, 0.005177671110, 0.006361514294])
    assert_vector(X / 2.0, [1.15, 0.6, 0.35], [0.005] * 3)
    assert_vector(2.0 / X,
                  [0.869565217391, 1.666666666667, 2.857142857143],
                  [0.003780718336, 0.013888888889, 0.040816326531])

    # unary minus is a no-op — same libalps mcdata bug as in the scalar
    # test above; flip to negated means once the C++ operator is fixed
    assert_vector(-X, [2.3, 1.2, 0.7], [0.01] * 3)
    assert_vector(abs(X), [2.3, 1.2, 0.7], [0.01] * 3)

    assert_vector(pow(X, 2.71),
                  [9.556138502711, 1.639008390308, 0.380378260851],
                  [0.112596240619, 0.037014272814, 0.014726072670])
    assert_vector(X.sq(), [5.29, 1.44, 0.49], [0.046, 0.024, 0.014])
    assert_vector(X.sqrt(),
                  [1.516575088810, 1.095445115010, 0.836660026534],
                  [0.003296902367, 0.004564354646, 0.005976143047])
    assert_vector(X.cb(), [12.167, 1.728, 0.343], [0.1587, 0.0432, 0.0147])
    assert_vector(X.cbrt(),
                  [1.320006121796, 1.062658569183, 0.887904001743],
                  [0.001913052350, 0.002951829359, 0.004228114294])
    assert_vector(X.exp(),
                  [9.974182454815, 3.320116922737, 2.013752707470],
                  [0.099741824548, 0.033201169227, 0.020137527075])
    assert_vector(X.log(),
                  [0.832909122935, 0.182321556794, -0.356674943939],
                  [0.004347826087, 0.008333333333, 0.014285714286])
    assert_vector(X.sin(),
                  [0.745705212177, 0.932039085967, 0.644217687238],
                  [0.006662760213, 0.003623577545, 0.007648421873])
    assert_vector(X.cos(),
                  [-0.666276021280, 0.362357754477, 0.764842187284],
                  [0.007457052122, 0.009320390860, 0.006442176872])
    assert_vector(X.tan(),
                  [-1.119213641734, 2.572151622126, 0.842288380463],
                  [0.022526391758, 0.076159639672, 0.017094497159])
    assert_vector(X.sinh(),
                  [4.936961805546, 1.509461355412, 0.758583701840],
                  [0.050372206493, 0.018106555673, 0.012551690056])
    assert_vector(X.cosh(),
                  [5.037220649269, 1.810655567324, 1.255169005631],
                  [0.049369618055, 0.015094613554, 0.007585837018])
    assert_vector(X.tanh(),
                  [0.980096396266, 0.833654607012, 0.604367777117],
                  [0.000394110540, 0.003050199962, 0.006347395900])


if __name__ == "__main__":
    test_mcdata_scalar()
    test_mcdata_vector()
    print("SUCCESS")
