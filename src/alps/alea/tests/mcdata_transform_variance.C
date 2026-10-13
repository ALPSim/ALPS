/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* Copyright (C) 1994-2025 by the ALPS collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Regression test for the variance_opt self-assignment bug in mcdata::transform_linear
// and mcdata::transform. Each method accepted a variance_opt parameter but assigned
// variance_opt_ to itself instead of to the parameter.

#include <alps/alea.h>
#include <alps/alea/mcdata.hpp>
#include <boost/optional.hpp>
#include <cmath>
#include <iostream>

static bool check(bool cond, const char* msg) {
    if (!cond) std::cerr << "FAIL: " << msg << "\n";
    return cond;
}

int main() {
    bool ok = true;

    alps::RealObservable obs("test");
    for (int i = 0; i < 1024; ++i)
        obs << double(i % 10);

    // transform_linear: supplied variance must be stored
    {
        alps::alea::mcdata<double> data(obs);
        const double expected = 42.0;
        data.transform_linear([](double x){ return 2.0 * x; }, 1.5,
                              boost::optional<double>(expected));
        ok &= check(data.has_variance(),
                    "transform_linear: has_variance() should be true after supplying variance");
        ok &= check(data.has_variance() && std::abs(data.variance() - expected) < 1e-10,
                    "transform_linear: variance() should equal the supplied value");
    }

    // transform_linear: boost::none must clear the variance (unconditional: set it first)
    {
        alps::alea::mcdata<double> data(obs);
        data.transform_linear([](double x){ return x; }, 1.5, boost::optional<double>(42.0));
        data.transform_linear([](double x){ return x; }, 1.5, boost::none);
        ok &= check(!data.has_variance(),
                    "transform_linear: has_variance() should be false after supplying boost::none");
    }

    // transform (single-observable overload): supplied variance must be stored
    {
        alps::alea::mcdata<double> data(obs);
        const double expected = 99.0;
        data.transform([](double x){ return x + 1.0; }, 2.0,
                       boost::optional<double>(expected));
        ok &= check(data.has_variance(),
                    "transform: has_variance() should be true after supplying variance");
        ok &= check(data.has_variance() && std::abs(data.variance() - expected) < 1e-10,
                    "transform: variance() should equal the supplied value");
    }

    // transform (binary overload): supplied variance must be stored
    {
        alps::alea::mcdata<double> lhs(obs);
        alps::alea::mcdata<double> rhs(obs);
        const double expected = 7.0;
        lhs.transform(rhs, [](double x, double y){ return x + y; }, 3.0,
                      boost::optional<double>(expected));
        ok &= check(lhs.has_variance(),
                    "binary transform: has_variance() should be true after supplying variance");
        ok &= check(lhs.has_variance() && std::abs(lhs.variance() - expected) < 1e-10,
                    "binary transform: variance() should equal the supplied value");
    }

    return ok ? 0 : 1;
}
