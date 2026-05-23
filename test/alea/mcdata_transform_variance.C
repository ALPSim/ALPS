/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* Copyright (C) 1994-2025 by the ALPS collaboration
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
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

    // transform_linear: boost::none must clear the variance
    {
        alps::alea::mcdata<double> data(obs);
        if (data.has_variance()) {
            data.transform_linear([](double x){ return 2.0 * x; }, 1.5, boost::none);
            ok &= check(!data.has_variance(),
                        "transform_linear: has_variance() should be false after supplying boost::none");
        }
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

    return ok ? 0 : 1;
}
