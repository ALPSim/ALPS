/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <vector>
#include <complex>

#include <alps/hdf5.hpp>
#include <alps/hdf5/vector.hpp>

int main() {
    alps::hdf5::archive ar("creal.h5", "a");
    {
        std::vector<double> a(1e6);
        ar << alps::make_pvp("a",
            std::make_pair(
                static_cast< double const *>(&a.front())
                , std::vector<std::size_t>(1,a.size())
            )
        );
    }
    {
        std::vector<std::complex<double> > a(1e6);
        ar << alps::make_pvp("a",
            std::make_pair(
                static_cast<std::complex<double> const *>(&a.front())
                , std::vector<std::size_t>(1,a.size())
            )
        );
    }
    return 0;
}
