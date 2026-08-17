/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>

namespace alps {
    namespace detail {

        void save_results_export(mcresults const & res, params const & par, alps::hdf5::archive & ar, std::string const & path) {
            ar["/parameters"] << par;
            if (res.size())
                ar[path] << res;
        }
    }
}

BOOST_PYTHON_MODULE(pyngsapi_c) {

    boost::python::def("collectResults", static_cast<alps::results_type<alps::mcbase>::type (*)(alps::mcbase const &)>(&alps::collect_results<alps::mcbase>));

    boost::python::def("saveResults", &alps::detail::save_results_export);

}
