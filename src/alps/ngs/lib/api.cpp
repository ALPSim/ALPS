/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/api.hpp>
#include <alps/hdf5/archive.hpp>

#include <boost/filesystem.hpp>

namespace alps {

    namespace detail {
        template<typename R, typename P> void save_results_impl(R const & results, P const & params, boost::filesystem::path const & filename, std::string const & path) {
            if (results.size()) {
                hdf5::archive ar(filename.string(), "w");
                ar["/parameters"] << params;
                ar[path] << results;
            }
        }
    }

    #ifdef ALPS_NGS_USE_NEW_ALEA

        void save_results(alps::accumulator::result_set const & results, params const & params, boost::filesystem::path const & filename, std::string const & path) {
            detail::save_results_impl(results, params, filename, path);
        }

        void save_results(alps::accumulator::accumulator_set const & observables, params const & params, boost::filesystem::path const & filename, std::string const & path) {
            detail::save_results_impl(observables, params, filename, path);
        }

    #endif

    void save_results(mcresults const & results, params const & params, boost::filesystem::path const & filename, std::string const & path) {
        detail::save_results_impl(results, params, filename, path);
    }

    void save_results(mcobservables const & observables, params const & params, boost::filesystem::path const & filename, std::string const & path) {
        detail::save_results_impl(observables, params, filename, path);
    }

}
