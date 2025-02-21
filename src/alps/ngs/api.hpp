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

#ifndef ALPS_NGS_API_HPP
#define ALPS_NGS_API_HPP

#include <alps/ngs/config.hpp>
#include <alps/ngs/params.hpp>
#include <alps/ngs/mcresults.hpp>
#include <alps/ngs/mcobservables.hpp>

#ifdef ALPS_NGS_USE_NEW_ALEA
    #include <alps/ngs/accumulator.hpp>
#endif

#include <boost/filesystem/path.hpp>

#include <string>

namespace alps {

    template<typename S> struct result_names_type {
        typedef typename S::result_names_type type;
    };

    template<typename S> struct results_type {
        typedef typename S::results_type type;
    };

    template<typename S> struct parameters_type {
        typedef typename S::parameters_type type;
    };

    template<typename S> typename result_names_type<S>::type result_names(S const & s) {
        return s.result_names();
    }

    template<typename S> typename result_names_type<S>::type unsaved_result_names(S const & s) {
        return s.unsaved_result_names();
    }

    template<typename S> typename results_type<S>::type collect_results(S const & s) {
        return s.collect_results();
    }

    template<typename S> typename results_type<S>::type collect_results(S const & s, typename result_names_type<S>::type const & names) {
        return s.collect_results(names);
    }

    template<typename S> typename results_type<S>::type collect_results(S const & s, std::string const & name) {
        return collect_results(s, typename result_names_type<S>::type(1, name));
    }

    template<typename S> double fraction_completed(S const & s) {
        return s.fraction_completed();
    }

    #ifdef ALPS_NGS_USE_NEW_ALEA
        ALPS_DECL void save_results(alps::accumulator::accumulator_set const & observables, params const & params, boost::filesystem::path const & filename, std::string const & path);
        ALPS_DECL void save_results(alps::accumulator::result_set const & results, params const & params, boost::filesystem::path const & filename, std::string const & path);
    #endif

    ALPS_DECL void save_results(mcresults const & results, params const & params, boost::filesystem::path const & filename, std::string const & path);

    ALPS_DECL void save_results(mcobservables const & observables, params const & params, boost::filesystem::path const & filename, std::string const & path);

    template<typename C, typename P> void broadcast(C const & c, P & p, int r = 0) {
        p.broadcast(c, r);
    }

}

#endif
