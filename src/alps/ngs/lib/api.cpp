/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/api.hpp>
#include <alps/hdf5/archive.hpp>

#include <boost/filesystem.hpp>

namespace alps {

    namespace detail {
        template<typename R, typename P> void save_results_impl(R const & results, P const & params, boost::filesystem::path const & filename, std::string const & path, citation_history const & citations) {
            if (results.size()) {
                hdf5::archive ar(filename.string(), "w");
                replace_citations(ar, citations);
                ar["/parameters"] << params;
                ar[path] << results;
            }
        }
    }

    #ifdef ALPS_NGS_USE_NEW_ALEA

        void save_results(alps::accumulator::result_set const & results, params const & params, boost::filesystem::path const & filename, std::string const & path) {
            save_results(results, params, filename, path, citation_history{make_citation_snapshot("framework", "unspecified")});
        }

        void save_results(alps::accumulator::result_set const & results, params const & params, boost::filesystem::path const & filename, std::string const & path, citation_history const & citations) {
            detail::save_results_impl(results, params, filename, path, citations);
        }

        void save_results(alps::accumulator::accumulator_set const & observables, params const & params, boost::filesystem::path const & filename, std::string const & path) {
            save_results(observables, params, filename, path, citation_history{make_citation_snapshot("framework", "unspecified")});
        }

        void save_results(alps::accumulator::accumulator_set const & observables, params const & params, boost::filesystem::path const & filename, std::string const & path, citation_history const & citations) {
            detail::save_results_impl(observables, params, filename, path, citations);
        }

    #endif

    void save_results(mcresults const & results, params const & params, boost::filesystem::path const & filename, std::string const & path) {
        save_results(results, params, filename, path, citation_history{make_citation_snapshot("framework", "unspecified")});
    }

    void save_results(mcresults const & results, params const & params, boost::filesystem::path const & filename, std::string const & path, citation_history const & citations) {
        detail::save_results_impl(results, params, filename, path, citations);
    }

    void save_results(mcobservables const & observables, params const & params, boost::filesystem::path const & filename, std::string const & path) {
        save_results(observables, params, filename, path, citation_history{make_citation_snapshot("framework", "unspecified")});
    }

    void save_results(mcobservables const & observables, params const & params, boost::filesystem::path const & filename, std::string const & path, citation_history const & citations) {
        detail::save_results_impl(observables, params, filename, path, citations);
    }

}
