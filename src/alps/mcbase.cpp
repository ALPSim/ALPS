/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/signal.hpp>
#include <alps/mcbase.hpp>

namespace alps {

    mcbase::mcbase(parameters_type const & parms, std::size_t seed_offset)
        : parameters(parms)
        , params(parameters) // TODO: remove, deprecated!
        , random((parameters["SEED"] | 42) + seed_offset)
    {
        alps::ngs::signal::listen();
    }

    void mcbase::save(boost::filesystem::path const & filename) const {
        alps::hdf5::archive ar(filename, "w");
        replace_citations(ar, citations());
        ar["/simulation/realizations/0/clones/0"] << *this;
    }

    void mcbase::load(boost::filesystem::path const & filename) {
        alps::hdf5::archive ar(filename);
        ar["/simulation/realizations/0/clones/0"] >> *this;
    }

    bool mcbase::run(boost::function<bool ()> const & stop_callback) {
        bool stopped = false;
        while(!(stopped = stop_callback()) && fraction_completed() < 1.) {
            update();
            measure();
        }
        return !stopped;
    }

    // implement a nice keys(m) function
    mcbase::result_names_type mcbase::result_names() const {
        result_names_type names;
        for(observable_collection_type::const_iterator it = measurements.begin(); it != measurements.end(); ++it)
            names.push_back(it->first);
        return names;
    }

    mcbase::result_names_type mcbase::unsaved_result_names() const {
        return result_names_type(); 
    }

    mcbase::results_type mcbase::collect_results() const {
        return collect_results(result_names());
    }

    mcbase::results_type mcbase::collect_results(result_names_type const & names) const {
        results_type partial_results;
        for(result_names_type::const_iterator it = names.begin(); it != names.end(); ++it)
            #ifdef ALPS_NGS_USE_NEW_ALEA
                partial_results.insert(*it, measurements[*it].result());
            #else
                partial_results.insert(*it, alps::mcresult(measurements[*it]));
            #endif
        return partial_results;
    }

    void mcbase::set_citation_component(const std::string& component, const std::string& activity) {
        make_citation_snapshot(component, activity);
        citation_component_ = component; citation_activity_ = activity;
    }
    void mcbase::inherit_citations(const citation_history& history) {
        merge_citations(citation_history_, history);
    }
    citation_history mcbase::citations() const {
        citation_history result = citation_history_;
        merge_citations(result, citation_history{make_citation_snapshot(citation_component_, citation_activity_)});
        return result;
    }
    void mcbase::save(alps::hdf5::archive & ar) const {
        write_citations(ar, citations());
        ar["/parameters"] << parameters;
        ar["measurements"] << measurements;
        ar["checkpoint/engine"] << random;
    }

    void mcbase::load(alps::hdf5::archive & ar) {
        citation_history_ = read_citations(ar);
        ar["/parameters"] >> parameters;
        ar["measurements"] >> measurements;
        ar["checkpoint/engine"] >> random;
    }

}
