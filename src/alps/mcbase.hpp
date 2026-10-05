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

#ifndef ALPS_NGS_MCBASE_HPP
#define ALPS_NGS_MCBASE_HPP

#include <alps/ngs.hpp>

#include <boost/function.hpp>
#include <boost/filesystem/path.hpp>

#include <vector>
#include <string>

namespace alps {

    class ALPS_DECL mcbase {

        public:

            #ifdef ALPS_NGS_USE_NEW_ALEA
                typedef alps::accumulator::accumulator_set observable_collection_type;
            #else
                typedef alps::mcobservables observable_collection_type;
            #endif

            typedef alps::params parameters_type;
            typedef std::vector<std::string> result_names_type;

            #ifdef ALPS_NGS_USE_NEW_ALEA
                typedef alps::accumulator::result_set results_type;
            #else
                typedef alps::mcresults results_type;
            #endif

            mcbase(parameters_type const & parms, std::size_t seed_offset = 0);
            virtual ~mcbase();

            virtual void update() = 0;
            virtual void measure() = 0;
            virtual double fraction_completed() const = 0;
            bool run(boost::function<bool ()> const & stop_callback);

            result_names_type result_names() const;
            result_names_type unsaved_result_names() const;
            results_type collect_results() const;
            results_type collect_results(result_names_type const & names) const;

            void save(boost::filesystem::path const & filename) const;
            void load(boost::filesystem::path const & filename);
            virtual void save(alps::hdf5::archive & ar) const;
            virtual void load(alps::hdf5::archive & ar);

            // Non-virtual accessors for language bindings and downstream
            // exporters. Keeping these on the actual base class avoids
            // assuming that every derived simulation is a Python trampoline.
            alps::random01 & get_random() { return random; }
            parameters_type & get_parameters() { return parameters; }
            observable_collection_type & get_measurements() { return measurements; }

        protected:

            parameters_type parameters;
            parameters_type & params; // TODO: deprecated, remove!
            alps::random01 mutable random;
            observable_collection_type measurements;
    };

}

#endif
