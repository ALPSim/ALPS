/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

#include <alps/mc/api.hpp>
#include <alps/hdf5/archive.hpp>

namespace alps {

    void save_results(alps::alea::result_set const & results, params const & params, std::string const & filename, std::string const & path) {
        if (results.names().empty())
            return;
        hdf5::archive ar(filename, "w");
        ar["/parameters"] << params;
        ar[path] << results;
    }

    void save_results(alps::alea::observable_set const & observables, params const & params, std::string const & filename, std::string const & path) {
        save_results(observables.results(), params, filename, path);
    }

}
