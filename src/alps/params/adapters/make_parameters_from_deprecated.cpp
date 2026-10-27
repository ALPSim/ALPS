// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/ngs/make_parameters_from_deprecated.hpp>

namespace alps {

    params make_parameters_from_deprecated(Parameters const & parameters) {
        params res;
        for (Parameters::const_iterator it = parameters.begin(); it != parameters.end(); ++it)
            res.supply(it->key(), it->value());
        return res;
    }

}
