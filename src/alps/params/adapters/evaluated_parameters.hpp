// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_PARAMS_ADAPTERS_EVALUATED_PARAMETERS_HPP
#define ALPS_PARAMS_ADAPTERS_EVALUATED_PARAMETERS_HPP

#include <alps/params.hpp>
#include <alps/parameter.h>

namespace alps {
    namespace detail {

        /// An alps::Parameters as params, with expressions evaluated
        /** A value that the expression library evaluates to a number against
            the other parameters, such as `J/4`, is supplied as that number;
            numbers and other text, such as `chain lattice`, are supplied as
            written. Used by make_parameters_from_file() and
            make_parameters_from_xml(). */
        params make_parameters_evaluated(Parameters const & parameters);

    }
}

#endif
