// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_NGS_MAKE_PARAMETERS_FROM_DEPRECATED_HPP
#define ALPS_NGS_MAKE_PARAMETERS_FROM_DEPRECATED_HPP

#include <alps/params.hpp>
#include <alps/config.h>
#include <alps/parameter.h>

namespace alps {

    /// An alps::Parameters as params, for programs that read the older Parameters
    /** Each parameter is supplied as its text, as make_parameters_from_file()
        supplies it; `define<T>()` parses it. The inverse of
        make_deprecated_parameters(). */
    ALPS_DECL params make_parameters_from_deprecated(Parameters const & parameters);

}

#endif
