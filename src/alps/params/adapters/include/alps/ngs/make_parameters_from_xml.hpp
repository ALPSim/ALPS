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

#ifndef ALPS_NGS_MAKE_PARAMETERS_FROM_XML_HPP
#define ALPS_NGS_MAKE_PARAMETERS_FROM_XML_HPP

#include <alps/params.hpp>
#include <boost/filesystem/path.hpp>

#include <alps/parameter.h>
#include <alps/config.h>

namespace alps {

    /// Parameters supplied as text from the <PARAMETERS> of an XML job file
    /** SEED is supplied as 0 if the file does not set it. Expressions such
        as `J/4` are supplied as their value, other values as their text. */
    ALPS_DECL params make_parameters_from_xml(boost::filesystem::path const & filename);

}

#endif
