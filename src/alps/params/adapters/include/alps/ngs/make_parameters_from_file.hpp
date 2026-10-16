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

#ifndef ALPS_NGS_MAKE_PARAMETERS_FROM_FILE_HPP
#define ALPS_NGS_MAKE_PARAMETERS_FROM_FILE_HPP

#include <alps/params.hpp>
#include <boost/filesystem/path.hpp>

#include <alps/parameter.h>
#include <alps/config.h>

namespace alps {

    /// Parameters supplied as text from a classic ALPS parameter file
    /** The file is read as alps::Parameters reads it, `KEY = value` per
        line; each value is supplied as its text. */
    ALPS_DECL params make_parameters_from_file(boost::filesystem::path const & filename);

}

#endif
