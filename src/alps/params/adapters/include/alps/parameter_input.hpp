// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_PARAMETER_INPUT_HPP
#define ALPS_PARAMETER_INPUT_HPP

#include <alps/params.hpp>

#include <alps/parameter.h>
#include <alps/config.h>

#include <boost/filesystem/path.hpp>

namespace alps {

    // Applications read and write parameter files only through these
    // functions; the file format is chosen here and nowhere else.
    //
    // A file is read as HDF5 (parameters in the /parameters group) when it
    // carries the HDF5 signature, and as ALPS parameter text otherwise.
    // A file is written as HDF5 when its name ends in ".h5" or it already is
    // an HDF5 file, and as ALPS parameter text otherwise.

    ALPS_DECL Parameters read_parameter_file(boost::filesystem::path const & filename);

    ALPS_DECL params read_params_file(boost::filesystem::path const & filename);

    ALPS_DECL void write_parameter_file(boost::filesystem::path const & filename,
                                        Parameters const & parameters);

}

#endif
