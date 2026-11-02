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
    // Typed text values are read by make_parameters_from_file: expressions
    // such as U/2 are evaluated, other values stay text.
    //
    // An existing file is rewritten in its own format, replacing all
    // previous parameters; a new file is HDF5 when its name ends in ".h5".
    // Writing text throws if a value cannot be read back as text.

    ALPS_DECL Parameters read_parameter_file(boost::filesystem::path const & filename);

    ALPS_DECL params read_params_file(boost::filesystem::path const & filename);

    ALPS_DECL void write_parameter_file(boost::filesystem::path const & filename,
                                        Parameters const & parameters);

}

#endif
