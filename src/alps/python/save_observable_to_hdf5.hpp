/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2010 by Matthias Troyer <troyer@comp-phys.org>,
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

/* $Id: make_copy.hpp 4059 2010-03-29 08:36:25Z troyer $ */

#ifndef ALPS_PYTHON_VERY_LONG_FILENAME_FOR_SAVE_OBSERVABLE_TO_HDF5_HPP
#define ALPS_PYTHON_VERY_LONG_FILENAME_FOR_SAVE_OBSERVABLE_TO_HDF5_HPP

#include <alps/hdf5.hpp>

namespace alps { namespace python {
    
    template <typename Obs> void save_observable_to_hdf5(Obs const & obs, std::string const & filename) {
        hdf5::archive ar(filename, "a");
        ar["/simulation/results/"+obs.representation()] << obs;
    }
        
} } // end namespace alps::python

#endif // ALPS_PYTHON_VERY_LONG_FILENAME_FOR_SAVE_OBSERVABLE_TO_HDF5_HPP
