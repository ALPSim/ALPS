// Copyright (C) 2010 by Matthias Troyer <troyer@comp-phys.org>,
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
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
