// Copyright (C) 2010 by Matthias Troyer
// SPDX-License-Identifier: MIT
#ifndef ALPS_PYTHON_SAVE_OBSERVABLE_TO_HDF5_HPP
#define ALPS_PYTHON_SAVE_OBSERVABLE_TO_HDF5_HPP

#include <alps/hdf5.hpp>

#include <string>

namespace alps {
namespace python {

// Implementation of the Python observable.save(filename) convenience method.
template <typename Observable>
void save_observable_to_hdf5(Observable const & observable,
                             std::string const & filename) {
    hdf5::archive archive(filename, "a");
    archive["/simulation/results/" + observable.representation()] << observable;
}

}  // namespace python
}  // namespace alps

#endif
