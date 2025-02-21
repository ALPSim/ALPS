/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2010 by Matthias Troyer <troyer@comp-phys.org>,
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
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
