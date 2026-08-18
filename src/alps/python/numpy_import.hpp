/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2016 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Jan Gukelberger <j.gukelberger@usherbrooke.ca>     *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_PYTHON_NUMPY_IMPORT_HPP
#define ALPS_PYTHON_NUMPY_IMPORT_HPP

#include <alps/config.h>

#if defined(ALPS_HAVE_BOOST_NUMPY)
    #include <boost/python/numpy.hpp>
#else
    #include <boost/python/numeric.hpp>
#endif

// Allow callers to pin an earlier API version; default to the latest we have tested.
#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_11_API_VERSION
#endif
#include <numpy/arrayobject.h>

namespace alps {
    namespace {

        // Initialize numpy. 
        // This function has to be called from each translation unit before any function from the 
        // numpy C API is used. This function must reside in an anonymous namespace in order to 
        // ensure that it has internal linkage and that each translation unit ends up with its own
        // import_numpy function.
        //
        // Some resources explaining the numpy madness can be found at the following URLs.
        // Synopsis: The numpy API consists of macros that call functions trough a static dispatch
        // table. This table needs to be set up by a call to import_array() in each translation
        // unit lest the numpy calls segfault.
        // https://docs.scipy.org/doc/numpy/reference/c-api.array.html#miscellaneous
        // http://stackoverflow.com/a/31973355
        // https://sourceforge.net/p/numpy/mailman/message/5700519/
        void import_numpy() {
            static bool inited = false;
            if (!inited) {
                import_array1((void)0);
                #if defined(ALPS_HAVE_BOOST_NUMPY)
                boost::python::numpy::initialize();
                #else
                boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
                #endif
                inited = true;
            }
        }
    }
}

#endif
