/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2016 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Jan Gukelberger <j.gukelberger@usherbrooke.ca>     *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_PYTHON_NUMPY_IMPORT_HPP
#define ALPS_PYTHON_NUMPY_IMPORT_HPP

#if BOOST_VERSION >= 106300
    #include <boost/python/numpy.hpp>
#else
    #include <boost/python/numeric.hpp>
#endif

// this should be set to the latest numpy version we have tested
#define NPY_NO_DEPRECATED_API NPY_1_11_API_VERSION
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
                #if BOOST_VERSION >= 106300
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
