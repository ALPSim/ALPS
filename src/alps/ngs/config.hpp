/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_CONFIG_HPP
#define ALPS_NGS_CONFIG_HPP

#include <alps/config.h>

// if defined, no threading libraries are included
// #define ALPS_NGS_SINGLE_THREAD

// do not throw an error on accessing a not existing paht in a hdf5 file
// #define ALPS_HDF5_READ_GREEDY

// do not throw an error if closing a hdf5 gets dirty (e.g in Python)
// #define ALPS_HDF5_CLOSE_GREEDY

// blocksize in compressed hdf5. Default: 32
#ifndef ALPS_HDF5_SZIP_BLOCK_SIZE
    #define ALPS_HDF5_SZIP_BLOCK_SIZE 32
#endif

// maximal number of stack frames displayed in stacktrace. Default 63
#ifndef ALPS_NGS_MAX_FRAMES
    #define ALPS_NGS_MAX_FRAMES 63
#endif

// prevent the signal object from registering signals
#ifdef BOOST_MSVC
    #define ALPS_NGS_NO_SIGNALS
#endif

// do not print a stacktrace in error messages
#ifndef __GNUC__
    #define ALPS_NGS_NO_STACKTRACE
#endif

// TODO: have_python
// TODO: have_mpi
// TODO: have_thread

#endif

