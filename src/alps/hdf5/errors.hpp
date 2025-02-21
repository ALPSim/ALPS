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

#ifndef ALPS_NGS_HDF5_ERROR_HPP
#define ALPS_NGS_HDF5_ERROR_HPP

#include <string>
#include <stdexcept>

namespace alps {
    namespace hdf5 {

        class archive_error : public std::runtime_error {
            public:
                archive_error(std::string const & what)
                    : std::runtime_error(what) 
                {}
        };

        #define DEFINE_ALPS_HDF5_EXCEPTION(name)                                    \
            class name : public archive_error {                                     \
                public:                                                             \
                    name (std::string const & what)                                 \
                        : archive_error(what)                                       \
                    {}                                                              \
            };
        DEFINE_ALPS_HDF5_EXCEPTION(archive_not_found)
        DEFINE_ALPS_HDF5_EXCEPTION(archive_closed)
        DEFINE_ALPS_HDF5_EXCEPTION(invalid_path)
        DEFINE_ALPS_HDF5_EXCEPTION(path_not_found)
        DEFINE_ALPS_HDF5_EXCEPTION(wrong_type)
        #undef DEFINE_ALPS_HDF5_EXCEPTION
    }
};

#endif
