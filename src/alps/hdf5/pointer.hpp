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

#ifndef ALPS_NGS_HDF5_POINTER_HPP
#define ALPS_NGS_HDF5_POINTER_HPP

#include <alps/hdf5/pair.hpp>

namespace alps {

    template <typename T> hdf5::detail::make_pvp_proxy<std::pair<T *, std::vector<std::size_t> > > make_pvp(
          std::string const & path
        , T * value
        , std::size_t size
    ) {
        return hdf5::detail::make_pvp_proxy<std::pair<T *, std::vector<std::size_t> > >(
              path
            , std::make_pair(value, size > 0 
                ? std::vector<std::size_t>(1, size)
                : std::vector<std::size_t>()
            )
        );
    }

    template <typename T> hdf5::detail::make_pvp_proxy<std::pair<T *, std::vector<std::size_t> > > make_pvp(
          std::string const & path
        , T * value
        , std::vector<std::size_t> const & size
    ) {
        return hdf5::detail::make_pvp_proxy<std::pair<T *, std::vector<std::size_t> > >(
              path
            , std::make_pair(value, size)
        );
    }

}

#endif
