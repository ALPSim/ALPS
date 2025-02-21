/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>                          *
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

#ifndef ALPS_HDF5_NUMERIC_VECTOR_HPP
#define ALPS_HDF5_NUMERIC_VECTOR_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/numeric/matrix/vector.hpp>

namespace alps {
namespace hdf5 {

        template <typename T, typename MemoryBlock>
        void save(
                  alps::hdf5::archive & ar
                  , std::string const & path
                  , alps::numeric::vector<T, MemoryBlock> const & value
                  , std::vector<std::size_t> size = std::vector<std::size_t>()
                  , std::vector<std::size_t> chunk = std::vector<std::size_t>()
                  , std::vector<std::size_t> offset = std::vector<std::size_t>()
                  ) {
            ar[path] << static_cast<MemoryBlock const&>(value);
        }
        template <typename T, typename MemoryBlock>
        void load(
                  alps::hdf5::archive & ar
                  , std::string const & path
                  , alps::numeric::vector<T, MemoryBlock> & value
                  , std::vector<std::size_t> chunk = std::vector<std::size_t>()
                  , std::vector<std::size_t> offset = std::vector<std::size_t>()
                  ) {
            MemoryBlock tmp;
            ar[path] >> tmp;
            value = alps::numeric::vector<T, MemoryBlock>(tmp.begin(), tmp.end());
        }
}
}
#endif // ALPS_HDF5_NUMERIC_VECTOR_HPP
