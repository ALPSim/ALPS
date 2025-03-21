/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_NGS_HDF5_STD_ARRAY_HPP
#define ALPS_NGS_HDF5_STD_ARRAY_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/cast.hpp>

#include <array>
#include <vector>
#include <iterator>
#include <algorithm>

namespace alps {
    namespace hdf5 {

        template<typename T, std::size_t N> struct scalar_type<std::array<T, N> > {
            typedef typename scalar_type<typename std::array<T, N>::value_type>::type type;
        };

        template<typename T, std::size_t N> struct is_continuous<std::array<T, N> >
            : public is_continuous<T>
        {};
        template<typename T, std::size_t N> struct is_continuous<std::array<T, N> const >
            : public is_continuous<T>
        {};

        template<typename T, std::size_t N> struct has_complex_elements<std::array<T, N> >
            : public has_complex_elements<typename alps::detail::remove_cvr<typename std::array<T, N>::value_type>::type>
        {};

        namespace detail {

            template<typename T, std::size_t N> struct get_extent<std::array<T, N> > {
                static std::vector<std::size_t> apply(std::array<T, N> const & value) {
                    using alps::hdf5::get_extent;
                    std::vector<std::size_t> result(1, value.size());
                    if (value.size()) {
                        std::vector<std::size_t> first(get_extent(value[0]));
                        std::copy(first.begin(), first.end(), std::back_inserter(result));
                    }
                    return result;
                }
            };

            template<typename T, std::size_t N> struct set_extent<std::array<T, N> > {
                static void apply(std::array<T, N> & value, std::vector<std::size_t> const & extent) {
                    using alps::hdf5::set_extent;
                    if (extent.size() > 1)
                        for(typename std::array<T, N>::iterator it = value.begin(); it != value.end(); ++it)
                            set_extent(*it, std::vector<std::size_t>(extent.begin() + 1, extent.end()));
                    else if (extent.size() == 0 && !std::is_same<typename scalar_type<T>::type, T>::value)
                        throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
                }
            };

            template<typename T, std::size_t N> struct is_vectorizable<std::array<T, N> > {
                static bool apply(std::array<T, N> const & value) {
                    using alps::hdf5::get_extent;
                    using alps::hdf5::is_vectorizable;
                    if (!is_continuous<std::array<T, N> >::value) {
                        if (!is_vectorizable(value[0]))
                            return false;
                        std::vector<std::size_t> first(get_extent(value[0]));
                        for(typename std::array<T, N>::const_iterator it = value.begin(); it != value.end(); ++it)
                            if (!is_vectorizable(*it))
                                return false;
                            else {
                                std::vector<std::size_t> size(get_extent(*it));
                                if (
                                       first.size() != size.size() 
                                    || !std::equal(first.begin(), first.end(), size.begin())
                                )
                                    return false;
                            }
                    }
                    return true;
                }
            };

            template<typename T, std::size_t N> struct get_pointer<std::array<T, N> > {
                static typename alps::hdf5::scalar_type<std::array<T, N> >::type * apply(std::array<T, N> & value) {
                    using alps::hdf5::get_pointer;
                    return get_pointer(value[0]);
                }
            };

            template<typename T, std::size_t N> struct get_pointer<std::array<T, N> const > {
                static typename alps::hdf5::scalar_type<std::array<T, N> >::type const * apply(std::array<T, N> const & value) {
                    using alps::hdf5::get_pointer;
                    return get_pointer(value[0]);
                }
            };
        }

        template<typename T, std::size_t N> void save(
              archive & ar
            , std::string const & path
            , std::array<T, N> const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            using alps::cast;
            if (ar.is_group(path))
                ar.delete_group(path);
            if (is_continuous<T>::value && value.size() == 0)
                ar.write(path, static_cast<typename scalar_type<std::array<T, N> >::type const *>(NULL), std::vector<std::size_t>());
            else if (is_continuous<T>::value) {
                std::vector<std::size_t> extent(get_extent(value));
                std::copy(extent.begin(), extent.end(), std::back_inserter(size));
                std::copy(extent.begin(), extent.end(), std::back_inserter(chunk));
                std::fill_n(std::back_inserter(offset), extent.size(), 0);
                ar.write(path, get_pointer(value), size, chunk, offset);
            } else if (value.size() == 0)
                ar.write(path, static_cast<int const *>(NULL), std::vector<std::size_t>());
            else if (is_vectorizable(value)) {
                size.push_back(value.size());
                chunk.push_back(1);
                offset.push_back(0);
                for(typename std::array<T, N>::const_iterator it = value.begin(); it != value.end(); ++it) {
                    offset.back() = it - value.begin();
                    save(ar, path, *it, size, chunk, offset);
                }
            } else {
                if (ar.is_data(path))
                    ar.delete_data(path);
                for(typename std::array<T, N>::const_iterator it = value.begin(); it != value.end(); ++it)
                    save(ar, ar.complete_path(path) + "/" + cast<std::string>(it - value.begin()), *it);
            }
        }

        template<typename T, std::size_t N> void load(
              archive & ar
            , std::string const & path
            , std::array<T, N> & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        ) {
            using alps::cast;
            if (ar.is_group(path)) {
                std::vector<std::string> children = ar.list_children(path);
                if (children.size() != N)
                    throw invalid_path("size does not match: " + path + ALPS_STACKTRACE);
                for (typename std::vector<std::string>::const_iterator it = children.begin(); it != children.end(); ++it)
                    load(ar, ar.complete_path(path) + "/" + *it, value[cast<std::size_t>(*it)]);
            } else {
                if (ar.is_complex(path) != has_complex_elements<T>::value)
                    throw archive_error("no complex value in archive" + ALPS_STACKTRACE);
                std::vector<std::size_t> size(ar.extent(path));
                if (size.size() > 0 && N != *(size.begin() + chunk.size()) && (is_continuous<T>::value || *(size.begin() + chunk.size()) > 0))
                    throw archive_error("dimensions do not match" + ALPS_STACKTRACE);
                if (is_continuous<T>::value) {
                    set_extent(value, std::vector<std::size_t>(size.begin() + chunk.size(), size.end()));
                    if (value.size()) {
                        std::copy(size.begin() + chunk.size(), size.end(), std::back_inserter(chunk));
                        std::fill_n(std::back_inserter(offset), size.size() - offset.size(), 0);
                        ar.read(path, get_pointer(value), chunk, offset);
                    }
                } else {
                    set_extent(value, std::vector<std::size_t>(1, *(size.begin() + chunk.size())));
                    chunk.push_back(1);
                    offset.push_back(0);
                    for(typename std::array<T, N>::iterator it = value.begin(); it != value.end(); ++it) {
                        offset.back() = it - value.begin();
                        load(ar, path, *it, chunk, offset);
                    }
                }
            }
        }
    }
}

#endif
