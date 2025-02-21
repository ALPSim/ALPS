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

#ifndef ALPS_NGS_HDF5_PYTHON_CPP
#define ALPS_NGS_HDF5_PYTHON_CPP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/cast.hpp>
#include <alps/hdf5/pair.hpp>
#include <alps/ngs/stacktrace.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>

#include <alps/ngs/boost_python.hpp>
#include <alps/python/numpy_array.hpp>

#include <boost/scoped_ptr.hpp>

#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/object.hpp>

#include <string>
#include <iterator>
#include <stdexcept>

namespace alps {
    namespace hdf5 {

        namespace detail {

            template<> struct is_vectorizable<boost::python::object> {
                static bool apply(boost::python::object const & value);
            };

            template<> struct get_extent<boost::python::object> {
                static std::vector<std::size_t> apply(boost::python::object const & value);
            };

            template<> struct set_extent<boost::python::object> {
                static void apply(boost::python::object & value, std::vector<std::size_t> const & extent);
            };
        }

        ALPS_DECL void save(
              archive & ar
            , std::string const & path
            , boost::python::object const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        );

        ALPS_DECL void load(
              archive & ar
            , std::string const & path
            , boost::python::object & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        );
        
        namespace detail {

            template<> struct is_vectorizable<boost::python::list> {
                static bool apply(boost::python::list const & value);
            };

            template<> struct get_extent<boost::python::list> {
                static std::vector<std::size_t> apply(boost::python::list const & value);
            };

            template<> struct set_extent<boost::python::list> {
                static void apply(boost::python::list & value, std::vector<std::size_t> const & extent);
            };
        }

        ALPS_DECL void save(
              archive & ar
            , std::string const & path
            , boost::python::list const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        );

        ALPS_DECL void load(
              archive & ar
            , std::string const & path
            , boost::python::list & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        );

        namespace detail {

            template<> struct is_vectorizable<boost::python::tuple> {
                static bool apply(boost::python::tuple const & value);
            };

            template<> struct get_extent<boost::python::tuple> {
                static std::vector<std::size_t> apply(boost::python::tuple const & value);
            };
        }

        ALPS_DECL void save(
              archive & ar
            , std::string const & path
            , boost::python::tuple const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        );

        namespace detail {

            template<> struct is_vectorizable<alps::python::numpy::array> {
                static bool apply(alps::python::numpy::array const & value);
            };

            template<>  struct get_extent<alps::python::numpy::array> {
                static std::vector<std::size_t> apply(alps::python::numpy::array const & value);
            };

            template<>  struct set_extent<alps::python::numpy::array> {
                // To set the extent of a numpy array, we need the type, extent is set in load
                static void apply(alps::python::numpy::array & value, std::vector<std::size_t> const & extent);
            };
        }

        ALPS_DECL void save(
              archive & ar
            , std::string const & path
            , alps::python::numpy::array const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        );

        ALPS_DECL void load(
              archive & ar
            , std::string const & path
            , alps::python::numpy::array & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        );
        
        ALPS_DECL void save(
              archive & ar
            , std::string const & path
            , boost::python::dict const & value
            , std::vector<std::size_t> size = std::vector<std::size_t>()
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        );
        
        ALPS_DECL void load(
              archive & ar
            , std::string const & path
            , boost::python::dict & value
            , std::vector<std::size_t> chunk = std::vector<std::size_t>()
            , std::vector<std::size_t> offset = std::vector<std::size_t>()
        );
    }
}

#endif
