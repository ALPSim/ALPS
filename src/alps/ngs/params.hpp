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

#ifndef ALPS_NGS_PARAMS_HPP
#define ALPS_NGS_PARAMS_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/config.hpp>
#include <alps/ngs/detail/paramvalue.hpp>
#include <alps/ngs/detail/paramproxy.hpp>
#include <alps/ngs/detail/paramiterator.hpp>

#ifdef ALPS_HAVE_PYTHON
    #include <alps/ngs/boost_python.hpp>
    #include <boost/python/dict.hpp>
#endif

#include <boost/filesystem.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp> 

#ifdef ALPS_HAVE_MPI
    #include <alps/ngs/boost_mpi.hpp>
#endif

#include <map>
#include <vector>
#include <string>

namespace alps {

    class ALPS_DECL params {

        typedef std::map<std::string, detail::paramvalue>::value_type iterator_value_type;

        friend class detail::paramiterator<params, iterator_value_type>;
        friend class detail::paramiterator<params const, iterator_value_type const>;

        public:

            typedef detail::paramiterator<params, iterator_value_type> iterator;
            typedef detail::paramiterator<params const, iterator_value_type const> const_iterator;
            typedef detail::paramproxy value_type;

            params() {}

            params(params const & arg)
                : keys(arg.keys)
                , values(arg.values)
            {}

            params(hdf5::archive ar, std::string const & path = "/parameters");

            params(boost::filesystem::path const &);

            #ifdef ALPS_HAVE_PYTHON
                params(boost::python::dict const & arg);
                params(boost::python::str const & arg);
            #endif

            std::size_t size() const;

            void erase(std::string const &);

            value_type operator[](std::string const &);

            value_type const operator[](std::string const &) const;

            bool defined(std::string const &) const;

            iterator begin();
            const_iterator begin() const;

            iterator end();
            const_iterator end() const;

            void save(hdf5::archive &) const;

            void load(hdf5::archive &);

            #ifdef ALPS_HAVE_MPI
                void broadcast(boost::mpi::communicator const &, int = 0);
            #endif

        private:

            friend class boost::serialization::access;
            
            template<class Archive> void serialize(Archive & ar, const unsigned int) {
                ar & keys
                   & values
                ;
            }

            void setter(std::string const &, detail::paramvalue const &);

            detail::paramvalue getter(std::string const &);

            std::vector<std::string> keys;
            std::map<std::string, detail::paramvalue> values;
    };

    ALPS_DECL std::ostream & operator<<(std::ostream & os, params const & arg);
}

#endif
