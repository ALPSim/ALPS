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

#include <alps/ngs/params.hpp>

#include <alps/parameter.h>

#include <boost/bind.hpp>

#include <algorithm>

namespace alps {

    params::params(hdf5::archive ar, std::string const & path) {
        std::string context = ar.get_context();
        ar.set_context(path);
        load(ar);
        ar.set_context(context);
    }

    params::params(boost::filesystem::path const & path) {
        boost::filesystem::ifstream ifs(path);
        Parameters par(ifs);
        for (Parameters::const_iterator it = par.begin(); it != par.end(); ++it) {
            detail::paramvalue val(it->value());
            setter(it->key(), val);
        }
    }

    #ifdef ALPS_HAVE_PYTHON
        params::params(boost::python::dict const & arg) {
            boost::python::extract<boost::python::dict> dict(arg);
            if (!dict.check())
                throw std::invalid_argument("parameters can only be created from a dict" + ALPS_STACKTRACE);
            const boost::python::list keys = dict().keys();
            for (std::size_t i = 0; i < boost::python::len(keys); ++i) {
                boost::python::object pyk = keys[i];
                std::string k = boost::python::call_method<std::string>(pyk.ptr(), "__str__");
                setter(k, dict().get(pyk));
            }
        }

        // TODO: merge with params::params(boost::filesystem::path const & path);
        params::params(boost::python::str const & arg) {
            std::string path = boost::python::extract<std::string>(arg)();
            boost::filesystem::ifstream ifs(path);
            Parameters par(ifs);
            for (Parameters::const_iterator it = par.begin(); it != par.end(); ++it) {
                detail::paramvalue val(it->value());
                setter(it->key(), val);
            }
        }

    #endif

    std::size_t params::size() const {
        return keys.size();
    }

    void params::erase(std::string const & key) {
        if (!defined(key))
            throw std::invalid_argument("the key " + key + " does not exists" + ALPS_STACKTRACE);
        keys.erase(find(keys.begin(), keys.end(), key));
        values.erase(key);
    }

    params::value_type params::operator[](std::string const & key) {
        return value_type(
            defined(key),
            boost::bind(&params::getter, boost::ref(*this), key),
            boost::bind(&params::setter, boost::ref(*this), key, _1),
            key
        );
    }

    params::value_type const params::operator[](std::string const & key) const {
        return defined(key)
            ? value_type(values.find(key)->second, key)
            : value_type(key)
        ;
    }

    bool params::defined(std::string const & key) const {
        return values.find(key) != values.end();
    }

    params::iterator params::begin() {
        return iterator(*this, keys.begin());
    }

    params::const_iterator params::begin() const {
        return const_iterator(*this, keys.begin());
    }

    params::iterator params::end() {
        return iterator(*this, keys.end());
    }

    params::const_iterator params::end() const {
        return const_iterator(*this, keys.end());
    }

    void params::save(hdf5::archive & ar) const {
        for (params::const_iterator it = begin(); it != end(); ++it)
            ar[it->first] << it->second;
    }

    void params::load(hdf5::archive & ar) {
        keys.clear();
        values.clear();
        std::vector<std::string> list = ar.list_children(ar.get_context());
        for (std::vector<std::string>::const_iterator it = list.begin(); it != list.end(); ++it) {
            detail::paramvalue value;
            ar[*it] >> value;
            setter(*it, value);
        }
    }

    #ifdef ALPS_HAVE_MPI
        void params::broadcast(boost::mpi::communicator const & comm, int root) {
            boost::mpi::broadcast(comm, *this, root);
        }
    #endif

    void params::setter(std::string const & key, detail::paramvalue const & value) {
        if (!defined(key))
            keys.push_back(key);
        values[key] = value;
    }

    detail::paramvalue params::getter(std::string const & key) {
        return values[key];
    }


    std::ostream & operator<<(std::ostream & os, params const & v) {
        for (params::const_iterator it = v.begin(); it != v.end(); ++it)
            os << it->first << " = " << it->second << std::endl;
        return os;
    }
}
