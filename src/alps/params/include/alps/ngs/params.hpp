/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_PARAMS_HPP
#define ALPS_NGS_PARAMS_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/config.hpp>
#include <alps/params_export.h>
#include <alps/ngs/detail/paramvalue.hpp>
#include <alps/ngs/detail/paramproxy.hpp>
#include <alps/ngs/detail/paramiterator.hpp>

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
#include <functional>

namespace alps {

    class params {

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
                , value_reader_(arg.value_reader_)
            {}

            ALPS_PARAMS_DECL params(hdf5::archive ar, std::string const & path = "/parameters");

            ALPS_PARAMS_DECL std::size_t size() const;

            ALPS_PARAMS_DECL void erase(std::string const &);

            ALPS_PARAMS_DECL value_type operator[](std::string const &);

            ALPS_PARAMS_DECL value_type const operator[](std::string const &) const;

            ALPS_PARAMS_DECL bool defined(std::string const &) const;

            // Direct native lookup for consumers that need to inspect the
            // stored variant. The returned pointer remains owned by params
            // and is null when the key is absent.
            ALPS_PARAMS_DECL detail::paramvalue const * find(std::string const &) const;

            ALPS_PARAMS_DECL iterator begin();
            ALPS_PARAMS_DECL const_iterator begin() const;

            ALPS_PARAMS_DECL iterator end();
            ALPS_PARAMS_DECL const_iterator end() const;

            ALPS_PARAMS_DECL void save(hdf5::archive &) const;

            ALPS_PARAMS_DECL void load(hdf5::archive &);

            // A binding-owned decoder, preserved when parameters are copied
            // into a native simulation. Native-only parameters need none.
            typedef std::function<detail::paramvalue(hdf5::archive &)> value_reader;
            void set_value_reader(value_reader reader) { value_reader_ = std::move(reader); }

            #ifdef ALPS_HAVE_MPI
                ALPS_PARAMS_DECL void broadcast(boost::mpi::communicator const &, int = 0);
            #endif

        private:

            friend class boost::serialization::access;
            
            template<class Archive> void serialize(Archive & ar, const unsigned int) {
                ar & keys
                   & values
                ;
            }

            ALPS_PARAMS_DECL void setter(std::string const &, detail::paramvalue const &);

            ALPS_PARAMS_DECL detail::paramvalue getter(std::string const &);

            std::vector<std::string> keys;
            std::map<std::string, detail::paramvalue> values;
            value_reader value_reader_;
    };

    ALPS_PARAMS_DECL std::ostream & operator<<(std::ostream & os, params const & arg);
}

#endif
