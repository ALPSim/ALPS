/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_DETAIL_PARAMPROXY_HPP
#define ALPS_NGS_DETAIL_PARAMPROXY_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/config.hpp>
#include <alps/ngs/detail/paramvalue.hpp>

#include <boost/function.hpp>
#include <boost/optional/optional.hpp>

#include <string>
#include <iostream>

namespace alps {
    
    namespace detail {

        class ALPS_DECL paramproxy {

            public:

                paramproxy(std::string const & k)
                    : defined(false)
                    , key(k)
                {}

                paramproxy(paramvalue const & v, std::string const & k)
                    : defined(true)
                    , key(k)
                    , value(v)
                {}

                paramproxy(
                      bool d
                    , boost::function<paramvalue()> const & g
                    , boost::function<void(paramvalue)> const & s
                    , std::string const & k
                )
                    : defined(d)
                    , key(k)
                    , getter(g)
                    , setter(s)
                {}

                paramproxy(paramproxy const & arg)
                    : defined(arg.defined)
                    , key(arg.key)
                    , value(arg.value)
                    , getter(arg.getter)
                    , setter(arg.setter)
                {}

                template<typename T> T cast() const {
                    if (!defined)
                        throw std::runtime_error("No parameter '" + key + "' available" + ALPS_STACKTRACE);
                    return (!value ? getter() : *value).cast<T>();
                }

                template<typename T> operator T () const {
                    return cast<T>();
                }

                template<typename T> paramproxy & operator=(T const & arg) {
                    if (!!value)
                        throw std::runtime_error("No reference to parameter '" + key + "' available" + ALPS_STACKTRACE);
                    setter(detail::paramvalue(arg));
                    return *this;
                }

                template<typename T> T or_default(T const & value) const {
                    return defined ? cast<T>() : value;
                }

                template<typename T> T operator|( T const & value) const {
                    return or_default(value);
                }

                std::string operator|(char const * value) const {
                    return or_default(std::string(value));
                }

                paramproxy const & operator|(paramproxy const & value) const {
                    return defined ? *this : value;
                }

                void save(hdf5::archive & ar) const;
                void load(hdf5::archive &);

				void print(std::ostream &) const;

            private:

                bool defined;
                std::string key;
                boost::optional<paramvalue> value;
                boost::function<paramvalue()> getter;
                boost::function<void(paramvalue)> setter;
        };

        ALPS_DECL std::ostream & operator<<(std::ostream & os, paramproxy const &);

        #define ALPS_NGS_PARAMPROXY_ADD_OPERATOR_DECL(T)                                 \
            ALPS_DECL T operator+(paramproxy const & p, T s);                            \
            ALPS_DECL T operator+(T s, paramproxy const & p);
        // vector<bool> is a native stored parameter type, but unlike the
        // historic numeric/string alternatives it has no meaningful or
        // portable element-wise operator+=.
        ALPS_NGS_FOREACH_PARAMETERVALUE_ADDABLE_TYPE(ALPS_NGS_PARAMPROXY_ADD_OPERATOR_DECL)
        #undef ALPS_NGS_PARAMPROXY_ADD_OPERATOR_DECL

        ALPS_DECL std::string operator+(paramproxy const & p, char const * s);
        ALPS_DECL std::string operator+(char const * s, paramproxy const & p);

    }
}
    
#endif
