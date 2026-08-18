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

#include <string>
#include <cassert>

namespace alps {
    namespace detail {

        template<typename params_type, typename value_type> class paramiterator
            : public boost::forward_iterator_helper<
                  paramiterator<params_type, value_type>
                , value_type
                , std::ptrdiff_t
                , value_type *
                , value_type &
            >
        {
            public:

                paramiterator(paramiterator const & arg)
                    : params(arg.params)
                    , it(arg.it)
                {}

                paramiterator(
                      params_type & p
                    , std::vector<std::string>::const_iterator i
                )
                    : params(p)
                    , it(i)
                {}

                operator paramiterator<const params_type, const value_type>() const {
                    return paramiterator<const params_type, const value_type>(params, it);
                }

                value_type & operator*() const {
                    assert(params.values.find(*it) != params.values.end());
                    return *params.values.find(*it);
                }

                void operator++() {
                    ++it;
                }

                bool operator==(paramiterator<params_type, value_type> const & arg) const {
                    return it == arg.it;
                }

            private:

                params_type & params;
                std::vector<std::string>::const_iterator it;
        };

    }
}
