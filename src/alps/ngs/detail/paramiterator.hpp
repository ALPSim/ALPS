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
