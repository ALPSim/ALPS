/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2013 by Mario Koenz <mkoenz@ethz.ch>                       *
 *                              Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_NGS_ACCUMULATOR_WEIGHT_HPP
#define ALPS_NGS_ACCUMULATOR_WEIGHT_HPP

#include <alps/ngs/accumulator/feature.hpp>
#include <alps/ngs/accumulator/parameter.hpp>
#include <alps/ngs/accumulator/feature/count.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/stacktrace.hpp>
#include <alps/ngs/short_print.hpp>

#include <boost/utility.hpp>

#include <stdexcept>

namespace alps {
    namespace accumulator {

        template<typename T> class base_wrapper;

        // this should be called namespace tag { struct weight; }
        // but gcc <= 4.4 has lookup error, so name it different
        struct weight_tag;

        template<typename T> struct has_feature<T, weight_tag> {
            template<typename R, typename C> static char helper(R(C::*)() const);
            template<typename C> static char check(boost::integral_constant<std::size_t, sizeof(helper(&C::owns_weight))>*);
            template<typename C> static double check(...);
            typedef boost::integral_constant<bool, sizeof(char) == sizeof(check<T>(0))> type;
        };

        namespace detail {
            struct no_weight_type {};
            template <bool, typename T> struct weight_type_impl {
                typedef no_weight_type type;
            };
            template <typename T> struct weight_type_impl<true, T> {
                typedef typename T::weight_type type;
            };
        }

        template<typename T> struct weight_type {
            typedef typename detail::weight_type_impl<has_feature<T, weight_tag>::type::value, T>::type type;
        };

        template<typename T> base_wrapper<typename value_type<T>::type> const * weight(T const & arg) {
            return arg.weight();
        }

        namespace detail {

            template<typename A> typename boost::enable_if<
                  typename has_feature<A, weight_tag>::type
                , base_wrapper<typename value_type<A>::type> const *
            >::type weight_impl(A const & acc) {
                return weight(acc);
            }

            template<typename A> typename boost::disable_if<
                  typename has_feature<A, weight_tag>::type
                , base_wrapper<typename value_type<A>::type> const *
            >::type weight_impl(A const & acc) {
                throw std::runtime_error(std::string(typeid(A).name()) + " has no weight-method" + ALPS_STACKTRACE);
                return NULL;
            }
        }

        namespace impl {

            template<typename T, typename B> class BaseWrapper<T, weight_tag, B> : public B {
                public:
                    virtual bool has_weight() const = 0;
                    virtual base_wrapper<T> const * weight() const = 0;
            };

            template<typename T, typename B> class DerivedWrapper<T, weight_tag, B> : public B {
                public:
                    DerivedWrapper(): B() {}
                    DerivedWrapper(T const & arg): B(arg) {}

                    bool has_weight() const { return has_feature<T, weight_tag>::type::value; }
                    base_wrapper<typename value_type<T>::type> const * weight() const { return detail::weight_impl(this->m_data); } 
            };

        }
    }
}

 #endif
