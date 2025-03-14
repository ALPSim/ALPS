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
#ifndef ALPS_NUMERIC_IS_BLAS_DISPATCHABLE_HPP
#define ALPS_NUMERIC_IS_BLAS_DISPATCHABLE_HPP

#include <boost/type_traits/remove_const.hpp>
#include <boost/mpl/bool_fwd.hpp>
#include <boost/mpl/and.hpp>
#include <alps/numeric/matrix/detail/blasmacros.hpp>

namespace alps {
namespace numeric {

template <typename T>
struct supports_blas : boost::mpl::false_
{
};

#define ALPS_FUNDAMENTAL_TYPES_BLAS_TRAITS(T) \
template <> \
struct supports_blas<T> : boost::mpl::true_ {};

ALPS_IMPLEMENT_FOR_ALL_BLAS_TYPES(ALPS_FUNDAMENTAL_TYPES_BLAS_TRAITS)

#undef ALPS_FUNDAMENTAL_TYPES_BLAS_TRAITS

template <typename T1, typename T2>
struct is_blas_dispatchable;

template <typename T1, typename T2, typename Tag1, typename Tag2>
struct is_blas_dispatchable_helper
:   boost::mpl::and_<
          typename boost::mpl::and_<
                supports_blas<typename boost::remove_const<T1>::type>
              , supports_blas<typename boost::remove_const<T2>::type>
          >::type
        , is_blas_dispatchable<typename T1::value_type, typename T2::value_type>
    >::type
{
};

template <typename T1, typename T2, typename Tag1>
struct is_blas_dispatchable_helper<T1,T2,Tag1,tag::scalar>
:   boost::mpl::and_<
          typename boost::mpl::and_<
                supports_blas<typename boost::remove_const<T1>::type>
              , supports_blas<typename boost::remove_const<T2>::type>
          >::type
        , is_blas_dispatchable<typename T1::value_type, T2>
    >::type
{
};

template <typename T1, typename T2, typename Tag2>
struct is_blas_dispatchable_helper<T1,T2,tag::scalar,Tag2>
:   boost::mpl::and_<
          typename boost::mpl::and_<
                supports_blas<typename boost::remove_const<T1>::type>
              , supports_blas<typename boost::remove_const<T2>::type>
          >::type
        , is_blas_dispatchable<T1, typename T2::value_type>
    >::type
{
};

template <typename T1, typename T2>
struct is_blas_dispatchable_helper<T1,T2,tag::scalar,tag::scalar>
: boost::mpl::false_
{
};

template <typename T>
struct is_blas_dispatchable_helper<T,T,tag::scalar,tag::scalar>
: supports_blas<typename boost::remove_const<T>::type>::type
{
};

template <typename T1,typename T2>
struct is_blas_dispatchable
: is_blas_dispatchable_helper<T1,T2, typename get_entity<T1>::type, typename get_entity<T2>::type>::type
{
};


} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_IS_BLAS_DISPATCHABLE_HPP
