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
#ifndef ALPS_NUMERIC_OPERATORS_OP_ASSIGN_VECTOR_HPP
#define ALPS_NUMERIC_OPERATORS_OP_ASSIGN_VECTOR_HPP

#include <alps/numeric/matrix/detail/debug_output.hpp>
#include <alps/numeric/matrix/is_blas_dispatchable.hpp>
#include <boost/numeric/bindings/blas/level1/axpy.hpp>
#include <boost/numeric/bindings/blas/level1/scal.hpp>

namespace alps {
namespace numeric {
    namespace detail {
            template <typename T, typename T2>
            struct multiplies : public std::binary_function<T,T2,T>
            {
                inline T operator()(T t, T2 const& t2) const
                {
                    return t*t2;
                }
            };
    } // end namespace detail

    namespace impl {
#if defined(__clang_major__) && __clang_major__ < 3 || (__clang_major__ == 3 && __clang_minor__ == 0) || defined(BOOST_MSVC)
// Workaround for a compiler bug in clang 3.0 (and maybe earlier versions)
    template <typename Vector1, typename Vector2>
    void plus_assign_impl(Vector1& lhs, Vector2 const& rhs, tag::vector, tag::vector, boost::mpl::false_)
    {
        for(std::size_t i=0; i< lhs.size(); ++i)
            lhs[i] += rhs[i];
    }

    template <typename Vector1, typename Vector2>
    void minus_assign_impl(Vector1& lhs, Vector2 const& rhs, tag::vector, tag::vector, boost::mpl::false_)
    {
        for(std::size_t i=0; i< lhs.size(); ++i)
            lhs[i] -= rhs[i];
    }

    template <typename Vector, typename T2>
    void multiplies_assign_impl(Vector& lhs, T2 lambda, tag::vector tag1, tag::scalar tag2, boost::mpl::false_)
    {
        for(std::size_t i=0; i< lhs.size(); ++i)
            lhs[i] *= lambda;
    }
#else // defined(__clang_major__) && __clang_major__ < 3 || (__clang_major__ == 3 && __clang_minor__ == 0)
    template <typename Vector1, typename Vector2>
    void plus_assign_impl(Vector1& lhs, Vector2 const& rhs, tag::vector, tag::vector, boost::mpl::false_)
    {
        std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), std::plus<typename Vector1::value_type>());
    }

    template <typename Vector1, typename Vector2>
    void minus_assign_impl(Vector1& lhs, Vector2 const& rhs, tag::vector, tag::vector, boost::mpl::false_)
    {
        std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), std::minus<typename Vector1::value_type>());
    }

    template <typename Vector, typename T2>
    void multiplies_assign_impl(Vector& lhs, T2 lambda, tag::vector tag1, tag::scalar tag2, boost::mpl::false_)
    {
        using detail::multiplies;
        std::transform(lhs.begin(), lhs.end(), lhs.begin(), std::bind2nd(multiplies<typename Vector::value_type, T2>(), lambda));
    }
#endif // defined(__clang_major__) && __clang_major__ < 3 || (__clang_major__ == 3 && __clang_minor__ == 0)

    // BLAS overloads
    template <typename Vector1, typename Vector2>
    void plus_assign_impl(Vector1& lhs, Vector2 const& rhs, tag::vector tag1, tag::vector tag2, boost::mpl::true_)
    {
        ALPS_NUMERIC_MATRIX_DEBUG_OUTPUT( "using blas axpy for " << typeid(lhs).name() << " " << typeid(rhs).name() );
        boost::numeric::bindings::blas::axpy(1., rhs, lhs);
    }

    template <typename Vector1, typename Vector2>
    void minus_assign_impl(Vector1& lhs, Vector2 const& rhs, tag::vector tag1, tag::vector tag2, boost::mpl::true_)
    {
        ALPS_NUMERIC_MATRIX_DEBUG_OUTPUT( "using blas axpy for " << typeid(lhs).name() << " " << typeid(rhs).name() );
        boost::numeric::bindings::blas::axpy(-1., rhs, lhs);
    }

    template <typename Vector, typename T2>
    void multiplies_assign_impl(Vector& lhs, T2 lambda, tag::vector tag1, tag::scalar tag2, boost::mpl::true_)
    {
        ALPS_NUMERIC_MATRIX_DEBUG_OUTPUT( "using blas scal for " << typeid(lhs).name() << " " << typeid(lambda).name() );
        boost::numeric::bindings::blas::scal(lambda, lhs);
    }

    } // end namespace impl


    template <typename Vector1, typename Vector2>
    void plus_assign(Vector1& lhs, Vector2 const& rhs, tag::vector tag1, tag::vector tag2)
    {
        using impl::plus_assign_impl;
        plus_assign_impl(lhs,rhs,tag1,tag2,typename is_blas_dispatchable<Vector1,Vector2>::type());
    }

    template <typename Vector1, typename Vector2>
    void minus_assign(Vector1& lhs, Vector2 const& rhs, tag::vector tag1, tag::vector tag2)
    {
        using impl::minus_assign_impl;
        minus_assign_impl(lhs,rhs,tag1,tag2,typename is_blas_dispatchable<Vector1,Vector2>::type());
    }

    template <typename Vector, typename T2>
    void multiplies_assign(Vector& lhs, T2 lambda, tag::vector tag1, tag::scalar tag2)
    {
        using impl::multiplies_assign_impl;
        multiplies_assign_impl(lhs, lambda, tag1, tag2, typename is_blas_dispatchable<Vector,T2>::type());
    }

} // end namespace numeric
} // end namespace alps
#endif // ALPS_NUMERIC_OPERATORS_OP_ASSIGN_VECTOR_HPP
