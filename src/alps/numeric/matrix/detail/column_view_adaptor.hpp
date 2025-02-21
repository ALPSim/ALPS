/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Andreas Hehn <hehn@phys.ethz.ch>                   *
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
#ifndef ALPS_NUMERIC_MATRIX_COLUMN_VIEW_ADAPTOR_HPP
#define ALPS_NUMERIC_MATRIX_COLUMN_VIEW_ADAPTOR_HPP

#include <boost/numeric/bindings/detail/adaptor.hpp>

namespace alps { namespace numeric {
    template <typename Matrix>
    class column_view;
} }

//
// An adaptor for the column_view of alps::numeric::matrix to the boost::numeric::bindings
//

namespace boost { namespace numeric { namespace bindings { namespace detail {

    template< typename T, typename MemoryBlock, typename Id, typename Enable >
    struct adaptor< ::alps::numeric::column_view< ::alps::numeric::matrix<T,MemoryBlock> >, Id, Enable>
    {
        typedef typename copy_const< Id, T >::type value_type;
        typedef std::ptrdiff_t  size_type;

        typedef mpl::map<
            mpl::pair< tag::value_type,     value_type >,
            mpl::pair< tag::entity,         tag::vector >,
            mpl::pair< tag::size_type<1>,   size_type >,
            mpl::pair< tag::data_structure, tag::linear_array >,
            mpl::pair< tag::stride_type<1>, tag::contiguous >
        > property_map;

        static std::ptrdiff_t size1( const Id& id ) {
            return id.size();
        }

        static value_type* begin_value( Id& id ) {
            return &(*id.begin());
        }

        static value_type* end_value( Id& id ) {
            return &(*id.end());
        }
    };
}}}}

#endif // ALPS_NUMERIC_MATRIX_COLUMN_VIEW_ADAPTOR_HPP
