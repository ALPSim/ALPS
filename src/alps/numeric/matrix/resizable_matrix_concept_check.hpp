/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Andreas Hehn <hehn@phys.ethz.ch>                   *
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

#ifndef ALPS_RESIZABLE_MATRIX_CONCEPT_CHECK_HPP
#define ALPS_RESIZABLE_MATRIX_CONCEPT_CHECK_HPP
#include <boost/concept_check.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <stdexcept>
#include <alps/numeric/matrix/matrix_concept_check.hpp>

namespace alps {
    namespace numeric {

    template <typename X>
    struct ResizableMatrix
            : Matrix<X>
    {
        public:
        BOOST_CONCEPT_USAGE(ResizableMatrix)
        {
            typename boost::remove_const<X>::type x(1,1);

            // Resize
            resize(x,2,2);

            // Append
            std::vector<typename X::value_type> dataA(2,typename X::value_type());
            std::vector<typename X::value_type> dataB(4,typename X::value_type());
            append_rows(x, std::make_pair(dataA.begin(),dataA.end()) );
            append_rows(x, std::make_pair(dataA.begin(),dataA.end()),1);
            append_cols(x, std::make_pair(dataB.begin(),dataB.end()) );
            append_cols(x, std::make_pair(dataB.begin(),dataB.end()),1);

            // Remove
            remove_rows(x,1);
            remove_rows(x,1,1);
            remove_cols(x,1);
            remove_cols(x,1,1);

            // Insert
            insert_rows(x,1, std::make_pair(dataA.begin(),dataA.end()) );
            insert_rows(x,1, std::make_pair(dataA.begin(),dataA.end()),1);
            insert_cols(x,1, std::make_pair(dataB.begin(),dataB.end()) );
            insert_cols(x,1, std::make_pair(dataB.begin(),dataB.end()),1); 
        }
    };

    }  // end namespace numeric
} // end namespace alps
#endif //ALPS_RESIZABLE_MATRIX_CONCEPT_CHECK_HPP
