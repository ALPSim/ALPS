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

#ifndef ALPS_MATRIX_INTERFACE_HPP
#define ALPS_MATRIX_INTERFACE_HPP

#include <alps/numeric/matrix/matrix_concept_check.hpp>


// BOOST_CONCEPT_ASSERT((MATRIX_TEMPLATE)); 
// This macro creates free functions that call member functions with the same
// name, e.g. swap_cols(A,i,j) -> A.swap_cols(i,j)
#define ALPS_IMPLEMENT_MATRIX_INTERFACE(MATRIX_TEMPLATE, TEMPLATE_PARAMETERS) \
    template TEMPLATE_PARAMETERS \
    typename MATRIX_TEMPLATE::size_type num_rows(MATRIX_TEMPLATE const& m) { \
        return m.num_rows(); \
    } \
    template TEMPLATE_PARAMETERS \
    typename MATRIX_TEMPLATE::size_type num_cols(MATRIX_TEMPLATE const& m) { \
        return m.num_cols(); \
    } \
    template TEMPLATE_PARAMETERS \
    void swap_rows(MATRIX_TEMPLATE & m, typename MATRIX_TEMPLATE::size_type i1, typename MATRIX_TEMPLATE::size_type i2) { \
        m.swap_rows(i1,i2); \
    } \
    template TEMPLATE_PARAMETERS \
    void swap_cols(MATRIX_TEMPLATE & m, typename MATRIX_TEMPLATE::size_type i1, typename MATRIX_TEMPLATE::size_type i2) { \
        m.swap_cols(i1,i2); \
    } \
    template TEMPLATE_PARAMETERS \
    std::pair<typename MATRIX_TEMPLATE::row_element_iterator, typename MATRIX_TEMPLATE::row_element_iterator> \
    row(MATRIX_TEMPLATE & m, typename MATRIX_TEMPLATE::size_type i) { \
        return m.row(i); \
    } \
    template TEMPLATE_PARAMETERS \
    std::pair<typename MATRIX_TEMPLATE::const_row_element_iterator, typename MATRIX_TEMPLATE::const_row_element_iterator> \
    row(MATRIX_TEMPLATE const& m, typename MATRIX_TEMPLATE::size_type i) { \
        return m.row(i); \
    } \
    template TEMPLATE_PARAMETERS \
    std::pair<typename MATRIX_TEMPLATE::col_element_iterator, typename MATRIX_TEMPLATE::col_element_iterator> \
    col(MATRIX_TEMPLATE & m, typename MATRIX_TEMPLATE::size_type j) { \
        return m.col(j); \
    } \
    template TEMPLATE_PARAMETERS \
    std::pair<typename MATRIX_TEMPLATE::const_col_element_iterator, typename MATRIX_TEMPLATE::const_col_element_iterator> \
    col(MATRIX_TEMPLATE const& m, typename MATRIX_TEMPLATE::size_type j) { \
        return m.col(j); \
    }

#define ALPS_IMPLEMENT_MATRIX_DIAGONAL_ITERATOR_INTERFACE(MATRIX_TEMPLATE, TEMPLATE_PARAMETERS) \
    template TEMPLATE_PARAMETERS \
    std::pair<typename MATRIX_TEMPLATE::diagonal_iterator, typename MATRIX_TEMPLATE::diagonal_iterator> \
    diagonal(MATRIX_TEMPLATE & m) { \
        return m.diagonal(); \
    } \
    template TEMPLATE_PARAMETERS \
    std::pair<typename MATRIX_TEMPLATE::const_diagonal_iterator, typename MATRIX_TEMPLATE::const_diagonal_iterator> \
    diagonal(MATRIX_TEMPLATE const& m)  { \
        return m.diagonal(); \
    }

#define ALPS_IMPLEMENT_MATRIX_ELEMENT_ITERATOR_INTERFACE(MATRIX_TEMPLATE, TEMPLATE_PARAMETERS) \
    template TEMPLATE_PARAMETERS \
    std::pair<typename MATRIX_TEMPLATE::element_iterator, typename MATRIX_TEMPLATE::element_iterator> \
    elements(MATRIX_TEMPLATE & m) { \
        return m.elements(); \
    } \
    template TEMPLATE_PARAMETERS \
    std::pair<typename MATRIX_TEMPLATE::const_element_iterator, typename MATRIX_TEMPLATE::const_element_iterator> \
    elements(MATRIX_TEMPLATE const& m)  { \
        return m.elements(); \
    }

#endif //ALPS_MATRIX_INTERFACE_HPP
