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
#include "matrix_unit_tests.hpp"

using alps::numeric::matrix;
using alps::numeric::transpose_view;

BOOST_AUTO_TEST_CASE_TEMPLATE( transpose_test, T, test_types )
{
    matrix<T> a(30,20);
    fill_matrix_with_numbers(a);

    matrix<T> b;

    b = transpose(a);

    matrix<T> const c(a);

    matrix<T> d = transpose(c);

    for(unsigned int j=0; j < num_cols(a); ++j)
        for(unsigned int i=0; i < num_rows(a); ++i)
        {
            BOOST_CHECK_EQUAL(a(i,j),b(j,i));
            BOOST_CHECK_EQUAL(c(i,j),d(j,i));
        }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( transpose_plus_test, T, test_types )
{
    matrix<T> a(30,20);
    matrix<T> b(20,30);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);
    matrix<T> ao(a);
    matrix<T> bo(b);

    matrix<T> c = transpose(a) + b;
    matrix<T> d = a + transpose(b);

    for(unsigned int j=0; j < num_cols(c); ++j)
        for(unsigned int i=0; i < num_rows(c); ++i)
            BOOST_CHECK_EQUAL(c(i,j), a(j,i)+b(i,j));

    for(unsigned int j=0; j < num_cols(d); ++j)
        for(unsigned int i=0; i < num_rows(d); ++i)
            BOOST_CHECK_EQUAL(d(i,j), a(i,j)+b(j,i));

    BOOST_CHECK_EQUAL(a,ao);
    BOOST_CHECK_EQUAL(b,bo);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( transpose_plus_transpose_test, T, test_types )
{
    matrix<T> a(17,24);
    matrix<T> b(17,24);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);
    matrix<T> ao(a);
    matrix<T> bo(b);

    matrix<T> c = transpose(a) + transpose(b);

    for(unsigned int j=0; j < num_cols(c); ++j)
        for(unsigned int i=0; i < num_rows(c); ++i)
            BOOST_CHECK_EQUAL(c(i,j), a(j,i)+b(j,i));

    BOOST_CHECK_EQUAL(a,ao);
    BOOST_CHECK_EQUAL(b,bo);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( transpose_minus_test, T, test_types )
{
    matrix<T> a(30,20);
    matrix<T> b(20,30);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b,T(4));
    matrix<T> ao(a);
    matrix<T> bo(b);

    matrix<T> c = transpose(a) - b;
    matrix<T> d = a - transpose(b);

    for(unsigned int j=0; j < num_cols(c); ++j)
        for(unsigned int i=0; i < num_rows(c); ++i)
            BOOST_CHECK_EQUAL(c(i,j), a(j,i)-b(i,j));

    for(unsigned int j=0; j < num_cols(d); ++j)
        for(unsigned int i=0; i < num_rows(d); ++i)
            BOOST_CHECK_EQUAL(d(i,j), a(i,j)-b(j,i));

    BOOST_CHECK_EQUAL(a,ao);
    BOOST_CHECK_EQUAL(b,bo);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( transpose_minus_transpose_test, T, test_types )
{
    matrix<T> a(17,24);
    matrix<T> b(17,24);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b,T(5));
    matrix<T> ao(a);
    matrix<T> bo(b);

    matrix<T> c = transpose(a) - transpose(b);

    for(unsigned int j=0; j < num_cols(c); ++j)
        for(unsigned int i=0; i < num_rows(c); ++i)
            BOOST_CHECK_EQUAL(c(i,j), a(j,i)-b(j,i));

    BOOST_CHECK_EQUAL(a,ao);
    BOOST_CHECK_EQUAL(b,bo);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( transpose_matrix_matrix_multiply_test, T, test_types )
{
    matrix<T> a(30,20);
    matrix<T> b(30,50);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);

    matrix<T> c = transpose(a) * b;

    BOOST_CHECK_EQUAL(num_rows(c), num_cols(a));
    BOOST_CHECK_EQUAL(num_cols(c), num_cols(b));

    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_cols(c); ++j)
        {
            T result(0);
            for(unsigned int k=0; k< num_rows(a); ++k)
                result += a(k,i) * b(k,j);
            BOOST_CHECK_EQUAL(c(i,j),result);
        }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( matrix_transpose_matrix_multiply_test, T, test_types )
{
    matrix<T> a(20,30);
    matrix<T> b(50,30);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);

    matrix<T> c = a * transpose(b);

    BOOST_CHECK_EQUAL(num_rows(c), num_rows(a));
    BOOST_CHECK_EQUAL(num_cols(c), num_rows(b));

    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_cols(c); ++j)
        {
            T result(0);
            for(unsigned int k=0; k< num_cols(a); ++k)
                result += a(i,k) * b(j,k);
            BOOST_CHECK_EQUAL(c(i,j),result);
        }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( transpose_conugate_matrix_matrix_multiply_test, T, test_types )
{
    using alps::numeric::conj;
    matrix<T> a(30,20);
    matrix<T> b(30,50);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);

    matrix<T> c = transpose(conj(a)) * b;

    BOOST_CHECK_EQUAL(num_rows(c), num_cols(a));
    BOOST_CHECK_EQUAL(num_cols(c), num_cols(b));

    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_cols(c); ++j)
        {
            T result(0);
            for(unsigned int k=0; k< num_rows(a); ++k)
                result += conj(a(k,i)) * b(k,j);
            BOOST_CHECK_EQUAL(c(i,j),result);
        }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( transpose_transpose_multiply_test, T, test_types)
{
    matrix<T> a(50,60);
    matrix<T> b(40,50);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);

    matrix<T> c = transpose(a) * transpose(b);

    BOOST_CHECK_EQUAL(num_rows(c), num_cols(a));
    BOOST_CHECK_EQUAL(num_cols(c), num_rows(b));

    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_cols(c); ++j)
        {
            T result(0);
            for(unsigned int k=0; k< num_rows(a); ++k)
                result += a(k,i) * b(j,k);
            BOOST_CHECK_EQUAL(c(i,j),result);
        }
}
