/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "matrix_unit_tests.hpp"

using alps::numeric::matrix;
using alps::numeric::transpose_view;

template<class T> class MatrixTransposeView : public ::testing::Test {};
TYPED_TEST_SUITE(MatrixTransposeView, test_types);

TYPED_TEST(MatrixTransposeView, transpose_test)
{
    using T = TypeParam;
    matrix<T> a(30,20);
    fill_matrix_with_numbers(a);

    matrix<T> b;

    b = transpose(a);

    matrix<T> const c(a);

    matrix<T> d = transpose(c);

    for(unsigned int j=0; j < num_cols(a); ++j)
        for(unsigned int i=0; i < num_rows(a); ++i)
        {
            EXPECT_EQ(a(i,j),b(j,i));
            EXPECT_EQ(c(i,j),d(j,i));
        }
}

TYPED_TEST(MatrixTransposeView, transpose_plus_test)
{
    using T = TypeParam;
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
            EXPECT_EQ(c(i,j), a(j,i)+b(i,j));

    for(unsigned int j=0; j < num_cols(d); ++j)
        for(unsigned int i=0; i < num_rows(d); ++i)
            EXPECT_EQ(d(i,j), a(i,j)+b(j,i));

    EXPECT_EQ(a,ao);
    EXPECT_EQ(b,bo);
}

TYPED_TEST(MatrixTransposeView, transpose_plus_transpose_test)
{
    using T = TypeParam;
    matrix<T> a(17,24);
    matrix<T> b(17,24);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);
    matrix<T> ao(a);
    matrix<T> bo(b);

    matrix<T> c = transpose(a) + transpose(b);

    for(unsigned int j=0; j < num_cols(c); ++j)
        for(unsigned int i=0; i < num_rows(c); ++i)
            EXPECT_EQ(c(i,j), a(j,i)+b(j,i));

    EXPECT_EQ(a,ao);
    EXPECT_EQ(b,bo);
}

TYPED_TEST(MatrixTransposeView, transpose_minus_test)
{
    using T = TypeParam;
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
            EXPECT_EQ(c(i,j), a(j,i)-b(i,j));

    for(unsigned int j=0; j < num_cols(d); ++j)
        for(unsigned int i=0; i < num_rows(d); ++i)
            EXPECT_EQ(d(i,j), a(i,j)-b(j,i));

    EXPECT_EQ(a,ao);
    EXPECT_EQ(b,bo);
}

TYPED_TEST(MatrixTransposeView, transpose_minus_transpose_test)
{
    using T = TypeParam;
    matrix<T> a(17,24);
    matrix<T> b(17,24);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b,T(5));
    matrix<T> ao(a);
    matrix<T> bo(b);

    matrix<T> c = transpose(a) - transpose(b);

    for(unsigned int j=0; j < num_cols(c); ++j)
        for(unsigned int i=0; i < num_rows(c); ++i)
            EXPECT_EQ(c(i,j), a(j,i)-b(j,i));

    EXPECT_EQ(a,ao);
    EXPECT_EQ(b,bo);
}

TYPED_TEST(MatrixTransposeView, transpose_matrix_matrix_multiply_test)
{
    using T = TypeParam;
    matrix<T> a(30,20);
    matrix<T> b(30,50);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);

    matrix<T> c = transpose(a) * b;

    EXPECT_EQ(num_rows(c), num_cols(a));
    EXPECT_EQ(num_cols(c), num_cols(b));

    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_cols(c); ++j)
        {
            T result(0);
            for(unsigned int k=0; k< num_rows(a); ++k)
                result += a(k,i) * b(k,j);
            EXPECT_EQ(c(i,j),result);
        }
}

TYPED_TEST(MatrixTransposeView, matrix_transpose_matrix_multiply_test)
{
    using T = TypeParam;
    matrix<T> a(20,30);
    matrix<T> b(50,30);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);

    matrix<T> c = a * transpose(b);

    EXPECT_EQ(num_rows(c), num_rows(a));
    EXPECT_EQ(num_cols(c), num_rows(b));

    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_cols(c); ++j)
        {
            T result(0);
            for(unsigned int k=0; k< num_cols(a); ++k)
                result += a(i,k) * b(j,k);
            EXPECT_EQ(c(i,j),result);
        }
}

TYPED_TEST(MatrixTransposeView, transpose_conugate_matrix_matrix_multiply_test)
{
    using T = TypeParam;
    using alps::numeric::conj;
    matrix<T> a(30,20);
    matrix<T> b(30,50);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);

    matrix<T> c = transpose(conj(a)) * b;

    EXPECT_EQ(num_rows(c), num_cols(a));
    EXPECT_EQ(num_cols(c), num_cols(b));

    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_cols(c); ++j)
        {
            T result(0);
            for(unsigned int k=0; k< num_rows(a); ++k)
                result += conj(a(k,i)) * b(k,j);
            EXPECT_EQ(c(i,j),result);
        }
}

TYPED_TEST(MatrixTransposeView, transpose_transpose_multiply_test)
{
    using T = TypeParam;
    matrix<T> a(50,60);
    matrix<T> b(40,50);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);

    matrix<T> c = transpose(a) * transpose(b);

    EXPECT_EQ(num_rows(c), num_cols(a));
    EXPECT_EQ(num_cols(c), num_rows(b));

    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_cols(c); ++j)
        {
            T result(0);
            for(unsigned int k=0; k< num_rows(a); ++k)
                result += a(k,i) * b(j,k);
            EXPECT_EQ(c(i,j),result);
        }
}
