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
#include <alps/numeric/matrix/column_view.hpp>

using alps::numeric::matrix;
using alps::numeric::vector;
using alps::numeric::column_view;

template<class T> class MatrixColumnView : public ::testing::Test {};
TYPED_TEST_SUITE(MatrixColumnView, test_types);

TYPED_TEST(MatrixColumnView, size_test)
{
    using T = TypeParam;
    matrix<T> c(15,5);
    column_view<matrix<T> > cv(c,3);
    EXPECT_EQ(cv.size(),15u);
}

TYPED_TEST(MatrixColumnView, const_element_access)
{
    using T = TypeParam;
    using std::distance;
    matrix<T> c(15,20);
    fill_matrix_with_numbers(c);
    column_view<matrix<T> > const cv(c,3);
    for(std::size_t i = 0; i < num_rows(c); ++i)
    {
        EXPECT_EQ( c(i,3), cv(i));
        EXPECT_EQ( c(i,3), cv[i]);
    }
}

TYPED_TEST(MatrixColumnView, const_iterator_test)
{
    using T = TypeParam;
    using std::distance;
    matrix<T> c(15,20);
    fill_matrix_with_numbers(c);
    column_view<matrix<T> > const cv(c,3);
    for(typename column_view<matrix<T> >::const_iterator it = cv.begin(); it != cv.end(); ++it)
        EXPECT_EQ( c(distance(cv.begin(),it),3), *it);
}

TYPED_TEST(MatrixColumnView, iterator_test)
{
    using T = TypeParam;
    using std::distance;
    matrix<T> c(15,20);
    fill_matrix_with_numbers(c);
    matrix<T> d(c);
    column_view<matrix<T> > cv(c,5);
    for(typename column_view<matrix<T> >::iterator it = cv.begin(); it != cv.end(); ++it)
        *it = 100*distance(cv.begin(),it);
    for(std::size_t i = 0; i < num_rows(c); ++i)
        EXPECT_EQ(c(i,5),T(100*i));

    for(std::size_t j=0; j < num_cols(c); ++j)
    {
        if(j != 5)
            for(std::size_t i = 0; i < num_rows(c); ++i)
                EXPECT_EQ(d(i,j),c(i,j));
    }
}

TYPED_TEST(MatrixColumnView, element_assign_test)
{
    using T = TypeParam;
    using std::distance;
    matrix<T> c(15,20);
    fill_matrix_with_numbers(c);
    matrix<T> d(c);
    column_view<matrix<T> > cv(c,5);
    for(std::size_t i = 0; i < num_rows(c); ++i)
        cv(i) = 100*i;

    for(std::size_t i = 0; i < num_rows(c); ++i)
        EXPECT_EQ(c(i,5),T(100*i));

    for(std::size_t i = 0; i < num_rows(c); ++i)
        cv[i] = 1000*i;

    for(std::size_t i = 0; i < num_rows(c); ++i)
        EXPECT_EQ(c(i,5),T(1000*i));

    for(std::size_t j=0; j < num_cols(c); ++j)
    {
        if(j != 5)
            for(std::size_t i = 0; i < num_rows(c); ++i)
                EXPECT_EQ(d(i,j),c(i,j));
    }
}

TYPED_TEST(MatrixColumnView, conversion_to_vector)
{
    using T = TypeParam;
    matrix<T> c(15,20);
    fill_matrix_with_numbers(c);
    matrix<T> d(c);
    column_view<matrix<T> > cv(c,5);
    vector<T> v(cv);
    EXPECT_EQ(num_rows(c),v.size());
    for(std::size_t i = 0; i < num_rows(c); ++i)
        EXPECT_EQ(c(i,5),v(i));
    EXPECT_EQ(c,d);
}

TYPED_TEST(MatrixColumnView, plus_assign)
{
    using T = TypeParam;
    matrix<T> c(15,20);
    fill_matrix_with_numbers(c);
    matrix<T> d(c);
    column_view<matrix<T> > cv(c,5);
    vector<T> v(15);
    fill_range_with_numbers(v.begin(),v.end(),0);
    cv += v;
    for(std::size_t i = 0; i < num_rows(c); ++i)
        EXPECT_EQ( c(i,5), d(i,5)+v[i]);

    for(std::size_t j=0; j < num_cols(c); ++j)
    {
        if(j != 5)
            for(std::size_t i = 0; i < num_rows(c); ++i)
                EXPECT_EQ(d(i,j),c(i,j));
    }

    cv += cv;

    for(std::size_t i = 0; i < num_rows(c); ++i)
        EXPECT_EQ( c(i,5), T(2)*(d(i,5)+v[i]));
}

TYPED_TEST(MatrixColumnView, minus_assign)
{
    using T = TypeParam;
    matrix<T> c(15,20);
    fill_matrix_with_numbers(c);
    matrix<T> d(c);
    column_view<matrix<T> > cv(c,5);
    vector<T> v(15);
    fill_range_with_numbers(v.begin(),v.end(),0);
    cv -= v;
    for(std::size_t i = 0; i < num_rows(c); ++i)
        EXPECT_EQ( c(i,5), d(i,5)-v[i]);

    for(std::size_t j=0; j < num_cols(c); ++j)
    {
        if(j != 5)
            for(std::size_t i = 0; i < num_rows(c); ++i)
                EXPECT_EQ(d(i,j),c(i,j));
    }

    cv -= cv;

    for(std::size_t i = 0; i < num_rows(c); ++i)
        EXPECT_EQ( c(i,5), T(0));
}

TYPED_TEST(MatrixColumnView, multiply_assign)
{
    using T = TypeParam;
    matrix<T> c(15,20);
    fill_matrix_with_numbers(c);
    matrix<T> d(c);
    column_view<matrix<T> > cv(c,19);

    cv *= T(5);

    for(std::size_t i=0; i < num_rows(c); ++i)
        EXPECT_EQ(c(i,19), T(5)*d(i,19));

    for(std::size_t j=0; j < num_cols(c); ++j)
    {
        if(j != 19)
            for(std::size_t i = 0; i < num_rows(c); ++i)
                EXPECT_EQ(c(i,j),d(i,j));
        else
            for(std::size_t i=0; i < num_rows(c); ++i)
                EXPECT_EQ(c(i,19), T(5)*d(i,19));
    }
}

TYPED_TEST(MatrixColumnView, scalar_product_column_view_vector)
{
    using T = TypeParam;
    using alps::numeric::conj;
    // We assume conj works properly
    ASSERT_EQ(conj(std::complex<double>(1,2)), std::complex<double>(1,-2));
    matrix<T> c(15,20);
    fill_matrix_with_numbers(c);
    matrix<T> d(c);
    column_view<matrix<T> > cv(c,5);
    vector<T> v(15);
    fill_range_with_numbers(v.begin(),v.end(),0);

    T r = scalar_product(cv, v);

    T ref(0);

    for(std::size_t i=0; i < v.size(); ++i)
        ref += conj(d(i,5))*v[i];

    EXPECT_EQ(ref, r);

    T r2 = scalar_product(v, cv);

    EXPECT_EQ(conj(ref) , r2);
    EXPECT_EQ(d,c);
}

TYPED_TEST(MatrixColumnView, scalar_product_column_view_column_view)
{
    using T = TypeParam;
    using alps::numeric::conj;
    // We assume conj works properly
    ASSERT_EQ(conj(std::complex<double>(1,2)), std::complex<double>(1,-2));
    matrix<T> c(15,20);
    fill_matrix_with_numbers(c);
    matrix<T> d(c);
    column_view<matrix<T> > cv(c,5);
    vector<T> v(15);
    fill_range_with_numbers(v.begin(),v.end(),0);

    T r = scalar_product(cv, cv);

    T ref(0);

    for(std::size_t i=0; i < v.size(); ++i)
        ref += conj(d(i,5))*d(i,5);

    EXPECT_EQ(ref, r);
    EXPECT_EQ(d,c);
}
