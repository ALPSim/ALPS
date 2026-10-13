/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012        by Michele Dolfi <dolfim@phys.ethz.ch>                *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <iterator>
#include <complex>
#include <vector>
#include <alps/numeric/real.hpp>

#include <gtest/gtest.h>
#include <algorithm>

//
// List of types T for which the real(T) is tested
//
typedef ::testing::Types<
      float
    , double
    , std::complex<float>
    , std::complex<double>
    , std::vector<float>
    , std::vector<double>
    , std::vector<std::complex<float> >
    , std::vector<std::complex<double> >
    , std::vector<std::vector<float> >
    , std::vector<std::vector<double> >
    , std::vector<std::vector<std::complex<float> > >
    , std::vector<std::vector<std::complex<double> > >
> test_types;


//
// Filling functions
//
float fill_val = 1;
const std::size_t vecsize = 2;
template <typename T>
void fill (T & v)
{
    v = (fill_val++);
}
template <typename T>
void fill (std::complex<T> & v) {
    T real_part = fill_val++;
    T imag_part = fill_val++;
    v = std::complex<T>(real_part,imag_part);
}
template <typename T>
void fill (std::vector<T> & v)
{
    v.resize(vecsize);
    std::for_each(v.begin(), v.end(), static_cast<void (*)(T &)>(&fill));
}



// Independent scalar/recursive reference, including nested vectors.
template<class T> T reference_real(const T& value) { return value; }
template<class T> T reference_real(const std::complex<T>& value) { return value.real(); }
template<class T>
auto reference_real(const std::vector<T>& values)
{
    using result_type = decltype(reference_real(T{}));
    std::vector<result_type> result;
    for (const auto& value : values) result.push_back(reference_real(value));
    return result;
}

template<class T> class Real : public ::testing::Test {
protected:
    void SetUp() override { fill_val = 1; }
};
TYPED_TEST_SUITE(Real, test_types);

TYPED_TEST(Real, QualifiedLookup)
{
    TypeParam value;
    fill(value);
    EXPECT_EQ(alps::numeric::real(value), reference_real(value));
}

TYPED_TEST(Real, UnqualifiedLookup)
{
    TypeParam value;
    fill(value);
    using alps::numeric::real;
    EXPECT_EQ(real(value), reference_real(value));
}
