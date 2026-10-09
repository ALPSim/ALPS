/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/integer_range.h>
#include <gtest/gtest.h>
#include <limits>

// Includes signed and unsigned 32/64-bit domains from the original test.
template <class T> class IntegerRangeBoundaries : public ::testing::Test {};
using IntegerTypes = ::testing::Types<int, unsigned int, long long, unsigned long long>;
TYPED_TEST_SUITE(IntegerRangeBoundaries, IntegerTypes);
TYPED_TEST(IntegerRangeBoundaries, RejectsUnrepresentableSizeWithoutOverflow) {
  using T = TypeParam;
  const T low = (std::numeric_limits<T>::min)(), high = (std::numeric_limits<T>::max)();
  const alps::integer_range<T> full(low, high), empty;
  EXPECT_FALSE(full.empty());
  EXPECT_TRUE(full.valid());
  EXPECT_TRUE(empty.empty());
  EXPECT_EQ(empty.size(), 0);
  EXPECT_EQ(alps::integer_range<T>(low).size(), 1);
  EXPECT_EQ(alps::integer_range<T>(high).size(), 1);
  EXPECT_EQ(alps::integer_range<T>(0, high - 1).size(), high);
  EXPECT_EQ(unify(full, full), full);
  EXPECT_THROW(full.size(), std::overflow_error);
}

TEST(IntegerRange, ParsesSignedEndpointsAndExpressions) {
  alps::Parameters params;
  params["L"] = 4;
  params["T"] = .5;
  const int low = (std::numeric_limits<int>::min)(), high = (std::numeric_limits<int>::max)();
  struct Sample {
    const char *input;
    int first, last;
  };
  const Sample samples[] = {{"[0:3]", 0, 3},      {"[3:9]", 3, 9}, {"[:3]", low, 3},
                            {"[-1:]", -1, high},  {"[3]", 3, 3},   {"6", 6, 6},
                            {"[1/T:L^2]", 2, 16}, {"[1/T]", 2, 2}, {"[L:]", 4, high},
                            {"[:L]", low, 4},     {"[2*L]", 8, 8}, {"exp(L)", 54, 54}};
  for (const auto &sample : samples) {
    SCOPED_TRACE(sample.input);
    const alps::integer_range<int> range(sample.input, params);
    EXPECT_EQ(range.min(), sample.first);
    EXPECT_EQ(range.max(), sample.last);
  }
  EXPECT_TRUE(alps::integer_range<int>("[]", params).empty());
  EXPECT_THROW(alps::integer_range<int>("[3:9:0]", params), std::exception);
}
TEST(IntegerRange, ParsesUnsignedEndpointsAndRejectsNegativeBounds) {
  using Range = alps::integer_range<unsigned int>;
  EXPECT_EQ(Range("[0:3]"), Range(0, 3));
  EXPECT_EQ(Range("[3:9]"), Range(3, 9));
  EXPECT_EQ(Range("[:3]"), Range(0, 3));
  EXPECT_EQ(Range("[1:]"), Range(1, (std::numeric_limits<unsigned int>::max)()));
  EXPECT_EQ(Range("6"), Range(6));
  EXPECT_THROW(Range("[-3:9]"), std::exception);
}
TEST(IntegerRange, AssignmentInclusionAndSaturatingScaling) {
  using Range = alps::integer_range<int>;
  Range range(0, 5);
  EXPECT_FALSE(range.is_included(7));
  range = 3;
  EXPECT_EQ(range, Range(3));
  range.include(8);
  EXPECT_EQ(range, Range(3, 8));
  EXPECT_TRUE(range.is_included(7));
  Range unbounded("[3:]");
  const int high = (std::numeric_limits<int>::max)();
  EXPECT_EQ(3.5 * unbounded, Range(10, high));
  EXPECT_EQ(unbounded * 2000000000, Range(high, high));
  unbounded *= .1;
  EXPECT_EQ(unbounded, Range(0, static_cast<int>(high * .1)));
}
TEST(IntegerRange, OverlapAndUnion) {
  using Range = alps::integer_range<int>;
  const Range range(0, 5), partial(-2, 3), full(-2, 10), disjoint(7, 10);
  EXPECT_EQ(overlap(range, partial), Range(0, 3));
  EXPECT_EQ(unify(range, partial), Range(-2, 5));
  EXPECT_EQ(overlap(range, full), range);
  EXPECT_EQ(unify(range, full), full);
  EXPECT_TRUE(overlap(range, disjoint).empty());
  EXPECT_THROW(unify(range, disjoint), std::exception);
}
