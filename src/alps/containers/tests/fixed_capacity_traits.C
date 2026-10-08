/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2003 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/fixed_capacity_traits.h>
#include <alps/fixed_capacity_vector.h>
#include <alps/fixed_capacity_deque.h>
#include <gtest/gtest.h>
#include <list>
#include <queue>
#include <stack>
#include <vector>

template<class T> class DynamicCapacityTraits : public ::testing::Test {};
using DynamicTypes = ::testing::Types<std::vector<int>, std::list<int>,
    std::stack<int>, std::queue<int>, std::priority_queue<int>, double>;
TYPED_TEST_SUITE(DynamicCapacityTraits, DynamicTypes);
TYPED_TEST(DynamicCapacityTraits, CapacityIsNotFixed)
{
    EXPECT_FALSE(alps::fixed_capacity_traits<TypeParam>::capacity_is_fixed);
}

template<class Container, std::size_t Capacity> struct FixedCase {
    using container_type = Container;
    static constexpr std::size_t capacity = Capacity;
};
using FixedTypes = ::testing::Types<
    FixedCase<alps::fixed_capacity_vector<int, 8>, 8>,
    FixedCase<alps::fixed_capacity_deque<int, 8>, 8>,
    FixedCase<std::stack<int, alps::fixed_capacity_vector<int, 4>>, 4>,
    FixedCase<std::queue<int, alps::fixed_capacity_deque<int, 6>>, 6>,
    FixedCase<std::priority_queue<int, alps::fixed_capacity_vector<int, 16>>, 16>>;
template<class T> class FixedCapacityTraits : public ::testing::Test {};
TYPED_TEST_SUITE(FixedCapacityTraits, FixedTypes);
TYPED_TEST(FixedCapacityTraits, ExposesFixedMaximum)
{
    using Traits = alps::fixed_capacity_traits<typename TypeParam::container_type>;
    EXPECT_TRUE(Traits::capacity_is_fixed);
    EXPECT_EQ(Traits::static_max_size, TypeParam::capacity);
}
