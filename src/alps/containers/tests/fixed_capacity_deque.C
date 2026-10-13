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

#include <alps/fixed_capacity_deque.h>
#include <gtest/gtest.h>
#include <list>
#include <vector>

namespace {
using Container = alps::fixed_capacity_deque<float, 6>;
std::vector<float> elements(const Container& value)
{
    return {value.begin(), value.end()};
}
}

TEST(FixedCapacityDeque, ModifiersAndCapacity)
{
    Container a(3);
    EXPECT_EQ(a.size(), 3u);
    a.assign(2);
    EXPECT_EQ(elements(a), (std::vector<float>{2, 2, 2}));
    a.push_back(4);
    a.erase(a.begin() + 1);
    a.insert(a.begin() + 2, 5);
    a.insert(a.begin() + 1, 2, 6);
    a.pop_back();
    EXPECT_EQ(a.size(), 5u);
    EXPECT_FALSE(a.empty());
    EXPECT_EQ(a.max_size(), 6u);
    EXPECT_EQ(a.front(), 2);
    EXPECT_EQ(a.back(), 5);
    EXPECT_EQ(elements(a), (std::vector<float>{2, 6, 6, 2, 5}));
    a.pop_back();
    a.pop_back();
    a.push_front(1);
    a.push_front(1);
    a.push_front(1);
    EXPECT_EQ(elements(a), (std::vector<float>{1, 1, 1, 2, 6, 6}));
    a.erase(a.begin(), a.end());
    EXPECT_TRUE(a.empty());
    EXPECT_EQ(a.begin(), a.end());
    a.resize(2, 1);
    const std::vector<float> vector{2, 3, 4, 5};
    a.insert(a.begin() + 1, vector.begin(), vector.end());
    EXPECT_EQ(elements(a), (std::vector<float>{1, 2, 3, 4, 5, 1}));
    a.erase(a.begin() + 3, a.begin() + 6);
    const std::list<float> list{6, 7, 8};
    a.insert(a.begin(), list.begin(), list.end());
    // Characterize the historical forward-iterator insertion order, recorded
    // in the former stdout fixture. Changing it is a separate compatibility decision.
    EXPECT_EQ(elements(a), (std::vector<float>{8, 7, 6, 1, 2, 3}));
    EXPECT_EQ((std::vector<float>{a.rbegin(), a.rend()}),
              (std::vector<float>{3, 2, 1, 6, 7, 8}));
}

TEST(FixedCapacityDeque, CopiesAndSwapsDifferentSizes)
{
    const std::vector<float> values{8, 7, 6, 1, 2, 3};
    Container a(values.begin(), values.end());
    Container b(a);
    Container c;
    c = a;
    EXPECT_EQ(a, b);
    EXPECT_EQ(a, c);
    Container d(a);
    Container e(2);
    e = a;
    EXPECT_EQ(a, d);
    EXPECT_EQ(a, e);
    b.clear();
    b.push_front(12);
    b.push_front(11);
    b.push_front(10);
    EXPECT_EQ(elements(a), values);
    EXPECT_EQ(elements(b), (std::vector<float>{10, 11, 12}));
    swap(a, b);
    EXPECT_EQ(elements(a), (std::vector<float>{10, 11, 12}));
    EXPECT_EQ(elements(b), values);
    swap(a, b);
    EXPECT_EQ(elements(a), values);
    EXPECT_EQ(elements(b), (std::vector<float>{10, 11, 12}));
}
