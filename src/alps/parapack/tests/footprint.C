/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/footprint.h>
#include <gtest/gtest.h>

class CustomObject {
public:
  CustomObject() : integer_(0), value_(0) {}
  std::size_t footprint() const { return sizeof(*this); }

private:
  int integer_;
  double value_;
};
TEST(Footprint, AccountsForValuesPointersAndAllocatedCapacity) {
  int value = 0;
  EXPECT_EQ(alps::footprint(value), sizeof(value));
  std::vector<int> vector(30);
  EXPECT_EQ(alps::footprint(vector), sizeof(vector) + vector.capacity() * sizeof(int));
  std::string text = "my string";
  EXPECT_EQ(alps::footprint(text), sizeof(text) + text.capacity());
  CustomObject object;
  EXPECT_EQ(alps::footprint(object), sizeof(object));
  double *pointer = nullptr;
  EXPECT_EQ(alps::footprint(pointer), sizeof(pointer));
}
