/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/util.h>
#include <gtest/gtest.h>
#include <ostream>

namespace {
struct IdentifierSpelling {
  int id;
  const char* spelling;
};
void PrintTo(const IdentifierSpelling& value, std::ostream* out) {
  *out << value.id;
}
class ProcessIdentifier : public ::testing::TestWithParam<IdentifierSpelling> {};
TEST_P(ProcessIdentifier, KeepsHistoricalSortableSpelling) {
  EXPECT_EQ(alps::id2string(GetParam().id), GetParam().spelling);
}
INSTANTIATE_TEST_SUITE_P(
    HistoricalValues, ProcessIdentifier,
    ::testing::Values(IdentifierSpelling{0, "0"}, IdentifierSpelling{9, "9"},
                      IdentifierSpelling{10, "_10"}, IdentifierSpelling{15, "_15"},
                      IdentifierSpelling{100, "__100"}, IdentifierSpelling{200, "__200"},
                      IdentifierSpelling{1000, "___1000"}, IdentifierSpelling{1001, "___1001"},
                      IdentifierSpelling{100000, "_____100000"}));
} // namespace
