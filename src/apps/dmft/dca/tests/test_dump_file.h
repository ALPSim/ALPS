// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

#pragma once
#include <cctype>
#include <string>
#include "gtest/gtest.h"

// CTest runs every test case as its own process, possibly in parallel, in one
// working directory. Give each case its own cluster dump file (dca.LATTICE_DUMP).
inline std::string test_dump_file() {
  const ::testing::TestInfo* info = ::testing::UnitTest::GetInstance()->current_test_info();
  std::string name = std::string(info->test_suite_name()) + "." + info->name();
  for (char& c : name)
    if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';
  return "./cluster_" + name + ".h5";
}
