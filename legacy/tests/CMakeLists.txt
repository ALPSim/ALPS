# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#
# Build-tree checks for the deprecated C++ interfaces. The installed-SDK
# consumer is exercised by tests/cmake/test_sdk_contracts.py.
add_executable(legacy_utilities_spellings utilities_spellings.cpp)
target_link_libraries(legacy_utilities_spellings PRIVATE ALPS::utilities)
alps_add_gtest(legacy_utilities_spellings COMPONENT utilities LABELS legacy)
