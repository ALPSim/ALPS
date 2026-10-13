# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#
# Build-tree checks for the deprecated C++ interfaces. The installed-SDK
# consumer is exercised by tests/cmake/test_sdk_contracts.py.
foreach(name utilities_spellings retired_headers bitops)
  add_executable(legacy_${name} ${name}.cpp)
  target_link_libraries(legacy_${name} PRIVATE ALPS::utilities)
  alps_add_gtest(legacy_${name} COMPONENT utilities LABELS legacy)
endforeach()

add_executable(legacy_module_spellings module_spellings.cpp)
target_link_libraries(legacy_module_spellings PRIVATE ALPS::cli ALPS::numerics)
alps_add_gtest(legacy_module_spellings COMPONENT cli LABELS legacy)
