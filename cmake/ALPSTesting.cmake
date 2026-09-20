# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
include_guard(GLOBAL)

# Golden-output tests need stdin and byte comparison in addition to CTest's
# ordinary exit-status checks. Keep that behavior in one small runner.
function(alps_add_test name)
  if(NOT BUILD_TESTING)
    return()
  endif()
  cmake_parse_arguments(PARSE_ARGV 1 TEST "" "TARGET;INPUT;OUTPUT" "")
  if(TEST_UNPARSED_ARGUMENTS OR TEST_KEYWORDS_MISSING_VALUES)
    message(FATAL_ERROR "Invalid arguments to alps_add_test(${name})")
  endif()
  foreach(argument IN ITEMS TARGET INPUT OUTPUT)
    if(NOT DEFINED TEST_${argument})
      set(TEST_${argument} "${name}")
    endif()
  endforeach()
  add_test(
    NAME "${name}"
    COMMAND
      "${CMAKE_COMMAND}" "-Dname=${name}" "-Dcmd_path=$<TARGET_FILE:${TEST_TARGET}>"
      "-Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}" "-Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}"
      "-Dinput=${TEST_INPUT}" "-Doutput=${TEST_OUTPUT}" -P
      "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/run_test.cmake")
  set_tests_properties("${name}" PROPERTIES TIMEOUT 600)
endfunction()
