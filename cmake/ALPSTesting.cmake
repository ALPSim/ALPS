##########################################################################
# Regression Testing Support for Boost                                   #
##########################################################################
# Copyright (C) 2007-8 Douglas Gregor <doug.gregor@gmail.com>            #
# Copyright (C) 2007-8 Troy D. Straszheim                                #
#                                                                        #
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
##########################################################################
# This file provides a set of CMake macros that support regression
# testing for Boost libraries. For each of the test macros below, the
# first argument, testname, states the name of the test that will be
# created. If no other arguments are provided, the source file
# testname.cpp will be used as the source file; otherwise, source
# files should be listed immediately after the name of the test.
#
# The macros for creating regression tests are:
#   boost_test_run: Builds an executable and runs it as a test. The test
#                   succeeds if it builds and returns 0 when executed.
#
#   boost_test_run_fail: Builds an executable and runs it as a test. The
#                        test succeeds if it builds but returns a non-zero
#                        exit code when executed.
#  
#   boost_test_compile: Tests that the given source file compiles without 
#                       any errors.
#
#   boost_test_compile_fail: Tests that the given source file produces 
#                            errors when compiled.

# User-controlled option that can be used to enable/disable regression
# testing. By default, we ena testing, because most users building from source will
# want to check whether the buils is correct
include_guard(GLOBAL)

# Golden-output tests need stdin and byte comparison in addition to CTest's
# ordinary exit-status checks. Keep that behavior in one small runner.
function(alps_add_test name)
  if(NOT ALPS_BUILD_TESTING)
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
  if(ALPS_DATA_DIR)
    set(xml_resources "${ALPS_DATA_DIR}/xml")
  else()
    set(xml_resources "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../src/alps/resources")
  endif()
  set_tests_properties("${name}" PROPERTIES TIMEOUT 600
    ENVIRONMENT_MODIFICATION "ALPS_XML_PATH=set:${xml_resources}")
  if(WIN32)
    foreach(directory IN LISTS ALPS_RUNTIME_LIBRARY_DIRS)
      set_property(TEST "${name}" APPEND PROPERTY ENVIRONMENT_MODIFICATION
        "PATH=path_list_prepend:${directory}")
    endforeach()
  endif()
endfunction()
