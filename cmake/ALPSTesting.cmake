##########################################################################
# ALPS Test Registration                                                 #
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
# alps_add_test(<name> [TARGET <target>] [INPUT <file>] [OUTPUT <file>])
#   Runs <target> (default: <name>) through run_test.cmake, feeding INPUT on
#   stdin and comparing stdout byte-for-byte with OUTPUT. Does nothing unless
#   ALPS_BUILD_TESTING is ON.
include_guard(GLOBAL)

function(alps_require_googletest)
  if(TARGET GTest::gtest_main)
    return()
  endif()
  # Prefer an installed provider for offline/HPC builds. A source override via
  # FETCHCONTENT_SOURCE_DIR_GOOGLETEST also avoids any configure-time download.
  # GLOBAL keeps the provider visible to sibling component directories.
  find_package(GTest 1.14 CONFIG QUIET GLOBAL)
  if(NOT TARGET GTest::gtest_main)
    include(FetchContent)
    set(INSTALL_GTEST OFF)
    set(BUILD_GMOCK OFF)
    set(gtest_force_shared_crt ON)
    FetchContent_Declare(googletest
      URL https://codeload.github.com/google/googletest/tar.gz/refs/tags/v1.18.0
      URL_HASH SHA256=6e3191c1455468b3fc35a417fb565c1c5071aee1b7e7f85e30cf48a98d37d8b5
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_MakeAvailable(googletest)
  endif()
endfunction()
# The installed tutorial examples also use this module's process runner. They
# should not acquire GoogleTest unless they actually register GoogleTest cases.
include(GoogleTest)

# Register real GoogleTest cases individually with CTest. Source files and SDK
# linkage remain ordinary CMake targets owned by their component.
function(alps_add_gtest target)
  if(NOT ALPS_BUILD_TESTING)
    return()
  endif()
  alps_require_googletest()
  cmake_parse_arguments(PARSE_ARGV 1 TEST "" "COMPONENT;TIMEOUT;THREADS" "LABELS;PROPERTIES")
  if(TEST_UNPARSED_ARGUMENTS OR TEST_KEYWORDS_MISSING_VALUES OR NOT TEST_COMPONENT)
    message(FATAL_ERROR "alps_add_gtest(${target}) requires COMPONENT and valid arguments")
  endif()
  if(NOT TEST_TIMEOUT)
    set(TEST_TIMEOUT 60)
  endif()
  if(NOT TEST_THREADS)
    set(TEST_THREADS 1)
  endif()
  target_link_libraries(${target} PRIVATE GTest::gtest_main)
  target_compile_features(${target} PRIVATE cxx_std_17)
  target_include_directories(${target} PRIVATE
    "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../tests/support")
  target_compile_definitions(${target} PRIVATE
    ALPS_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    ALPS_TEST_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
  set(workdir "${CMAKE_CURRENT_BINARY_DIR}/test-work/${target}")
  file(MAKE_DIRECTORY "${workdir}")
  set(environment
    "ALPS_XML_PATH=set:${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../src/alps/resources"
    "OMP_NUM_THREADS=set:${TEST_THREADS}" "OPENBLAS_NUM_THREADS=set:1")
  if(WIN32)
    foreach(directory IN LISTS ALPS_RUNTIME_LIBRARY_DIRS)
      list(APPEND environment "PATH=path_list_prepend:${directory}")
    endforeach()
  endif()
  set(labels "${TEST_COMPONENT};${TEST_LABELS}")
  if(NOT TEST_LABELS)
    list(APPEND labels unit)
  endif()
  list(FILTER labels EXCLUDE REGEX "^$")
  gtest_discover_tests(${target}
    TEST_PREFIX "${TEST_COMPONENT}.${target}." DISCOVERY_MODE PRE_TEST
    DISCOVERY_TIMEOUT 60 WORKING_DIRECTORY "${workdir}"
    PROPERTIES TIMEOUT "${TEST_TIMEOUT}" PROCESSORS "${TEST_THREADS}"
      ${TEST_PROPERTIES})
  # Set list-valued properties after discovery. CMake 3.27's GoogleTest
  # discovery script flattens semicolons passed through PROPERTIES.
  set(properties_file "${CMAKE_CURRENT_BINARY_DIR}/${target}_properties.cmake")
  get_property(multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
  if(multi_config)
    set(properties_config "\${CTEST_CONFIGURATION_TYPE}")
  else()
    set(properties_config "${CMAKE_BUILD_TYPE}")
  endif()
  file(WRITE "${properties_file}"
    "include(\"${CMAKE_CURRENT_BINARY_DIR}/${target}_properties-${properties_config}.cmake\")\n")
  string(CONCAT properties_content
    "if(EXISTS [==[$<TARGET_FILE:${target}>]==] AND NOT ${target}_TESTS)\n"
    "  message(FATAL_ERROR [==[${target} built successfully but discovered zero GoogleTest cases]==])\n"
    "endif()\n"
    "if(${target}_TESTS)\n"
    "  set_tests_properties(\${${target}_TESTS} PROPERTIES\n"
    "    LABELS [==[${labels}]==]\n"
    "    ENVIRONMENT_MODIFICATION [==[${environment}]==])\n"
    "endif()\n")
  file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${target}_properties-$<CONFIG>.cmake"
    CONTENT "${properties_content}")
  set_property(DIRECTORY APPEND PROPERTY TEST_INCLUDE_FILES "${properties_file}")
endfunction()

# MPI cases execute as a topology: every rank must enter the same collectives.
# Each rank writes a distinct GoogleTest report; CTest records the whole run.
function(alps_add_mpi_gtest target)
  if(NOT ALPS_BUILD_TESTING)
    return()
  endif()
  alps_require_googletest()
  cmake_parse_arguments(PARSE_ARGV 1 TEST "" "COMPONENT;RANKS;TIMEOUT" "LABELS")
  if(TEST_UNPARSED_ARGUMENTS OR NOT TEST_COMPONENT OR NOT TEST_RANKS)
    message(FATAL_ERROR "alps_add_mpi_gtest requires COMPONENT and RANKS")
  endif()
  if(NOT ALPS_ENABLE_MPI OR NOT MPIEXEC_EXECUTABLE)
    message(FATAL_ERROR "MPI tests require ALPS_ENABLE_MPI and an MPI launcher")
  endif()
  if(NOT TEST_TIMEOUT)
    set(TEST_TIMEOUT 120)
  endif()
  if(NOT TARGET alps_mpi_test_main)
    add_library(alps_mpi_test_main STATIC
      "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../tests/support/mpi_main.cpp")
    target_link_libraries(alps_mpi_test_main PUBLIC GTest::gtest MPI::MPI_CXX)
  endif()
  target_include_directories(${target} PRIVATE
    "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../tests/support")
  target_link_libraries(${target} PRIVATE alps_mpi_test_main)
  target_compile_features(${target} PRIVATE cxx_std_17)
  set(name "${TEST_COMPONENT}.${target}.np${TEST_RANKS}")
  set(workdir "${CMAKE_CURRENT_BINARY_DIR}/test-work/${name}")
  file(MAKE_DIRECTORY "${workdir}")
  add_test(NAME "${name}" COMMAND "${MPIEXEC_EXECUTABLE}"
    ${MPIEXEC_NUMPROC_FLAG} ${TEST_RANKS} ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
  set_tests_properties("${name}" PROPERTIES
    WORKING_DIRECTORY "${workdir}" TIMEOUT "${TEST_TIMEOUT}"
    PROCESSORS "${TEST_RANKS}" LABELS "${TEST_COMPONENT};mpi;integration;${TEST_LABELS}"
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endfunction()

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
  # Transitional text/CLI fixtures retain their own runner and label. Keep
  # generated files away from other tests and retain output on failure.
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-work/${name}")
  set_tests_properties("${name}" PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-work/${name}"
    LABELS "legacy" PROCESSORS 1)
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
