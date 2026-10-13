# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#
# Deprecated compatibility layer for projects written against ALPS 2.x:
#
#   find_package(ALPS REQUIRED)
#   include(${ALPS_USE_FILE})
#   target_link_libraries(app ${ALPS_LIBRARIES})
#   add_alps_test(app)
#
# The imported ALPS::alps target carries the include directories, compile
# settings and dependencies. Unlike the 2.x file, this one changes no compilers,
# global flags, include/link directories or cache entries. New projects link
# ALPS::alps directly.

include_guard(GLOBAL)
if(NOT TARGET ALPS::alps)
  message(FATAL_ERROR "Include \${ALPS_USE_FILE} after find_package(ALPS)")
endif()
message(DEPRECATION "UseALPS.cmake is deprecated: link ALPS::alps instead of \${ALPS_LIBRARIES}")

cmake_path(GET CMAKE_CURRENT_LIST_DIR PARENT_PATH _alps_datadir)
cmake_path(GET _alps_datadir PARENT_PATH ALPS_ROOT_DIR)
set(ALPS_INCLUDE_DIRS "${ALPS_ROOT_DIR}/include")
set(ALPS_LIBRARY_DIRS "")
foreach(_alps_config IN ITEMS "" _RELEASE _RELWITHDEBINFO _MINSIZEREL _DEBUG)
  get_target_property(_alps_location ALPS::alps IMPORTED_LOCATION${_alps_config})
  if(_alps_location)
    cmake_path(GET _alps_location PARENT_PATH ALPS_LIBRARY_DIRS)
    break()
  endif()
endforeach()
set(ALPS_LIBRARIES ALPS::alps)
set(ALPS_FORTRAN_LIBRARIES ALPS::fortran)
# Carried by ALPS::alps; kept empty for scripts that still expand them.
set(ALPS_EXTRA_DEFINITIONS "")
set(ALPS_EXTRA_INCLUDE_DIRS "")
set(ALPS_EXTRA_LIBRARIES "")
set(ALPS_EXTRA_LINKER_FLAGS "")

# Enable testing in the including directory; function scope does not generate
# CTestTestfile.cmake for legacy projects that rely on add_alps_test alone.
enable_testing()

# add_alps_test(name), (name io), (name input output) or (name cmd input output):
# run <cmd> with <input>.input on stdin and compare stdout with <output>.output.
function(add_alps_test name)
  set(cmd "${name}")
  set(input "${name}")
  set(output "${name}")
  if(ARGC EQUAL 2)
    set(input "${ARGV1}")
    set(output "${ARGV1}")
  elseif(ARGC EQUAL 3)
    set(input "${ARGV1}")
    set(output "${ARGV2}")
  elseif(ARGC GREATER 3)
    set(cmd "${ARGV1}")
    set(input "${ARGV2}")
    set(output "${ARGV3}")
  endif()
  add_test(
    NAME "${name}"
    COMMAND
      "${CMAKE_COMMAND}" "-Dname=${name}" "-Dcmd_path=$<TARGET_FILE:${cmd}>"
      "-Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}" "-Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}"
      "-Dinput=${input}" "-Doutput=${output}" -P
      "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/run_test.cmake")
  set_tests_properties("${name}" PROPERTIES TIMEOUT 600
    ENVIRONMENT_MODIFICATION "ALPS_XML_PATH=set:${ALPS_DATA_DIR}/xml")
  if(WIN32)
    foreach(directory IN LISTS ALPS_RUNTIME_LIBRARY_DIRS)
      set_property(TEST "${name}" APPEND PROPERTY ENVIRONMENT_MODIFICATION
        "PATH=path_list_prepend:${directory}")
    endforeach()
  endif()
endfunction()
