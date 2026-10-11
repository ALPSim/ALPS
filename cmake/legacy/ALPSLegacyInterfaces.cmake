# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#
# Deprecated ALPS 2.x entry points, installed next to ALPSConfig.cmake:
# UseALPS.cmake, UsePyALPS.cmake and the Makefile fragment include.mk.
# find_package(ALPS) only records their locations; projects opt in by
# including them. Everything for these interfaces lives in this directory and
# tests/cmake/legacy/; removing them later means deleting both, the include in
# the top-level CMakeLists.txt and the two USE_FILE lines in ALPSConfig.cmake.in.
include_guard(GLOBAL)

install(
  FILES "${CMAKE_CURRENT_LIST_DIR}/UseALPS.cmake" "${CMAKE_CURRENT_LIST_DIR}/UsePyALPS.cmake"
        "${PROJECT_SOURCE_DIR}/cmake/run_test.cmake"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps
  COMPONENT build)

if(WIN32)
  return()
endif()

# include.mk needs the external libraries as plain linker arguments. Flatten
# the usage requirements of the exported runtime targets into flags.
set(_alps_mk_cppflags "")
set(_alps_mk_libs "")
set(_alps_mk_seen "")
function(_alps_mk_add_item item)
  if(TARGET "${item}")
    get_target_property(aliased "${item}" ALIASED_TARGET)
    if(aliased)
      set(item "${aliased}")
    endif()
    if(item IN_LIST _alps_mk_seen)
      return()
    endif()
    list(APPEND _alps_mk_seen "${item}")
    get_target_property(imported "${item}" IMPORTED)
    get_target_property(type "${item}" TYPE)
    if(imported)
      get_target_property(includes "${item}" INTERFACE_INCLUDE_DIRECTORIES)
      foreach(directory IN LISTS includes)
        if(directory AND NOT directory MATCHES "^\\$<")
          list(APPEND _alps_mk_cppflags "-I${directory}")
        endif()
      endforeach()
      get_target_property(definitions "${item}" INTERFACE_COMPILE_DEFINITIONS)
      foreach(definition IN LISTS definitions)
        if(definition AND NOT definition MATCHES "^\\$<")
          list(APPEND _alps_mk_cppflags "-D${definition}")
        endif()
      endforeach()
      if(NOT type STREQUAL "INTERFACE_LIBRARY")
        set(location "")
        string(TOUPPER "${CMAKE_BUILD_TYPE}" config)
        foreach(suffix IN ITEMS "_${config}" "" _RELEASE _RELWITHDEBINFO _MINSIZEREL _DEBUG _NOCONFIG)
          get_target_property(location "${item}" IMPORTED_LOCATION${suffix})
          if(location)
            break()
          endif()
        endforeach()
        if(location)
          list(APPEND _alps_mk_libs "${location}")
        endif()
      endif()
    endif()
    get_target_property(links "${item}" INTERFACE_LINK_LIBRARIES)
    # Keep static link dependencies; drop build-tree, compile-only and
    # configuration-specific entries before splitting the list.
    string(REGEX REPLACE "\\$<LINK_ONLY:([^<>]*)>" "\\1" links "${links}")
    string(GENEX_STRIP "${links}" links)
    foreach(link IN LISTS links)
      if(link)
        _alps_mk_add_item("${link}")
      endif()
    endforeach()
  elseif(item MATCHES "^(.*/)?([^/]+)\\.framework$")
    list(APPEND _alps_mk_libs "-framework ${CMAKE_MATCH_2}")
  elseif(IS_ABSOLUTE "${item}" OR item MATCHES "^-")
    list(APPEND _alps_mk_libs "${item}")
  else()
    list(APPEND _alps_mk_libs "-l${item}")
  endif()
  set(_alps_mk_seen "${_alps_mk_seen}" PARENT_SCOPE)
  set(_alps_mk_cppflags "${_alps_mk_cppflags}" PARENT_SCOPE)
  set(_alps_mk_libs "${_alps_mk_libs}" PARENT_SCOPE)
endfunction()

set(_alps_mk_alps_libs "")
foreach(target IN LISTS ALPS_SDK_RUNTIME_TARGETS)
  get_target_property(name "${target}" OUTPUT_NAME)
  if(NOT name)
    get_target_property(name "${target}" ALIASED_TARGET)
  endif()
  list(APPEND _alps_mk_alps_libs "-l${name}")
  _alps_mk_add_item("${target}")
endforeach()
list(REMOVE_DUPLICATES _alps_mk_cppflags)
list(REMOVE_DUPLICATES _alps_mk_libs)
if(NOT BUILD_SHARED_LIBS)
  list(APPEND _alps_mk_cppflags -DALPS_STATIC_DEFINE)
endif()
list(JOIN _alps_mk_cppflags " " ALPS_MK_CPPFLAGS)
list(JOIN _alps_mk_alps_libs " " ALPS_MK_ALPS_LIBS)
list(JOIN _alps_mk_libs " " ALPS_MK_LIBS)
string(TOUPPER "${CMAKE_BUILD_TYPE}" _alps_mk_config)
set(ALPS_MK_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_alps_mk_config}} ${CMAKE_CXX17_STANDARD_COMPILE_OPTION}")
if(ALPS_ENABLE_OPENMP)
  string(APPEND ALPS_MK_CXXFLAGS " ${OpenMP_CXX_FLAGS}")
endif()
configure_file("${CMAKE_CURRENT_LIST_DIR}/include.mk.in" "${PROJECT_BINARY_DIR}/cmake/include.mk" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/cmake/include.mk"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps
  COMPONENT build)
