#  Copyright Matthias Troyer 2009 - 2010.
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

#
# libalps
#

add_library(alps)
add_library(ALPS::alps ALIAS alps)
add_library(alps_headers INTERFACE)
add_library(ALPS::headers ALIAS alps_headers)
set_target_properties(alps_headers PROPERTIES EXPORT_NAME headers)
# Foundations own their compile interfaces. The aggregate below is for the
# remaining simulation modules and for callers explicitly requesting all headers.
add_library(alps_configuration INTERFACE)
add_library(ALPS::configuration ALIAS alps_configuration)
set_target_properties(alps_configuration PROPERTIES EXPORT_NAME configuration)
target_compile_features(alps_configuration INTERFACE cxx_std_17)
target_compile_features(alps PRIVATE c_std_11)
target_link_libraries(alps_configuration INTERFACE Boost::headers)
target_include_directories(alps_configuration INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third_party/boost_numeric_bindings/include>)
target_link_libraries(alps PUBLIC alps_headers HDF5::HDF5 Threads::Threads
                                  ${ALPS_BLAS_TARGET} ${ALPS_LAPACK_TARGET})
foreach(component IN LISTS ALPS_BOOST_COMPONENTS)
  target_link_libraries(alps PUBLIC Boost::${component})
endforeach()
if(ALPS_ENABLE_MPI)
  target_link_libraries(alps PUBLIC MPI::MPI_CXX)
endif()
if(ALPS_ENABLE_OPENMP)
  target_link_libraries(alps PUBLIC OpenMP::OpenMP_CXX)
endif()
include(GenerateExportHeader)
generate_export_header(alps EXPORT_MACRO_NAME ALPS_DECL EXPORT_FILE_NAME
                       "${PROJECT_BINARY_DIR}/generated/include/alps/export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_configuration INTERFACE ALPS_STATIC_DEFINE)
endif()
set_target_properties(
  alps
  PROPERTIES VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}
             INSTALL_RPATH "${ALPS_LOADER_ORIGIN}")

function(alps_dependency_rpath target)
  if(NOT ALPS_LOADER_ORIGIN)
    return()
  endif()
  # CMake's automatic install RPATH omits dependency prefixes inside the source
  # or build tree. Imported shared libraries still need to be found after the
  # SDK is installed, including when CI extracts its dependencies in _build/.
  get_target_property(dependencies ${target} LINK_LIBRARIES)
  list(APPEND dependencies ${ARGN})
  foreach(dependency IN LISTS dependencies)
    if(TARGET "${dependency}")
      get_target_property(dependency_type "${dependency}" TYPE)
      get_target_property(dependency_imported "${dependency}" IMPORTED)
      if(dependency_imported AND dependency_type MATCHES "^(SHARED|UNKNOWN)_LIBRARY$")
        set_property(TARGET ${target} APPEND PROPERTY INSTALL_RPATH "$<TARGET_FILE_DIR:${dependency}>")
      endif()
    elseif(IS_ABSOLUTE "${dependency}")
      get_filename_component(dependency_directory "${dependency}" DIRECTORY)
      set_property(TARGET ${target} APPEND PROPERTY INSTALL_RPATH "${dependency_directory}")
    endif()
  endforeach()
endfunction()
# Keep paths for transitive libraries present on the core/program link lines.
get_target_property(hdf5_dependencies HDF5::HDF5 INTERFACE_LINK_LIBRARIES)
alps_dependency_rpath(alps ${hdf5_dependencies})

function(alps_copy_runtime_dependencies target)
  if(WIN32 AND BUILD_SHARED_LIBS)
    # Follow imported configurations: a Debug consumer may use a Release C-ABI
    # dependency that vcpkg's configuration-specific DLL search does not find.
    set(runtime_copy_script "${CMAKE_CURRENT_BINARY_DIR}/copy-runtime-$<CONFIG>.cmake")
    file(GENERATE OUTPUT "${runtime_copy_script}" CONTENT
      "set(dlls \"$<TARGET_RUNTIME_DLLS:${target}>\")\nif(dlls)\nfile(COPY \${dlls} DESTINATION \"$<TARGET_FILE_DIR:${target}>\")\nendif()\n")
    add_custom_command(TARGET ${target} POST_BUILD
      COMMAND "${CMAKE_COMMAND}" -P "${runtime_copy_script}"
      VERBATIM)
  endif()
endfunction()
alps_copy_runtime_dependencies(alps)

# Shared export/install policy; sources and dependencies stay with each component.
function(alps_configure_component target component)
  add_library(ALPS::${component} ALIAS ${target})
  set_target_properties(${target} PROPERTIES
    EXPORT_NAME ${component} VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR} INSTALL_RPATH "${ALPS_LOADER_ORIGIN}")
  alps_dependency_rpath(${target} ${ARGN})
  alps_copy_runtime_dependencies(${target})
  string(TOUPPER "${component}" upper)
  set(export_header "${PROJECT_BINARY_DIR}/generated/include/alps/${component}_export.h")
  generate_export_header(${target} EXPORT_MACRO_NAME ALPS_${upper}_DECL
    EXPORT_FILE_NAME "${export_header}")
  if(NOT BUILD_SHARED_LIBS)
    target_compile_definitions(${target} PUBLIC ALPS_${upper}_STATIC_DEFINE)
  endif()
  target_sources(${target} PUBLIC FILE_SET HEADERS
    BASE_DIRS "${PROJECT_BINARY_DIR}/generated/include" FILES "${export_header}")
endfunction()
target_sources(alps_configuration INTERFACE
  FILE_SET HEADERS
  BASE_DIRS "${PROJECT_BINARY_DIR}/generated/include"
  FILES "${PROJECT_BINARY_DIR}/generated/include/alps/config.h"
        "${PROJECT_BINARY_DIR}/generated/include/alps/version.h"
        "${PROJECT_BINARY_DIR}/generated/include/alps/export.h")

# Filesystem ownership follows scientific responsibilities. Some modules still
# contribute to libalps until their dependency cycles have been disentangled.
foreach(module IN ITEMS containers utilities cli numerics ietl solvers xml plotting
                        legacy_parameters expression graph lattice model random
                        osiris hdf5 numeric_io numeric_xml params alea accumulators mc scheduler parapack fortran)
  add_subdirectory(${module})
endforeach()
target_link_libraries(alps PUBLIC ALPS::utilities ALPS::osiris ALPS::hdf5 ALPS::params
                                  ALPS::xml ALPS::cli)

target_link_libraries(alps_headers INTERFACE alps_configuration
  "$<COMPILE_ONLY:ALPS::utilities;ALPS::osiris;ALPS::hdf5;ALPS::params;ALPS::xml;ALPS::cli;ALPS::numerics;ALPS::numeric_io;ALPS::numeric_xml;ALPS::solver_headers>")
if(ALPS_ENABLE_OPENMP)
  target_link_libraries(alps_headers INTERFACE "$<COMPILE_ONLY:OpenMP::OpenMP_CXX>")
endif()

# The one list of SDK targets. ALPSConfig.cmake derives its component names and
# runtime targets from it, so adding a module here is the only edit needed.
set(ALPS_SDK_TARGETS alps alps_utilities alps_hdf5 alps_params alps_osiris
                     alps_xml alps_cli alps_headers alps_configuration
                     alps_containers alps_numerics alps_numeric_io alps_numeric_xml alps_solver_headers)
set(sdk_components "")
set(sdk_runtime_targets "")
foreach(target IN LISTS ALPS_SDK_TARGETS ITEMS alps_fortran)
  get_target_property(component ${target} EXPORT_NAME)
  if(NOT component)
    set(component ${target})
  endif()
  list(APPEND sdk_components ${component})
  get_target_property(type ${target} TYPE)
  if(NOT type STREQUAL "INTERFACE_LIBRARY" AND NOT target STREQUAL "alps_fortran")
    list(APPEND sdk_runtime_targets ALPS::${component})
  endif()
endforeach()
set_property(GLOBAL PROPERTY ALPS_SDK_COMPONENTS ${sdk_components})
set_property(GLOBAL PROPERTY ALPS_SDK_RUNTIME_TARGETS ${sdk_runtime_targets})

# Install each target's actual public file sets, including generated exports.
# This keeps header installation policy identical for runtimes and interfaces.
foreach(target IN LISTS ALPS_SDK_TARGETS)
  get_target_property(header_sets ${target} INTERFACE_HEADER_SETS)
  set(header_install_args "")
  if(header_sets)
    foreach(header_set IN LISTS header_sets)
      list(APPEND header_install_args FILE_SET ${header_set}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers)
    endforeach()
  endif()
  install(TARGETS ${target} EXPORT ALPSTargets
    RUNTIME_DEPENDENCY_SET alps_runtime
    COMPONENT libraries
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ${header_install_args})
endforeach()

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/third_party/boost_numeric_bindings/include/boost"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT headers
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
  PATTERN "resources" EXCLUDE)

install(
  FILES "${PROJECT_SOURCE_DIR}/third_party/xdr/LICENSE"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/licenses/xdr
  COMPONENT libraries)

install(
  FILES "${PROJECT_SOURCE_DIR}/third_party/boost_numeric_bindings/LICENSE_1_0.txt"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/licenses/boost_numeric_bindings
  COMPONENT headers)

install(
  DIRECTORY resources/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/xml
  COMPONENT xml
  FILES_MATCHING PATTERN "*.xml" PATTERN "*.xsl")
