# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

# Eigen 5 retains the Eigen3 package name but rejects a version-3 request.
find_package(Eigen3 REQUIRED NO_MODULE)
if(Eigen3_VERSION VERSION_LESS 3.3)
  message(FATAL_ERROR "MaxEnt requires Eigen 3.3 or newer (found ${Eigen3_VERSION})")
endif()

add_library(alps_maxent STATIC
  src/maxent_run.cpp src/maxent_helper.cpp src/maxent_simulation.cpp
  src/maxent_params.cpp src/maxent_grid.cpp src/maxent_kernel.cpp
  src/maxent_backcont.cpp src/default_model.cpp)
add_library(ALPS::maxent ALIAS alps_maxent)
set_target_properties(alps_maxent PROPERTIES EXPORT_NAME maxent
  CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
target_link_libraries(alps_maxent PUBLIC ALPS::solver_headers ALPS::params ALPS::hdf5 ALPS::utilities ALPS::osiris
  Boost::filesystem Boost::date_time
  PRIVATE $<BUILD_INTERFACE:Eigen3::Eigen> ${ALPS_BLAS_TARGET} ${ALPS_LAPACK_TARGET})
target_compile_definitions(alps_maxent PRIVATE HAVE_LAPACK)
install(TARGETS alps_maxent EXPORT ALPSTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries)
install(FILES LICENSE.TXT ACKNOWLEDGE.TXT
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/licenses/maxent COMPONENT libraries)

add_executable(maxent cli/maxent.cpp)
target_link_libraries(maxent PRIVATE ALPS::maxent ALPS::cli)
set_target_properties(maxent PROPERTIES EXPORT_NAME maxent_executable)
set_property(TARGET alps APPEND PROPERTY ALPS_APPLICATION_TARGETS ALPS::maxent_executable)
install(TARGETS maxent EXPORT ALPSApplicationTargets RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
