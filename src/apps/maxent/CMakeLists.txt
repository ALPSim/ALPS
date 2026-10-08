# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

add_library(alps_maxent STATIC
  src/maxent_run.cpp src/maxent_helper.cpp src/maxent_simulation.cpp src/maxent_parms.cpp)
add_library(ALPS::maxent ALIAS alps_maxent)
set_target_properties(alps_maxent PROPERTIES EXPORT_NAME maxent
  CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
target_link_libraries(alps_maxent PUBLIC ALPS::solver_headers ALPS::params ALPS::hdf5 ALPS::utilities ALPS::osiris
  Boost::filesystem Boost::date_time ${ALPS_BLAS_TARGET} ${ALPS_LAPACK_TARGET})
install(TARGETS alps_maxent EXPORT ALPSTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries)

add_executable(maxent cli/maxent.cpp)
target_link_libraries(maxent PRIVATE ALPS::maxent ALPS::cli)
set_target_properties(maxent PROPERTIES EXPORT_NAME maxent_executable)
set_property(TARGET alps APPEND PROPERTY ALPS_APPLICATION_TARGETS ALPS::maxent_executable)
install(TARGETS maxent EXPORT ALPSApplicationTargets RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
