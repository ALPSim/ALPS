# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

add_library(alps_maxent STATIC maxent_run.cpp maxent_helper.cpp maxent_simulation.cpp maxent_parms.cpp)
add_library(ALPS::maxent ALIAS alps_maxent)
set_target_properties(alps_maxent PROPERTIES EXPORT_NAME maxent
  CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
target_link_libraries(alps_maxent PUBLIC ALPS::alps)
install(TARGETS alps_maxent EXPORT ALPSTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries)

add_executable(maxent maxent.cpp)
target_link_libraries(maxent PRIVATE ALPS::maxent)
install(TARGETS maxent RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
