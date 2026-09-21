#  Copyright Matthias Troyer 2009 - 2010.
#  Copyright ALPS collaboration 2026.
#   SPDX-License-Identifier: MIT

# Export applications with their configuration-specific executable locations.
function(alps_install_applications)
  install(
    TARGETS ${ARGN}
    EXPORT ALPSApplicationTargets
    RUNTIME_DEPENDENCY_SET alps_runtime
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT applications)
  list(TRANSFORM ARGN PREPEND "ALPS::")
  set_property(TARGET alps APPEND PROPERTY ALPS_APPLICATION_TARGETS ${ARGN})
endfunction()

add_subdirectory(mc)
add_subdirectory(qmc)
add_subdirectory(diag)
add_subdirectory(dmrg)
add_subdirectory(dmft)
add_subdirectory(maxent)

install(
  EXPORT ALPSApplicationTargets
  NAMESPACE ALPS::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps
  COMPONENT applications)
