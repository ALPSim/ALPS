# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT

# The SDK exports the same application targets that the wheel bundles. Target
# locations preserve custom install directories and executable suffixes.
function(alps_install_applications)
  install(
    TARGETS ${ARGN}
    EXPORT ALPSApplicationTargets
    RUNTIME_DEPENDENCY_SET alps_runtime
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT applications)
  foreach(application IN LISTS ARGN)
    set_property(
      TARGET alps
      APPEND
      PROPERTY ALPS_APPLICATION_TARGETS "ALPS::${application}")
  endforeach()
endfunction()
