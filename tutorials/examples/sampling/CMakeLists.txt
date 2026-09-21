#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

#
# examples of MC error analysis
#

set(EXAMPLES fleas_correlated fleas_direct fleas_independent fleas_simpleminded fleas_uncorrelated)

foreach(program ${EXAMPLES})
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps)
  install(
    TARGETS ${program}
    COMPONENT examples
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach()
