#  Copyright Matthias Troyer and Synge Todo 2009 - 2015.
#   SPDX-License-Identifier: MIT

if(ALPS_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

install(
  DIRECTORY .
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/tutorials
  COMPONENT tutorials
  EXCLUDE_FROM_ALL
  FILES_MATCHING
  PATTERN "*.md"
  PATTERN "*.jpg"
  PATTERN "*.tex"
  PATTERN "*.bib"
  PATTERN "*.py"
  PATTERN "*.ipynb"
  PATTERN "*.sh"
  PATTERN "parm*"
  PATTERN "params*"
  PATTERN "*params"
  PATTERN "*.ip"
  PATTERN "*.op"
  PATTERN "*.op-*"
  PATTERN "*.dat"
  PATTERN "*.parm"
  PATTERN "*.xml"
  PATTERN "*input*"
  PATTERN "*.pvsm"
  PATTERN "*.h"
  PATTERN "*.hpp"
  PATTERN "*.ipp"
  PATTERN "*.c"
  PATTERN "*.cpp"
  PATTERN "*.C"
  PATTERN "*.f"
  PATTERN "*.f90"
  PATTERN "Makefile"
  PATTERN "CMakeLists.txt"
  PATTERN "*.h5" EXCLUDE
  PATTERN "parm*xml" EXCLUDE
  PATTERN "*.chkp" EXCLUDE
  PATTERN "*.xsl" EXCLUDE
  PATTERN "*.run*" EXCLUDE
  PATTERN "*clone*" EXCLUDE
  PATTERN "dump.xml" EXCLUDE)

# The examples share the golden-output runner with the repository tests.
install(
  FILES ../cmake/ALPSTesting.cmake ../cmake/run_test.cmake
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/cmake
  COMPONENT tutorials EXCLUDE_FROM_ALL)
# This is a fixed example input, unlike generated tutorial result archives.
install(
  FILES examples/alea/testfile.h5
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/tutorials/examples/alea
  COMPONENT tutorials EXCLUDE_FROM_ALL)
