#  Copyright Matthias Troyer and Synge Todo 2009-2013.
#   SPDX-License-Identifier: MIT

foreach(
  program IN
  ITEMS parameter2xml
        parameter2hdf5
        convert2xml
        compactrun
        lattice2xml
        pevaluate
        pconfig
        poutput
        printgraph
        snap2vtk)
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps)
  install(TARGETS ${program} RUNTIME_DEPENDENCY_SET alps_runtime
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
endforeach()
add_executable(p2h5 p2h5.cpp)
target_link_libraries(p2h5 PRIVATE ALPS::alps)
install(TARGETS p2h5 RUNTIME_DEPENDENCY_SET alps_runtime RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                                                                 COMPONENT tools)

if(UNIX)
  file(RELATIVE_PATH ALPS_XML_FROM_BINDIR
    "${CMAKE_INSTALL_FULL_BINDIR}" "${CMAKE_INSTALL_FULL_DATADIR}/alps/xml")
  configure_file(alps-xml.in alps-xml @ONLY NEWLINE_STYLE UNIX)
  install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/alps-xml"
    DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT xml)
endif()
