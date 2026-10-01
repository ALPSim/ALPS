#  Copyright Matthias Troyer and Synge Todo 2009-2013.
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
# basic C++ tools
#


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
