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

add_subdirectory(alea)

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
        txt2archive
        xml2archive
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

if(ALPS_BUILD_ARCHIVE)
  find_package(SQLite3 REQUIRED)
  add_executable(archive archive.cpp archive_index.cpp archive_node.cpp archive_plot.cpp
                         archive_sqlite.cpp archive_xml.cpp)
  target_link_libraries(archive PRIVATE ALPS::alps SQLite::SQLite3)
  install(
    TARGETS
    archive RUNTIME_DEPENDENCY_SET alps_runtime
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
endif()

if(UNIX)
  install(
    PROGRAMS archivecat
             convert2html
             convert2text
             extractgp
             extracthtml
             extractmpl
             extracttext
             extractxmgr
             plot2gp
             plot2html
             plot2mpl
             plot2text
             plot2xmgr
             use_local_stylesheet
             xslttransform
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT tools)
endif()

include(maxent_sources.cmake)
add_library(maxent_impl STATIC ${ALPS_MAXENT_SOURCES})
target_link_libraries(maxent_impl PUBLIC ALPS::alps)
add_executable(maxent maxent.cpp)
target_link_libraries(maxent PRIVATE maxent_impl)
install(TARGETS maxent RUNTIME_DEPENDENCY_SET alps_runtime
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
if(BUILD_TESTING)
  add_executable(maxent_linear_grid_numeric maxent_linear_grid_numeric.cpp)
  target_link_libraries(maxent_linear_grid_numeric PRIVATE maxent_impl)
  alps_add_test(maxent_linear_grid_numeric)
endif()
